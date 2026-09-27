# AC-OPF as a complex polynomial optimization problem, following Section 4 of Josz & Molzahn.
#
#   min  sum_k c_k * (v^H H_k v + pd_k)
#   s.t. pmin_k - pd_k <= v^H H_k  v <= pmax_k - pd_k      (real power injection)
#        qmin_k - qd_k <= v^H H~_k v <= qmax_k - qd_k      (reactive power injection)
#        vmin_k^2      <= v^H e_k e_k^T v <= vmax_k^2      (voltage magnitude)
#
# with S_k = v_k conj((Yv)_k) = v^H (Y^H e_k e_k^T) v, so H_k and H~_k are the Hermitian and
# skew-Hermitian parts of Y^H e_k e_k^T:
#
#   H_k = (Y^H e_k e_k^T + e_k e_k^T Y) / 2,   H~_k = (Y^H e_k e_k^T - e_k e_k^T Y) / (2i).
#
# Every one of these is v^H A v with A Hermitian, i.e. a CPoly whose terms all have |alpha| =
# |beta| = 1. So the problem is T-INVARIANT and the order-1 relaxation block-diagonalises into
# [1] and the n x n matrix W = v v^H -- which IS the standard SDP relaxation of AC-OPF. That is
# the correctness gate.
#
# LINEAR COSTS ONLY, deliberately. A quadratic cost in p is quartic in v and the paper handles it
# by introducing real epigraph variables t_k, which mixes real variables into a complex POP and
# breaks T-invariance. Nothing about the hierarchy prevents that, but it changes the order-1
# structure and so cannot be part of the correctness gate; it is the next step, not this one.

"Bus admittance matrix, including transformer taps, line charging and bus shunts."
function complex_ybus(data::Dict{String,Any})
    ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]
    buses = sort(collect(keys(ref[:bus])))
    idx = Dict(b => i for (i, b) in enumerate(buses))
    n = length(buses)
    Y = zeros(ComplexF64, n, n)
    for (_, br) in ref[:branch]
        f, t = idx[br["f_bus"]], idx[br["t_bus"]]
        gs, bs = PowerModels.calc_branch_y(br)
        trr, tii = PowerModels.calc_branch_t(br)
        y = complex(gs, bs); T = complex(trr, tii); tm2 = abs2(T)
        Y[f, f] += (y + complex(br["g_fr"], br["b_fr"])) / tm2
        Y[f, t] += -y / conj(T)
        Y[t, f] += -y / T
        Y[t, t] += y + complex(br["g_to"], br["b_to"])
    end
    for (_, sh) in ref[:shunt]
        k = idx[sh["shunt_bus"]]
        Y[k, k] += complex(sh["gs"], sh["bs"])
    end
    return Y, buses, idx
end

"v^H A v as a CPoly over variables 1..n."
function quadform(A::AbstractMatrix{ComplexF64})
    p = CPoly()
    for i in axes(A, 1), j in axes(A, 2)
        A[i, j] == 0 && continue
        caddterm!(p, ([i], [j]), A[i, j])
    end
    return p
end

"""
    build_complex_power_pop(data) -> CPOP

AC-OPF in complex voltages. Requires at most one in-service generator per bus and linear
generation costs; both are checked, because silently aggregating generators or dropping a
quadratic cost term would make the comparison against the real hierarchy meaningless.
"""
function build_complex_power_pop(data::Dict{String,Any})
    ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]
    Y, buses, idx = complex_ybus(data)
    n = length(buses)
    pop = CPOP()
    for b in buses
        bus = ref[:bus][b]
        add_cvar!(pop, "v[$b]"; absmin = Float64(bus["vmin"]), absmax = Float64(bus["vmax"]))
    end
    gens = Dict{Int,Any}()
    for (_, g) in ref[:gen]
        Int(get(g, "gen_status", 1)) == 0 && continue
        k = Int(g["gen_bus"])
        haskey(gens, k) && error("bus $k has more than one in-service generator; the complex " *
                                 "formulation here writes generation as the bus injection")
        length(g["cost"]) <= 2 || (g["cost"][1] == 0 ||
            error("generator at bus $k has a quadratic cost; the complex POP supports linear " *
                  "costs only (see the note in src/cpower.jl)"))
        gens[k] = g
    end
    pd = Dict(b => 0.0 for b in buses); qd = Dict(b => 0.0 for b in buses)
    for (_, l) in ref[:load]
        pd[Int(l["load_bus"])] += Float64(l["pd"]); qd[Int(l["load_bus"])] += Float64(l["qd"])
    end

    obj = CPoly()
    for b in buses
        k = idx[b]
        ek = zeros(ComplexF64, n, n); ek[k, k] = 1
        YHe = adjoint(Y) * ekcol(n, k)
        Hk = (YHe + adjoint(YHe)) / 2
        Hqk = (YHe - adjoint(YHe)) / (2im)
        Pk = quadform(Hk); Qk = quadform(Hqk); Wk = quadform(ek)
        bus = ref[:bus][b]
        add_cineq!(pop, Wk - (Float64(bus["vmin"])^2), "vmin[$b]")
        add_cineq!(pop, (Float64(bus["vmax"])^2) - Wk, "vmax[$b]")
        if haskey(gens, b)
            g = gens[b]
            add_cineq!(pop, Pk - (Float64(g["pmin"]) - pd[b]), "pg_min[$b]")
            add_cineq!(pop, (Float64(g["pmax"]) - pd[b]) - Pk, "pg_max[$b]")
            add_cineq!(pop, Qk - (Float64(g["qmin"]) - qd[b]), "qg_min[$b]")
            add_cineq!(pop, (Float64(g["qmax"]) - qd[b]) - Qk, "qg_max[$b]")
            c = g["cost"]
            c1 = length(c) >= 2 ? Float64(c[end-1]) : 0.0
            obj = obj + c1 * (Pk + pd[b])
        else
            # no generator: the injection is exactly minus the load
            add_ceq!(pop, Pk + pd[b], "p_balance[$b]")
            add_ceq!(pop, Qk + qd[b], "q_balance[$b]")
        end
    end
    pop.obj = obj
    pop.meta["buses"] = buses
    pop.meta["bus_index"] = idx
    pop.meta["const_cost"] = sum((length(g["cost"]) >= 1 ? Float64(g["cost"][end]) : 0.0)
                                 for g in values(gens); init = 0.0)
    return pop
end

"n x n matrix with a single 1 in column k of row k, i.e. e_k e_k^T, as a full matrix product helper."
ekcol(n::Int, k::Int) = (M = zeros(ComplexF64, n, n); M[k, k] = 1; M)
