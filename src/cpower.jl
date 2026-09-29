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

AC-OPF in complex voltages.

GENERATORS. A bus with zero or one generator needs no auxiliary variable: its generation IS the
bus injection plus the load, so the limits and the cost go straight onto the injection polynomial.
A bus with SEVERAL generators gets one real auxiliary variable per generator, summed into the
balance. Those auxiliaries are ordinary reals that never enter a moment matrix, so they add no
monomials and cost nothing in the hierarchy's size -- only buses that actually need them pay.

THERMAL LIMITS go in through the cone channel. |S_lm|^2 <= rate^2 is degree 4 in v, so it is out
of reach at order 1 as a polynomial inequality, but P_lm and Q_lm are each degree (1,1) and the
cone is convex, so (rate, L(P), L(Q)) in SOC is valid at EVERY order.

COSTS are linear only; a quadratic cost is quartic in v and the paper lifts it with real epigraph
variables, which breaks T-invariance and with it the order-1 = SDP structure.
"""
function build_complex_power_pop(data::Dict{String,Any}; thermal_limits::Bool = true)
    ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]
    Y, buses, idx = complex_ybus(data)
    n = length(buses)
    pop = CPOP()
    for b in buses
        bus = ref[:bus][b]
        add_cvar!(pop, "v[$b]"; absmin = Float64(bus["vmin"]), absmax = Float64(bus["vmax"]))
    end
    gens = Dict{Int,Vector{Any}}()
    for (_, g) in ref[:gen]
        Int(get(g, "gen_status", 1)) == 0 && continue
        length(g["cost"]) <= 3 ||
            error("generator at bus $(g["gen_bus"]) has a cost of degree > 2; unsupported")
        push!(get!(gens, Int(g["gen_bus"]), Any[]), g)
    end
    pd = Dict(b => 0.0 for b in buses); qd = Dict(b => 0.0 for b in buses)
    for (_, l) in ref[:load]
        pd[Int(l["load_bus"])] += Float64(l["pd"]); qd[Int(l["load_bus"])] += Float64(l["qd"])
    end
    cost1(g) = (c = g["cost"]; length(c) >= 2 ? Float64(c[end-1]) : 0.0)
    cost2(g) = (c = g["cost"]; length(c) >= 3 ? Float64(c[end-2]) : 0.0)
    isquad(g) = cost2(g) > 0

    obj = CPoly()
    obj_aux = Dict{Int,Float64}()
    n_aux_buses = 0
    for b in buses
        k = idx[b]
        ek = ekcol(n, k)
        YHe = adjoint(Y) * ek
        Pk = quadform((YHe + adjoint(YHe)) / 2)
        Qk = quadform((YHe - adjoint(YHe)) / (2im))
        bus = ref[:bus][b]
        add_cineq!(pop, quadform(ek) - (Float64(bus["vmin"])^2), "vmin[$b]")
        add_cineq!(pop, (Float64(bus["vmax"])^2) - quadform(ek), "vmax[$b]")
        gl = get(gens, b, Any[])
        if isempty(gl)
            add_ceq!(pop, Pk + pd[b], "p_balance[$b]")
            add_ceq!(pop, Qk + qd[b], "q_balance[$b]")
        elseif length(gl) == 1 && !isquad(gl[1])
            g = gl[1]
            add_cineq!(pop, Pk - (Float64(g["pmin"]) - pd[b]), "pg_min[$b]")
            add_cineq!(pop, (Float64(g["pmax"]) - pd[b]) - Pk, "pg_max[$b]")
            add_cineq!(pop, Qk - (Float64(g["qmin"]) - qd[b]), "qg_min[$b]")
            add_cineq!(pop, (Float64(g["qmax"]) - qd[b]) - Qk, "qg_max[$b]")
            obj = obj + cost1(g) * (Pk + pd[b])
        else
            # several generators, or one with a quadratic cost: a real auxiliary per generator,
            # summed into the balance. A quadratic cost is then quadratic in an ORDINARY REAL
            # variable and goes in as a rotated cone, so it never reaches the moment matrix.
            n_aux_buses += 1
            pidx = Int[]; qidx = Int[]
            for g in gl
                ip = add_caux!(pop, "pg[$(g["index"])]"; lo = Float64(g["pmin"]), hi = Float64(g["pmax"]))
                iq = add_caux!(pop, "qg[$(g["index"])]"; lo = Float64(g["qmin"]), hi = Float64(g["qmax"]))
                push!(pidx, ip); push!(qidx, iq)
                obj_aux[ip] = get(obj_aux, ip, 0.0) + cost1(g)
                if isquad(g)
                    it = add_caux!(pop, "tcost[$(g["index"])]"; lo = 0.0)
                    add_caux_quad!(pop, it, cost2(g), ip)
                    obj_aux[it] = get(obj_aux, it, 0.0) + 1.0
                end
            end
            # sum_g pg_g - pd_b - P_k(v) = 0, and likewise for reactive
            add_cmixed!(pop, Dict(i => 1.0 for i in pidx), (-1.0) * Pk - pd[b], :eq, "p_balance[$b]")
            add_cmixed!(pop, Dict(i => 1.0 for i in qidx), (-1.0) * Qk - qd[b], :eq, "q_balance[$b]")
        end
    end

    n_thermal = 0
    if thermal_limits
        for (l, br) in ref[:branch]
            rate = Float64(get(br, "rate_a", Inf))
            (isfinite(rate) && rate > 0) || continue
            f, t = idx[br["f_bus"]], idx[br["t_bus"]]
            gs, bs = PowerModels.calc_branch_y(br)
            trr, tii = PowerModels.calc_branch_t(br)
            y = complex(gs, bs); T = complex(trr, tii); tm2 = abs2(T)
            # S_fr = v_f conj(i_fr), i_fr = Yff v_f + Yft v_t
            Yff = (y + complex(br["g_fr"], br["b_fr"])) / tm2
            Yft = -y / conj(T)
            Ytt = y + complex(br["g_to"], br["b_to"])
            Ytf = -y / T
            for (a, c, d, tagend) in ((f, conj(Yff), conj(Yft), "fr"), (t, conj(Ytt), conj(Ytf), "to"))
                o = a == f ? t : f
                S = CPoly(); caddterm!(S, ([a], [a]), c); caddterm!(S, ([o], [a]), d)
                P = 0.5 * (S + cconj(S))
                Q = (-0.5im) * (S - cconj(S))
                add_csoc!(pop, CPoly(rate), [P, Q], "thermal[$l _$tagend]")
                n_thermal += 1
            end
        end
    end

    pop.obj = obj
    pop.meta["obj_aux"] = obj_aux
    pop.meta["ybus"] = Y
    pop.meta["buses"] = buses
    pop.meta["bus_index"] = idx
    pop.meta["n_aux_buses"] = n_aux_buses
    pop.meta["n_thermal_cones"] = n_thermal
    pop.meta["const_cost"] = sum(sum((length(g["cost"]) >= 1 ? Float64(g["cost"][end]) : 0.0)
                                     for g in gl; init = 0.0) for gl in values(gens); init = 0.0)
    return pop
end

"n x n matrix with a single 1 at (k, k), i.e. e_k e_k^T."
ekcol(n::Int, k::Int) = (M = zeros(ComplexF64, n, n); M[k, k] = 1; M)

"""
    complex_bus_cliques(data, bus_index; order2_buses = Int[]) -> Vector{Vector{Int}}

Clique decomposition for the complex hierarchy, from the BUS ADJACENCY. This is where the complex
formulation is structurally simpler than the real one: a bus is ONE complex variable, so the
interaction graph is the network graph itself, with no e/f pairing to undo.

Buses promoted to order 2 have their closed neighbourhood clique-ified first, exactly as in
Molzahn & Hiskens and as the paper's E^con prescribes -- a high-order constraint couples all the
variables it touches, so its support must sit inside a single clique. The rest of the graph is
left alone and a min-degree chordal extension supplies the maximal cliques.
"""
function complex_bus_cliques(data::Dict{String,Any}, bus_index::Dict{Int,Int};
        order2_buses = Int[])
    supports = Vector{Int}[]
    nbrs = Dict{Int,Set{Int}}()
    for (_, br) in get(data, "branch", Dict{String,Any}())
        Int(get(br, "br_status", 1)) == 0 && continue
        f, t = Int(br["f_bus"]), Int(br["t_bus"])
        f == t && continue
        (haskey(bus_index, f) && haskey(bus_index, t)) || continue
        push!(supports, sort!([bus_index[f], bus_index[t]]))
        push!(get!(nbrs, f, Set{Int}([f])), t)
        push!(get!(nbrs, t, Set{Int}([t])), f)
    end
    for b in order2_buses
        haskey(nbrs, b) || continue
        cl = sort!([bus_index[j] for j in nbrs[b] if haskey(bus_index, j)])
        length(cl) >= 2 && push!(supports, cl)
    end
    for (b, i) in bus_index
        push!(supports, [i])          # isolated buses still need a block
    end
    return chordal_cliques(supports)
end
