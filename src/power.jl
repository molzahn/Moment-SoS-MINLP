# AC-OPF / AC-OTS / AC-UC (single period) as POPs in rectangular voltage coordinates,
# built from PowerModels data so that NLP recovery with PowerModels' ACPPowerModel solves
# exactly the same problem for fixed binaries. See docs/formulations.md.

"""
    build_power_pop(data; switchable = Int[], commitable = Int[],
                    exact_switching = true, bigM_switching = true, capacity_cut = true)

* `switchable`: branch ids with a binary status z_l (AC-OTS).
* `commitable`: generator ids with a binary status u_g (AC-UC).
* `exact_switching`: add the degree-3 constraints p_l = z_l·P_l(V) and z_l·(angle limit) >= 0
  (enforced only in cliques of order >= 2).
* `bigM_switching`: add degree-2 on/off big-M constraints (usable at order 1).
* `capacity_cut`: valid inequality Σ_g u_g pmax_g >= Σ pd (only if losses are nonnegative).
"""
function build_power_pop(data::Dict{String,Any}; switchable = Int[], commitable = Int[],
    exact_switching::Bool = true, bigM_switching::Bool = true, capacity_cut::Bool = true, name::String = "")
    ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]
    pop = POP()
    binaries = Tuple{Symbol,Int}[]
    buses = sort(collect(keys(ref[:bus])))
    refbus = first(sort(collect(keys(ref[:ref_buses]))))

    E = Dict{Int,Poly}()
    F = Dict{Int,Poly}()
    for i in buses
        vmax = ref[:bus][i]["vmax"]
        E[i] = pvar(add_var!(pop, "e[$i]"; lb = i == refbus ? 0.0 : -vmax, ub = vmax, start = 1.0))
        F[i] = i == refbus ? Poly() : pvar(add_var!(pop, "f[$i]"; lb = -vmax, ub = vmax, start = 0.0))
    end
    Vsq(i) = E[i] * E[i] + F[i] * F[i]
    for i in buses
        b = ref[:bus][i]
        add_ineq!(pop, Vsq(i) - b["vmin"]^2, "vmin[$i]")
        add_ineq!(pop, b["vmax"]^2 - Vsq(i), "vmax[$i]")
    end
    add_ineq!(pop, E[refbus], "eref[$refbus]")

    # generators
    PG = Dict{Int,Poly}()
    QG = Dict{Int,Poly}()
    obj = Poly()
    for g in sort(collect(keys(ref[:gen])))
        gen = ref[:gen][g]
        pmin, pmax, qmin, qmax = gen["pmin"], gen["pmax"], gen["qmin"], gen["qmax"]
        on = g in commitable
        U = Poly(1.0)
        if on
            U = pvar(add_var!(pop, "u[$g]"; binary = true, start = 1.0))
            push!(binaries, (:gen, g))
        end
        PG[g] = pvar(add_var!(pop, "pg[$g]"; lb = on ? min(0.0, pmin) : pmin, ub = on ? max(0.0, pmax) : pmax,
            start = (pmin + pmax) / 2))
        QG[g] = pvar(add_var!(pop, "qg[$g]"; lb = on ? min(0.0, qmin) : qmin, ub = on ? max(0.0, qmax) : qmax,
            start = (qmin + qmax) / 2))
        add_ineq!(pop, PG[g] - pmin * U, "pg_min[$g]")
        add_ineq!(pop, pmax * U - PG[g], "pg_max[$g]")
        add_ineq!(pop, QG[g] - qmin * U, "qg_min[$g]")
        add_ineq!(pop, qmax * U - QG[g], "qg_max[$g]")
        gen["model"] == 2 || error("only polynomial generator costs are supported")
        c = gen["cost"]
        n = length(c)
        for (k, ck) in enumerate(c)
            p = n - k
            obj = obj + (p == 0 ? ck * U : ck * PG[g]^p)
        end
    end
    pop.obj = obj

    # branches
    Pout = Dict(i => Poly() for i in buses)
    Qout = Dict(i => Poly() for i in buses)
    for l in sort(collect(keys(ref[:branch])))
        br = ref[:branch][l]
        f, t = br["f_bus"], br["t_bus"]
        gs, bs = PowerModels.calc_branch_y(br)
        tr, ti = PowerModels.calc_branch_t(br)
        y = complex(gs, bs)
        T = complex(tr, ti)
        A = (y + complex(br["g_fr"], br["b_fr"])) / abs2(T)
        B = -y / conj(T)
        C = -y / T
        D = y + complex(br["g_to"], br["b_to"])
        c = E[f] * E[t] + F[f] * F[t]      # |Vf||Vt| cos(θf-θt)
        s = F[f] * E[t] - E[f] * F[t]      # |Vf||Vt| sin(θf-θt)
        Pfr = real(A) * Vsq(f) + real(B) * c + imag(B) * s
        Qfr = -imag(A) * Vsq(f) + real(B) * s - imag(B) * c
        Pto = real(D) * Vsq(t) + real(C) * c - imag(C) * s
        Qto = -imag(D) * Vsq(t) - real(C) * s - imag(C) * c
        rate = get(br, "rate_a", Inf)
        amax = tan(br["angmax"]) * c - s
        amin = s - tan(br["angmin"]) * c
        vmf, vmt = ref[:bus][f]["vmax"], ref[:bus][t]["vmax"]

        if l in switchable
            Z = pvar(add_var!(pop, "z[$l]"; binary = true, start = 1.0))
            push!(binaries, (:branch, l))
            Mfr = abs(A) * vmf^2 + abs(B) * vmf * vmt     # bound on |Pfr(V)|, |Qfr(V)|
            Mto = abs(D) * vmt^2 + abs(C) * vmf * vmt
            bfr = isfinite(rate) ? min(rate, Mfr) : Mfr
            bto = isfinite(rate) ? min(rate, Mto) : Mto
            pfr = pvar(add_var!(pop, "p_fr[$l]"; lb = -bfr, ub = bfr))
            qfr = pvar(add_var!(pop, "q_fr[$l]"; lb = -bfr, ub = bfr))
            pto = pvar(add_var!(pop, "p_to[$l]"; lb = -bto, ub = bto))
            qto = pvar(add_var!(pop, "q_to[$l]"; lb = -bto, ub = bto))
            Pout[f] += pfr; Qout[f] += qfr
            Pout[t] += pto; Qout[t] += qto
            if exact_switching
                add_eq!(pop, pfr - Z * Pfr, "sw_flow_exact[$l]")
                add_eq!(pop, qfr - Z * Qfr, "sw_flow_exact[$l]")
                add_eq!(pop, pto - Z * Pto, "sw_flow_exact[$l]")
                add_eq!(pop, qto - Z * Qto, "sw_flow_exact[$l]")
                add_ineq!(pop, Z * amax, "sw_angle_exact[$l]")
                add_ineq!(pop, Z * amin, "sw_angle_exact[$l]")
            end
            if bigM_switching || !exact_switching
                for (v, bnd) in ((pfr, bfr), (qfr, bfr), (pto, bto), (qto, bto))
                    add_ineq!(pop, bnd * Z - v, "sw_flow_onoff[$l]")
                    add_ineq!(pop, bnd * Z + v, "sw_flow_onoff[$l]")
                end
                for (v, V, M) in ((pfr, Pfr, Mfr), (qfr, Qfr, Mfr), (pto, Pto, Mto), (qto, Qto, Mto))
                    add_ineq!(pop, M * (1 - Z) - V + v, "sw_flow_bigM[$l]")
                    add_ineq!(pop, M * (1 - Z) + V - v, "sw_flow_bigM[$l]")
                end
                Ma = (max(abs(tan(br["angmax"])), abs(tan(br["angmin"]))) + 1) * vmf * vmt
                add_ineq!(pop, amax + Ma * (1 - Z), "sw_angle_bigM[$l]")
                add_ineq!(pop, amin + Ma * (1 - Z), "sw_angle_bigM[$l]")
            end
            if isfinite(rate)
                add_ineq!(pop, rate^2 - pfr * pfr - qfr * qfr, "sw_thermal[$l]")
                add_ineq!(pop, rate^2 - pto * pto - qto * qto, "sw_thermal[$l]")
            end
        else
            Pout[f] += Pfr; Qout[f] += Qfr
            Pout[t] += Pto; Qout[t] += Qto
            if isfinite(rate)
                for (P, Q) in ((Pfr, Qfr), (Pto, Qto))
                    G = [Poly(rate^2) P Q; P Poly(1.0) Poly(); Q Poly() Poly(1.0)]
                    add_pmi!(pop, G, rate^2 - P * P - Q * Q, "thermal[$l]")
                end
            end
            add_ineq!(pop, amax, "angle[$l]")
            add_ineq!(pop, amin, "angle[$l]")
        end
    end

    # power balance
    total_pd = 0.0
    for i in buses
        pd = sum((ref[:load][d]["pd"] for d in ref[:bus_loads][i]); init = 0.0)
        qd = sum((ref[:load][d]["qd"] for d in ref[:bus_loads][i]); init = 0.0)
        gsh = sum((ref[:shunt][k]["gs"] for k in ref[:bus_shunts][i]); init = 0.0)
        bsh = sum((ref[:shunt][k]["bs"] for k in ref[:bus_shunts][i]); init = 0.0)
        total_pd += pd
        pg = sum((PG[g] for g in ref[:bus_gens][i]); init = Poly())
        qg = sum((QG[g] for g in ref[:bus_gens][i]); init = Poly())
        add_eq!(pop, pg - pd - gsh * Vsq(i) - Pout[i], "p_balance[$i]")
        add_eq!(pop, qg - qd + bsh * Vsq(i) - Qout[i], "q_balance[$i]")
    end

    lossless_ok = all(br["br_r"] >= 0 && br["g_fr"] >= 0 && br["g_to"] >= 0 for br in values(ref[:branch])) &&
                  all(sh["gs"] >= 0 for sh in values(ref[:shunt]))
    if capacity_cut && !isempty(commitable) && lossless_ok
        cap = Poly()
        for g in keys(ref[:gen])
            k = findfirst(==((:gen, g)), binaries)
            pmax = ref[:gen][g]["pmax"]
            cap = cap + (k === nothing ? Poly(pmax) : pmax * pvar(findfirst(==("u[$g]"), pop.names)))
        end
        add_ineq!(pop, cap - total_pd, "capacity_cut")
    end

    # Optional clique augmentation (pass to solve_moment_relaxation via extra_supports):
    # "bus_binaries": for each bus, the binaries of incident switchable branches / local units together
    # with that bus's voltage variables, so that joint moments of neighboring binaries exist.
    bus_bins = Vector{Vector{Int}}()
    for i in buses
        vs = Int[]
        for l in keys(ref[:branch])
            br = ref[:branch][l]
            (l in switchable && (br["f_bus"] == i || br["t_bus"] == i)) &&
                push!(vs, findfirst(==("z[$l]"), pop.names))
        end
        for g in ref[:bus_gens][i]
            g in commitable && push!(vs, findfirst(==("u[$g]"), pop.names))
        end
        if length(vs) > 1
            append!(vs, [k for k in (findfirst(==("e[$i]"), pop.names), findfirst(==("f[$i]"), pop.names)) if k !== nothing])
            push!(bus_bins, sort(vs))
        end
    end
    pop.meta["bus_binaries"] = bus_bins
    pop.meta["all_binaries"] = [sort([findfirst(==(k == :gen ? "u[$id]" : "z[$id]"), pop.names) for (k, id) in binaries])]

    pop.meta["name"] = name
    pop.meta["data"] = data
    pop.meta["binaries"] = binaries
    pop.meta["binary_vars"] = [findfirst(==(k == :gen ? "u[$id]" : "z[$id]"), pop.names) for (k, id) in binaries]
    pop.meta["refbus"] = refbus
    return pop
end

# ---------------------------------------------------------------------------------------------
# Evaluation of binary configurations with PowerModels (polar AC-OPF, Ipopt)

function config_data(data::Dict{String,Any}, binaries, bits::AbstractVector{Bool})
    d = deepcopy(data)
    for (k, (kind, id)) in enumerate(binaries)
        if kind == :branch
            d["branch"][string(id)]["br_status"] = bits[k] ? 1 : 0
        else
            d["gen"][string(id)]["gen_status"] = bits[k] ? 1 : 0
        end
    end
    return d
end

"Cheap necessary conditions: network connected, enough active capacity."
function screen_config(d::Dict{String,Any})
    comps = PowerModels.calc_connected_components(d)
    length(comps) == 1 || return (false, "islanded")
    pmax = sum((g["pmax"] for g in values(d["gen"]) if g["gen_status"] != 0); init = 0.0)
    pd = sum((l["pd"] for l in values(d["load"]) if l["status"] != 0); init = 0.0)
    pmax >= pd || return (false, "capacity")
    return (true, "")
end

struct ConfigResult
    feasible::Bool
    cost::Float64
    reason::String      # "", "islanded", "capacity", or the Ipopt termination status
    time::Float64
end

mutable struct ConfigEvaluator
    data::Dict{String,Any}
    binaries::Vector{Tuple{Symbol,Int}}
    cache::Dict{BitVector,ConfigResult}
    optimizer
end
ConfigEvaluator(pop::POP; optimizer = ipopt_optimizer()) =
    ConfigEvaluator(pop.meta["data"], pop.meta["binaries"], Dict{BitVector,ConfigResult}(), optimizer)

function evaluate!(ev::ConfigEvaluator, bits::AbstractVector{Bool})
    key = BitVector(bits)
    haskey(ev.cache, key) && return ev.cache[key]
    t0 = time()
    d = config_data(ev.data, ev.binaries, key)
    ok, why = screen_config(d)
    r = if !ok
        ConfigResult(false, Inf, why, time() - t0)
    else
        res = PowerModels.solve_opf(d, PowerModels.ACPPowerModel, ev.optimizer)
        st = res["termination_status"]
        feas = st in (MOI.LOCALLY_SOLVED, MOI.ALMOST_LOCALLY_SOLVED, MOI.OPTIMAL)
        ConfigResult(feas, feas ? res["objective"] : Inf, feas ? "" : string(st), time() - t0)
    end
    ev.cache[key] = r
    return r
end

"Evaluate all 2^n configurations (small instances only)."
function enumerate_configs!(ev::ConfigEvaluator)
    n = length(ev.binaries)
    n <= 16 || error("too many binaries to enumerate ($n)")
    for k in 0:(2^n-1)
        evaluate!(ev, BitVector(digits(k; base = 2, pad = n) .== 1))
    end
    return ev
end

function best_config(ev::ConfigEvaluator)
    best = nothing
    for (b, r) in ev.cache
        if r.feasible && (best === nothing || r.cost < best[2].cost)
            best = (b, r)
        end
    end
    return best
end
