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
* `merge_zmax`: buses joined by non-transformer branches ("ties") with |r + jx| < `merge_zmax` share one
  voltage (see `low_impedance_groups`). Each tie becomes a lossless flow (p, q) with its thermal limit and
  is never switchable; every bus keeps its own power balance. This is the zero-impedance limit of the tie
  (as for PowerModels switches) and removes admittances of order 1e3–1e4 that break the SDP numerics.
  Unlike `merge_low_impedance` on the data, tie thermal limits are kept (they bind on 89-pegase).
* `merge_safe` (default `true`): keep out of the merge any tie that would make the merged model infeasible
  or shift its AC-OPF value by more than `merge_tol` (see `safe_merge_exclusions`). Without this,
  1354-pegase merges into an infeasible model and its relaxation returns a "bound" above the true optimum.
  `merge_exclude` names further ties to leave unmerged. The ties actually dropped from the merge are
  recorded in `pop.meta["merge_excluded"]`.
* `fix_radial`: keep closed (not switchable) every switchable bridge whose removal leaves an island that
  provably cannot balance active power (`radial_fixings`).
* `angle_sign`: keep `e_refbus >= 0`. It only selects V over -V, and it is the only constraint
  here with an odd power of the voltage, so it must be off for parity block splitting.
* `conn_cuts`: add cuts Σ_{l∈δ(S)} z_l >= 1 for bus sets |S| <= `conn_cuts` that cannot balance power on
  their own (`connectivity_cuts`; 0 = none). Cuts with at most `cut_graph_max` lines are tagged
  "conn_cut" (they may enter the sparsity graph), wider ones "conn_cut_wide". Putting cuts into the graph
  merges cliques badly (118-ieee: max PSD block 54 -> 1029 at 4 lines), so the default is 0 and
  `pop.meta["conn_cut_cliques"]` lists, per cut, its binary variables: pass the small ones as
  `extra_cliques` (order-2 blocks on just those binaries) and the tags in `out_of_graph_tags`.
Both are valid for every configuration that passes `island_screen`, so bounds stay valid.
"""
function build_power_pop(data::Dict{String,Any}; switchable = Int[], commitable = Int[],
    exact_switching::Bool = true, bigM_switching::Bool = true, capacity_cut::Bool = true, name::String = "",
    merge_zmax::Float64 = 0.0, merge_safe::Bool = true, merge_tol::Float64 = 0.01, merge_exclude = Int[],
    fix_radial::Bool = false, conn_cuts::Int = 0, cut_graph_max::Int = 0,
    angle_sign::Bool = false)
    ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]
    pop = POP()
    binaries = Tuple{Symbol,Int}[]
    buses = sort(collect(keys(ref[:bus])))
    refbus = first(sort(collect(keys(ref[:ref_buses]))))
    excluded = sort(collect(merge_exclude))
    if merge_zmax > 0 && merge_safe
        excluded = sort(union(excluded, cached_merge_exclusions(data; zmax = merge_zmax, tol = merge_tol, name = name)))
    end
    groups, ties = low_impedance_groups(data; zmax = merge_zmax, exclude = excluded)
    rep(i) = get(groups, i, i)
    switchable = setdiff(switchable, ties)
    radial = fix_radial ? radial_fixings(data, switchable) : Int[]
    switchable = setdiff(switchable, radial)

    E = Dict{Int,Poly}()
    F = Dict{Int,Poly}()
    for i in buses
        rep(i) == i || continue
        mem = [j for j in buses if rep(j) == i]
        vmax = minimum(ref[:bus][j]["vmax"] for j in mem)
        vmin = maximum(ref[:bus][j]["vmin"] for j in mem)
        E[i] = pvar(add_var!(pop, "e[$i]"; lb = i == refbus ? 0.0 : -vmax, ub = vmax, start = 1.0))
        F[i] = i == refbus ? Poly() : pvar(add_var!(pop, "f[$i]"; lb = -vmax, ub = vmax, start = 0.0))
        Vi = E[i] * E[i] + F[i] * F[i]
        add_ineq!(pop, Vi - vmin^2, "vmin[$i]")
        add_ineq!(pop, max(vmax, vmin)^2 - Vi, "vmax[$i]")
    end
    for i in buses
        E[i], F[i] = E[rep(i)], F[rep(i)]
    end
    Vsq(i) = E[i] * E[i] + F[i] * F[i]
    # e_refbus >= 0 only picks V over -V (reference angle 0 rather than 180 degrees); the angle
    # reference itself is already imposed by eliminating f_refbus above. It is the one constraint
    # in this model with an odd power of the voltage, so it has to go before the moment matrices
    # can be split by parity -- dropping it cannot change the optimal value.
    angle_sign && add_ineq!(pop, E[refbus], "eref[$refbus]")

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
    cap = sum((max(abs(g["pmax"]), abs(g["qmax"]), abs(g["qmin"])) for g in values(ref[:gen])); init = 0.0) +
          sum((abs(ld["pd"]) + abs(ld["qd"]) for ld in values(ref[:load])); init = 0.0)
    for l in sort(collect(keys(ref[:branch])))
        br = ref[:branch][l]
        f, t = br["f_bus"], br["t_bus"]
        if l in ties                                  # zero-impedance tie: lossless flow variable
            rate = get(br, "rate_a", Inf)
            bnd = isfinite(rate) ? min(rate, cap) : cap
            pt = pvar(add_var!(pop, "p_tie[$l]"; lb = -bnd, ub = bnd))
            qt = pvar(add_var!(pop, "q_tie[$l]"; lb = -bnd, ub = bnd))
            Pout[f] += pt + br["g_fr"] * Vsq(f); Qout[f] += qt - br["b_fr"] * Vsq(f)
            Pout[t] += br["g_to"] * Vsq(t) - pt; Qout[t] += -br["b_to"] * Vsq(t) - qt
            isfinite(rate) && add_ineq!(pop, rate^2 - pt * pt - qt * qt, "tie_thermal[$l]")
            continue
        end
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
            append!(vs, [k for k in (findfirst(==("e[$(rep(i))]"), pop.names), findfirst(==("f[$(rep(i))]"), pop.names)) if k !== nothing])
            push!(bus_bins, sort(vs))
        end
    end
    pop.meta["bus_binaries"] = bus_bins
    # "binary_neighborhoods": each binary with the voltage variables of its bus(es), so that cliques
    # containing binaries also contain the network variables they interact with
    vidx(i) = [k for k in (findfirst(==("e[$(rep(i))]"), pop.names), findfirst(==("f[$(rep(i))]"), pop.names)) if k !== nothing]
    nbhd = Vector{Vector{Int}}()
    for (kind, id) in binaries
        if kind == :gen
            push!(nbhd, sort([findfirst(==("u[$id]"), pop.names); vidx(ref[:gen][id]["gen_bus"])]))
        else
            br = ref[:branch][id]
            push!(nbhd, sort([findfirst(==("z[$id]"), pop.names); vidx(br["f_bus"]); vidx(br["t_bus"])]))
        end
    end
    pop.meta["binary_neighborhoods"] = nbhd
    # "pair_cliques": for each bus and each pair of binaries acting at that bus, a small moment block with
    # both binaries, their local flow / output variables at that bus, and the bus voltage
    idx(nm) = findfirst(==(nm), pop.names)
    pair_cl = Vector{Vector{Int}}()
    for i in buses
        local_bins = Tuple{Int,Vector{Int}}[]
        for l in sort(collect(keys(ref[:branch])))
            l in switchable || continue
            br = ref[:branch][l]
            side = br["f_bus"] == i ? "fr" : br["t_bus"] == i ? "to" : nothing
            side === nothing && continue
            push!(local_bins, (idx("z[$l]"), [idx("p_$(side)[$l]"), idx("q_$(side)[$l]")]))
        end
        for g in ref[:bus_gens][i]
            g in commitable && push!(local_bins, (idx("u[$g]"), [idx("pg[$g]"), idx("qg[$g]")]))
        end
        for a in eachindex(local_bins), b in a+1:length(local_bins)
            (ba, va), (bb, vb) = local_bins[a], local_bins[b]
            push!(pair_cl, sort([ba; bb; va; vb; vidx(i)]))
        end
    end
    pop.meta["pair_cliques"] = pair_cl
    pop.meta["all_binaries"] = [sort([findfirst(==(k == :gen ? "u[$id]" : "z[$id]"), pop.names) for (k, id) in binaries])]

    pop.meta["name"] = name
    pop.meta["data"] = data
    pop.meta["binaries"] = binaries
    pop.meta["binary_vars"] = [findfirst(==(k == :gen ? "u[$id]" : "z[$id]"), pop.names) for (k, id) in binaries]
    pop.meta["refbus"] = refbus
    # the voltage variables, i.e. the set the problem is even in -- used for parity splitting
    pop.meta["voltage_vars"] = [i for (i, n) in enumerate(pop.names) if occursin(r"^[ef]\[", n)]
    pop.meta["merged_ties"] = ties
    pop.meta["merge_excluded"] = excluded
    pop.meta["radial_fixed"] = radial
    ncuts = 0
    cut_cliques = Vector{Vector{Int}}()
    if conn_cuts > 0
        zidx = Dict(id => findfirst(==("z[$id]"), pop.names) for (k, id) in binaries if k == :branch)
        for cut in connectivity_cuts(data, switchable; maxset = conn_cuts)
            all(l -> haskey(zidx, l), cut) || continue
            g = sum(pvar(zidx[l]) for l in cut) - 1.0
            add_ineq!(pop, g, length(cut) <= cut_graph_max ? "conn_cut[$(ncuts + 1)]" : "conn_cut_wide[$(ncuts + 1)]")
            push!(cut_cliques, sort([zidx[l] for l in cut]))
            ncuts += 1
        end
    end
    pop.meta["n_conn_cuts"] = ncuts
    pop.meta["conn_cut_cliques"] = cut_cliques
    pop.meta["bus_groups"] = groups
    return pop
end

# ---------------------------------------------------------------------------------------------
# Network preprocessing

"""
    low_impedance_groups(data; zmax, exclude = Int[]) -> (groups = Dict(bus id => representative bus id),
                                                          ties = branch ids)

Ties: in-service, non-transformer branches (tap = 1, shift = 0) with |r + jx| < `zmax`, except the branch
ids in `exclude`, which stay ordinary (switchable, loss-carrying) branches. Buses connected by ties form a
group; its representative is the reference bus if the group has one, else the smallest id.
"""
function low_impedance_groups(data::Dict{String,Any}; zmax::Float64, exclude = Int[])
    skip = Set(exclude)
    istie(b) = get(b, "br_status", 1) != 0 && abs(complex(b["br_r"], b["br_x"])) < zmax && b["tap"] == 1 &&
               b["shift"] == 0 && !(b["index"] in skip)
    ties = sort([b["index"] for b in values(data["branch"]) if istie(b)])
    groups = Dict(b["index"] => b["index"] for b in values(data["bus"]))
    isempty(ties) && return (groups, ties)
    parent = copy(groups)
    findr(x) = (while parent[x] != x; parent[x] = parent[parent[x]]; x = parent[x]; end; x)
    for l in ties
        br = data["branch"][string(l)]
        ra, rb = findr(br["f_bus"]), findr(br["t_bus"])
        ra != rb && (parent[max(ra, rb)] = min(ra, rb))
    end
    members = Dict{Int,Vector{Int}}()
    for i in keys(parent)
        push!(get!(members, findr(i), Int[]), i)
    end
    for (_, mem) in members
        refs = [i for i in mem if data["bus"][string(i)]["bus_type"] == 3]
        r = isempty(refs) ? minimum(mem) : first(refs)
        for i in mem
            groups[i] = r
        end
    end
    return (groups, ties)
end

"""
    safe_merge_exclusions(data; zmax, tol = 0.01, maxiter = 25, name = "", verbose = false)
        -> (exclude = branch ids, ok, reference, merged, iters)

Branch ids that must be kept *out* of the low-impedance merge for the merged model to remain a faithful
model of `data`.

Merging a tie forces its end buses to share a voltage and remodels its flow as lossless. That is harmless
when the tie has slack, but a tie sitting exactly *at* its thermal limit has none. On 1354-pegase (184
merged ties) the all-closed merged model is infeasible while the unmerged one solves to the paper's AC-OPF
cost, and scaling only the merged-tie ratings by 1.10 restores feasibility. A relaxation built on such a
model relaxes the wrong problem, and its "bound" can exceed the true optimum: 1354-pegase reported
2.58e6–3.77e6 against a known feasible 1.498e6.

Rather than guessing a smaller `zmax`, the rule validates itself: compare the all-closed AC-OPF of the
merged model with the unmerged one, and while it is infeasible or off by more than `tol` (relative), drop
the most heavily loaded merged tie (|S| / rate_a at the unmerged solution) from the merge and retry. It
terminates in the worst case with nothing merged. `ok = false` means no validated merge was found within
`maxiter`; the caller should treat bounds from the resulting model as unverified.
"""
function safe_merge_exclusions(data::Dict{String,Any}; zmax::Float64, tol::Float64 = 0.01,
    maxiter::Int = 25, name::String = "", verbose::Bool = false)
    exclude = Int[]
    zmax > 0 || return (exclude = exclude, ok = true, reference = NaN, merged = NaN, iters = 0)
    _, ties0 = low_impedance_groups(data; zmax = zmax)
    isempty(ties0) && return (exclude = exclude, ok = true, reference = NaN, merged = NaN, iters = 0)

    sol = PowerModels.solve_ac_opf(data, ipopt_optimizer())
    reference = get(sol, "objective", NaN)
    if !(string(sol["termination_status"]) in ("LOCALLY_SOLVED", "OPTIMAL")) || !isfinite(reference)
        @warn "safe_merge_exclusions: unmerged AC-OPF did not solve; merge left unvalidated" name
        return (exclude = exclude, ok = false, reference = reference, merged = NaN, iters = 0)
    end
    loading = Dict{Int,Float64}()                    # |S| / rate_a at the unmerged solution
    for l in ties0
        rate = get(data["branch"][string(l)], "rate_a", 0.0)
        s = get(sol["solution"]["branch"], string(l), nothing)
        loading[l] = (s === nothing || !(rate > 0)) ? 0.0 : hypot(s["pf"], s["qf"]) / rate
    end

    merged = NaN
    for it in 1:maxiter
        _, ties = low_impedance_groups(data; zmax = zmax, exclude = exclude)
        isempty(ties) && return (exclude = exclude, ok = true, reference = reference, merged = reference, iters = it)
        pop = build_power_pop(data; name = name, merge_zmax = zmax, merge_exclude = exclude, merge_safe = false)
        r = solve_nlp(pop)
        merged = r.objective
        if r.feasible && isfinite(merged) && abs(merged - reference) <= tol * max(abs(reference), 1.0)
            return (exclude = exclude, ok = true, reference = reference, merged = merged, iters = it)
        end
        worst = argmax(l -> get(loading, l, 0.0), ties)
        push!(exclude, worst)
        verbose && @info "safe_merge_exclusions: un-merging tie $worst" loading = get(loading, worst, 0.0) status = string(r.status)
    end
    @warn "safe_merge_exclusions: no validated merge after $maxiter un-merges" name nexcluded = length(exclude)
    return (exclude = sort(exclude), ok = false, reference = reference, merged = merged, iters = maxiter)
end

"Cache for `safe_merge_exclusions`, which costs one AC-OPF per un-merge and is hit once per POP build."
const MERGE_EXCLUSION_CACHE = Dict{Any,Vector{Int}}()

function cached_merge_exclusions(data::Dict{String,Any}; zmax::Float64, tol::Float64, name::String = "")
    key = (zmax, tol, length(data["bus"]), length(data["branch"]),
        hash(sort([(b["index"], Float64(get(b, "rate_a", 0.0)), Float64(b["br_r"]), Float64(b["br_x"]),
            Int(get(b, "br_status", 1))) for b in values(data["branch"])])),
        hash(sort([(l["index"], Float64(l["pd"]), Float64(l["qd"])) for l in values(get(data, "load", Dict{String,Any}()))])))
    get!(MERGE_EXCLUSION_CACHE, key) do
        safe_merge_exclusions(data; zmax = zmax, tol = tol, name = name).exclude
    end
end

"""
    merge_low_impedance(data; zmax = 1e-3) -> (data = merged, groups = Dict(bus id => representative id),
                                                ties = merged branch ids, dropped = all removed branch ids)

Collapse buses joined by in-service branches with |r + jx| < `zmax` (p.u.) into one bus, as in the
LCOTS project (`merge_zero_impedance` in parameter_optimized_LCOTS/LCOPF.jl). Such branches have
admittances of order 1e3–1e4 (89-pegase: 19 branches with |z| ≈ 2.2e-4), which wreck the conditioning
of the moment/SOS relaxations (MOSEK stalls within a few iterations).

* Transformers (tap ≠ 1 or shift ≠ 0) are never merged, since their endpoints do not share a voltage.
* Representative bus: the reference bus if the group has one, else the smallest id. Bus type is the most
  specific in the group (ref > PV > PQ); the voltage boxes are intersected.
* Loads, generators, shunts and storage are moved to the representative. The merged ties are removed,
  together with any other branch that becomes a self-loop; the charging of removed branches is kept as a
  bus shunt.
* Surviving branches, generators and loads keep their original ids, so switching decisions still refer
  to the original network.

This is a modelling approximation (merged buses share one voltage, and removed branches can no longer be
switched), so a bound from the merged model is not a rigorous bound for the original network. It also drops
the thermal limits of the ties, which bind on 89-pegase (AC-OPF 125456 merged vs 130175 original); prefer
`build_power_pop(...; merge_zmax)`, which merges voltages but keeps tie flows and limits.
"""
function merge_low_impedance(data::Dict{String,Any}; zmax::Float64 = 1e-3)
    d = deepcopy(data)
    isactive(b) = get(b, "br_status", 1) != 0
    groups, ties = low_impedance_groups(d; zmax = zmax)
    isempty(ties) && return (data = d, groups = groups, ties = ties, dropped = Int[])
    members = Dict{Int,Vector{Int}}()
    for (i, r) in groups
        push!(get!(members, r, Int[]), i)
    end
    for (rep, mem) in members
        length(mem) == 1 && continue
        rb = d["bus"][string(rep)]
        types = [d["bus"][string(i)]["bus_type"] for i in mem]
        rb["bus_type"] = 3 in types ? 3 : (2 in types ? 2 : minimum(types))
        rb["vmin"] = maximum(d["bus"][string(i)]["vmin"] for i in mem)
        rb["vmax"] = max(rb["vmin"], minimum(d["bus"][string(i)]["vmax"] for i in mem))
        for i in mem
            i == rep || delete!(d["bus"], string(i))
        end
    end
    for (comp, key) in (("gen", "gen_bus"), ("load", "load_bus"), ("shunt", "shunt_bus"), ("storage", "storage_bus"))
        for c in values(get(d, comp, Dict{String,Any}()))
            c[key] = groups[c[key]]
        end
    end
    for c in values(get(d, "dcline", Dict{String,Any}()))
        c["f_bus"], c["t_bus"] = groups[c["f_bus"]], groups[c["t_bus"]]
    end
    dropped = Int[]
    nshunt = maximum(parse.(Int, collect(keys(d["shunt"]))); init = 0)
    for (k, br) in collect(d["branch"])
        f, t = groups[br["f_bus"]], groups[br["t_bus"]]
        br["f_bus"], br["t_bus"] = f, t
        f == t || continue
        push!(dropped, br["index"])
        delete!(d["branch"], k)
        isactive(br) || continue
        gs = br["g_fr"] / br["tap"]^2 + br["g_to"]
        bs = br["b_fr"] / br["tap"]^2 + br["b_to"]
        if gs != 0 || bs != 0
            nshunt += 1
            d["shunt"][string(nshunt)] = Dict{String,Any}("index" => nshunt, "shunt_bus" => f, "gs" => gs, "bs" => bs,
                "status" => 1, "source_id" => Any["merged_branch", br["index"]])
        end
    end
    return (data = d, groups = groups, ties = ties, dropped = sort(dropped))
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

"""
Cheap necessary conditions for a configuration.

* `islands = :allow` (default): islands are allowed, as in PowerModels' AC-OTS model; an island is rejected
  only if it provably cannot balance active power (`island_screen`).
* `islands = :connected`: the network must be connected (the rule used in experiments 1-3).
"""
function screen_config(d::Dict{String,Any}; islands::Symbol = :allow)
    if islands == :connected
        comps = PowerModels.calc_connected_components(d)
        length(comps) == 1 || return (false, "islanded")
    else
        ok, why, _ = island_screen(d)
        ok || return (false, why)
    end
    pmax = sum((g["pmax"] for g in values(d["gen"]) if g["gen_status"] != 0); init = 0.0)
    pd = sum((l["pd"] for l in values(d["load"]) if l["status"] != 0); init = 0.0)
    pmax >= pd || return (false, "capacity")
    return (true, "")
end

struct ConfigResult
    feasible::Bool
    cost::Float64
    reason::String      # "", "islanded", "island_deficit", "capacity", or the Ipopt termination status
    time::Float64
end

mutable struct ConfigEvaluator
    data::Dict{String,Any}
    binaries::Vector{Tuple{Symbol,Int}}
    cache::Dict{BitVector,ConfigResult}
    optimizer
    islands::Symbol      # :allow or :connected, see `screen_config`
end
ConfigEvaluator(data, binaries, cache, optimizer) = ConfigEvaluator(data, binaries, cache, optimizer, :allow)
ConfigEvaluator(pop::POP; optimizer = ipopt_optimizer(), islands::Symbol = :allow) =
    ConfigEvaluator(pop.meta["data"], pop.meta["binaries"], Dict{BitVector,ConfigResult}(), optimizer, islands)

function evaluate!(ev::ConfigEvaluator, bits::AbstractVector{Bool})
    key = BitVector(bits)
    haskey(ev.cache, key) && return ev.cache[key]
    t0 = time()
    d = config_data(ev.data, ev.binaries, key)
    ok, why = screen_config(d; islands = ev.islands)
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
