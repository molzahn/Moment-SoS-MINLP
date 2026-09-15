# Experiment 4: rounding that avoids infeasible islands (DC-OTS paper cases, all branches switchable).
#   julia --project=. scripts/experiment4_islanding.jl <case name, e.g. 118-ieee-api>
#
# The four steps (src/islands.jl, src/rounding.jl):
#   1. evaluation allows islands (as PowerModels' AC-OTS model does) and rejects only islands that provably
#      cannot balance active power; the strict "connected" screen of experiments 1-3 is reported for reference
#   2. repair: re-close opened lines (largest marginal first) until every island can balance power
#   3. guarded conditional sampling: open a line only if that cannot create such an island
#   4. relaxation with radial lines fixed closed and cut constraints Σ_{δ(S)} z >= 1 (|S| <= CUT_SET); cuts stay
#      out of the sparsity graph, and each cut with <= CUT_BLOCK_MAX lines gets an order-2 moment block on its binaries
#
# Relaxations: "base" = mixed order policy without big-M (experiment 3's best at >= 73 buses), with low-impedance
# merging; "conn" = the same plus step 4. Samples from both are evaluated on the same original network.
include("dcots_cases.jl")
using Random, Printf, JSON, Statistics

const SEED = 20260913
const N = parse(Int, get(ENV, "N_SAMPLES", "100"))
const MAX_SDP_TIME = parse(Float64, get(ENV, "MAX_SDP_TIME", "3600"))
const POLISH_EVALS = parse(Int, get(ENV, "POLISH_EVALS", "40"))
const CUT_SET = parse(Int, get(ENV, "CUT_SET", "2"))
const CUT_BLOCK_MAX = parse(Int, get(ENV, "CUT_BLOCK_MAX", "6"))   # cuts with <= this many lines get an order-2 block on their binaries
const RELAX = split(get(ENV, "RELAX", "base,conn"), ",")

jsonsafe(x::AbstractFloat) = isfinite(x) ? x : nothing
jsonsafe(x::AbstractDict) = Dict(k => jsonsafe(v) for (k, v) in x)
jsonsafe(x::AbstractVector) = map(jsonsafe, x)
jsonsafe(x) = x

function polish(ev, start::BitVector, μ::Vector{Float64}, budget::Int)
    best, bestc = copy(start), evaluate!(ev, start).cost
    order = sortperm(abs.(μ .- 0.5))
    evals, improved = 0, true
    while improved && evals < budget
        improved = false
        for j in order
            evals >= budget && break
            cand = copy(best)
            cand[j] = !cand[j]
            haskey(ev.cache, cand) || (evals += 1)
            c = evaluate!(ev, cand).cost
            if c < bestc - 1e-6 * abs(bestc)
                best, bestc, improved = cand, c, true
                break
            end
        end
    end
    return best, bestc, evals
end

function run_case(name)
    row = paper_row(name)
    base_pop, data = load_dcots_instance(name; bigM_switching = false)
    base_ids = [id for (_, id) in base_pop.meta["binaries"]]
    ev = ConfigEvaluator(base_pop)                       # islands allowed (step 1)
    allclosed = evaluate!(ev, trues(length(base_ids)))
    best_known = minimum(filter(!isnothing, [row[4], row[5]]))
    res = Dict{String,Any}("case" => name, "buses" => length(data["bus"]), "n_binaries_base" => length(base_ids),
        "paper_odcots" => row[4], "paper_acots" => row[5], "acopf" => allclosed.cost, "N" => N, "cut_set" => CUT_SET,
        "relaxations" => Dict{String,Any}())
    @printf("\n=== %s: %d binaries; AC-OPF %.2f; paper O-DC-OTS %d, AC-OTS %s\n", name, length(base_ids), allclosed.cost,
        row[4], something(row[5], "–"))
    flush(stdout)
    params = Dict{String,Any}("MSK_DPAR_OPTIMIZER_MAX_TIME" => MAX_SDP_TIME)

    for lbl in RELAX
        t_build = @elapsed pop = lbl == "base" ? base_pop :
            first(load_dcots_instance(name; bigM_switching = false, fix_radial = true, conn_cuts = CUT_SET))
        ids = [id for (_, id) in pop.meta["binaries"]]
        pos = Dict(id => k for (k, id) in enumerate(base_ids))
        tobase(b) = (x = trues(length(base_ids)); for (k, id) in enumerate(ids); x[pos[id]] = b[k]; end; x)
        bv = binvars(pop)
        cut_blocks = filter(c -> length(c) <= CUT_BLOCK_MAX, get(pop.meta, "conn_cut_cliques", Vector{Vector{Int}}()))
        t = @elapsed rel = solve_moment_relaxation(pop; order = binary_clique_order(pop), global_linear = 8,
            solver_params = params, out_of_graph_tags = ["conn_cut", "conn_cut_wide"], extra_cliques = cut_blocks)
        μ = marginals(rel, bv)
        μbase = ones(length(base_ids))
        for (k, id) in enumerate(ids)
            μbase[pos[id]] = μ[k]
        end
        o = Dict{String,Any}("bound" => rel.bound, "raw_bound" => get(rel.info, "raw_bound", nothing), "status" => string(rel.status),
            "time" => t, "build_time" => t_build, "n_binaries" => length(ids), "radial_fixed" => pop.meta["radial_fixed"],
            "n_conn_cuts" => pop.meta["n_conn_cuts"], "n_cut_blocks" => length(cut_blocks), "max_psd" => maximum(rel.psd_sizes), "n_cliques" => length(rel.cliques),
            "marginals_below_half" => count(<(0.5), μ), "expected_opened" => sum(1 .- μ), "schemes" => Dict{String,Any}())
        @printf("  %-5s bound %.2f (raw %.2f) %s  %.1fs  binaries %d (radial fixed %d, cuts %d)  max PSD %d  E[#opened] %.1f\n",
            lbl, rel.bound, something(o["raw_bound"], NaN), rel.status, t, length(ids), length(pop.meta["radial_fixed"]),
            pop.meta["n_conn_cuts"], o["max_psd"], o["expected_opened"])
        flush(stdout)

        guard = IslandGuard(pop)
        stats = Dict{String,Any}()
        raw = Dict(
            "independent" => sample_independent(rel, bv, N; rng = MersenneTwister(SEED)),
            "gaussian" => sample_gaussian(rel, bv, N; rng = MersenneTwister(SEED)),
            "conditional" => first(sample_conditional(rel, bv, N; rng = MersenneTwister(SEED))))
        sets = Tuple{String,Vector{BitVector},Float64}[]
        for s in ("independent", "gaussian", "conditional")
            push!(sets, (s, raw[s], 0.0))
            rep, nclosed = repair_samples(guard, raw[s], μ)
            push!(sets, (s * "+repair", rep, mean(nclosed)))
        end
        tg = @elapsed g = first(sample_conditional(rel, bv, N; rng = MersenneTwister(SEED), guard = guard, order = :open_first, stats = stats))
        push!(sets, ("conditional_guarded", g, get(stats, "forced_closed", 0) / N))

        bestb, bestc = trues(length(base_ids)), allclosed.cost
        for (s, samples, fixes) in sets
            bs = tobase.(samples)
            strict = mean(first(screen_config(config_data(data, base_pop.meta["binaries"], b); islands = :connected)) for b in bs)
            te = @elapsed rs = [evaluate!(ev, b) for b in bs]
            costs = [r.cost for r in rs]
            k = argmin(costs)
            costs[k] < bestc && ((bestb, bestc) = (bs[k], costs[k]))
            reasons = Dict{String,Int}()
            for r in rs
                r.feasible || (reasons[r.reason] = get(reasons, r.reason, 0) + 1)
            end
            fc = filter(isfinite, costs)
            o["schemes"][s] = Dict("best" => isfinite(costs[k]) ? costs[k] : nothing, "feas" => mean(isfinite.(costs)),
                "connected" => strict, "unique" => length(unique(bs)), "mean_fixes" => fixes, "eval_time" => te,
                "median_feasible" => isempty(fc) ? nothing : median(fc), "infeasible_reasons" => reasons,
                "within_1pct_best_known" => mean(costs .<= best_known * 1.01),
                "best_opened" => [base_ids[i] for i in findall(.!bs[k])])
            @printf("    %-22s best %-12s feas %3.0f%%  connected %3.0f%%  ≤1%% of best known %3.0f%%  unique %3d  fixes/sample %.1f  (%.0fs)  %s\n",
                s, isfinite(costs[k]) ? @sprintf("%.2f", costs[k]) : "infeasible", 100mean(isfinite.(costs)), 100strict,
                100o["schemes"][s]["within_1pct_best_known"], length(unique(bs)), fixes, te, string(reasons))
            flush(stdout)
        end
        o["guarded_sampling_time"] = tg
        tp = @elapsed (pb, pc, pe) = polish(ev, bestb, μbase, POLISH_EVALS)
        o["best_rounded"] = isfinite(bestc) ? bestc : nothing
        o["best_rounded_opened"] = [base_ids[i] for i in findall(.!bestb)]
        o["polished"] = isfinite(pc) ? pc : nothing
        o["polished_opened"] = [base_ids[i] for i in findall(.!pb)]
        o["polish_evals"] = pe
        @printf("    best rounded %.2f (%d lines opened); after 1-flip polish %.2f (%d evals, %.0fs)\n", bestc,
            length(o["best_rounded_opened"]), pc, pe, tp)
        flush(stdout)
        res["relaxations"][lbl] = o
        open(joinpath(@__DIR__, "..", "results", "experiment4_$(name).json"), "w") do io
            JSON.print(io, jsonsafe(res), 1)
        end
    end
    feas = [r.cost for r in values(ev.cache) if r.feasible]
    res["best_found"] = minimum(feas)
    res["n_configs_evaluated"] = length(ev.cache)
    bounds = [o["bound"] for o in values(res["relaxations"]) if isfinite(o["bound"])]
    res["best_bound"] = isempty(bounds) ? nothing : maximum(bounds)
    bk = min(res["best_found"], best_known)
    res["certified_gap_best_known"] = isempty(bounds) ? nothing : (bk - maximum(bounds)) / bk
    @printf("  SUMMARY %s: best found %.2f (paper AC-OTS %s, O-DC-OTS %d); best bound %s; gap %s; %d configs\n", name,
        res["best_found"], something(row[5], "–"), row[4], something(res["best_bound"], "–"),
        res["certified_gap_best_known"] === nothing ? "–" : @sprintf("%.2f%%", 100res["certified_gap_best_known"]), length(ev.cache))
    open(joinpath(@__DIR__, "..", "results", "experiment4_$(name).json"), "w") do io
        JSON.print(io, jsonsafe(res), 1)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    foreach(run_case, ARGS)
end
