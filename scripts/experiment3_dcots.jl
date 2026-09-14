# Experiment 3: moment relaxation + rounding on the AC-OTS test cases of Taheri & Molzahn (DC-OTS paper).
#   julia --project=. scripts/experiment3_dcots.jl <case name, e.g. 24-ieee-rts-api>
#
# All branches switchable. For each case:
#   * relaxations (certified bounds): mixed (order 2 on cliques with a binary); adj16 (order 2 on
#     size-capped cliques adjacent to binaries; ≤ 118 buses); adj16_pairs (+ pair moment blocks; ≤ 60 buses)
#   * rounding from each relaxation: threshold, independent, Gaussian, conditional (N samples)
#   * "polish": 1-flip local search around the best rounded configuration, flipping the lines whose
#     marginals are least certain first (reported separately from pure rounding)
#   * every configuration is screened (connected network) and solved with PowerModels AC-OPF + Ipopt
include("dcots_cases.jl")
using Random, Printf, JSON, Statistics

const SEED = 20260913
const MAX_SDP_TIME = parse(Float64, get(ENV, "MAX_SDP_TIME", "3600"))
const POLISH_EVALS = parse(Int, get(ENV, "POLISH_EVALS", "40"))
# VARIANTS (comma-separated labels) restricts the relaxation variants, e.g. VARIANTS=mixed

"Replace non-finite floats (e.g. a NaN bound from a failed SDP) by `nothing` so the results are valid JSON."
jsonsafe(x::AbstractFloat) = isfinite(x) ? x : nothing
jsonsafe(x::AbstractDict) = Dict(k => jsonsafe(v) for (k, v) in x)
jsonsafe(x::AbstractVector) = map(jsonsafe, x)
jsonsafe(x) = x

function polish(ev, start::BitVector, μ::Vector{Float64}, budget::Int)
    best, bestc = copy(start), evaluate!(ev, start).cost
    order = sortperm(abs.(μ .- 0.5))           # least certain first
    evals = 0
    improved = true
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
    t_build = @elapsed (pop, data) = load_dcots_instance(name)
    nbus, nbr = length(data["bus"]), length(data["branch"])
    bv = binvars(pop)
    N = nbus <= 60 ? 200 : 100
    ev = ConfigEvaluator(pop)
    base = evaluate!(ev, trues(length(bv)))
    res = Dict{String,Any}("case" => name, "buses" => nbus, "branches" => nbr, "n_vars" => length(pop.names),
        "paper_acopf" => row[3], "paper_odcots" => row[4], "paper_acots" => row[5], "paper_acots_time" => row[6],
        "acopf" => base.cost, "pop_build_time" => t_build, "relaxations" => Dict{String,Any}(),
        "merge_zmax" => MERGE_ZMAX, "merged_ties" => pop.meta["merged_ties"], "n_binaries" => length(bv))
    @printf("\n=== %s: %d buses, %d branches (%d merged ties, %d switchable), %d POP variables (build %.1fs); AC-OPF %.2f (paper %d); paper O-DC-OTS %s, AC-OTS %s\n",
        name, nbus, nbr, length(pop.meta["merged_ties"]), length(bv), length(pop.names), t_build, base.cost, row[3], row[4], something(row[5], "–"))
    flush(stdout)

    params = Dict{String,Any}("MSK_DPAR_OPTIMIZER_MAX_TIME" => MAX_SDP_TIME)
    variants = Any[("mixed", (order = binary_clique_order(pop), global_linear = 8, solver_params = params))]
    # same POP without the degree-2 big-M switching constraints (same variables, so marginals line up);
    # used where big-M makes the SDP numerically fail (89-pegase)
    push!(variants, ("mixed_nobigM", (order = binary_clique_order(pop), global_linear = 8, solver_params = params)))
    if nbus <= 118
        ord = adjacent_clique_order(pop; maxsize = 16)
        push!(variants, ("adj16", (order = ord, global_linear = 8, solver_params = params)))
        nbus <= 60 && push!(variants, ("adj16_pairs", (order = ord, global_linear = 8, solver_params = params,
            extra_cliques = pop.meta["pair_cliques"])))
    end
    keep = split(get(ENV, "VARIANTS", "mixed,adj16,adj16_pairs"), ",")   # mixed_nobigM only on request
    filter!(v -> v[1] in keep, variants)
    pop_nobigM = nothing
    for (lbl, kw) in variants
        t = @elapsed rel = try
            if lbl == "mixed_nobigM"
                pop_nobigM = first(load_dcots_instance(name; bigM_switching = false))
                @assert pop_nobigM.names == pop.names
                solve_moment_relaxation(pop_nobigM; kw...)
            else
                solve_moment_relaxation(pop; kw...)
            end
        catch err
            @printf("  %-12s FAILED: %s\n", lbl, sprint(showerror, err)[1:min(200, end)])
            nothing
        end
        rel === nothing && (res["relaxations"][lbl] = Dict("failed" => true, "time" => t); continue)
        μ = marginals(rel, bv)
        o = Dict{String,Any}("bound" => rel.bound, "raw_bound" => get(rel.info, "raw_bound", nothing),
            "status" => string(rel.status), "time" => t, "solve_time" => rel.solve_time,
            "n_cliques" => length(rel.cliques), "n_order2" => count(==(2), rel.orders), "max_psd" => maximum(rel.psd_sizes),
            "n_moments" => rel.n_moments, "marginals" => μ, "frac_integral" => mean(min.(μ, 1 .- μ) .< 1e-3),
            "n_marg_below_half" => count(<(0.5), μ), "schemes" => Dict{String,Any}())
        @printf("  %-12s bound %.2f (raw %.2f) %s  total %.1fs (solve %.1fs)  order-2 cliques %d/%d  max PSD %d  integral marginals %.0f%%, marginals<0.5: %d\n",
            lbl, rel.bound, get(rel.info, "raw_bound", NaN), rel.status, t, rel.solve_time, o["n_order2"], o["n_cliques"],
            o["max_psd"], 100o["frac_integral"], o["n_marg_below_half"])
        flush(stdout)
        bestb, bestc = trues(length(bv)), base.cost
        for (sch, f) in (("threshold", () -> sample_threshold(rel, bv)),
                         ("independent", () -> sample_independent(rel, bv, N; rng = MersenneTwister(SEED))),
                         ("gaussian", () -> sample_gaussian(rel, bv, N; rng = MersenneTwister(SEED))),
                         ("conditional", () -> first(sample_conditional(rel, bv, N; rng = MersenneTwister(SEED)))))
            s = f()
            te = @elapsed costs = [evaluate!(ev, b).cost for b in s]
            k = argmin(costs)
            if costs[k] < bestc
                bestb, bestc = s[k], costs[k]
            end
            o["schemes"][sch] = Dict("best" => isfinite(costs[k]) ? costs[k] : nothing, "feas" => mean(isfinite.(costs)),
                "unique" => length(unique(s)), "eval_time" => te,
                "best_opened" => [pop.meta["binaries"][i][2] for i in findall(.!s[k])])
            @printf("    %-12s best %s  feas %.0f%%  unique %d  (eval %.1fs)\n", sch,
                isfinite(costs[k]) ? @sprintf("%.2f", costs[k]) : "infeasible", 100mean(isfinite.(costs)), length(unique(s)), te)
            flush(stdout)
        end
        tp = @elapsed (pb, pc, pe) = polish(ev, bestb, μ, POLISH_EVALS)
        o["best_rounded"] = isfinite(bestc) ? bestc : nothing
        o["best_rounded_opened"] = [pop.meta["binaries"][i][2] for i in findall(.!bestb)]
        o["polished"] = isfinite(pc) ? pc : nothing
        o["polished_opened"] = [pop.meta["binaries"][i][2] for i in findall(.!pb)]
        o["polish_evals"] = pe
        o["polish_time"] = tp
        @printf("    best rounded %.2f (opens %s); after 1-flip polish %.2f (opens %s; %d evals, %.1fs)\n", bestc,
            string(o["best_rounded_opened"]), pc, string(o["polished_opened"]), pe, tp)
        flush(stdout)
        res["relaxations"][lbl] = o
        open(joinpath(@__DIR__, "..", "results", "experiment3_$(name).json"), "w") do io
            JSON.print(io, jsonsafe(res), 1)
        end
    end
    feas = [r.cost for r in values(ev.cache) if r.feasible]
    res["best_found"] = minimum(feas)
    res["n_configs_evaluated"] = length(ev.cache)
    bounds = [o["bound"] for o in values(res["relaxations"]) if !get(o, "failed", false) && isfinite(o["bound"])]
    res["best_bound_variant"] = isempty(bounds) ? nothing :
        first(k for (k, o) in res["relaxations"] if !get(o, "failed", false) && o["bound"] === maximum(bounds))
    res["best_bound"] = isempty(bounds) ? nothing : maximum(bounds)
    ref = minimum(filter(!isnothing, [res["best_found"], row[5], row[4]]))   # ours, paper AC-OTS, paper O-DC-OTS (AC cost)
    res["certified_gap_best_known"] = res["best_bound"] === nothing ? nothing : (ref - res["best_bound"]) / ref
    @printf("  SUMMARY %s: best found %.2f vs paper AC-OTS %s / O-DC-OTS %d; best bound %.2f; certified gap of best known %.2f%%; %d configs\n",
        name, res["best_found"], something(row[5], "–"), row[4], something(res["best_bound"], NaN),
        100something(res["certified_gap_best_known"], NaN), length(ev.cache))
    open(joinpath(@__DIR__, "..", "results", "experiment3_$(name).json"), "w") do io
        JSON.print(io, jsonsafe(res), 1)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    foreach(run_case, ARGS)
end
