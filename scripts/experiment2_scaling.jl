# Experiment 2: sparse (mixed-order) moment relaxations and rounding on larger instances.
#   julia --project=. scripts/experiment2_scaling.jl [instance ...]
#
# For each instance:
#   * baselines: all binaries on; NLP relax-and-round (Ipopt on the POP with binaries in [0,1],
#     then threshold and independent rounding)
#   * relaxations: order 1; mixed order (order 2 on cliques containing binaries), without and with capped
#     clique augmentation (MAXCLIQUE); large linear constraints handled on first moments (global_linear = 8)
#   * rounding from each relaxation: threshold, independent, gaussian, conditional
#   * every distinct configuration is evaluated with PowerModels AC-OPF (Ipopt)
#   * certified gap = (best found cost − best lower bound) / best found cost
#   * instances with ≤ 12 binaries are also enumerated for ground truth
include("enumerate_instances.jl")
using Random, Printf, JSON, Statistics

const N = 200
const SEED = 20260913
const ENUM_MAX = 12
const MAXCLIQUE = parse(Int, get(ENV, "MAXCLIQUE", "14"))

function nlp_relax_round(pop, ev, rng)
    t = @elapsed r = solve_nlp(pop)
    bv = binvars(pop)
    if !r.feasible
        return Dict("status" => string(r.status), "time" => t)
    end
    μ = clamp.(r.x[bv], 0.0, 1.0)
    thr = evaluate!(ev, BitVector(μ .>= 0.5))
    samples = [BitVector(rand(rng, length(μ)) .< μ) for _ in 1:N]
    costs = [evaluate!(ev, s).cost for s in samples]
    return Dict("status" => string(r.status), "time" => t, "objective" => r.objective, "marginals" => μ,
        "threshold_cost" => isfinite(thr.cost) ? thr.cost : nothing,
        "independent_best" => isfinite(minimum(costs)) ? minimum(costs) : nothing,
        "independent_feas" => mean(isfinite.(costs)), "independent_unique" => length(unique(samples)))
end


"Default variant set: order 1, mixed (order 2 on binary cliques), mixed with capped clique augmentation."
function base_variants(pop, res)
    # augmentation candidates: each binary with its bus voltages, then pairs of binaries at a common bus;
    # accepted greedily while binary-containing cliques stay within MAXCLIQUE variables
    pairs = unique([sort([s[i], s[j]]) for s in pop.meta["bus_binaries"] for i in eachindex(s) for j in i+1:length(s)
                    if pop.isbin[s[i]] && pop.isbin[s[j]]])
    t_aug = @elapsed aug = capped_augmentation(pop, [pop.meta["binary_neighborhoods"]; pairs];
        maxclique = MAXCLIQUE, global_linear = 8)
    res["augmentation"] = Dict("candidates" => length(pop.meta["binary_neighborhoods"]) + length(pairs),
        "accepted" => length(aug), "time" => t_aug, "maxclique" => MAXCLIQUE)
    @printf("  augmentation: accepted %d of %d candidate supports (%.1fs)\n", length(aug),
        length(pop.meta["binary_neighborhoods"]) + length(pairs), t_aug)
    variants = (("order1", (order = 1, global_linear = 8)),
        ("mixed", (order = binary_clique_order(pop), global_linear = 8)),
        ("mixed_aug", (order = binary_clique_order(pop), global_linear = 8, extra_supports = aug)))
    return variants
end

"Order-elevation variants: order 2 on cliques adjacent to binaries (size-capped), with/without pair moment blocks."
function adjacent_variants(pop, res)
    cap = parse(Int, get(ENV, "ADJ_MAXSIZE", "16"))
    res["adj_maxsize"] = cap
    ord = adjacent_clique_order(pop; maxsize = cap)
    return (("adj$(cap)", (order = ord, global_linear = 8)),
        ("adj$(cap)_pairs", (order = ord, global_linear = 8, extra_cliques = pop.meta["pair_cliques"])))
end

const VARIANT_SETS = Dict("base" => base_variants, "adjacent" => adjacent_variants)
const VARIANT_SET = get(ENV, "VARIANTS", "base")

function run_scaling(name)
    pop, desc = load_instance(name)
    bv = binvars(pop)
    nb = length(bv)
    ev = ConfigEvaluator(pop)
    res = Dict{String,Any}("instance" => name, "description" => desc, "n_binaries" => nb, "n_vars" => length(pop.names))
    @printf("\n=== %s: %s (%d vars, %d binaries)\n", name, desc, length(pop.names), nb)
    allon = evaluate!(ev, trues(nb))
    res["all_on"] = isfinite(allon.cost) ? allon.cost : nothing
    @printf("  all-on: %s\n", isfinite(allon.cost) ? @sprintf("%.2f", allon.cost) : allon.reason)
    if nb <= ENUM_MAX
        loaded = load_enumeration!(ev, name)
        t = @elapsed loaded || enumerate_configs!(ev)
        loaded || save_enumeration(name, pop, desc, ev, t)
        b, r = best_config(ev)
        res["enumerated_opt"] = r.cost
        res["enumerated_opt_config"] = join(Int.(b))
        @printf("  enumerated optimum: %.2f at %s (%.0fs)\n", r.cost, join(Int.(b)), t)
    end
    flush(stdout)
    nlp = nlp_relax_round(pop, ev, MersenneTwister(SEED))
    res["nlp_relax_round"] = nlp
    @printf("  NLP relaxation: %s obj %s; threshold %s; independent best %s (feas %s)\n", nlp["status"],
        get(nlp, "objective", "–"), get(nlp, "threshold_cost", "–"), get(nlp, "independent_best", "–"),
        get(nlp, "independent_feas", "–"))
    flush(stdout)

    variants = VARIANT_SETS[VARIANT_SET](pop, res)
    res["relaxations"] = Dict{String,Any}()
    for (lbl, kw) in variants
        t = @elapsed rel = solve_moment_relaxation(pop; kw...)
        μ, _, R = binary_correlation(rel, bv)
        npairs = count(!ismissing(moment(rel, sort([bv[i], bv[j]]))) for i in 1:nb for j in i+1:nb)
        o = Dict{String,Any}("bound" => rel.bound, "status" => string(rel.status), "time" => t, "solve_time" => rel.solve_time,
            "build_time" => rel.build_time, "n_cliques" => length(rel.cliques), "n_order2" => count(==(2), rel.orders),
            "max_clique" => maximum(length.(rel.cliques)), "max_psd" => maximum(rel.psd_sizes), "n_moments" => rel.n_moments,
            "skipped" => rel.skipped, "marginals" => μ, "pair_moments" => [npairs, nb * (nb - 1) ÷ 2],
            "frac_integral" => mean(min.(μ, 1 .- μ) .< 1e-3), "schemes" => Dict{String,Any}())
        @printf("  %-12s bound %.2f %s total %.1fs (solve %.1fs) cliques %d (order2 %d) max %d maxPSD %d, pair moments %d/%d, integral marginals %.0f%%\n",
            lbl, rel.bound, rel.status, t, rel.solve_time, length(rel.cliques), o["n_order2"], o["max_clique"], o["max_psd"],
            npairs, nb * (nb - 1) ÷ 2, 100o["frac_integral"])
        flush(stdout)
        samplers = (("threshold", () -> sample_threshold(rel, bv)),
            ("independent", () -> sample_independent(rel, bv, N; rng = MersenneTwister(SEED))),
            ("gaussian", () -> sample_gaussian(rel, bv, N; rng = MersenneTwister(SEED))),
            ("conditional", () -> first(sample_conditional(rel, bv, N; rng = MersenneTwister(SEED)))))
        for (sch, f) in samplers
            s = f()
            t_eval = @elapsed costs = [evaluate!(ev, b).cost for b in s]
            best = minimum(costs)
            counts = Dict{String,Int}()
            for b in s
                counts[join(Int.(b))] = get(counts, join(Int.(b)), 0) + 1
            end
            o["schemes"][sch] = Dict("best" => isfinite(best) ? best : nothing, "feas" => mean(isfinite.(costs)),
                "unique" => length(unique(s)), "eval_time" => t_eval, "config_counts" => counts,
                "costs_sorted" => sort([isfinite(c) ? c : nothing for c in costs]; by = x -> x === nothing ? Inf : x))
            @printf("    %-12s best %s  feas %.0f%%  unique %d  (eval %.1fs)\n", sch,
                isfinite(best) ? @sprintf("%.2f", best) : "infeasible", 100mean(isfinite.(costs)), length(unique(s)), t_eval)
            flush(stdout)
        end
        res["relaxations"][lbl] = o
    end

    feas = [(b, r) for (b, r) in ev.cache if r.feasible]
    best_found = isempty(feas) ? Inf : minimum(r.cost for (_, r) in feas)
    best_bound = maximum(o["bound"] for o in values(res["relaxations"]) if isfinite(o["bound"]))
    res["best_found"] = best_found
    res["best_bound"] = best_bound
    res["certified_gap"] = (best_found - best_bound) / best_found
    res["n_configs_evaluated"] = length(ev.cache)
    @printf("  best found %.2f, best bound %.2f, certified gap %.3f%%, %d configurations evaluated\n",
        best_found, best_bound, 100res["certified_gap"], length(ev.cache))
    return res
end

if abspath(PROGRAM_FILE) == @__FILE__
    names = isempty(ARGS) ? ["case14_ots", "case30_uc", "case30_ots", "case24_uc", "case24_ots"] : ARGS
    for name in names
        r = run_scaling(name)
        suffix = VARIANT_SET == "base" ? "" : "_" * VARIANT_SET
        open(joinpath(@__DIR__, "..", "results", "experiment2$(suffix)_$(name).json"), "w") do io
            JSON.print(io, r, 1)
        end
    end
end
