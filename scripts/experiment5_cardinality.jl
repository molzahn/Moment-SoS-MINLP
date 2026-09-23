# Experiment 5: how many lines should a rounded sample open at once?
#   julia --project=. scripts/experiment5_cardinality.jl <case ...>
#   env: N (samples/scheme, default 50), KMAX_SET (default "1,2,3,4,5,6,8"), POLISH_EVALS (default 0),
#        SEED (default 1), RELAX_LABEL (marginal source; default mixed_nobigM then base)
#
# Rounding collapses above 118 buses (0-2% AC-feasible) not because samples island the network but
# because they open far too many lines at once: independent rounding opens ~Sum(1 - y_l), which is 25 on
# 300-ieee and 724 of 1807 binaries on 1354-pegase. `sample_cardinality` caps the count while still
# letting the relaxation choose which lines (Gumbel top-k weighted by the opening marginal 1 - y_l).
#
# This sweeps the cap. It reuses the marginals already stored by experiments 3/4, so no SDP is re-solved
# and the whole study is cheap -- which also means it fits comfortably in the free `embers` QOS.
# Writes results/experiment5_<case>.json.
include("dcots_cases.jl")
using JSON, Printf, Random, Statistics

const N            = parse(Int, get(ENV, "N", "50"))
const SEED         = parse(Int, get(ENV, "SEED", "1"))
const POLISH_EVALS = parse(Int, get(ENV, "POLISH_EVALS", "0"))
const KMAX_SET     = [parse(Int, x) for x in split(get(ENV, "KMAX_SET", "1,2,3,4,5,6,8"), ",")]

jsonsafe(x::AbstractFloat) = isfinite(x) ? x : nothing
jsonsafe(x::AbstractDict) = Dict(k => jsonsafe(v) for (k, v) in x)
jsonsafe(x::AbstractVector) = map(jsonsafe, x)
jsonsafe(x) = x

"Marginals saved by experiment 3 or 4, preferring the no-big-M relaxation."
function load_marginals(name)
    for (f, labels) in (("experiment3_$(name).json", ("mixed_nobigM", "mixed")),
                        ("experiment4_$(name).json", ("base", "conn")))
        p = joinpath(@__DIR__, "..", "results", f)
        isfile(p) || continue
        rel = get(JSON.parsefile(p), "relaxations", Dict())
        for l in labels
            haskey(rel, l) && haskey(rel[l], "marginals") || continue
            return Float64.(rel[l]["marginals"]), "$(f):$(l)"
        end
    end
    return nothing, ""
end

"One-flip local search from `bits`, at most `budget` evaluations."
function polish(ev, bits, budget)
    best = copy(bits); r = evaluate!(ev, best)
    r.feasible || return (best, Inf, 0)
    bestcost = r.cost; used = 0
    improved = true
    while improved && used < budget
        improved = false
        for i in eachindex(best)
            used >= budget && break
            trial = copy(best); trial[i] = !trial[i]
            rr = evaluate!(ev, trial); used += 1
            if rr.feasible && rr.cost < bestcost - 1e-6
                best, bestcost, improved = trial, rr.cost, true
            end
        end
    end
    return (best, bestcost, used)
end

function run_case(name)
    μ, src = load_marginals(name)
    if μ === nothing
        @warn "no saved marginals; run experiment 3 or 4 first" case = name
        return
    end
    pop, _ = load_dcots_instance(name)
    # The marginals belong to the relaxation that produced them, so the POP must be built with the SAME
    # merge setting; otherwise the binary vectors have different lengths and meanings. (This is only a
    # bookkeeping constraint: candidate configurations are always evaluated on the original, unmerged
    # network via pop.meta["data"], so a feasible solution found here is valid regardless.)
    nb = length(binvars(pop))
    if length(μ) != nb
        error("marginal/POP mismatch for $name: $(length(μ)) saved marginals from $src but the POP has " *
              "$nb binaries (MERGE_ZMAX=$(MERGE_ZMAX)). Either set MERGE_ZMAX to the value used for that " *
              "run, or the marginals predate a model change and must be regenerated. On 1354-pegase the " *
              "latter applies: its stored marginals come from the pre-fix merged model, which was " *
              "infeasible (see results/merging_findings.md), so re-run experiment 3 or 4 first.")
    end
    guard = IslandGuard(pop)
    ev = ConfigEvaluator(pop)
    bk = best_known_cost(name)
    expected = sum(1 .- μ)
    res = Dict{String,Any}("case" => name, "N" => N, "seed" => SEED, "marginal_source" => src,
        "n_binaries" => length(μ), "expected_opened" => expected, "best_known" => bk,
        "schemes" => Dict{String,Any}())
    @printf("%s: %d binaries, E[#opened] = %.1f, best known %.2f (marginals from %s)\n",
        name, length(μ), expected, bk, src)

    bins = pop.meta["binaries"]
    opened_lines(bits) = Int[bins[i][2] for i in eachindex(bits) if !bits[i]]

    function report(tag, samples)
        t = @elapsed rs = [evaluate!(ev, s) for s in samples]
        feas = [r for r in rs if r.feasible]
        best = isempty(feas) ? nothing : minimum(r.cost for r in feas)
        bestbits = best === nothing ? nothing : samples[argmin([r.feasible ? r.cost : Inf for r in rs])]
        o = Dict{String,Any}("mean_opened" => mean(count(!, s) for s in samples),
            "unique" => length(unique(samples)), "feasible" => length(feas) / length(samples),
            "n_feasible" => length(feas), "best" => best, "eval_time" => t,
            "gap_vs_best_known" => best === nothing ? nothing : (best - bk) / bk,
            "best_opened" => bestbits === nothing ? nothing : opened_lines(bestbits),
            "within_1pct" => count(r -> r.cost <= bk * 1.01, feas) / length(samples))
        res["schemes"][tag] = o
        @printf("  %-22s opened %6.1f  unique %3d  feasible %3d/%3d (%3.0f%%)  best %-24s %5.0fs\n",
            tag, o["mean_opened"], o["unique"], length(feas), length(samples),
            100 * o["feasible"], best === nothing ? "–" : @sprintf("%.2f (%+.2f%%)", best, 100 * (best - bk) / bk), t)
        flush(stdout)
        return feas
    end

    rng = MersenneTwister(SEED)
    ind = [BitVector(rand(rng, length(μ)) .< μ) for _ in 1:N]
    report("independent", ind)
    report("independent+repair", first(repair_samples(guard, ind, μ)))
    report("cardinality_E", sample_cardinality(μ, N; rng = MersenneTwister(SEED), guard = guard))
    bestfeas = nothing; besttag = ""
    for k in sort(KMAX_SET)
        f = report("cardinality_k$k", sample_cardinality(μ, N; rng = MersenneTwister(SEED), kmax = k, guard = guard))
        if !isempty(f) && (bestfeas === nothing || minimum(r.cost for r in f) < bestfeas)
            bestfeas = minimum(r.cost for r in f); besttag = "cardinality_k$k"
        end
    end
    res["best_scheme"] = besttag
    res["best_rounded"] = bestfeas
    if besttag != ""
        res["best_opened"] = res["schemes"][besttag]["best_opened"]
        res["improves_best_known"] = bestfeas < bk - 1e-6
    end

    if POLISH_EVALS > 0 && bestfeas !== nothing
        cand = argmin(r -> r.cost, [r for r in values(ev.cache) if r.feasible])
        bits = first(k for (k, v) in ev.cache if v === cand)
        _, pc, used = polish(ev, BitVector(bits), POLISH_EVALS)
        res["polished"] = pc; res["polish_evals"] = used
        @printf("  polish: %.2f (%+.2f%%) after %d evals\n", pc, 100 * (pc - bk) / bk, used)
    end
    res["n_configs_evaluated"] = length(ev.cache)
    open(joinpath(@__DIR__, "..", "results", "experiment5_$(name).json"), "w") do io
        JSON.print(io, jsonsafe(res), 1)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    foreach(run_case, ARGS)
end
