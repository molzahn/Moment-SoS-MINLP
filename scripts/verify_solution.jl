# Independently re-check a rounded solution: reproduce the sample, then re-solve the AC-OPF for that
# switching configuration on the ORIGINAL (unmerged) network with PowerModels, from scratch.
#   julia --project=. scripts/verify_solution.jl <case> <N> <kmax ...>
include("dcots_cases.jl")
using JSON, Printf, Random, PowerModels
PowerModels.silence()

function main(name, N, kmaxes)
    mu = Float64.(JSON.parsefile(joinpath(@__DIR__, "..", "results", "experiment3_$(name).json"))["relaxations"]["mixed_nobigM"]["marginals"])
    pop, _ = load_dcots_instance(name)
    guard = IslandGuard(pop); ev = ConfigEvaluator(pop); bk = best_known_cost(name)
    best = Inf; bestbits = nothing; bestk = 0
    for k in kmaxes
        for s in sample_cardinality(mu, N; rng = MersenneTwister(1), kmax = k, guard = guard)
            r = evaluate!(ev, s)
            if r.feasible && r.cost < best
                best = r.cost; bestbits = copy(s); bestk = k
            end
        end
    end
    bestbits === nothing && (println("no feasible sample"); return)
    bins = pop.meta["binaries"]
    opened = Int[bins[i][2] for i in eachindex(bestbits) if !bestbits[i]]
    @printf("%s: best rounded %.2f (kmax=%d), opens %d line(s): %s\n", name, best, bestk, length(opened), string(opened))
    @printf("  previous best known %.2f -> improvement %.2f (%.3f%%)\n", bk, bk - best, 100 * (bk - best) / bk)
    d = deepcopy(pop.meta["data"])
    for l in opened
        d["branch"][string(l)]["br_status"] = 0
    end
    sol = PowerModels.solve_ac_opf(d, MomentMINLP.ipopt_optimizer())
    @printf("  INDEPENDENT PowerModels AC-OPF on the original network: %s  obj=%.2f  (match: %s)\n",
        string(sol["termination_status"]), sol["objective"], abs(sol["objective"] - best) < 1.0 ? "YES" : "NO")
end

main(ARGS[1], parse(Int, ARGS[2]), [parse(Int, x) for x in ARGS[3:end]])
