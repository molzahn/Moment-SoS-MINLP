# Experiment 1b: effect of clique augmentation on joint binary moments and correlated rounding.
#   julia --project=. scripts/experiment1b_cliques.jl [instance ...]
include("enumerate_instances.jl")
using Random, Printf, JSON

const N = 200
const SEED = 20260913

function pair_coverage(rel, bv)
    n = length(bv)
    have = count(!ismissing(moment(rel, sort([bv[i], bv[j]]))) for i in 1:n for j in i+1:n)
    return have, n * (n - 1) ÷ 2
end

names = isempty(ARGS) ? ["case5_ots", "case9_ots", "case5_uc_sym"] : ARGS
out = Dict{String,Any}()
for name in names
    pop, desc = load_instance(name)
    ev = ConfigEvaluator(pop)
    load_enumeration!(ev, name) || enumerate_configs!(ev)
    bopt, ropt = best_config(ev)
    opt = ropt.cost
    bv = binvars(pop)
    @printf("\n=== %s (opt %.2f at %s)\n", name, opt, join(Int.(bopt)))
    out[name] = Dict{String,Any}()
    for (variant, extra) in (("base", Vector{Vector{Int}}()), ("bus_binaries", pop.meta["bus_binaries"]),
                             ("all_binaries", pop.meta["all_binaries"]))
        for ord in (1, 2)
            rel = solve_moment_relaxation(pop; order = ord, extra_supports = extra)
            have, tot = pair_coverage(rel, bv)
            μ, _, R = binary_correlation(rel, bv)
            @printf("  %-12s order %d: bound %.2f (gap %.3f%%) %s solve %.1fs  max clique %d  pair moments %d/%d\n", variant, ord,
                rel.bound, 100(opt - rel.bound) / abs(opt), rel.status, rel.solve_time, maximum(length.(rel.cliques)), have, tot)
            println("      marginals: ", join([@sprintf("%.3f", m) for m in μ], " "))
            offdiag = [R[i, j] for i in axes(R, 1) for j in axes(R, 2) if i < j && abs(R[i, j]) > 0.05]
            println("      |corr|>0.05 pairs: ", join([@sprintf("(%d,%d)=%.2f", i, j, R[i, j]) for i in axes(R, 1) for j in axes(R, 2) if i < j && abs(R[i, j]) > 0.05], " "))
            rec = Dict{String,Any}("bound" => rel.bound, "solve_time" => rel.solve_time, "status" => string(rel.status),
                "max_clique" => maximum(length.(rel.cliques)), "pair_moments" => [have, tot], "marginals" => μ,
                "corr" => [R[i, :] for i in axes(R, 1)])
            for (sch, f) in (("independent", () -> sample_independent(rel, bv, N; rng = MersenneTwister(SEED))),
                             ("gaussian", () -> sample_gaussian(rel, bv, N; rng = MersenneTwister(SEED))),
                             ("conditional", () -> first(sample_conditional(rel, bv, N; rng = MersenneTwister(SEED)))))
                s = f()
                st = summarize_samples(ev, s, opt; tol = 1e-5)
                @printf("      %-12s uniq=%3d feas=%4s opt=%4s ≤1%%=%4s P(opt in 10)=%4s\n", sch, st.unique,
                    @sprintf("%.0f%%", 100st.feas), @sprintf("%.0f%%", 100st.opt), @sprintf("%.0f%%", 100st.within1),
                    @sprintf("%.0f%%", 100st.opt_k[10]))
                rec[sch] = Dict("unique" => st.unique, "feas" => st.feas, "opt" => st.opt, "within1" => st.within1,
                    "opt_k10" => st.opt_k[10])
            end
            flush(stdout)
            out[name]["$(variant)_order$(ord)"] = rec
        end
    end
end
open(joinpath(@__DIR__, "..", "results", "experiment1b_cliques.json"), "w") do io
    JSON.print(io, out, 1)
end
