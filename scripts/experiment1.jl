# Experiment 1: order-1 vs order-2 moment relaxations × rounding schemes on enumerable instances.
#   julia --project=. scripts/experiment1.jl [instance ...]
# Requires results/enum_<instance>.json (scripts/enumerate_instances.jl); missing configurations
# are evaluated on demand.
include("enumerate_instances.jl")
using Random, Printf, JSON, Statistics, LinearAlgebra

const N_SAMPLES = 200
const N_DIVE = Dict(1 => 20, 2 => 5)
const SEED = 20260913

fmt(x) = x === nothing ? "–" : @sprintf("%.2f", x)
pct(x) = @sprintf("%.0f%%", 100x)

function run_instance(name)
    pop, desc = load_instance(name)
    ev = ConfigEvaluator(pop)
    load_enumeration!(ev, name) || enumerate_configs!(ev)
    bopt, ropt = best_config(ev)
    opt = ropt.cost
    bv = binvars(pop)
    out = Dict{String,Any}("instance" => name, "description" => desc, "opt" => opt, "opt_config" => join(Int.(bopt)),
        "binaries" => [string(k, ":", id) for (k, id) in pop.meta["binaries"]], "orders" => Dict{String,Any}())
    @printf("\n=== %s (opt %.2f at %s)\n", name, opt, join(Int.(bopt)))
    for ord in (1, 2)
        rel = solve_moment_relaxation(pop; order = ord)
        μ, Σ, R = binary_correlation(rel, bv)
        gap = (opt - rel.bound) / abs(opt)
        @printf("  order %d: bound %.2f gap %.3f%% status %s solve %.1fs, max clique %d, skipped %s\n", ord, rel.bound,
            100gap, rel.status, rel.solve_time, maximum(length.(rel.cliques)), rel.skipped)
        println("    marginals: ", join([@sprintf("%.3f", m) for m in μ], " "))
        o = Dict{String,Any}("bound" => rel.bound, "gap" => gap, "status" => string(rel.status),
            "solve_time" => rel.solve_time, "build_time" => rel.build_time, "n_moments" => rel.n_moments,
            "max_clique" => maximum(length.(rel.cliques)), "n_cliques" => length(rel.cliques),
            "skipped" => rel.skipped, "marginals" => μ, "corr" => [R[i, :] for i in axes(R, 1)],
            "schemes" => Dict{String,Any}())
        schemes = Dict{String,Any}()
        rng = MersenneTwister(SEED)
        t = @elapsed s = sample_threshold(rel, bv)
        schemes["threshold"] = (s, t, Dict())
        t = @elapsed s = sample_independent(rel, bv, N_SAMPLES; rng = rng)
        schemes["independent"] = (s, t, Dict())
        t = @elapsed s = sample_gaussian(rel, bv, N_SAMPLES; rng = rng)
        schemes["gaussian"] = (s, t, Dict())
        t = @elapsed (s, neg) = sample_conditional(rel, bv, N_SAMPLES; rng = rng)
        schemes["conditional"] = (s, t, Dict("neg_mass_per_draw" => neg))
        t = @elapsed (s, ns) = sample_dive(pop, rel, N_DIVE[ord]; rng = rng, order = ord)
        schemes["dive"] = (s, t, Dict("sdp_solves" => ns))
        for sch in ("threshold", "independent", "gaussian", "conditional", "dive")
            s, t, extra = schemes[sch]
            st = summarize_samples(ev, s, opt)
            @printf("    %-12s n=%3d uniq=%3d feas=%5s opt=%5s ≤1%%=%5s  P(opt in 10)=%5s  best gap=%s  (%.1fs) %s\n",
                sch, st.n, st.unique, pct(st.feas), pct(st.opt), pct(st.within1), pct(st.opt_k[10]),
                isfinite(st.best) ? @sprintf("%.3f%%", 100(st.best - opt) / abs(opt)) : "infeasible", t,
                isempty(extra) ? "" : string(extra))
            flush(stdout)
            counts = Dict{String,Int}()
            for b in s
                k = join(Int.(b))
                counts[k] = get(counts, k, 0) + 1
            end
            o["schemes"][sch] = merge(Dict("n" => st.n, "unique" => st.unique, "feas" => st.feas, "opt" => st.opt,
                    "within1" => st.within1, "best" => isfinite(st.best) ? st.best : nothing,
                    "opt_k" => st.opt_k, "feas_k" => st.feas_k, "time" => t,
                    "infeasible_reasons" => st.infeasible_reasons, "config_counts" => counts), extra)
        end
        out["orders"][string(ord)] = o
    end
    return out
end

if abspath(PROGRAM_FILE) == @__FILE__
    names = isempty(ARGS) ? INSTANCES : ARGS
    results = Dict{String,Any}()
    for name in names
        results[name] = run_instance(name)
        open(joinpath(@__DIR__, "..", "results", "experiment1_$(name).json"), "w") do io
            JSON.print(io, results[name], 1)
        end
    end
end
