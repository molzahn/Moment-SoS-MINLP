# Re-solve the relaxations of experiments 1 and 2 with the numerics defaults (normalize + scale_vars,
# certified bounds) and compare with the previously reported (raw, unnormalized) bounds.
#   julia --project=. scripts/recheck_bounds.jl
include("instances.jl")
using Printf, JSON

const RES = joinpath(@__DIR__, "..", "results")
rows = Any[]
function run!(inst, variant, pop, kw, old)
    rel = solve_moment_relaxation(pop; kw...)
    i = rel.info
    @printf("%-12s %-12s %-14s %8.1fs  certified %12.4f  raw %12.4f  old %s  SOS res %.1e\n", inst, variant, rel.status,
        rel.solve_time, rel.bound, get(i, "raw_bound", NaN), old === nothing ? "–" : @sprintf("%.4f", old),
        get(i, "sos_residual_max", NaN))
    flush(stdout)
    push!(rows, Dict("instance" => inst, "variant" => variant, "status" => string(rel.status), "solve_time" => rel.solve_time,
        "bound" => rel.bound, "raw_bound" => get(i, "raw_bound", nothing), "old_bound" => old,
        "sos_residual_max" => get(i, "sos_residual_max", nothing)))
    open(joinpath(RES, "recheck_bounds.json"), "w") do io
        JSON.print(io, rows, 1)
    end
end

for name in ["case5_uc", "case5_uc_sym", "case14_uc", "case5_ots", "case9_ots"]
    ex = JSON.parsefile(joinpath(RES, "experiment1_$(name).json"))
    pop, _ = load_instance(name)
    for ord in (1, 2)
        run!(name, "order$(ord)", pop, (order = ord,), ex["orders"][string(ord)]["bound"])
    end
end
for name in ["case30_uc", "case24_uc", "case30_ots", "case24_ots"]
    pop, _ = load_instance(name)
    b = JSON.parsefile(joinpath(RES, "experiment2_$(name).json"))
    a = JSON.parsefile(joinpath(RES, "experiment2_adjacent_$(name).json"))
    run!(name, "mixed", pop, (order = binary_clique_order(pop), global_linear = 8), b["relaxations"]["mixed"]["bound"])
    ord = adjacent_clique_order(pop; maxsize = 16)
    run!(name, "adj16", pop, (order = ord, global_linear = 8), a["relaxations"]["adj16"]["bound"])
    run!(name, "adj16_pairs", pop, (order = ord, global_linear = 8, extra_cliques = pop.meta["pair_cliques"]),
        a["relaxations"]["adj16_pairs"]["bound"])
end
