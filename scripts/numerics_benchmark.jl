# Numerics benchmark: effect of formulation / solver options on MOSEK behaviour for order-2 relaxations.
#   julia --project=. scripts/numerics_benchmark.jl [case_id ...]
include("instances.jl")
using Printf, JSON

const BIGM_TAGS = ["sw_flow_bigM", "sw_angle_bigM", "sw_flow_onoff"]

"Benchmark relaxations: (id, pop builder, relaxation keyword arguments)."
function bench_cases()
    opf5() = (build_power_pop(parse_case("case5")), (order = 2,))
    function inst(name, kind)
        pop, _ = load_instance(name)
        kw = kind == :full ? (order = 2,) :
             kind == :mixed ? (order = binary_clique_order(pop), global_linear = 8) :
             kind == :adj16 ? (order = adjacent_clique_order(pop; maxsize = 16), global_linear = 8) :
             (order = adjacent_clique_order(pop; maxsize = 16), global_linear = 8, extra_cliques = pop.meta["pair_cliques"])
        return pop, kw
    end
    return [("case5_opf_o2", opf5),
        ("case14_uc_o2", () -> inst("case14_uc", :full)),
        ("case5_ots_o2", () -> inst("case5_ots", :full)),
        ("case30_ots_mixed", () -> inst("case30_ots", :mixed)),
        ("case24_ots_adj16", () -> inst("case24_ots", :adj16)),
        ("case24_uc_adj16pairs", () -> inst("case24_uc", :adj16pairs))]
end

const CONFIGS = [
    ("baseline", (;)),
    ("normalize", (normalize = true,)),
    ("norm+drop", (normalize = true, skip_high_order_tags = BIGM_TAGS)),
    ("norm+quotient", (normalize = true, quotient_basis = true)),
    ("norm+scale", (normalize = true, scale_vars = true)),
    ("norm+scale+drop", (normalize = true, scale_vars = true, skip_high_order_tags = BIGM_TAGS)),
    ("norm+scale+drop+quot", (normalize = true, scale_vars = true, skip_high_order_tags = BIGM_TAGS, quotient_basis = true)),
]

getn(d, k) = get(d, k, NaN)

if abspath(PROGRAM_FILE) == @__FILE__
    sel = isempty(ARGS) ? nothing : Set(ARGS)
    cfgsel = haskey(ENV, "CONFIGS") ? Set(split(ENV["CONFIGS"], ",")) : nothing
    out = Dict{String,Any}()
    outfile = joinpath(@__DIR__, "..", "results", "numerics_benchmark.json")
    isfile(outfile) && (out = JSON.parsefile(outfile))
    for (id, builder) in bench_cases()
        (sel === nothing || id in sel) || continue
        pop, kw = builder()
        @printf("\n=== %s\n", id)
        @printf("  %-20s %-14s %7s %5s %13s %13s %9s %9s %9s %9s  %s\n", "config", "status", "time", "iter", "raw bound",
            "certified", "SOS res", "λmin(X)", "relgap", "λmin(M)", "dropped")
        for (cfg, ckw) in CONFIGS
            (cfgsel === nothing || cfg in cfgsel) || continue
            rel = solve_moment_relaxation(pop; kw..., ckw...)
            i = rel.info
            @printf("  %-20s %-14s %7.1f %5s %13.4f %13.4f %9.1e %9.1e %9.1e %9.1e  %s\n", cfg, rel.status, rel.solve_time,
                string(get(i, "iterations", "–")), getn(i, "raw_bound"), rel.bound, getn(i, "sos_residual_max"),
                getn(i, "gram_min_eig"), getn(i, "rel_gap"), getn(i, "min_eig_moment"), i["dropped"])
            flush(stdout)
            out["$(id)/$(cfg)"] = Dict("status" => string(rel.status), "solve_time" => rel.solve_time,
                "bound" => rel.bound, "info" => i)
            open(outfile, "w") do io
                JSON.print(io, out, 1)
            end
        end
    end
end
