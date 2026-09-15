# Re-certify the relaxations of experiments 3 and 4 with the improved certificates (src/certify.jl).
#   julia --project=. scripts/recertify.jl <case name, e.g. 118-ieee-api>
#
# For every relaxation variant saved in results/experiment3_<case>.json (and experiment4_<case>.json), the
# SDP is re-solved with the same settings and certified with:
#   box (original), box with tight trace bounds, implicit Gram, per-block hybrid (step 1; also available when
#   MOSEK returns no feasible point, step 2), and the optimized certificate (step 3; smoothed L-BFGS).
# MOSEK's stopping point can differ between runs, so the re-solve's raw objective is reported next to the
# saved one. Writes results/recertify_<case>.json.
include("dcots_cases.jl")
using JSON, Printf

const MAX_SDP_TIME = parse(Float64, get(ENV, "MAX_SDP_TIME", "43200"))
const OPT_TIME = parse(Float64, get(ENV, "CERT_OPT_TIME", "1800"))
const CUT_SET = parse(Int, get(ENV, "CUT_SET", "2"))
const CUT_BLOCK_MAX = parse(Int, get(ENV, "CUT_BLOCK_MAX", "6"))

jsonsafe(x::AbstractFloat) = isfinite(x) ? x : nothing
jsonsafe(x::AbstractDict) = Dict(k => jsonsafe(v) for (k, v) in x if !startswith(string(k), "_"))
jsonsafe(x::AbstractVector) = map(jsonsafe, x)
jsonsafe(x) = x

"(label, pop, solve kwargs, saved bound) for every saved relaxation of the case."
function variants(name)
    out = Any[]
    params = Dict{String,Any}("MSK_DPAR_OPTIMIZER_MAX_TIME" => MAX_SDP_TIME)
    f3 = joinpath(@__DIR__, "..", "results", "experiment3_$(name).json")
    if isfile(f3)
        saved = JSON.parsefile(f3)["relaxations"]
        popM = first(load_dcots_instance(name))
        popN = haskey(saved, "mixed_nobigM") ? first(load_dcots_instance(name; bigM_switching = false)) : nothing
        for lbl in ("mixed", "mixed_nobigM", "adj16", "adj16_pairs")
            haskey(saved, lbl) || continue
            pop = lbl == "mixed_nobigM" ? popN : popM
            kw = Dict{Symbol,Any}(:global_linear => 8, :solver_params => params)
            if lbl in ("mixed", "mixed_nobigM")
                kw[:order] = binary_clique_order(pop)
            else
                kw[:order] = adjacent_clique_order(pop; maxsize = 16)
                lbl == "adj16_pairs" && (kw[:extra_cliques] = pop.meta["pair_cliques"])
            end
            push!(out, ("exp3/" * lbl, pop, kw, get(saved[lbl], "bound", nothing), get(saved[lbl], "raw_bound", nothing)))
        end
    end
    f4 = joinpath(@__DIR__, "..", "results", "experiment4_$(name).json")
    if isfile(f4)
        saved = JSON.parsefile(f4)["relaxations"]
        for lbl in ("base", "conn")
            haskey(saved, lbl) || continue
            pop = lbl == "base" ? first(load_dcots_instance(name; bigM_switching = false)) :
                first(load_dcots_instance(name; bigM_switching = false, fix_radial = true, conn_cuts = CUT_SET))
            kw = Dict{Symbol,Any}(:global_linear => 8, :solver_params => params, :order => binary_clique_order(pop),
                :out_of_graph_tags => ["conn_cut", "conn_cut_wide"],
                :extra_cliques => filter(c -> length(c) <= CUT_BLOCK_MAX, get(pop.meta, "conn_cut_cliques", Vector{Vector{Int}}())))
            push!(out, ("exp4/" * lbl, pop, kw, get(saved[lbl], "bound", nothing), get(saved[lbl], "raw_bound", nothing)))
        end
    end
    return out
end

function run_case(name)
    row = paper_row(name)
    res = Dict{String,Any}("case" => name, "paper_acots" => row[5], "paper_odcots" => row[4], "variants" => Dict{String,Any}())
    out = joinpath(@__DIR__, "..", "results", "recertify_$(name).json")
    for (lbl, pop, kw, saved_bound, saved_raw) in variants(name)
        t = @elapsed rel = solve_moment_relaxation(pop; kw..., certify = :optimize,
            bundle_params = Dict(:method => :smooth, :maxtime => OPT_TIME))
        i = rel.info
        g(k) = get(i, k, nothing)
        o = Dict{String,Any}("status" => string(rel.status), "primal_status" => string(rel.primal_status), "time" => t,
            "solve_time" => rel.solve_time, "saved_bound" => saved_bound, "saved_raw" => saved_raw, "raw" => g("raw_bound"),
            "box" => g("certified_box"), "box_tight_rho" => g("certified_box_tight_rho"), "implicit" => g("certified_implicit"),
            "hybrid" => g("certified_hybrid"), "optimized" => g("certified_optimized"), "bound" => rel.bound,
            "optimize_time" => g("optimize_time"), "optimize_evals" => g("optimize_evals"),
            "unrepresented_monomials" => g("n_unrepresented_monomials"), "square_groups" => g("n_square_groups"))
        fmt(x) = (x === nothing || (x isa Real && !isfinite(x))) ? "–" : @sprintf("%.2f", x)
        @printf("%-16s %-18s %-14s raw %s (saved %s) | box %s | box ρ* %s | implicit %s | hybrid %s | optimized %s | saved certified %s | %.0fs (opt %.0fs)\n",
            name, lbl, o["status"], fmt(o["raw"]), fmt(saved_raw), fmt(o["box"]), fmt(o["box_tight_rho"]), fmt(o["implicit"]),
            fmt(o["hybrid"]), fmt(o["optimized"]), fmt(saved_bound), t, something(o["optimize_time"], 0.0))
        flush(stdout)
        res["variants"][lbl] = o
        open(out, "w") do io
            JSON.print(io, jsonsafe(res), 1)
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    foreach(run_case, ARGS)
end
