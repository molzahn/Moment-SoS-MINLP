# Validate POP formulations against PowerModels (polar AC-OPF) and report relaxation bounds.
#   julia --project=. scripts/validate_formulations.jl
using MomentMINLP, PowerModels, JuMP, Printf
const DATA = joinpath(@__DIR__, "..", "data")

for case in ["case5", "case9", "case14"]
    data = PowerModels.parse_file(joinpath(DATA, case * ".m"))
    pm = PowerModels.solve_opf(data, ACPPowerModel, ipopt_optimizer())
    pop = build_power_pop(data; name = case)
    nlp = solve_nlp(pop)
    @printf("\n%s  PowerModels ACP: %.4f (%s)   rectangular POP NLP: %.4f (%s)\n", case,
            pm["objective"], pm["termination_status"], nlp.objective, nlp.status)
    # all-switchable OTS and all-commitable UC models with binaries fixed to 1 must agree with OPF
    ots = build_power_pop(data; switchable = parse.(Int, collect(keys(data["branch"]))))
    uc = build_power_pop(data; commitable = parse.(Int, collect(keys(data["gen"]))))
    for (lbl, q) in (("OTS(z=1)", ots), ("UC(u=1)", uc))
        r = solve_nlp(q; fixed = Dict(v => 1.0 for v in binvars(q)))
        @printf("   %-9s NLP: %.4f (%s)\n", lbl, r.objective, r.status)
    end
    for ord in (1, 2)
        case == "case14" && ord == 2 && continue
        rel = solve_moment_relaxation(pop; order = ord)
        @printf("   order %d: bound %.4f  gap %.3f%%  status %s  cliques %d (max %d)  moments %d  time %.1fs  λ2/λ1 max %.2e\n",
                ord, rel.bound, 100 * (pm["objective"] - rel.bound) / pm["objective"], rel.status,
                length(rel.cliques), maximum(length.(rel.cliques)), rel.n_moments, rel.solve_time,
                maximum(rel.rank_ratio))
    end
end
