# Validate the low-impedance merge on every DC-OTS paper case.
#   julia --project=. scripts/validate_merging.jl [case ...]
#
# Merging forces the end buses of a low-impedance tie to share a voltage and remodels its flow as a
# lossless one. Where a tie sits at its thermal limit that is not harmless: with a plain |z| threshold the
# 1354-pegase model comes out infeasible, and the relaxation built on it reported "lower bounds" of
# 2.58e6-3.77e6 against a known feasible 1.498e6. `safe_merge_exclusions` (src/power.jl) now keeps such
# ties out of the merge, and this script is the regression test that would have caught the bug.
#
# For each case the merged POP with every switch closed must solve, and must match the paper's AC-OPF cost
# to within MERGE_TOL (default 1%). Note that agreement here does not make a bound from the merged model a
# rigorous lower bound for the original network: merging both restricts (shared voltages) and relaxes
# (tie losses dropped). See results/merging_findings.md.

include("dcots_cases.jl")
using Test, Printf

const TOL = parse(Float64, get(ENV, "MERGE_TOL", "0.01"))

cases = isempty(ARGS) ? [r[1] for r in DCOTS_PAPER] : ARGS

@testset "merged model reproduces the AC-OPF of every paper case" begin
    for nm in cases
        row = paper_row(nm)
        pop, _ = load_dcots_instance(nm)
        ties, excl = pop.meta["merged_ties"], pop.meta["merge_excluded"]
        res = solve_nlp(pop; fixed = Dict(v => 1.0 for v in pop.meta["binary_vars"]))
        rel = (res.objective - row[3]) / row[3]
        @printf("%-17s ties %3d (unmerged by validation %2d)  %-18s %14s  paper %10d  %s\n",
            nm, length(ties), length(excl), string(res.status),
            isfinite(res.objective) ? @sprintf("%.2f", res.objective) : "Inf", row[3],
            isfinite(rel) ? @sprintf("%+.3f%%", 100 * rel) : "--")
        flush(stdout)
        @test res.feasible
        @test isfinite(res.objective) && abs(rel) <= TOL
    end
end
