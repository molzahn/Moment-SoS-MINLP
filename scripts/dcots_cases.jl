# Test cases and reported results from Taheri & Molzahn, "AC-Informed DC Optimal Transmission Switching
# Problems via Parameter Optimization" (docs/taheri_molzahn-optimal_dcots.pdf), Tables I and II.
# Data: PGLib-OPF v21.07 "api" cases (data/pglib_api/, copied from the authors' OTS project data).
# All branches are switchable. Costs in $/h; AC-OTS = PowerModels AC-OTS solved with Juniper (local).

using MomentMINLP, PowerModels, JSON

const API_DIR = joinpath(@__DIR__, "..", "data", "pglib_api")

"name => (file stem, AC-OPF cost, O-DC-OTS cost, AC-OTS (Juniper) cost, AC-OTS time [s], AC-OTS opened lines)"
const DCOTS_PAPER = [
    ("3-lmbd-api", "case3_lmbd", 11236, 10636, 10636, 0.971, [3]),
    ("5-pjm-api", "case5_pjm", 76377, 75190, 75190, 0.521, [3]),
    ("14-ieee-api", "case14_ieee", 5999, 5999, 5999, 0.100, Int[]),
    ("24-ieee-rts-api", "case24_ieee_rts", 134944, 122283, 119743, 4.446, [8, 9, 24, 30, 34, 35]),
    ("30-as-api", "case30_as", 4996, 2797, 2797, 7.141, [11, 22]),
    ("30-ieee-api", "case30_ieee", 18044, 17939, 17936, 3.230, [12, 23, 26]),
    ("39-epri-api", "case39_epri", 249672, 246850, 246723, 5.234, [6, 31]),
    ("57-ieee-api", "case57_ieee", 49290, 49290, 49274, 14.17, [13, 34, 52, 62, 73]),
    ("60-c-api", "case60_c", 185239, 182028, 182028, 15.79, [3, 4, 16]),
    ("73-ieee-rts-api", "case73_ieee_rts", 422627, 413133, 385194, 611.0,
        [1, 24, 27, 34, 35, 36, 37, 50, 65, 73, 74, 75, 76, 80, 113, 114, 119]),
    ("89-pegase-api", "case89_pegase", 130175, 100702, 100344, 1867.0, Int[]),     # 49 lines opened (see paper)
    ("118-ieee-api", "case118_ieee", 242237, 195918, 180312, 3098.0, Int[]),      # 37 lines opened (see paper)
    ("179-goc-api", "case179_goc", 1932044, 1931004, nothing, nothing, Int[]),
    ("200-activ-api", "case200_active", 35701, 35701, 35701, 1.000, Int[]),
    ("240-pserc-api", "case240_pserc", 4640589, 4627155, nothing, nothing, Int[]),
    ("300-ieee-api", "case300_ieee", 684985, 684985, 683968, 35248.0, Int[]),
    ("500-goc-api", "case500_goc", 692407, 692271, nothing, nothing, Int[]),
    ("1354-pegase-api", "case1354_pegase", 1498271, 1496750, nothing, nothing, Int[]),
]

paper_row(name) = DCOTS_PAPER[findfirst(r -> r[1] == name, DCOTS_PAPER)]

"""
    best_known_cost(name) -> Float64

Lowest AC-feasible cost known for `name`: the paper's AC-OPF (all lines closed, feasible by construction),
its O-DC-OTS and AC-OTS values, and any `best_found` recorded by experiments 3 and 4.

A certified lower bound must never exceed this. When one does, the relaxation is not a relaxation of this
problem — 1354-pegase reported 2.58e6-3.77e6 against a known feasible 1.498e6 because low-impedance merging
had made its model infeasible (see `safe_merge_exclusions` and results/merging_findings.md).
"""
function best_known_cost(name)
    row = paper_row(name)
    vals = Float64[Float64(v) for v in (row[3], row[4], row[5]) if v !== nothing]
    for f in ("experiment3_$(name).json", "experiment4_$(name).json")
        p = joinpath(@__DIR__, "..", "results", f)
        isfile(p) || continue
        v = get(JSON.parsefile(p), "best_found", nothing)
        (v isa Real && isfinite(v)) && push!(vals, Float64(v))
    end
    return minimum(vals)
end

"""
    check_bound(value, best_known; tol = 1e-6) -> Bool

`false` when `value` is a lower bound that exceeds a known feasible cost, which is impossible. Callers warn
and record the flag rather than throwing: these bounds cost hours of PACE time, so the run should finish and
save, and the tabulation step is where an invalid row is refused.
"""
check_bound(value, best_known; tol::Float64 = 1e-6) =
    !(value isa Real) || !isfinite(value) || value <= best_known + tol * max(abs(best_known), 1.0)
parse_api(stem) = PowerModels.parse_file(joinpath(API_DIR, "pglib_opf_$(stem)__api.m"))

"Buses joined by non-transformer branches with |z| below this (p.u.) share a voltage (env MERGE_ZMAX; 0 = off)."
const MERGE_ZMAX = parse(Float64, get(ENV, "MERGE_ZMAX", "1e-3"))

"""
AC-OTS POP with every branch switchable (as in the paper), except low-impedance ties (|z| < `merge_zmax`),
whose end buses share a voltage and which stay closed (see `build_power_pop`). Up to 118 buses this only
affects 89-pegase (19 ties). Configurations are still evaluated on the original network (`pop.meta["data"]`).
"""
function load_dcots_instance(name; merge_zmax = MERGE_ZMAX, kwargs...)
    row = paper_row(name)
    d = parse_api(row[2])
    return build_power_pop(d; switchable = sort(parse.(Int, collect(keys(d["branch"])))), name = name,
        merge_zmax = merge_zmax, kwargs...), d
end
