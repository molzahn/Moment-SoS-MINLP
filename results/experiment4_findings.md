# Experiment 4: islanding-aware rounding, repair, and connectivity cuts

*Script `scripts/experiment4_islanding.jl`; raw data `results/experiment4_<case>.json`. N = 100 samples
per scheme, configurations evaluated on the original network with PowerModels ACP + Ipopt/MA97.*

Experiment 4 asked whether the 0–2% AC-feasible sample rates seen above 118 buses in experiment 3 were
caused by **islanding**, and whether they could be fixed by (a) evaluating with islands allowed rather
than requiring connectivity, (b) repairing island-creating samples, (c) guarding the sampler so it never
opens a line that would create an unservable island, and (d) strengthening the relaxation with radial
fixings and connectivity cuts.

Two relaxations per case: `base` (no big-M) and `conn` (no big-M + radial lines fixed closed + cuts
Σ_{l∈δ(S)} z_l ≥ 1 for bus sets that cannot balance power alone).

## AC-feasible fraction of 100 samples

| case | variant | binaries | radial fixed | cuts | E[#opened] | indep | +repair | gauss | cond | cond+rep | guarded |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 89-pegase | base | 191 | 0 | 0 | 35.8 | 39% | 39% | 38% | 49% | 54% | 43% |
| 89-pegase | conn | 189 | 2 | 113 | 41.1 | **50%** | **53%** | 42% | **55%** | **55%** | **54%** |
| 118-ieee | base | 186 | 0 | 0 | 10.1 | 0% | 0% | 0% | 1% | 1% | 0% |
| 118-ieee | conn | 182 | 4 | 186 | 10.4 | 0% | 0% | 0% | 0% | 0% | 0% |
| 179-goc | base | 261 | 0 | 0 | 39.3 | 0% | 0% | 0% | 0% | 1% | 0% |
| 179-goc | conn | 258 | 3 | 175 | 42.2 | 0% | 1% | 0% | 0% | 0% | 0% |
| 200-activ | base | 245 | 0 | 0 | 19.2 | 4% | 4% | 2% | 2% | 3% | 2% |
| 200-activ | conn | 222 | 23 | 217 | 21.0 | 1% | 1% | 1% | 2% | 2% | 3% |
| 240-pserc | base | 393 | 0 | 0 | 33.0 | 27% | 27% | 24% | 31% | 32% | **37%** |
| 240-pserc | conn | 392 | 1 | 199 | 36.4 | 24% | 25% | 33% | 27% | 30% | 32% |
| 300-ieee | base | 409 | 0 | 0 | 25.0 | 0% | 0% | 0% | 2% | 2% | 1% |
| 300-ieee | conn | 365 | 44 | 315 | 28.4 | 0% | 0% | 0% | 0% | 0% | 1% |
| 500-goc | base | 728 | 0 | 0 | 50.5 | 0% | 0% | 1% | 0% | 0% | 0% |
| 500-goc | conn | 708 | 20 | 598 | 50.8 | 1% | 1% | 0% | 0% | 0% | 0% |
| 1354-pegase | base | 1807 | 0 | 0 | 724.1 | 0% | 0% | 0% | 0% | 0% | 0% |
| 1354-pegase | conn | 1607 | 200 | 737 | 660.4 | 0% | 0% | 0% | 0% | 0% | 0% |

## Findings

**1. The islanding hypothesis is only partly right, and the machinery does not rescue the hard cases.**
Repair and guarding help where rounding already works (89-pegase 39% → 53-55%, 240-pserc 27% → 37%) and
do essentially nothing where it does not (118, 179, 300, 500, 1354 stay at 0-2% under every scheme).

**2. Why: the failure mode differs by case.** From `infeasible_reasons` on 100 independent `base` samples:

| case | AC locally infeasible / solver failure | rejected by island screen |
|---|---:|---:|
| 118-ieee | 78 | 22 |
| 200-activ | 78 | 18 |
| 500-goc | 95 | 5 |
| 240-pserc | 72 | 1 |
| 89-pegase | 58 | 3 |
| 179-goc | 18 | 82 |
| 300-ieee | 46 | 54 |
| 1354-pegase | 0 | 100 |

Islanding dominates on 179, 300 and 1354; AC infeasibility dominates on 118, 200 and 500. But repairing
islands is not sufficient in either regime — repaired samples still fail the AC problem — which is why
`+repair` barely moves the columns above.

**3. Connectivity cuts and radial fixings help the bound where the network has real bridges.** They fix
44 radial lines on 300-ieee and 200 on 1354-pegase, and improve the 89-pegase bound and its feasible
fraction. They did not improve rounding on 118/179/300. On 1354-pegase the `conn` variant gives the best
certified bound (see `results/recertify_findings.md`).

**4. The real driver is how many lines are opened at once, not connectivity per se.** E[#opened] is 724
of 1807 binaries on 1354-pegase — the relaxation's marginals want to open 40% of the network, and every
such sample is rejected. This observation motivated `sample_cardinality`, which caps the number of
simultaneous openings and lifts 300-ieee from 0% to 93% AC-feasible; see
**`results/rounding_scale_findings.md`**, which supersedes this experiment's conclusions about how to fix
rounding at scale.

**5. Stressed "api" cases are genuinely tight.** On 118-ieee, 100 of the 186 single-line openings are
AC-infeasible on their own, and feasibility is **not monotone** in the set of opened lines — opening a
second line can restore feasibility lost by the first. Any screen based on single-line reasoning will
therefore be incomplete.

## Caveats

Ipopt returning `LOCALLY_INFEASIBLE` is not a proof of infeasibility, so the "AC locally infeasible"
column is an upper bound on genuine infeasibility; some of those samples may be feasible points Ipopt
failed to find. The island screen, by contrast, rejects only islands that provably cannot balance active
power (Σpmax < Σ(pd + gs·vmin²)), so its rejections are sound.
