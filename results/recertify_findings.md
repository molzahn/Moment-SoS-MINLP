# Re-certification of experiments 3 and 4 (September 15, 2026)

`scripts/recertify.jl` re-solves every saved relaxation and certifies it with the box certificate, the box certificate
with tight trace bounds, the implicit-Gram certificate, the per-block hybrid, and the optimized certificate (smoothed
soft-min with L-BFGS, starting from the MOSEK dual). Jobs ran on PACE with MOSEK, 4 threads.

"Loss" is the raw SDP objective minus the certified bound. The old loss is measured against the saved certified bound;
the new loss against the optimized certificate.

| case | variant | raw | old certified | new certified | loss old → new |
|---|---|---:|---:|---:|---|
| 24-ieee-rts | adj16_pairs | 111819.69 | 111220.25 | 111740.58 | 591 → 79 |
| 30-ieee | adj16_pairs | 17291.55 | 17101.85 | 17286.70 | 188 → 5 |
| 39-epri | adj16_pairs | 245369.03 | 244874.05 | 245368.87 | 495 → 0.2 |
| 60-c | adj16_pairs | 180902.27 | 180023.10 | 180836.29 | 876 → 66 |
| 73-ieee-rts | mixed (big-M) | 368345.10 | 367152.40 | 368338.12 | 1193 → 7 |
| 73-ieee-rts | mixed_nobigM | 368244.24 | 368109.94 | 368243.67 | 134 → 0.6 |
| 89-pegase | mixed_nobigM | 99914.69 | 99293.18 | 99890.19 | 622 → 25 |
| 118-ieee | mixed (big-M) | 169579.98 | 168241.54 | 169371.10 | 1338 → 209 |
| 118-ieee | mixed_nobigM | 169535.35 | 169350.89 | 169523.05 | 185 → 12 |
| 179-goc | mixed (big-M) | 1841292.30 | 1809725.20 | 1836682.93 | 31567 → 4609 |
| 200-activ | mixed_nobigM | 35677.99 | 35594.66 | 35670.27 | 83 → 8 |
| 300-ieee | mixed (big-M) | 678537.48 | none (NaN) | 678179.49 | – → 358 |
| 300-ieee | exp4/conn | 678537.93 | 677648.24 | 678466.75 | 890 → 71 |

Across all 49 re-solved relaxations, the optimized certificate recovers roughly 80–99.9% of the gap between the box certificate
and the raw objective. Tight trace bounds alone never changed the box certificate: the residual term dominates.
Implicit/hybrid recovers most of the loss, and optimization recovers most of the rest (1–5 min per relaxation).

## Updated certified gaps (best bound over all variants vs. best known = min(ours, paper AC-OTS, paper O-DC-OTS))

| case | old gap | new gap | raw-objective gap (not certified) |
|---|---:|---:|---:|
| 5-pjm | 2.42% | 2.41% | 2.41% |
| 14-ieee | 4.49% | 4.39% | 4.39% |
| 24-ieee-rts | 6.86% | 6.68% | 6.62% |
| 30-as | 1.31% | 1.06% | 1.05% |
| 30-ieee | 4.22% | 3.62% | 3.59% |
| 39-epri | 0.51% | 0.48% | 0.48% |
| 57-ieee | 0.17% | 0.10% | 0.09% |
| 60-c | 1.08% | 0.65% | 0.62% |
| 73-ieee-rts | 4.44% | 4.38% | 4.37% |
| 89-pegase* | 0.71% | 0.45% | 0.35% |
| 118-ieee | 6.08% | 5.98% | 5.95% |
| 179-goc | 6.28% | 4.88% | 4.64% |
| 200-activ | 0.30% | 0.09% | 0.06% |
| 300-ieee | 0.92% | 0.80% | 0.79% |

\* The 89-pegase bound is for the low-impedance-merged model, while the best known value is the paper's unmerged AC-OTS.

Notes:
- Certification losses were large mainly for the big-M formulation and for adj16/pairs. Once certified tightly, big-M
  becomes the best bound on 179-goc (1836683 vs. 1810514 without big-M) and the 300-ieee big-M variant, which had no bound before, now gets 678179.
  So the earlier conclusion that no-big-M gives the better bound was partly a certification artifact. No-big-M is
  still 1.5–2× faster.
- The certificate is now within 0.25% of the raw objective on every relaxation (below 0.03% on most). The remaining gaps come from the
  relaxation, not from the inexact solve.
- 57-ieee adj16_pairs ran out of memory at 16G. It was resubmitted at 48G, together with 240-pserc, 500-goc and 1354-pegase.

---

## Addendum, September 20, 2026: the four late cases, and an invalid bound on 1354-pegase

The tables above were written on September 15 at 10:18, before the resubmitted 48–160 GB jobs
(57-ieee at 48G, plus 240-pserc, 500-goc and 1354-pegase) finished at 10:45–12:14. Their results were
never folded in. Doing that now:

| case | variant | raw | optimized certificate | comment |
|---|---|---:|---:|---|
| 57-ieee | adj16_pairs | 49229.35 | 49220.93 | completes the 57-bus row |
| 240-pserc | mixed (big-M) | 4607420 | 4605168 | big-M beats no-big-M here too |
| 240-pserc | mixed_nobigM | 4607220 | 4606997 | best 240-bus bound; gap ≈ 0.49% |
| 500-goc | mixed_nobigM | 667245 | 667141 | first bound for 500-goc; vs best found 685131 → gap ≈ 2.6% |
| 1354-pegase | all variants | 2.58M–3.77M | 2.58M–3.77M | **invalid — see below** |

### 1354-pegase: invalid bound, diagnosed and now fixed

The re-certification reported 2580996 (mixed), 2612454 (no big-M) and 3772084 (conn) against a known
feasible cost of 1496750. A lower bound cannot exceed a feasible solution. **These numbers must not be
reported.**

The cause has been found and fixed: low-impedance merging at the default `merge_zmax = 1e-3` merged a
tie that sits exactly at its thermal limit, which made the 1354-pegase model *infeasible*, so the SDP
was relaxing the wrong problem. The certificate machinery was not at fault. Full diagnosis, evidence
and fix are in **`results/merging_findings.md`**.

**Re-run completed** (PACE job 13409920, 3 h 24 m, 47 GB peak). All four variants now return valid
bounds, the best being **1486086.12** from `exp4/conn`, a certified gap of **0.712%** against the best
known 1496750. This is the first valid bound ever obtained for this case — the original run gave `NaN`.
The 1354-pegase row in the table above should read:

| case | gap before re-certification | gap after |
|---|---|---|
| 1354-pegase | – (NaN / invalid) | **0.71%** |

MOSEK still stops at `SLOW_PROGRESS`, but `rel_gap` is 1.2e-11 with residuals ~1e-10, and the
pathological `box_tight_rho` of -2.6e6 is now a sane +1.46e6. The certification loss (1488307 raw →
1486086 certified, 0.15%) is in line with the other cases.

### Caveat that applies to the tables above

A bound computed on a low-impedance-merged model is **not a rigorous lower bound** for the original
network: merging both restricts (merged buses share a voltage) and relaxes (tie losses are dropped).
This affects every case with merged ties — **89-pegase (19), 179-goc (2), 240-pserc (55), 300-ieee (2)
and 1354-pegase (178)** — not only 89-pegase, which is the only one flagged above. See
`results/merging_findings.md` for what would be needed to remove the caveat.
