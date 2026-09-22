# Low-impedance merging: an invalid bound on 1354-pegase, its cause, and the fix

*2026-09-20. Code: `safe_merge_exclusions` in `src/power.jl`; regression test
`scripts/validate_merging.jl`.*

## Symptom

Re-certification reported these **lower bounds** for `1354-pegase-api`:

| variant | reported bound |
|---|---:|
| exp3/mixed | 2580996.38 |
| exp3/mixed_nobigM | 2612453.62 |
| exp4/base | 2612453.62 |
| exp4/conn | 3772084.28 |

The best known feasible cost is **1496750** (the paper's O-DC-OTS; our own best found is 1498116.94,
opening lines [646, 770, 1718, 1921]). A lower bound cannot exceed a feasible cost, so every one of
these numbers was wrong — by 72% to 152%.

## The certificate machinery was not at fault

`certified_value` (`src/certify.jl:156`) forms `θ[tcol]` and then subtracts only non-positive terms
(residual penalties and `ρ_k · min(λ_min, 0)`). It therefore returns a valid lower bound **for
whatever POP produced the certificate**, for any `θ`, no matter how badly MOSEK converged. "MOSEK
stalled" cannot by itself produce a bound above the optimum. The POP itself had to be wrong.

## Root cause: merging made the model infeasible

Low-impedance merging forces the end buses of a tie to share one voltage and remodels the tie's flow
as lossless. Rebuilding 1354-pegase at several thresholds, with every switch closed:

| `merge_zmax` | merged ties | all-closed AC-OPF of the POP |
|---|---:|---|
| 0 | 0 | `LOCALLY_SOLVED` 1498271.03 (= the paper's AC-OPF exactly) |
| 1e-4 | 0 | `LOCALLY_SOLVED` 1498271.03 |
| **1e-3 (the default)** | **184** | **`LOCALLY_INFEASIBLE`** |

An SDP relaxation of an infeasible POP can take any value, which is exactly what was observed: three
different variants produced three wildly different numbers, all above the incumbent, all at
`SLOW_PROGRESS`, and the *original* (pre-recertification) run returned `NaN` outright.

### Mechanism: a merged tie at its thermal limit

At the true (unmerged) AC-OPF solution, measured across the 184 merged ties:

- max |ΔV_m| across a tie = **0.00194 pu** (median 1.33e-4), max |Δθ| = 0.00647 rad — so the
  shared-voltage restriction is very nearly harmless;
- **max flow / rate_a = 1.000** — a merged tie sits *exactly* at its thermal limit.

A binding constraint has no slack to absorb the remodelling. Scaling **only the merged-tie ratings**
confirms it is the thermal limits and nothing else:

| merged-tie `rate_a` | all-closed AC-OPF of the merged POP |
|---|---|
| ×1.00 | `LOCALLY_INFEASIBLE` |
| ×1.02 | `LOCALLY_INFEASIBLE` |
| ×1.10 | `LOCALLY_SOLVED` 1497037.28 |

Merging is not trivially at fault in the obvious way: no merged group anywhere (89-pegase,
240-pserc, 1354-pegase) has an empty shared voltage window, i.e. merging never forces together buses
whose [vmin, vmax] ranges are disjoint. And the incumbent's opened lines are not merged ties, so the
incumbent *switching pattern* was never excluded. The failure is purely the binding tie ratings.

Per this project's own standard, Ipopt `LOCALLY_INFEASIBLE` is not a proof of infeasibility. The
rating experiment localises the cause regardless, and the fix is validated by construction below.

## Fix: a merge rule that validates itself

`safe_merge_exclusions(data; zmax, tol)` replaces the bare |z| threshold. It solves the all-closed
AC-OPF of the unmerged network once, then, while the merged model is infeasible or off by more than
`tol` (default 1%), un-merges the most heavily loaded merged tie (|S| / rate_a at the unmerged
solution) and retries. It terminates in the worst case with nothing merged, and needs no new
threshold to be guessed. `build_power_pop(...; merge_safe = true)` is the default; the ties dropped
from the merge are recorded in `pop.meta["merge_excluded"]`.

On 1354-pegase it un-merges **6 of 184** ties, in this order:

| tie | 1708 | 1709 | 780 | 1028 | 643 | 1411 |
|---|---|---|---|---|---|---|
| loading |S|/rate_a | 1.000 | 0.999 | 0.982 | 0.958 | 0.731 | 0.694 |

leaving 178 merged. The all-closed AC-OPF then gives **1497994.29** against the paper's 1498271
(−0.018%). The case gains 6 binaries (1807 → 1813), since an un-merged tie is an ordinary switchable
branch again. On every other case the rule excludes nothing, so previous models are unchanged — in
particular 89-pegase keeps all 19 of its merged ties, which it needs: they carry |z| ≈ 2.2e-4 and
merging them is what cured its SDP stall.

## Regression test: `scripts/validate_merging.jl`

For each paper case, the merged POP with every switch closed must solve and match the paper's AC-OPF
cost to within 1%. All 18 now pass (36/36 assertions, ~2 min).

| case | merged ties | un-merged by validation | all-closed AC-OPF | vs paper |
|---|---:|---:|---:|---:|
| 3-lmbd / 5-pjm / 14-ieee / 24-ieee-rts / 30-as / 30-ieee / 39-epri / 57-ieee / 60-c / 73-ieee-rts / 118-ieee / 200-activ / 500-goc | 0 | 0 | — | ≤ 0.006% |
| 89-pegase | 19 | 0 | 129734.63 | −0.338% |
| 179-goc | 2 | 0 | 1932012.17 | −0.002% |
| 240-pserc | 55 | 0 | 4640357.48 | −0.005% |
| 300-ieee | 2 | 0 | 684963.51 | −0.003% |
| **1354-pegase** | **178** | **6** | **1497994.29** | **−0.018%** |

Before the fix 1354-pegase was the only failure, and it failed outright (infeasible).

## Other guards added

- `best_known_cost(name)` and `check_bound` (`scripts/dcots_cases.jl`) — the lowest feasible cost
  known for a case, from the paper's AC-OPF / O-DC-OTS / AC-OTS columns and from `best_found` in the
  saved experiment 3 and 4 results. Note the paper's O-DC-OTS column is populated for **all 18**
  cases (1354-pegase: 1496750); only the AC-OTS column is missing for the four largest, so a
  reference cost is always available.
- `scripts/recertify.jl`, `scripts/experiment3_dcots.jl`, `scripts/experiment4_islanding.jl` now
  record `best_known` and a per-variant `bound_valid`, warn loudly, and keep an invalid bound out of
  the summary and out of the reported gap. They **warn rather than throw**: these bounds cost hours
  of PACE time, so the run should still finish and save.
- `scripts/analyze_experiment3.jl` writes `**INVALID**` in place of such a bound and then **fails**,
  so an invalid number cannot reach a table silently.
- `scripts/recertify.jl` now also persists the solver and certificate diagnostics it was discarding
  (`gram_min_eig`, `residual_correction`, `psd_correction`, `sos_residual_max`, `min_eig_moment`,
  `min_eig_localizing`, `max_eq_residual`, `iterations`, `primal_obj`, `dual_obj`, `rel_gap`) plus
  the merge state. Their absence was why this bug could not be diagnosed from the saved JSON at all.

## Open caveat: a bound from a merged model is still not a rigorous lower bound

This is separate from the bug above and is **not** fixed by the new rule. Merging simultaneously

- **restricts** the feasible set (the merged buses are forced to share one voltage), and
- **relaxes** it (the ties' losses are dropped; their flows become lossless),

so the merged optimum is neither provably above nor provably below the original one. Every bound
reported for a case with merged ties therefore carries this caveat: **89-pegase (19 ties), 179-goc
(2), 240-pserc (55), 300-ieee (2) and 1354-pegase (178)**. The 89-pegase headline gap of 0.45% rests
on it. Empirically merging looks net-relaxing — the merged all-closed AC-OPF comes out at or below the
unmerged one on all five cases (89-pegase −0.338%, the rest ≤ 0.02%) — but that is evidence, not proof.

There is a **third** non-rigorous effect, easy to miss: `build_power_pop` also removes merged ties from
`switchable` (`src/power.jl:50`), forcing them closed. Minimising over a subset of switching
configurations can only raise the optimum, so this is a restriction too — independent of the voltage
sharing, and it applies to 178 branches on 1354-pegase.

### Resolved empirically: rigour is free, and merging now looks obsolete

The unmerged model *is* the exact problem, so it is rigorous by construction; and since
`certified_value` is valid at any dual point, even a stalled unmerged solve yields a rigorous bound.
Merging was only ever bought for bound *quality*. Like-for-like re-runs (same variant
`exp3/mixed_nobigM`, `MAX_SDP_TIME=5400`, `CERT_OPT_TIME=600`, same machine):

| case | ties | merged raw / certified | time | unmerged raw / certified | time | Δ certified |
|---|---:|---:|---:|---:|---:|---:|
| 179-goc | 2 | 1813862.13 / 1807472.50 | 162 s | **1813901.24 / 1807973.47** | 180 s | **+500.97 better** |
| 300-ieee | 2 | 678372.22 / 678213.35 | 1785 s | **678379.30 / 678249.20** | **350 s** | **+35.85 better** |
| 89-pegase | 19 | 99914.97 / 99883.17 | 168 s | **99919.85 / 99858.29** | **146 s** | −24.88 (−0.025%) |

On two of three the rigorous bound is *better*, and on 300-ieee it is 5× faster. Critically, **the
89-pegase stall did not reproduce**: the unmerged solve finished in 146 s, faster than merged. Merging
was introduced precisely because that case stalled on |z| ≈ 2.2e-4 branches, so the most likely
explanation is that the later normalization and variable-scaling work (`results/numerics_findings.md`,
default since 2026-09-13) cured the conditioning problem merging was working around.

### All-variant result: rigour costs at most 0.04 pp

Re-running **every** variant unmerged on the three cases tractable without PACE:

| case | merged best | rigorous best | certified gap merged → rigorous |
|---|---:|---:|---|
| 179-goc | 1836682.93 | 1835908.57 | 4.883% → 4.923% |
| 89-pegase | 99890.19 | 99858.29 | 0.452% → 0.484% |
| 300-ieee | 678466.75 | 678283.45 | 0.804% → 0.831% |
| 240-pserc | 4606997.21 | 4606386.34 | 0.436% → 0.449% |
| 1354-pegase | 1486086.12 | **1487078.55** | 0.712% → **0.646%** (better) |

**All 18 cases now have rigorous bounds** (13 never merged anything; 5 re-run on the exact model). The
change in certified gap spans +0.040 pp to −0.066 pp — free on four cases, and on 1354-pegase the exact
model gives a strictly *better* bound than the merged one, which is the cleanest possible refutation of
the idea that merging was buying anything.

**Correction to the single-variant conclusion above:** merging is *not* obsolete in general. The
89-pegase **big-M** variant fails numerically when unmerged — MOSEK returns raw 861.03 and the
certificates collapse to ≈ −1.6e8 (still a valid bound, just useless). The no-big-M variants are
unaffected and give the best bound anyway, so this does not change any reported number, but it means
`merge_zmax` cannot simply be defaulted to 0 and the merge machinery deleted. Merging remains a
workaround for big-M conditioning on 89-pegase specifically.

Use `MERGE_ZMAX=0 RESULT_TAG=_unmerged` to produce rigorous results without overwriting the merged ones.

## Status: resolved

`1354-pegase-api` was re-run on PACE with the repaired merge rule (job 13409920, 3 h 24 m, 47 GB peak,
`gts-dmolzahn6-fy20phase3`). **All four variants now return valid bounds**, and the bound-validity guard
reported no violations:

| variant | raw | certified bound | gap vs best known (1496750) |
|---|---:|---:|---:|
| exp3/mixed (big-M) | 1488573.44 | 1482716.77 | 0.938% |
| exp3/mixed_nobigM | 1488338.07 | 1485994.32 | 0.719% |
| exp4/base | 1488338.07 | 1485994.32 | 0.719% |
| **exp4/conn** | 1488307.30 | **1486086.12** | **0.712%** |

This is the **first valid bound ever obtained for 1354-pegase**: the original run produced `NaN`, and the
re-certification produced the impossible 2.58e6-3.77e6 above. The four numbers at the top of this note
remain wrong and must never be reported.

Notes on the repaired run:
- 178 ties merged, the same 6 un-merged by validation ([643, 780, 1028, 1411, 1708, 1709]) as locally.
- MOSEK still reports `SLOW_PROGRESS` with `primal_status = UNKNOWN_RESULT_STATUS`, but `rel_gap` is now
  1.2e-11 and the equality/PSD residuals are ~1e-10 — i.e. the solve is accurate, and the earlier
  wildly-negative `box_tight_rho` is gone (now +1.457e6, a sane value close to raw). The remaining ~0.1%
  certification loss (raw 1488307 → certified 1486086) is in line with every other case.
- The connectivity-cut variant gives the best bound here, marginally ahead of no-big-M; big-M is worst.
  That matches 89-pegase, where cuts also helped, and contrasts with 118/179/300, where they did not.
