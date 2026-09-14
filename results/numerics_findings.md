# Numerics of order-2 sparse moment relaxations: diagnosis and fixes

*2026-09-13. Code: `scripts/numerics_benchmark.jl`, `scripts/recheck_bounds.jl`. Logs and data:*
- *round 1: `results/numerics_benchmark_log_round1.txt`, `results/numerics_benchmark_round1.json`*
- *rounds 2–3: `results/numerics_benchmark_log.txt`, `results/numerics_benchmark_log_round3.txt`, `results/numerics_benchmark.json`*
- *recheck: `results/recheck_bounds_log.txt`, `results/recheck_bounds.json`*

## Problem

In experiments 1–2, almost every order-2 solve ended with MOSEK status SLOW_PROGRESS.
- Solve times varied erratically (56–1083 s for similar models).
- In one case (case5_uc_sym), the reported "lower bound" exceeded a feasible cost.

## Diagnostics added (`rel.info`)

- barrier iterations; primal and dual objective and their gap
- feasibility of the recovered pseudo-moments (smallest eigenvalue of moment and localizing matrices, largest equality residual)
- **residuals of the SOS problem itself:** coefficient-matching residuals and the smallest Gram-matrix eigenvalue
- **a certified lower bound** that corrects the SOS objective for those residuals using the variable box (derivation in `docs/formulations.md` §5.1)

## What was tested

Benchmark relaxations:
- case5 AC-OPF, order 2
- case14_uc, order 2
- case5_ots, order 2
- case30_ots, `mixed`
- case24_ots, `adj16`

| option | effect |
|---|---|
| `normalize` (each constraint divided by its max coefficient) | SOS residuals drop 100–1000×. case5 OPF and case14_uc become OPTIMAL. Certified bounds become meaningful. |
| `scale_vars` (all variables scaled to [−1, 1]) | Residual correction shrinks further; certified bound within ≤ 0.03% of the raw bound on every case. case5 OPF: 1.2 s → 0.5 s. |
| `skip_high_order_tags` (drop big-M on/off constraints in order-2 cliques) | **Weakens** the bound (case5_ots 15026 → 15019; case24_ots 74375 → 74007), with no consistent speed-up. Not adopted. |
| `quotient_basis` (remove pivot monomials of degree-2 equalities from order-2 bases; exact reformulation) | No measurable effect on status, time, or bound. Not adopted. The near-singularity of the moment matrices at optimality comes mostly from low rank (tight relaxations), not from the equality constraints. |
| MOSEK `INTPNT_SOLVE_FORM` primal/dual | No change (identical iterates). |

**New defaults:** `normalize = true`, `scale_vars = true`, certified `rel.bound`.

Switching instances still end with SLOW_PROGRESS, but their certified bounds are now within 0.03% of the raw SOS objective. The status is harmless once the bound is certified.

## Corrected bounds (old = raw unnormalized bound reported in experiments 1–2)

| instance | relaxation | old bound | certified bound (new) | solve time old → new | comment |
|---|---|---|---|---|---|
| case5_uc | order 2 | 21865.80 | 21865.71 | 4 s → 3.5 s | old bound was 0.0002% above the optimum 21865.76 |
| case5_uc_sym | order 2 | **21891.22** | **21888.35** | 4 s → 7 s | old value exceeded the optimum 21889.54 (invalid); new gap 0.005% |
| case14_uc | order 2 | 9042.48 | 9042.46 | 38 s → 28 s | |
| case5_ots | order 2 | **15071.65** | **15025.97** | 10–60 s → 11 s | gap vs 15174.03 is **0.98%**, not 0.67% |
| case9_ots | order 2 | 638.18 | 638.18 | 5 s → 5 s | exact |
| case24_uc | adj16 | 75728.19 | 75716.19 | **1083 s → 55 s** | |
| case24_uc | adj16_pairs | 75747.02 | 75736.79 | 68 s → 67 s | certified gap 3.63% |
| case30_ots | adj16 | 174.38 | 173.80 | 218 s → 85 s | |
| case30_ots | adj16_pairs | 175.07 | 174.46 | **1035 s → 91 s** | certified gap 10.2% |
| case24_ots | adj16 | 74415.51 | 74374.89 | 56 s → 51 s | |
| case24_ots | adj16_pairs | 74476.28 | 74425.89 | 72 s → 53 s | certified gap 5.04% |

Order-1 and `mixed` bounds were essentially unaffected (changes < 0.02%).

## Takeaways

1. **Raw SOS objectives from SLOW_PROGRESS solves must not be reported as bounds.**
   - Before normalization, the correction needed to certify them exceeded 20% (case5_ots: raw 15072, certified 11713).
   - After normalization and scaling, it is below 0.03%.
2. **Conditioning, not size, caused the erratic 1000-s solves.** With normalization and scaling, the worst cases drop to 55–91 s.
3. **Experiment 1's qualitative conclusions stand.** Order 2 is exact or nearly exact on case5_uc, case5_uc_sym, case14_uc and case9_ots, and correlations help on case5_uc_sym. The case5_ots order-2 gap is about 1%.
4. **Experiment 2's certified gaps are essentially unchanged** (3.6% / 10.2% / 5.0% on case24_uc / case30_ots / case24_ots).
5. **Rounding results in experiments 1–2 used moments from the old formulation.** Re-running rounding with the new defaults is still to do. Marginals can move, because the optimal face is not unique (experiment 1, F4).

## Remaining numerical issues / ideas

- **OTS relaxations still stall.** Likely causes are the degree-3 switching equalities together with big-M constraints, and a non-unique optimal face.
- **Candidates to try:**
  - Replace big-M with tighter on/off bounds derived from angle-difference limits.
  - Term sparsity to shrink the order-2 blocks.
  - Exploit that $z\,p = p$ for switched flows, i.e. add the redundant equalities $z\,p - p = 0$, which hold exactly.
