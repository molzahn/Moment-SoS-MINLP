# Experiment 2: scaling sparse moment relaxations + rounding to case14/24/30

*2026-09-13. Code: `scripts/experiment2_scaling.jl` (`VARIANTS=base` or `VARIANTS=adjacent`), tables: `results/experiment2_tables.md`, raw data: `results/experiment2*_*.json`.*

## ⚠️ Update after the numerics fix (rerun 2026-09-13)

Rerun with `normalize = true`, `scale_vars = true` and certified bounds; same seeds. Side-by-side: `results/rerun_comparison.md`. Updated tables: `results/experiment2_tables.md`.

**Certified bounds and gaps:**

| instance | best certified bound | best known | certified gap |
|---|---|---|---|
| case14_ots | 8077.59 | 8081.42 | 0.05% |
| case30_uc | 959.37 | 1005.74 | 4.6% |
| case24_uc | 75736.79 | 78586.97 | 3.6% |
| case30_ots | 174.46 | 194.37 | 10.2% |
| case24_ots | 74425.89 | 78377.09 | 5.0% |

The best-known solutions are unchanged.

**Solve times:**
- The pathological solves are gone: case24_uc `adj16` 1083 s → 85 s; case30_ots `adj16_pairs` 1035 s → 133 s.
- Some small solves are slightly slower (e.g. case30_ots `mixed` 5 s → 12 s). This includes the scaled-model construction and certification. S3 is therefore resolved in its severe form.

**Rounding:**
- Essentially unchanged at this scale. S4 stands: correlated schemes give no consistent gain over independent rounding, and relax-and-round from the NLP relaxation is competitive.
- Small improvements on case24_ots from order-2 cliques around binaries: `mixed` best rounded gap 0.45% → 0.19%; `mixed_aug` 1.64% → 0.91%.
- Order-1 OTS rounding got worse (case30_ots: 38% → 0% feasible; case24_ots: 10% → 0%), consistent with the arbitrariness of order-1 OTS marginals (experiment 1 update, item 4).

---

# Original analysis (before the numerics fix)

## Setup

- **Instances** (`scripts/instances.jl`; case24/case30 are the PGLib versions shipped with PowerModels):

  | instance | variables | binaries | description |
  |---|---|---|---|
  | case14_ots | 87 | 10 | case14 with 10 switchable branches |
  | case30_uc | 77 | 6 | 2 units + 4 synchronous condensers, with no-load costs |
  | case24_uc | 146 | 33 | RTS-96, all units commitable, original Pmin and no-load costs |
  | case30_ots | 146 | 15 | 15 switchable non-bridge branches |
  | case24_ots | 188 | 15 | 15 switchable non-bridge branches |

  Only case14_ots and case30_uc are small enough to enumerate. For the others, "best known" is the best configuration found by any method below.
- **Baselines:**
  - all binaries = 1;
  - **NLP relax-and-round**: Ipopt on the POP with binaries in [0,1], then threshold rounding plus 200 independent samples.
- **Relaxation variants.** All are sparse (chordal cliques, binary reduction, SOS form, MOSEK). Linear constraints with more than 8 variables (the capacity cut) are enforced on first moments only (`global_linear = 8`).
  - `order1`
  - `mixed`: order 2 on cliques containing a binary, order 1 elsewhere (`binary_clique_order`).
  - `mixed_aug`: `mixed` plus clique augmentation (binary + bus voltages, and pairs of binaries at a bus). Augmentations are accepted greedily while every binary clique stays at ≤ 14 variables (`capped_augmentation`).
  - `adj16`: order 2 on cliques of ≤ 16 variables that contain a binary *or* a variable sharing a constraint with a binary (`adjacent_clique_order`).
  - `adj16_pairs`: `adj16` plus small extra moment blocks outside the chordal structure. There is one block per pair of binaries at a common bus: {z_i, z_j, local flows/outputs, bus voltage} (`extra_cliques`, `meta["pair_cliques"]`).
- **Rounding:** threshold, independent, Gaussian, and conditional, with 200 samples each. Every distinct configuration is solved with PowerModels AC-OPF + Ipopt (60 s CPU limit).

## Summary

| instance | all-on | NLP relax-and-round | best known | best bound (variant) | certified gap | order1 → mixed → adj16 bound gap |
|---|---|---|---|---|---|---|
| case14_ots | 8081.52 | 8081.52 | 8081.42 (enum.) | 8077.91 (mixed) | 0.04% | 100% → 0.04% → – |
| case30_uc | 1204.97 | 1005.74 | 1005.74 (enum.) | 959.37 (adj16) | 4.6% | 24.5% → 5.2% → 4.6% |
| case24_uc | 79804.96 | 78586.97 | 78586.97 | 75747.02 (adj16_pairs) | 3.6% | 8.1% → 7.8% → 3.6% |
| case30_ots | 204.97 | 194.38 | 194.37 | 175.07 (adj16_pairs) | 9.9% | 100% → 11.0% → 10.3% |
| case24_ots | 79804.96 | 78522.89 | **78377.09** | 74476.28 (adj16_pairs) | 5.0% | 49% → 7.7% → 5.1% |

## Findings

### S1 (positive): mixed-order sparse relaxations run on 150–190-variable MINLPs in seconds
- `mixed` (order 2 only on binary cliques) solves in 0.2–5 s on every instance.
- It closes the huge order-1 gaps: 100% → 11% on case30_ots, 49% → 7.7% on case24_ots, 24% → 5% on case30_uc.

### S2 (positive but costly): pushing order 2 onto cliques *adjacent* to binaries is what tightens bounds
- **Why `mixed` alone leaves the network coupling at order 1.** A binary sits in a tiny clique with its own p_g/q_g or its line flows. The bus power-balance equations that couple it to the network are assigned to large voltage cliques that stay at order 1.
- **Effect of `adj16`:**
  - case24_uc: gap 7.8% → 3.6%
  - case24_ots: 7.7% → 5.0%
  - case30_uc: 5.2% → 4.6%
  - case30_ots: 11.0% → 10.3%
- **Cost:** solve times become large and erratic, from 56 s to 1083 s. This is thread 1's core trade-off: which cliques to elevate.

### S3 (negative): solve time is dominated by MOSEK numerical difficulty, not size
- Nearly every order-2 solve ends with SLOW_PROGRESS.
- **Evidence from case24_uc:** `adj16_pairs` has *more* PSD blocks (130 vs 79) than `adj16`, yet solved in 68 s vs 1083 s, with the same bound to 0.03%. On case30_ots the pattern reversed: `adj16` took 218 s and `adj16_pairs` 1035 s.
- **Consequences:**
  - Timing comparisons between variants are unreliable until conditioning is addressed.
  - Bounds from SLOW_PROGRESS solves are only approximately valid. In experiment 1, one bound exceeded a feasible cost by 0.008%.
- **Likely culprits:**
  - explicit big-M on/off constraints with large M;
  - many redundant constraints at order 2 (big-M plus exact switching);
  - poorly scaled cost coefficients (for case24, costs in $/h with pu powers).

### S4 (negative): with non-tight relaxations at this scale, correlated rounding shows no clear advantage
- **No consistent advantage for Gaussian/conditional over independent rounding**, measured by best-of-200.
  - Gaussian was better in one case: case24_ots `mixed` (0.45% vs 1.62% gap).
  - Independent was better or equal elsewhere.
- **Pair moments rarely add information.** Clique augmentation and pair blocks create some pair moments (up to 51/528 on case24_uc), but not much correlation.
- **Versus the NLP relax-and-round baseline:**
  - *Tied* on case24_uc (78586.97) and case30_ots (194.37 vs 194.38).
  - *Better* on case24_ots: 78377.09 (from `adj16` and `adj16_pairs`, independent/Gaussian/conditional) vs 78522.89. That is 0.19% better, and 1.8% better than all lines on.
- **Moment marginals mostly add value through better feasibility rates:**
  - case30_uc: 80–86% feasible with `mixed` or `adj16` vs 25% with the NLP relaxation.
  - case24_ots: up to 40% with `adj16_pairs` vs 36% with the NLP relaxation.
  - case24_uc: only 6–15%, i.e., most sampled commitments are infeasible or islanded, which wastes AC-OPF solves.
- **Deterministic threshold rounding is poor:** it is often infeasible on UC, and on OTS it often returns all lines on.

### S5: the best-known solutions are probably not optimal on the 15–33-binary instances
- Certified gaps of 3.6–10% remain. We cannot tell how much comes from the relaxation and how much from the heuristics.
- For case30_ots in particular, the ≈10% gap is consistent across all variants, and every method converges to 194.37. That suggests a relaxation weakness rather than an unfound better configuration, but only a global MINLP solve or a tighter relaxation can decide.

## Implications / next steps

1. **Numerics first (S3).**
   - Rescale costs and variables.
   - Drop the big-M constraints where exact order-2 constraints exist.
   - Try MOSEK tolerance settings, or explicit facial reduction of the binary-reduced equalities.
   - Consider term sparsity (TSSOS-style) to shrink the order-2 blocks.
2. **Thread 1 now has a concrete target.** Choose which cliques to elevate (and whether to add pair blocks) to maximize bound improvement per second. The `adjacent_clique_order` / `capped_augmentation` / `extra_cliques` knobs are the action space. Solve-time variability means the reward must be measured carefully (e.g., over several solves).
3. **Rounding at scale needs more than one-shot sampling (S4):**
   - feasibility-aware sampling (reject islanding / capacity violations *before* AC-OPF, resample);
   - diving with the cheaper mixed relaxation;
   - local search (1- and 2-flip neighborhoods) around the best rounded configuration.
4. **Stronger ground truth** for the 15–33-binary instances: a global MINLP solve of a QC/SOC relaxation-based MISOCP (Juniper/Gurobi are installed), or long branch-and-bound runs.
5. **Larger standard test sets** (PGLib case57/118 and PGLib-UC) would require downloading data; not done yet.
