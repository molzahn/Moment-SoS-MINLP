# Experiment 1: moment relaxations + rounding on small enumerable AC-UC / AC-OTS instances

*2026-09-13. Code: `scripts/experiment1.jl`, `scripts/experiment1b_cliques.jl`. Full tables: `results/experiment1_tables.md`. Raw data: `results/experiment1_*.json`, `results/experiment1b_cliques.json`, `results/enum_*.json`.*

## ⚠️ Update after the numerics fix (rerun 2026-09-13; supersedes parts of F1–F5)

All experiments were rerun with `normalize = true`, `scale_vars = true` and certified bounds (`results/numerics_findings.md`), using the same seeds and sample counts.
- Side-by-side tables: `results/rerun_comparison.md`
- Updated full tables: `results/experiment1_tables.md`
- The findings below the line are the **original** analysis; this section lists what changed.

| instance, order 2 | certified bound (old raw) | marginals of interest | independent opt | Gaussian opt | conditional opt | dive opt |
|---|---|---|---|---|---|---|
| case5_uc | 21865.71 (21865.80) | integral | 100% | 100% | 100% | 100% |
| case5_uc_sym | 21888.35 (21891.22, invalid) | u1, u2 = 0.30, 0.40; corr −0.54 (was −0.20) | 46% | 60% (50%) | **73%** (57%) | 100% |
| case14_uc | 9042.46 | integral | 100% | 100% | 100% | 100% |
| case5_ots | 15025.97 (15071.65), gap 0.98% | z5, z6 = **0.17, 0.14** (were 0.66, 0.67) | **66%** (10%) | 68% (8%) | **73%** (13%) | 40% (20%) |
| case9_ots | 638.18 | integral | 100% | 100% | 100% | 100% |

In parentheses: values from the original run.

**What changes in the conclusions:**

1. **F4 was largely a numerical artifact.**
   - With accurate solves, the *base* order-2 relaxation of case5_ots already puts low marginals on the two phase shifters that the optimum switches off.
   - Independent rounding finds the optimum in 66% of samples (was 10%), and threshold rounding returns it directly.
   - The old "marginals 0.66 vs 0.08 depending on cliques" effect was mostly inaccurate moments, not optimal-face non-uniqueness. The general caveat about non-unique optimal faces remains valid, but it was not the main driver here.

2. **Clique augmentation (F3) still helps, now on top of good marginals.**

   | case5_ots, order 2 | independent | Gaussian | conditional |
   |---|---|---|---|
   | base | 66% | 64% | 70% |
   | bus_binaries | 83% | 87% | **92%** |
   | all_binaries | 84% | 88% | **94%** |

   With augmentation, marginals drop to 0.05 and correlated schemes add 8–10 points over independent rounding.

3. **The correlation effect (F5) is stronger and cleaner.**
   - case5_uc_sym, base cliques: corr(u1,u2) = −0.54; conditional 73–74% (experiments 1 and 1b) vs independent 46%.
   - With bus_binaries: corr(u1,u2) = **−1.00**, marginals 0.52/0.48. Gaussian and conditional rounding are optimal in **100%** of samples, vs 44% for independent rounding.
   - The pseudo-moments encode "exactly one of the two identical units" exactly.

4. **Order-1 OTS rounding is arbitrary (new negative finding).**
   - The order-1 big-M relaxations have bound 0 (case5_ots) or 75 (case9_ots) and a huge optimal face.
   - With the new scaling, the solver returns a different point on it: marginals ≈ 0.4–0.5 instead of 0.8–0.95.
   - Rounding collapsed: case9_ots threshold 100% → 0%, independent 32% → 0%.
   - Order-1 OTS marginals carry no reliable information, which reinforces F2.

5. **Diving at order 2** is now optimal in 100% of samples (5 of 5) on case5_uc and case5_uc_sym, and 40% on case5_ots.

---

# Original analysis (before the numerics fix)

## Setup

- **Instances** (definitions in `scripts/instances.jl`), all derived from PowerModels test cases:
  - `case5_uc`: 5 units
  - `case5_uc_sym`: two identical units at bus 1
  - `case14_uc`: 5 units; the three synchronous condensers carry no-load costs
  - `case5_ots`: 7 switchable branches
  - `case9_ots`: 9 switchable branches, 100 MVA ratings
- **Ground truth:** all 2^n configurations solved with PowerModels AC-OPF (polar) + Ipopt, after screening for islanding and capacity.
- **Relaxations:** sparse moment relaxation (chordal cliques, binary reduction, SOS dual form, MOSEK 11) at order 1 and order 2. Bound = SOS objective.
- **Rounding:** 200 samples each for threshold, independent, Gaussian copula, and pseudo-moment conditional. Diving uses 20 samples at order 1 and 5 at order 2. "Optimal" means within 1e-5 of the best enumerated cost.

## Relaxation bounds

| instance | opt | order 1 (gap) | order 2 (gap) | order-2 solve |
|---|---|---|---|---|
| case5_uc | 21865.76 | 21172.77 (3.17%) | 21865.80 (0.000%) | 4 s |
| case5_uc_sym | 21889.54 | 21094.75 (3.63%) | 21891.22 (−0.008%)* | 4 s |
| case14_uc | 9042.48 | 8626.42 (4.60%) | 9042.48 (0.000%) | 38 s |
| case5_ots | 15174.03 | **0.00** (100%) | 15071.65 (0.675%) | 10–60 s |
| case9_ots | 638.18 | 75.00 (88%) | 638.18 (0.000%) | 5 s |

\*The bound exceeds the enumerated optimum by 0.008%. MOSEK stopped with SLOW_PROGRESS, so this is treated as numerical inaccuracy; it needs a higher-accuracy re-solve.

## Findings

### F1 (positive): order 2 is (nearly) exact on 4 of 5 instances
On case5_uc, case14_uc and case9_ots, the order-2 marginals are integral and equal the optimal configuration, so every rounding scheme returns the optimum. This matches the AC-OPF experience in the Molzahn–Hiskens papers: low-order moment relaxations are often exact. It carries over to these small UC and OTS instances.

### F2 (negative): order-1 AC-OTS relaxations are useless with big-M on/off constraints
The exact switching constraints p = z·P(V) have degree 3, so they cannot be enforced at order 1. The degree-2 big-M replacement lets fractional z decouple flows from voltages. The resulting bound is 0 (case5_ots) or 75 vs 638 (case9_ots). Any first-order OTS approach needs a stronger lifted on/off formulation, such as QC/SOC on/off envelopes, or must go straight to order 2.

### F3 (negative, and important for thread 2): correlative sparsity removes binary–binary moments in OTS
With the default chordal cliques, **none** of the 21 (case5) or 36 (case9) pairs z_i z_j share a clique. The reason is structural: each z_l appears only in its own line's constraints. So correlated rounding reduces exactly to independent rounding.

UC instances do not have this problem, because the capacity cut Σ u_g P̄_g ≥ Σ P_d puts all u_g in one clique.

**Remedy tested (experiment 1b):** add "coupling supports" (`extra_supports`) that force binaries into a common clique:
- `bus_binaries`: binaries of branches/units incident to the same bus, plus that bus's voltage
- `all_binaries`: every binary together

Cost is modest on these cases (max clique size 13–17).

| case5_ots, order 2 | pair moments | marginals z5, z6 | independent opt | Gaussian opt | conditional opt |
|---|---|---|---|---|---|
| base cliques | 0/21 | 0.66, 0.67 | 10% | 8% | 8% |
| bus_binaries | 16/21 | 0.08, 0.08 | 75% | 82% | 83% |
| all_binaries | 21/21 | 0.07, 0.08 | 75% | 80% | 86% |

The optimum turns off both parallel phase shifters (z5 = z6 = 0).

### F4 (caveat): most of the F3 gain comes from different *marginals*, not from correlations
- **Symptom:** changing the clique structure changed the order-2 marginals of z5 and z6 from 0.66 to 0.08.
- **Why this can happen:**
  - The chordal extensions are not nested; greedy min-degree produced a *smaller* max clique.
  - The bound changed slightly (15071.65 → 15060.70). Nested cliques could only increase the bound, so this confirms the relaxations differ, and SLOW_PROGRESS tolerances add noise.
  - More importantly, when the relaxation is not exact, the optimal face is not a single point. The interior-point solver returns a point near its analytic center, a "maximum-entropy-like" pseudo-distribution, and that point depends on the formulation.
- **Correlations still help, but by less:** on top of better marginals, conditional sampling raises the optimal rate from 75% to 83–86%.
- **Research question this raises:** *which* optimal pseudo-moment vector should be sampled from? Options include selecting within the optimal face (max log-det, entropy-like regularization, or objective penalization in the spirit of Molzahn–Josz–Hiskens–Panciatici, CDC 2015), or using several solves.

### F5 (positive): cleanest evidence that correlations help comes from the symmetric UC instance
In case5_uc_sym, units 1 and 2 are identical, and exactly one of them should be committed (01101 or 10101).

| order 2 | corr(u1,u2) | marginals u1, u2 | independent opt | Gaussian opt | conditional opt |
|---|---|---|---|---|---|
| base cliques | −0.20 | 0.28, 0.41 | 46% | 51% | 54% |
| bus_binaries | **−0.73** | 0.50, 0.51 | 44% | 76% | **86%** |

- With symmetric marginals of about 0.5, independent rounding commits both units or neither half the time.
- The pseudo-moments encode "exactly one", and conditional sampling recovers it.
- All four configurations sampled are within 0.06% of optimal, so this is about *exact* optimality in a degenerate instance. It is illustrative rather than practically large.

### F6: order-1 correlations are essentially zero
At order 1, all binary correlations have magnitude ≤ 0.08 on every instance, so Gaussian and conditional rounding behave like independent rounding. Higher order is needed for informative correlations, not only for tighter bounds.

### F7: randomized diving is the strongest scheme at order 1, but expensive at order 2

| instance | order-1 dive: optimal rate | SDP solves |
|---|---|---|
| case5_uc | 40% | 34 |
| case14_uc | 25% | 50 |

On both instances, one-shot sampling from order-1 moments found the optimum in 0–14% of samples. At order 2, diving costs 30–260 s for 5 samples (each order-2 SDP takes 4–60 s, often with SLOW_PROGRESS). Closed-form conditioning costs nothing extra and was about as good whenever the joint moments existed.

## Caveats

- **Tiny instances** (5–9 binaries). The results show mechanisms, not performance.
- **Ipopt is local.** The enumerated optimum is the best local solution. It is certified only where the order-2 bound matches: case5_uc, case14_uc, case9_ots, and case5_uc_sym up to numerical tolerance.
- **case5_ots is not certified.** Its order-2 gap is 0.67%, so either the relaxation is not exact or Ipopt misses a better local solution for some configuration.
- **Seeds are fixed** (20260913). No variance across seeds has been measured yet.
- **MOSEK often terminates with SLOW_PROGRESS at order 2.** Bounds may be off by ~1e-4 relative, and returned moments may be less accurate. This needs better scaling and/or tighter tolerances.

## Suggested next steps

1. **Clique design for rounding:** augment cliques with neighboring binaries (bus- or corridor-based, or chosen from order-1 information), and measure pair coverage vs. cost on larger cases (case14_ots with 10 switchable lines, case30/case57).
2. **Choice within the optimal face:** compare IPM analytic-center moments with log-det or entropy-regularized moments and objective-penalized relaxations, judged by rounding quality.
3. **Scaling:** term sparsity (TSSOS-style) and better numerics (scaling, tolerances) to remove SLOW_PROGRESS.
4. **Stronger order-1 OTS formulations** (QC/SOC on/off), since F2 makes plain Shor + big-M useless.
5. **Baselines:** relax-and-round from the NLP relaxation (Gómez et al. 2026) and MISOCP/QC-based heuristics, on larger instances where enumeration is impossible (use a best-known solution from a MINLP solver or long runs).
6. **Report variability** across seeds and load perturbations.
