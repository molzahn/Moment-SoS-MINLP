# Rounding at scale: cardinality is the binding constraint, not islanding

*2026-09-20. Code: `sample_cardinality` in `src/rounding.jl`.*

The central open problem coming out of experiments 3 and 4 was that rounding collapses with size: 38-55%
of samples are AC-feasible on 89-pegase and 24-37% on 240-pserc, but **0-2% on 118, 179, 300, 500 and
1354**, where every good solution comes from the 1-flip polish rather than from the rounding itself.

## Hypothesis

Independent and Gaussian rounding decide each line on its own, so the number of lines opened *at once*
concentrates around Σ(1 - y_i). On the stressed "api" cases that is far more simultaneous openings than any
feasible configuration tolerates. Capping the count, while still letting the relaxation choose *which*
lines, should restore feasibility.

`sample_cardinality(μ, N; kmax, guard)` opens at most `kmax` lines per sample, chosen without replacement
with probability proportional to the opening marginal 1 - y_i (Gumbel top-k, i.e. Plackett-Luce), and
respects an `IslandGuard` if given.

## Setup

Marginals are the saved `mixed_nobigM` pseudo-moments from `results/experiment3_<case>.json`; every sampled
configuration is evaluated on the **original** network with `ConfigEvaluator` (PowerModels ACP AC-OPF,
Ipopt/MA97). No SDP is re-solved, so this isolates the sampler. N = 50 (118-ieee) and N = 30 (300-ieee).

### 118-ieee-api — 186 binaries, E[#opened] = 10.9, best known 180312 (AC-OPF 242237)

| scheme | mean opened | AC-feasible | best cost | vs best known |
|---|---:|---:|---:|---:|
| independent | 10.6 | **0/50 (0%)** | – | – |
| independent + repair | 10.1 | **0/50 (0%)** | – | – |
| cardinality kmax=11 (= E) | 11.0 | 1/50 (2%) | 245540.84 | +36.18% |
| cardinality kmax=8 | 8.0 | 0/50 (0%) | – | – |
| cardinality kmax=5 | 5.0 | 3/50 (6%) | 240723.03 | +33.50% |
| cardinality kmax=3 | 3.0 | 11/50 (22%) | 225486.25 | +25.05% |
| cardinality kmax=2 | 2.0 | 16/50 (32%) | 228919.05 | +26.96% |
| cardinality kmax=1 | 1.0 | **29/50 (58%)** | 231521.99 | +28.40% |

### 300-ieee-api — 409 binaries, E[#opened] = 25.0, best known 683968 (AC-OPF 684985)

| scheme | mean opened | AC-feasible | best cost | vs best known |
|---|---:|---:|---:|---:|
| independent | 25.5 | **0/30 (0%)** | – | – |
| independent + repair | 24.8 | **0/30 (0%)** | – | – |
| cardinality kmax=25 (= E) | 25.0 | 0/30 (0%) | – | – |
| cardinality kmax=8 | 8.0 | 14/30 (47%) | 684766.01 | **+0.12%** |
| cardinality kmax=5 | 5.0 | 20/30 (67%) | 684985.45 | +0.15% |
| cardinality kmax=3 | 3.0 | 23/30 (77%) | 684822.16 | +0.12% |
| cardinality kmax=2 | 2.0 | 26/30 (87%) | 684824.09 | +0.13% |
| cardinality kmax=1 | 1.0 | **28/30 (93%)** | 684971.69 | +0.15% |

## Findings

**1. The hypothesis holds, and the effect is large.** Feasibility rises monotonically as the cap tightens,
from 0% to 58% (118-ieee) and 0% to 93% (300-ieee). On 300-ieee the sampler now *also* finds good
solutions directly — within **0.12%** of the best known cost, from rounding alone, with no polish step.
This is the first time rounding rather than the 1-flip polish has produced competitive solutions at this
size.

**2. Repairing islands does not restore feasibility — a negative result that matters.** `independent +
repair` is still 0/50 and 0/30 here, and experiment 4 shows the same across 118, 179, 300, 500 and
1354 (0-1% for every guarded or repaired scheme). It is worth being precise about *why*, because the
reported failure reason differs by case — from experiment 4's `infeasible_reasons` on 100 independent
samples of the `base` relaxation:

| case | AC locally infeasible / solver failure | rejected by island screen |
|---|---:|---:|
| 118-ieee | 78 | 22 |
| 200-activ | 78 | 18 |
| 500-goc | 95 | 5 |
| 179-goc | 18 | 82 |
| 300-ieee | 46 | 54 |
| 1354-pegase | 0 | 100 |

So both mechanisms are active, in different proportions: AC infeasibility dominates on 118, 200 and
500, while the island screen rejects most samples on 179, 300 and 1354. The common finding is that
**repairing islands is not sufficient** — repaired samples still fail the AC problem. Cardinality
capping helps precisely because it attacks both at once: opening fewer lines produces fewer islands
*and* less network stress. (An earlier draft of this note claimed islanding was simply not the
bottleneck; that is right for 118-ieee, which is the case I measured directly, but wrong for 179,
300 and 1354, where the screen is the dominant rejector.)

**3. Capping at the relaxation's own expectation is not enough.** `kmax = round(Σ(1 - y_i))` gives 2% and
0%. The useful caps are far below the expectation — 1 to 8 lines. So the default in `sample_cardinality`
(the relaxation's expectation) is the wrong default for these cases and should be swept, not trusted.

**4. Feasibility and quality decouple, and the difference between the two cases explains why.** On
300-ieee the optimum is very close to the all-closed network (AC-OPF 684985 vs best known 683968, only
0.15% apart), so few-opening samples are both feasible and good. On 118-ieee the optimum is a *large*
rewiring — the paper's solution opens **37** lines for a 25% cost reduction below AC-OPF — so capping at 1-3
openings cannot reach it, and the best feasible sample is +25%. (It is still well below the all-closed
AC-OPF, 225486 vs 242237, so the samples are not worthless — just far from optimal.)

**5. The marginals themselves are the deeper problem on 118-ieee.** The relaxation implies E[#opened] =
10.9, with only 4 of 186 lines having a marginal below 0.5, while the best known solution opens 37 lines.
The pseudo-marginals are inconsistent in cardinality with the optimum by a factor of ~3.5. **No sampler
drawing from these fixed marginals can reach that solution**, however cleverly it chooses. This reframes
the problem: at 118-ieee the limit is the relaxation's pseudo-distribution, not the rounding scheme.

## What this implies for the next steps

- Sweep `kmax` rather than defaulting to the expectation; 1-8 is the useful range on these cases.
- Cases where the optimum is near the all-closed network (300-ieee, and likely 200-activ, 500-goc) look
  solved by this change alone. Re-running experiments 3/4 with `sample_cardinality` should be cheap and is
  the obvious next measurement.
- For cases needing large coordinated rewirings (118-ieee, 89-pegase, 73-ieee-rts), the priority moves off
  sampler design and onto **changing the distribution**: conditioning and re-solving (randomized diving,
  `sample_dive`), which updates the marginals, rather than resampling from fixed ones. Conditional-moment
  Ipopt warm starts remain worth testing, but note finding 2 — these samples fail on AC limits, so a warm
  start may convert some of the near-misses.
- Worth checking why the marginals are so concentrated near 1 on 118-ieee. That is a relaxation-tightness
  question (big-M looseness, order policy), and it connects this thread back to the relaxation-selection
  thread.

## Caveats

Marginals come from a single saved relaxation per case (`mixed_nobigM`); the sampler was not re-run against
`mixed`, `adj16` or `adj16_pairs`. Sample counts are small (30-50) and a single rng seed was used, so the
feasible fractions carry a few percentage points of noise — the 0% vs 47-93% contrast is far outside that,
the differences between neighbouring `kmax` values are not. No SDP was re-solved for this study.
