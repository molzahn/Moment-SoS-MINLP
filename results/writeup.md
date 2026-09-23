# Moment/SOS relaxations for AC optimal transmission switching — consolidated write-up

*2026-09-21. Consolidates `results/experiment{1,2,3,4}_findings.md`, `results/numerics_findings.md`,
`results/recertify_findings.md`, `results/merging_findings.md` and `results/rounding_scale_findings.md`.
Draft for review — not yet synced to Overleaf.*

## 1. Scope

Sparse moment/sum-of-squares (Lasserre) relaxations applied to **AC optimal transmission switching**
(AC-OTS) and single-period AC unit commitment, with two threads:

1. **Certified lower bounds** from per-clique mixed-order relaxations of the switching MINLP.
2. **Correlated randomized rounding** of the switching binaries from the relaxation's pseudo-moments,
   with continuous recovery by Ipopt.

Benchmark: the 18 PGLib-OPF v21.07 "api" cases of Taheri & Molzahn (`docs/taheri_molzahn-optimal_dcots.pdf`),
every branch switchable. Reference values are that paper's AC-OPF, O-DC-OTS and AC-OTS (Juniper) columns.

## 2. Method

AC-OTS in rectangular voltage coordinates (`docs/formulations.md`, implemented in `src/power.jl`).
Key modelling choices, all load-bearing:

- **Binaries by basis reduction.** z ∈ {0,1} is imposed by substituting z^k → z in the monomial basis, not
  by adding z² = z equality constraints, which destroy Slater's condition and wreck the order-2 interior-point solve.
- **Switching is degree 3** (p_l = z_l·P_l(V)). Order 1 cannot represent it and the bound collapses
  (73-bus 78116, 118-bus 0, 200-bus 20079 — i.e. useless).
- **Order policy, not a global order.** Only cliques containing a binary are raised to order 2
  (`binary_clique_order`); the rest stay at order 1. Roughly 40–48% of cliques go to order 2 and the
  largest PSD block is 54–75, which is what makes degree-3 constraints affordable.
- **Variants compared:** `mixed` (big-M), `mixed_nobigM`, `adj16` (adjacent-clique merging, max size 16),
  `adj16_pairs` (plus order-2 blocks on binary pairs), and from experiment 4 `base` / `conn` (radial
  fixings + connectivity cuts).
- **Certified bounds from inexact solves.** MOSEK stops at `SLOW_PROGRESS` on these problems, so the raw
  dual objective is *not* a bound. Following Oustry, D'Ambrosio, Liberti & Ruiz, *Certified and accurate
  SDP bounds for the ACOPF problem* (PSCC 2022) — `docs/Certified and Accurate SDP Bounds for ACOPF
  Problem.pdf` — `src/certify.jl` produces a bound valid at *any* dual point: box, box with tight trace
  bounds, implicit-Gram, per-block hybrid, and an L-BFGS-optimized certificate (`:smooth`, the one to
  use). We extend their method to sparse moment/SOS relaxations with polynomial equality multipliers and
  order-2 cliques; their bundle method is implemented as `bundle_certify` but is weaker in practice than
  the smoothed L-BFGS scheme. That this bound holds at any dual point is what makes every number in §3
  rigorous despite no solve reaching optimality.

## 3. Main results

Costs in $/h. "Best known" = min(ours, paper AC-OTS, paper O-DC-OTS, paper AC-OPF). Gap = (best known −
certified bound) / best known. All bounds below are **valid** — each is checked against the best known
feasible cost and no violations remain.

| case | AC-OPF | paper O-DC-OTS | paper AC-OTS | best found (ours) | best known | certified bound | gap | best variant | merged ties |
|---|---:|---:|---:|---:|---:|---:|---:|---|---:|
| 3-lmbd-api | 11,236 | 10,636 | 10,636 | 10,636 | 10,636 | 10,636 | -0.00% | exp3/adj16_pairs |  |
| 5-pjm-api | 76,377 | 75,190 | 75,190 | 75,190 | 75,190 | 73,376 | 2.41% | exp3/adj16_pairs |  |
| 14-ieee-api | 5,999 | 5,999 | 5,999 | 5,999 | 5,999 | 5,735 | 4.39% | exp3/adj16_pairs |  |
| 24-ieee-rts-api | 134,944 | 122,283 | 119,743 | 119,743 | 119,743 | 111,741 | 6.68% | exp3/adj16_pairs |  |
| 30-as-api | 4,996 | 2,797 | 2,797 | 2,797 | 2,797 | 2,767 | 1.06% | exp3/adj16 |  |
| 30-ieee-api | 18,044 | 17,939 | 17,936 | 17,940 | 17,936 | 17,287 | 3.62% | exp3/adj16_pairs |  |
| 39-epri-api | 249,672 | 246,850 | 246,723 | 246,561 | 246,561 | 245,369 | 0.48% | exp3/adj16_pairs |  |
| 57-ieee-api | 49,290 | 49,290 | 49,274 | 49,272 | 49,272 | 49,224 | 0.10% | exp3/adj16 |  |
| 60-c-api | 185,239 | 182,028 | 182,028 | 182,028 | 182,028 | 180,836 | 0.65% | exp3/adj16_pairs |  |
| 73-ieee-rts-api | 422,627 | 413,133 | 385,194 | 389,142 | 385,194 | 368,338 | 4.38% | exp3/mixed |  |
| 89-pegase-api | 130,175 | 100,702 | 100,344 | 100,503 | 100,344 | 99,858 | 0.48% | exp4/base | 19 |
| 118-ieee-api | 242,237 | 195,918 | 180,312 | 222,376 | 180,312 | 169,523 | 5.98% | exp4/base |  |
| 179-goc-api | 1,932,044 | 1,931,004 | – | 1,928,125 | 1,928,125 | 1,835,909 | 4.78% | exp3/mixed | 2 |
| 200-activ-api | 35,701 | 35,701 | 35,701 | 35,701 | 35,701 | 35,670 | 0.09% | exp4/base |  |
| 240-pserc-api | 4,640,589 | 4,627,155 | – | 4,635,816 | 4,627,155 | 4,606,386 | 0.45% | exp4/base | 55 |
| 300-ieee-api | 684,985 | 684,985 | 683,968 | 684,401 | 683,968 | 678,283 | 0.83% | exp4/conn | 2 |
| 500-goc-api | 692,407 | 692,271 | – | 676,183 | 676,183 | 667,141 | 1.34% | exp4/base |  |
| 1354-pegase-api | 1,498,271 | 1,496,750 | – | 1,498,117 | 1,496,750 | 1,487,079 | 0.65% | exp4/base | 178 |

> **Rigour note.** The five rows with merged ties (89-pegase, 179-goc, 240-pserc, 300-ieee,
> 1354-pegase) were originally bounds on a *merged* model, which is not a relaxation of the problem. All
> All five have since been re-run on the exact unmerged model and the table above already shows those
> rigorous values, so every bound here is valid for the original problem. See §7.

Notable solutions we found that beat the reference: **179-goc 1928125** and **500-goc 676183** (both
new, from the experiment-5 cardinality sweep, §5b), **39-epri 246561** (opening lines {4,6}, better than
the paper's Juniper AC-OTS 246723), **500-goc 685131** (better than the paper's O-DC-OTS 692271),
**57-ieee 49272**, **300-ieee 684401**.

Certified gaps are below 1% on 8 of 18 cases and below 0.5% on 5. The hardest are the small, highly
combinatorial ones (24-ieee-rts 6.68%, 118-ieee 5.98%, 179-goc 4.88%, 14-ieee 4.39%) — small cases are
*not* easy here, because a larger fraction of the network is switchable relative to its redundancy.

## 4. Certification matters as much as the relaxation

Re-certifying every saved relaxation with the optimized certificate tightened every case, in several
cases substantially (179-goc 6.28% → 4.88%, 60-c 1.08% → 0.65%, 200-activ 0.30% → 0.09%). Two findings:

- **The big-M vs no-big-M ranking was a certification artifact.** Big-M looked worse because its
  certificates were looser, not its relaxation. Certified properly, big-M gives the *better* bound on
  179-goc (1836683 vs 1810514). No-big-M still solves 1.5–2× faster.
- **The tight-trace-bound refinement never helped**, because the residual term dominates the box
  certificate. Negative result; the code is kept but `:smooth` is what to use.

After re-certification the certified bound is within 0.25% of the raw SDP objective on every relaxation
and within 0.03% on most — so the remaining gaps come from the *relaxations*, not from inexact solves.

## 5. Rounding

Experiments 1–2 (small enumerable cases, ground truth by enumeration) established that correlated
schemes beat independent rounding, and that order-2 pseudo-moments carry usable joint information.

At scale the picture is worse and more interesting. From experiment 4 (100 samples/scheme, evaluated on
the original network), AC-feasible fractions are 39–55% on 89-pegase and 24–37% on 240-pserc, but
**0–2% on 118, 179, 300, 500 and 1354** under every scheme including guarded sampling and island repair.
Full table in `results/experiment4_findings.md`.

**The cause is cardinality.** Independent and Gaussian rounding decide each line separately, so the
number of lines opened at once concentrates around Σ(1 − y_i) — which is 25 on 300-ieee and **724 of 1807
binaries** on 1354-pegase. Capping it (`sample_cardinality`, Gumbel top-k weighted by 1 − y_i):

| case | independent | cap 8 | cap 5 | cap 3 | cap 1 | best found by capped sampling |
|---|---:|---:|---:|---:|---:|---|
| 118-ieee | 0/50 | 0% | 6% | 22% | **58%** | 225486 (+25.05% vs best known) |
| 300-ieee | 0/30 | **47%** | 67% | 77% | **93%** | 684766 (**+0.12%**) |

On 300-ieee this is the first time rounding alone — with no 1-flip polish — has produced competitive
solutions at that size. On 118-ieee feasibility is restored but quality is not, and the reason is
diagnostic: the relaxation implies ~11 openings while the best known solution opens **37**, so no
sampler drawing from those fixed marginals can reach it. Details in `results/rounding_scale_findings.md`.

## 5b. Experiment 5: how many lines should a sample open?

`scripts/experiment5_cardinality.jl` sweeps the opening cap directly, reusing the stored marginals so no
SDP is re-solved. N = 100 (N = 60 on 500-goc), evaluated on the original network.

| case | E[#opened] | independent | best cap | feasible at that cap | best cost | vs best known |
|---|---:|---:|---|---:|---:|---:|
| 118-ieee | 10.9 | 0% | k3 | 17% | 225486.25 | +25.05% |
| **179-goc** | 39.3 | 0% | **k2** | 78% | **1928124.99** | **−0.15%** |
| 200-activ | 19.2 | 5% | k2 | 77% | 35700.87 | +0.00% |
| 240-pserc | 33.0 | 26% | k5 | 89% | 4638056.84 | +0.24% |
| 300-ieee | 25.0 | 0% | k1 | 91% | 684584.26 | +0.09% |
| **500-goc** | 50.5 | 0% | **k6** | 52% | **676182.83** | **−1.31%** |
| 89-pegase | 37.2 | 38% | k8 | 74% | 101016.95 | +0.67% |

**Two new best-known solutions, both independently verified** by re-solving the AC-OPF for that switching
configuration from scratch on the original network (`scripts/verify_solution.jl`):

- **179-goc: 1928124.99**, opening just two lines **{147, 158}** — 2851 (0.148%) below our previous
  incumbent, and below the paper's O-DC-OTS (1931004) and AC-OPF (1932044).
- **500-goc: 676182.83**, opening six lines **{131, 228, 464, 607, 629, 698}** — 8949 (1.306%) below our
  incumbent and 2.3% below the paper's O-DC-OTS (692271).

Both tighten their certified gaps, because the best known cost falls: 179-goc 4.923% → **4.783%**, and
500-goc 2.626% → **1.337%**, roughly halving the latter.

**Findings.** The effect is large and consistent: independent rounding lands at 0–38% AC-feasible,
capping at 52–91%. The useful caps are **1–8**, always far below the relaxation's own expectation, which
remains useless everywhere (0–2%). Quality is flat across neighbouring caps, so the result does not
depend on tuning the cap precisely. 118-ieee remains the exception: feasibility is restored but cost
stays +25%, because its optimum opens 37 lines while the marginals imply about 11 — the limitation there
is the pseudo-distribution, not the sampler.

**1354-pegase is excluded, for a real reason.** Its stored marginals come from the pre-fix merged model,
which was infeasible (§7), so they cannot be used for rounding; the dimension guard in the script catches
this. Regenerating them requires re-running experiment 3 or 4 on the corrected model.

## 6. Negative results worth keeping

- Tight trace bounds in the box certificate: never changed the answer.
- Island repair and guarded sampling: help where rounding already works (89-pegase 39% → 53%), do nothing
  where it does not (0–2% stays 0–2%).
- Connectivity cuts and radial fixings: improve the bound on 89-pegase and 1354-pegase, but not rounding
  on 118/179/300.
- Capping openings at the relaxation's *own* expectation is not enough (2% and 0%); the useful caps are
  far tighter, 1–8 lines.
- The "api" cases are genuinely tight: on 118-ieee, 100 of 186 single-line openings are AC-infeasible on
  their own, and feasibility is **not monotone** in the set of opened lines.

## 7. Rigour of the bounds — resolved for all 18 cases

### The problem

Low-impedance merging (`merge_zmax`, default 1e-3) was used on 5 of the 18 cases. It is **not a
relaxation**, for two independent reasons:

1. **Restriction (voltages).** Buses joined by a merged tie are forced to share one (e, f) pair. The
   original problem does not require that.
2. **Restriction (switching).** `build_power_pop` removes merged ties from `switchable`
   (`src/power.jl:50`), i.e. it forces them closed. Minimising over a *subset* of switching
   configurations can only raise the optimum.
3. **Relaxation (losses).** Merged ties are modelled as lossless, which cuts the other way.

Because the effects have opposite signs, the merged optimum is neither provably above nor below the true
one — so a "lower bound" computed on a merged model is not a valid lower bound for the original problem.
This affected **89-pegase (19 ties), 179-goc (2), 240-pserc (55), 300-ieee (2), 1354-pegase (178)**;
the other 13 cases merge nothing and were always rigorous.

Separately, merging at the default threshold made the 1354-pegase model *infeasible* by merging a tie at
its thermal limit, which is what produced the impossible 2.58e6–3.77e6 "bounds". That is fixed
independently by the self-validating merge rule (`safe_merge_exclusions`); see
`results/merging_findings.md`.

### The fix, and what it costs

The unmerged model **is** the exact problem, so it is rigorous by construction. And because
`certified_value` yields a valid bound from *any* dual point, a stalled unmerged solve still produces a
rigorous bound — merging was only ever bought for bound *quality*. So the question is purely empirical:
what does rigour cost?

All variants re-run unmerged on the three cases that are tractable without PACE
(`MAX_SDP_TIME=7200`, `CERT_OPT_TIME=900`). Best valid bound per case:

| case | merged ties | merged best bound | **rigorous best bound** | certified gap: merged → **rigorous** |
|---|---:|---:|---:|---|
| 179-goc | 2 | 1836682.93 | **1835908.57** | 4.883% → **4.923%** |
| 89-pegase | 19 | 99890.19 | **99858.29** | 0.452% → **0.484%** |
| 300-ieee | 2 | 678466.75 | **678283.45** | 0.804% → **0.831%** |
| 240-pserc | 55 | 4606997.21 | **4606386.34** | 0.436% → **0.449%** |
| 1354-pegase | 178 | 1486086.12 | **1487078.55** | 0.712% → **0.646%** (better) |

**Rigour is free.** The change in certified gap ranges from +0.040 pp to **−0.066 pp**: negligible on
four cases, and an *improvement* on 1354-pegase, where the exact model beats the merged one outright. That is the headline: the bounds were never
meaningfully dependent on merging, so the sound model can simply be used instead.

Part of even that small deficit is budget, not modelling: the merged numbers come from PACE runs with a
1800 s certificate budget, the rigorous ones from local runs with 900 s. In a *controlled*
single-variant comparison at equal budget, unmerged came out **ahead** on 179-goc (1807973.47 vs
1807472.50) and 300-ieee (678249.20 vs 678213.35, and 5× faster), and 0.025% behind on 89-pegase.

**One case where merging still earns its keep.** The 89-pegase **big-M** variant fails numerically when
unmerged: MOSEK returns a garbage point (raw 861.03) and the certificates collapse to ≈ −1.6e8. That is
still a *valid* lower bound — the guard correctly does not flag it, since a very negative bound violates
nothing — but it is useless. The no-big-M variants are unaffected and supply the best bound anyway. So
merging is **not** obsolete in general, as an earlier draft of this note suggested on one variant's
evidence; it remains a workaround for big-M conditioning on 89-pegase specifically. Everywhere else
tested, the unmerged solve is as good and often faster.

### Status per case

| case | rigorous bound | status |
|---|---|---|
| the 13 cases with 0 merged ties | as in §3 | **rigorous already** — nothing was ever merged |
| 179-goc | **1835908.57** (gap 4.923%) | **done**, all variants re-run unmerged |
| 89-pegase | **99858.29** (gap 0.484%) | **done**, all variants (big-M variant unusable, see above) |
| 300-ieee | **678283.45** (gap 0.831%) | **done**, all variants re-run unmerged |
| 240-pserc | **4606386.34** (gap 0.449%) | **done**, all variants re-run unmerged |
| 1354-pegase | **1487078.55** (gap 0.646%) | **done** — rigorous bound *better* than merged |

So **all 18 cases now have rigorous bounds** — every number in §3 is a valid lower bound for the
original problem; results files for the finished
ones are `results/recertify_<case>_unmerged.json`.

Mechanics: `MERGE_ZMAX=0` gives the rigorous model, and `RESULT_TAG` writes to a separate file so the two
sets of results do not overwrite each other, e.g.
`MERGE_ZMAX=0 RESULT_TAG=_unmerged julia --project=. scripts/recertify.jl 240-pserc-api`.

## 8. Open problems and next steps

2. **Consider deleting merging entirely.** If the unmerged solves hold up across all variants and cases,
   `merge_zmax` should default to 0 and the merge machinery becomes dead weight — a simplification that
   also removes a whole class of soundness bug.
3. **Rounding at scale.** Sweep `kmax` (1–8 is the useful range, not the relaxation's expectation).
   Re-run experiments 3/4 with `sample_cardinality`. For cases needing large coordinated rewirings
   (118-ieee, 89-pegase, 73-ieee-rts), move from resampling fixed marginals to **conditioning and
   re-solving** (randomized diving, `sample_dive`), which changes the distribution; and test conditional
   moments E[V | z] as Ipopt warm starts.
4. **Why are the 118-ieee marginals so concentrated?** E[#opened] = 11 against an optimum that opens 37.
   This is a relaxation-tightness question (big-M looseness, order policy) and links this thread to the
   untouched relaxation-selection thread (OBBT, QC envelopes, dual-/eigenvalue-guided order elevation).
5. **Still untouched:** the ML-guided relaxation-selection thread. Non-ML baselines should come first.

## 9. Infrastructure note

**PACE scratch currently has a storage fault.** 26 of 227 files under `~/scratch/Moment-SoS-MINLP` are
persistently unreadable (`Input/output error` on repeated reads, Lustre `/coda1/scratch1`), including
`Project.toml` and several `results/*.json`. Deleting and rewriting repaired some but not all — new
writes into `results/` also land unreadable, which points at a failed OST rather than anything local.
25 of the 26 have intact copies on the workstation; the only loss is one stderr log
(`hpc/logs/mm-240-pserc-api_13190482.err`). **This needs a PACE support ticket**, and the 240-pserc and
1354-pegase rigour re-runs are blocked until it is resolved.

## 10. Reproduction

Julia 1.12.7 is not on PATH:

```
~/julia-1.12.7/bin/julia --project=. scripts/<script>.jl <case>
```

- `scripts/validate_merging.jl` — merge regression test, all 18 cases (~2 min).
- `scripts/test_islands.jl` — island/bridge/cut unit tests (96/96).
- `scripts/recertify.jl <case>` — re-solve and certify; env `MERGE_ZMAX`, `RESULT_TAG`, `VARIANTS`,
  `CERT_OPT_TIME`, `MAX_SDP_TIME`.
- `scripts/experiment3_dcots.jl`, `scripts/experiment4_islanding.jl` — the main experiments.
- `scripts/analyze_experiment3.jl` — regenerates `results/experiment3_tables.md`; **fails** if any bound
  exceeds a known feasible cost.
- PACE: `hpc/submit_recertify.sh <case>`. Use `PACE_ACCOUNT=gts-dmolzahn6-fy20phase3` — the default
  account is out of allocation. Never rsync with `--delete` (the Julia depot lives inside the tree).
