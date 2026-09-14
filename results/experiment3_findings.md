# Experiment 3: AC-OTS on the test cases of Taheri & Molzahn (DC-OTS paper)

*Scripts: `scripts/dcots_cases.jl`, `scripts/experiment3_dcots.jl`, `scripts/analyze_experiment3.jl`. Full tables: `results/experiment3_tables.md`. Log: `results/experiment3_log.txt` (3-lmbd was run separately in the foreground). Status: nine small cases (3–60 buses) done; 73/89/118 and larger not yet run.*

## Setup

- **Data:** PGLib-OPF v21.07 "api" cases, the same files the paper used (`data/pglib_api/`).
  - AC-OPF costs reproduce Table II for all 18 cases.
  - The paper's line numbers are PowerModels branch ids; its reported AC-OTS topologies reproduce its costs.
- **Problem:** AC-OTS with **every branch switchable**, as in the paper. Rectangular-coordinate POP with exact degree-3 switching constraints plus big-M constraints (`build_power_pop`).
- **Relaxations** (SOS form, MOSEK, normalized and variable-scaled, certified bounds; see `docs/formulations.md` §5.1):
  - `mixed`: order 2 on cliques that contain a binary, order 1 elsewhere.
  - `adj16`: order 2 on cliques up to 16 variables adjacent to binaries.
  - `adj16_pairs`: adj16 plus extra moment blocks for pairs of binaries.
  - Linear constraints enter through first moments only (`global_linear = 8`).
- **Rounding from each relaxation:**
  - threshold
  - independent (N = 200)
  - Gaussian copula (N = 200)
  - pseudo-moment conditional (N = 200)
- **Evaluating configurations:** each is screened for connectivity, then solved with PowerModels ACP AC-OPF + Ipopt (local).
- **"Polish":** a separately reported 1-flip local search around the best rounded configuration. It flips the least certain lines first, with a budget of 40 new AC-OPF solves per relaxation.
- **Reference numbers:** the paper's AC-OTS costs come from Juniper, a *local* MINLP solver, so they are not proven optimal. Our certified bounds bound how far they can be from optimal.

## Results

| case | AC-OPF | paper O-DC-OTS | paper AC-OTS | best found (ours) | vs paper AC-OTS | best certified bound | certified gap | peak memory | wall time |
|---|---|---|---|---|---|---|---|---|---|
| 3-lmbd | 11235.68 | 10636 | 10636 | 10635.95 | 0 | 10635.95 | **0 (proven optimal)** | – | 0.5 min |
| 5-pjm | 76377.41 | 75190 | 75190 | 75190.29 | 0 | 73369.15 | 2.42% | 2.1 GB | 1 min |
| 14-ieee | 5999.36 | 5999 | 5999 | 5999.36 | 0 | 5729.78 | 4.49% | 7.7 GB | 5.5 min |
| 24-ieee-rts | 134943.63 | 122283 | 119743 | 119743.13 | 0 | 111523.24 | 6.86% | 10.3 GB | 8 min |
| 30-as | 4996.20 | 2797 | 2797 | 2797.47 | 0 | 2760.26 | 1.31% | 10.0 GB | 7 min |
| 30-ieee | 18043.92 | 17939 | 17936 | 17939.51 | +0.02% | 17179.82 | 4.22% | 10.1 GB | 7 min |
| 39-epri | 249672.32 | 246850 | 246723 | **246561.20** | **−0.07%** | 245299.51 | **0.51%** | 8.7 GB | 5.5 min |
| 57-ieee | 49290.36 | 49290 | 49274 | **49271.85** | −0.004% | 49188.32 | **0.17%** | 13.5 GB | 17 min |
| 60-c | 185239.00 | 182028 | 182028 | 182027.63 | 0 | 180064.72 | 1.08% | 15.5 GB | 13 min |

The certified gap is (best known − best bound) / best known. For 30-ieee the best known solution is the paper's Juniper cost.

### Positive findings

1. **Optimality certificates for the paper's AC-OTS solutions.**
   - **3-lmbd:** the relaxation is exact, so the paper's solution is globally optimal.
   - **Within 1.4% of optimal:** 57-ieee (0.17%), 39-epri (0.51%), 60-c (1.08%) and 30-as (1.31%).
   - **Larger gaps:** 5-pjm (2.4%), 30-ieee (4.2%), 14-ieee (4.5%) and 24-ieee-rts (6.9%). These are much looser.
2. **Better solutions than the paper on two cases.**
   - **39-epri:** opening lines {4, 6} costs 246561.20, compared with 246722.96 for the paper's {6, 31}. Both costs were re-verified with an independent PowerModels AC-OPF solve.
   - **Pure rounding finds the 39-epri topology.** Independent and conditional rounding from `mixed`, Gaussian from `adj16`, and conditional from `adj16_pairs` all reach it without local search.
   - **57-ieee:** opening {13, 34, 62, 73, 77} costs 49271.85. The improvement is only 0.004%, and it comes after polishing.
3. **The best known solution is matched or beaten on 8 of 9 cases.** The exception is 30-ieee, 3.5 $/h (0.02%) above Juniper; that topology needs two extra lines opened simultaneously, which 1-flip search cannot reach.
4. **Pure rounding without polish reaches the best known cost on 3-lmbd, 5-pjm, 14-ieee, 30-as (adj16), 39-epri and 60-c (adj16, adj16_pairs).** On 24-ieee-rts the best rounded solution is +0.31% (adj16_pairs), which already beats the paper's O-DC-OTS (+2.12%). Polishing then reaches the paper's Juniper topology.
5. **The mixed relaxation is the best value.** It takes 27–49 s per SDP. Its certified bound is within 0.2% of the best bound on 7 of 9 cases (all except 14-ieee at 0.7% and 24-ieee-rts at 0.6%), and it *is* the best bound on 3-lmbd, 30-as, 57-ieee and 60-c. Its rounding is weaker than adj16 on 30-as and 60-c, however.

### Negative findings

1. **The all-switchable problem gives much weaker bounds than instances with a few switchable lines.**
   - case14 has a 4.5% gap here vs 0.05% in experiment 2 (few switchable lines).
   - 24-ieee-rts has a 6.9% gap even though the paper's solution is plausibly optimal.
   - Order 2 on the cliques containing binaries is not enough to close these gaps.
2. **Larger relaxations are held back by solver accuracy, not relaxation strength.**
   - Every solve ends in MOSEK `SLOW_PROGRESS` except the three 3-lmbd solves and 5-pjm `mixed`.
   - The raw SOS objective increases monotonically from mixed to adj16 to adj16_pairs on every case.
   - The certification correction usually grows faster than the raw objective: 24-bus 13 → 198 → 591, 39-bus 73 → 52 → 495, 60-bus 270 → 1329 → 876.
   - As a result, `adj16_pairs` has the *lowest* certified bound on 30-as, 30-ieee, 39-epri and 57-ieee. On 60-c, `adj16` is the lowest.
   - Improving SDP accuracy (solver tolerances, better conditioning, facial reduction of the binary blocks) should be worth more than larger cliques.
3. **Correlated rounding gives no systematic gain over independent rounding.**
   - Gaussian/conditional rounding helps on 24-ieee-rts adj16_pairs (+0.31–0.35% vs +2.12%) and on 30-as mixed (Gaussian +9.6% vs independent +22.8%).
   - It is worse on 30-as adj16 (independent reaches the optimum, Gaussian +9.6%) and on 60-c adj16 (conditional +0.12%).
   - This mirrors experiment 2.
4. **Threshold rounding is poor.** It is up to 79% above best (30-as), and it is infeasible for all 57-ieee variants and for 24-ieee-rts `adj16`. Relaxation marginals are rarely below 0.5, and at most 4 lines have a marginal below 0.5 in any case.
5. **Many samples are infeasible or disconnected.** Only 3–44% of samples are feasible on the 30–60 bus cases, so most of the evaluation budget is spent on islanded or failed configurations.
6. **Memory is the scaling bottleneck.**
   - Peak memory reaches 13.5 GB at 57 buses and 15.5 GB at 60 buses, dominated by the `adj16`/`adj16_pairs` SDPs.
   - The 73–118 bus cases should run `mixed` only, or reduce clique size.

## Next steps

- Run 73-ieee-rts, 89-pegase and 118-ieee, where Juniper took 10–52 min and O-DC-OTS is notably worse than AC-OTS (73: 413133 vs 385194). These are the cases where certificates and alternative solutions are most informative.
- Try `mixed` alone on 179–1354 buses to find the scaling limit.
- Investigate the accuracy loss: tighter MOSEK tolerances, fewer redundant localizing constraints, and whether big-M plus exact switching constraints together hurt conditioning.
