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

## Update: 73-ieee-rts, 89-pegase, 118-ieee

*Logs: `results/experiment3_log_large.txt` (73 and 118; 89 before merging), `results/experiment3_log_89merged.txt`.*

These runs use only the cheap relaxations. The adj16 relaxation was also run on 89-pegase, where it needed 21 GB of memory.
- **`mixed`:** as above.
- **`mixed_nobigM`:** the same order policy on the POP without the degree-2 big-M switching constraints (same variables).

| case | paper O-DC-OTS | paper AC-OTS (Juniper time) | best found (ours) | vs paper AC-OTS | best bound (variant) | certified gap |
|---|---|---|---|---|---|---|
| 73-ieee-rts | 413133 | 385194 (611 s) | 389141.92 | +1.02% | 368109.94 (mixed_nobigM) | 4.44% |
| 89-pegase (merged ties) | 100702 | 100344 (1867 s) | 100539.74 | +0.20% | 99517.69 (mixed) | 0.82%† |
| 118-ieee | 195918 | 180312 (3098 s) | 230418.58 | +27.8% | 169350.89 (mixed_nobigM) | 6.08% |

† Bound of the merged model; see the caveat below.

### 89-pegase: big-M failure and low-impedance merging

- **Without merging, `mixed` fails.**
  - MOSEK stalls after 3 iterations (PFEAS ≈ 2.3e3), reports `SLOW_PROGRESS` with no feasible SOS point, and gives no bound (NaN).
  - Its "marginals" are noise: all 210 are below 0.5, and every sample is disconnected.
  - Turning off variable scaling or constraint normalization does not help. Order 1 solves (bound 9138, useless).
  - Dropping big-M solves (bound 99445.94).
- **Cause:** 19 non-transformer branches with |z| ≈ 2.2e-4, i.e. admittances up to 4.5e3, next to thermal limits as small as 0.04 p.u. The big-M constants M = |y|·Vmax² then multiply terms that must cancel.
- **Plain merge (`merge_low_impedance`):** collapse the buses as in the LCOTS project (`LCOPF.jl`, `merge_zero_impedance`, |z| < 1e-3).
  - This **drops the ties' thermal limits, and three of them bind**: ties 73, 96, 160 at 3.58/3.58, 3.28/3.32, 5.45/5.66 p.u.
  - AC-OPF on the merged data is 125456 vs 130175 on the original (−3.6%). Removing only the tie limits on the original network gives 126100.
- **Voltage-level merge (`build_power_pop(...; merge_zmax = 1e-3)`, used here):**
  - Buses joined by a tie share one voltage (E, F); every bus keeps its own power balance.
  - Each tie is a lossless flow (p, q) with its thermal limit and is not switchable. This is the zero-impedance limit, as for PowerModels switches.
  - The POP shrinks from 1251 to 1156 variables and 210 to 191 binaries. With all lines closed, its NLP gives 129734.63 vs 130174.85 on the original (−0.34%).
  - Configurations are still evaluated with PowerModels on the original network, with the ties closed.
  - Merging changes nothing else up to 118 buses. Larger paper cases with ties: 179-goc (2), 240-pserc (55), 300-ieee (2 non-transformer), 1354-pegase (184).
- **Effect of merging on 89-pegase:**
  - `mixed` now solves (bound 99517.69, 70 s).
  - 34–52% of rounded samples are feasible (0% before).
  - Rounding finds 100545.72 (Gaussian from `mixed_nobigM`), polished to 100539.74. That is +0.20% vs Juniper and beats the paper's O-DC-OTS (100702).
- **Caveat on the 0.82% gap.** The bound is for the merged model: ties are zero-impedance and cannot be switched, so it is not a rigorous bound for the original AC-OTS.
  - The paper's Juniper topology opens 4 ties (99, 129, 138, 155); the 56 opened lines were parsed from Table I and reproduce 100343.59.
  - With those ties closed instead, the same topology costs **100321.78** on the original network, slightly better than the paper. On the merged model it costs 100312.69.
  - So keeping ties closed loses nothing here, but that is evidence, not proof.

### Big-M on the larger cases

Dropping the big-M constraints gives **higher certified bounds and 2–3× faster solves** on 73-ieee-rts (368110 vs 367152; 19 s vs 53 s) and 118-ieee (169351 vs 168242; 42 s vs 87 s).
- The raw objectives are nearly identical (73: 368244 vs 368346; 118: 169536 vs 169580).
- The difference is the inexact-solve correction: 73: 134 vs 1193; 118: 185 vs 1338. Big-M hurts conditioning more than it tightens the relaxation.
- On merged 89-pegase the order reverses (99518 with big-M vs 99293 without).
- This qualifies the earlier numerics finding that dropping big-M weakens bounds, which came from the small instances.

### Rounding at 73–118 buses

- **73-ieee-rts:** best rounded +3.7% (Gaussian); 1-flip polish reaches +1.02%. The polish hit its 40-evaluation budget while still improving. The paper's solution opens 17 lines; ours opens 10.
- **118-ieee: rounding fails.**
  - Only 2–4 marginals are below 0.5, and 0–1% of the 100 samples per scheme are feasible (almost all disconnect the network or fail in Ipopt).
  - The best configuration (+27.8%) comes from 1-flip polish starting at the all-closed topology.
  - The paper's solution opens 33 lines; the relaxation marginals give no sign of that.
- **89-pegase (merged):** Gaussian/adj16 rounding does best (+0.20–0.22%), with independent at +0.36–0.38%. As in the small cases, correlated rounding has no systematic edge.
- **Memory:** 3.3–5.8 GB for `mixed` at 73–118 buses; 21 GB peak on 89-pegase with `adj16`, whose certified bound is also the weakest (97009, correction 2992).

## Next steps

- Run 73-ieee-rts, 89-pegase and 118-ieee, where Juniper took 10–52 min and O-DC-OTS is notably worse than AC-OTS (73: 413133 vs 385194). These are the cases where certificates and alternative solutions are most informative.
- Try `mixed` alone on 179–1354 buses to find the scaling limit.
- Investigate the accuracy loss: tighter MOSEK tolerances, fewer redundant localizing constraints, and whether big-M plus exact switching constraints together hurt conditioning.
