# Project blueprint — Moment/SOS relaxations for power-system MINLPs

Written 2026-09-19 for hand-off to a new Claude account. It records the project's instructions, its
reference material, everything decided and built so far, and what to do next.

---

## 1. Project name

**Moment-SoS-MINLP** — moment/sum-of-squares (Lasserre) relaxations for power-system mixed-integer
nonlinear programs: AC optimal transmission switching (AC-OTS) and AC unit commitment (AC-UC).

- Local path: `/Users/dmolzahn6/Documents/Dan's Files/Moment-SoS-MINLP`
- GitHub: https://github.com/molzahn/Moment-SoS-MINLP
- Owner: Daniel Molzahn (dmolzahn6@gatech.edu), Georgia Tech
- Collaborator: Sam Chevalier (schevali@uvm.edu), who maintains the Overleaf notes repo
  `SamChevalier/Moment-SoS-MINLP-Notes`, which syncs one-way into `notes/notes.tex`

---

## 2. Custom instructions (verbatim)

This project ran in Claude Code, not in a claude.ai Project, so there is no "custom instructions"
text box. There is **no `CLAUDE.md` and no `.claude/` directory in the repository**. The standing
instructions came from the user in conversation and from a saved plan file. Both are reproduced
verbatim below; treat them as the project's instructions.

### 2a. Standing working rules given by the user in conversation

- Commit and push **only when explicitly asked**. Every push so far was individually requested.
- PACE (the Georgia Tech HPC cluster) requires Duo two-factor authentication. **Never enter
  passwords.** When the SSH ControlMaster socket expires, ask the user to run, in their own terminal:
  `ssh -o ServerAliveInterval=60 phoenix`
- Record significant findings, **positive and negative**, in `results/*.md` and in the shared
  Overleaf document.

### 2b. Saved plan file (verbatim)

Stored at `~/.claude/plans/many-important-optimization-problems-vast-blanket.md`. This is the
research plan the whole project follows.

> # Moment/SOS relaxations for AC-OTS / AC-UC: conceptual perspective and next steps
>
> ## Context
> The goal is to extend moment/SOS relaxations, which work well for AC-OPF (see `docs/`), to MINLPs with switching and commitment binaries (AC-OTS, AC-UC). There are two research threads: (1) choosing a mix of relaxation-strengthening techniques that is tractable and still tight, possibly with ML; and (2) using the higher-order pseudo-moments of the binaries for correlated randomized rounding, then recovering the continuous variables with Ipopt. This is a discussion document. No code or repo changes yet.
>
> ---
>
> ## 0. Formulation issues that shape both ideas
>
> - **Binaries in the moment framework.** Encode z ∈ {0,1} as z² = z. **Substitute the reduction z^k → z directly in the monomial basis** instead of adding equality constraints. Kept as equalities, they make the SDP feasible set lower-dimensional, so Slater's condition fails and interior-point solvers run into numerical trouble. This matters in practice at order 2.
> - **Degree mismatch in OTS.** Flow equations become z_ij · (quadratic in V), which has degree 3. The same happens to z-gated angle-difference and thermal limits. So the order-1 (Shor) relaxation cannot represent AC-OTS directly. Two options:
>   (a) Go straight to order ≥ 2. This is the "pure" Lasserre approach, but expensive.
>   (b) Reformulate with on/off lifted line variables (W_ij^z = z_ij W_ij plus McCormick/big-M, or per-line voltage copies as in Kocuk–Dey–Sun), then apply moments.
>   Choosing between (a) and (b) is an early design decision, and it affects what "correlations" the relaxation can even express.
> - **UC is easier structurally.** u_g P_min ≤ P_g ≤ u_g P_max is degree 2, so order 1 already couples u with P and Q. Multi-period constraints (min up/down, ramping) are linear. Suggest starting single-period (or 2–3 periods) so the combinatorics of interest are commitment vs. network, not time coupling.
> - **Sparsity is mandatory.** z_ij involves only buses i and j, so it fits naturally into the chordal cliques from the sparse moment papers. Adding binaries enlarges cliques, though. **Term sparsity** (TSSOS / CS-TSSOS, Wang–Magron–Lasserre) looks well suited, because z·V products are very sparse. This is worth testing early.
> - **OTS-specific pitfalls.** Islanding (angle reference per island, or connectivity constraints). The rectangular-coordinate rotational symmetry per island. Many near-equivalent switching configurations (flat landscape), which matters for rounding (below).
>
> ## 1. Relaxation selection (ML-guided or otherwise)
>
> **Perspective:** the idea is sound. The main risk is that the ML becomes the story before strong non-ML baselines exist. Recommendations:
>
> - **Establish baselines first.** Your earlier sparse-moment work already uses an iterative, mismatch-based heuristic for which buses to elevate to order 2. Extend it with dual-/rank-based signals: SDP dual multipliers, and nonzero eigenvalues of the clique moment matrices beyond the first, which localize where the relaxation fails to be tight. ML only earns its place if it beats these baselines, measured as gap closed per second.
> - **The action space mixes discrete choices with very different costs:**
>   - per-clique order elevation
>   - OBBT on selected variables (|V_i|², angle differences, flows)
>   - QC envelopes (trig/product) added to the moment relaxation
>   - possibly cheaper SOCP/DD/SDD surrogates for parts of the moment matrix (your PowerTech 2015 "mixing SDP and SOCP" idea)
>
>   This is naturally a sequential decision problem (imitation learning or RL, like learning-to-cut: Tang–Agrawal–Faenza; Baltean-Lugojan et al. score SDP cuts for QCQP with neural nets). The practical bottleneck is **label cost**: every label is an SDP solve. Prefer learning from cheap features (mismatch, bound widths, local eigenvalue ratios) and use many small or perturbed instances.
> - **Synergies between OBBT and moments.** Tighter bounds do more than add linear cuts. They enable redundant localizing products like (|V_i|² − l)(u − |V_i|²) ≥ 0 and z-gated bound products at order 2, and those products strengthen the moment relaxation substantially. So OBBT and order elevation are complements, not alternatives.
> - **Key reframing for MINLPs:** for thread 2, the right quality metric is **not only the objective gap**. It is the quality of the rounded solutions produced from the relaxation (feasible-solution cost, success rate). A relaxation with a small bound gap may still spread its pseudo-distribution uselessly. So the selection policy should be evaluated, and perhaps trained, against the rounding outcome. This ties the two threads into one pipeline.
>
> ## 2. Correlated randomized rounding from pseudo-moments
>
> **Perspective:** this is the more distinctive contribution. The theory supports it more strongly than the proposal text states:
>
> - **Local consistency (Rothvoß notes, §2.4).** For binaries, an order-t Lasserre solution restricted to any set S of ≲ t binaries is a genuine probability distribution on {0,1}^S. Order 2 therefore gives exact joint distributions over small groups of switching or commitment decisions. That is exactly the "correlations" in the proposal, now with a rigorous meaning.
> - **Refinement:** even order 1 carries pairwise correlations E[z_i z_j] (and E[z_i V_k]). The marginal-vs-correlated distinction starts at order 1. Higher orders add higher-order, consistent joint distributions. The experiments should separate these effects.
> - **Candidate rounding schemes, from simplest to most principled:**
>   1. Independent rounding of the marginals y_i (baseline).
>   2. **Gaussian/hyperplane rounding:** sample from N(mean, covariance) built from the order-1 moments, then threshold. This is Goemans–Williamson style and standard SDR randomization (Luo et al.). Also in Barak–Kelner–Steurer, "Rounding SoS relaxations".
>   3. **Sequential conditioning** (Karlin–Mathieu–Nguyen decomposition theorem; global correlation rounding of Barak–Raghavendra–Steurer and Guruswami–Sinop, Rothvoß §2.6, §3.5):
>      - Sample z_i from its marginal.
>      - Condition the moments: y'_α = y_{α+e_i}/y_i, or (y_α − y_{α+e_i})/(1 − y_i). This yields a valid order-(t−1) moment vector, including for the continuous variables.
>      - Repeat, or **re-solve the relaxation with z_i fixed** (randomized diving).
>      - Bonus: the conditional moments E[V | z] give warm starts for Ipopt.
>   4. Fit a max-entropy/Ising model to the pairwise or 4-wise binary moments and sample it with Gibbs sampling.
>   5. Approximate atom extraction (Henrion–Lasserre / Christoffel–Darboux) on the binary block when the moment matrix is nearly flat.
> - **Where correlations should visibly help:** symmetric or degenerate structure such as parallel lines or identical units. Marginals of 0.5 are uninformative, but E[z_i z_j] ≈ 0 reveals "exactly one on". Build small illustrative examples around this, in the spirit of your "illustrative example" paper.
> - **Practical issues:**
>   - Rounded OTS configurations may island the network, and UC samples may lack capacity. Use cheap pre-screens (connectivity, capacity ≥ load + reserve) before Ipopt.
>   - Ipopt failure is not a proof of infeasibility.
>   - Deduplicate samples.
>   - The number of samples needed is a result to report.
> - **Relation to prior work (initial scan; needs a careful review):**
>   - Moment-SOS plus rounding for MIPOPs exists as *granularity/feasible rounding* (Neumann/Stein et al., JOTA 2025). That is deterministic rounding of inner-parallel-set solutions, **not** correlation-based sampling.
>   - The TCS literature (SoS rounding, global correlation rounding) is mature but rarely applied to physics-constrained MINLPs.
>   - Power systems examples: a low-rank rounding heuristic for SDR of hydro UC (arXiv 1711.06194); relax-and-round for UC-ACOPF (arXiv 2501.11355, based on the continuous relaxation); SDP-based Benders for UC (arXiv 1903.02628); Fattahi–Lavaei–Atamtürk OTS bound strengthening; Kocuk et al. MISOCP OTS.
>   - I did not find work that samples from **Lasserre pseudo-distributions over switching/commitment binaries** for AC problems. That looks like a genuine gap, pending a thorough search.
>
> ## 3. Proposed next steps (for later execution, after discussion)
>
> 1. **Literature review.** Record it in the repo as `docs/literature_notes.md` and summarize significant items in the Overleaf file. Topics: SoS rounding/conditioning; MINLP rounding heuristics (RENS, diving); SDP/moment work on OTS and UC; ML for relaxation and cut selection; OBBT for OPF (Coffrin–Hijazi–Van Hentenryck; Chen–Atamtürk–Oren; Sundar et al.).
> 2. **Formulation note.** Write out AC-OPF, single-period AC-UC, and AC-OTS as polynomial programs in rectangular coordinates, with the degree, sparsity, and z-reduction choices above and the OTS option (a) vs (b).
> 3. **Toolchain (decided): Julia + MOSEK.** Use JuMP, PowerModels.jl for PGLib data parsing, and TSSOS.jl for the sparse moment-SOS construction, with custom JuMP-based moment code where finer control is needed (z-reduction, per-clique orders, conditioning). MosekTools.jl is the SDP solver and Ipopt.jl handles NLP recovery.
>    - **First action:** set up a Julia project environment in the repo (`Project.toml`). Then verify the MOSEK license with a tiny JuMP SDP. MOSEK is believed installed and licensed, but this has not been tested yet.
> 4. **Small, enumerable test cases** (for example, 5–30 bus PGLib cases with a handful of switchable lines or units) so global optima are known by enumeration + Ipopt/BARON. Only then scale up.
> 5. **First experiment:** order-1 vs order-2 relaxations × rounding schemes 1–3, reporting bound gap, best feasible cost, and success rate vs. number of samples.
>
> ## Verification (for the eventual experiments)
> - Enumeration gives ground-truth optima on small instances.
> - Check local consistency numerically: the moments on small binary subsets form valid distributions.
> - Confirm the conditioned moment matrices stay PSD.
> - Compare rounding results against independent rounding and against a MISOCP/QC-based heuristic baseline.

### 2c. Saved memory (verbatim)

`~/.claude/projects/-Users-dmolzahn6-Documents-Dan-s-Files-Moment-SoS-MINLP/memory/toolchain-julia-mosek.md`:

> Moment/SOS MINLP project (AC-OPF/OTS/UC relaxations + correlated rounding) uses a Julia toolchain with MOSEK as the SDP solver (Ipopt for NLP recovery).
>
> **Why:** User chose Julia + MOSEK explicitly on 2026-09-12. MOSEK is believed installed and licensed, but that hadn't been tested yet.
> **How to apply:** Build code in Julia/JuMP. Verify the MOSEK license with a tiny SDP before relying on it. Collaboration happens through a shared GitHub repo and an Overleaf file that tracks positive and negative findings.

(The license has since been verified and used on both the Mac and PACE.)

---

## 3. Knowledge files

Nothing was uploaded to a claude.ai project. The equivalent reference material lives in `docs/` in
the repository.

### Papers (PDF)

| file | description |
|---|---|
| `docs/taheri_molzahn-optimal_dcots.pdf` | **The benchmark paper.** Taheri & Molzahn, optimal DC-OTS. Defines the 18 PGLib "api" test cases used in experiments 3 and 4 and reports reference values (AC-OTS, and O-DC-OTS = a DC-optimal switching solution evaluated in AC). |
| `docs/Certified and Accurate SDP Bounds for ACOPF Problem.pdf` | **Oustry, D'Ambrosio, Liberti & Ruiz**, *Certified and accurate SDP bounds for the ACOPF problem*, PSCC 2022. (Not Josz — an earlier version of this table misattributed it.) Certified bounds from inexact SDP solves for AC-OPF. Our certification work extends this to per-clique higher orders and binaries. |
| `docs/molzahn_hiskens-sparse_moment_opf.pdf` | Sparse moment relaxations for OPF; source of the per-clique order-elevation heuristic. |
| `docs/molzahn_hiskens-Moment_OPF.pdf` | Original moment relaxation of OPF. |
| `docs/molzahn_hiskens-fnt2019.pdf` | Survey monograph on OPF relaxations. |
| `docs/molzahn_hiskens-illustrative_example.pdf` | Small examples where low-order relaxations fail; template for the "correlations matter" illustrative examples. |
| `docs/molzahn_hiskens-powertech2015.pdf` | Mixing SDP and SOCP per-bus; the cheaper-surrogate idea in the plan. |
| `docs/molzahn_josz_hiskens_panciatici-cdc2015.pdf`, `docs/molzahn_josz_hiskens_pantiatici-pscc2016.pdf` | Laplacian-based and mixed SDP/SOCP approaches for large systems. |
| `docs/lasserre survey.pdf` | Lasserre's moment-SOS survey; background theory. |

### Notes written in this project

| file | description |
|---|---|
| `docs/formulations.md` | AC-OPF, AC-OTS and AC-UC as polynomial programs in rectangular coordinates, exactly as implemented in `src/power.jl`, including the z² = z reduction and the degree-3 switching constraints. |
| `docs/literature_notes.md` | Literature review: SoS rounding and conditioning, MINLP rounding heuristics, SDP/moment work on OTS and UC, ML for cut/relaxation selection, OBBT for OPF. |
| `notes/notes.tex` | The shared Overleaf document (one-way sync from Sam's Overleaf repo into `notes/`). Contains a section "Status as of September 14, 2026". **Anything written here can be overwritten by the next Overleaf sync**; coordinate with Sam before editing. |
| `results/*.md` | Per-experiment findings (see below). |
| `data/` | MATPOWER cases 5/9/14/24/30 plus `data/pglib_api/` for the DC-OTS paper cases. All data modifications are documented in `scripts/instances.jl`. |

---

## 4. What has been done

### 4.1 Research threads

1. **Relaxation strengthening selection** — which cliques get order 2, which extra blocks and cuts
   to add, and how to certify the resulting bound.
2. **Correlated randomized rounding** — sampling binaries from the relaxation's pseudo-moments, then
   recovering the continuous variables with Ipopt through PowerModels.

### 4.2 Key modeling decisions (and why)

- **Binaries:** z ∈ {0,1} is handled by substituting z^k → z in the monomial basis, not by adding
  z² = z equality constraints. Equalities destroy Slater's condition and wreck the interior-point
  solve at order 2.
- **Switching constraints are degree 3** (p = z · P(V)). Two formulations are implemented and both
  have the same feasible set over binary z:
  - **big-M** (`bigM_switching=true`, the default): adds degree-2 big-M, on/off and angle big-M
    constraints;
  - **no-big-M** (`bigM_switching=false`): keeps only the exact degree-3 constraints.
  A third option, `scalar_tags`, enforces chosen constraints only as L(g) ≥ 0 with no localizing
  matrix, which is how the big-M constraints can be kept at order 1 in an order-2 model.
- **Order policy, not a global order.** Only cliques containing a binary get order 2; the rest stay
  at order 1 (`binary_clique_order`). That is what "mixed" means, and it is what makes the degree-3
  constraints affordable — roughly 40–48% of cliques are order 2 and the largest PSD block is 54–75.
  Order 1 everywhere cannot represent the degree-3 constraints and the bound collapses (73-bus
  78116, 118-bus 0, 200-bus 20079).
- **Relaxation variants compared:** `mixed`, `mixed_nobigM`, `adj16` (adjacent-clique merging with
  maximum size 16), `adj16_pairs` (plus extra order-2 blocks on binary pairs).
- **Low-impedance bus merging** (`merge_zmax`, default 1e-3, non-transformer branches only). Merging
  is done **inside the POP**, not in the data: the merged buses share one (E,F) voltage pair, but each
  bus keeps its own power balance and the tie lines keep their thermal limits as lossless flows and
  are not switchable. A plain data-level merge drops binding tie limits and changes the answer
  (89-pegase AC-OPF: 125456 vs. the correct 130175). Merging fixed an SDP stall on 89-pegase that was
  caused by lines with |z| ≈ 2.2e-4.
- **Islanding.** PowerModels' AC-OTS allows islands, and the benchmark paper's own 89- and 118-bus
  topologies contain islands, so requiring connectivity is wrong. The screen rejects only islands
  that provably cannot balance active power: Σpmax < Σ(pd + gs·vmin²), which is valid because losses
  are nonnegative. Ipopt "locally infeasible" is not a proof of infeasibility, and neither is DC
  infeasibility.
- **Ipopt uses HSL MA97 by default** through a licensed local `HSL_jll`, falling back to MUMPS with a
  warning. MA97 and MUMPS were verified to agree on objectives (118, 300, 1354) and on single-line
  feasibility (118). ⚠️ `Manifest.toml` pins a **Mac-specific local path** to the licensed HSL_jll;
  a collaborator needs their own licensed build or the registry HSL_jll (which falls back to MUMPS).
- **Certified bounds from inexact SDP solves.** MOSEK stops early on these problems
  (`SLOW_PROGRESS`), so the raw dual objective is not a valid bound. Four certificates are
  implemented in `src/certify.jl`:
  - **box**: t − Σ|r_α|·max|x^α| − Σ ρ_k λ_min⁻ (the residual term dominates);
  - **box with tight trace bounds** from quadratic group constraints such as e² + f² ≤ V̄²
    (in practice this never changed the answer);
  - **implicit Gram**: residuals absorbed into representative Gram entries of the moment blocks;
  - **per-block hybrid** of the two, and an **optimized** certificate: a smoothed soft-min surrogate
    maximized with L-BFGS from the MOSEK dual. A proximal bundle method was also implemented but is
    weak in practice; `:smooth` is the method to use.

### 4.3 Code

Julia package `MomentMINLP` (Julia 1.12.7 at `~/julia-1.12.7/bin`, not on PATH; MOSEK, Ipopt+HSL).

| file | contents |
|---|---|
| `src/MomentMINLP.jl` | module; `__init__` finds the MOSEK license, sets up HSL, silences PowerModels |
| `src/polynomials.jl`, `src/pop.jl` | sparse polynomials, POP container, Ipopt solve of a POP |
| `src/sparsity.jl` | chordal extension and maximal cliques |
| `src/relaxation.jl` | the sparse moment relaxation: per-clique orders, binary reduction, `out_of_graph_tags`, `scalar_tags`, `extra_cliques`, `certify=:implicit/:optimize/:bundle` |
| `src/certify.jl` | `SOSCertificate`, `certified_value`, `smooth_certify` (L-BFGS), `bundle_certify` |
| `src/power.jl` | `build_power_pop` (switching/commitment, big-M toggle, merging, radial fixing, connectivity cuts), `merge_low_impedance`, `ConfigEvaluator`, `screen_config` |
| `src/islands.jl` | `island_screen`, `IslandGuard`, `repair_islands!`, `can_open`, `radial_fixings`, `connectivity_cuts` |
| `src/rounding.jl` | independent, Gaussian, conditional (with guard), and diving rounding |
| `scripts/experiment1.jl` … `experiment4_islanding.jl`, `recertify.jl` | the experiments (see below) |
| `scripts/test_islands.jl` | island unit tests — **all pass**; run these after touching `src/islands.jl` |
| `hpc/` | PACE SLURM scripts: `pace_env.sh`, `setup_pace.sh`, `experiment3.sbatch`, `experiment4.sbatch`, `recertify.sbatch`, and the `submit_*.sh` wrappers with per-case memory/time maps |

### 4.4 Experiments and results

- **Experiment 1** (`results/experiment1_findings.md`): order 1 vs. order 2 × rounding schemes on
  small enumerable cases, with ground truth by enumeration.
- **Experiment 2** (`results/experiment2_findings.md`): scaling on cases 14/24/30, mixed-order and
  augmented sparse relaxations vs. NLP relax-and-round.
- **Numerics study** (`results/numerics_findings.md`): normalization and variable scaling.
- **Experiment 3** (`results/experiment3_findings.md`, `results/experiment3_tables.md`): the 18
  DC-OTS paper cases with every branch switchable, all four relaxation variants.
- **Experiment 4** (`results/experiment4_*.json`): islanding-aware rounding — evaluation that allows
  islands, repair, guarded conditional sampling, and relaxations with radial lines fixed and
  connectivity cuts.
- **Re-certification** (`results/recertify_findings.md`): every saved relaxation re-solved and
  certified with all certificate variants.

**Certified gaps after re-certification** (best bound over all variants vs. the best known value,
which is the minimum of ours, the paper's AC-OTS and the paper's O-DC-OTS):

| case | gap before re-certification | gap after |
|---|---|---|
| 5-pjm | 2.42% | 2.41% |
| 14-ieee | 4.49% | 4.39% |
| 24-ieee-rts | 6.86% | 6.68% |
| 30-as | 1.31% | 1.06% |
| 30-ieee | 4.22% | 3.62% |
| 39-epri | 0.51% | 0.48% |
| 57-ieee | 0.17% | 0.10% |
| 60-c | 1.08% | 0.65% |
| 73-ieee-rts | 4.44% | 4.38% |
| 89-pegase (merged model) | 0.71% | 0.45% |
| 118-ieee | 6.08% | 5.98% |
| 179-goc | 6.28% | 4.88% |
| 200-activ | 0.30% | 0.09% |
| 300-ieee | 0.92% | 0.80% |

The certified bound is now within 0.25% of the raw SDP objective on every relaxation, and within
0.03% on most: the remaining gaps come from the relaxations, not from the inexact solves.

**Best solutions found**, with notable wins: 39-epri 246561.20 (better than Juniper, opening lines
{4,6}), 500-goc 685131 (better than the paper's O-DC-OTS 692271), 200-activ 35700.74, 240-pserc
4635816, 300-ieee 684431, 1354-pegase 1498117.

**Findings worth keeping**

- Certification matters as much as the relaxation. The big-M formulation looked worse than no-big-M
  largely because its certificates were loose; certified properly, big-M gives the better bound on
  179-goc (1836683 vs. 1810514). No-big-M still solves 1.5–2× faster.
- The tight-trace-bound refinement never helped, because the residual term dominates the box
  certificate.
- Rounding quality falls off a cliff with size. 89-pegase reaches 38–55% AC-feasible samples and
  240-pserc 24–37%, but 118, 179, 300, 500 and 1354 give 0–2%. At those sizes the good solutions come
  from the one-flip polishing step, not from the rounding itself. **Making rounding work at scale is
  the central open problem.**
- The connectivity cuts and radial fixings help on 89-pegase (bound 99486 → 99627 before
  re-certification, and more feasible samples), but not on 118/179/300.
- These are stressed "api" cases: on 118-ieee, 100 of the 186 single-line openings are AC-infeasible,
  and feasibility is **not monotone** in the set of opened lines.

### 4.5 Infrastructure

- **PACE Phoenix.** Account `gts-dmolzahn6`, scratch directory `~/scratch/Moment-SoS-MINLP`, Julia at
  `~/scratch/parameter_optimized_OTS/julia-1.12.7/bin/julia`, depot on scratch,
  `JULIA_CPU_TARGET="generic;x86-64-v3,clone_all"`, MOSEK license at `~/mosek/mosek.lic`.
  Sync with rsync, excluding `.git`, `docs/*.pdf`, `Manifest.toml` (Mac-specific HSL path),
  `hpc/logs` and results. Duo 2FA is required; when the connection dies, ask the user to run
  `ssh -o ServerAliveInterval=60 phoenix`.
- **Git.** Work is on branch `experiment3-dcots`. Open PR: https://github.com/molzahn/Moment-SoS-MINLP/pull/1
  (`experiment3-dcots` → `main`). Last commit is `0dd8367`; commits `c9539bc`, `cd00c3a` and `0dd8367`
  are local and unpushed, and the re-certification and experiment-4 results below are uncommitted.

---

## 5. Current status (2026-09-19)

**All PACE jobs have finished.** Experiment 4 completed for every case, and all 18 re-certification
jobs completed (57-ieee needed 48G after an out-of-memory failure at 16G).

Uncommitted in the working tree: `results/recertify_*.json`, `results/recertify_findings.md`,
`results/experiment4_{240-pserc,500-goc,1354-pegase}-api.json`, `results/pace_logs/*`, the
`hpc/submit_recertify.sh` memory-map edit, and `docs/Certified and Accurate SDP Bounds for ACOPF Problem.pdf`.

### Results that arrived last and are not yet written up

| case | variant | raw | optimized certificate | comment |
|---|---|---:|---:|---|
| 57-ieee | adj16_pairs | 49229.35 | 49220.93 | completes the 57-bus row |
| 240-pserc | mixed (big-M) | 4607420 | 4605168 | big-M now beats no-big-M here too |
| 240-pserc | mixed_nobigM | 4607220 | 4606997 | best 240-bus bound; gap ≈ 0.49% |
| 500-goc | mixed_nobigM | 667245 | 667141 | **first bound for 500-goc**; vs. best found 685131 → gap ≈ 2.6% |
| 1354-pegase | all variants | 2.58M–3.77M | 2.58M–3.77M | ⚠️ **invalid, see below** |

### ⚠️ Open bug: the 1354-pegase "bounds" exceed a known feasible solution

The re-certification reports 2580996 (mixed), 2612454 (no big-M) and 3772084 (conn) for
1354-pegase, but a feasible solution of cost **1498117** is known. A lower bound above a feasible
cost is impossible, so something is wrong and **these numbers must not be reported**. Likely
suspects, in order:

1. MOSEK is far from converged at this size (`SLOW_PROGRESS`), so the dual objective is meaningless
   and the certificate machinery may be operating on a garbage Gram matrix — check
   `gram_min_eig` and `residual_correction` in the saved info.
2. A modeling difference in the 1354 POP itself: low-impedance merging, radial fixings or
   connectivity cuts may be cutting off the incumbent. Evaluate the incumbent's switching
   configuration inside the POP used for the bound and check that it is feasible and its cost
   matches.
3. A scaling or normalization overflow specific to this case.

A good regression test: the certified bound must never exceed the best known feasible cost. Add that
assertion to `scripts/recertify.jl` so this cannot slip through again.

---

## 6. Suggested next steps

1. **Diagnose the 1354-pegase bound** (above). Until it is resolved, report "no valid bound" for that
   case. Add the bound-vs-incumbent assertion to the scripts.
2. **Finish the write-up.** Update `results/experiment3_tables.md` and `results/experiment3_findings.md`
   with the re-certified gaps, write `results/experiment4_findings.md` (rounding feasibility by case,
   the effect of cuts and repair, and the stressed-case infeasibility observation), and fold both into
   `results/recertify_findings.md`'s tables, which currently cover only the first 14 cases.
3. **Update `notes/notes.tex`** with sections on certification and on experiment 4, then update
   PR #1. Two things the PR description must say: `notes/` will be overwritten by the next Overleaf
   sync, and the `Manifest.toml` HSL_jll path is Mac-specific.
4. **Write the email summary for Sam Chevalier** (schevali@uvm.edu). This was requested earlier and
   has never been delivered. It should cover: the mixed per-clique order policy, low-impedance
   merging, big-M vs. no-big-M, the certification results, and the rounding-at-scale problem.
   **Draft it for the user to send; do not send email.**
5. **Attack the central research problem: rounding at scale.** Ideas in rough priority order:
   - sample fewer lines at once (limit the number of simultaneous openings to match E[#opened]);
   - condition and re-solve (randomized diving) rather than sampling all binaries at once;
   - use the conditional moments E[V | z] as Ipopt warm starts;
   - fit a max-entropy/Ising model to the pairwise moments and sample it with Gibbs sampling;
   - use a local AC feasibility screen stronger than the power-balance island screen, since the api
     cases are so tightly constrained.
6. **The relaxation-selection thread** (plan §1) is still mostly untouched: OBBT, QC envelopes and a
   dual-/eigenvalue-guided order-elevation policy, with the ML layer only after those baselines exist.

## 7. Working preferences to carry over

- Commit and push only when asked. Never enter passwords; for PACE, ask the user to re-establish the
  SSH connection themselves.
- Record negative results as carefully as positive ones, in `results/*.md` and in Overleaf.
- Prefer running long jobs on PACE, one SLURM job per case, with per-case memory and time settings.
- Be precise about what is and is not a proof: an inexact SDP objective is not a bound, an Ipopt
  failure is not an infeasibility certificate, and a DC infeasibility does not imply AC infeasibility.
- Julia is at `~/julia-1.12.7/bin/julia` and is not on PATH; run scripts as
  `~/julia-1.12.7/bin/julia --project=. scripts/<script>.jl <case>`.
