# Literature notes: moment/SOS relaxations and rounding for AC-OTS / AC-UC

*Started 2026-09-13. Bibliographic details were checked by web search where marked ✓. Everything else is from memory and should be verified before citing. "Relevance" notes are our own assessment.*

Contents
1. [Moment/SOS hierarchy: foundations and sparsity](#1-momentsos-hierarchy-foundations-and-sparsity)
2. [Moment/SOS relaxations of AC-OPF](#2-momentsos-relaxations-of-ac-opf)
3. [Rounding from Lasserre / SoS solutions (TCS)](#3-rounding-from-lasserre--sos-solutions-theoretical-cs)
4. [Rounding and solution extraction in continuous/MINLP optimization](#4-rounding-and-solution-extraction-in-continuous--minlp-optimization)
5. [AC optimal transmission switching](#5-ac-optimal-transmission-switching)
6. [Unit commitment with AC power flow](#6-unit-commitment-with-ac-power-flow)
7. [Bound tightening, QC, and strengthening relaxations](#7-bound-tightening-qc-and-strengthening-relaxations)
8. [Machine learning for relaxation / cut / bound selection](#8-machine-learning-for-relaxation--cut--bound-selection)
9. [Takeaways and open gaps](#9-takeaways-and-open-gaps)

---

## 1. Moment/SOS hierarchy: foundations and sparsity

- **Lasserre (2001)**, *Global optimization with polynomials and the problem of moments*, SIAM J. Optim. 11(3). The original hierarchy.
- **Lasserre (2001)**, *An explicit exact SDP relaxation for nonlinear 0-1 programs*, IPCO. For n binary variables, the hierarchy is exact at order n. This justifies the reduction x² = x in the monomial basis.
- **Laurent (2003)**, *A comparison of the Sherali–Adams, Lovász–Schrijver, and Lasserre relaxations for 0-1 programming*, Math. Oper. Res. 28(3).
  - Lasserre dominates the other two.
  - For binaries, order t implies Sherali–Adams-type consistency on small subsets (§3).
- **Rothvoß (2013)**, *The Lasserre hierarchy in approximation algorithms*, MAPSP tutorial notes (`docs/lasserre survey.pdf`). The most useful single reference for thread 2:
  - §2.4: locally consistent probability distributions.
  - §2.6: decomposition theorem.
  - §3.5: global correlation rounding.
- **Waki, Kim, Kojima, Muramatsu (2006)**, SIAM J. Optim. Correlative sparsity via chordal extension, which is what our code uses.
- ✓ **Josz & Molzahn (2018)**, *Lasserre hierarchy for large scale polynomial optimization in real and complex variables*, SIAM J. Optim. 28(2):1017–1048 ([arXiv:1709.04376](https://arxiv.org/abs/1709.04376)).
  - Complex hierarchy and multi-ordered hierarchy (per-clique orders, with global convergence).
  - The multi-ordered idea is the formal backbone of thread 1's order selection.
- **Wang, Magron, Lasserre (2021)**, *TSSOS: a moment-SOS hierarchy that exploits term sparsity*, SIAM J. Optim. **CS-TSSOS** (Wang, Magron, Lasserre, Mai; ACM TOMS 2022) combines correlative and term sparsity.
- **Lasserre, Pauwels, Putinar**, *The Christoffel–Darboux Kernel for Data Analysis* (Cambridge Univ. Press, 2022).
  - ✓ **Lasserre**, *The Moment-SOS hierarchy and the Christoffel–Darboux kernel*, Optim. Lett. (2021) ([arXiv:2011.08566](https://arxiv.org/abs/2011.08566)).
  - ✓ *Leveraging Christoffel–Darboux kernels to strengthen moment-SOS relaxations* ([arXiv:2501.14281](https://arxiv.org/abs/2501.14281)). Uses CD kernels from low-order relaxations to tighten bounds and help extract minimizers. This is a possible alternative to rounding for extracting candidate points.

## 2. Moment/SOS relaxations of AC-OPF

Papers in `docs/`:
- **Molzahn & Hiskens**, *Moment-based relaxation of the optimal power flow problem*, PSCC 2014.
- **Molzahn & Hiskens**, *Sparsity-exploiting moment-based relaxations of the optimal power flow problem*, IEEE TPWRS 2015. Iterative, mismatch-based selection of buses whose cliques get order 2. This is the baseline to beat in thread 1.
- **Molzahn & Hiskens**, *Mixed SDP/SOCP moment relaxations of the optimal power flow problem*, PowerTech 2015. Cheaper cones for parts of the higher-order constraints.
- **Molzahn, Josz, Hiskens, Panciatici**, *Solution of OPF problems using moment relaxations augmented with objective function penalization*, CDC 2015.
- **Molzahn, Josz, Hiskens, Panciatici**, *Computational analysis of sparsity-exploiting moment relaxations of the OPF problem*, PSCC 2016.
- **Molzahn & Hiskens**, *Convex relaxations of optimal power flow problems: an illustrative example*, IEEE TCAS-I 2016.
- **Molzahn & Hiskens**, *A survey of relaxations and approximations of the power flow equations*, Foundations and Trends in Electric Energy Systems, 2019.

Other:
- ✓ **Ghaddar, Marecek, Mevissen**, *Optimal power flow as a polynomial optimization problem*, IEEE TPWRS 31(1):539–546, 2016 ([arXiv:1404.3626](https://arxiv.org/abs/1404.3626)). Sparse Lasserre hierarchy plus valid inequalities, and global solutions up to 39 buses.
- ✓ **Wang, Magron, Lasserre (and Mai)**, *Certifying global optimality of AC-OPF solutions via sparse polynomial optimization*, Electric Power Systems Research, 2022 ([arXiv:2109.10005](https://arxiv.org/abs/2109.10005)). CS-TSSOS certifies 1% optimality on PGLib cases up to 24,464 buses. **Implication:** term sparsity is the most promising scalability route for our OTS/UC formulations.
- ✓ **Molzahn, Josz, Hiskens**, *Moment relaxations of optimal power flow problems: beyond the convex hull* ([arXiv:1612.02519](https://arxiv.org/pdf/1612.02519)).

## 3. Rounding from Lasserre / SoS solutions (theoretical CS)

This is the theory behind thread 2. For binary problems, an order-t Lasserre solution defines pseudo-distributions that are **genuine distributions on every subset of ≤ t variables**. Conditioning on a variable's value produces a valid order-(t−1) solution.

- ✓ **Karlin, Mathieu, Nguyen**, *Integrality gaps of linear and semi-definite programming relaxations for Knapsack*, IPCO 2011.
  - Decomposition theorem: a Lasserre solution is a convex combination of solutions that are integral on a chosen small set.
  - The Knapsack integrality gap drops to t/(t−1) after t rounds.
- ✓ **Barak, Raghavendra, Steurer**, *Rounding semidefinite programming hierarchies via global correlation*, FOCS 2011 ([arXiv:1104.4680](https://arxiv.org/abs/1104.4680)). Condition on a few randomly chosen variables to reduce global correlation, then round independently.
- ✓ **Guruswami & Sinop**, *Lasserre hierarchy, higher eigenvalues, and approximation schemes for graph partitioning and quadratic integer programming with PSD objectives*, FOCS 2011 ([arXiv:1104.4746](https://arxiv.org/abs/1104.4746)). Chooses which variables to condition on using column selection / low-rank approximation of the moment matrix.
- ✓ **Barak, Kelner, Steurer**, *Rounding sum-of-squares relaxations*, STOC 2014 ([arXiv:1312.6652](https://arxiv.org/abs/1312.6652)). Uses pseudo-moments as moments of a Gaussian or other distribution for rounding.
- **Goemans & Williamson (1995)**, JACM. Hyperplane rounding: the prototype of correlated rounding from the second-moment matrix.

**Relevance:** these results give approximation guarantees for structured combinatorial problems. None of them treats continuous physics-based constraints (power flow) coupled to the binaries. Our "conditional" and "dive" schemes are direct implementations of conditioning. Our "gaussian" scheme is the Goemans–Williamson / Barak–Kelner–Steurer idea with a copula so the marginals are preserved.

## 4. Rounding and solution extraction in continuous / MINLP optimization

- **Henrion & Lasserre (2005)**, *Detecting global optimality and extracting solutions in GloptiPoly*. Flat-extension extraction of minimizers from exact relaxations.
- **Luo, Ma, So, Ye, Zhang (2010)**, *Semidefinite relaxation of quadratic optimization problems*, IEEE Signal Processing Magazine. Gaussian randomization from the SDR solution, standard in signal processing.
- ✓ **Eggen, Stein, Volkwein**, *Granularity for mixed-integer polynomial optimization problems*, JOTA 205:22, 2025 ([link](https://link.springer.com/article/10.1007/s10957-025-02631-6)).
  - Feasible rounding based on inner parallel sets.
  - The moment-SOS hierarchy is used to solve the continuous problems whose optima are rounded.
  - Tested on MINLPLib.
  - **Closest "moment-SOS + rounding for MIPOP" work found**, but it uses deterministic rounding of a point with a feasibility guarantee. It does **not** sample from pseudo-moment correlations. Power-flow equality constraints make granularity (non-empty inner parallel set) fail, so their approach does not apply directly to AC-OTS/UC.
- **Berthold (2014)**, *RENS: the optimal rounding*, Math. Prog. Comp. Also diving and other relaxation-based primal heuristics in MIP/MINLP solvers. Our "dive" scheme is a randomized diving heuristic driven by an SDP relaxation.
- ✓ **Gally, Pfetsch, Ulbrich**, *A framework for solving mixed-integer semidefinite programs*, Optim. Methods Softw. 33(3):594–632, 2018 (SCIP-SDP).
  - Branch-and-bound over SDP relaxations, with SDP-based rounding heuristics.
  - This is a natural host if we move from heuristics to exact solution of the mixed-integer SDP.

## 5. AC optimal transmission switching

- **Fisher, O'Neill, Ferris (2008)**, *Optimal transmission switching*, IEEE TPWRS 23(3). DC model.
- ✓ **Hijazi, Coffrin, Van Hentenryck**, *Convex quadratic relaxations for mixed-integer nonlinear programs in power systems*, Math. Prog. Comp. 9:321–367, 2017. QC relaxations with on/off constraints, including OTS.
- ✓ **Kocuk, Dey, Sun**, *New formulation and strong MISOCP relaxations for AC optimal transmission switching problem*, IEEE TPWRS 32(6):4161–4170, 2017 ([arXiv:1510.02064](https://arxiv.org/abs/1510.02064)).
- **Bestuzheva, Hijazi, Coffrin**, *Convex relaxations for quadratic on/off constraints and applications to optimal transmission switching*, INFORMS J. Comput. (2020). ✓ title (OSTI).
- ✓ **Fattahi, Lavaei, Atamtürk**, *A bound strengthening method for optimal transmission switching in power systems*, IEEE TPWRS 34(1):280–291, 2019.
- ✓ **Aigner, Burlacu, Liers, Martin**, *Solving AC optimal power flow with discrete decisions to global optimality*, INFORMS J. Comput. 35(2):458–474, 2023.
- ✓ **Pineda, Morales, Jiménez-Cordero**, *Learning-assisted optimization for transmission switching*, TOP 2024 ([arXiv:2304.07269](https://arxiv.org/abs/2304.07269)). DC-OTS.
- ✓ **Pineda & Morales**, *The sweet spot of bound tightening for topology optimization* ([arXiv:2507.16496](https://arxiv.org/abs/2507.16496)).

**Search result:** no work found that applies moment/Lasserre relaxations of order ≥ 2 to AC-OTS. SDP (Shor-level) relaxations for OTS exist, e.g., Fattahi et al. and others.

## 6. Unit commitment with AC power flow

- ✓ **Castillo, Laird, Silva-Monroy, Watson, O'Neill**, *The unit commitment problem with AC optimal power flow constraints*, IEEE TPWRS 31(6), 2016. Outer approximation.
- ✓ **Madani, Atamtürk, Davoudi**, *Scalable unit commitment with AC power flow via semidefinite programming relaxation* ([arXiv:1707.03541](https://arxiv.org/abs/1707.03541), BCOL report 17.03, 2017).
  - Proposes a "third-order SDP" (TSDP) relaxation for UC with AC power flow.
  - **Must read closely.** It is the nearest thing to a higher-order moment relaxation for AC-UC.
- **Fattahi, Ashraphijuo, Lavaei, Atamtürk (2017)**, *Conic relaxations of the unit commitment problem*, Energy 134.
- ✓ **Paredes, Martins, Soares**, *A low-rank rounding heuristic for semidefinite relaxation of hydro unit commitment problems* ([arXiv:1711.06194](https://arxiv.org/abs/1711.06194)). **Closest "SDP + rounding for UC" work**, but it uses a hydro model and rounds from a low-rank approximation of the Shor SDP, not from higher-order moments.
- ✓ **Paredes, Martins, Soares, Ye**, *Benders' decomposition of the unit commitment problem with semidefinite relaxation of AC power flow constraints*, EPSR 192:106965, 2021 ([arXiv:1903.02628](https://arxiv.org/abs/1903.02628)).
- ✓ **Gómez, Göttlich, Ríos, Salgado**, *Relax-and-round strategies for solving the unit commitment problem with AC power flow constraints*, Optimization and Engineering, 2026 ([arXiv:2501.11355](https://arxiv.org/abs/2501.11355)). Rounds the NLP relaxation (after rescaling), with rounding formulas that enforce combinatorial constraints. **This is the natural baseline for our UC rounding experiments.**
- ✓ *Sequential relaxation of unit commitment with AC transmission constraints* ([arXiv:1809.09812](https://arxiv.org/pdf/1809.09812)).

## 7. Bound tightening, QC, and strengthening relaxations

- **Coffrin, Hijazi, Van Hentenryck (2015)**, *Strengthening convex relaxations with bound tightening for power network optimization*, CP 2015.
- **Chen, Atamtürk, Oren (2016)**, *Bound tightening for the alternating current optimal power flow problem*, IEEE TPWRS 31(5).
- **Coffrin, Hijazi, Van Hentenryck (2016)**, *The QC relaxation: a theoretical and computational study on optimal power flow*, IEEE TPWRS 31(4).
- ✓ **Coffrin, Hijazi, Van Hentenryck**, *Strengthening the SDP relaxation of AC power flows with convex envelopes, bound tightening, and lifted nonlinear cuts* ([arXiv:1512.04644](https://arxiv.org/abs/1512.04644)); IEEE TPWRS 2017.
- **Kocuk, Dey, Sun (2016)**, *Strong SOCP relaxations for the optimal power flow problem*, Oper. Res. 64(6).
- ✓ **Sundar, Nagarajan, Misra, Lu, Coffrin, Bent**, *Optimization-based bound tightening using a strengthened QC-relaxation of the optimal power flow problem* ([arXiv:1809.04565](https://arxiv.org/abs/1809.04565)).

## 8. Machine learning for relaxation / cut / bound selection

- ✓ **Baltean-Lugojan, Bonami, Misener, Tramontani (2018)**, *Scoring positive semidefinite cutting planes for quadratic optimization via trained neural networks*, Optimization Online 6943 ([code](https://github.com/rb2309/SDPCutSel-via-NN)). A neural net scores which small SDP cuts to add. **This is the closest analogue to learning which cliques get higher order.**
- **Tang, Agrawal, Faenza (2020)**, *Reinforcement learning for integer programming: learning to cut*, ICML.
- **Bengio, Lodi, Prouvost (2021)**, *Machine learning for combinatorial optimization: a methodological tour d'horizon*, EJOR.
- ✓ **Cengil, Nagarajan, Bent, et al.**, *Learning to accelerate tightening of convex relaxations of the AC optimal power flow problem*, Comput. Optim. Appl., 2025. A learned dynamic policy chooses subsets of voltage-magnitude and angle-difference variables for OBBT. **This is directly relevant to thread 1, and is a prior on ML-OBBT for OPF.**
- ✓ **Deza & Khalil (2023)**, *Machine learning for cutting planes in integer programming: a survey*, IJCAI.

## 9. Takeaways and open gaps

1. **No prior work found that samples binary decisions from Lasserre pseudo-distributions for AC-OTS/AC-UC.**
   - Nearest neighbors: Paredes et al. (low-rank rounding of the Shor SDP for hydro UC), Gómez et al. (NLP relax-and-round for UC-ACOPF), Madani et al. (third-order SDP for UC), Eggen–Stein–Volkwein (moment-SOS + granularity rounding).
   - The TCS rounding literature supplies the theory (local consistency, conditioning) but not the continuous coupling.
2. **Scalability:** CS-TSSOS certified AC-OPF on very large systems, so combining term sparsity with binary reduction is the most promising route beyond small test cases.
3. **Thread 1 has clear priors:** Josz–Molzahn (multi-ordered hierarchy), Molzahn–Hiskens (mismatch-based order selection), Cengil et al. (ML-OBBT), Baltean-Lugojan et al. (ML-scored SDP cuts). The novelty would be (a) learning over a *mixed* action space (order / OBBT / QC / cone type) and (b) using rounding quality as the reward.
4. **To read in detail next:**
   - Madani–Atamtürk–Davoudi (is their third-order SDP a partial order-2 moment relaxation?)
   - Gómez et al. (baseline design)
   - Aigner et al. (global solution methods for comparison)
   - arXiv:2501.14281 (CD kernels)
   - Guruswami–Sinop (choosing which variables to condition on)
