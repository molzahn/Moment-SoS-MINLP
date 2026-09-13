# Polynomial formulations of AC-OPF, AC-UC, and AC-OTS, and their moment relaxations

This note documents exactly what `src/power.jl`, `src/relaxation.jl`, and `src/rounding.jl` implement. The math is written so it can be pasted into Overleaf.

## 1. Data and notation

- Data come from PowerModels (`PowerModels.parse_file`, per-unit, with PowerModels' data corrections).
- Formulations are built from `PowerModels.build_ref`. For fixed binaries, the rectangular POP and PowerModels' polar `ACPPowerModel` are therefore the *same* problem. This is verified in `scripts/validate_formulations.jl`, where objectives agree to 1e-4 on case5/9/14.
- Buses $\mathcal N$, reference bus $r$, generators $\mathcal G$ ($\mathcal G_i$ at bus $i$), branches $\mathcal L$ with from/to buses $(f,t)$.
- Voltage in rectangular coordinates: $V_i = e_i + \mathrm{j} f_i$, with $f_r = 0$ substituted (not an equality constraint) and $e_r \ge 0$.
- Squared magnitude $|V_i|^2 = e_i^2 + f_i^2$.
- For branch $l=(f,t)$:
  $$c_l = e_f e_t + f_f f_t = |V_f||V_t|\cos(\theta_f-\theta_t),\qquad s_l = f_f e_t - e_f f_t = |V_f||V_t|\sin(\theta_f-\theta_t).$$

### Branch flows (π-model with complex tap $T=\tau e^{\mathrm j\sigma}$)

Let $y = (r+\mathrm j x)^{-1}$ and
$$A = \frac{y + y^{sh}_{fr}}{|T|^2},\quad B = -\frac{y}{\bar T},\quad C = -\frac{y}{T},\quad D = y + y^{sh}_{to}.$$
Then $I_{fr} = A V_f + B V_t$ and $I_{to} = C V_f + D V_t$, so $S = V\bar I$ gives the quadratic polynomials
$$\begin{aligned}
P^{fr}_l(V) &= \Re A\,|V_f|^2 + \Re B\, c_l + \Im B\, s_l, &
Q^{fr}_l(V) &= -\Im A\,|V_f|^2 + \Re B\, s_l - \Im B\, c_l,\\
P^{to}_l(V) &= \Re D\,|V_t|^2 + \Re C\, c_l - \Im C\, s_l, &
Q^{to}_l(V) &= -\Im D\,|V_t|^2 - \Re C\, s_l - \Im C\, c_l .
\end{aligned}$$
These match PowerModels' `constraint_ohms_yt_from/to` term by term.

## 2. AC-OPF (degree 2 in $(e,f,p^g,q^g)$)

Generator outputs are kept as explicit variables. This is uniform with UC and keeps cost terms at degree 2.

$$\begin{aligned}
\min\ & \textstyle\sum_{g} \big(c_{2,g} (p^g_g)^2 + c_{1,g} p^g_g + c_{0,g}\big)\\
\text{s.t. }& \underline v_i^2 \le e_i^2+f_i^2 \le \overline v_i^2 && i\in\mathcal N\\
& \underline p_g \le p^g_g \le \overline p_g,\ \ \underline q_g \le q^g_g \le \overline q_g && g\in\mathcal G\\
& \textstyle\sum_{g\in\mathcal G_i} p^g_g - p^d_i - g^{sh}_i |V_i|^2 = \sum_{l\in\mathcal L_i} P_l^{\cdot}(V) && i\in\mathcal N\\
& \textstyle\sum_{g\in\mathcal G_i} q^g_g - q^d_i + b^{sh}_i |V_i|^2 = \sum_{l\in\mathcal L_i} Q_l^{\cdot}(V) && i\in\mathcal N\\
& \begin{bmatrix} \overline s_l^2 & P^{\cdot}_l(V) & Q^{\cdot}_l(V)\\ P^{\cdot}_l(V) & 1 & 0\\ Q^{\cdot}_l(V) & 0 & 1\end{bmatrix} \succeq 0 && l\in\mathcal L,\ \text{both ends}\\
& \tan(\underline\theta_l)\, c_l \le s_l \le \tan(\overline\theta_l)\, c_l && l\in\mathcal L
\end{aligned}$$

- The thermal limit is written as a **polynomial matrix inequality** (PMI), the Schur complement of $\overline s^2 - P^2 - Q^2 \ge 0$. Its entries have degree 2, so it can be enforced at order 1 as an LMI in the moments. The scalar form has degree 4 and would require order 2.
- The angle limits imply $c_l\ge 0$ when $\underline\theta_l < 0 < \overline\theta_l$.

## 3. Single-period AC-UC

Add a binary $u_g\in\{0,1\}$ for each commitable unit, and replace the generator bounds and costs with
$$\underline p_g u_g \le p^g_g \le \overline p_g u_g,\qquad \underline q_g u_g \le q^g_g \le \overline q_g u_g,\qquad \text{cost}_g = c_{2,g}(p^g_g)^2 + c_{1,g}p^g_g + c_{0,g} u_g .$$

- All constraints remain degree 2, so **order 1 already couples $u$ with $p^g, q^g$**.
- Optional valid inequality (on by default; valid because all branch resistances and shunt conductances are nonnegative, so losses are nonnegative):
  $$\textstyle\sum_{g} \overline p_g u_g \ge \sum_i p^d_i \qquad(\texttt{capacity\_cut}).$$
- Not yet included: startup/shutdown costs, min up/down times, ramping, reserves (multi-period).

## 4. AC-OTS

For each switchable branch $l$, add a binary $z_l$ and explicit flow variables $p^{fr}_l, q^{fr}_l, p^{to}_l, q^{to}_l$. These replace $P_l(V), Q_l(V)$ in the balance equations. Fixed branches keep the §2 formulation.

**(a) Exact polynomial model (degree 3):**
$$p^{fr}_l = z_l P^{fr}_l(V),\ \ q^{fr}_l = z_l Q^{fr}_l(V),\ \ p^{to}_l = z_l P^{to}_l(V),\ \ q^{to}_l = z_l Q^{to}_l(V),$$
$$z_l\big(\tan(\overline\theta_l)c_l - s_l\big) \ge 0,\qquad z_l\big(s_l - \tan(\underline\theta_l)c_l\big)\ge 0,\qquad (p^{\cdot}_l)^2 + (q^{\cdot}_l)^2 \le \overline s_l^2 .$$
When $z_l=0$, all flows vanish and the angle limit is dropped. This is exactly PowerModels with `br_status = 0`.

**(b) Big-M on/off model (degree ≤ 2):** with $M^{fr}_l = |A|\overline v_f^2 + |B|\overline v_f\overline v_t \ge |P^{fr}_l(V)|, |Q^{fr}_l(V)|$ (similarly for $M^{to}_l$), $b_l = \min(\overline s_l, M_l)$, and $M^\theta_l = (\max|\tan\theta| + 1)\overline v_f\overline v_t$:
$$|p^{\cdot}_l| \le b_l z_l,\qquad |p^{\cdot}_l - P^{\cdot}_l(V)| \le M_l(1-z_l),\qquad \tan(\overline\theta_l)c_l - s_l + M^\theta_l(1-z_l)\ge 0,\ \text{etc.}$$

**Implementation choice.** Both (a) and (b) are included in the POP. The relaxation builder enforces a constraint of degree $d$ only in cliques of order $t \ge \lceil d/2\rceil$. So:
- **Order 1** sees only (b), a big-M "Shor + on/off" relaxation. The skipped (a) constraints are reported in `rel.skipped`.
- **Order 2** sees (a) and (b). At order 2, (b) is redundant in principle but may still help numerically.

**Not modeled:**
- Connectivity / no-islanding constraints in the relaxation. Islanded configurations are rejected during rounding and enumeration, so the relaxation is a relaxation of a *larger* set and still gives a valid bound.
- Per-island angle references.

## 5. Moment relaxation (Lasserre hierarchy) as implemented

POP: $\min f(x)$ s.t. $g_i(x)\ge 0$, $h_j(x)=0$, $G_k(x)\succeq 0$, with $x_b\in\{0,1\}$ for $b\in\mathcal B$.

1. **Binary reduction.** Every monomial is reduced with $x_b^2 = x_b$ *before* building the relaxation. Basis monomials with a repeated binary are dropped. There are no constraints $y_{x_b^2} = y_{x_b}$, so the moment matrix keeps strictly feasible points. Keeping those equalities explicitly destroys Slater's condition.
2. **Correlative sparsity.**
   - Build the interaction graph (edges between variables that appear together in the objective, any constraint, or any PMI).
   - Take a chordal extension by greedy min-degree elimination and keep its maximal cliques $C_1,\dots,C_p$.
   - Give each clique its own order $t_k$ (`order` may be a function of the clique). This is the multi-ordered hierarchy of Josz & Molzahn, and the main knob for thread 1.
   - **Clique augmentation** (`extra_supports`): additional variable sets can be added to the interaction graph, so the chordal extension places them in a common clique.
   - `build_power_pop` provides `meta["bus_binaries"]` (binaries incident to each bus, plus the bus voltage) and `meta["all_binaries"]`.
   - This matters for rounding. Without it, AC-OTS relaxations contain *no* joint moments $y_{z_iz_j}$ (see `results/experiment1_findings.md`, F3).
   - **Scaling options** (experiment 2):
     - `global_linear = k`: linear constraints with more than $k$ variables are kept out of the interaction graph and, if no clique contains them, enforced only as $L_y(g)\ge 0$ / $L_y(h)=0$.
     - `capped_augmentation`: accept augmentation supports only while binary cliques stay within a size cap.
     - `extra_cliques`: additional moment blocks added *without* re-chordalizing. The relaxation stays valid, but the running intersection property is lost.
     - Order policies: `binary_clique_order` (order 2 on cliques with a binary) and `adjacent_clique_order` (order 2 on size-capped cliques containing a binary or a variable that shares a constraint with one).
3. **Moment matrices.** $M_{t_k}(y; C_k) = \big[y_{\alpha\beta}\big]_{\alpha,\beta\in \mathbb B_{t_k}(C_k)} \succeq 0$, where $\mathbb B_t(C)$ is the set of binary-reduced monomials of degree $\le t$ in the variables of $C$.
4. **Assignment.** Each constraint is assigned to one clique that contains its support, choosing the highest order.
5. **Localizing matrices.** $M_{t_k-\lceil \deg g/2\rceil}(g\,y;C_k)\succeq 0$. Skipped if $t_k < \lceil \deg g/2\rceil$.
6. **Equalities.** $L_y(h\cdot m) = 0$ for all $m\in \mathbb B_{2t_k-\deg h}(C_k)$.
7. **PMIs.** $\big[L_y(G_{ab}\, m_i m_j)\big]_{(a,i),(b,j)}\succeq 0$ with $m_i\in\mathbb B_{t_k-\lceil \deg G/2\rceil}(C_k)$.
8. **Objective and bound.** $\min L_y(f)$, scaled by $\max|f_\alpha|$.
9. **SDP form.** The code solves the equivalent **SOS (dual) form** by default (`form = :sos`):
   $$\max\ \lambda\ \ \text{s.t.}\ \ f-\lambda = \sum_k \langle X_k, \text{(weighted) monomial matrix}_k\rangle + \sum_{h,m}\lambda_{h,m}\, h\, m\ \ (\text{coefficient-wise, after binary reduction}),\ \ X_k\succeq 0 .$$
   - There is one equality per moment. The pseudo-moments are the normalized equality duals, $y_\alpha = \text{dual}_\alpha/\text{dual}_\emptyset$.
   - Compared with the primal moment form (`form = :moment`, affine PSD constraints bridged by JuMP), this is much smaller for MOSEK: case5 AC-OPF at order 2 takes 1.0 s vs 36.6 s, with identical bounds.

**Size check** (`scripts/validate_formulations.jl`): case5 AC-OPF has a 3.84% gap at order 1 and 0% at order 2 (7 s). case9 and case14 are tight at order 1.

## 6. Rounding schemes (`src/rounding.jl`)

Let $\mu_i = y_{z_i}$ (marginal) and $Y_{ij} = y_{z_iz_j}$ (available if $z_i,z_j$ share a clique).

| scheme | description | uses |
|---|---|---|
| `threshold` | $z_i = [\mu_i\ge 1/2]$, deterministic | $\mu$ |
| `independent` | $z_i\sim\mathrm{Bernoulli}(\mu_i)$ | $\mu$ |
| `gaussian` | Gaussian copula: $g\sim N(0,R)$ with $R_{ij} = (Y_{ij}-\mu_i\mu_j)/\sqrt{\mu_i(1-\mu_i)\mu_j(1-\mu_j)}$ (PSD-projected), $z_i = [g_i < \Phi^{-1}(\mu_i)]$. Marginals are preserved exactly. | $\mu, Y$ |
| `conditional` | Random variable order. For each $j$, choose up to $2t-1$ already-sampled $B$ (largest $|R_{ij}|$, with $y_{B\cup\{j\}}$ available), then $P(z_j=1\mid z_B=a) = \frac{\sum_{T\subseteq B_0}(-1)^{|T|}y_{B_1\cup T\cup\{j\}}}{\sum_{T\subseteq B_0}(-1)^{|T|}y_{B_1\cup T}}$ ($B_1/B_0$: ones/zeros in $a$). Negative pseudo-probabilities are clipped, and the clipped mass is reported. | all binary moments up to degree $2t$ |
| `dive` | Random order: sample $z_j\sim\mathrm{Bernoulli}(\mu_j)$ from the *current* relaxation, fix $z_j$ (substitute into the POP), re-solve. If infeasible, flip. | re-solved relaxations |

**Theory.**
- For pure binary problems, the pseudo-moments of any $\le t$ binaries at order $t$ form a genuine distribution (Rothvoß §2.4).
- Conditioning on $z_j = 1$, i.e. $y'_\alpha = y_{\alpha\cup j}/y_j$, gives a feasible order-$(t-1)$ moment sequence. With the reduction $z_j^2 = z_j$, the localizing matrix of $z_j$ is a principal submatrix of $M_t(y)$. The same argument goes through with continuous variables present.
- `conditional` uses moments up to degree $2t$, beyond the guaranteed range, so it clips.
- `dive` re-solves instead of conditioning in closed form. It is more expensive but uses all constraints.

**Recovery.**
- Each distinct sampled configuration is screened (network connected; active capacity ≥ load), then solved with PowerModels `ACPPowerModel` + Ipopt.
- Results are cached per configuration. For small instances the cache is pre-filled by full enumeration (`scripts/enumerate_instances.jl`).

## 7. Known limitations (to revisit)

- **Ipopt is local.** "Enumerated optimum" means best local solution from a flat start, per configuration. Where the order-2 bound matches it, global optimality is certified.
- **Explicit $p^g, q^g$ variables enlarge cliques.** Eliminating them (as in the original moment-OPF papers) is possible for OPF/OTS but not with $u_g$ coupling.
- **Moment matrices are built as dense JuMP affine matrices.** This is fine for ≤ ~30 variables per clique, but term sparsity (TSSOS-style) will be needed to scale.
