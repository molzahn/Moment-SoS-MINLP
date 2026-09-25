# Sparse moment (Lasserre) relaxation of a POP with binary variables.
#
# * Binary reduction x^2 = x is applied directly to the monomial basis (no equality constraints),
#   which keeps the SDP strictly feasible.
# * Correlative sparsity: one moment matrix per maximal clique of a chordal extension.
# * Each clique can have its own relaxation order (knob for adaptive/learned order selection).
# * A constraint of degree d is enforced in a clique of order t only if ceil(d/2) <= t;
#   otherwise it is skipped and reported (e.g. degree-3 switching constraints at order 1).
#
# Two equivalent SDP forms are available:
#   form = :sos    (default) SOS/dual form: PSD Gram matrices + one coefficient-matching equality per
#                  moment. Pseudo-moments are recovered from the equality duals. Much smaller Schur
#                  complement for interior-point solvers.
#   form = :moment primal moment form with affine PSD constraints (bridged by JuMP; memory-hungry).

struct MomentRelaxation
    status::MOI.TerminationStatusCode
    primal_status::MOI.ResultStatusCode   # status of the pseudo-moment vector y
    bound::Float64                         # lower bound (SOS form: certified from the inexact solution, see info)
    primal_objective::Float64              # L_y(f) at the recovered pseudo-moments
    build_time::Float64
    solve_time::Float64
    cliques::Vector{Vector{Int}}
    orders::Vector{Int}
    y::Dict{Vector{Int},Float64}           # pseudo-moments (binary-reduced monomials)
    skipped::Dict{String,Int}              # constraint tags skipped because of insufficient order
    n_moments::Int
    psd_sizes::Vector{Int}
    rank_ratio::Vector{Float64}            # λ2/λ1 of each clique's order-1 moment matrix
    info::Dict{String,Any}                 # solver diagnostics (iterations, objectives, feasibility of y)
end

function monomial_basis(vars::Vector{Int}, d::Int, isbin::AbstractVector{Bool})
    out = [Int[]]
    function rec(prefix, start, remaining)
        for k in start:length(vars)
            v = vars[k]
            (isbin[v] && !isempty(prefix) && prefix[end] == v) && continue
            m = [prefix; v]
            push!(out, m)
            remaining > 1 && rec(m, k, remaining - 1)
        end
    end
    d > 0 && rec(Int[], 1, d)
    sort!(out; by = m -> (length(m), m))
    return out
end

monoprod(a::Vector{Int}, b::Vector{Int}, isbin) = reduce_mono(sort!(vcat(a, b)), isbin)

_ceil_half(d) = (d + 1) ÷ 2

"""
Supports (variable sets) defining the interaction graph. Linear constraints with more than
`global_linear` variables are kept out of the graph; if no clique contains their support they are
enforced on first moments only (L(g) >= 0, L(h) = 0).
"""
function interaction_supports(obj, ineqs, eqs, pmis; global_linear = typemax(Int))
    supports = Vector{Vector{Int}}()
    for m in keys(obj.terms)
        isempty(m) || push!(supports, unique(m))
    end
    in_graph(p) = !(degree(p) <= 1 && length(support(p)) > global_linear)
    append!(supports, filter(!isempty, support.(filter(in_graph, ineqs))))
    append!(supports, filter(!isempty, support.(filter(in_graph, eqs))))
    for G in pmis
        s = sort!(unique(reduce(vcat, support.(G))))
        isempty(s) || push!(supports, s)
    end
    return supports, sort!(unique(reduce(vcat, supports)))
end

function interaction_supports(pop::POP; global_linear = typemax(Int))
    r(p) = reduce_poly(p, pop.isbin)
    return interaction_supports(r(pop.obj), r.(pop.ineqs), r.(pop.eqs), [map(r, G) for G in pop.pmis];
        global_linear = global_linear)
end

"""
    capped_augmentation(pop, candidates; maxclique, global_linear) -> accepted supports

Greedily add candidate supports (in the given order) as long as every clique of the chordal extension
that contains a binary variable has at most `maxclique` variables.
"""
function capped_augmentation(pop::POP, candidates::Vector{Vector{Int}}; maxclique::Int = 14,
    global_linear = typemax(Int))
    supports, active = interaction_supports(pop; global_linear = global_linear)
    act = Set(active)
    maxbin(cl) = maximum((length(c) for c in cl if any(pop.isbin[v] for v in c)); init = 0)
    accepted = Vector{Vector{Int}}()
    for cand in candidates
        s = filter(in(act), cand)
        length(s) > 1 || continue
        trial = [supports; accepted; [s]]
        if maxbin(chordal_cliques(trial)) <= maxclique
            push!(accepted, s)
        end
    end
    return accepted
end

"""
    solve_moment_relaxation(pop; order = 1, sparse = true, form = :sos, optimizer = mosek_optimizer())

`order` is either an Int or a function `clique::Vector{Int} -> Int` (see `binary_clique_order`).
"""
function solve_moment_relaxation(pop::POP; order = 1, sparse::Bool = true, cliques = nothing, form::Symbol = :sos,
    extra_supports = Vector{Vector{Int}}(), extra_cliques = Vector{Vector{Int}}(), global_linear::Int = typemax(Int),
    optimizer = mosek_optimizer(), silent::Bool = true, normalize::Bool = true,
    skip_high_order_tags = String[], solver_params = Dict{String,Any}(), diagnostics::Bool = true,
    quotient_basis::Bool = false, scale_vars::Bool = true, out_of_graph_tags = String[], scalar_tags = String[],
    certify::Symbol = :implicit, bundle_params = Dict{Symbol,Any}(), parity_vars = nothing)
    t0 = time()
    if scale_vars
        # substitute x_i = s_i x̂_i with s_i = max(|lb_i|, |ub_i|) so that every variable lies in [-1, 1];
        # pseudo-moments are mapped back to the original variables afterwards
        svec = [pop.isbin[i] ? 1.0 : (m = max(abs(pop.lb[i]), abs(pop.ub[i])); isfinite(m) && m > 0 ? m : 1.0)
                for i in eachindex(pop.lb)]
        sp(p) = Poly(Dict(m => c * prod((svec[v] for v in m); init = 1.0) for (m, c) in p.terms))
        q = deepcopy(pop)
        q.obj = sp(pop.obj)
        q.ineqs = sp.(pop.ineqs)
        q.eqs = sp.(pop.eqs)
        q.pmis = [map(sp, G) for G in pop.pmis]
        # the cones must be substituted too, or they would still reference unscaled variables
        q.socs = [(sp(a), sp(b), [sp(x) for x in xs]) for (a, b, xs) in pop.socs]
        q.lb = pop.lb ./ svec
        q.ub = pop.ub ./ svec
        r = solve_moment_relaxation(q; order = order, sparse = sparse, cliques = cliques, form = form,
            extra_supports = extra_supports, extra_cliques = extra_cliques, global_linear = global_linear,
            optimizer = optimizer, silent = silent, normalize = normalize, skip_high_order_tags = skip_high_order_tags,
            solver_params = solver_params, diagnostics = diagnostics, quotient_basis = quotient_basis, scale_vars = false,
            out_of_graph_tags = out_of_graph_tags, scalar_tags = scalar_tags, certify = certify, bundle_params = bundle_params, parity_vars = parity_vars)
        y = Dict(m => v * prod((svec[i] for i in m); init = 1.0) for (m, v) in r.y)
        r.info["scale_vars"] = true
        return MomentRelaxation(r.status, r.primal_status, r.bound, r.primal_objective, r.build_time, r.solve_time,
            r.cliques, r.orders, y, r.skipped, r.n_moments, r.psd_sizes, r.rank_ratio, r.info)
    end
    isbin = pop.isbin
    obj = reduce_poly(pop.obj, isbin)
    ineqs = [reduce_poly(g, isbin) for g in pop.ineqs]
    eqs = [reduce_poly(h, isbin) for h in pop.eqs]
    pmis = [map(p -> reduce_poly(p, isbin), G) for G in pop.pmis]
    if normalize   # scale every constraint to unit max-coefficient (same feasible set)
        nrm(p) = (m = maximum(abs, values(p.terms); init = 0.0); m > 0 ? (1 / m) * p : p)
        ineqs = nrm.(ineqs)
        eqs = nrm.(eqs)
        pmis = [(m = maximum(p -> maximum(abs, values(p.terms); init = 0.0), G); m > 0 ? map(p -> (1 / m) * p, G) : G) for G in pmis]
    end
    tagkey(tag) = first(split(tag, '['))
    drop_high(tag, k) = k > 0 && orders[k] >= 2 && tagkey(tag) in skip_high_order_tags

    # linear constraints whose tag is in `out_of_graph_tags` do not shape the cliques; they are still localized
    # in a clique that happens to contain them, otherwise enforced on first moments
    graph_ineqs = [g for (g, tag) in zip(ineqs, pop.ineq_tags) if !(degree(g) <= 1 && tagkey(tag) in out_of_graph_tags)]
    supports, active = interaction_supports(obj, graph_ineqs, eqs, pmis; global_linear = global_linear)
    # extra supports force sets of variables into a common clique (e.g. to create joint binary moments)
    for s in extra_supports
        s = filter(in(Set(active)), s)
        length(s) > 1 && push!(supports, s)
    end
    if cliques === nothing
        cliques = sparse ? chordal_cliques(supports) : [active]
    end
    # extra moment blocks added without re-chordalizing (valid relaxation; running intersection not required)
    for c in extra_cliques
        c = sort(filter(in(Set(active)), c))
        length(c) > 1 && !any(issubset(c, k) for k in cliques) && push!(cliques, c)
    end
    orders = order isa Integer ? fill(Int(order), length(cliques)) : [Int(order(c)) for c in cliques]

    function assign(s::Vector{Int}; allow_global = false)
        best = 0
        for (k, c) in enumerate(cliques)
            if issubset(s, c) && (best == 0 || orders[k] > orders[best])
                best = k
            end
        end
        best == 0 && !allow_global && error("support $(s) not contained in any clique")
        return best
    end

    # A constraint whose *full support* lies in no clique can still be imposed on the shared moment
    # vector as L_y(g) >= 0 / L_y(h) == 0, provided every monomial it uses exists -- i.e. each
    # monomial's own support lies in some clique. This is the matrix-completion reading of the
    # sparse relaxation (Fukuda et al.; Jabr; Molzahn et al.) rather than the stricter Waki et al.
    # condition that each constraint live inside one clique, and it is what lets the clique graph be
    # the bus adjacency instead of the constraint-support graph. On case200 that is the difference
    # between 23-bus and 9-bus cliques, i.e. an order-2 block of 1128 against 190.
    #
    # Inert unless `cliques` is supplied: the default support-derived decomposition already covers
    # every constraint, so `covered_by_monomials` is never reached.
    function covered_by_monomials(p::Poly)
        for m in keys(p.terms)
            isempty(m) && continue
            s = sort(unique(m))
            any(issubset(s, c) for c in cliques) ||
                error("monomial $(s) of a constraint is in no clique; the clique decomposition " *
                      "does not cover this problem")
        end
        return true
    end
    # ---------------------------------------------------------------------------------------
    # Parity block splitting.
    #
    # When every polynomial in the problem has even total degree in a designated variable set S,
    # the problem is invariant under x_S -> -x_S, so there is an optimal moment vector with every
    # odd-in-S moment zero. Each moment and localizing matrix then block-diagonalises by the
    # parity of its basis monomials, and the two halves can be imposed as independent PSD
    # constraints. This is the scheme in Molzahn's msos reference (firstRow2mon_moment.m), gated
    # there by `anyFirstOrderMon`.
    #
    # For rectangular AC-OPF, S is the set of voltage variables e[i], f[i]: the power flow
    # equations carry only even powers of the voltage. Measured saving at order 2 is 21-34% of the
    # flops for the clique sizes this code produces, and at order 1 the even half collapses to the
    # constant so the block is just the W matrix.
    #
    # The precondition is *checked*, never assumed: one odd constraint makes the split invalid.
    # On by default: the voltage variables recorded by build_power_pop. The precondition below is
    # the interlock -- if any constraint is odd in them (e.g. `eref` is still present) nothing is
    # split and `info["parity_blocked_by"]` says which family objected.
    Sp = parity_vars === nothing ? collect(get(pop.meta, "voltage_vars", Int[])) : collect(parity_vars)
    Spset = Set(Sp)
    vparity(m) = isodd(count(v -> v in Spset, m))
    function poly_parity(p::Poly)          # :even, :odd, or :mixed
        e = o = false
        for m in keys(p.terms)
            vparity(m) ? (o = true) : (e = true)
        end
        o && e && return :mixed
        return o ? :odd : :even
    end
    parity_ok = !isempty(Sp)
    parity_block = String[]
    if parity_ok
        for (p, tag) in vcat(collect(zip(ineqs, pop.ineq_tags)), collect(zip(eqs, pop.eq_tags)))
            poly_parity(p) == :even || (parity_ok = false; push!(parity_block, tagkey(tag)))
        end
        for (G, tag) in zip(pmis, pop.pmi_tags), q in G
            poly_parity(q) == :even || (parity_ok = false; push!(parity_block, tagkey(tag)))
        end
        poly_parity(obj) == :even || (parity_ok = false; push!(parity_block, "OBJECTIVE"))
    end
    "Split a basis by monomial parity; returns the halves that are non-empty."
    split_basis(B) = parity_ok ? filter(!isempty, [filter(m -> !vparity(m), B),
                                                   filter(vparity, B)]) : [B]

    bases = Dict{Tuple{Int,Int},Vector{Vector{Int}}}()
    basis(k, d) = get!(() -> monomial_basis(cliques[k], d, isbin), bases, (k, d))
    skipped = Dict{String,Int}()
    skip!(tag) = (key = tagkey(tag); skipped[key] = get(skipped, key, 0) + 1)
    dropped = Dict{String,Int}()
    drop!(tag) = (key = tagkey(tag); dropped[key] = get(dropped, key, 0) + 1)

    # Rotated second-order cones, imposed on first moments. Only the binary reduction is applied:
    # 2ab >= ||x||^2 is invariant under scaling a, b and x by a COMMON factor but not under the
    # per-polynomial normalisation used for the other families, so normalising them separately
    # would change the cone.
    socs = [(reduce_poly(a, isbin), reduce_poly(b, isbin), [reduce_poly(x, isbin) for x in xs])
            for (a, b, xs) in pop.socs]

    # PSD blocks: (G, B) means [L_y(G_ab m_i m_j)] ⪰ 0; equality blocks: (h, multipliers)
    one = fill(Poly(1.0), 1, 1)
    blocks = Tuple{Matrix{Poly},Vector{Vector{Int}}}[]
    eqblocks = Tuple{Poly,Vector{Vector{Int}}}[]
    n_pivots = 0
    for k in eachindex(cliques)
        B = basis(k, orders[k])
        if quotient_basis && orders[k] >= 2
            # Degree-2 equalities h with support in the clique force L(h·h) = 0, i.e. a singular moment
            # matrix. Remove one pivot monomial per independent equality from the basis and impose
            # L(h·m) = 0 for all m in the full basis, so that the full moment matrix is a congruence of
            # the reduced one (exact reformulation, strictly feasible reduced matrix).
            Ek = [h for h in eqs if degree(h) == 2 && !is_constant(h) && issubset(support(h), cliques[k])]
            piv = _pivot_monomials(Ek)
            if !isempty(piv)
                pset = Set(piv)
                B = filter(m -> !(m in pset), B)
                n_pivots += length(piv)
                full = basis(k, 2 * orders[k] - 2)
                for h in Ek
                    push!(eqblocks, (h, full))
                end
            end
        end
        for Bp in split_basis(B)
            push!(blocks, (one, Bp))
        end
    end
    n_moment_blocks = length(blocks)   # localizing blocks follow; no longer one per clique
    for (g, tag) in zip(ineqs, pop.ineq_tags)
        is_constant(g) && continue
        k = assign(support(g); allow_global = true)
        if k == 0
            degree(g) <= 1 || covered_by_monomials(g)
            push!(blocks, (fill(g, 1, 1), [Int[]]))
            continue
        end
        if drop_high(tag, k)
            drop!(tag)
            continue
        end
        dd = orders[k] - _ceil_half(degree(g))
        # inequalities tagged in `scalar_tags` are enforced as L(g) >= 0 only (no localizing matrix), as at order 1
        tagkey(tag) in scalar_tags && (dd = min(dd, 0))
        if dd < 0
            skip!(tag)
        else
            # An even weight preserves the parity of the multiplier basis, so the localizing
            # matrix splits like the moment matrix. An odd weight *swaps* the halves (the block is
            # anti-block-diagonal, not block-diagonal) and a mixed weight mixes them; neither
            # splits. `parity_ok` already guarantees every weight here is even, but the check is
            # local so this stays correct if that precondition is ever relaxed.
            for Bp in (poly_parity(g) == :even ? split_basis(basis(k, dd)) : [basis(k, dd)])
                push!(blocks, (fill(g, 1, 1), Bp))
            end
        end
    end
    for (G, tag) in zip(pmis, pop.pmi_tags)
        k = assign(sort!(unique(reduce(vcat, support.(G)))))
        dd = orders[k] - _ceil_half(maximum(degree.(G)))
        if dd < 0
            skip!(tag)
        else
            for Bp in (all(poly_parity(q) == :even for q in G) ? split_basis(basis(k, dd)) :
                       [basis(k, dd)])
                push!(blocks, (G, Bp))
            end
        end
    end
    for (h, tag) in zip(eqs, pop.eq_tags)
        is_constant(h) && continue
        k = assign(support(h); allow_global = true)
        if k == 0
            degree(h) <= 1 || covered_by_monomials(h)
            push!(eqblocks, (h, [Int[]]))
            continue
        end
        if drop_high(tag, k)
            drop!(tag)
            continue
        end
        dd = 2 * orders[k] - degree(h)
        dd < 0 ? skip!(tag) : push!(eqblocks, (h, basis(k, dd)))
    end
    for m in keys(obj.terms)
        isempty(m) && continue
        k = assign(unique(m))
        length(m) <= 2 * orders[k] || error("objective monomial $(m) exceeds relaxation degree")
        # A monomial no block can produce becomes `0 == coefficient` further down, i.e. a silently
        # infeasible SDP. Parity is the way that can now happen.
        !parity_ok || !vparity(m) ||
            error("objective monomial $(m) is odd in the parity variables; parity splitting " *
                  "would make the relaxation infeasible")
    end
    scale = maximum(abs, values(obj.terms); init = 1.0)
    psd_sizes = [size(G, 1) * length(B) for (G, B) in blocks]

    model = Model(optimizer)
    silent && set_silent(model)
    for (k, v) in solver_params
        set_attribute(model, k, v)
    end
    if form == :sos
        maxabs = [max(abs(pop.lb[i]), abs(pop.ub[i])) for i in eachindex(pop.lb)]
        yv, bound, pobj, st, ps, nmom, cert = _solve_sos_form(model, obj, blocks, eqblocks, socs, isbin, scale, maxabs;
            certify = certify, bundle_params = bundle_params, groups = square_groups(ineqs))
    elseif form == :moment
        yv, bound, pobj, st, ps, nmom = _solve_moment_form(model, obj, blocks, eqblocks, socs, isbin, scale)
        cert = Dict{String,Any}()
    else
        error("unknown form $form")
    end
    build_time = time() - t0 - solve_time(model)

    ratios = Float64[]
    if !isempty(yv)
        for k in eachindex(cliques)
            # With parity splitting the odd moments are gone, and `get(yv, mm, 0.0)` would turn
            # every missing entry into a silent zero -- an eigenvalue ratio computed on the wrong
            # matrix. Use the odd half of the basis alone, whose Gram matrix is exactly the W
            # block (what msos extracts as W{i}); without splitting, keep the full basis.
            B = basis(k, 1)
            Bk = parity_ok ? filter(vparity, B) : B
            isempty(Bk) && (push!(ratios, 0.0); continue)
            M1 = [(mm = monoprod(Bk[i], Bk[j], isbin); isempty(mm) ? 1.0 : get(yv, mm, 0.0))
                  for i in eachindex(Bk), j in eachindex(Bk)]
            ev = sort(eigvals(Symmetric(M1)); rev = true)
            push!(ratios, length(ev) > 1 ? ev[2] / ev[1] : 0.0)
        end
    end
    info = Dict{String,Any}("dropped" => dropped, "normalize" => normalize, "n_pivots" => n_pivots,
        "parity_split" => parity_ok, "n_moment_blocks" => n_moment_blocks)
    isempty(Sp) || parity_ok ||
        (info["parity_blocked_by"] = sort(unique(parity_block)))
    merge!(info, cert)
    try
        info["iterations"] = MOI.get(model, MOI.BarrierIterations())
    catch
    end
    try
        info["primal_obj"] = objective_value(model) * scale
        info["dual_obj"] = dual_objective_value(model) * scale
        info["rel_gap"] = abs(info["primal_obj"] - info["dual_obj"]) / max(1.0, abs(info["primal_obj"]))
    catch
    end
    if diagnostics && !isempty(yv)
        info["min_eig_moment"], info["min_eig_localizing"], info["max_eq_residual"] =
            _moment_feasibility(yv, blocks, eqblocks, isbin, n_moment_blocks)
    end
    return MomentRelaxation(st, ps, bound, pobj, build_time, solve_time(model), cliques, orders, yv,
        skipped, nmom, psd_sizes, ratios, info)
end

function _solve_sos_form(model, obj, blocks, eqblocks, socs, isbin, scale, maxabs; certify::Symbol = :implicit,
    bundle_params = Dict{Symbol,Any}(), groups = Tuple{Vector{Int},Float64}[])
    coef = Dict{Vector{Int},AffExpr}()
    addc!(m, c, v) = add_to_expression!(get!(() -> AffExpr(0.0), coef, m), c, v)
    Xs = Any[]
    for (G, B) in blocks
        ng, nb = size(G, 1), length(B)
        n = ng * nb
        if n == 1
            X = reshape([@variable(model, lower_bound = 0.0)], 1, 1)
        else
            X = @variable(model, [1:n, 1:n], PSD)
        end
        push!(Xs, X)
        for a in 1:ng, b in 1:ng, i in 1:nb, j in 1:nb
            v = X[(a-1)*nb+i, (b-1)*nb+j]
            mij = monoprod(B[i], B[j], isbin)
            for (mono, c) in G[a, b].terms
                addc!(monoprod(mono, mij, isbin), c, v)
            end
        end
    end
    for (h, mults) in eqblocks
        for m in mults
            λ = @variable(model)
            for (mono, c) in h.terms
                addc!(monoprod(mono, m, isbin), c, λ)
            end
        end
    end
    # A rotated SOC on the moments (a, b, x) contributes a dual multiplier from the same cone --
    # RotatedSecondOrderCone is self-dual, so <mu, (a,b,x)> >= 0 on the feasible set and the
    # certificate stays a valid lower bound. This is the conic analogue of the nonnegative
    # multiplier used for a scalar inequality above.
    for (a, b, xs) in socs
        mu = @variable(model, [1:(2 + length(xs))])
        @constraint(model, mu in RotatedSecondOrderCone())
        for (q, v) in zip(vcat([a, b], xs), mu)
            for (mono, c) in q.terms
                addc!(reduce_mono(mono, isbin), c, v)
            end
        end
    end
    t = @variable(model)
    add_to_expression!(get!(() -> AffExpr(0.0), coef, Int[]), 1.0, t)
    monos = union(keys(coef), keys(obj.terms))
    cons = Dict{Vector{Int},ConstraintRef}()
    for m in monos
        cons[m] = @constraint(model, get(coef, m, AffExpr(0.0)) == get(obj.terms, m, 0.0) / scale)
    end
    @objective(model, Max, t)
    optimize!(model)
    st = termination_status(model)
    ok = primal_status(model) in (MOI.FEASIBLE_POINT, MOI.NEARLY_FEASIBLE_POINT)
    # Certified lower bound from an inexact SOS solution. With residuals r_α = C_α(X,λ,t) − f_α/scale and
    # Gram matrices X_k with smallest eigenvalue λ_k, for every feasible x in the variable box
    #   f(x)/scale ≥ t − Σ_α |r_α| max|x^α| − Σ_k max(−λ_k, 0) · max tr W_k(x),
    # where W_k(x) is the (weighted) monomial matrix of block k.
    cert = Dict{String,Any}()
    bound = NaN
    if ok
        mb(m) = prod((maxabs[v] for v in m); init = 1.0)
        tval = value(t)
        rsum, rmax = 0.0, 0.0
        for m in monos
            r = value(get(coef, m, AffExpr(0.0))) - get(obj.terms, m, 0.0) / scale
            rsum += abs(r) * mb(m)
            rmax = max(rmax, abs(r))
        end
        esum, emin = 0.0, Inf
        for ((G, B), X) in zip(blocks, Xs)
            Xv = size(X, 1) == 1 ? fill(value(X[1, 1]), 1, 1) : value.(X)
            λ = size(Xv, 1) == 1 ? Xv[1, 1] : eigmin(Symmetric(Xv))
            emin = min(emin, λ)
            if λ < 0
                tr = 0.0
                for a in axes(G, 1), mi in B
                    tr += sum((abs(c) * mb(mono) for (mono, c) in G[a, a].terms); init = 0.0) * mb(mi)^2
                end
                esum += -λ * tr
            end
        end
        bound = (tval - rsum - esum) * scale
        cert = Dict{String,Any}("raw_bound" => tval * scale, "sos_residual_max" => rmax,
            "residual_correction" => rsum * scale, "psd_correction" => esum * scale, "gram_min_eig" => emin,
            "certified_box" => bound)
    end
    # Certificates valid for any point (so also when MOSEK returns no feasible point), see src/certify.jl:
    # :implicit = best of box / implicit / per-block hybrid at MOSEK's point; :optimize additionally maximizes
    # the certificate (bundle_params[:method] = :smooth (default, L-BFGS on a smoothed surrogate) or :bundle).
    if certify in (:implicit, :optimize, :bundle) && has_values(model)
        tc = time()
        sc = SOSCertificate(model, coef, obj, monos, blocks, Xs, t, isbin, scale, maxabs; groups = groups)
        θ0 = value.(all_variables(model))
        haskey(cert, "raw_bound") || (cert["raw_bound"] = θ0[sc.tcol] * scale)
        cert["certified_box_tight_rho"] = first(certified_value(sc, θ0; mode = :box, gradient = false))
        cert["certified_implicit"] = first(certified_value(sc, θ0; gradient = false))
        cert["certified_hybrid"] = first(certified_value(sc, θ0; mode = :hybrid, gradient = false))
        cert["n_unrepresented_monomials"] = count(==(0), sc.repcol)
        cert["n_square_groups"] = length(groups)
        for k in ("certified_box_tight_rho", "certified_implicit", "certified_hybrid")
            isfinite(cert[k]) && (!isfinite(bound) || cert[k] > bound) && (bound = cert[k])
        end
        cert["implicit_time"] = time() - tc
        bp = Dict{Symbol,Any}(bundle_params)
        get(bp, :keep, false) && (cert["_certificate"] = (sc, θ0))
        delete!(bp, :keep)
        method = pop!(bp, :method, certify == :bundle ? :bundle : :smooth)
        if certify in (:optimize, :bundle)
            tb = time()
            Fb, _, hist = method == :bundle ? bundle_certify(sc, θ0; bp...) : smooth_certify(sc, θ0; bp...)
            cert["certified_optimized"] = Fb
            cert["optimize_method"] = string(method)
            cert["optimize_evals"] = length(hist)
            cert["optimize_time"] = time() - tb
            isfinite(Fb) && Fb > bound && (bound = Fb)
        end
    end
    yv = Dict{Vector{Int},Float64}()
    ps = dual_status(model)
    if has_duals(model)
        d0 = dual(cons[Int[]])
        for (m, c) in cons
            isempty(m) || (yv[m] = dual(c) / d0)
        end
    end
    pobj = isempty(yv) ? NaN : sum(c * (isempty(m) ? 1.0 : yv[m]) for (m, c) in obj.terms)
    return yv, bound, pobj, st, ps, length(monos) - 1, cert
end

function _solve_moment_form(model, obj, blocks, eqblocks, socs, isbin, scale)
    y = Dict{Vector{Int},VariableRef}()
    getY(m) = get!(() -> @variable(model), y, m)
    function lin(p::Poly, mult::Vector{Int})
        e = AffExpr(0.0)
        for (m, c) in p.terms
            mm = monoprod(m, mult, isbin)
            isempty(mm) ? add_to_expression!(e, c) : add_to_expression!(e, c, getY(mm))
        end
        return e
    end
    for (G, B) in blocks
        ng, nb = size(G, 1), length(B)
        M = Matrix{AffExpr}(undef, ng * nb, ng * nb)
        for a in 1:ng, b in 1:ng, i in 1:nb, j in 1:nb
            M[(a-1)*nb+i, (b-1)*nb+j] = lin(G[a, b], monoprod(B[i], B[j], isbin))
        end
        size(M, 1) == 1 ? @constraint(model, M[1, 1] >= 0) : @constraint(model, Symmetric(M) in PSDCone())
    end
    for (h, mults) in eqblocks, m in mults
        @constraint(model, lin(h, m) == 0)
    end
    # rotated second-order cones on the first moments (JuMP: 2*u*v >= ||w||^2)
    for (a, b, xs) in socs
        @constraint(model, vcat(lin(a, Int[]), lin(b, Int[]), [lin(x, Int[]) for x in xs])
                    in RotatedSecondOrderCone())
    end
    @objective(model, Min, lin(obj, Int[]) / scale)
    optimize!(model)
    st = termination_status(model)
    ps = primal_status(model)
    have = ps in (MOI.FEASIBLE_POINT, MOI.NEARLY_FEASIBLE_POINT)
    bound = try
        dual_objective_value(model) * scale
    catch
        NaN
    end
    pobj = have ? objective_value(model) * scale : NaN
    yv = have ? Dict(m => value(v) for (m, v) in y) : Dict{Vector{Int},Float64}()
    return yv, bound, pobj, st, ps, length(y)
end

"Pseudo-moment of a monomial (1 for the empty monomial, `missing` if not in the relaxation)."
function moment(rel::MomentRelaxation, m::Vector{Int})
    isempty(m) && return 1.0
    return get(rel.y, sort(m), missing)
end

"Order function: `high` for cliques containing a binary variable, `low` otherwise."
binary_clique_order(pop::POP; high = 2, low = 1) = c -> any(pop.isbin[v] for v in c) ? high : low

"""
Order function: `high` for cliques (of at most `maxsize` variables) that contain a binary or a variable
appearing in some constraint together with a binary; `low` otherwise.
"""
function adjacent_clique_order(pop::POP; high = 2, low = 1, maxsize = typemax(Int))
    S = Set(findall(pop.isbin))
    for p in [pop.ineqs; pop.eqs]
        sp = support(reduce_poly(p, pop.isbin))
        any(pop.isbin[v] for v in sp) && union!(S, sp)
    end
    return c -> (length(c) <= maxsize && any(in(S), c)) ? high : low
end

"""
Feasibility of recovered pseudo-moments: smallest eigenvalue of the moment blocks (first `nmom` blocks)
and of the localizing / PMI blocks (each normalized by max(1, largest |eigenvalue|)), and the largest
equality residual |L_y(h·m)|.
"""
function _moment_feasibility(yv, blocks, eqblocks, isbin, nmom)
    Ly(p::Poly, mult) = sum((c * (isempty(mm) ? 1.0 : get(yv, mm, 0.0)) for (m, c) in p.terms for mm in (monoprod(m, mult, isbin),)); init = 0.0)
    emin_m, emin_l = Inf, Inf
    for (k, (G, B)) in enumerate(blocks)
        ng, nb = size(G, 1), length(B)
        M = zeros(ng * nb, ng * nb)
        for a in 1:ng, b in 1:ng, i in 1:nb, j in 1:nb
            M[(a-1)*nb+i, (b-1)*nb+j] = Ly(G[a, b], monoprod(B[i], B[j], isbin))
        end
        ev = eigvals(Symmetric(M))
        r = ev[1] / max(1.0, maximum(abs, ev))
        k <= nmom ? (emin_m = min(emin_m, r)) : (emin_l = min(emin_l, r))
    end
    res = 0.0
    for (h, mults) in eqblocks, m in mults
        res = max(res, abs(Ly(h, m)))
    end
    return emin_m, emin_l, res
end

"Degree-2 pivot monomials of a set of degree-2 polynomials (Gaussian elimination with partial pivoting)."
function _pivot_monomials(E::Vector{Poly}; tol = 1e-9)
    isempty(E) && return Vector{Int}[]
    cols = sort!(unique([m for h in E for m in keys(h.terms) if length(m) == 2]))
    isempty(cols) && return Vector{Int}[]
    cidx = Dict(m => j for (j, m) in enumerate(cols))
    A = zeros(length(E), length(cols))
    for (i, h) in enumerate(E), (m, c) in h.terms
        length(m) == 2 && (A[i, cidx[m]] = c)
    end
    for i in axes(A, 1)                       # row scaling
        r = maximum(abs, A[i, :]; init = 0.0)
        r > 0 && (A[i, :] ./= r)
    end
    piv = Vector{Int}[]
    row = 1
    for j in axes(A, 2)
        row > size(A, 1) && break
        i = argmax(abs.(A[row:end, j])) + row - 1
        abs(A[i, j]) < tol && continue
        A[[row, i], :] = A[[i, row], :]
        for r in row+1:size(A, 1)
            A[r, :] .-= (A[r, j] / A[row, j]) .* A[row, :]
        end
        push!(piv, cols[j])
        row += 1
    end
    return piv
end
