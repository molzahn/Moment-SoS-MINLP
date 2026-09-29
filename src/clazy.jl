# Lazy determinant (minor) relaxation.
#
# The minor relaxation's FEASIBLE SET is convex: enforcing every principal minor of order <= k on a
# subset S is exactly `W[S,S] >= 0`, an intersection of convex cones. A local minimiser of a LINEAR
# objective over a convex set is global, so a CONVERGED Ipopt solution is a valid lower bound --
# provided the constraints are written as convex FUNCTIONS, which is what
# `post_minor_hermitian_nlp!` now does for the 1x1/2x2 base.
#
# 3x3 Hermitian PSD is not SOC-representable, so a lazily added 3x3 determinant is unavoidably a
# nonconvex description of a convex constraint. Adding only the most-violated ones keeps that set
# small. Crucially, because each subset's 1x1 and 2x2 minors are already posted, adding just its
# 3x3 determinant is EXACTLY `W[S,S] >= 0` -- so every round is a relaxation of the SDP, the
# converged values are monotone non-decreasing, and the loop may stop anywhere without invalidating
# the bound it has.

"""
Candidate principal submatrices of a Hermitian block, scored by their OWN violation.

The eigenvector for `lambda_min` names the rows carrying the violation, so candidates are drawn from
its `width` largest-magnitude entries. Each candidate subset is then scored by `-lambda_min(H[S,S])`
-- its own violation, not the whole block's -- because that is what posting it actually removes.

`exclude` lists subsets already posted. Scoring per subset and enumerating within the eigenvector's
support is what lets the loop keep making progress after the single worst triple is already in: the
first version returned only that one triple and stalled with `added = 0` while the block was still
violated by 6.4e-05.
"""
function candidate_minors(H::AbstractMatrix{ComplexF64}, k::Int;
        width::Int = 6, exclude = Vector{Int}[], tol::Float64 = 0.0)
    n = size(H, 1)
    (n < k || k < 2) && return Tuple{Float64,Vector{Int}}[]
    E = eigen(Hermitian(Matrix(H)))
    E.values[1] >= -tol && return Tuple{Float64,Vector{Int}}[]
    v = abs.(@view E.vectors[:, 1])
    pool = sort(collect(partialsortperm(v, 1:min(width, n); rev = true)))
    out = Tuple{Float64,Vector{Int}}[]
    for S in combinations_k(pool, k)
        S in exclude && continue
        lam = eigvals!(Hermitian(Matrix(H[S, S])), 1:1)[1]
        lam < -tol && push!(out, (-lam, S))
    end
    sort!(out; by = t -> -t[1])
    return out
end

"All sorted k-subsets of `xs` (k is 3 here, so the count stays small)."
function combinations_k(xs::Vector{Int}, k::Int)
    out = Vector{Int}[]
    n = length(xs)
    k > n && return out
    idx = collect(1:k)
    while true
        push!(out, [xs[i] for i in idx])
        i = k
        while i >= 1 && idx[i] == n - k + i
            i -= 1
        end
        i == 0 && break
        idx[i] += 1
        for j in (i + 1):k
            idx[j] = idx[j - 1] + 1
        end
    end
    return out
end

"""
    solve_complex_lazy_minors(pop; cliques, clique_orders, adjacency, ...) -> NamedTuple

Iteratively strengthen a 2x2-minor relaxation by posting the most violated 3x3 determinants, solving
each round with Ipopt.

`crosscheck` re-solves the IDENTICAL cover conically with Mosek each round and compares. It also
supplies Ipopt's warm start, which costs nothing extra and matters a lot -- cold-started Ipopt spent
3000 iterations getting nowhere on case14. After `crosscheck_streak` consecutive agreeing rounds the
conic solve is dropped and Ipopt runs alone, which is the point of the exercise.

Returns the per-round history rather than only the final value, so monotonicity is inspectable.
"""
function solve_complex_lazy_minors(pop::CPOP;
        cliques::Vector{Vector{Int}}, clique_orders::Vector{Int},
        adjacency::Union{Nothing,Set{Tuple{Int,Int}}} = nothing,
        subset_size::Int = 3, n_add::Int = 8, max_rounds::Int = 12,
        tol::Float64 = 1e-6, max_seconds::Float64 = 1800.0, round_seconds::Float64 = 600.0,
        crosscheck::Bool = true, crosscheck_streak::Int = 3, verbose::Bool = true)
    t0 = time()
    c0 = Float64(get(pop.meta, "const_cost", 0.0))
    # block index -> posted subsets. Populated on the first build from the 2x2 cover, then grown.
    cover = Dict{Int,Vector{Vector{Int}}}()
    coverfn = (bidx, basis) -> get!(cover, bidx) do
        basis_minor_subsets(basis; kmax = 2, core_degree = 0, adjacency = adjacency)
    end
    # minor_kmax = subset_size throughout: the COVER decides strength (which subsets exist), while
    # kmax only has to be large enough that a posted 3-subset gets its 3x3 determinant and not just
    # the 2x2s it already had.
    common = (cliques = cliques, clique_orders = clique_orders, certify = false,
              psd_mode = :nlp, minor_kmax = subset_size, minor_core_degree = 0,
              adjacency = adjacency, minor_cover = coverfn)
    rounds = Any[]
    best = -Inf
    warm = nothing
    agree = 0
    checking = crosscheck
    for it in 1:max_rounds
        cbound, cstat = NaN, "-"
        if checking
            rc = solve_complex_moment_relaxation(pop; common..., psd_mode = :minors,
                    optimizer = mosek_optimizer(),
                    solver_params = Dict{String,Any}("MSK_DPAR_OPTIMIZER_MAX_TIME" => round_seconds))
            cbound, cstat = rc.bound + c0, string(rc.status)
            isempty(rc.y) || (warm = rc.y)
        end
        ref = Ref{Any}(nothing)
        r = solve_complex_moment_relaxation(pop; common..., optimizer = ipopt_optimizer(),
                retain = ref, warm_start = warm,
                solver_params = Dict{String,Any}("max_iter" => 5000, "tol" => 1e-8,
                                                 "mu_strategy" => "adaptive",
                                                 "max_cpu_time" => round_seconds))
        ib = r.bound + c0
        ipr = r.primal_objective + c0
        converged = r.status in (MOI.LOCALLY_SOLVED, MOI.ALMOST_LOCALLY_SOLVED, MOI.OPTIMAL)
        # Only a CONVERGED round carries the local-minimum-of-a-convex-set argument, so only a
        # converged round may raise the reported bound.
        converged && isfinite(ib) && ib > best && (best = ib)
        isempty(r.y) || (warm = r.y)
        d = (checking && isfinite(cbound) && isfinite(ib)) ? abs(ib - cbound) : NaN
        if checking && isfinite(d)
            agree = d < 1e-4 * max(1.0, abs(cbound)) ? agree + 1 : 0
            if agree >= crosscheck_streak
                checking = false
                verbose && println("  cross-check agreed $agree rounds running; Ipopt alone from here")
            end
        end
        # score every block at the Ipopt point and add the most violated subsets
        st = ref[]
        worst, added = 0.0, 0
        if st !== nothing && converged
            cand = Tuple{Float64,Int,Vector{Int}}[]
            for (bi, (A, C, _)) in enumerate(st.blocks_ab)
                size(A, 1) >= subset_size || continue
                H = try complex.(value.(A), value.(C)) catch; continue end
                posted = get(cover, bi, Vector{Int}[])
                cs = candidate_minors(H, subset_size; exclude = posted, tol = tol)
                isempty(cs) || (worst = max(worst, cs[1][1]))
                for (v, S) in cs
                    push!(cand, (v, bi, S))
                end
            end
            sort!(cand; by = t -> -t[1])
            for (_, bi, rows) in first(cand, n_add)
                haskey(cover, bi) || continue
                rows in cover[bi] && continue
                push!(cover[bi], rows); added += 1
            end
        end
        push!(rounds, (iter = it, ipopt = ib, ipopt_primal = ipr, mosek = cbound, diff = d, best = best,
                       status = string(r.status), worst_violation = worst, added = added,
                       n_subsets = sum(length, values(cover); init = 0), time = time() - t0))
        verbose && @printf("  r%-2d ipopt d %13.4f p %13.4f %-22s mosek %13.4f  |d| %8.1e  worst %9.2e  +%d subs (%d)  %5.1fs\n",
            it, ib, ipr, string(r.status), cbound, d, worst, added,
            sum(length, values(cover); init = 0), time() - t0)
        flush(stdout)
        (worst <= tol || added == 0) && break
        time() - t0 > max_seconds && break
    end
    return (bound = best, rounds = rounds, cover = cover, time = time() - t0)
end
