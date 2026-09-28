# Determinant / minor relaxation of the semidefinite constraints.
#
# A symmetric or Hermitian matrix is PSD iff every principal minor is nonnegative, so enforcing
# only the minors up to order k is NECESSARY but not sufficient: the feasible set grows and the
# optimum remains a valid LOWER bound. Two references prescribe the two variants:
#
#   Molzahn & Hiskens, "Mixed SDP/SOCP Moment Relaxations of the OPF Problem", PowerTech 2015.
#   Its (15) is A_ii*A_kk >= |A_ik|^2 for k > i; its (16) keeps a full SDP constraint on the
#   FIRST-ORDER sub-block of each clique's moment matrix and applies cones only above it.
#   Reported 1.13x-18.70x faster than the SDP-based relaxation at equal iteration count.
#
#   Hijazi, Coffrin & Van Hentenryck, "Polynomial SDP Cuts for OPF", PSCC 2016. Writes det >= 0 as
#   explicit polynomial constraints for a general NLP solver. Its Remark 1 is load-bearing:
#   convexity of the region holds only if all LOWER-order determinants are imposed together.
#
# THE IDENTITY THAT MAKES THIS CHEAP. Enforcing every principal minor of order <= k on an index
# subset S is exactly "A[S,S] is PSD". So a determinant relaxation is fully described by a family
# of index subsets, the relaxed region is an intersection of small PSD cones -- convex by
# construction, no determinant polynomial required -- and for |S| = 2 that cone is a rotated
# second-order cone. Mosek therefore carries 3x3 and larger minors perfectly well, as small PSD
# cones rather than scalar determinants; a general NLP solver is NOT required to go beyond k = 2.
#
# ROW GROUPS, AND WHY THE REAL HIERARCHY NEEDS THEM. In the complex hierarchy a moment block at
# order 1 IS the Hermitian W, so a 2x2 principal minor is exactly the Jabr/SOCP constraint
# W_ii*W_jj >= |W_ij|^2. In the REAL rectangular hierarchy it is not: there
#
#     W_ij = L(e_i e_j + f_i f_j) + i * L(f_i e_j - e_i f_j)
#
# so W is a linear IMAGE of the real moment matrix M and a 2x2 minor of W is not a principal minor
# of M at all -- it is implied by the 4x4 principal submatrix on rows (e_i, f_i, e_j, f_j). Taking
# row-wise pairs in the real hierarchy therefore gives something strictly WEAKER than SOCP. The
# cover unit must be the BUS (its two rows), so k buses give a 2k x 2k cone. `row_groups` carries
# that; row-wise pairs remain available only as an explicitly labelled ablation.

"""
    minor_subsets(groups; kmax = 2, core = Int[]) -> Vector{Vector{Int}}

Index subsets whose principal submatrices replace one PSD block.

`groups` is a vector of row-index groups -- one entry per bus in the real hierarchy (its `e` and
`f` rows), one entry per row in the complex hierarchy. `core` is a set of row indices kept together
in one full PSD sub-block, which is PowerTech (16): the first-order part stays exact while the
higher-order part is relaxed. Every combination of up to `kmax` groups is emitted, so `kmax = 1`
gives the diagonal blocks, `kmax = 2` the pairwise ones, and so on.

Singletons are always emitted, so the diagonal is covered without separate constraints.
"""
function minor_subsets(groups::Vector{Vector{Int}}; kmax::Int = 2, core::Vector{Int} = Int[])
    subs = Vector{Int}[]
    isempty(core) || push!(subs, sort(unique(core)))
    ng = length(groups)
    for a in 1:ng
        push!(subs, sort(groups[a]))
        kmax >= 2 || continue
        for b in (a + 1):ng
            push!(subs, sort(vcat(groups[a], groups[b])))
            kmax >= 3 || continue
            for c in (b + 1):ng
                push!(subs, sort(vcat(groups[a], groups[b], groups[c])))
            end
        end
    end
    # a subset contained in another adds nothing: only the maximal ones need posting
    sort!(subs; by = length, rev = true)
    maximal = Vector{Int}[]
    for s in subs
        any(m -> issubset(s, m), maximal) || push!(maximal, s)
    end
    return maximal
end

"""
    post_minor_real!(model, M, S) -> ConstraintRef

`M[S,S] >= 0` for a real symmetric `M`, in the cheapest cone that expresses it.

The two-index case uses the sqrt(2) scaling `[M_ii, M_jj, sqrt(2) M_ij] in RotatedSecondOrderCone`
rather than the more obvious `[M_ii/2, M_jj, M_ij]`. Both encode the same minor, but this one makes
the Euclidean metric on the cone vector equal the Frobenius metric on `[M_ii M_ij; M_ij M_jj]`, so
the eigen-projection in `certify.jl`'s `project_dual` is literally the correct cone projection
rather than an approximation of it. Getting the scaling wrong degrades the bound silently.
"""
function post_minor_real!(model, M::AbstractMatrix, S::Vector{Int})
    k = length(S)
    if k == 1
        i = S[1]
        return @constraint(model, M[i, i] >= 0)
    elseif k == 2
        i, j = S[1], S[2]
        return @constraint(model, [M[i, i], M[j, j], sqrt(2) * M[i, j]] in RotatedSecondOrderCone())
    else
        return @constraint(model, Symmetric([M[i, j] for i in S, j in S]) in PSDCone())
    end
end

"""
    post_minor_hermitian!(model, A, B, S) -> ConstraintRef

`H[S,S] >= 0` for a Hermitian `H = A + iB` given its real and imaginary parts entry-wise.

For two indices this is a rotated cone directly on the four real quantities,
`2*A_ii*A_jj >= 2(A_ij^2 + B_ij^2)`, i.e. `A_ii*A_jj >= |H_ij|^2` -- the Jabr constraint, and
strictly cheaper than lifting a 2x2 Hermitian block to a real 4x4 PSD cone. Three or more indices
use the same real lift `[A -B; B A]` as the full blocks, which the certificate's de-lift depends on.
"""
function post_minor_hermitian!(model, A::AbstractMatrix, B::AbstractMatrix, S::Vector{Int})
    k = length(S)
    if k == 1
        i = S[1]
        return @constraint(model, A[i, i] >= 0)
    elseif k == 2
        i, j = S[1], S[2]
        return @constraint(model, [A[i, i], A[j, j], sqrt(2) * A[i, j], sqrt(2) * B[i, j]]
                           in RotatedSecondOrderCone())
    else
        X = [p <= k ? (q <= k ? A[S[p], S[q]] : -B[S[p], S[q - k]]) :
                      (q <= k ? B[S[p - k], S[q]] : A[S[p - k], S[q - k]]) for p in 1:2k, q in 1:2k]
        return @constraint(model, Symmetric(X) in PSDCone())
    end
end

"Row groups for a complex moment block: one row per group, since the block already IS W."
complex_row_groups(nb::Int) = [[i] for i in 1:nb]

"""
    real_bus_row_groups(B, varbus) -> Vector{Vector{Int}}

Row groups for a real moment block: the rows of one bus travel together, because a Jabr constraint
on buses (i, j) needs the 4x4 submatrix on (e_i, f_i, e_j, f_j) and not a 2x2 minor. `varbus[v]`
gives the bus owning POP variable `v`, or 0 for a variable that is not a voltage coordinate.
"""
function real_bus_row_groups(B::Vector{Vector{Int}}, varbus::Dict{Int,Int})
    bybus = Dict{Int,Vector{Int}}()
    other = Vector{Int}[]
    for (i, m) in enumerate(B)
        bs = unique([get(varbus, v, 0) for v in m])
        if length(m) == 1 && length(bs) == 1 && bs[1] != 0
            push!(get!(bybus, bs[1], Int[]), i)
        else
            push!(other, [i])
        end
    end
    return vcat([bybus[k] for k in sort(collect(keys(bybus)))], other)
end

"""
    basis_minor_subsets(basis; kmax = 2, core_degree = 1, adjacency = nothing) -> Vector{Vector{Int}}

Cover for a block whose rows are MONOMIALS rather than single variables, which is what orders 2 and
3 produce. This is the order-1 branch-pair rule generalised.

Two parts, following PowerTech (16) and then (15):

* every row of degree <= `core_degree` goes into ONE full PSD sub-block. At order 1 that is the
  whole block; at orders 2 and 3 it is the cheap, strong part that should stay exact.
* above it, two rows are paired when their variable supports OVERLAP or are joined by an edge of
  `adjacency`. Pairing densely instead would cost C(165,2) = 13530 cones for a single 9-bus
  order-3 block and C(364,2) = 66066 for a 12-bus one; support-overlap pairing is a small multiple
  of the clique's edge count, which is what makes order 3 reachable at all.

With single-variable rows and an `adjacency` of branches this reduces exactly to the branch-pair
cover that reproduces SOCWR, so the order-1 gate still applies unchanged.
"""
function basis_minor_subsets(basis::Vector{Vector{Int}}; kmax::Int = 2, core_degree::Int = 1,
        adjacency::Union{Nothing,Set{Tuple{Int,Int}}} = nothing)
    n = length(basis)
    n == 0 && return Vector{Int}[]
    sup = [Set(m) for m in basis]
    core = [i for i in 1:n if length(basis[i]) <= core_degree]
    rest = [i for i in 1:n if length(basis[i]) > core_degree]
    linked(i, j) = begin
        isempty(intersect(sup[i], sup[j])) || return true
        adjacency === nothing && return true
        for u in sup[i], v in sup[j]
            u == v && return true
            ((min(u, v), max(u, v)) in adjacency) && return true
        end
        return false
    end
    subs = Vector{Int}[]
    isempty(core) || push!(subs, sort(core))
    for i in rest
        push!(subs, [i])
    end
    if kmax >= 2
        for a in eachindex(rest)
            i = rest[a]
            for j in core
                linked(i, j) && push!(subs, sort([i, j]))
            end
            for b in (a + 1):length(rest)
                j = rest[b]
                linked(i, j) && push!(subs, sort([i, j]))
            end
        end
    end
    if kmax >= 3
        pairs = [s for s in subs if length(s) == 2]
        for p in pairs, i in rest
            i in p && continue
            (linked(i, p[1]) && linked(i, p[2])) || continue
            push!(subs, sort(vcat(p, i)))
        end
    end
    sort!(subs; by = length, rev = true)
    maximal = Vector{Int}[]
    for s in subs
        any(m -> issubset(s, m), maximal) || push!(maximal, s)
    end
    return maximal
end

"Number of scalar cone entries a cover implies, as a cheap cost proxy for the promotion schedule."
cover_cost(subs::Vector{Vector{Int}}) =
    sum(length(S) == 1 ? 1 : (length(S) == 2 ? 4 : (2 * length(S))^2) for S in subs; init = 0)
