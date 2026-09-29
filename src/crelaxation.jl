# Complex moment relaxation (Josz & Molzahn, (3.9) and the multi-ordered (3.31)).
#
#   rho_d = inf_y  L_y(f)   s.t.  y_{0,0} = 1,  M_{d-k_i}(g_i y) >= 0  (Hermitian PSD)
#
# with M_d(phi y)(alpha, beta) = sum_{gamma,delta} phi_{gamma,delta} y_{alpha+gamma, beta+delta}.
#
# HERMITIAN TO REAL. Mosek takes real symmetric cones, and a Hermitian H = A + iB (A symmetric,
# B skew-symmetric) is positive semidefinite exactly when the real 2k x 2k matrix
# [A -B; B A] is. That is the only conversion needed, and it is applied per block.
#
# T-INVARIANCE (the paper's Section 3.5) is the complex analogue of the parity splitting already
# used in the real hierarchy. The AC-OPF objective and constraints depend on z only through
# conj(z_i) z_j pairs, so they are invariant under z -> exp(i*theta) z, and there is an optimal
# moment vector with y_{alpha,beta} = 0 whenever |alpha| != |beta|. Every moment and localizing
# matrix then block-diagonalises by |alpha|, and at order 1 the blocks are exactly [1] and the
# n x n matrix W = v v^H. THAT IS THE CORRECTNESS GATE: complex order 1 with T-invariance is the
# standard SDP (Shor) relaxation of AC-OPF, so it must reproduce `sdp`/`moment1` to solver
# tolerance. The precondition is checked, never assumed.

using LinearAlgebra

"Cached (off-diagonal weight, dual sign) for the SOS identity; see the certificate assembly."
const _DUAL_CONVENTION = Ref{Union{Nothing,Tuple{Float64,Float64}}}(nothing)

struct ComplexMomentRelaxation
    status::Any
    bound::Float64
    primal_objective::Float64
    build_time::Float64
    solve_time::Float64
    cliques::Vector{Vector{Int}}
    psd_sizes::Vector{Int}          # sizes of the REAL blocks actually posted
    n_moments::Int
    t_invariant::Bool
    y::Dict{CMono,ComplexF64}       # the moment vector, for rank-one recovery
    info::Dict{String,Any}
end

"Canonical key for y_{alpha,beta}: the pair itself if alpha <= beta lexicographically, else swapped."
function _ckey(a::Vector{Int}, b::Vector{Int})
    return (a < b || a == b) ? ((a, b), false) : ((b, a), true)
end

"""
    solve_complex_moment_relaxation(pop; order, cliques, t_invariant, optimizer) -> ComplexMomentRelaxation

Moment form of the complex hierarchy. Returns MOSEK's bound; the rigorous certificate is a separate
step and is not applied here.
"""
function solve_complex_moment_relaxation(pop::CPOP; order::Int = 1,
        cliques::Union{Nothing,Vector{Vector{Int}}} = nothing,
        clique_orders::Union{Nothing,Vector{Int}} = nothing,
        t_invariant::Union{Nothing,Bool} = nothing,
        optimizer = mosek_optimizer(), silent::Bool = true, certify::Bool = false,
        psd_mode::Symbol = :full, minor_kmax::Int = 2, psd_threshold::Int = 0,
        minor_core_degree::Int = 1,
        adjacency::Union{Nothing,Set{Tuple{Int,Int}}} = nothing,
        minor_cover = nothing, retain = nothing, warm_start = nothing,
        solver_params = Dict{String,Any}())
    psd_mode === :nlp && certify &&
        error("psd_mode = :nlp posts determinant inequalities, not PSD cones, so the certificate's " *
              "Gram recovery does not apply; solve with certify = false and verify the returned " *
              "moments against the true blocks instead")
    t0 = time()
    n = cnvars(pop)
    cls = cliques === nothing ? [collect(1:n)] : cliques
    # the paper's multi-ordered relaxation (3.31): one order per clique, so second-order blocks
    # are paid for only where they are wanted
    ords = clique_orders === nothing ? fill(order, length(cls)) : clique_orders
    length(ords) == length(cls) || error("clique_orders must have one entry per clique")

    # T-invariance precondition: every term of every polynomial must have |alpha| == |beta|
    tinv_ok = all(length(a) == length(b)
                  for p in vcat([pop.obj], pop.ineqs, pop.eqs) for (a, b) in keys(p.terms))
    tinv = t_invariant === nothing ? tinv_ok : (t_invariant && tinv_ok)
    t_invariant === true && !tinv_ok &&
        error("T-invariance requested but some term has |alpha| != |beta|; the hierarchy is " *
              "still valid without it, but the block splitting is not")

    model = Model(optimizer)
    silent && set_silent(model)
    for (k, v) in solver_params; set_attribute(model, k, v); end

    # one real pair per canonical moment; a diagonal moment y_{alpha,alpha} is real
    yre = Dict{CMono,VariableRef}()
    # a diagonal moment y_{alpha,alpha} is real, so it has no imaginary variable
    yim = Dict{CMono,Union{Nothing,VariableRef}}()
    function ypair(a::Vector{Int}, b::Vector{Int})
        (key, swapped) = _ckey(a, b)
        if !haskey(yre, key)
            yre[key] = @variable(model)
            yim[key] = key[1] == key[2] ? nothing : @variable(model)
        end
        re = yre[key]
        im = yim[key]
        im === nothing && return (re, nothing, false)
        return (re, im, swapped)          # y_{b,a} = conj(y_{a,b}) when swapped
    end
    "Real and imaginary parts of y_{a,b} as affine expressions."
    function ymom(a::Vector{Int}, b::Vector{Int})
        (re, im, swapped) = ypair(a, b)
        # AffExpr(::VariableRef) is not a constructor; build the affine expression explicitly
        rex = one(Float64) * re
        im === nothing && return (rex, AffExpr(0.0))
        imx = one(Float64) * im
        return (rex, swapped ? -imx : imx)
    end

    one_poly = CPoly(1)
    psd_sizes = Int[]
    blockmeta = Tuple{CPoly,Vector{Vector{Int}}}[]   # (weight, basis) per posted block
    # (A, C, basis) per hermitian_psd! CALL, i.e. per moment block rather than per posted cone.
    # The lazy minor loop needs these to rebuild each block's Hermitian matrix at the current
    # solution, score its violation and post further minors into the same model; without them the
    # matrices built by `loc_matrix` are local to the build loop and thrown away.
    blocks_ab = Tuple{Matrix,Matrix,Vector{Vector{Int}}}[]
    psd_refs = Any[]                 # (constraint, k) for Gram recovery in the certificate
    """
    Post a Hermitian PSD constraint from its entry-wise (real, imag) affine parts.

    With `psd_mode = :minors` a block wider than `psd_threshold` is replaced by PSD constraints on
    the principal submatrices of up to `minor_kmax` rows -- a valid relaxation, since a PSD matrix
    has PSD principal submatrices. `blockmeta` is pushed HERE rather than by the callers, because
    one call now emits many constraints and the three registries would otherwise drift apart.
    """
    function hermitian_psd!(A, B, weight = one_poly, basis = Vector{Int}[])
        k = size(A, 1)
        k == 0 && return
        push!(blocks_ab, (A, B, basis))
        bidx = length(blocks_ab)
        # An explicit cover overrides the generated one: this is what lets the NLP and the conic
        # model be built over the SAME posted subsets, which is the only way the Ipopt-vs-Mosek
        # comparison means anything.
        cover(kmx) = minor_cover === nothing ?
            basis_minor_subsets(basis; kmax = kmx, core_degree = minor_core_degree,
                                adjacency = adjacency) :
            minor_cover(bidx, basis)
        if psd_mode === :nlp && !isempty(basis)
            # Variant 2: no PSD cone at all, so EVERY block becomes principal minors -- including
            # the core block that :minors keeps exact. The model is therefore weaker than :minors,
            # which is the price of reaching order 3 with a solver that has no conic support.
            posted = 0
            for S in cover(minor_kmax)
                post_minor_hermitian_nlp!(model, A, B, S; kmax = minor_kmax) > 0 || continue
                push!(psd_sizes, length(S))
                push!(blockmeta, (weight, [basis[i] for i in S]))
                posted += 1
            end
            posted > 0 && return
        end
        if psd_mode === :minors && k > psd_threshold && !isempty(basis)
            posted = 0
            # With T-invariance each block holds monomials of a single |alpha|, so
            # `minor_core_degree = 1` keeps the |alpha| = 1 block (the Hermitian W) EXACT and
            # relaxes only the higher-order blocks -- PowerTech (16) generalised, and exactly what
            # makes order 3 reachable without giving up the order-1 strength that is nearly free.
            for S in cover(minor_kmax)
                cref = post_minor_hermitian!(model, A, B, S)
                push!(psd_sizes, length(S) <= 2 ? length(S) : 2 * length(S))
                push!(psd_refs, (cref, length(S)))
                push!(blockmeta, (weight, [basis[i] for i in S]))
                posted += 1
            end
            posted > 0 && return
            # nothing survived the adjacency filter: fall through to the full block rather than
            # silently dropping the constraint entirely
        end
        if k == 1
            c1 = @constraint(model, A[1, 1] >= 0)
            push!(psd_sizes, 1); push!(psd_refs, (c1, 1)); push!(blockmeta, (weight, basis))
            return
        end
        X = [i <= k ? (j <= k ? A[i, j] : -B[i, j - k]) :
                      (j <= k ? B[i - k, j] : A[i - k, j - k]) for i in 1:2k, j in 1:2k]
        cc = @constraint(model, Symmetric(X) in PSDCone())
        push!(psd_sizes, 2k); push!(psd_refs, (cc, k)); push!(blockmeta, (weight, basis))
    end

    "Blocks of a basis: split by |alpha| when T-invariant, else one block."
    function split_basis(B::Vector{Vector{Int}})
        tinv || return [B]
        d = Dict{Int,Vector{Vector{Int}}}()
        for m in B; push!(get!(d, length(m), Vector{Int}[]), m); end
        return [d[k] for k in sort(collect(keys(d)))]
    end

    "Localizing matrix of weight phi over basis B, as (Re, Im) matrices."
    function loc_matrix(phi::CPoly, B::Vector{Vector{Int}})
        k = length(B)
        A = Matrix{AffExpr}(undef, k, k); C = Matrix{AffExpr}(undef, k, k)
        for i in 1:k, j in 1:k
            re, im = AffExpr(0.0), AffExpr(0.0)
            for ((g, d), c) in phi.terms
                r, m = ymom(sort!(vcat(B[i], g)), sort!(vcat(B[j], d)))
                add_to_expression!(re, real(c), r); add_to_expression!(re, -imag(c), m)
                add_to_expression!(im, real(c), m); add_to_expression!(im, imag(c), r)
            end
            A[i, j] = re; C[i, j] = im
        end
        return A, C
    end

    nblocks = 0
    for (ci, cl) in enumerate(cls)
        B = cmonomial_basis(cl, ords[ci])
        for Bp in split_basis(B)
            A, C = loc_matrix(one_poly, Bp)
            hermitian_psd!(A, C, one_poly, Bp); nblocks += 1
        end
    end
    skipped = Dict{String,Int}()
    # A monomial y_{alpha,beta} is available whenever alpha and beta together sit inside SOME
    # clique -- that is all a LINEAR functional of the moments needs. Requiring a constraint's
    # whole support to lie in one clique is the right condition for a localizing MATRIX and quite
    # wrong for a scalar multiplier, and getting that backwards silently dropped every injection
    # constraint spanning two cliques: case14's order-1 bound collapsed from 1147.59 to -12701.95.
    produced(g) = all(any(issubset(vcat(a, b), c) for c in cls) for (a, b) in keys(g.terms))
    for (g, tag) in zip(pop.ineqs, pop.ineq_tags)
        cis_constant(g) && continue
        sup = csupport(g)
        k = findfirst(c -> issubset(sup, c), cls)
        dd = k === nothing ? -1 : ords[k] - cld(cdegree(g), 2)
        if k !== nothing && dd >= 1
            B = cmonomial_basis(cls[k], dd)
            for Bp in split_basis(B)
                A, C = loc_matrix(g, Bp)
                hermitian_psd!(A, C, g, Bp); nblocks += 1
            end
        elseif produced(g)
            # scalar multiplier: impose L(g) >= 0 directly
            A, C = loc_matrix(g, [Int[]])
            hermitian_psd!(A, C, g, [Int[]]); nblocks += 1
        else
            skipped[tag] = get(skipped, tag, 0) + 1
        end
    end
    for (h, tag) in zip(pop.eqs, pop.eq_tags)
        cis_constant(h) && continue
        sup = csupport(h)
        k = findfirst(c -> issubset(sup, c), cls)
        dd = k === nothing ? -1 : 2 * ords[k] - cdegree(h)
        if dd < 0
            if produced(h)
                # L(h) = 0 alone, the scalar-multiplier analogue for equalities
                re, im = AffExpr(0.0), AffExpr(0.0)
                for ((g, d), c) in h.terms
                    r, i2 = ymom(g, d)
                    add_to_expression!(re, real(c), r); add_to_expression!(re, -imag(c), i2)
                    add_to_expression!(im, real(c), i2); add_to_expression!(im, imag(c), r)
                end
                @constraint(model, re == 0); @constraint(model, im == 0)
            else
                skipped[tag] = get(skipped, tag, 0) + 1
            end
            continue
        end
        for m in cmonomial_basis(cls[k], dd ÷ 2), m2 in cmonomial_basis(cls[k], dd ÷ 2)
            re, im = AffExpr(0.0), AffExpr(0.0)
            for ((g, d), c) in h.terms
                r, i2 = ymom(sort!(vcat(m, g)), sort!(vcat(m2, d)))
                add_to_expression!(re, real(c), r); add_to_expression!(re, -imag(c), i2)
                add_to_expression!(im, real(c), i2); add_to_expression!(im, imag(c), r)
            end
            @constraint(model, re == 0)
            @constraint(model, im == 0)
        end
    end

    # normalisation y_{0,0} = 1
    r0, _ = ymom(Int[], Int[])
    @constraint(model, r0 == 1)

    "Real part of L(p) for a real-valued p (the imaginary part vanishes by Hermitian symmetry)."
    function lin_real(p::CPoly)
        e = AffExpr(0.0)
        for ((a, b), c) in p.terms
            re, im = ymom(a, b)
            add_to_expression!(e, real(c), re)
            add_to_expression!(e, -imag(c), im)
        end
        return e
    end

    # auxiliary REAL variables: generator outputs at buses carrying more than one generator. They
    # are ordinary JuMP variables and never enter a moment matrix, which is the point -- they add
    # no monomials and so cost nothing in the hierarchy's size.
    auxv = VariableRef[]
    for (nm, lo, hi) in pop.aux
        v = @variable(model, base_name = nm)
        isfinite(lo) && set_lower_bound(v, lo)
        isfinite(hi) && set_upper_bound(v, hi)
        push!(auxv, v)
    end
    for ((coefs, p, sense), tag) in zip(pop.mixed, pop.mixed_tags)
        e = lin_real(p)
        for (i, c) in coefs; add_to_expression!(e, c, auxv[i]); end
        sense === :eq ? @constraint(model, e == 0) : @constraint(model, e >= 0)
    end
    # convex quadratic epigraphs on the auxiliaries, as rotated cones: 2*(t/(2a))*1 >= x^2
    for (ti, a, xi) in pop.aux_quad
        a > 0 || continue
        if psd_mode === :nlp
            # 2*u*v >= w^2 with u = t/(2a), v = 1 is just t/a >= x^2, already smooth
            @constraint(model, auxv[ti] / a - auxv[xi] * auxv[xi] >= 0)
        else
            @constraint(model, [auxv[ti] / (2a), 1.0, auxv[xi]] in RotatedSecondOrderCone())
        end
    end
    # second-order cones on the moment image (thermal limits); valid at every order
    for ((a, xs), tag) in zip(pop.socs, pop.soc_tags)
        if psd_mode === :nlp
            # Same convexity trap as the 2x2 minor: `a^2 - sum(x^2) >= 0` describes the cone with a
            # nonconcave function. `a >= ||x||` is convex as written, smoothed so the set only grows.
            ea = lin_real(a)
            ex = [lin_real(x) for x in xs]
            ce = NLP_CONE_EPS[]
            sq = sum(e * e for e in ex)
            @constraint(model, ea - (ce > 0 ? sqrt(sq + ce) - sqrt(ce) : sqrt(sq)) >= 0)
        else
            @constraint(model, vcat(lin_real(a), [lin_real(x) for x in xs]) in SecondOrderCone())
        end
    end

    # Normalise the objective. Cost coefficients are O(1e3) while every moment is O(1), and handing
    # MOSEK that spread is what made the first gate run return SLOW_PROGRESS. The real hierarchy
    # does the same thing via `scale`.
    oscale = max(maximum((abs(c) for c in values(pop.obj.terms)); init = 1.0),
                 maximum((abs(c) for c in values(get(pop.meta, "obj_aux", Dict{Int,Float64}())));
                         init = 1.0))
    oscale > 0 || (oscale = 1.0)
    obj = AffExpr(0.0)
    for (i, c) in get(pop.meta, "obj_aux", Dict{Int,Float64}())
        add_to_expression!(obj, c / oscale, auxv[i])
    end
    for ((a, b), c) in pop.obj.terms
        re, im = ymom(a, b)
        add_to_expression!(obj, real(c) / oscale, re); add_to_expression!(obj, -imag(c) / oscale, im)
    end
    @objective(model, Min, obj)
    # Ipopt is strongly start-dependent, and cold-starting it on thousands of moment variables was
    # worth 3000 iterations of nothing (ITERATION_LIMIT at 2e7 on case14 order 1). `warm_start` is a
    # moment dictionary from a previous solve -- typically the conic solve of the SAME cover, which
    # the lazy loop computes anyway for its cross-check, so the start costs nothing extra.
    if warm_start !== nothing
        # Sanitise. A start value that is NaN/Inf makes Ipopt reject the model outright with
        # INVALID_MODEL, and the source here is another solver's output -- an unconverged or weakly
        # constrained conic solve can hand back enormous moments. Skipping the bad entries leaves
        # those variables at their default start rather than poisoning the whole solve.
        nbad = 0
        for (key, rv) in yre
            v = get(warm_start, key, nothing)
            v === nothing && continue
            re, im = real(v), imag(v)
            if !isfinite(re) || !isfinite(im) || abs(re) > 1e8 || abs(im) > 1e8
                nbad += 1
                continue
            end
            set_start_value(rv, re)
            iv = yim[key]
            iv === nothing || set_start_value(iv, im)
        end
        nbad > 0 && @warn "warm start: skipped $nbad non-finite or oversized moment values"
    end
    if retain !== nothing
        retain[] = (model = model, yre = yre, yim = yim, blocks_ab = blocks_ab,
                    oscale = oscale, auxv = auxv)
    end
    build_time = time() - t0
    optimize!(model)
    st = termination_status(model)
    bound = oscale * (try
        dual_objective_value(model)
    catch
        try objective_value(model) catch; NaN end
    end)
    pobj = oscale * (try objective_value(model) catch; NaN end)
    # ---- certificate data -------------------------------------------------------------------
    # The residual of the SOS identity is the DUAL-FEASIBILITY residual of the moment form: each
    # moment y is a free variable, so at an exact solution its reduced cost is zero, and what
    # MOSEK actually leaves is what the certificate has to pay for. Computed by one pass over
    # every constraint, so every channel -- PSD blocks, equalities, cones, the auxiliary
    # variables -- is covered without reconstructing any of them by hand.
    certdata = Dict{String,Any}()
    if certify && has_duals(model)
        rc = Dict{VariableRef,Float64}()
        addrc!(v, c) = (rc[v] = get(rc, v, 0.0) + c)
        # The real-part variable of y_{a,b} is fed by BOTH y_{a,b} and its conjugate partner
        # y_{b,a}, since Re(y_{b,a}) = Re(y_{a,b}). Counting only one of them left a residual the
        # size of the objective's own coefficients.
        function objcoef_re(key::CMono)
            (a, b) = key
            t = 0.0
            for ((p1, p2), cc) in pop.obj.terms
                ((p1, p2) == (a, b) || (p1, p2) == (b, a)) && (t += real(cc))
            end
            return t / oscale
        end
        function walk!(fexpr, dv)
            for (t, d) in zip(fexpr, dv)
                for (v, c) in linear_terms_of(t)
                    addrc!(v, d * c)
                end
            end
        end
        # A PSD constraint posted as `Symmetric(X) in PSDCone()` has a MATRIX function and a
        # matrix dual, not vectors. Sending it down the vector path silently contributed nothing,
        # which left a residual of 21.2 where it should be ~1e-8 -- the whole semidefinite channel
        # was missing. Matrices are paired entry by entry; `offdiag` is the weight for i != j and
        # is settled by measurement below rather than by assuming a convention.
        function sweep!(offdiag, sgn)
            empty!(rc)
            for (F, S) in list_of_constraint_types(model)
                for c in all_constraints(model, F, S)
                    d = try dual(c) catch; continue end
                    fn = constraint_object(c).func
                    if fn isa AbstractVector && d isa AbstractMatrix
                        # PSD: JuMP hands back the VECTORISED UPPER TRIANGLE as the function and a
                        # Symmetric MATRIX as the dual. zip-ing the two pairs them in different
                        # orders -- column-major over the matrix against triangle order over the
                        # vector -- which is what scrambled the identity and left a residual of
                        # 21.2. Walk the triangle explicitly instead, weighting off-diagonals by
                        # `offdiag` because <D, X> counts each of them twice.
                        kdim = size(d, 1)
                        t = 0
                        for j in 1:kdim, i in 1:j
                            t += 1
                            t <= length(fn) || break
                            w = i == j ? 1.0 : offdiag
                            for (v, cf) in linear_terms_of(fn[t])
                                addrc!(v, sgn * w * d[i, j] * cf)
                            end
                        end
                    elseif fn isa AbstractVector
                        walk!(fn, sgn .* d)
                    else
                        walk!([fn], [sgn * d])
                    end
                end
            end
        end
        # the correct off-diagonal weight is the one that makes the identity hold
        # JuMP's dual sign convention varies with the cone and the problem sense, and the
        # off-diagonal weight of a matrix pairing is a convention too. Rather than assume either,
        # try all four and keep whichever actually makes the identity hold -- the residual is a
        # sharp test, since the right combination gives ~1e-8 and every wrong one gives O(1).
        # The convention is a property of JuMP/MOI, not of this model, so it is determined once
        # and cached. Re-deriving it means four full passes over every constraint, which on a
        # 200-bus model is the dominant cost of the certificate.
        cand = _DUAL_CONVENTION[] === nothing ? ((1.0, 1.0), (1.0, -1.0), (2.0, 1.0), (2.0, -1.0)) :
               (_DUAL_CONVENTION[]::Tuple{Float64,Float64},)
        best_off, best_sgn, best_err = 1.0, 1.0, Inf
        for (off, sg) in cand
            sweep!(off, sg)
            err = 0.0
            for (key, rv) in yre
                fc = objcoef_re(key)
                err = max(err, abs(fc - get(rc, rv, 0.0)))
            end
            err < best_err && ((best_off, best_sgn, best_err) = (off, sg, err))
        end
        # only trust the cache when it actually reproduced the identity
        best_err < 1e-6 ? (_DUAL_CONVENTION[] = (best_off, best_sgn)) :
                          (_DUAL_CONVENTION[] = nothing)
        sweep!(best_off, best_sgn)
        certdata["psd_offdiag_weight"] = best_off
        certdata["dual_sign"] = best_sgn
        certdata["identity_err"] = best_err
        # r for the real and imaginary part of each moment; the identity row is the real part
        rows = CMono[]; fvec = Float64[]; resid = Float64[]
        for (key, rv) in yre
            push!(rows, key)
            fc = objcoef_re(key)
            push!(fvec, fc)
            push!(resid, fc - get(rc, rv, 0.0))
        end
        blocksH = Matrix{ComplexF64}[]
        rhos = Float64[]
        absmax = [isfinite(a) ? a : 1.0 for a in pop.absmax]
        for (bi, (cref, k)) in enumerate(psd_refs)
            X = zeros(ComplexF64, k, k)
            try
                D = dual(cref)
                if k == 1
                    X[1, 1] = ComplexF64(D isa Number ? D : D[1])
                else
                    M = Matrix(D)                      # 2k x 2k real dual of the lift
                    # adjoint of the lift: the structure-preserving projection back to Hermitian
                    D11 = M[1:k, 1:k]; D22 = M[k+1:2k, k+1:2k]
                    D12 = M[1:k, k+1:2k]; D21 = M[k+1:2k, 1:k]
                    X = (D11 .+ D22) ./ 2 .+ im .* ((D21 .- D12) ./ 2)
                end
            catch
            end
            push!(blocksH, X)
            g, B = bi <= length(blockmeta) ? blockmeta[bi] : (CPoly(1), [Int[]])
            # rho >= max over the box of ||p(z)||^2 * g(z)
            gmax = sum((abs(c) * prod((absmax[v] for v in vcat(a, b)); init = 1.0)
                        for ((a, b), c) in g.terms); init = 0.0)
            psq = sum((prod((absmax[v] for v in m); init = 1.0)^2 for m in B); init = 0.0)
            push!(rhos, gmax * psq)
        end
        certdata["rows"] = rows; certdata["f"] = fvec; certdata["resid"] = resid
        certdata["blocks"] = blocksH; certdata["rho"] = rhos
        certdata["scale"] = oscale
        certdata["lambda"] = try objective_value(model) catch; NaN end
    end

    yval = Dict{CMono,ComplexF64}()
    if has_values(model)
        for (key, rv) in yre
            iv = yim[key]
            yval[key] = complex(value(rv), iv === nothing ? 0.0 : value(iv))
        end
    end

    info = Dict{String,Any}("skipped" => skipped, "n_blocks" => nblocks,
        "t_invariant_possible" => tinv_ok, "obj_scale" => oscale,
        "clique_orders" => ords, "psd_mode" => string(psd_mode),
        "minor_kmax" => minor_kmax,
        "minor_core_degree" => minor_core_degree)
    isempty(certdata) || (info["cert_data"] = certdata)
    return ComplexMomentRelaxation(st, bound, pobj, build_time, solve_time(model), cls,
        psd_sizes, length(yre), tinv, yval, info)
end

"""
    cmoment(y, a, b) -> ComplexF64

Read y_{alpha,beta} from a solved moment vector, honouring y_{beta,alpha} = conj(y_{alpha,beta}).
"""
function cmoment(y::Dict{CMono,ComplexF64}, a::Vector{Int}, b::Vector{Int})
    (key, swapped) = _ckey(a, b)
    v = get(y, key, ComplexF64(0))
    return swapped ? conj(v) : v
end

"""
    clique_W(y, clique) -> Matrix{ComplexF64}

The Hermitian W = v v^H block of a clique, W[i,j] = y_{(c_i),(c_j)}. Under T-invariance this is
the |alpha| = 1 diagonal block of the clique's moment matrix, and at order 1 it IS the standard
SDP relaxation's W restricted to the clique.
"""
function clique_W(y::Dict{CMono,ComplexF64}, clique::Vector{Int})
    k = length(clique)
    W = zeros(ComplexF64, k, k)
    for i in 1:k, j in 1:k
        W[i, j] = cmoment(y, [clique[i]], [clique[j]])
    end
    return W
end

"""
    rank_one_point(W) -> (v, lambda2_over_lambda1)

Closest rank-one Hermitian factor of W: the leading eigenvector scaled by sqrt(lambda_1). The
eigenvalue ratio is returned alongside because it is the standard exactness indicator -- a ratio
near zero means W is essentially rank one and the relaxation is tight on that clique.
"""
function rank_one_point(W::Matrix{ComplexF64})
    k = size(W, 1)
    k == 0 && return (ComplexF64[], 0.0)
    E = eigen(Hermitian(W))
    lam = E.values
    i = argmax(lam)
    v = E.vectors[:, i] * sqrt(max(lam[i], 0.0))
    ratio = length(lam) > 1 ? (sort(lam; rev = true)[2] / max(lam[i], 1e-300)) : 0.0
    return (v, ratio)
end

"Linear (variable, coefficient) pairs of an affine expression, as a plain iterator."
linear_terms_of(e::AffExpr) = ((v, c) for (v, c) in e.terms)
linear_terms_of(v::VariableRef) = ((v, 1.0),)
linear_terms_of(x) = ()
