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
        t_invariant::Union{Nothing,Bool} = nothing,
        optimizer = mosek_optimizer(), silent::Bool = true,
        solver_params = Dict{String,Any}())
    t0 = time()
    n = cnvars(pop)
    cls = cliques === nothing ? [collect(1:n)] : cliques

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

    psd_sizes = Int[]
    "Post a Hermitian PSD constraint given entry-wise (real, imag) affine expressions."
    function hermitian_psd!(A, B)
        k = size(A, 1)
        k == 0 && return
        if k == 1
            @constraint(model, A[1, 1] >= 0)
            push!(psd_sizes, 1)
            return
        end
        X = [i <= k ? (j <= k ? A[i, j] : -B[i, j - k]) :
                      (j <= k ? B[i - k, j] : A[i - k, j - k]) for i in 1:2k, j in 1:2k]
        @constraint(model, Symmetric(X) in PSDCone())
        push!(psd_sizes, 2k)
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

    one_poly = CPoly(1)
    nblocks = 0
    for cl in cls
        B = cmonomial_basis(cl, order)
        for Bp in split_basis(B)
            A, C = loc_matrix(one_poly, Bp)
            hermitian_psd!(A, C); nblocks += 1
        end
    end
    skipped = Dict{String,Int}()
    for (g, tag) in zip(pop.ineqs, pop.ineq_tags)
        cis_constant(g) && continue
        sup = csupport(g)
        k = findfirst(c -> issubset(sup, c), cls)
        dd = order - cld(cdegree(g), 2)
        if k === nothing || dd < 0
            skipped[tag] = get(skipped, tag, 0) + 1
            continue
        end
        B = cmonomial_basis(cls[k], dd)
        for Bp in split_basis(B)
            A, C = loc_matrix(g, Bp)
            hermitian_psd!(A, C); nblocks += 1
        end
    end
    for (h, tag) in zip(pop.eqs, pop.eq_tags)
        cis_constant(h) && continue
        sup = csupport(h)
        k = findfirst(c -> issubset(sup, c), cls)
        dd = k === nothing ? -1 : 2 * order - cdegree(h)
        if dd < 0
            skipped[tag] = get(skipped, tag, 0) + 1
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
    # second-order cones on the moment image (thermal limits); valid at every order
    for ((a, xs), tag) in zip(pop.socs, pop.soc_tags)
        @constraint(model, vcat(lin_real(a), [lin_real(x) for x in xs]) in SecondOrderCone())
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
    build_time = time() - t0
    optimize!(model)
    st = termination_status(model)
    bound = oscale * (try
        dual_objective_value(model)
    catch
        try objective_value(model) catch; NaN end
    end)
    pobj = oscale * (try objective_value(model) catch; NaN end)
    info = Dict{String,Any}("skipped" => skipped, "n_blocks" => nblocks,
        "t_invariant_possible" => tinv_ok, "obj_scale" => oscale)
    return ComplexMomentRelaxation(st, bound, pobj, build_time, solve_time(model), cls,
        psd_sizes, length(yre), tinv, info)
end
