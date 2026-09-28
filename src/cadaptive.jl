# Selective complex hierarchy: the multi-ordered relaxation (3.31) driven by the same
# power-injection-mismatch heuristic the real selective algorithm uses (Molzahn & Hiskens,
# Algorithm 1).
#
# The premise this exists to test: the complex hierarchy is WEAKER than the real one at a given
# order but much SMALLER, because the order-d basis over n complex variables has binomial(n+d, d)
# monomials against binomial(2n+d, d) over 2n real ones. If the saving is large enough to promote
# more buses to order 2 inside the same budget, the weaker-per-order hierarchy can still end up
# with the better bound.
#
# THE HEURISTIC. Solve, take each clique's Hermitian W block, factor out its closest rank-one point
# v (leading eigenvector scaled by sqrt(lambda_1)), and compare the power injections that v implies
# with the injections the relaxation itself reports. Buses where those disagree most are where the
# moment matrix is furthest from rank one, so they are the ones promoted. This is the same rule as
# the real version, as agreed, so that the comparison isolates the hierarchy rather than the
# selection.

"""
    complex_rank_one_point(y, cliques, nvars) -> (vhat, ratios)

A GLOBAL rank-one voltage estimate stitched from the per-clique Hermitian W blocks, plus each
clique's lambda2/lambda1.

Stitching is what makes every bus usable. Taking each clique's rank-one point in isolation only
lets you evaluate an injection when that bus's whole electrical neighbourhood happens to sit inside
one clique -- true for 5 of case14's 14 buses and 100 of case200's 200 -- so half the network could
never be scored, and therefore never selected for promotion. The real hierarchy has always built a
global point (`xhat` in moment_adaptive_worker.jl); this is the complex counterpart.

The complex twist is PHASE. T-invariance means v and exp(i*theta)*v are equally valid rank-one
factors, so each clique's point arrives with an arbitrary global phase. Cliques are visited in
order and each is rotated by the phase that best matches the buses already fixed, obtained in
closed form as angle(sum over the overlap of conj(v_local) * vhat). The first clique is anchored so
that its largest-magnitude entry is real and positive.
"""
function complex_rank_one_point(y::Dict{CMono,ComplexF64}, cliques::Vector{Vector{Int}}, n::Int)
    vhat = fill(ComplexF64(NaN), n)
    ratios = zeros(length(cliques))
    order = sortperm(cliques; by = c -> -length(c))       # biggest first: a better phase anchor
    for ci in order
        c = cliques[ci]
        W = clique_W(y, c)
        v, r = rank_one_point(W)
        ratios[ci] = r
        isempty(v) && continue
        ov = [(a, c[a]) for a in eachindex(c) if isfinite(real(vhat[c[a]]))]
        if isempty(ov)
            k = argmax(abs.(v))
            abs(v[k]) > 1e-12 && (v = v .* conj(v[k] / abs(v[k])))
        else
            acc = sum(conj(v[a]) * vhat[g] for (a, g) in ov)
            abs(acc) > 1e-12 && (v = v .* (acc / abs(acc)))
        end
        for a in eachindex(c)
            isfinite(real(vhat[c[a]])) || (vhat[c[a]] = v[a])
        end
    end
    for i in 1:n
        isfinite(real(vhat[i])) || (vhat[i] = ComplexF64(1))
    end
    return vhat, ratios
end

"""
    complex_injection_mismatch(cpop, data, y, cliques) -> Dict{Int,Float64}

Per-bus |S| mismatch in MVA between the injection implied by the global rank-one point and the
injection the relaxation reports. EVERY bus gets a score, so every bus is a promotion candidate.
"""
function complex_injection_mismatch(pop::CPOP, data::Dict{String,Any},
        y::Dict{CMono,ComplexF64}, cliques::Vector{Vector{Int}})
    Y = pop.meta["ybus"]::Matrix{ComplexF64}
    buses = pop.meta["buses"]::Vector{Int}
    baseMVA = Float64(get(data, "baseMVA", 100.0))
    n = length(buses)
    vhat, _ = complex_rank_one_point(y, cliques, n)
    mism = Dict{Int,Float64}()
    for c in 1:n
        nb = [j for j in axes(Y, 2) if Y[c, j] != 0]
        isempty(nb) && continue
        s_rank  = vhat[c] * conj(sum(Y[c, j] * vhat[j] for j in nb))
        s_relax = sum(conj(Y[c, j]) * cmoment(y, [j], [c]) for j in nb)
        mism[buses[c]] = abs(s_rank - s_relax) * baseMVA
    end
    return mism
end

"""
    solve_complex_adaptive(cpop, data; h, max_iter, ...) -> NamedTuple

Iteratively promote the `h` worst-mismatch buses to order 2 and re-solve, keeping the best bound
seen. The bound is the running MAXIMUM, as in the real version: each iteration's relaxation is
strictly tighter than the last, so any decrease is solver error rather than a real regression.
"""
function solve_complex_adaptive(pop::CPOP, data::Dict{String,Any};
        h::Int = 3, max_iter::Int = 10, mismatch_tol::Float64 = 1e-3,
        max_order::Int = 2, total_seconds::Float64 = 3600.0,
        max_seconds::Float64 = 900.0, verbose::Bool = true, certify::Bool = true,
        psd_mode::Symbol = :full, minor_kmax::Int = 2, psd_threshold::Int = 0,
        minor_core_degree::Int = 1,
        adjacency::Union{Nothing,Set{Tuple{Int,Int}}} = nothing)
    t0 = time()
    idx = pop.meta["bus_index"]::Dict{Int,Int}
    buses = pop.meta["buses"]::Vector{Int}
    c0 = Float64(get(pop.meta, "const_cost", 0.0))
    promoted = Int[]
    best = -Inf
    iters = Any[]
    for it in 1:max_iter
        cls = complex_bus_cliques(data, idx; order2_buses = promoted)
        # a clique is order 2 exactly when it covers the closed neighbourhood of a promoted bus
        ords = ones(Int, length(cls))
        for b in promoted
            nb = [idx[j] for j in _closed_nbrs(data, b) if haskey(idx, j)]
            for (ci, c) in enumerate(cls)
                issubset(nb, c) && (ords[ci] = max_order)
            end
        end
        # The determinant relaxation is WEAKER than full PSD at a fixed order, so it cannot tighten
        # anything directly. The premise is tractability: where a full order-2 block of 156-240 rows
        # sends MOSEK to TIME_LIMIT and the bound degenerates, small cones may solve, and a solved
        # weaker relaxation at order 2 can still beat a solved exact one at order 1.
        rel = solve_complex_moment_relaxation(pop; cliques = cls, clique_orders = ords,
            certify = certify, psd_mode = psd_mode, minor_kmax = minor_kmax,
            psd_threshold = psd_threshold, adjacency = adjacency,
            minor_core_degree = minor_core_degree,
            solver_params = Dict{String,Any}("MSK_DPAR_OPTIMIZER_MAX_TIME" => max_seconds))
        raw = rel.bound + c0
        # The RIGOROUS bound is what counts, and it is what the running maximum is taken over. A
        # promotion makes the relaxation strictly tighter, so a fall is solver error, not a real
        # regression -- the same reasoning as the real selective algorithm.
        cb, cinfo = certify ? certify_complex(pop, rel) : (NaN, Dict{String,Any}())
        bound = certify && isfinite(cb) ? cb + c0 : raw
        isfinite(bound) && bound > best && (best = bound)
        mism = isempty(rel.y) ? Dict{Int,Float64}() :
               complex_injection_mismatch(pop, data, rel.y, cls)
        # only UNPROMOTED buses count: a promoted bus keeps its mismatch, and including it made
        # the reported maximum constant across iterations and the stopping rule unreachable
        rem = [v for (b, v) in mism if !(b in promoted)]
        worst = isempty(rem) ? 0.0 : maximum(rem)
        push!(iters, (iter = it, bound = bound, raw = raw, best = best,
                      psd_corr = get(cinfo, "psd_correction", NaN),
                      box_corr = get(cinfo, "box_correction", NaN),
                      max_resid = get(cinfo, "max_resid", NaN), status = string(rel.status),
                      n_order2 = count(==(max_order), ords), psd_max = maximum(rel.psd_sizes; init = 0),
                      n_cliques = length(cls), max_mismatch_MVA = worst,
                      solve_time = rel.solve_time))
        verbose && println("  iter ", lpad(it, 2),
            "  raw ", round(raw, digits = 4),
            "  certified ", round(bound, digits = 4),
            "  best ", round(best, digits = 4),
            "  |r| ", round(get(cinfo, "max_resid", NaN), sigdigits = 2),
            "  max_mismatch ", round(worst, sigdigits = 4), " MVA",
            "  order2 ", count(==(max_order), ords), "/", length(cls),
            "  psd_max ", maximum(rel.psd_sizes; init = 0),
            "  ", rel.status, "  ", round(rel.solve_time, digits = 1), "s")
        flush(stdout)
        worst < mismatch_tol && break
        time() - t0 > total_seconds && break
        cand = sort([b for b in buses if !(b in promoted)]; by = b -> -get(mism, b, 0.0))
        isempty(cand) && break
        append!(promoted, cand[1:min(h, length(cand))])
    end
    return (bound = best, iterations = iters, promoted = promoted, time = time() - t0)
end

"Closed neighbourhood of a bus over in-service branches."
function _closed_nbrs(data::Dict{String,Any}, b::Int)
    s = Set{Int}([b])
    for (_, br) in get(data, "branch", Dict{String,Any}())
        Int(get(br, "br_status", 1)) == 0 && continue
        f, t = Int(br["f_bus"]), Int(br["t_bus"])
        f == b && push!(s, t)
        t == b && push!(s, f)
    end
    return collect(s)
end
