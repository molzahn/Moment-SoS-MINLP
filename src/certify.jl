# Certified lower bounds from inexact SOS solutions, following the idea of Oustry, D'Ambrosio, Liberti and
# Ruiz, "Certified and accurate SDP bounds for the ACOPF problem" (PSCC 2022), extended to our sparse
# moment/SOS relaxations with polynomial equality multipliers and order-2 cliques.
#
# SOS form (objective normalized by `scale`, variables scaled to the box): maximize t subject to, for every
# monomial α,   C_α(θ) := Σ_k ⟨X_k, C_k^α⟩ + Σ_j λ_j h_j^α + t·[α = ∅] = f_α,   X_k ⪰ 0,
# with θ = (Gram entries X_k, equality multipliers λ, t). For any θ, write r = f − C(θ). For every x feasible
# for the POP (h(x) = 0, W_k(x) ⪰ 0, x in the box):
#     f(x) = t + Σ_k ⟨X_k, W_k(x)⟩ + Σ_α r_α x^α.
#
# * "box" certificate (the original one):  f(x) ≥ t − Σ_α |r_α| max|x^α| − Σ_k ρ_k λ_min(X_k)⁻.
# * "implicit" certificate (step 1): every monomial α that is produced by an entry of some *moment* block
#   (G = 1) gets a representative Gram entry; adding r_α / c to that entry makes the identity exact, i.e.
#   the Gram matrices become implicit affine functions A_k(θ) of θ, as in the paper's unconstrained dual:
#     F(θ) = t + Σ_k ρ_k min(λ_min(A_k(θ)), 0) − Σ_{α without representative} |r_α| max|x^α|,
#   with ρ_k ≥ max tr W_k(x) over the box. F is concave and F(θ) ≤ val(POP) for EVERY θ, so it also gives a
#   bound when MOSEK stops without a feasible point (step 2).
# * "bundle" (step 3): maximize F with a proximal bundle method started at MOSEK's point.
#
# Floating-point caveat: eigenvalues and sums are computed in Float64 (the paper uses rational arithmetic).

using SparseArrays

struct SOSCertificate
    C::SparseMatrixCSC{Float64,Int}        # rows = monomials, cols = SOS variables
    f::Vector{Float64}                     # objective coefficients / scale, per row
    tcol::Int                              # column of t
    blocks::Vector{Matrix{Int}}            # Gram blocks: n×n matrix of column indices (symmetric)
    rho::Vector{Float64}                   # trace bound of each block over the box
    repcol::Vector{Int}                    # per row: representative column (0 = none)
    repcoef::Vector{Float64}               # per row: coefficient of the representative column in that row
    mb::Vector{Float64}                    # per row: max |x^α| over the box
    blockrows::Vector{Vector{Int}}         # rows whose representative lies in each block
    scale::Float64
end

"""
Groups of variables with a quadratic bound Σ_{i∈S} x_i² <= β, detected from inequalities of the form
c − Σ_{i∈S} a_i x_i² >= 0 with a_i > 0 (e.g. voltage magnitude and thermal limits). Returns (S, β) pairs.
"""
function square_groups(ineqs)
    groups = Tuple{Vector{Int},Float64}[]
    for g in ineqs
        c = get(g.terms, Int[], 0.0)
        c > 0 || continue
        S, amin, ok = Int[], Inf, true
        for (m, a) in g.terms
            isempty(m) && continue
            if length(m) == 2 && m[1] == m[2] && a < 0
                push!(S, m[1]); amin = min(amin, -a)
            else
                ok = false; break
            end
        end
        ok && length(S) >= 2 && push!(groups, (sort(S), c / amin))
    end
    return groups
end

"Upper bound on Σ_{v∈vars} x_v² over the box and the square groups."
function _sumsq_bound(vars::Vector{Int}, maxabs, groups)
    vs = Set(vars)
    covered = Set{Int}()
    T = 0.0
    for (S, β) in groups
        inter = [v for v in S if v in vs && !(v in covered)]
        length(inter) >= 2 || continue
        box = sum(maxabs[v]^2 for v in inter)
        β < box && (T += β; union!(covered, inter))
    end
    return T + sum((maxabs[v]^2 for v in vars if !(v in covered)); init = 0.0)
end

"Upper bound on Σ_{m∈B} m(x)² using degree-wise bounds (T = bound on Σ x_i², valid when |x_i| <= 1)."
function _basis_sq_bound(B, maxabs, groups)
    mbf(m) = prod((maxabs[v] for v in m); init = 1.0)
    box = Dict{Int,Float64}()
    for m in B
        box[length(m)] = get(box, length(m), 0.0) + mbf(m)^2
    end
    vars = sort(unique(reduce(vcat, B; init = Int[])))
    all(v -> maxabs[v] <= 1 + 1e-12, vars) || return sum(values(box))
    T = _sumsq_bound(vars, maxabs, groups)
    tot = 0.0
    for (d, bx) in box
        tot += d == 0 ? bx : d == 1 ? min(bx, T) : d == 2 ? min(bx, (T^2 + T) / 2) : bx
    end
    return tot
end

"Build the certificate data from the SOS model (after `optimize!`)."
function SOSCertificate(model, coef::Dict{Vector{Int},AffExpr}, obj::Poly, monos, blocks, Xs, t, isbin, scale, maxabs;
    groups = Tuple{Vector{Int},Float64}[])
    vars = all_variables(model)
    col = Dict(v => k for (k, v) in enumerate(vars))
    rows = collect(monos)
    rowof = Dict(m => i for (i, m) in enumerate(rows))
    I, J, V = Int[], Int[], Float64[]
    for (i, m) in enumerate(rows)
        haskey(coef, m) || continue
        for (v, c) in coef[m].terms
            push!(I, i); push!(J, col[v]); push!(V, c)
        end
    end
    C = sparse(I, J, V, length(rows), length(vars))
    f = [get(obj.terms, m, 0.0) / scale for m in rows]
    mbf(m) = prod((maxabs[v] for v in m); init = 1.0)
    mb = mbf.(rows)
    bmats = Matrix{Int}[]
    rho = Float64[]
    repcol, repcoef, repscore = zeros(Int, length(rows)), zeros(length(rows)), fill(Inf, length(rows))
    for ((G, B), X) in zip(blocks, Xs)
        n = size(X, 1)
        push!(bmats, [col[X[i, j]] for i in 1:n, j in 1:n])
        bsq = _basis_sq_bound(B, maxabs, groups)
        tr = sum(sum((abs(c) * mbf(mono) for (mono, c) in G[a, a].terms); init = 0.0) for a in axes(G, 1)) * bsq
        push!(rho, tr)
        # representatives: entries of moment blocks (G = [1]) map to exactly one monomial
        (size(G, 1) == 1 && length(G[1, 1].terms) == 1 && get(G[1, 1].terms, Int[], 0.0) == 1.0) || continue
        nb = length(B)
        for i in 1:nb, j in i:nb
            m = monoprod(B[i], B[j], isbin)
            r = get(rowof, m, 0)
            r == 0 && continue
            c = i == j ? 1.0 : 2.0
            score = tr / c                     # eigenvalue shift per unit residual is 1/c; penalty weight ρ
            if score < repscore[r]
                repscore[r] = score
                repcol[r] = col[X[i, j]]
                repcoef[r] = c
            end
        end
    end
    colblock = zeros(Int, length(vars))
    for (k, B) in enumerate(bmats), c in B
        colblock[c] = k
    end
    blockrows = [Int[] for _ in bmats]
    for (i, c) in enumerate(repcol)
        c > 0 && push!(blockrows[colblock[c]], i)
    end
    return SOSCertificate(C, f, col[t], bmats, rho, repcol, repcoef, mb, blockrows, scale)
end

"""
    certified_value(cert, θ; mode = :implicit, mu = 0.0, gradient = true) -> (value, supergradient, exact)

Lower bound on the POP value (original objective units) valid for ANY θ:
* `mode = :box`: residuals bounded over the box, PSD violation of the explicit Gram matrices penalized;
* `mode = :implicit`: every represented residual absorbed into its representative Gram entry (concave in θ);
* `mode = :hybrid`: per block, the better of the two (value only).
With `mu > 0` (implicit mode) the smooth concave surrogate ρ_k·softmin(softmin_i λ_i(A_k), 0) with temperature
`mu` is returned instead; it is <= the exact implicit value, which is returned as the third output.
Gradients are with respect to θ in normalized objective units.
"""
function certified_value(cert::SOSCertificate, θ::Vector{Float64}; mode::Symbol = :implicit, mu::Float64 = 0.0,
    gradient::Bool = true, implicit = nothing)
    implicit === false && (mode = :box)
    r = cert.f .- cert.C * θ
    if mode == :box
        F = θ[cert.tcol] - sum(abs.(r) .* cert.mb)
        for (k, B) in enumerate(cert.blocks)
            F += cert.rho[k] * min(_mineig(θ, B)[1], 0.0)
        end
        return F * cert.scale, nothing, F * cert.scale
    end
    θa = copy(θ)
    unrep = 0.0
    w = zeros(length(r))
    for i in eachindex(r)
        if cert.repcol[i] > 0
            θa[cert.repcol[i]] += r[i] / cert.repcoef[i]
        else
            unrep += abs(r[i]) * cert.mb[i]
            w[i] = -sign(r[i]) * cert.mb[i]
        end
    end
    if mode == :hybrid
        F = θ[cert.tcol] - unrep
        for (k, B) in enumerate(cert.blocks)
            absorbed = cert.rho[k] * min(_mineig(θa, B)[1], 0.0)
            boxed = cert.rho[k] * min(_mineig(θ, B)[1], 0.0) - sum((abs(r[i]) * cert.mb[i] for i in cert.blockrows[k]); init = 0.0)
            F += max(absorbed, boxed)
        end
        return F * cert.scale, nothing, F * cert.scale
    end
    F = θ[cert.tcol] - unrep
    Fexact = F
    Ga = gradient ? zeros(length(θ)) : Float64[]
    for (k, B) in enumerate(cert.blocks)
        n = size(B, 1)
        if mu <= 0 || n == 1
            λ, u = _mineig(θa, B)
            Fexact += cert.rho[k] * min(λ, 0.0)
            if mu <= 0
                λ >= 0 && continue
                F += cert.rho[k] * λ
                gradient && _addouter!(Ga, B, u, cert.rho[k])
            else                                   # 1×1 block, smooth softmin(λ, 0)
                h, dh = _softmin0(λ, mu)
                F += cert.rho[k] * h
                gradient && (Ga[B[1, 1]] += cert.rho[k] * dh)
            end
        else
            M = Symmetric([θa[B[i, j]] for i in 1:n, j in 1:n])
            E = eigen(M)
            λs = E.values
            Fexact += cert.rho[k] * min(λs[1], 0.0)
            z = exp.(-(λs .- λs[1]) ./ mu)
            fs = λs[1] - mu * log(sum(z))           # softmin of the eigenvalues (<= λ_min)
            h, dh = _softmin0(fs, mu)
            F += cert.rho[k] * h
            if gradient
                wts = z ./ sum(z)
                for q in eachindex(λs)
                    wts[q] < 1e-12 && continue
                    _addouter!(Ga, B, E.vectors[:, q], cert.rho[k] * dh * wts[q])
                end
            end
        end
    end
    gradient || return F * cert.scale, nothing, Fexact * cert.scale
    for i in eachindex(r)
        cert.repcol[i] > 0 && (w[i] = Ga[cert.repcol[i]] / cert.repcoef[i])
    end
    g = Ga .- cert.C' * w
    g[cert.tcol] += 1.0
    return F * cert.scale, g, Fexact * cert.scale
end

"softmin(a, 0) = −μ log(1 + exp(−a/μ)) (<= min(a, 0)) and its derivative."
function _softmin0(a::Float64, mu::Float64)
    if a >= 0
        e = exp(-a / mu)
        return -mu * log1p(e), e / (1 + e)
    else
        e = exp(a / mu)
        return a - mu * log1p(e), 1 / (1 + e)
    end
end

function _addouter!(G::Vector{Float64}, B::Matrix{Int}, u::AbstractVector, c::Float64)
    n = size(B, 1)
    for i in 1:n, j in i:n
        G[B[i, j]] += c * (i == j ? u[i]^2 : 2 * u[i] * u[j])
    end
end

function _mineig(θ::Vector{Float64}, B::Matrix{Int})
    n = size(B, 1)
    n == 1 && return θ[B[1, 1]], [1.0]
    M = Symmetric([θ[B[i, j]] for i in 1:n, j in 1:n])
    E = eigen(M, 1:1)
    return E.values[1], E.vectors[:, 1]
end

"Minimize ½ αᵀQα − aᵀα over the unit simplex (small dense problem; accelerated projected gradient)."
function _simplex_qp(Q::Matrix{Float64}, a::Vector{Float64}; iters::Int = 3000)
    n = length(a)
    L = max(opnorm(Q), 1e-12)
    α = fill(1.0 / n, n)
    z, tk = copy(α), 1.0
    proj(v) = (u = sort(v; rev = true); css = cumsum(u); ρ = findlast(k -> u[k] - (css[k] - 1) / k > 0, 1:n);
               τ = (css[ρ] - 1) / ρ; max.(v .- τ, 0.0))
    for _ in 1:iters
        αn = proj(z .- (Q * z .- a) ./ L)
        tn = (1 + sqrt(1 + 4tk^2)) / 2
        z = αn .+ ((tk - 1) / tn) .* (αn .- α)
        α, tk = αn, tn
    end
    return α
end

"""
    bundle_certify(cert, θ0; maxiter = 200, maxtime = 600.0, tol = 1e-4, m = 0.01, maxcuts = 20)
        -> (best bound, θ, history)

Proximal bundle method (maximization of the concave F) with a bounded, aggregated cutting-plane model,
started at θ0 (MOSEK's point). Every evaluated point yields a valid bound; the best one is returned. Stops when
the predicted gain falls below `tol` times the initial gap between MOSEK's raw objective and its certified value.
"""
function bundle_certify(cert::SOSCertificate, θ0::Vector{Float64}; maxiter::Int = 200, maxtime::Float64 = 600.0,
    tol::Float64 = 1e-4, m::Float64 = 0.01, maxcuts::Int = 20, verbose::Bool = false)
    t0 = time()
    s = cert.scale
    θc = copy(θ0)
    Fc, gc, _ = certified_value(cert, θc)
    Fc /= s
    best, bestθ = Fc, copy(θc)
    hist = [Fc * s]
    gs = [gc]                                  # cut gradients
    as = [Fc]                                  # cut values at the current center: a_i = F_i + g_i·(θc − θ_i)
    gn = norm(gc)
    gap = max(θ0[cert.tcol] - Fc, 1e-9 * max(1.0, abs(Fc)))    # raw t minus the certified value: the room to gain
    κ = gn > 0 ? gn^2 / gap : 1.0                                # first step predicts a gain of about `gap`
    κmin, κmax = κ * 1e-6, κ * 1e6
    nserious, nnull, lastserious = 0, 0, 0
    for it in 1:maxiter
        time() - t0 > maxtime && break
        k = length(gs)
        Gm = reduce(hcat, gs)
        Q = (Gm' * Gm) ./ κ
        α = _simplex_qp(Q, as)
        d = (Gm * α) ./ κ
        pred = minimum(as[i] + dot(gs[i], d) for i in 1:k)
        δ = pred - Fc
        if δ <= tol * gap                       # model sees little to gain: first try longer steps
            κ > 10κmin ? (κ /= 10; continue) : break
        end
        θn = θc .+ d
        Fn, gnw, _ = certified_value(cert, θn)
        Fn /= s
        push!(hist, Fn * s)
        Fn > best && ((best, bestθ) = (Fn, copy(θn)))
        # aggregate the active cuts, then add the new one (linearization at θn, expressed at the center)
        agg_g = Gm * α
        agg_a = dot(α, as)
        if Fn >= Fc + m * δ                    # serious step: move the center
            shiftd = d
            as = [as[i] + dot(gs[i], shiftd) for i in 1:k]
            agg_a += dot(agg_g, shiftd)
            (Fn - Fc >= 0.5δ) && (κ = max(0.5κ, κmin))      # model was reliable: allow longer steps
            θc, Fc = θn, Fn
            nserious += 1
            lastserious = it
            anew = Fn
        else
            nnull += 1
            Fn < Fc && (κ = min(2κ, κmax))                   # the step made things worse: shorten
            anew = Fn + dot(gnw, θc .- θn)
        end
        keep = [i for i in 1:k if α[i] > 1e-9]
        gs, as = gs[keep], as[keep]
        if length(gs) >= maxcuts
            gs, as = [agg_g], [agg_a]
        end
        push!(gs, gnw); push!(as, anew)
        verbose && @printf("  bundle it %3d  F %.6f  center %.6f  best %.6f  δ %.2e  κ %.2e  cuts %d\n", it, Fn * s, Fc * s, best * s, δ * s, κ, length(gs))
    end
    return best * s, bestθ, hist
end

"""
    smooth_certify(cert, θ0; mus = nothing, bias = 0.1, stages = 3, iters = 60, maxtime = 600.0, memory = 10)
        -> (best bound, θ, history)

Step 3 alternative to the bundle method: maximize the smooth concave surrogate F_μ (every value of which is
itself a valid bound, since F_μ <= F) with L-BFGS and backtracking, for a decreasing sequence of temperatures
`mus` (normalized objective units; by default chosen so that the smoothing bias is `bias` times the gap between
MOSEK's raw objective and the hybrid certificate, then divided by 10 per stage). The best exact hybrid
certificate over all iterates is returned.
"""
function smooth_certify(cert::SOSCertificate, θ0::Vector{Float64}; mus = nothing, bias::Float64 = 0.1, stages::Int = 3,
    iters::Int = 60, maxtime::Float64 = 600.0, memory::Int = 10, verbose::Bool = false)
    t0 = time()
    s = cert.scale
    θ = copy(θ0)
    best, bestθ = first(certified_value(cert, θ; mode = :hybrid, gradient = false)), copy(θ)
    hist = [best]
    if mus === nothing
        # temperature so that the smoothing bias Σ_k ρ_k μ log(n_k) is `bias` times the remaining gap
        gapn = max(θ0[cert.tcol] - best / s, 1e-12)
        wsum = sum(cert.rho[k] * log(size(B, 1) + 1.0) for (k, B) in enumerate(cert.blocks))
        mus = [bias * gapn / wsum / 10.0^(i - 1) for i in 1:stages]
    end
    for mu in mus
        F, g, Fx = certified_value(cert, θ; mu = mu)
        Sv, Yv, ρv = Vector{Float64}[], Vector{Float64}[], Float64[]
        step0 = 1.0 / max(norm(g), 1e-12)
        for it in 1:iters
            time() - t0 > maxtime && break
            # two-loop recursion (ascent direction for maximization)
            q = copy(g)
            al = zeros(length(Sv))
            for i in length(Sv):-1:1
                al[i] = ρv[i] * dot(Sv[i], q)
                q .-= al[i] .* Yv[i]
            end
            γ = isempty(Sv) ? step0 : dot(Sv[end], Yv[end]) / dot(Yv[end], Yv[end])
            q .*= abs(γ)
            for i in eachindex(Sv)
                b = ρv[i] * dot(Yv[i], q)
                q .+= (al[i] - b) .* Sv[i]
            end
            d = q
            slope = dot(g, d)
            slope <= 0 && (d = g .* step0; slope = dot(g, d); empty!(Sv); empty!(Yv); empty!(ρv))
            τ, accepted = 1.0, false
            local Fn, gn, Fxn, θn
            for _ in 1:30
                θn = θ .+ τ .* d
                Fn, gn, Fxn = certified_value(cert, θn; mu = mu)
                if Fn >= F + 1e-4 * τ * slope * s
                    accepted = true
                    break
                end
                τ /= 2
            end
            accepted || break
            sk, yk = θn .- θ, g .- gn                # yk = −(∇F(θn) − ∇F(θ)) for the concave objective
            if dot(sk, yk) > 1e-16
                push!(Sv, sk); push!(Yv, yk); push!(ρv, 1 / dot(yk, sk))
                length(Sv) > memory && (popfirst!(Sv); popfirst!(Yv); popfirst!(ρv))
            end
            gain = Fn - F
            θ, F, g = θn, Fn, gn
            Fh = first(certified_value(cert, θ; mode = :hybrid, gradient = false))
            push!(hist, Fh)
            Fh > best && ((best, bestθ) = (Fh, copy(θ)))
            verbose && @printf("  smooth μ %.0e it %3d  Fμ %.4f  exact implicit %.4f  hybrid %.4f  best %.4f\n", mu, it, F, Fxn, Fh, best)
            abs(gain) < 1e-9 * max(1.0, abs(F)) && break
        end
    end
    return best, bestθ, hist
end
