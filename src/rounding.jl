# Randomized rounding of binary variables using pseudo-moments from a moment relaxation.
#
# Schemes (see docs/formulations.md §5):
#   :threshold    deterministic, z_i = [y_i >= 1/2]
#   :independent  z_i ~ Bernoulli(y_i) independently
#   :gaussian     Gaussian copula matching marginals y_i and pairwise correlations from y_ij
#   :conditional  sequential sampling; P(z_j | already-sampled z_B) from pseudo-moments of B∪{j}
#                 via inclusion–exclusion (|B| <= 2t-1, restricted to moments present in a clique)
#   :dive         sequential sampling with the relaxation re-solved after fixing each binary

binvars(pop::POP) = pop.meta["binary_vars"]

function marginals(rel::MomentRelaxation, bvars::Vector{Int})
    return [clamp(coalesce(moment(rel, [v]), 0.5), 0.0, 1.0) for v in bvars]
end

"Correlation matrix of the binaries implied by first/second pseudo-moments (0 when y_ij is absent)."
function binary_correlation(rel::MomentRelaxation, bvars::Vector{Int})
    n = length(bvars)
    μ = marginals(rel, bvars)
    Σ = zeros(n, n)
    for i in 1:n, j in 1:n
        if i == j
            Σ[i, i] = μ[i] * (1 - μ[i])
        else
            yij = moment(rel, sort([bvars[i], bvars[j]]))
            Σ[i, j] = ismissing(yij) ? 0.0 : yij - μ[i] * μ[j]
        end
    end
    d = sqrt.(max.(diag(Σ), 0.0))
    R = Matrix{Float64}(I, n, n)
    for i in 1:n, j in 1:n
        if i != j && d[i] > 1e-6 && d[j] > 1e-6
            R[i, j] = clamp(Σ[i, j] / (d[i] * d[j]), -1.0, 1.0)
        end
    end
    return μ, Σ, R
end

sample_threshold(rel, bvars) = [BitVector(marginals(rel, bvars) .>= 0.5)]

function sample_independent(rel, bvars, N; rng = Random.default_rng())
    μ = marginals(rel, bvars)
    return [BitVector(rand(rng, length(μ)) .< μ) for _ in 1:N]
end

function sample_gaussian(rel, bvars, N; rng = Random.default_rng())
    μ, _, R = binary_correlation(rel, bvars)
    F = eigen(Symmetric(R))
    λ = max.(F.values, 0.0)
    L = F.vectors * Diagonal(sqrt.(λ))
    scale = sqrt.(max.(sum(abs2, L; dims = 2)[:], 1e-12))   # renormalize to unit variances
    L = L ./ scale
    thr = [p <= 0 ? -Inf : p >= 1 ? Inf : quantile(Normal(), p) for p in μ]
    return [BitVector((L * randn(rng, length(μ))) .< thr) for _ in 1:N]
end

"Σ_{T ⊆ zeros} (-1)^|T| y(ones ∪ T ∪ extra): pseudo-probability that z_ones = 1, z_zeros = 0 (and z_extra = 1)."
function _joint(rel, ones_vars::Vector{Int}, zero_vars::Vector{Int}, extra::Vector{Int})
    s = 0.0
    nz = length(zero_vars)
    for mask in 0:(2^nz-1)
        T = [zero_vars[k] for k in 1:nz if (mask >> (k - 1)) & 1 == 1]
        s += (-1)^length(T) * moment(rel, sort([ones_vars; T; extra]))
    end
    return s
end

function sample_conditional(rel, bvars, N; rng = Random.default_rng(), kmax = 2 * maximum(rel.orders) - 1)
    n = length(bvars)
    μ, _, R = binary_correlation(rel, bvars)
    samples = BitVector[]
    neg_mass = 0.0
    for _ in 1:N
        z = falses(n)
        assigned = Int[]
        for j in randperm(rng, n)
            cand = sort(assigned; by = i -> -abs(R[i, j]))
            B = Int[]
            for i in cand
                (length(B) >= kmax || abs(R[i, j]) < 1e-6) && break
                ismissing(moment(rel, sort(bvars[[B; i; j]]))) && continue
                push!(B, i)
            end
            p = μ[j]
            if !isempty(B)
                ones_v = [bvars[i] for i in B if z[i]]
                zero_v = [bvars[i] for i in B if !z[i]]
                p1 = _joint(rel, ones_v, zero_v, [bvars[j]])
                pB = _joint(rel, ones_v, zero_v, Int[])
                p0 = pB - p1
                neg_mass += max(-p1, 0.0) + max(-p0, 0.0)
                p1, p0 = max(p1, 0.0), max(p0, 0.0)
                p = p1 + p0 > 1e-9 ? p1 / (p1 + p0) : μ[j]
            end
            z[j] = rand(rng) < p
            push!(assigned, j)
        end
        push!(samples, z)
    end
    return samples, neg_mass / max(N * n, 1)
end

_rel_ok(r) = r.primal_status in (MOI.FEASIBLE_POINT, MOI.NEARLY_FEASIBLE_POINT) &&
             r.status in (MOI.OPTIMAL, MOI.ALMOST_OPTIMAL, MOI.SLOW_PROGRESS)

"""
Randomized diving: sample a binary from the current relaxation's marginal, fix it, re-solve.
If the relaxation with the sampled value is infeasible, the opposite value is used.
"""
function sample_dive(pop::POP, rel0::MomentRelaxation, N; rng = Random.default_rng(), order = 1, kwargs...)
    bvars = binvars(pop)
    n = length(bvars)
    cache = Dict{Vector{Pair{Int,Float64}},Union{Nothing,MomentRelaxation}}()
    function relax(fixed)
        key = sort(collect(fixed))
        get!(cache, key) do
            q, ok = fix_variables(pop, fixed)
            ok || return nothing
            r = solve_moment_relaxation(q; order = order, kwargs...)
            _rel_ok(r) ? r : nothing
        end
    end
    samples = BitVector[]
    nsolves = 0
    for _ in 1:N
        fixed = Dict{Int,Float64}()
        rel = rel0
        z = falses(n)
        perm = randperm(rng, n)
        for (step, j) in enumerate(perm)
            p = clamp(coalesce(moment(rel, [bvars[j]]), 0.5), 0.0, 1.0)
            v = rand(rng) < p
            if step == n
                z[j] = v
                break
            end
            before = length(cache)
            r = relax(merge(fixed, Dict(bvars[j] => Float64(v))))
            if r === nothing
                v = !v
                r = relax(merge(fixed, Dict(bvars[j] => Float64(v))))
            end
            nsolves += length(cache) - before
            z[j] = v
            fixed[bvars[j]] = Float64(v)
            if r === nothing          # both branches infeasible: finish by independent rounding
                for jj in perm[step+1:end]
                    z[jj] = rand(rng) < clamp(coalesce(moment(rel, [bvars[jj]]), 0.5), 0.0, 1.0)
                end
                break
            end
            rel = r
        end
        push!(samples, z)
    end
    return samples, nsolves
end

"""
Summary statistics of a sample set given the optimal cost. Best-of-k probabilities assume
i.i.d. draws from the empirical sample distribution: P = 1 - (1 - q)^k.
"""
function summarize_samples(ev::ConfigEvaluator, samples::Vector{BitVector}, opt::Float64; ks = (1, 5, 10, 20), tol = 1e-4)
    res = [evaluate!(ev, s) for s in samples]
    costs = [r.cost for r in res]
    q_feas = mean(isfinite.(costs))
    q_opt = mean(costs .<= opt + tol * abs(opt))
    q_1pct = mean(costs .<= opt + 0.01 * abs(opt))
    reasons = Dict{String,Int}()
    for r in res
        r.feasible || (reasons[r.reason] = get(reasons, r.reason, 0) + 1)
    end
    return (n = length(samples), unique = length(unique(samples)), feas = q_feas, opt = q_opt, within1 = q_1pct,
        best = minimum(costs), opt_k = Dict(k => 1 - (1 - q_opt)^k for k in ks),
        feas_k = Dict(k => 1 - (1 - q_feas)^k for k in ks), infeasible_reasons = reasons)
end
