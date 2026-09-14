# Polynomial optimization problem:  min f(x)  s.t.  g_i(x) >= 0,  h_j(x) = 0,  G_k(x) ⪰ 0,
# with some variables binary (x_i ∈ {0,1}, handled via x_i^2 = x_i).

mutable struct POP
    names::Vector{String}
    isbin::Vector{Bool}
    lb::Vector{Float64}
    ub::Vector{Float64}
    start::Vector{Float64}
    obj::Poly
    ineqs::Vector{Poly}
    ineq_tags::Vector{String}
    eqs::Vector{Poly}
    eq_tags::Vector{String}
    pmis::Vector{Matrix{Poly}}
    pmi_scalar::Vector{Poly}   # equivalent scalar inequality (>= 0) used by the NLP solver
    pmi_tags::Vector{String}
    meta::Dict{String,Any}
end

POP() = POP(String[], Bool[], Float64[], Float64[], Float64[], Poly(), Poly[], String[], Poly[], String[],
    Matrix{Poly}[], Poly[], String[], Dict{String,Any}())

nvars(pop::POP) = length(pop.names)
binary_indices(pop::POP) = findall(pop.isbin)

function add_var!(pop::POP, name::String; binary = false, lb = -Inf, ub = Inf, start = 0.0)
    push!(pop.names, name)
    push!(pop.isbin, binary)
    push!(pop.lb, binary ? 0.0 : lb)
    push!(pop.ub, binary ? 1.0 : ub)
    push!(pop.start, start)
    return length(pop.names)
end

add_ineq!(pop::POP, g::Poly, tag::String) = (push!(pop.ineqs, g); push!(pop.ineq_tags, tag); pop)
add_eq!(pop::POP, h::Poly, tag::String) = (push!(pop.eqs, h); push!(pop.eq_tags, tag); pop)
function add_pmi!(pop::POP, G::Matrix{Poly}, scalar::Poly, tag::String)
    push!(pop.pmis, G)
    push!(pop.pmi_scalar, scalar)
    push!(pop.pmi_tags, tag)
    return pop
end

"""
    fix_variables(pop, fixed) -> (newpop, consistent)

Substitute fixed values. Constraints that become constant are dropped; `consistent`
is false if any of them is violated.
"""
function fix_variables(pop::POP, fixed::AbstractDict{Int,<:Real}; tol = 1e-9)
    q = deepcopy(pop)
    ok = true
    for (i, v) in fixed
        q.lb[i] = v
        q.ub[i] = v
        q.start[i] = v
    end
    q.obj = substitute(pop.obj, fixed)
    keep_i = Int[]
    q.ineqs = [substitute(g, fixed) for g in pop.ineqs]
    for (k, g) in enumerate(q.ineqs)
        if is_constant(g)
            constant_term(g) < -tol && (ok = false)
        else
            push!(keep_i, k)
        end
    end
    q.ineqs = q.ineqs[keep_i]
    q.ineq_tags = q.ineq_tags[keep_i]
    keep_e = Int[]
    q.eqs = [substitute(h, fixed) for h in pop.eqs]
    for (k, h) in enumerate(q.eqs)
        if is_constant(h)
            abs(constant_term(h)) > tol && (ok = false)
        else
            push!(keep_e, k)
        end
    end
    q.eqs = q.eqs[keep_e]
    q.eq_tags = q.eq_tags[keep_e]
    q.pmis = [map(p -> substitute(p, fixed), G) for G in pop.pmis]
    q.pmi_scalar = [substitute(p, fixed) for p in pop.pmi_scalar]
    q.meta["fixed"] = merge(get(pop.meta, "fixed", Dict{Int,Float64}()), Dict{Int,Float64}(fixed))
    return q, ok
end

function _jump_poly(model, x, p::Poly)
    terms = Any[]
    for (m, c) in p.terms
        push!(terms, isempty(m) ? c : c * prod(x[v] for v in m))
    end
    isempty(terms) && return 0.0
    return sum(terms)
end

"""
    solve_nlp(pop; fixed, optimizer) -> NamedTuple

Local solution of the POP (binaries relaxed to [0,1] unless fixed) with Ipopt.
Used to validate formulations and to recover continuous variables.
"""
function solve_nlp(pop::POP; fixed = Dict{Int,Float64}(), optimizer = ipopt_optimizer(), start = pop.start)
    n = nvars(pop)
    model = Model(optimizer)
    @variable(model, x[1:n])
    for i in 1:n
        isfinite(pop.lb[i]) && set_lower_bound(x[i], pop.lb[i])
        isfinite(pop.ub[i]) && set_upper_bound(x[i], pop.ub[i])
        set_start_value(x[i], start[i])
    end
    for (i, v) in fixed
        fix(x[i], v; force = true)
    end
    for g in pop.ineqs
        @constraint(model, _jump_poly(model, x, g) >= 0)
    end
    for h in pop.eqs
        @constraint(model, _jump_poly(model, x, h) == 0)
    end
    for s in pop.pmi_scalar
        @constraint(model, _jump_poly(model, x, s) >= 0)
    end
    @objective(model, Min, _jump_poly(model, x, pop.obj))
    optimize!(model)
    st = termination_status(model)
    ok = st in (MOI.LOCALLY_SOLVED, MOI.ALMOST_LOCALLY_SOLVED, MOI.OPTIMAL)
    return (status = st, feasible = ok, objective = ok ? objective_value(model) : Inf,
        x = ok ? value.(x) : fill(NaN, n))
end
