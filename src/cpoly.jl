# Complex polynomials and complex polynomial optimization problems, for the complex
# moment/sum-of-squares hierarchy of Josz & Molzahn, "Moment/Sum-of-Squares Hierarchy for Complex
# Polynomial Optimization" (arXiv:1508.02068).
#
# WHY A PARALLEL TYPE RATHER THAN A REAL ENCODING. Converting complex to real and relaxing do not
# commute: converting first and relaxing second is what the existing real hierarchy does, and it
# gives a TIGHTER relaxation at a given order but a much larger one. Over n complex variables the
# order-d basis has binomial(n+d, d) monomials against binomial(2n+d, d) for the real formulation
# in 2n real variables -- on a 6-bus clique at order 2 that is 28 against 105. A Hermitian block of
# size k costs a real PSD block of size 2k, so 56 against 105: still about half. The hierarchy is
# weaker per order and cheaper per order, which is the trade the selective algorithm exploits, by
# promoting more constraints to order 2 inside the same budget.
#
# REPRESENTATION. A term is a pair of sorted multi-indices (alpha, beta) meaning conj(z)^alpha
# z^beta, following the paper's y_{alpha,beta} = integral of conj(z)^alpha z^beta. A polynomial is
# REAL-VALUED exactly when its coefficients satisfy c[(beta,alpha)] = conj(c[(alpha,beta)]); every
# objective and constraint here is of that kind, and `is_real_valued` checks it rather than
# assuming it.

const CMono = Tuple{Vector{Int},Vector{Int}}     # (alpha, beta) -> conj(z)^alpha * z^beta

struct CPoly
    terms::Dict{CMono,ComplexF64}
end

CPoly() = CPoly(Dict{CMono,ComplexF64}())
CPoly(c::Number) = c == 0 ? CPoly() : CPoly(Dict((Int[], Int[]) => ComplexF64(c)))
"z_i, as a complex polynomial."
cvar(i::Int) = CPoly(Dict((Int[], [i]) => ComplexF64(1)))
"conj(z_i), as a complex polynomial."
cbar(i::Int) = CPoly(Dict(([i], Int[]) => ComplexF64(1)))

Base.zero(::Type{CPoly}) = CPoly()
Base.copy(p::CPoly) = CPoly(copy(p.terms))
Base.isempty(p::CPoly) = isempty(p.terms)

function caddterm!(p::CPoly, m::CMono, c::Number)
    c == 0 && return p
    v = get(p.terms, m, ComplexF64(0)) + c
    abs(v) < 1e-14 ? delete!(p.terms, m) : (p.terms[m] = v)
    return p
end

function Base.:+(a::CPoly, b::CPoly)
    r = copy(a)
    for (m, c) in b.terms; caddterm!(r, m, c); end
    return r
end
Base.:+(a::CPoly, b::Number) = a + CPoly(b)
Base.:+(a::Number, b::CPoly) = b + CPoly(a)
Base.:-(a::CPoly) = CPoly(Dict(m => -c for (m, c) in a.terms))
Base.:-(a::CPoly, b::CPoly) = a + (-b)
Base.:-(a::CPoly, b::Number) = a + CPoly(-b)
Base.:-(a::Number, b::CPoly) = CPoly(a) + (-b)
Base.:*(a::Number, b::CPoly) = a == 0 ? CPoly() : CPoly(Dict(m => a * c for (m, c) in b.terms))
Base.:*(b::CPoly, a::Number) = a * b
function Base.:*(a::CPoly, b::CPoly)
    r = CPoly()
    for ((a1, b1), c1) in a.terms, ((a2, b2), c2) in b.terms
        caddterm!(r, (sort!(vcat(a1, a2)), sort!(vcat(b1, b2))), c1 * c2)
    end
    return r
end

"Complex conjugate: conj(z)^a z^b -> conj(z)^b z^a with the conjugated coefficient."
cconj(p::CPoly) = CPoly(Dict((b, a) => conj(c) for ((a, b), c) in p.terms))

"True when the polynomial takes real values for every z, i.e. c[(b,a)] == conj(c[(a,b)])."
function is_real_valued(p::CPoly; tol = 1e-10)
    for ((a, b), c) in p.terms
        abs(get(p.terms, (b, a), ComplexF64(0)) - conj(c)) <= tol || return false
    end
    return true
end

"Total degree, max over terms of |alpha| + |beta| (the paper's notion for R_l[zbar, z])."
cdegree(p::CPoly) = isempty(p.terms) ? 0 : maximum(length(a) + length(b) for (a, b) in keys(p.terms))
"Variables the polynomial touches."
function csupport(p::CPoly)
    v = Int[]
    for (a, b) in keys(p.terms); append!(v, a); append!(v, b); end
    return sort!(unique(v))
end
cis_constant(p::CPoly) = all(isempty(a) && isempty(b) for (a, b) in keys(p.terms))

"""
Complex polynomial optimization problem: min f(z) s.t. g_i(z) >= 0, h_j(z) = 0, with f, g, h
real-valued complex polynomials. Bounds are on |z_i|, which is what a voltage magnitude limit is.
"""
mutable struct CPOP
    names::Vector{String}
    absmin::Vector{Float64}       # lower bound on |z_i|
    absmax::Vector{Float64}       # upper bound on |z_i|
    obj::CPoly
    ineqs::Vector{CPoly}
    ineq_tags::Vector{String}
    eqs::Vector{CPoly}
    eq_tags::Vector{String}
    # Second-order cones on the moment image: (a, xs) means L(a) >= ||L(xs)||_2, with a and every
    # x a REAL-VALUED CPoly. This is how thermal limits enter. |S_lm|^2 <= rate^2 is degree 4 in v
    # and so unreachable at order 1 as a polynomial inequality, but P_lm and Q_lm are each degree
    # (1,1), and for any feasible v the triple (rate, P(v), Q(v)) lies in the cone. The cone is
    # convex and L is linear, so the moment image lies in it too -- valid at EVERY order, order 1
    # included, and with no lifting.
    socs::Vector{Tuple{CPoly,Vector{CPoly}}}
    soc_tags::Vector{String}
    # Auxiliary REAL variables that are not complex coordinates and never enter a moment matrix:
    # one per generator at a bus carrying more than one. Stored as (name, lo, hi).
    aux::Vector{Tuple{String,Float64,Float64}}
    # Constraints mixing the auxiliary reals with the moments:
    #   sum_i coef_i * aux_i + L(p)  {== 0 | >= 0}
    # Linear in both, so the relaxation stays a valid lower bound for the same reason the ordinary
    # moment constraints do: L is linear and the relation holds pointwise on the feasible set.
    mixed::Vector{Tuple{Dict{Int,Float64},CPoly,Symbol}}
    mixed_tags::Vector{String}
    meta::Dict{String,Any}
end

CPOP() = CPOP(String[], Float64[], Float64[], CPoly(), CPoly[], String[], CPoly[], String[],
              Tuple{CPoly,Vector{CPoly}}[], String[], Tuple{String,Float64,Float64}[],
              Tuple{Dict{Int,Float64},CPoly,Symbol}[], String[], Dict{String,Any}())

"""Add L(a) >= ||L(xs)||_2. Every argument must be real-valued; the cone is on the moments."""
function add_csoc!(pop::CPOP, a::CPoly, xs::Vector{CPoly}, tag::String = "")
    is_real_valued(a) || error("SOC bound $tag is not real-valued")
    for x in xs
        is_real_valued(x) || error("SOC component of $tag is not real-valued")
    end
    push!(pop.socs, (a, xs)); push!(pop.soc_tags, tag); return pop
end

"""Add a real auxiliary variable (not a complex coordinate, never in a moment matrix)."""
function add_caux!(pop::CPOP, name::String; lo = -Inf, hi = Inf)
    push!(pop.aux, (name, Float64(lo), Float64(hi)))
    return length(pop.aux)
end

"""Add `sum_i coef_i * aux_i + L(p) {== 0 | >= 0}`, with `sense` one of `:eq`, `:geq`."""
function add_cmixed!(pop::CPOP, coefs::Dict{Int,Float64}, p::CPoly, sense::Symbol, tag::String = "")
    sense in (:eq, :geq) || error("sense must be :eq or :geq")
    is_real_valued(p) || error("mixed constraint $tag has a non-real-valued polynomial part")
    push!(pop.mixed, (coefs, p, sense)); push!(pop.mixed_tags, tag); return pop
end

cnvars(pop::CPOP) = length(pop.names)

function add_cvar!(pop::CPOP, name::String; absmin = 0.0, absmax = Inf)
    push!(pop.names, name); push!(pop.absmin, Float64(absmin)); push!(pop.absmax, Float64(absmax))
    return length(pop.names)
end

function add_cineq!(pop::CPOP, g::CPoly, tag::String = "")
    is_real_valued(g) || error("inequality $tag is not real-valued: g and conj(g) must agree")
    push!(pop.ineqs, g); push!(pop.ineq_tags, tag); return pop
end

function add_ceq!(pop::CPOP, h::CPoly, tag::String = "")
    is_real_valued(h) || error("equality $tag is not real-valued: h and conj(h) must agree")
    push!(pop.eqs, h); push!(pop.eq_tags, tag); return pop
end

"""
    cmonomial_basis(vars, d) -> Vector{Vector{Int}}

Monomials z^alpha with |alpha| <= d over `vars`, as sorted index vectors. This is the complex
basis: it indexes by alpha ALONE, where the real formulation over the same buses would index by
pairs. That single fact is where the size advantage comes from.
"""
function cmonomial_basis(vars::Vector{Int}, d::Int)
    basis = [Int[]]
    cur = [Int[]]
    for _ in 1:d
        nxt = Vector{Int}[]
        for m in cur, v in vars
            (isempty(m) || v >= m[end]) && push!(nxt, vcat(m, v))
        end
        append!(basis, nxt)
        cur = nxt
    end
    return basis
end
