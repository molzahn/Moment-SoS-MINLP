# Minimal sparse polynomial type used to state polynomial optimization problems (POPs).
# A monomial is a sorted Vector{Int} of variable indices with multiplicity, e.g. x1^2*x3 -> [1,1,3].

struct Poly
    terms::Dict{Vector{Int},Float64}
end

Poly() = Poly(Dict{Vector{Int},Float64}())
Poly(c::Real) = c == 0 ? Poly() : Poly(Dict(Int[] => Float64(c)))
pvar(i::Int) = Poly(Dict([i] => 1.0))

Base.zero(::Type{Poly}) = Poly()
Base.zero(::Poly) = Poly()
Base.copy(p::Poly) = Poly(copy(p.terms))

function addterm!(p::Poly, m::Vector{Int}, c::Real)
    c == 0 && return p
    v = get(p.terms, m, 0.0) + c
    if abs(v) < 1e-14
        delete!(p.terms, m)
    else
        p.terms[m] = v
    end
    return p
end

function Base.:+(a::Poly, b::Poly)
    r = copy(a)
    for (m, c) in b.terms
        addterm!(r, m, c)
    end
    return r
end
Base.:+(a::Poly, b::Real) = a + Poly(b)
Base.:+(a::Real, b::Poly) = b + Poly(a)
Base.:-(a::Poly) = Poly(Dict(m => -c for (m, c) in a.terms))
Base.:-(a::Poly, b::Poly) = a + (-b)
Base.:-(a::Poly, b::Real) = a + Poly(-b)
Base.:-(a::Real, b::Poly) = Poly(a) + (-b)
Base.:*(a::Real, b::Poly) = a == 0 ? Poly() : Poly(Dict(m => a * c for (m, c) in b.terms))
Base.:*(b::Poly, a::Real) = a * b
function Base.:*(a::Poly, b::Poly)
    r = Poly()
    for (m1, c1) in a.terms, (m2, c2) in b.terms
        addterm!(r, sort!(vcat(m1, m2)), c1 * c2)
    end
    return r
end
Base.:^(a::Poly, k::Integer) = k == 0 ? Poly(1.0) : prod(fill(a, k))

degree(p::Poly) = isempty(p.terms) ? 0 : maximum(length(m) for m in keys(p.terms))
support(p::Poly) = sort!(unique(reduce(vcat, keys(p.terms); init = Int[])))
constant_term(p::Poly) = get(p.terms, Int[], 0.0)
is_constant(p::Poly) = all(isempty(m) for m in keys(p.terms))

"Apply x^2 = x for binary variables (monomial must be sorted)."
function reduce_mono(m::Vector{Int}, isbin::AbstractVector{Bool})
    r = Int[]
    for v in m
        if isbin[v] && !isempty(r) && r[end] == v
            continue
        end
        push!(r, v)
    end
    return r
end

function reduce_poly(p::Poly, isbin::AbstractVector{Bool})
    r = Poly()
    for (m, c) in p.terms
        addterm!(r, reduce_mono(m, isbin), c)
    end
    return r
end

"Substitute fixed numerical values for some variables."
function substitute(p::Poly, fixed::AbstractDict{Int,<:Real})
    r = Poly()
    for (m, c) in p.terms
        coef = c
        rest = Int[]
        for v in m
            if haskey(fixed, v)
                coef *= fixed[v]
            else
                push!(rest, v)
            end
        end
        addterm!(r, rest, coef)
    end
    return r
end

"Evaluate a polynomial at a point."
function evaluate(p::Poly, x::AbstractVector{<:Real})
    s = 0.0
    for (m, c) in p.terms
        t = c
        for v in m
            t *= x[v]
        end
        s += t
    end
    return s
end
