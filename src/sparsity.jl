# Correlative sparsity: chordal extension of the variable interaction graph via greedy
# minimum-degree elimination, returning the maximal cliques.

function chordal_cliques(supports::Vector{Vector{Int}})
    adj = Dict{Int,Set{Int}}()
    for s in supports
        for v in s
            haskey(adj, v) || (adj[v] = Set{Int}())
        end
        for a in s, b in s
            a != b && push!(adj[a], b)
        end
    end
    remaining = Set(keys(adj))
    cliques = Vector{Vector{Int}}()
    while !isempty(remaining)
        # min-degree vertex (ties broken by index for determinism)
        v = minimum(u -> (length(adj[u]), u), remaining)[2]
        nb = collect(adj[v])
        push!(cliques, sort!([v; nb]))
        for a in nb, b in nb
            a != b && push!(adj[a], b)
        end
        for a in nb
            delete!(adj[a], v)
        end
        delete!(remaining, v)
        delete!(adj, v)
    end
    # keep maximal cliques only
    sort!(cliques; by = length, rev = true)
    maximal = Vector{Vector{Int}}()
    for c in cliques
        any(issubset(c, m) for m in maximal) || push!(maximal, c)
    end
    return maximal
end
