# Islanding: which switching configurations can possibly be feasible, and how to keep rounded samples
# (and the relaxation) away from the ones that cannot.
#
# PowerModels' AC-OTS model has no connectivity requirement: every bus stays energized, each island
# must balance its own power, and islands need no angle reference. The paper's own AC-OTS topologies
# contain islands (89-pegase: three buses without load/generation; 118-ieee: a generator-only island).
# We therefore do not reject islands per se. We reject an island only when it provably cannot balance
# active power: with nonnegative branch losses (r, g_fr, g_to >= 0) and nonnegative shunt conductances,
# an island S needs
#     Σ_{g in S} pmax_g  >=  Σ_{i in S} (pd_i + gs_i · vmin_i²).
# A configuration that violates this for some island is infeasible for the POP as well, so fixing lines
# or adding cuts that exclude such configurations keeps the relaxation bound valid.

"Active-power island data on the in-service part of a network (PowerModels data dict)."
struct IslandData
    buses::Vector{Int}                       # bus ids
    idx::Dict{Int,Int}                       # bus id -> index
    demand::Vector{Float64}                  # pd + gs * vmin^2 at each bus
    cap::Vector{Float64}                     # Σ pmax of in-service generators at each bus
    edges::Vector{Tuple{Int,Int,Int}}        # (from idx, to idx, branch id) for in-service branches
    adj::Vector{Vector{Tuple{Int,Int}}}      # bus idx -> (neighbor idx, edge idx)
    lossless_ok::Bool                        # the balance test above is valid
end

function IslandData(d::Dict{String,Any})
    buses = sort([b["index"] for b in values(d["bus"]) if b["bus_type"] != 4])
    idx = Dict(b => k for (k, b) in enumerate(buses))
    n = length(buses)
    demand, cap = zeros(n), zeros(n)
    for ld in values(d["load"])
        get(ld, "status", 1) != 0 && haskey(idx, ld["load_bus"]) && (demand[idx[ld["load_bus"]]] += ld["pd"])
    end
    ok = true
    for sh in values(d["shunt"])
        (get(sh, "status", 1) != 0 && haskey(idx, sh["shunt_bus"])) || continue
        k = idx[sh["shunt_bus"]]
        demand[k] += sh["gs"] * d["bus"][string(sh["shunt_bus"])]["vmin"]^2
        sh["gs"] < 0 && (ok = false)
    end
    for g in values(d["gen"])
        get(g, "gen_status", 1) != 0 && haskey(idx, g["gen_bus"]) && (cap[idx[g["gen_bus"]]] += g["pmax"])
    end
    edges = Tuple{Int,Int,Int}[]
    adj = [Tuple{Int,Int}[] for _ in 1:n]
    for br in values(d["branch"])
        (get(br, "br_status", 1) != 0 && haskey(idx, br["f_bus"]) && haskey(idx, br["t_bus"])) || continue
        (br["br_r"] < 0 || br["g_fr"] < 0 || br["g_to"] < 0) && (ok = false)
        u, v = idx[br["f_bus"]], idx[br["t_bus"]]
        push!(edges, (u, v, br["index"]))
        e = length(edges)
        push!(adj[u], (v, e))
        u != v && push!(adj[v], (u, e))
    end
    isempty(get(d, "dcline", Dict())) || (ok = false)   # dc lines transfer power between islands
    return IslandData(buses, idx, demand, cap, edges, adj, ok)
end

"True if the island with bus indices `members` provably cannot balance active power."
function unservable(isd::IslandData, members)
    isd.lossless_ok || return false
    dem = sum((isd.demand[k] for k in members); init = 0.0)
    cap = sum((isd.cap[k] for k in members); init = 0.0)
    return cap < dem - 1e-6 * max(1.0, abs(dem))
end

"Component label of every bus index, using only edges with `enabled[e]`."
function island_labels(isd::IslandData, enabled::AbstractVector{Bool})
    n = length(isd.buses)
    parent = collect(1:n)
    findr(x) = (while parent[x] != x; parent[x] = parent[parent[x]]; x = parent[x]; end; x)
    for (e, (u, v, _)) in enumerate(isd.edges)
        enabled[e] || continue
        ru, rv = findr(u), findr(v)
        ru != rv && (parent[ru] = rv)
    end
    return [findr(i) for i in 1:n]
end

function _unservable_labels(isd::IslandData, lab::Vector{Int})
    dem, cap = Dict{Int,Float64}(), Dict{Int,Float64}()
    for (k, r) in enumerate(lab)
        dem[r] = get(dem, r, 0.0) + isd.demand[k]
        cap[r] = get(cap, r, 0.0) + isd.cap[k]
    end
    isd.lossless_ok || return Set{Int}()
    return Set(r for r in keys(dem) if cap[r] < dem[r] - 1e-6 * max(1.0, abs(dem[r])))
end

"""
    island_screen(d) -> (ok, reason, n_islands)

`reason` is "island_deficit" when some island provably cannot balance active power.
"""
function island_screen(d::Dict{String,Any})
    isd = IslandData(d)
    lab = island_labels(isd, trues(length(isd.edges)))
    bad = _unservable_labels(isd, lab)
    return (isempty(bad), isempty(bad) ? "" : "island_deficit", length(unique(lab)))
end

# ---------------------------------------------------------------------------------------------
# Guard for sampling: the network with every switchable line closed, plus the binary -> edge map

"""
    IslandGuard(data, binaries)

Island data of `data` with all switchable branches closed, and for each binary the index of its
edge (0 for generator binaries). Used by `repair_islands!` and by the connectivity-aware conditional
sampler (`sample_conditional(...; guard)`).
"""
struct IslandGuard
    isd::IslandData
    edge_of::Vector{Int}                 # binary k -> edge index (0 if not a branch binary)
    binary_of::Vector{Int}               # edge index -> binary k (0 if not switchable)
end

function IslandGuard(data::Dict{String,Any}, binaries)
    d = deepcopy(data)
    for (kind, id) in binaries
        kind == :branch && (d["branch"][string(id)]["br_status"] = 1)
    end
    isd = IslandData(d)
    eidx = Dict(id => e for (e, (_, _, id)) in enumerate(isd.edges))
    edge_of = [kind == :branch ? get(eidx, id, 0) : 0 for (kind, id) in binaries]
    binary_of = zeros(Int, length(isd.edges))
    for (k, e) in enumerate(edge_of)
        e > 0 && (binary_of[e] = k)
    end
    return IslandGuard(isd, edge_of, binary_of)
end

IslandGuard(pop::POP) = IslandGuard(pop.meta["data"], pop.meta["binaries"])

_enabled(g::IslandGuard, closed::AbstractVector{Bool}) = [g.binary_of[e] == 0 || closed[g.binary_of[e]] for e in eachindex(g.isd.edges)]

"""
    repair_islands!(guard, bits, μ) -> number of lines closed

Step 2 of the islanding plan: while some island provably cannot balance active power, close the opened
line with the largest marginal μ among those joining such an island to another island.
"""
function repair_islands!(g::IslandGuard, bits::AbstractVector{Bool}, μ::AbstractVector{<:Real})
    nclosed = 0
    while true
        lab = island_labels(g.isd, _enabled(g, bits))
        bad = _unservable_labels(g.isd, lab)
        isempty(bad) && return nclosed
        best, bestμ = 0, -Inf
        for (k, e) in enumerate(g.edge_of)
            (e == 0 || bits[k]) && continue
            u, v, _ = g.isd.edges[e]
            lab[u] != lab[v] && (lab[u] in bad || lab[v] in bad) && μ[k] > bestμ && ((best, bestμ) = (k, μ[k]))
        end
        best == 0 && return nclosed            # cannot be repaired by closing switchable lines
        bits[best] = true
        nclosed += 1
    end
end

"Buses reachable from bus index `s` over enabled edges, never using edge `skip`."
function _reach(isd::IslandData, enabled, s::Int, skip::Int)
    seen = falses(length(isd.buses))
    seen[s] = true
    stack = [s]
    while !isempty(stack)
        u = pop!(stack)
        for (v, e) in isd.adj[u]
            (e == skip || !enabled[e] || seen[v]) && continue
            seen[v] = true
            push!(stack, v)
        end
    end
    return seen
end

"""
    can_open(guard, closed, k) -> Bool

Step 3: can binary `k` be opened given the current state `closed` (undecided lines count as closed)
without creating an island that provably cannot balance active power? Generator binaries always pass.
"""
function can_open(g::IslandGuard, closed::AbstractVector{Bool}, k::Int)
    e = g.edge_of[k]
    e == 0 && return true
    en = _enabled(g, closed)
    u, v, _ = g.isd.edges[e]
    side_u = _reach(g.isd, en, u, e)
    side_u[v] && return true                      # not a bridge: no new island
    side_v = _reach(g.isd, en, v, e)
    return !unservable(g.isd, findall(side_u)) && !unservable(g.isd, findall(side_v))
end

# ---------------------------------------------------------------------------------------------
# Step 4: structure in the relaxation

"Bridges (branch ids) of the in-service network (Tarjan; parallel branches are not bridges)."
function _bridges(isd::IslandData)
    n = length(isd.buses)
    disc, low = zeros(Int, n), zeros(Int, n)
    t = 0
    br = Int[]
    for s in 1:n
        disc[s] == 0 || continue
        t += 1; disc[s] = low[s] = t
        stack = [(s, 0, 1)]                      # (vertex, edge used to reach it, next adjacency position)
        while !isempty(stack)
            u, pe, pos = stack[end]
            if pos <= length(isd.adj[u])
                stack[end] = (u, pe, pos + 1)
                v, e = isd.adj[u][pos]
                e == pe && continue
                if disc[v] == 0
                    t += 1; disc[v] = low[v] = t
                    push!(stack, (v, e, 1))
                else
                    low[u] = min(low[u], disc[v])
                end
            else
                pop!(stack)
                if !isempty(stack)
                    p = stack[end][1]
                    low[p] = min(low[p], low[u])
                    low[u] > disc[p] && push!(br, isd.edges[pe][3])
                end
            end
        end
    end
    return br
end

"""
    radial_fixings(data, switchable) -> branch ids

Switchable bridges whose removal would leave an island that provably cannot balance active power.
Every feasible configuration keeps them closed, so they can be fixed (removed from the binaries).
"""
function radial_fixings(data::Dict{String,Any}, switchable)
    d = deepcopy(data)
    for l in switchable
        d["branch"][string(l)]["br_status"] = 1
    end
    isd = IslandData(d)
    isd.lossless_ok || return Int[]
    sw = Set(switchable)
    eidx = Dict(id => e for (e, (_, _, id)) in enumerate(isd.edges))
    en = trues(length(isd.edges))
    fixed = Int[]
    for l in _bridges(isd)
        l in sw || continue
        e = eidx[l]
        u, v, _ = isd.edges[e]
        side_u = _reach(isd, en, u, e)
        side_v = _reach(isd, en, v, e)
        (unservable(isd, findall(side_u)) || unservable(isd, findall(side_v))) && push!(fixed, l)
    end
    return sort(fixed)
end

"""
    connectivity_cuts(data, switchable; maxset = 2) -> Vector{Vector{Int}} (branch ids)

Cut constraints Σ_{l ∈ δ(S)} z_l >= 1 for connected bus sets S with |S| <= `maxset` (1, 2 or 3) that
provably cannot balance active power on their own and whose boundary δ(S) consists of switchable
branches only. Valid for every feasible configuration. Duplicate boundaries are removed.
"""
function connectivity_cuts(data::Dict{String,Any}, switchable; maxset::Int = 2)
    d = deepcopy(data)
    for l in switchable
        d["branch"][string(l)]["br_status"] = 1
    end
    isd = IslandData(d)
    isd.lossless_ok || return Vector{Vector{Int}}()
    sw = Set(switchable)
    sets = Set{Vector{Int}}()
    for u in eachindex(isd.buses)
        push!(sets, [u])
        maxset >= 2 || continue
        for (v, _) in isd.adj[u]
            v == u && continue
            push!(sets, sort([u, v]))
            maxset >= 3 || continue
            for (w, _) in isd.adj[v]
                w in (u, v) || push!(sets, sort([u, v, w]))
            end
        end
    end
    cuts = Set{Vector{Int}}()
    for S in sets
        unservable(isd, S) || continue
        inS = Set(S)
        δ = Int[]
        ok = true
        for (a, b, id) in isd.edges
            ((a in inS) == (b in inS)) && continue
            id in sw || (ok = false; break)
            push!(δ, id)
        end
        ok && !isempty(δ) && push!(cuts, sort(δ))
    end
    return sort(collect(cuts); by = c -> (length(c), c))
end
