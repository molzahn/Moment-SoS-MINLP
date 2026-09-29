# Checks for the islanding tools (src/islands.jl), without SDP solves.
#   julia --project=. scripts/test_islands.jl
include("dcots_cases.jl")
using JSON, Printf, Random, Statistics, Test

expand(s) = reduce(vcat, [occursin('–', p) ? collect(parse(Int, split(p, '–')[1]):parse(Int, split(p, '–')[2])) : [parse(Int, p)]
                          for p in strip.(split(s, ","))])
# Paper AC-OTS topologies (Table I; 89 and 118 parsed from the table, they reproduce Table II costs)
const PAPER_TOPO = Dict(
    "24-ieee-rts-api" => [8, 9, 24, 30, 34, 35],
    "73-ieee-rts-api" => [1, 24, 27, 34, 35, 36, 37, 50, 65, 73, 74, 75, 76, 80, 113, 114, 119],
    "89-pegase-api" => expand("11, 13, 20, 22, 27, 28, 33, 38–41, 43, 45–46, 53, 76–77, 88–89, 92–93, 99, 101, 107, 109, 111, 116, 118–125, 127–129, 134, 138–139, 141–142, 149, 151, 155, 159, 169, 173, 184, 187–188, 190, 192, 196, 200"),
    "118-ieee-api" => expand("1, 12, 14–15, 19, 24, 26, 43–44, 57–58, 65, 69, 72, 79–80, 82, 84–87, 100, 102–103, 106, 117, 127, 148–149, 157, 176, 180, 186"))

function opened_data(d, opened)
    dd = deepcopy(d)
    for l in opened
        dd["branch"][string(l)]["br_status"] = 0
    end
    return dd
end

@testset "bridges agree with brute force" begin
    for name in ("30-as-api", "57-ieee-api", "118-ieee-api")
        d = parse_api(paper_row(name)[2])
        isd = MomentMINLP.IslandData(d)
        br = Set(MomentMINLP._bridges(isd))
        ncomp(dd) = length(PowerModels.calc_connected_components(dd))
        base = ncomp(d)
        brute = Set(l for l in parse.(Int, collect(keys(d["branch"]))) if ncomp(opened_data(d, [l])) > base)
        @test br == brute
    end
end

@testset "paper topologies pass the island screen (and 89/118 fail the strict one)" begin
    for (name, opened) in PAPER_TOPO
        d = parse_api(paper_row(name)[2])
        dd = opened_data(d, opened)
        ok, why, nisl = island_screen(dd)
        @test ok
        strict = length(PowerModels.calc_connected_components(dd)) == 1
        @printf("  %-16s islands %d, island screen %s, strict screen %s\n", name, nisl, ok, strict)
        name in ("89-pegase-api", "118-ieee-api") && @test !strict
    end
    # the 118 topology through the evaluator (islands allowed) reproduces the paper's cost
    pop, d = load_dcots_instance("118-ieee-api")
    ev = ConfigEvaluator(pop)
    bits = BitVector([!(id in PAPER_TOPO["118-ieee-api"]) for (_, id) in pop.meta["binaries"]])
    r = evaluate!(ev, bits)
    @printf("  118-ieee paper topology via evaluator: feasible %s, cost %.2f (paper 180312)\n", r.feasible, r.cost)
    @test r.feasible && abs(r.cost - 180312) < 1.0
    @test !evaluate!(ConfigEvaluator(pop; islands = :connected), bits).feasible
end

@testset "radial fixings and cuts are satisfied by known feasible topologies" begin
    for name in ("24-ieee-rts-api", "73-ieee-rts-api", "89-pegase-api", "118-ieee-api", "179-goc-api", "200-activ-api",
                 "240-pserc-api", "300-ieee-api", "500-goc-api", "1354-pegase-api")
        pop, d = load_dcots_instance(name)
        sw = [id for (_, id) in pop.meta["binaries"]]
        t = @elapsed begin
            fixed = radial_fixings(d, sw)
            cuts2 = connectivity_cuts(d, setdiff(sw, fixed); maxset = 2)
            cuts3 = connectivity_cuts(d, setdiff(sw, fixed); maxset = 3)
        end
        topos = Vector{Vector{Int}}()
        haskey(PAPER_TOPO, name) && push!(topos, PAPER_TOPO[name])
        f = joinpath(@__DIR__, "..", "results", "experiment3_$(name).json")
        if isfile(f)
            js = JSON.parsefile(f)
            for o in values(js["relaxations"])
                get(o, "failed", false) && continue
                for k in ("best_rounded_opened", "polished_opened")
                    haskey(o, k) && push!(topos, Int.(o[k]))
                end
            end
        end
        for opened in topos
            ok, _, _ = island_screen(opened_data(d, opened))
            ok || continue                                      # only feasible-by-screen topologies are relevant
            @test isempty(intersect(opened, fixed))
            os = Set(opened)
            @test all(c -> any(l -> !(l in os), c), cuts3)
        end
        @printf("  %-16s binaries %4d  radial fixed %3d  cuts |S|<=2: %4d  |S|<=3: %4d  (%.1fs; %d topologies checked)\n",
            name, length(sw), length(fixed), length(cuts2), length(cuts3), t, length(topos))
    end
end

@testset "repair and guard produce screen-feasible samples" begin
    for name in ("118-ieee-api", "200-activ-api", "300-ieee-api")
        js = JSON.parsefile(joinpath(@__DIR__, "..", "results", "experiment3_$(name).json"))
        pop, d = load_dcots_instance(name)
        μ = Float64.(js["relaxations"]["mixed_nobigM"]["marginals"])
        g = IslandGuard(pop)
        rng = MersenneTwister(1)
        raw = [BitVector(rand(rng, length(μ)) .< μ) for _ in 1:50]
        screen(b) = first(island_screen(config_data(d, pop.meta["binaries"], b)))
        rep, nclosed = repair_samples(g, raw, μ)
        # guard: sequential independent draws in open-first order
        guarded = BitVector[]
        forced = 0
        for _ in 1:50
            closed = trues(length(μ))
            for j in sortperm(μ)
                if rand(rng) >= μ[j]
                    can_open(g, closed, j) ? (closed[j] = false) : (forced += 1)
                end
            end
            push!(guarded, closed)
        end
        @printf("  %-14s raw pass %2d/50 | repaired pass %2d/50 (mean %.1f lines re-closed) | guarded pass %2d/50 (mean %.1f forced)\n",
            name, count(screen, raw), count(screen, rep), mean(nclosed), count(screen, guarded), forced / 50)
        @test all(screen, rep)
        @test all(screen, guarded)
    end
end
