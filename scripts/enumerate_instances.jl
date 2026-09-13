# Enumerate all binary configurations of each instance with PowerModels AC-OPF (Ipopt) and
# store the table in results/enum_<instance>.json.
#   julia --project=. scripts/enumerate_instances.jl [instance ...]
include("instances.jl")
using JSON, Printf

function save_enumeration(name, pop, desc, ev, elapsed)
    rows = [Dict("bits" => join(Int.(b)), "feasible" => r.feasible, "cost" => r.feasible ? r.cost : nothing,
                 "reason" => r.reason) for (b, r) in ev.cache]
    sort!(rows; by = r -> r["cost"] === nothing ? Inf : r["cost"])
    out = Dict("instance" => name, "description" => desc,
               "binaries" => [string(k, ":", id) for (k, id) in pop.meta["binaries"]],
               "bit_order_note" => "bits[k] is the status of binaries[k]",
               "enumeration_seconds" => elapsed, "configs" => rows)
    open(joinpath(@__DIR__, "..", "results", "enum_$(name).json"), "w") do io
        JSON.print(io, out, 1)
    end
end

function load_enumeration!(ev::ConfigEvaluator, name)
    f = joinpath(@__DIR__, "..", "results", "enum_$(name).json")
    isfile(f) || return false
    js = JSON.parsefile(f)
    for r in js["configs"]
        b = BitVector([c == '1' for c in r["bits"]])
        ev.cache[b] = MomentMINLP.ConfigResult(r["feasible"], r["cost"] === nothing ? Inf : r["cost"], r["reason"], 0.0)
    end
    return true
end

if abspath(PROGRAM_FILE) == @__FILE__
    names = isempty(ARGS) ? INSTANCES : ARGS
    for name in names
        pop, desc = load_instance(name)
        ev = ConfigEvaluator(pop)
        t = @elapsed enumerate_configs!(ev)
        save_enumeration(name, pop, desc, ev, t)
        feas = [(b, r) for (b, r) in ev.cache if r.feasible]
        sort!(feas; by = x -> x[2].cost)
        @printf("\n%s: %s\n  %d binaries, %d configs, %d feasible, %.1fs\n", name, desc,
                length(ev.binaries), length(ev.cache), length(feas), t)
        println("  binaries: ", join([string(k, ":", id) for (k, id) in ev.binaries], " "))
        allon = ev.cache[trues(length(ev.binaries))]
        @printf("  all-on cost: %s\n", allon.feasible ? @sprintf("%.2f", allon.cost) : allon.reason)
        for (b, r) in feas[1:min(8, end)]
            @printf("  %s  %.2f  (+%.3f%%)\n", join(Int.(b)), r.cost, 100 * (r.cost / feas[1][2].cost - 1))
        end
        within1 = count(x -> x[2].cost <= 1.01 * feas[1][2].cost, feas)
        println("  configs within 1% of best: ", within1)
    end
end
