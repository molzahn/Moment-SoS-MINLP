# Small, enumerable AC-OTS / AC-UC test instances derived from PowerModels test cases (data/).
# All modifications to the original data are made here, in code, so they are documented.

using MomentMINLP, PowerModels

const DATA_DIR = joinpath(@__DIR__, "..", "data")
parse_case(case) = PowerModels.parse_file(joinpath(DATA_DIR, case * ".m"))
ids(d, comp) = sort(parse.(Int, collect(keys(d[comp]))))

"UC modifications: Pmin = frac·Pmax and no-load cost c0 = nl·c1·Pmax (c1 = linear cost coefficient, pu)."
function uc_modify!(data; pmin_frac = 0.4, noload_frac = 0.15, c0 = nothing)
    for (k, g) in data["gen"]
        g["pmin"] = pmin_frac * g["pmax"]
        c = g["cost"]
        while length(c) < 2            # PowerModels may store zero-cost units with a truncated cost vector
            pushfirst!(c, 0.0)
        end
        g["ncost"] = length(c)
        c1 = c[end-1]
        c[end] = c0 === nothing ? noload_frac * c1 * g["pmax"] : c0[parse(Int, k)]
    end
    return data
end

"Branch ids whose removal alone does not disconnect the network."
function nonbridge_branches(d)
    out = Int[]
    for l in ids(d, "branch")
        dd = deepcopy(d)
        dd["branch"][string(l)]["br_status"] = 0
        length(PowerModels.calc_connected_components(dd)) == 1 && push!(out, l)
    end
    return out
end

function load_instance(name::String)
    if name == "case5_ots"
        # PJM 5-bus: every branch switchable (7 binaries, 128 configurations)
        d = parse_case("case5")
        return build_power_pop(d; switchable = ids(d, "branch"), name = name),
            "case5, all 7 branches switchable"
    elseif name == "case9_ots"
        # WSCC 9-bus with branch ratings reduced to 1.0 pu to create congestion (9 binaries)
        d = parse_case("case9")
        for br in values(d["branch"])
            br["rate_a"] = br["rate_b"] = br["rate_c"] = 1.0
        end
        return build_power_pop(d; switchable = ids(d, "branch"), name = name),
            "case9, rate_a = 100 MVA on all branches, all 9 branches switchable"
    elseif name == "case5_uc"
        d = uc_modify!(parse_case("case5"))
        return build_power_pop(d; commitable = ids(d, "gen"), name = name),
            "case5, Pmin = 0.4 Pmax, no-load cost = 0.15 c1 Pmax, all 5 units commitable"
    elseif name == "case5_uc_sym"
        # two identical units at bus 1 (symmetric commitment decisions)
        d = parse_case("case5")
        for k in ("1", "2")
            g = d["gen"][k]
            g["pmax"] = 1.05
            g["qmax"], g["qmin"] = 0.8, -0.8
            g["cost"][end-1] = 1450.0
        end
        uc_modify!(d)
        return build_power_pop(d; commitable = ids(d, "gen"), name = name),
            "case5_uc with units 1 and 2 (bus 1) made identical: Pmax 105 MW, Q ±80 MVAr, c1 = 14.5 \$/MWh"
    elseif name == "case14_uc"
        # reactive-power-driven commitment: units 3–5 are synchronous condensers in the original case
        d = parse_case("case14")
        uc_modify!(d; pmin_frac = 0.2, c0 = Dict(1 => 500.0, 2 => 300.0, 3 => 150.0, 4 => 150.0, 5 => 150.0))
        return build_power_pop(d; commitable = ids(d, "gen"), name = name),
            "case14, Pmin = 0.2 Pmax, no-load costs (500, 300, 150, 150, 150) \$/h, all 5 units commitable"
    elseif name == "case14_ots"
        d = parse_case("case14")
        sw = [1, 2, 3, 4, 5, 6, 7, 10, 13, 16]
        return build_power_pop(d; switchable = sw, name = name),
            "case14, 10 switchable branches $(sw)"
    elseif name == "case30_ots"
        # PGLib IEEE 30-bus: the 15 lowest-numbered non-bridge branches are switchable
        d = parse_case("case30")
        sw = nonbridge_branches(d)[1:15]
        return build_power_pop(d; switchable = sw, name = name),
            "case30 (pglib), 15 switchable non-bridge branches $(sw)"
    elseif name == "case24_ots"
        d = parse_case("case24")
        sw = nonbridge_branches(d)[1:15]
        return build_power_pop(d; switchable = sw, name = name),
            "case24 (pglib RTS), 15 switchable non-bridge branches $(sw)"
    elseif name == "case24_uc"
        # PGLib RTS-96 already has Pmin > 0 and no-load costs; all 33 units commitable (startup costs ignored)
        d = parse_case("case24")
        return build_power_pop(d; commitable = ids(d, "gen"), name = name),
            "case24 (pglib RTS), all 33 units commitable, original Pmin and no-load costs"
    elseif name == "case30_uc"
        d = parse_case("case30")
        uc_modify!(d; pmin_frac = 0.2, c0 = Dict(1 => 300.0, 2 => 300.0, 3 => 100.0, 4 => 100.0, 5 => 100.0, 6 => 100.0))
        return build_power_pop(d; commitable = ids(d, "gen"), name = name),
            "case30 (pglib), Pmin = 0.2 Pmax, no-load costs (300, 300, 100, 100, 100, 100) \$/h, all 6 units commitable"
    end
    error("unknown instance $name")
end

const INSTANCES = ["case5_ots", "case9_ots", "case5_uc", "case5_uc_sym", "case14_uc"]
