# Re-validate every stored certified bound against the current best-known costs.
# Run after anything that LOWERS a best-known cost (e.g. a new incumbent from experiment 5), since that
# can retroactively invalidate a bound that was fine before.
#   julia --project=. scripts/check_all_bounds.jl
include("dcots_cases.jl")
using JSON, Printf

function main()
    bad = 0; checked = 0
    for r in DCOTS_PAPER
        nm = r[1]; bk = best_known_cost(nm)
        for suf in ("", "_unmerged")
            p = joinpath(@__DIR__, "..", "results", "recertify_$(nm)$(suf).json")
            isfile(p) || continue
            for (k, o) in JSON.parsefile(p)["variants"]
                b = get(o, "bound", nothing)
                (b isa Real && isfinite(b)) || continue
                checked += 1
                if !check_bound(b, bk)
                    bad += 1
                    @printf("INVALID %s%s %s: bound %.2f > best known %.2f\n", nm, suf, k, b, bk)
                end
            end
        end
    end
    @printf("checked %d stored bounds against current best-known costs: %d invalid\n", checked, bad)
    bad == 0 || error("$bad stored bound(s) now exceed a known feasible cost")
end
main()
