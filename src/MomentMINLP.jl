module MomentMINLP

using JuMP
using LinearAlgebra
using Random
using Statistics
using Printf
using Distributions: Normal, quantile
import MathOptInterface as MOI
import Mosek
import MosekTools
import Ipopt
import HSL_jll
import PowerModels

export Poly, pvar, POP, CPoly, CPOP, cvar, cbar, cconj, cmonomial_basis,
    build_complex_power_pop, solve_complex_moment_relaxation, complex_ybus, complex_bus_cliques, cmoment, clique_W, rank_one_point, solve_complex_adaptive, complex_injection_mismatch, complex_rank_one_point, certify_complex, complex_certified_value, add_var!, add_ineq!, add_eq!, add_pmi!, add_soc!, fix_variables, solve_nlp,
    solve_moment_relaxation, MomentRelaxation, moment, chordal_cliques, binary_clique_order, adjacent_clique_order, capped_augmentation,
    SOSCertificate, certified_value, free_absorb, project_dual, bundle_certify, smooth_certify,
    merge_low_impedance, low_impedance_groups, safe_merge_exclusions, IslandGuard, island_screen, repair_islands!, can_open, radial_fixings, connectivity_cuts,
    build_power_pop, ConfigEvaluator, evaluate!, enumerate_configs!, best_config, config_data, screen_config,
    binvars, marginals, binary_correlation, sample_threshold, sample_independent, sample_gaussian,
    sample_conditional, sample_cardinality, repair_samples, sample_dive, summarize_samples, mosek_optimizer, ipopt_optimizer

include("polynomials.jl")
include("cpoly.jl")
include("crelaxation.jl")
include("cpower.jl")
include("ccertify.jl")
include("cadaptive.jl")
include("pop.jl")
include("sparsity.jl")
include("relaxation.jl")
include("certify.jl")
include("power.jl")
include("islands.jl")
include("rounding.jl")

"""
MOSEK looks for ~/mosek/mosek.lic by default; also accept ~/mosek/<version>/mosek.lic.
Mosek.jl creates its global environment when it loads, so the path is passed to that
environment directly (setting ENV at this point is too late).
"""
function _find_mosek_license()
    haskey(ENV, "MOSEKLM_LICENSE_FILE") && return
    root = joinpath(homedir(), "mosek")
    (isfile(joinpath(root, "mosek.lic")) || !isdir(root)) && return
    for d in sort(readdir(root); rev = true)
        f = joinpath(root, d, "mosek.lic")
        if isfile(f)
            ENV["MOSEKLM_LICENSE_FILE"] = f
            Mosek.putlicensepath(Mosek.msk_global_env, f)
            return
        end
    end
end

function __init__()
    _find_mosek_license()
    _setup_hsl()
    PowerModels.silence()
end

"Linear solver used by `ipopt_optimizer` (set in `__init__`; \"mumps\" if HSL is unavailable)."
const IPOPT_LINEAR_SOLVER = Ref("mumps")

"""
Use an HSL linear solver in Ipopt when the licensed HSL_jll is installed (its `override` directory
holds libhsl): `ma97` by default, or the solver named in env `IPOPT_LINEAR_SOLVER` (ma27, ma57, ma77,
ma86, ma97, or mumps). The public registry HSL_jll has no libhsl, in which case Ipopt keeps MUMPS.
"""
function _setup_hsl()
    want = lowercase(get(ENV, "IPOPT_LINEAR_SOLVER", "ma97"))
    if want == "mumps"
        IPOPT_LINEAR_SOLVER[] = "mumps"
    elseif HSL_jll.is_available() && isfile(HSL_jll.libhsl_path)
        IPOPT_LINEAR_SOLVER[] = want
    else
        IPOPT_LINEAR_SOLVER[] = "mumps"
        @warn "HSL library not available (licensed HSL_jll not installed); Ipopt will use MUMPS"
    end
end

"MOSEK optimizer; thread count from env `MOSEK_THREADS` (default 4)."
mosek_optimizer() = optimizer_with_attributes(MosekTools.Optimizer,
    "MSK_IPAR_NUM_THREADS" => parse(Int, get(ENV, "MOSEK_THREADS", "4")))
"Ipopt for AC-OPF evaluations and local NLP solves; HSL `IPOPT_LINEAR_SOLVER[]` (default ma97) when available."
function ipopt_optimizer(; print_level::Int = 0)
    attrs = Any["print_level" => print_level, "sb" => "yes", "max_iter" => 3000, "max_cpu_time" => 60.0]
    if IPOPT_LINEAR_SOLVER[] != "mumps"
        push!(attrs, "linear_solver" => IPOPT_LINEAR_SOLVER[], "hsllib" => HSL_jll.libhsl_path)
    end
    return optimizer_with_attributes(Ipopt.Optimizer, attrs...)
end

end # module
