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
import PowerModels

export Poly, pvar, POP, add_var!, add_ineq!, add_eq!, add_pmi!, fix_variables, solve_nlp,
    solve_moment_relaxation, MomentRelaxation, moment, chordal_cliques,
    build_power_pop, ConfigEvaluator, evaluate!, enumerate_configs!, best_config, config_data,
    binvars, marginals, binary_correlation, sample_threshold, sample_independent, sample_gaussian,
    sample_conditional, sample_dive, summarize_samples, mosek_optimizer, ipopt_optimizer

include("polynomials.jl")
include("pop.jl")
include("sparsity.jl")
include("relaxation.jl")
include("power.jl")
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
    PowerModels.silence()
end

mosek_optimizer() = optimizer_with_attributes(MosekTools.Optimizer, "MSK_IPAR_NUM_THREADS" => 4)
ipopt_optimizer() = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0, "sb" => "yes",
    "max_iter" => 3000)

end # module
