#!/bin/bash
# Environment for running the moment/SOS experiments on GT PACE (Phoenix). Sourced by setup and jobs.
#   source hpc/pace_env.sh
export PACE_ACCOUNT="${PACE_ACCOUNT:-gts-dmolzahn6}"
export MM_ROOT="${MM_ROOT:-$HOME/scratch/Moment-SoS-MINLP}"
# Julia depot on scratch (home quota is small); reuse the Julia 1.12.7 binary installed for the OTS project
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-$MM_ROOT/.julia}"
export JULIA="${JULIA:-$HOME/scratch/parameter_optimized_OTS/julia-1.12.7/bin/julia}"
# Precompile once on the login node for a portable CPU target, identical in every job
export JULIA_CPU_TARGET="${JULIA_CPU_TARGET:-generic;x86-64-v3,clone_all}"
# MOSEK finds its license at ~/mosek/mosek.lic (default location); Ipopt/BLAS single-threaded
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
module purge 2>/dev/null || true
echo "[pace_env] root=$MM_ROOT julia=$JULIA depot=$JULIA_DEPOT_PATH account=$PACE_ACCOUNT"
