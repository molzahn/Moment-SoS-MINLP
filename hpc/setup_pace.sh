#!/bin/bash
# One-time setup on the PACE login node (compute nodes have no internet):
#   cd ~/scratch/Moment-SoS-MINLP && bash hpc/setup_pace.sh
set -euo pipefail
cd "$(dirname "$0")/.."
source hpc/pace_env.sh
mkdir -p "$JULIA_DEPOT_PATH" hpc/logs
# Licensed HSL (MA97 for Ipopt): develop the local HSL_jll (the workstation Manifest points at a Mac path)
export HSL_JLL_PATH="${HSL_JLL_PATH:-$HOME/scratch/parameter_optimized_OTS/HSL/HSL_jll.jl.v2026.8.4}"
"$JULIA" --project=. -e 'using Pkg; Pkg.develop(path = ENV["HSL_JLL_PATH"]); Pkg.instantiate(); Pkg.precompile()'
"$JULIA" --project=. -e 'using MomentMINLP; println("Ipopt linear solver: ", MomentMINLP.IPOPT_LINEAR_SOLVER[])'
"$JULIA" --project=. scripts/check_mosek.jl
