#!/bin/bash
# One-time setup on the PACE login node (compute nodes have no internet):
#   cd ~/scratch/Moment-SoS-MINLP && bash hpc/setup_pace.sh
set -euo pipefail
cd "$(dirname "$0")/.."
source hpc/pace_env.sh
mkdir -p "$JULIA_DEPOT_PATH" hpc/logs
"$JULIA" --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
"$JULIA" --project=. scripts/check_mosek.jl
