#!/bin/bash
# Submit the larger DC-OTS paper cases (one job each), with memory/time scaled to case size.
#   cd ~/scratch/Moment-SoS-MINLP && bash hpc/submit_large.sh
set -euo pipefail
cd "$(dirname "$0")/.."
source hpc/pace_env.sh
pace_check_account
submit() {  # case mem time
    local q; q="$(pace_qos_for "$3")"
    echo "[submit] $1  mem=$2 time=$3 qos=$q account=$PACE_ACCOUNT"
    CASE="$1" sbatch --account="$PACE_ACCOUNT" --qos="$q" --job-name="mm-$1" --mem="$2" --time="$3" --export=ALL hpc/experiment3.sbatch
}
submit 179-goc-api     64G  36:00:00
submit 200-activ-api   64G  36:00:00
submit 240-pserc-api   64G  36:00:00
submit 300-ieee-api   128G  48:00:00
submit 500-goc-api     48G  36:00:00
submit 1354-pegase-api 128G 72:00:00
