#!/bin/bash
# Submit experiment 4 (islanding-aware rounding) for the DC-OTS paper cases, one job each.
#   cd ~/scratch/Moment-SoS-MINLP && bash hpc/submit_exp4.sh
set -euo pipefail
cd "$(dirname "$0")/.."
source hpc/pace_env.sh
pace_check_account
submit() {  # case mem time
    local q; q="$(pace_qos_for "$3")"
    echo "[submit] $1  mem=$2 time=$3 qos=$q account=$PACE_ACCOUNT"
    CASE="$1" sbatch --account="$PACE_ACCOUNT" --qos="$q" --job-name="mm4-$1" --mem="$2" --time="$3" --export=ALL hpc/experiment4.sbatch
}
submit 89-pegase-api    32G 12:00:00
submit 118-ieee-api     32G 12:00:00
submit 179-goc-api      32G 12:00:00
submit 200-activ-api    32G 12:00:00
submit 240-pserc-api    32G 16:00:00
submit 300-ieee-api     48G 16:00:00
submit 500-goc-api      64G 24:00:00
submit 1354-pegase-api 128G 48:00:00
