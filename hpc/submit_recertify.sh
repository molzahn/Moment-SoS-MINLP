#!/bin/bash
# Re-certify the relaxations of experiments 3 and 4 (scripts/recertify.jl), one job per case.
#   cd ~/scratch/Moment-SoS-MINLP && bash hpc/submit_recertify.sh [case ...]
set -euo pipefail
cd "$(dirname "$0")/.."
source hpc/pace_env.sh
declare -A MEM=( [73-ieee-rts-api]=32G [89-pegase-api]=48G [118-ieee-api]=32G [179-goc-api]=32G [200-activ-api]=32G
                 [240-pserc-api]=48G [300-ieee-api]=64G [500-goc-api]=96G [1354-pegase-api]=160G )
declare -A TIM=( [240-pserc-api]=24:00:00 [300-ieee-api]=24:00:00 [500-goc-api]=36:00:00 [1354-pegase-api]=72:00:00 )
CASES=("$@")
[ ${#CASES[@]} -gt 0 ] || CASES=(3-lmbd-api 5-pjm-api 14-ieee-api 24-ieee-rts-api 30-as-api 30-ieee-api 39-epri-api 57-ieee-api
    60-c-api 73-ieee-rts-api 89-pegase-api 118-ieee-api 179-goc-api 200-activ-api 300-ieee-api)
for c in "${CASES[@]}"; do
    CASE="$c" sbatch --account="$PACE_ACCOUNT" --job-name="mmrc-$c" --mem="${MEM[$c]:-16G}" --time="${TIM[$c]:-16:00:00}" \
        --export=ALL hpc/recertify.sbatch
done
