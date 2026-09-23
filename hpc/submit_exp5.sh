#!/bin/bash
# Submit the experiment-5 cardinality sweep, one job per case.
#   cd <repo> && bash hpc/submit_exp5.sh [case ...]
# Every case fits well under the 8 h embers cap, so this whole sweep costs no credits.
set -euo pipefail
cd "$(dirname "$0")/.."
source hpc/pace_env.sh
declare -A MEM=( [500-goc-api]=48G [1354-pegase-api]=64G )
declare -A TIM=( [500-goc-api]=06:00:00 [1354-pegase-api]=08:00:00 )
declare -A NS=(  [500-goc-api]=60      [1354-pegase-api]=30 )
CASES=("$@")
[ ${#CASES[@]} -gt 0 ] || CASES=(89-pegase-api 118-ieee-api 179-goc-api 200-activ-api 240-pserc-api
                                 300-ieee-api 500-goc-api 1354-pegase-api)
pace_check_account
for c in "${CASES[@]}"; do
    m="${MEM[$c]:-32G}"; t="${TIM[$c]:-04:00:00}"; n="${NS[$c]:-100}"; q="$(pace_qos_for "$t")"
    echo "[submit] $c  mem=$m time=$t N=$n qos=$q account=$PACE_ACCOUNT"
    CASE="$c" N="$n" sbatch --account="$PACE_ACCOUNT" --qos="$q" --job-name="mm5-$c" \
        --mem="$m" --time="$t" --export=ALL hpc/experiment5.sbatch
done
