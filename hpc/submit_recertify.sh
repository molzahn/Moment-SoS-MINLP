#!/bin/bash
# Re-certify the relaxations of experiments 3 and 4 (scripts/recertify.jl), one job per case.
#   cd ~/scratch/Moment-SoS-MINLP && bash hpc/submit_recertify.sh [case ...]
set -euo pipefail
cd "$(dirname "$0")/.."
source hpc/pace_env.sh
# Memory and walltime are sized from MEASURED runs (see results/writeup.md). Observed totals over all
# variants: every case except 1354-pegase finishes in under 1 h (peak RSS 14.5 GB on 240-pserc); the
# unmerged 1354-pegase run took 9 h 50 m and 48.5 GB. Keeping the limits realistic matters because
# anything at or under 8 h runs on the FREE embers QOS (see pace_qos_for in hpc/pace_env.sh) -- so
# 17 of the 18 cases cost no credits at all. Only 1354-pegase needs the paid inferno QOS.
declare -A MEM=( [57-ieee-api]=48G [73-ieee-rts-api]=32G [89-pegase-api]=48G [118-ieee-api]=32G [179-goc-api]=32G [200-activ-api]=32G
                 [240-pserc-api]=48G [300-ieee-api]=64G [500-goc-api]=64G [1354-pegase-api]=96G )
declare -A TIM=( [1354-pegase-api]=24:00:00 )
: "${DEFAULT_TIME:=04:00:00}"     # ~4x the slowest observed non-1354 case; stays under the embers cap
CASES=("$@")
[ ${#CASES[@]} -gt 0 ] || CASES=(3-lmbd-api 5-pjm-api 14-ieee-api 24-ieee-rts-api 30-as-api 30-ieee-api 39-epri-api 57-ieee-api
    60-c-api 73-ieee-rts-api 89-pegase-api 118-ieee-api 179-goc-api 200-activ-api 300-ieee-api)
pace_check_account
for c in "${CASES[@]}"; do
    t="${TIM[$c]:-$DEFAULT_TIME}"; q="$(pace_qos_for "$t")"
    echo "[submit] $c  mem=${MEM[$c]:-16G} time=$t qos=$q account=$PACE_ACCOUNT"
    CASE="$c" sbatch --account="$PACE_ACCOUNT" --qos="$q" --job-name="mmrc-$c" --mem="${MEM[$c]:-16G}" --time="$t" \
        --export=ALL hpc/recertify.sbatch
done
