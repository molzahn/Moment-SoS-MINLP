#!/bin/bash
# Submit the larger DC-OTS paper cases (one job each), with memory/time scaled to case size.
#   cd ~/scratch/Moment-SoS-MINLP && bash hpc/submit_large.sh
set -euo pipefail
cd "$(dirname "$0")/.."
source hpc/pace_env.sh
submit() {  # case mem time
    CASE="$1" sbatch --account="$PACE_ACCOUNT" --job-name="mm-$1" --mem="$2" --time="$3" --export=ALL hpc/experiment3.sbatch
}
submit 179-goc-api     64G  36:00:00
submit 200-activ-api   64G  36:00:00
submit 240-pserc-api   64G  36:00:00
submit 300-ieee-api   128G  48:00:00
submit 500-goc-api    192G  72:00:00
submit 1354-pegase-api 384G 96:00:00
