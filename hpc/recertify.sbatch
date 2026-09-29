#!/bin/bash
#SBATCH --job-name=mm-recert
#SBATCH --output=hpc/logs/%x_%j.out
#SBATCH --error=hpc/logs/%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#
# One job per case (memory and walltime differ by case). Submit with hpc/submit_recertify.sh, or e.g.
#   CASE=300-ieee-api sbatch --account=$PACE_ACCOUNT --mem=48G --time=16:00:00 hpc/recertify.sbatch
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
source hpc/pace_env.sh
: "${CASE:?set CASE, e.g. CASE=179-goc-api}"
export MOSEK_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export CERT_OPT_TIME="${CERT_OPT_TIME:-1800}"
export MAX_SDP_TIME="${MAX_SDP_TIME:-43200}"      # 12 h per SDP
export MERGE_ZMAX="${MERGE_ZMAX:-1e-3}"
echo "case $CASE on $(hostname), cpus $MOSEK_THREADS, mem ${SLURM_MEM_PER_NODE:-?} MB, certificate optimization time $CERT_OPT_TIME s, start $(date)"
/usr/bin/time -v "$JULIA" --project=. scripts/recertify.jl "$CASE"
echo "finished $(date)"
