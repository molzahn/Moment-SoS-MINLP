#!/bin/bash
#SBATCH --job-name=mm5-cardinality
#SBATCH --output=hpc/logs/%x_%j.out
#SBATCH --error=hpc/logs/%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#
# Experiment 5 (cardinality sweep for rounding). No SDP is solved -- it reuses the marginals stored by
# experiments 3/4 and only evaluates configurations with Ipopt -- so it is cheap and fits the free
# `embers` QOS. Submit with hpc/submit_exp5.sh.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
source hpc/pace_env.sh
: "${CASE:?set CASE, e.g. CASE=300-ieee-api}"
export MOSEK_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export N="${N:-100}"
export KMAX_SET="${KMAX_SET:-1,2,3,4,5,6,8}"
export POLISH_EVALS="${POLISH_EVALS:-0}"
# Must match the run that produced the saved marginals (experiments 3/4 used the 1e-3 default).
# Rounding does not need the rigorous model: candidates are evaluated on the original network either way.
export MERGE_ZMAX="${MERGE_ZMAX:-1e-3}"
echo "case $CASE on $(hostname), cpus $MOSEK_THREADS, N=$N, kmax set $KMAX_SET, merge_zmax $MERGE_ZMAX, start $(date)"
/usr/bin/time -v "$JULIA" --project=. scripts/experiment5_cardinality.jl "$CASE"
echo "finished $(date)"
