#!/bin/bash
# Environment for running the moment/SOS experiments on GT PACE (Phoenix). Sourced by setup and jobs.
#   source hpc/pace_env.sh
export PACE_ACCOUNT="${PACE_ACCOUNT:-gts-dmolzahn6}"
export MM_ROOT="${MM_ROOT:-$HOME/scratch/Moment-SoS-MINLP}"
# Julia depot on scratch (home quota is small); reuse the Julia 1.12.7 binary installed for the OTS project
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-$MM_ROOT/.julia}"
export JULIA="${JULIA:-$HOME/scratch/parameter_optimized_OTS/julia-1.12.7/bin/julia}"
# Precompile once on the login node for a portable CPU target, identical in every job
export JULIA_CPU_TARGET="${JULIA_CPU_TARGET:-generic;x86-64-v3,clone_all}"
# MOSEK finds its license at ~/mosek/mosek.lic (default location); Ipopt uses HSL MA97 (see setup_pace.sh);
# sync code with rsync --exclude Manifest.toml (the workstation Manifest has a Mac path for HSL_jll).
# NEVER rsync with --delete: the Julia depot lives at $MM_ROOT/.julia, inside the synced tree, and --delete
# wipes its artifacts/compiled/registries. Recovery is `bash hpc/setup_pace.sh` on the login node.
# BLAS single-threaded
# ---------------------------------------------------------------------------------------------
# Account and QOS policy.
#
# Accounts, in priority order:
#   1. gts-dmolzahn6  -- free monthly allocation; ALWAYS spend this first (it is the default above).
#      It refreshes every month, so a job held with pending reason AssocGrpBillingMinutes means
#      "drained for this month", NOT "exhausted" and NOT in need of a top-up request. Overflow for
#      the rest of the month, then switch back at the start of the next one.
#   2. gts-dmolzahn6-ece / gts-dmolzahn6-fy20phase3 -- paid; overflow only.
#
# QOS:
#   embers  -- FREE (UsageFactor 0.000000), but MaxWall 08:00:00, preemptible by inferno with
#              GraceTime 01:00:00, and NoReserve (backfill only, so queue wait is unpredictable).
#   inferno -- paid (UsageFactor 1.000000); no 8 h cap and not preempted.
# PACE_QOS=auto (default) picks embers whenever the requested walltime fits under the 8 h cap.
export PACE_QOS="${PACE_QOS:-auto}"          # auto | embers | inferno

# Echo the QOS to use for a walltime "HH:MM:SS" or "D-HH:MM:SS".
pace_qos_for() {
    if [ "${PACE_QOS}" != "auto" ]; then printf '%s' "${PACE_QOS}"; return; fi
    local t="$1" d=0 hh mm ss secs
    case "$t" in *-*) d="${t%%-*}"; t="${t#*-}" ;; esac
    IFS=: read -r hh mm ss <<<"$t"
    secs=$(( 10#${d:-0} * 86400 + 10#${hh:-0} * 3600 + 10#${mm:-0} * 60 + 10#${ss:-0} ))
    if [ "$secs" -le 28800 ]; then printf 'embers'; else printf 'inferno'; fi
}

# Credits still available on an account, from pace-quota (empty if unavailable, e.g. off-cluster).
pace_available() { pace-quota 2>/dev/null | awk -v a="$1" '$1 == a { print $NF }'; }

# Warn if the account being charged is drained for the month.
pace_check_account() {
    local avail; avail=$(pace_available "$PACE_ACCOUNT")
    [ -n "$avail" ] || return 0
    if awk -v v="$avail" 'BEGIN { exit !(v + 0 < 1) }'; then
        echo "[pace_env] WARNING: $PACE_ACCOUNT has only $avail credits available (drained for this month)." >&2
        echo "[pace_env]   Jobs will hold with reason AssocGrpBillingMinutes. For the rest of the month use" >&2
        echo "[pace_env]   PACE_ACCOUNT=gts-dmolzahn6-fy20phase3 (or -ece); switch back next month." >&2
    else
        echo "[pace_env] account $PACE_ACCOUNT: $avail credits available"
    fi
}

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
module purge 2>/dev/null || true
echo "[pace_env] root=$MM_ROOT julia=$JULIA depot=$JULIA_DEPOT_PATH account=$PACE_ACCOUNT"
