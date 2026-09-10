#!/usr/bin/env bash
# One command, from nothing to a verified sample on bwUniCluster 3.0.
#
#   bash cluster_setup.sh --dry-run     # check everything, submit nothing
#   bash cluster_setup.sh               # do it
#
# Run this ON THE LOGIN NODE. It is idempotent: every step checks whether it is
# already done, so re-running after a failure resumes rather than restarts.
#
# What it does, in order:
#   1. preflight   sbatch/enroot/ws_* present, partition names, python module
#   2. workspace   ws_allocate if missing
#   3. repo        clone, and restore the committed fixtures
#   4. python      module load + venv + numpy/scipy/meshio/gmsh  (NOT torch)
#   5. sources     build_algohex_enroot.sh -- downloads, needs this login node
#   6. submit      build -> smoke -> pipeline, chained with --dependency=afterok
#
# Nothing is copied from another machine. AlgoHex is built from its pinned
# public commit, the base image comes from Docker Hub, the one test input comes
# out of the fixtures archive in git.
#
# After step 6 you can log out. Three jobs run in sequence; the last one writes
# a sample and round-trips it.

set -uo pipefail

DRY=0
[ "${1:-}" = "--dry-run" ] && DRY=1

WS_NAME="${WS_NAME:-hex3d}"
WS_DAYS="${WS_DAYS:-60}"
REPO_URL="${REPO_URL:-git@github.com:RentschlerTobias/domain_partition_3D.git}"
# The default branch is master and none of this exists there.
REPO_BRANCH="${REPO_BRANCH:-data_generation}"
PYTHON_MODULE="${PYTHON_MODULE:-devel/python/3.13.3-gnu-14.2}"
# eigenfrequencies/docs/cluster-resource-sizing.md measured the three
# partitions -- cpu_il 64 cores/256 GiB, cpu 96/384, dev_cpu_il 64/256 with
# 30-minute slots -- and recommends cpu_il on throughput.
#
# EVERYTHING HERE RUNS ON dev_cpu_il ANYWAY, and that overrides the doc:
# cpu_il measures about a WEEK of queue time in practice, which makes it the
# slowest route to a first result no matter how many cores it has. So the build
# is cut into three resumable stages that each aim at a 30-minute slot, and
# the tests are small by construction.
#
# Escape hatch if a stage will not fit: BUILD_PARTITION=cpu_il with
# STAGE=all and a 4 h walltime, then wait out the queue.
DEV_PARTITION="${DEV_PARTITION:-dev_cpu_il}"
BUILD_PARTITION="${BUILD_PARTITION:-$DEV_PARTITION}"
PIPE_PARTITION="${PIPE_PARTITION:-$DEV_PARTITION}"

say()  { printf '\n\033[1m[setup] %s\033[0m\n' "$*"; }
info() { printf '        %s\n' "$*"; }
die()  { printf '\n[setup] STOP: %s\n' "$*" >&2; exit 1; }
run()  { if [ "$DRY" = 1 ]; then info "would run: $*"; else "$@"; fi; }

# ---------------------------------------------------------------- 1 preflight
say "1/6 preflight"
for c in sbatch squeue enroot; do
    command -v "$c" >/dev/null || die "$c not on PATH -- is this a login node?"
done
info "sbatch, squeue, enroot: present"

# sinfo is denied on this cluster (measured 2026-09-11: "slurm_load_node:
# Access/permission denied"), so partitions are probed differently and, if that
# is denied too, left to `sbatch --test-only` to reject.
if scontrol show partition 2>/dev/null | grep -q PartitionName; then
    info "partitions visible:"
    scontrol show partition 2>/dev/null \
        | grep -oE "PartitionName=[^ ]+|MaxTime=[^ ]+" | paste - - | sed 's/^/        /'
    for p in "$BUILD_PARTITION" "$DEV_PARTITION"; do
        scontrol show partition "$p" >/dev/null 2>&1 \
            || info "WARNING: partition '$p' not found -- set BUILD_PARTITION / DEV_PARTITION"
    done
else
    info "scontrol denied too. Partition names stay UNVERIFIED:"
    info "  build='$BUILD_PARTITION'  dev='$DEV_PARTITION'"
    info "  --test-only below will reject them if they are wrong."
fi

module load "$PYTHON_MODULE" 2>/dev/null \
    && info "module $PYTHON_MODULE loaded" \
    || info "WARNING: could not load $PYTHON_MODULE (module avail python)"

# --------------------------------------------------------------- 2 workspace
say "2/6 workspace"
WS="$(ws_find "$WS_NAME" 2>/dev/null)"
if [ -z "$WS" ]; then
    run ws_allocate "$WS_NAME" "$WS_DAYS"
    WS="$(ws_find "$WS_NAME" 2>/dev/null)"
    [ "$DRY" = 1 ] && WS="<workspace>"
fi
[ -n "$WS" ] || die "no workspace '$WS_NAME'"
info "WS=$WS"
export WS
run mkdir -p "$WS/enroot-images" "$WS/enroot-data" "$WS/enroot-cache"

export ENROOT_DATA_PATH="$WS/enroot-data"
export ENROOT_CACHE_PATH="$WS/enroot-cache"
export ENROOT_TEMP_PATH="/tmp/${USER:-$(id -un)}-enroot"   # LOCAL disk, not the workspace
export ENROOT_SQUASH_OPTIONS="-comp zstd -noD"   # never lzo
run mkdir -p "$ENROOT_TEMP_PATH"

# -------------------------------------------------------------------- 3 repo
say "3/6 repo and fixtures"
REPO="$WS/domain_partition_3D"
if [ -d "$REPO/.git" ]; then
    info "already cloned: $REPO"
else
    # --branch: the repo's default is master, and NONE of these scripts exist
    # there. A plain clone lands on master and every path below is missing.
    run git clone --branch "$REPO_BRANCH" "$REPO_URL" "$REPO" \
        || die "clone failed. SSH key on this node? Or set REPO_URL to the https form."
fi
E="experimentell/hex3d_algohex"
if [ "$DRY" = 0 ]; then
    cd "$REPO" || die "cannot cd $REPO"
    # Checked, not just printed. Two ways to be on the wrong commit: a clone
    # that defaulted to master, or a detached HEAD from
    # `git checkout origin/<branch>` -- which prints "HEAD" here and makes a
    # later `git pull` fail.
    BR="$(git rev-parse --abbrev-ref HEAD)"
    info "branch: $BR"
    if [ "$BR" = "HEAD" ]; then
        info "detached HEAD -- attaching to $REPO_BRANCH"
        git switch "$REPO_BRANCH" || die "cannot switch to $REPO_BRANCH"
    elif [ "$BR" != "$REPO_BRANCH" ]; then
        info "on '$BR', switching to '$REPO_BRANCH'"
        git switch "$REPO_BRANCH" || die "cannot switch to $REPO_BRANCH"
    fi
    [ -f "$E/scripts/build_algohex.slurm" ] \
        || die "$E/scripts/build_algohex.slurm missing -- wrong branch?"
    # 22.4 MiB of block structures and the one tet input, 5.5 MiB compressed
    [ -f data/T1_9/T1_9_tet_v5.vtk ] || tar -I zstd -xf "$E/fixtures/blocks_core.tar.zst"
    info "fixtures: $(ls data/T1_9/T1_9_tet_v5.vtk 2>/dev/null || echo MISSING)"
fi

# ------------------------------------------------------------------ 4 python
say "4/6 python environment"
if [ "$DRY" = 0 ] && [ ! -x "$REPO/.venv/bin/python" ]; then
    python -m venv "$REPO/.venv" || die "venv failed"
    "$REPO/.venv/bin/pip" -q install --upgrade pip
    # Only what the dataset chain imports. torch is deliberately absent:
    # reattach.py is the only stage that needs it, it makes the full.* CFD
    # variant that decision J deletes, and SKIP_REATTACH=1 skips it.
    "$REPO/.venv/bin/pip" -q install \
        numpy==2.4.4 scipy==1.17.1 meshio==5.3.5 gmsh==4.15.2 \
        || die "pip install failed"
fi
if [ "$DRY" = 0 ]; then
    "$REPO/.venv/bin/python" -c \
        "import numpy,scipy,meshio;print('        numpy',numpy.__version__,
'scipy',scipy.__version__,'meshio',meshio.__version__)" || die "env incomplete"
else
    info "would create $REPO/.venv with numpy, scipy, meshio, gmsh"
fi

# ----------------------------------------------------------------- 5 sources
say "5/6 AlgoHex sources into the build container (downloads, ~10-20 min)"
if enroot list 2>/dev/null | grep -qxF algohex-build; then
    info "container 'algohex-build' exists, skipping"
elif enroot list 2>/dev/null | grep -qxF algohex; then
    info "container 'algohex' already built, skipping the build entirely"
else
    run bash "$REPO/$E/scripts/build_algohex_enroot.sh"
fi

# ------------------------------------------------------------------ 6 submit
say "6/6 submit the chain"
S="$REPO/$E/scripts"
# stdout carries ONLY the job id, because the caller captures it. Everything
# human-readable goes to stderr -- otherwise --test-only's output ends up
# inside the next job's --dependency argument, which is what the first dry run
# of this script did.
sub() {  # sub <script> <extra sbatch args...> -> job id on stdout
    local script="$1"; shift
    if [ "$DRY" = 1 ]; then
        { sbatch --test-only "$@" "$script" 2>&1 | sed 's/^/        /'; } >&2
        return 0
    fi
    sbatch --parsable "$@" "$script"
}

# Three resumable stages instead of one 4-hour job, so the whole build fits
# dev_cpu_il's 30-minute slots. Chained afterok; each stage exits in seconds if
# its artifact already exists, so a re-run of this script is cheap.
J1=""
if enroot list 2>/dev/null | grep -qxF algohex; then
    info "image already built; skipping the build jobs"
else
    for st in ipopt bonmin algohex; do
        JN=$(STAGE="$st" sub "$S/build_algohex.slurm" \
                --partition="$BUILD_PARTITION" --export=ALL,STAGE="$st" \
                --job-name="algohex_$st" ${J1:+--dependency=afterok:$J1})
        info "build:$(printf '%-8s' "$st") job ${JN:-(dry run)}   (30 min slot, resumable)"
        J1="$JN"
    done
fi

# A dependency on an empty id would be malformed, and in a dry run there are no
# ids at all.
J2=$(sub "$S/smoke_algohex.slurm" --partition="$DEV_PARTITION" \
        ${J1:+--dependency=afterok:$J1})
info "smoke      job ${J2:-(dry run)}   (5 min: image readable on a compute node)"

J3=$(sub "$S/smoke_pipeline.slurm" --partition="$PIPE_PARTITION" \
        ${J2:+--dependency=afterok:$J2})
info "pipeline   job ${J3:-(dry run)}   (30 min: one sample, round-trip gate)"

if [ "$DRY" = 1 ]; then
    say "dry run done. Nothing was submitted, nothing was downloaded."
    exit 0
fi

cat <<EOF

[setup] submitted. You can log out.

  watch:    squeue --me
  logs:     tail -f $REPO/algohex_build_${J1:-\*}.out
  result:   $REPO/output/hex3d_algohex/smoke2000/sample.npz

The last job prints "PASS -- the cluster produced a verified sample" on
success. If a job fails, re-running this script resumes from where it stopped.
EOF
