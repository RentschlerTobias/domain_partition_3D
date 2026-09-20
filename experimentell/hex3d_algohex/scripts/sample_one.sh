#!/usr/bin/env bash
# One sample, end to end. The unit of work that batch_samples.slurm runs many
# of in parallel.
#
#   sample_one.sh <tet.vtk> <n> <outdir> [params.json]
#
# AlgoHex -> clean_blocks -> tfi -> export_sample, all logs inside <outdir>.
# Writes <outdir>/sample.npz and nothing outside <outdir>.
#
# IDEMPOTENT: if <outdir>/sample.npz already exists, it exits 0 immediately.
# That is what makes a killed batch job cheap to resubmit -- the survivors are
# skipped in milliseconds.
#
# Exit codes are the batch summary's data: 0 ok, 2 algohex, 3 clean_blocks,
# 4 tfi, 5 export, 6 bad arguments, 7 pruned (sample budget or drain expired).
#
# Every stage chain of one sample runs under SAMPLE_TIMEOUT (default 2400 s,
# measured samples take 893-1113 s): when the budget or the drain deadline is
# gone, the whole process group is killed and the lane moves to the next
# sample -- a hung sample must not pin a lane for the rest of the slot. The
# pruned sample is retried by a later run; nothing checkpoints mid-sample.
# One row per finished attempt (ok/fail/prune, per-stage seconds) is appended
# to <OUT_ROOT>/sample_history.tsv.

set -u

TET="${1:?usage: sample_one.sh <tet.vtk> <n> <outdir> [params.json]}"
N="${2:?missing -n}"
D="${3:?missing outdir}"
PARAMS="${4:-}"

REPO="${REPO:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}"
E="$REPO/experimentell/hex3d_algohex"
PY="${PY:-$REPO/.venv/bin/python}"
H="${TARGET_H:-0.05}"
IMG="${ALGOHEX_IMAGE:-algohex}"
BACKEND="${BACKEND:-enroot}"
# Keep the AlgoHex .ovm? Decision J says no at dataset scale: 13 MB per sample
# at n=60000 is 130 GB across 10 000 samples. Worth keeping for a showcase set.
KEEP_OVM="${KEEP_OVM:-0}"

# One aggregate budget per sample (measured 893-1113 s at 32 lanes). When it
# -- or the drain deadline exported by batch_samples.slurm -- is gone, the
# watchdog kills the whole stage process group and the lane moves on.
SAMPLE_TIMEOUT="${SAMPLE_TIMEOUT:-2400}"

# One thread per worker. Measured: a single sample with 8 cores allocated used
# 1.45 cores on average, but three samples on a 64-core node used 6.4 each --
# numpy/scipy expand to whatever the node offers. Unpinned, 32 concurrent
# workers would ask for 205 cores of 64 and spend the difference on context
# switching. The parallelism that matters here is across samples, not inside
# one: AlgoHex itself is largely serial.
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"

[ -f "$TET" ] || { echo "[$D] no such tet mesh: $TET"; exit 6; }
mkdir -p "$D" || exit 6

if [ -f "$D/sample.npz" ]; then
    [ "${VERBOSE:-0}" = 1 ] && echo "[$(basename "$D")] already done, skipping"
    exit 0
fi

t0=$(date +%s)
name=$(basename "$D")
log() { echo "[$name] $*"; }
HIST="${OUT_ROOT:-$(dirname "$D")}/sample_history.tsv"
STAGE_TIMES=""

# The lane must not start a sample it cannot finish: batch_samples.slurm
# exports DRAIN_DEADLINE (walltime minus DRAIN_SECONDS); past it the lane
# stops claiming work so the batch can still print its summary instead of
# being killed mid-sample by the walltime.
if [ "${DRAIN_DEADLINE:-0}" -gt 0 ] && [ "$(date +%s)" -ge "${DRAIN_DEADLINE:-0}" ]; then
    log "drain: not starting, $(( ${JOB_END:-0} - $(date +%s) ))s to walltime"
    exit 0
fi

# Remaining budget for this sample: the sample-wide cap, clamped by the job
# end so a sample killed at the walltime still leaves a history row.
stage_left() {
    left=$(( SAMPLE_TIMEOUT - ( $(date +%s) - t0 ) ))
    if [ "${JOB_END:-0}" -gt 0 ]; then
        jleft=$(( JOB_END - $(date +%s) - 15 ))
        [ "$jleft" -lt "$left" ] && left=$jleft
    fi
    echo "$left"
}

# run_stage <label> <logfile> <cmd...>. setsid puts the stage in its own
# process group so the watchdog kills python AND the enroot/HexMeshing
# children together -- killing only the direct child would leave the real
# compute running orphaned. Sets RC; appends "<label>=<n>s" to STAGE_TIMES.
run_stage() {
    label="$1"; logfile="$2"; shift 2
    left=$(stage_left)
    s0=$(date +%s)
    if [ "$left" -le 0 ]; then
        RC=124; STAGE_TIMES="$STAGE_TIMES $label=0s"; return
    fi
    setsid "$@" > "$logfile" 2>&1 < /dev/null &
    pid=$!
    # The watchdog runs in its own process group too: killing the guard shell
    # alone would orphan its pending `sleep`, which inherits (and holds) the
    # log pipe -- the batch's final `wait` would then block on a tee that never
    # sees EOF (the 181ee81 deadlock, measured here). Kill the group.
    setsid bash -c 'sleep "$1"; kill -TERM -- "-$2" 2>/dev/null; sleep 10; kill -KILL -- "-$2" 2>/dev/null' _ "$left" "$pid" &
    guard=$!
    wait "$pid"; RC=$?
    kill -TERM -- "-$guard" 2>/dev/null || true
    wait "$guard" 2>/dev/null || true
    STAGE_TIMES="$STAGE_TIMES $label=$(( $(date +%s) - s0 ))s"
}

# 124 = budget spent before the stage started, >=128 = killed by the watchdog.
pruned() { [ "$1" -eq 124 ] || [ "$1" -ge 128 ]; }

# record <status> <stage> <rc>: one TSV row per finished attempt. Single
# O_APPEND writes from the lanes do not interleave.
record() {
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$(date -Is)" "$name" "$N" "$1" "$2" "$3" "$STAGE_TIMES" \
        >> "$HIST" 2>/dev/null || true
}

# stage_fail <label> <script-exit-code> <logfile>: logs, records, exits.
stage_fail() {
    label="$1"; code="$2"; logfile="$3"
    if pruned "$RC"; then
        log "PRUNED $label $STAGE_TIMES (rc=$RC), see $logfile"
        record prune "$label" "$RC"; exit 7
    fi
    log "FAIL $label (rc=$RC), see $logfile"
    record fail "$label" "$RC"; exit "$code"
}

log "start  tet=$(basename "$TET") n=$N h=$H"

run_stage algohex "$D/algohex.log" \
    "$PY" "$E/run_algohex.py" --backend "$BACKEND" --image "$IMG" \
    --prefix hex --tag "$name" --in-vtk "$TET" --out-dir "$D" \
    -- -n "$N"
[ "$RC" -eq 0 ] || stage_fail algohex 2 "$D/algohex.log"

run_stage clean_blocks "$D/clean_blocks.log" \
    "$PY" "$E/clean_blocks.py" "$D/hex_hex_${name}.ovm" --input-vtk "$TET" \
    --collapse-rounds 5 --untangle --untangle-rounds 6 \
    --out "$D/blocks.vtk"
[ "$RC" -eq 0 ] || stage_fail clean_blocks 3 "$D/clean_blocks.log"

# --require-lattices: a block without a lattice pins whole direction classes
# through the complex, and the refill then comes out distorted. For training
# data a skipped sample is cheaper than a distorted one (generate_dataset.sh
# has the measurement on cand_000).
run_stage tfi "$D/tfi.log" \
    "$PY" "$E/tfi.py" "$D/blocks.vtk" --input-vtk "$TET" --target-h "$H" \
    --require-lattices --solve-divisions
[ "$RC" -eq 0 ] || stage_fail tfi 4 "$D/tfi.log"

run_stage export "$D/export_sample.log" \
    "$PY" "$E/export_sample.py" "$D/blocks.vtk" --tet "$TET" --target-h "$H" \
    --n "$N" ${PARAMS:+--params "$PARAMS"} \
    --out "$D/sample.npz" --check
[ "$RC" -eq 0 ] || stage_fail export 5 "$D/export_sample.log"

# Retention, decision J: the sample plus blocks.vtk are the recovery boundary,
# everything else is derivable from them in seconds.
[ "$KEEP_OVM" = 1 ] || rm -f "$D"/*.ovm
rm -f "$D"/blocks_nfaces.vtk "$D"/blocks_quality.vtk "$D"/blocks.msh

sz=$(stat -c %s "$D/sample.npz" 2>/dev/null || echo 0)
total=$(( $(date +%s) - t0 ))
STAGE_TIMES="$STAGE_TIMES total=${total}s"
record ok done 0
log "OK     ${total}s, sample $(( sz / 1024 )) KB, kept $(du -sk "$D" | cut -f1) KB total [$STAGE_TIMES]"
