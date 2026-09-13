#!/bin/bash
# One production job for the whole dataset: dtOO meshes -> hex samples.
#
#   sbatch experimentell/hex3d_algohex/scripts/batch_generate.slurm
#   sbatch -t 24:00:00 .../batch_generate.slurm         # refill wants a long slot
#   GROW=0 sbatch -t 04:00:00 .../batch_generate.slurm  # single pass, no refill
#
# WHY cpu_il. The dev queues cap at 30 minutes, which fits exactly one wave of
# samples; every batch so far was a resubmit chain. This job is the escape
# hatch the cluster notes describe: cpu_il allows 3 days, so ONE queue wait
# covers the whole generation instead of a dozen 30-minute slots. The price is
# queue time (about a week, measured) -- so this is the job for a full or big
# generation, not for finishing a handful of missing samples. Those are still
# fastest as dev_cpu_il resubmits.
#
# WHY ONE NODE, for now. The measured optimum on 64 cores / 128 GB is 32
# workers with OMP_NUM_THREADS=2 (see the batch_samples.slurm header): roughly
# 80-100 samples per hour, so one 8-hour slot covers ~600-800 samples. More
# nodes would parallelize the samples linearly, but the algohex container runs
# with --root --rw from a shared unpacked rootfs, so concurrent nodes need one
# container copy each (plus a per-node share of the manifest). Deliberately not
# built yet; revisit when the target dataset outgrows a single node.
#
# REFILL LOOP. The job does not stop when the current samples are done: every
# round grows the Sobol prefix to twice its size (generate_machines.py --grow,
# new machines), exports their meshes, and then samples everything pending,
# until fewer than MIN_ROUND_SECONDS + DRAIN_SECONDS remain -- the queue never
# runs dry, so the slot keeps producing into its drain window. A round boundary
# costs the mesh export of the fresh machines (~10-25 min per 64 at PAR=32);
# that is the one idle window, left as is for now.
#
# The phase scripts are chained unchanged. Both are idempotent (a finished
# mesh.msh / sample.npz is skipped), each appends its per-job log under logs/,
# and their summary blocks -- the dtOO yield and the sample throughput -- are
# the measurements DATASET_PIPELINE.md hangs off.
#SBATCH --job-name=hex3d_generate
#SBATCH --output=hex3d_generate_%j.out
#SBATCH --partition=cpu_il
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G

set -u

WS="$(ws_find hex3d 2>/dev/null)"
REPO="${REPO:-${SLURM_SUBMIT_DIR:-$WS/domain_partition_3D}}"
[ -n "$WS" ] && [ -d "$REPO/experimentell/hex3d_algohex" ] || {
    echo "[gen] not the hex3d checkout: WS=$WS REPO=$REPO"; exit 1; }

LOGDIR="$REPO/experimentell/hex3d_algohex/logs"
mkdir -p "$LOGDIR"
LOG="$LOGDIR/batch_generate_${SLURM_JOB_ID:-manual}.log"
exec > >(tee -a "$LOG") 2>&1

cd "$REPO" || exit 1
E=experimentell/hex3d_algohex

# Thread budget and concurrency; both phase scripts read these from the
# environment. 32 x 2 is the measured optimum (48 x 1 produced fewer samples
# per slot), and sample_one.sh already forwards the pin into the container.
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
export PARALLEL="${PARALLEL:-32}"   # sample phase: xargs workers
export PAR="${PAR:-32}"             # dtOO phase: parallel exports
export NS="${NS:-2000 8000}"
export CANDS="${CANDS:-data/dataset/sobol/machine_*}"

# Refill controls. GROW=0 restores the old single pass; GROW_MAX_MACHINES caps
# the machine set (0 = unbounded) for anyone who wants a bounded dataset.
GROW="${GROW:-1}"
MIN_ROUND_SECONDS="${MIN_ROUND_SECONDS:-3600}"
DRAIN_SECONDS="${DRAIN_SECONDS:-1500}"
GROW_MAX_MACHINES="${GROW_MAX_MACHINES:-0}"
SAMPLER_JSON="data/dataset/sobol/sampler.json"

echo "[gen] host $(hostname), job ${SLURM_JOB_ID:-none}"
echo "[gen] NS='$NS' PAR=$PAR PARALLEL=$PARALLEL OMP_NUM_THREADS=$OMP_NUM_THREADS GROW=$GROW"

parse_wall() {
    case "$1" in
        [0-9]*-*:*)
            d=${1%%-*}; t=${1#*-}
            h=${t%%:*}; m=${t#*:}; m=${m%%:*}; s=${t##*:}
            echo $(( 10#${d:-0} * 86400 + 10#${h:-0} * 3600 + 10#${m:-0} * 60 + 10#${s:-0} )) ;;
        [0-9]*:*:*)
            h=${1%%:*}; r=${1#*:}; m=${r%%:*}; s=${r##*:}
            echo $(( 10#${h:-0} * 3600 + 10#${m:-0} * 60 + 10#${s:-0} )) ;;
        [0-9]*:*)
            m=${1%%:*}; s=${1##*:}
            echo $(( 10#${m:-0} * 60 + 10#${s:-0} )) ;;
        *) echo "" ;;
    esac
}

# The draw runs on the host, not in a container, so unlike the dtOO pass it
# needs scipy.stats -- probe candidates instead of assuming one.
PY_DRAW="${PY_DRAW:-}"
if [ -z "$PY_DRAW" ]; then
    for c in "$HOME/venvs/hydroflow/bin/python" python3; do
        if command -v "$c" >/dev/null 2>&1 && "$c" -c 'import scipy, yaml, numpy' >/dev/null 2>&1; then
            PY_DRAW="$c"; break
        fi
    done
fi
if [ -z "$PY_DRAW" ]; then
    echo "[gen] no python with scipy+yaml+numpy found -- refill disabled (set PY_DRAW)"
    GROW=0
else
    echo "[gen] python for draws: $PY_DRAW"
fi

round=0
while :; do
    LEFT="$(parse_wall "$(squeue -h -j "${SLURM_JOB_ID:-0}" -o %L 2>/dev/null || true)")"
    if [ -z "$LEFT" ]; then
        echo "[gen] no walltime from squeue -- running a single pass"
        GROW=0; LEFT=0
    fi

    if [ "$round" -gt 0 ]; then
        [ "$GROW" = 1 ] || break
        if [ "$LEFT" -le $(( MIN_ROUND_SECONDS + DRAIN_SECONDS )) ]; then
            echo "[gen] ${LEFT}s left (< $(( MIN_ROUND_SECONDS + DRAIN_SECONDS ))s) -- no next round"
            break
        fi
    fi

    round=$(( round + 1 ))
    echo "[gen] ================ round $round (${LEFT}s to walltime) ================"

    DREW=0
    if [ "$GROW" = 1 ] && [ "$LEFT" -gt $(( MIN_ROUND_SECONDS + DRAIN_SECONDS )) ]; then
        TARGET=0
        [ -f "$SAMPLER_JSON" ] && TARGET="$("$PY_DRAW" -c 'import json,sys; print(json.load(open(sys.argv[1])).get("count_target", 0))' "$SAMPLER_JSON" 2>/dev/null || echo 0)"
        if [ "$GROW_MAX_MACHINES" -gt 0 ] && [ "$TARGET" -ge "$GROW_MAX_MACHINES" ]; then
            echo "[gen] machine cap reached ($TARGET >= $GROW_MAX_MACHINES) -- sampling only"
        else
            echo "[gen] growing the machine set"
            if "$PY_DRAW" "$E/scripts/generate_machines.py" --preview --grow; then
                DREW=1
            else
                echo "[gen] draw failed -- refill disabled, sampling the pending queue"
                GROW=0
            fi
        fi
    elif [ "$GROW" = 1 ]; then
        echo "[gen] ${LEFT}s left -- sampling only, no new machines"
    fi

    echo "[gen] ---------------- round $round: dtOO meshes ----------------"
    bash "$E/scripts/batch_dtoo_export.slurm"
    RC_DTOO=$?
    echo "[gen] round $round: dtOO rc=$RC_DTOO"

    echo "[gen] ---------------- round $round: hex samples ----------------"
    bash "$E/scripts/batch_samples.slurm"
    RC_SAMPLES=$?
    echo "[gen] round $round: samples rc=$RC_SAMPLES"

    if [ "$RC_SAMPLES" -eq 0 ] && [ "$DREW" = 0 ]; then
        echo "[gen] nothing pending and no new machines -- done"
        break
    fi
done

echo "[gen] ================ generation loop complete ================"
echo "[gen] rounds run:     $round"
READY=$(find output/hex3d_algohex/batch -name sample.npz -path '*machine_*' 2>/dev/null | wc -l)
echo "[gen] samples ready:  $READY"
echo "[gen] leftovers at the drain are normal; the phase summaries carry the numbers"
echo "[gen] log:            $LOG"
# Let the tee drain before SLURM sees the job as done. A bare wait deadlocks
# here (shell holds the process-substitution pipe open, tee never sees EOF --
# job 6890146 idled 3.5 h on this line); closing stdout/stderr first fixes it.
exec >/dev/null 2>&1
wait
exit 0
