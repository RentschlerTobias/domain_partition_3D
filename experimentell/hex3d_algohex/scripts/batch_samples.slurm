#!/bin/bash
# T10: many samples in ONE job, concurrently. The throughput measurement.
#
#   sbatch --test-only .../batch_samples.slurm
#   sbatch            .../batch_samples.slurm                 # the -n sweep test
#   NS="2000 8000" CANDS="data/random_tistos/investigated/cand_00*" \
#       sbatch .../batch_samples.slurm                        # the real thing
#
# ---------------------------------------------------------------------------
# WHY ONE JOB AND NOT A SLURM ARRAY
#
# Two measurements decide the shape, both from 2026-09-11:
#
#   1. dev_cpu_il allows ~4 QUEUED jobs per user, shared with everything else
#      you run (QOSMaxSubmitJobPerUserLimit). An array of 30 tasks does not fit
#      a budget of 4.
#   2. AlgoHex ran at 18 % efficiency on 8 cores -- about 1.5 cores actually in
#      use -- and 1.87 GB at -n 2000 (job 6860177). One task per node would
#      leave 62 of 64 cores idle.
#
# CORES_PER_SAMPLE, measured and then corrected. With 8 cores allocated a
# single sample averaged 1.45; with 64 cores three samples averaged 6.4 each
# (job 6860772), because numpy/scipy expand into whatever is free.
# sample_one.sh now pins OMP/OPENBLAS/MKL/NUMEXPR to 1 thread, which is what
# makes the per-sample cost predictable enough to compute a concurrency from.
#
# So: one job, one queue slot, and `xargs -P` inside it. On a 64-core/256 GiB
# node that is ~30 concurrent samples, which is the whole point.
#
# CONCURRENCY IS BOUND BY MEMORY, not by cores. Two concurrent AlgoHex runs at
# -n 60000 (~6 GB each) are what killed v7 and v8 on a 7.7 GiB box. The
# defaults below are per-sample measurements; -n changes them a lot:
#
#   -n 2000    1.87 GB, 12:25   <- measured, job 6860177
#   -n 60000   ~6 GB,   ~55 min <- measured on the VPS, does NOT fit a 30 min slot
#
# WALLTIME. 30 minutes is dev_cpu_il's cap and one wave at -n 2000 takes 12-13,
# so ONE wave fits and two do not. sample_one.sh is idempotent, so a job killed
# mid-wave loses only the samples that were still running: resubmit and the
# finished ones are skipped in milliseconds.
# ---------------------------------------------------------------------------
#SBATCH --job-name=hex3d_batch
#SBATCH --output=hex3d_batch_%j.out
#SBATCH --partition=dev_cpu_il
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
# 128G, not 64: at 2.5 GB per sample, 64G caps concurrency at 26 while the
# cores allow 32. 128G moves the binding constraint to cores, which is where
# it belongs, and dev_cpu_il has 256 GiB.
#SBATCH --mem=128G

set -u

# ws_find, NOT the environment: an interactive shell that sourced the
# eigenfrequencies cluster_env.sh exports WS=<eigenfreq workspace> and every
# sbatch inherits it -- job 6877668 died on cd "<eigenfreq>/domain_partition_3D".
WS="$(ws_find hex3d 2>/dev/null)"
REPO="${REPO:-${SLURM_SUBMIT_DIR:-$WS/domain_partition_3D}}"
[ -n "$WS" ] && [ -d "$REPO/experimentell/hex3d_algohex" ] || {
    echo "[batch] not the hex3d checkout: WS=$WS REPO=$REPO"; exit 1; }

# Keep the batch log in the repository: the SLURM .out in the submit directory
# is not durable, and the throughput block at the end is a T10 measurement we
# want committed. Mirrors batch_dtoo_export.slurm.
LOGDIR="$REPO/experimentell/hex3d_algohex/logs"
mkdir -p "$LOGDIR"
LOG="$LOGDIR/batch_samples_${SLURM_JOB_ID:-manual}.log"
exec > >(tee -a "$LOG") 2>&1

PYTHON_MODULE="${PYTHON_MODULE:-devel/python/3.13.3-gnu-14.2}"
NS="${NS:-500 2000 8000}"          # the -n values; each becomes its own sample
CANDS="${CANDS:-}"                 # candidate dirs; empty = the fixtures tet
OUT_ROOT="${OUT_ROOT:-output/hex3d_algohex/batch}"

module load "$PYTHON_MODULE" 2>/dev/null || true
export PATH="$REPO/.venv/bin:$PATH"
export REPO
export ENROOT_DATA_PATH="${ENROOT_DATA_PATH:-$WS/enroot-data-$(hostname -s)}"
export ENROOT_TEMP_PATH="${ENROOT_TEMP_PATH:-/tmp/${USER:-$(id -un)}-enroot}"
export ENROOT_RUNTIME_PATH="${ENROOT_RUNTIME_PATH:-/tmp/${USER:-$(id -un)}-enroot-run}"
mkdir -p "$ENROOT_TEMP_PATH" "$ENROOT_RUNTIME_PATH"

# enroot start uses the NAME form ('algohex', see run_algohex.py): the
# container must exist unpacked under this node's ENROOT_DATA_PATH. Race-safe
# install (private name + mv -T), mirrors batch_dtoo_export.slurm.
CONTAINER="${CONTAINER:-algohex}"
SQSH="${SQSH:-$WS/enroot-images/algohex.sqsh}"
if ! enroot list 2>/dev/null | grep -qxF "$CONTAINER"; then
    tmp_name="$CONTAINER.${SLURM_JOB_ID:-$$}"
    echo "[batch] unpacking $SQSH -> $ENROOT_DATA_PATH/$CONTAINER"
    enroot create --name "$tmp_name" "$SQSH" \
        || { echo "[batch] enroot create failed: $SQSH"; exit 1; }
    if mv -T "$ENROOT_DATA_PATH/$tmp_name" "$ENROOT_DATA_PATH/$CONTAINER" 2>/dev/null; then
        echo "[batch] container $CONTAINER installed"
    else
        rm -rf "$ENROOT_DATA_PATH/$tmp_name"
        enroot list 2>/dev/null | grep -qxF "$CONTAINER" \
            || { echo "[batch] container install mismatch: $ENROOT_DATA_PATH/$CONTAINER"; exit 1; }
    fi
else
    echo "[batch] container $CONTAINER already unpacked"
fi

cd "$REPO" || exit 1
E=experimentell/hex3d_algohex
mkdir -p "$OUT_ROOT"

# --------------------------------------------------------------- concurrency
CORES="${SLURM_CPUS_PER_TASK:-$(nproc)}"
MEM_MB="${SLURM_MEM_PER_NODE:-131072}"
# Per-sample figures. MEM_PER_SAMPLE_MB defaults to the -n 2000 measurement
# plus headroom; override it upward for large -n or the node will OOM.
MEM_PER_SAMPLE_MB="${MEM_PER_SAMPLE_MB:-2500}"
# 2 with threads pinned to 1: AlgoHex is largely serial, plus some overlap.
# Without the pinning in sample_one.sh this would have to be 7.
CORES_PER_SAMPLE="${CORES_PER_SAMPLE:-2}"
P_MEM=$(( MEM_MB / MEM_PER_SAMPLE_MB ))
P_CPU=$(( CORES / CORES_PER_SAMPLE ))
P="${PARALLEL:-$(( P_MEM < P_CPU ? P_MEM : P_CPU ))}"
[ "$P" -ge 1 ] || P=1

echo "[batch] host $(hostname), job ${SLURM_JOB_ID:-none}"
echo "[batch] $CORES cores, ${MEM_MB} MB"
echo "[batch] concurrency $P   (memory allows $P_MEM, cores allow $P_CPU)"

# ------------------------------------------------- walltime budget & drain
# Lanes must not start work they cannot finish: at 32 lanes a sample takes
# 893-1113 s, so the last minutes of a slot would only produce samples killed
# at the walltime. DRAIN_SECONDS before the end every lane stops claiming work
# (sample_one.sh checks DRAIN_DEADLINE); two minutes before the end a snapshot
# is written into the log, because the summary below never prints when SLURM
# kills the job at the limit (job 6890146 idled into exactly that kill).
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
            # squeue prints minutes:seconds for a remaining time below one
            # hour -- i.e. always on the 30-minute dev partition.
            m=${1%%:*}; s=${1##*:}
            echo $(( 10#${m:-0} * 60 + 10#${s:-0} )) ;;
        *) echo "" ;;
    esac
}
DRAIN_SECONDS="${DRAIN_SECONDS:-1500}"
JOB_END=0; DRAIN_DEADLINE=0; REPORT_AT=0; LEFT=""
if [ -n "${SLURM_JOB_ID:-}" ]; then
    LEFT=$(parse_wall "$(squeue -h -j "$SLURM_JOB_ID" -o %L 2>/dev/null || true)")
fi
if [ -n "$LEFT" ] && [ "$LEFT" -gt 0 ]; then
    JOB_END=$(( $(date +%s) + LEFT ))
    DRAIN_DEADLINE=$(( JOB_END - DRAIN_SECONDS ))
    REPORT_AT=$(( JOB_END - 120 ))
    echo "[batch] walltime left ${LEFT}s: no new samples in the last ${DRAIN_SECONDS}s, snapshot 120s before the end"
else
    echo "[batch] no walltime from squeue -- drain and early snapshot disabled"
fi
export SAMPLE_TIMEOUT="${SAMPLE_TIMEOUT:-2400}" JOB_END DRAIN_DEADLINE

# ------------------------------------------------------------------ manifest
# One line per sample: <tet> <n> <outdir> [params]
# Parallel-safe: append-only (no truncation race) and every candidate is
# claimed atomically through the shared pool (mkdir), so several
# batch_generate jobs never sample the same machine. A claim older than
# CLAIM_TTL is treated as stale (crashed job) and taken over.
MAN="$OUT_ROOT/manifest.txt"
touch "$MAN"
GEN_POOL="${GEN_POOL:-$HOME/hex3d_pool}"
CLAIM_TTL_SECONDS="${CLAIM_TTL_SECONDS:-14400}"
mkdir -p "$GEN_POOL/claimed"
claim_machine() {
    local key="$1"
    local dir="$GEN_POOL/claimed/$key"
    if ! mkdir "$dir" 2>/dev/null; then
        # owned: stale only if the claim is older than the TTL
        if [ -n "$(find "$dir" -maxdepth 0 -mmin +$(( CLAIM_TTL_SECONDS / 60 )))" ]; then
            return 1
        fi
        rm -rf "$dir" 2>/dev/null
        mkdir "$dir" 2>/dev/null || return 1
    fi
    date -Is > "$dir/claim"
    return 0
}

if [ -z "$CANDS" ]; then
    # The -n sweep on the committed fixture: no geometry transfer needed, and
    # T6 established that different -n give genuinely different structures
    # rather than duplicates, so these are real samples and not a rehearsal.
    TET=data/T1_9/T1_9_tet_v5.vtk
    [ -f "$TET" ] || tar -I zstd -xf "$E/fixtures/blocks_core.tar.zst"
    for n in $NS; do
        echo "$TET $n $OUT_ROOT/T1_9_n$n" >> "$MAN"
    done
else
    for cand in $CANDS; do
        name=$(basename "$cand")
        if ! claim_machine "$name"; then
            echo "[batch] skip $name: claimed by another job"
            continue
        fi
        msh="$cand/mesh.msh"
        [ -f "$msh" ] || { echo "[batch] skip $name: no mesh.msh"; continue; }
        d0="$OUT_ROOT/$name"
        mkdir -p "$d0"
        # tet_prep is cheap and serial per candidate; do it once here rather
        # than once per -n, since every -n of a candidate shares the tet mesh.
        if [ ! -f "$d0/${name}_tet.vtk" ]; then
            python $E/tet_prep_v5.py --msh "$msh" --out-dir "$d0" \
                --out "${name}_tet.vtk" > "$d0/tet_prep.log" 2>&1 \
              || { echo "[batch] skip $name: tet_prep failed"; continue; }
        fi
        for n in $NS; do
            echo "$d0/${name}_tet.vtk $n $OUT_ROOT/${name}_n$n $cand/params.json" >> "$MAN"
        done
    done
fi

TOTAL=$(wc -l < "$MAN")
[ "$TOTAL" -gt 0 ] || { echo "[batch] nothing to do"; exit 1; }
echo "[batch] $TOTAL samples queued:"
sed 's/^/          /' "$MAN"

# ------------------------------------------------------------ state history
# One TSV row per finished sample attempt (ok/fail/prune, per-stage seconds),
# appended by sample_one.sh; the snapshot survives the walltime kill.
HIST="$OUT_ROOT/sample_history.tsv"
[ -f "$HIST" ] || printf '# ts\tmachine\tn\tstatus\tstage\trc\tstages\n' > "$HIST"
echo "[batch] history         $HIST"

if [ "$REPORT_AT" -gt 0 ] && [ "$(date +%s)" -lt "$REPORT_AT" ]; then
    (
        # Poll in 5 s ticks: kill only reaches this subshell, so one long sleep
        # would survive as an orphan and hold the tee pipe open -- job 6897408
        # idled 2 h 40 min at the final wait after its samples phase that way,
        # until it was cancelled. The whole subshell writes into the log file
        # so no surviving tick can hold the pipe.
        while [ "$(date +%s)" -lt "$REPORT_AT" ]; do sleep 5; done
        echo "[batch] ===== snapshot 2 min before the walltime ====="
        echo "[batch] samples ready  $(find "$OUT_ROOT" -name sample.npz 2>/dev/null | wc -l) of $TOTAL"
        echo "[batch] load           $(cut -d' ' -f1-3 /proc/loadavg)"
        echo "[batch] hm_threads     $(ps -eL -o comm= 2>/dev/null | grep -c '^HexMeshing')"
    ) >> "$LOG" 2>&1 &
    REPORTER=$!
else
    REPORTER=""
fi

# ------------------------------------------------------------- host monitor
# SLURM's per-job CPU accounting reports 0.00 % for these enroot payloads (job
# 6861458 finished six samples at 0.00 %, verified on disk), so it cannot tell
# a busy node from a stuck one -- the exact question job 6882450 left open.
# Sample the host once a minute into this log instead: load says whether the
# workers were doing work, hm_threads shows how many threads HexMeshing really
# uses (i.e. whether the OMP_NUM_THREADS=1 pin reaches the container). Appends
# straight to the file so its sleep ticks never hold the tee pipe either.
(
    while :; do
        printf '[mon] %s load=%s run=%s hm_threads=%s memMB=%s\n' \
            "$(date +%H:%M:%S)" \
            "$(cut -d' ' -f1-3 /proc/loadavg)" \
            "$(ps -eo stat= 2>/dev/null | grep -c '^R')" \
            "$(ps -eL -o comm= 2>/dev/null | grep -c '^HexMeshing')" \
            "$(free -m | awk '/Mem:/{print $3}')"
        sleep 60
    done
) >> "$LOG" 2>&1 &
MON=$!

# -------------------------------------------------------------------- run it
# Snapshot the samples that already exist: the end-of-run OK count includes
# them (a skip costs milliseconds) and would inflate the rate -- job 6897408
# started with 124 done and reported 141/hour for 114 new samples.
PRE=0
while read -r _tet _n d _rest; do
    if [ -f "$d/sample.npz" ]; then PRE=$(( PRE + 1 )); fi
done < "$MAN"
t0=$(date +%s)
# -P $P workers, one line each. --line means a whole manifest line becomes the
# argument list, so the worker gets <tet> <n> <outdir> [params] positionally.
xargs -a "$MAN" -P "$P" -L 1 bash "$E/scripts/sample_one.sh"
DT=$(( $(date +%s) - t0 ))
kill "$MON" 2>/dev/null || true
wait "$MON" 2>/dev/null || true
if [ -n "$REPORTER" ]; then
    kill "$REPORTER" 2>/dev/null || true
    wait "$REPORTER" 2>/dev/null || true
fi

# ------------------------------------------------------------------ summary
OK=0; MISS=0
while read -r _tet _n d _rest; do
    if [ -f "$d/sample.npz" ]; then OK=$(( OK + 1 )); else MISS=$(( MISS + 1 )); fi
done < "$MAN"

echo
echo "[batch] ================ throughput ================"
echo "[batch] wall            ${DT}s"
NEW=$(( OK - PRE ))
echo "[batch] samples         $OK ok ($NEW new this run), $MISS missing, of $TOTAL"
[ "$OK" -gt 0 ] && {
    echo "[batch] retained        $(du -sh "$OUT_ROOT" | cut -f1) total, $(( $(du -sk "$OUT_ROOT" | cut -f1) / OK )) KB per sample"
    echo "[batch] measured rate   $(( NEW * 3600 / (DT > 0 ? DT : 1) )) samples/hour this run"
    # The rate is only a capacity number if the workers were actually busy.
    # With fewer NEW samples than P the wall time is the duration of the
    # SLOWEST sample and the rate says nothing about the node -- the first run
    # of this script reported 8/hour for 3 samples on 32 slots, where the
    # capacity was about 89/hour. Say which of the two this run measured.
    if [ "$NEW" -lt "$P" ]; then
        echo "[batch] ^^ NOT capacity: only $NEW new samples for $P slots, so the"
        echo "[batch]    wall time is just the slowest single sample."
        echo "[batch]    capacity would be ~$(( P * 3600 / (DT > 0 ? DT : 1) ))/hour -- UNVERIFIED,"
        echo "[batch]    queue at least $P new samples to measure it."
    else
        echo "[batch] ^^ saturated ($NEW new >= $P), so this IS a capacity number."
    fi
}
[ "$MISS" -gt 0 ] && {
    echo "[batch] failures, newest log first:"
    while read -r _tet _n d _rest; do
        [ -f "$d/sample.npz" ] || echo "          $d  -> $(ls -t "$d"/*.log 2>/dev/null | head -1)"
    done < "$MAN"
}
echo "[batch] ============================================"
echo "[batch] Resubmit to retry only the missing ones -- sample_one.sh skips"
echo "[batch] anything that already has a sample.npz."
# Let the tee drain before the script exits and SLURM sees the job as done.
# A bare wait deadlocks here: the shell holds the process-substitution pipe
# open, the tee never sees EOF -- job 6890146 idled 3.5 h this way until the
# walltime killed it. Closing stdout/stderr first closes the pipe, the tee
# exits, then wait returns.
exec >/dev/null 2>&1
wait
[ "$MISS" -eq 0 ]
