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
PYTHON_MODULE="${PYTHON_MODULE:-devel/python/3.13.3-gnu-14.2}"
NS="${NS:-500 2000 8000}"          # the -n values; each becomes its own sample
CANDS="${CANDS:-}"                 # candidate dirs; empty = the fixtures tet
OUT_ROOT="${OUT_ROOT:-output/hex3d_algohex/batch}"

module load "$PYTHON_MODULE" 2>/dev/null || true
export PATH="$REPO/.venv/bin:$PATH"
export REPO
export ENROOT_DATA_PATH="${ENROOT_DATA_PATH:-$WS/enroot-data}"
export ENROOT_TEMP_PATH="${ENROOT_TEMP_PATH:-/tmp/${USER:-$(id -un)}-enroot}"
export ENROOT_RUNTIME_PATH="${ENROOT_RUNTIME_PATH:-/tmp/${USER:-$(id -un)}-enroot-run}"
mkdir -p "$ENROOT_TEMP_PATH" "$ENROOT_RUNTIME_PATH"

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

# ------------------------------------------------------------------ manifest
# One line per sample: <tet> <n> <outdir> [params]
MAN="$OUT_ROOT/manifest.txt"
: > "$MAN"

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

# -------------------------------------------------------------------- run it
t0=$(date +%s)
# -P $P workers, one line each. --line means a whole manifest line becomes the
# argument list, so the worker gets <tet> <n> <outdir> [params] positionally.
xargs -a "$MAN" -P "$P" -L 1 bash "$E/scripts/sample_one.sh"
DT=$(( $(date +%s) - t0 ))

# ------------------------------------------------------------------ summary
OK=0; MISS=0
while read -r _tet _n d _rest; do
    if [ -f "$d/sample.npz" ]; then OK=$(( OK + 1 )); else MISS=$(( MISS + 1 )); fi
done < "$MAN"

echo
echo "[batch] ================ throughput ================"
echo "[batch] wall            ${DT}s"
echo "[batch] samples         $OK ok, $MISS missing, of $TOTAL"
[ "$OK" -gt 0 ] && {
    echo "[batch] retained        $(du -sh "$OUT_ROOT" | cut -f1) total, $(( $(du -sk "$OUT_ROOT" | cut -f1) / OK )) KB per sample"
    echo "[batch] measured rate   $(( OK * 3600 / (DT > 0 ? DT : 1) )) samples/hour"
    # The measured rate is only a capacity number if the workers were actually
    # busy. With TOTAL < P the wall time is the duration of the SLOWEST single
    # sample and the rate says nothing about the node -- the first run of this
    # script reported 8/hour for 3 samples on 32 slots, where the capacity was
    # about 89/hour. Say which of the two this run measured.
    if [ "$TOTAL" -lt "$P" ]; then
        echo "[batch] ^^ NOT capacity: $TOTAL samples for $P slots, so the wall"
        echo "[batch]    time is just the slowest single sample."
        echo "[batch]    capacity would be ~$(( P * 3600 / (DT > 0 ? DT : 1) ))/hour -- UNVERIFIED,"
        echo "[batch]    queue at least $P samples to measure it."
    else
        echo "[batch] ^^ saturated ($TOTAL >= $P), so this IS a capacity number."
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
[ "$MISS" -eq 0 ]
