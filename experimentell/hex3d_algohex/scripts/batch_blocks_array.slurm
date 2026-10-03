#!/bin/bash
# batch_blocks_array.slurm: N sharded sample jobs, one node each.
#
# Same unit of work, environment and idempotency conventions as
# batch_samples.slurm (sample_one.sh: algohex -> clean_blocks -> tfi ->
# export_sample). The only difference is the manifest: the candidate
# directories are split across shards by index modulo, so N array jobs
# partition the machine set without overlapping or locking.
#
#   sbatch -a 0-2 .../batch_blocks_array.slurm                # 3 nodes
#   CANDS="output/hex3d_algohex/batch/*" NS="8000" \
#       sbatch -a 0-2 -t 10:00:00 .../batch_blocks_array.slurm
#
# WHY SHARDS INSTEAD OF ONE BIG MANIFEST. batch_samples.sh builds one shared
# manifest and truncates it (": > MAN"), and two concurrent jobs truncating
# the same manifest race. Sharding the candidates by index gives each job a
# disjoint manifest by construction, so sample_one.sh idempotency is only
# needed against a retry, not against another node.
#
# WHY A PER-SHARD ENROOT DATA PATH. The algohex enroot container unpacks
# with --root --rw on a shared unpacked rootfs (see batch_samples.slurm
# header); concurrent nodes need one container copy each. ENROOT_DATA_PATH
# therefore points to <base>/shard<SHARD> by default; existing unpacked
# containers are reused per shard automatically after the first run.
#SBATCH --job-name=hex3d_blocks
#SBATCH --output=experimentell/hex3d_algohex/logs/hex3d_blocks_%a_%j.out
#SBATCH --partition=cpu_il
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G

set -u

SHARD="${SHARD:-${SLURM_ARRAY_TASK_ID:-0}}"
NSHARDS="${NSHARDS:-1}"

WS="$(ws_find hex3d 2>/dev/null)"
REPO="${REPO:-${SLURM_SUBMIT_DIR:-$WS/domain_partition_3D}}"
[ -n "$WS" ] && [ -d "$REPO/experimentell/hex3d_algohex" ] || {
    echo "[blocks] not the hex3d checkout: WS=$WS REPO=$REPO"; exit 1; }

LOGDIR="$REPO/experimentell/hex3d_algohex/logs"
mkdir -p "$LOGDIR"
LOG="$LOGDIR/batch_blocks_s${SHARD}_of_${NSHARDS}_${SLURM_JOB_ID:-manual}.log"
exec > >(tee -a "$LOG") 2>&1

PYTHON_MODULE="${PYTHON_MODULE:-devel/python/3.13.3-gnu-14.2}"
CANDS="${CANDS:?CANDS=<machine-dirs|exported-tets> required}"
NS="${NS:-8000}"
OUT_ROOT="${OUT_ROOT:-output/hex3d_algohex/batch_blocks}"
MESH_BASE="${MESH_BASE:-}"          # where mesh.msh lives, only needed for CANDS without tets

module load "$PYTHON_MODULE" 2>/dev/null || true
export PATH="$REPO/.venv/bin:$PATH"
export REPO
export ENROOT_DATA_PATH="${ENROOT_DATA_PATH:-$WS/enroot-data-bs${SHARD}}"
export ENROOT_TEMP_PATH="${ENROOT_TEMP_PATH:-/tmp/${USER:-$(id -un)}-enroot-bs${SHARD}}"
export ENROOT_RUNTIME_PATH="${ENROOT_RUNTIME_PATH:-/tmp/${USER:-$(id -un)}-enroot-run-bs${SHARD}}"
mkdir -p "$ENROOT_TEMP_PATH" "$ENROOT_RUNTIME_PATH"

cd "$REPO" || exit 1
E=experimentell/hex3d_algohex
mkdir -p "$OUT_ROOT"

# --------------------------------------------------- shard the candidate set
# Discovery order is fixed by ls -1, so every shard enumerates the same list
# and takes its own residue classes exactly once.
MAN="$OUT_ROOT/manifest_s${SHARD}_of_${NSHARDS}.txt"
: > "$MAN"

i=0
for cand in $CANDS; do
    if [ $(( i % NSHARDS )) -ne "$SHARD" ]; then i=$(( i + 1 )); continue; fi
    i=$(( i + 1 ))
    name=$(basename "$cand")
    [ -d "$cand" ] && [ -f "$cand/sample.npz" ] && { echo "skip $name (sample.npz exists)"; continue; }
    if [ -f "$cand/${name}_tet.vtk" ] && [ ! -f "$cand/mesh.msh" ]; then
        # Pre-exported tet: skip the tet_prep step entirely.
        for n in $NS; do
            echo "$cand/${name}_tet.vtk $n $OUT_ROOT/${name}_n$n ${cand}/params.json" >> "$MAN"
        done
        continue
    fi
    msh="${MESH_BASE:+$MESH_BASE/}$cand/mesh.msh"
    [ -f "$msh" ] || msh="$cand/mesh.msh"
    [ -f "$msh" ] || { echo "skip $name: no mesh.msh and no tet"; continue; }
    d0="$OUT_ROOT/$name"
    mkdir -p "$d0"
    if [ ! -f "$d0/${name}_tet.vtk" ]; then
        python "$E/tet_prep_v5.py" --msh "$msh" --out-dir "$d0" \
            --out "${name}_tet.vtk" > "$d0/tet_prep.log" 2>&1 \
          || { echo "skip $name: tet_prep failed"; continue; }
    fi
    for n in $NS; do
        echo "$d0/${name}_tet.vtk $n $OUT_ROOT/${name}_n$n ${cand}/params.json" >> "$MAN"
    done
done

TOTAL=$(wc -l < "$MAN")
[ "$TOTAL" -gt 0 ] || { echo "[blocks] shard $SHARD/$NSHARDS: nothing to do"; exit 0; }
echo "[blocks] shard $SHARD/$NSHARDS: $TOTAL samples queued"
sed 's/^/          /' "$MAN"

# shared history, single O_APPEND rows are safe across nodes
HIST="$OUT_ROOT/sample_history.tsv"
[ -f "$HIST" ] || printf '# ts\tmachine\tn\tstatus\tstage\trc\tstages\n' > "$HIST"

# ------------------------------------------------ concurrency (same as batch_samples)
CORES="${SLURM_CPUS_PER_TASK:-$(nproc)}"
MEM_MB="${SLURM_MEM_PER_NODE:-131072}"
MEM_PER_SAMPLE_MB="${MEM_PER_SAMPLE_MB:-2500}"
CORES_PER_SAMPLE="${CORES_PER_SAMPLE:-2}"
P_MEM=$(( MEM_MB / MEM_PER_SAMPLE_MB ))
P_CPU=$(( CORES / CORES_PER_SAMPLE ))
P="${PARALLEL:-$(( P_MEM < P_CPU ? P_MEM : P_CPU ))}"
[ "$P" -ge 1 ] || P=1
echo "[blocks] host $(hostname), task ${SLURM_ARRAY_TASK_ID:-none}, concurrency $P"

# squeue prints minutes:seconds below one hour (the dev partitions).
parse_wall() {
    case "$1" in
        [0-9]*-*:*:*)
            d=${1%%-*}; t=${1#*-}; h=${t%%:*}; m=${t#*:}; echo $(( 10#${d} * 86400 + 10#$h * 3600 )) ;;
        [0-9]*-*:*)
            h=${1%%-*}; m=${1#*-}; echo $(( 10#$h * 3600 )) ;;
        [0-9]*:*:*)
            h=${1%%:*}; echo $(( 10#$h * 3600 )) ;;
        [0-9]*:*)
            m=${1%%:*}; echo $(( 10#$m * 60 )) ;;
        *) echo "" ;;
    esac
}
DRAIN_SECONDS="${DRAIN_SECONDS:-1500}"
LEFT=$(parse_wall "$(squeue -h -j "${SLURM_JOB_ID:-0}" -o %L 2>/dev/null || true)" 2>/dev/null || true)
JOB_END=0; DRAIN_DEADLINE=0
if [ -n "${LEFT:-}" ] && [ "$LEFT" -gt 0 ]; then
    JOB_END=$(( $(date +%s) + LEFT ))
    DRAIN_DEADLINE=$(( JOB_END - DRAIN_SECONDS ))
    echo "[blocks] walltime ${LEFT}s, no new samples after the last ${DRAIN_SECONDS}s"
fi
export SAMPLE_TIMEOUT="${SAMPLE_TIMEOUT:-2400}" JOB_END DRAIN_DEADLINE

PRE=0
while read -r _tet _n d _rest; do
    [ -f "$d/sample.npz" ] && PRE=$(( PRE + 1 ))
done < "$MAN"

MON=/dev/null
t0=$(date +%s)
xargs -a "$MAN" -P "$P" -L 1 bash "$E/scripts/sample_one.sh"
DT=$(( $(date +%s) - t0 ))

OK=0; MISS=0
while read -r _tet _n d _rest; do
    if [ -f "$d/sample.npz" ]; then OK=$(( OK + 1 )); else MISS=$(( MISS + 1 )); fi
done < "$MAN"
NEW=$(( OK - PRE ))
echo
echo "[blocks] =========== shard $SHARD/$NSHARDS throughput ==========="
echo "[blocks] wall            ${DT}s"
echo "[blocks] samples         $OK ok ($NEW new this run), $MISS missing, of $TOTAL"
[ "$OK" -gt 0 ] && echo "[blocks] measured rate   $(( NEW * 3600 / (DT > 0 ? DT : 1) )) samples/hour this run"
[ "$MISS" -gt 0 ] && {
    echo "[blocks] failures, newest log first:"
    while read -r _tet _n d _rest; do
        [ -f "$d/sample.npz" ] || echo "           $d  -> $(ls -t "$d"/*.log 2>/dev/null | head -1)"
    done < "$MAN"
}
echo "[blocks] ============================================"
echo "[blocks] Resubmit to retry only the missing ones."
exec >/dev/null 2>&1
wait
[ "$MISS" -eq 0 ]
