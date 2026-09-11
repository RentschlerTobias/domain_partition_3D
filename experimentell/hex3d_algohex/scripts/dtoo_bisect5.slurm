#!/bin/bash
#SBATCH --job-name=dtoo_bisect5
#SBATCH --partition=dev_cpu_il
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=dtoo_bisect5_%j.out
#SBATCH --exclude=uc2n601

# bisect5: is the silent container-start stall node-side state or the shared
# container store?
#
# History: bisect3 v2 started the shared-store container on uc2n601 in 2-3 s
# (job 6870239, 20:23). After the export run 6870688 was SIGKILLed at its
# 30-min time limit, every later start on uc2n601 stalled with zero output
# and zero CPU (bisect4 job 6871118: step RUNNING 15:52, TotalCPU 00:00:00,
# no marker file written). enroot takes an flock on the container rootfs
# during setup (src/runtime.sh); a stuck Lustre flock would match a silent
# zero-CPU hang. The shared store itself could also have been left poisoned
# by the SIGKILL.
#
# This probe therefore runs with --exclude=uc2n601, i.e. on a different node
# on purpose (remove the --exclude line to probe the suspect node instead).
#   S0  host diagnostics (runtime dirs, enroot/fuse processes, mounts)
#   S1  flock probe on the container rootfs lock
#   S2  dtOO sanity start (exec echo, 90 s bound)
#   S3  dolfinx control start (enroot-wide vs dtOO-specific, 90 s bound)
#   S4  one bounded real export (machine_0001) with a marker trace
#
# Read the summary line first: S2/S3/S4 rc values.
#   S2 OK + S4 export OK       -> uc2n601 is the problem; the store is fine.
#   S2 OK + S4 stalls          -> clean node stalls too: look at the store.
#   S2 stalls on this node too -> shared container store poisoned:
#                                 enroot remove -f dtOO; enroot create --name dtOO $SQSH
#   S3 stalls while S2 OK      -> enroot-wide / node bootstrap issue.
#   S1 rc=1 or 124             -> rootfs lock busy or stuck.

set -u

module load devel/python/3.13.3-gnu-14.2

WS=$(ws_find hex3d)
REPO="${REPO:-${SLURM_SUBMIT_DIR:-$WS/domain_partition_3D}}"
EIG_WS=$(ws_find eigenfreq)

export ENROOT_TEMP_PATH="${ENROOT_TEMP_PATH:-/tmp/$USER-enroot}"
export ENROOT_RUNTIME_PATH="${ENROOT_RUNTIME_PATH:-/tmp/$USER-enroot-run}"
export ENROOT_DATA_PATH="${ENROOT_DATA_PATH:-$EIG_WS/enroot-data}"
SQSH="$EIG_WS/enroot-images/dtOO.sqsh"
CONTAINER=dtOO

EF=""
for cand in "$EIG_WS/eigenfrequencies" "$HOME/eigen/eigenfrequencies" "$HOME/eigenfrequencies"; do
    if [ -f "$cand/src/eigenfrequencies/adapters/dtoo/adapter.py" ]; then EF="$cand"; break; fi
done
[ -n "$EF" ] || { echo "[bisect5] no eigenfrequencies checkout"; exit 1; }
RC="$EF/cluster/enroot_rc.sh"
[ -f "$RC" ] || { echo "[bisect5] no --rc command script: $RC"; exit 1; }
[ -f "$SQSH" ] || { echo "[bisect5] no dtOO image: $SQSH"; exit 1; }

LOGDIR="$REPO/experimentell/hex3d_algohex/logs"
mkdir -p "$LOGDIR"
LOG="$LOGDIR/dtoo_bisect5_${SLURM_JOB_ID:-manual}.log"
exec > >(tee "$LOG") 2>&1

echo "[bisect5] host $(hostname), job ${SLURM_JOB_ID:-manual}, $(date -Is)"
echo "[bisect5] repo $REPO"
echo "[bisect5] ef   $EF"
echo "[bisect5] img  $SQSH"
echo "[bisect5] data $ENROOT_DATA_PATH"
enroot list || true
enroot list | grep -qxF "$CONTAINER" || echo "[bisect5] warning: container $CONTAINER not in enroot list"

run_step() {
    local label=$1 t=$2
    shift 2
    echo "[bisect5] --- $label --- $(date -Is)"
    timeout "$t" srun -n1 -N1 enroot start --root "$@" < /dev/null
    local rc=$?
    echo "[bisect5] --- $label rc=$rc --- $(date -Is)"
    return "$rc"
}

echo "[bisect5] --- S0 host diagnostics ---"
ls -la "/tmp/$USER-enroot" 2>&1 | head -20
ls -la "/tmp/$USER-enroot-run" 2>&1 | head -20
ls -la "$ENROOT_DATA_PATH/$CONTAINER" 2>&1 | head -10
ps -ef | grep -E 'enroot|fuse|squashfuse|unsquashfs' | grep -v grep || echo "[bisect5] no enroot/fuse processes"
findmnt | grep -iE 'enroot|fuse|/repo|/ef' || echo "[bisect5] no relevant mounts"
ls -l /sys/fs/fuse/connections 2>&1 | head -5 || true

echo "[bisect5] --- S1 rootfs lock probe ---"
LOCK="$ENROOT_DATA_PATH/$CONTAINER/.lock"
if [ -e "$LOCK" ]; then
    ls -la "$LOCK"
    timeout 5 flock -n "$LOCK" -c 'echo LOCK-FREE'
    echo "[bisect5] S1 flock rc=$? (0 free, 1 busy, 124 stuck)"
else
    echo "[bisect5] no lock file $LOCK (normal between runs)"
fi

R2=n/a R3=n/a R4=n/a
MS=(--mount "$REPO:/repo" --mount "$EF:/ef")

run_step "S2 dtOO sanity start" 90 --rc "$RC" "$CONTAINER" bash -c 'exec echo SANITY-OK'
R2=$?

if enroot list | grep -qxF dolfinx; then
    run_step "S3 dolfinx control start" 90 --rc "$RC" dolfinx bash -c 'exec echo DOLFINX-OK'
    R3=$?
else
    echo "[bisect5] dolfinx container not present, S3 skipped"
fi

rm -f "$LOGDIR/probe_markers5.txt"
run_step "S4 one-machine export with markers" 600 --rc "$RC" "${MS[@]}" "$CONTAINER" bash -c 'source /usr/lib/openfoam/openfoam2606/etc/bashrc; source /dtOO-install/bin/env.sh; export EIGENFREQUENCIES_ROOT=/ef PYTHONUNBUFFERED=1; exec bash /repo/experimentell/hex3d_algohex/scripts/probe_markers.sh machine_0001'
R4=$?

echo "[bisect5] --- harvest ---"
echo "[bisect5] summary: S2_dtoo=$R2 S3_dolfinx=$R3 S4_export=$R4"
ls -la "$REPO/data/dataset/sobol/machine_0001/"
if [ -f "$LOGDIR/probe_markers5.txt" ]; then
    echo "[bisect5] --- marker trace ---"
    cat "$LOGDIR/probe_markers5.txt"
else
    echo "[bisect5] no marker file -- the container never wrote into /repo"
fi
echo "[bisect5] log: $LOG"
echo "[bisect5] commit this log: git add experimentell/hex3d_algohex/logs && git commit -m 'cluster: dtOO bisect5 log'"
wait
