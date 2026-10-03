#!/bin/bash
#SBATCH --job-name=dtoo_bisect4
#SBATCH --partition=dev_cpu_il
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=experimentell/hex3d_algohex/logs/dtoo_bisect4_%j.out
#
# dtoo bisect4 -- why does the batch export burn ~17 CPU-min and 1 GiB RAM
# yet write NOTHING (job 6870688: machine dirs stayed with params.json only,
# no mesh.msh.part, no dtoo_build.log, zero payload output)?
#
# The 6870688 evidence says the container step was ALIVE for the full 30 min
# (CPU Utilized 00:16:55, Memory 1.09 GB) while every filesystem trace the
# export would create (mesh.msh.part from writeMSH, dtoo_build.log from
# logMe.initLog, design_*.json temp from run_dtoo_export) is absent. Two
# candidates: (a) the lazy `from dtOOPythonSWIG import ...` (never exercised
# in-container; bisect3 P4d only imported the adapter, which is SWIG-free)
# or (b) the dtOO compute path itself (loadStateToConst/destroyAndCreate/
# makeGrid) with output lost in srun buffering on SIGKILL.
#
# Everything is appended to $REPO/probe_markers.txt INSIDE the container --
# a file in the repo is immune to the srun output buffering that swallowed
# every stdout line of job 6870688. Each stage is timeout-bounded so the
# step ends NORMALLY and flushes whatever came before.
set -u

module load devel/python/3.13.3-gnu-14.2

WS=$(ws_find hex3d)
EIG_WS=$(ws_find eigenfreq)
REPO="${REPO:-${SLURM_SUBMIT_DIR:-$WS/domain_partition_3D}}"
DATASET="${DATASET:-data/dataset/sobol}"
CONTAINER=dtOO

export ENROOT_TEMP_PATH="${ENROOT_TEMP_PATH:-/tmp/$USER-enroot}"
export ENROOT_RUNTIME_PATH="${ENROOT_RUNTIME_PATH:-/tmp/$USER-enroot-run}"
export ENROOT_DATA_PATH="${ENROOT_DATA_PATH:-$EIG_WS/enroot-data}"
export PYTHONUNBUFFERED=1

SQSH="$EIG_WS/enroot-images/dtOO.sqsh"
[ -f "$SQSH" ] || { echo "[bisect4] no dtOO image: $SQSH"; exit 1; }
[ -d "$REPO/$DATASET" ] || { echo "[bisect4] no dataset root: $REPO/$DATASET"; exit 1; }

EF=""
for cand in "$EIG_WS/eigenfrequencies" "$HOME/eigen/eigenfrequencies" "$HOME/eigenfrequencies"; do
    if [ -f "$cand/src/eigenfrequencies/adapters/dtoo/adapter.py" ]; then EF="$cand"; break; fi
done
[ -n "$EF" ] || { echo "[bisect4] no eigenfrequencies checkout"; exit 1; }
RC="$EF/cluster/enroot_rc.sh"
[ -f "$RC" ] || { echo "[bisect4] no --rc script: $RC"; exit 1; }

LOGDIR="$REPO/experimentell/hex3d_algohex/logs"
mkdir -p "$LOGDIR"
exec > >(tee "$LOGDIR/dtoo_bisect4_${SLURM_JOB_ID:-manual}.log") 2>&1

echo "[bisect4] host $(hostname -f), job ${SLURM_JOB_ID:-manual}, $(date -Is)"
echo "[bisect4] repo $REPO  ef $EF"
echo "[bisect4] img  $SQSH  data $ENROOT_DATA_PATH"

# stale markers out of the way (host side; the container sees the same dir)
rm -f "$REPO/probe_markers.txt"

enroot list 2>/dev/null | grep -qxF "$CONTAINER" \
    || { echo "[bisect4] container $CONTAINER not unpacked"; exit 1; }

echo "[bisect4] --- single bounded step ---"
timeout 1500 srun -n1 -N1 enroot start --root \
    --rc "$RC" \
    --mount "$REPO:/repo" \
    --mount "$EF:/ef" \
    "$CONTAINER" bash -c '
    source /usr/lib/openfoam/openfoam2606/etc/bashrc
    source /dtOO-install/bin/env.sh
    export EIGENFREQUENCIES_ROOT=/ef
    export PYTHONUNBUFFERED=1
    cd /repo || exit 3
    M=/repo/probe_markers.txt
    { date -Is; echo "stage=start pwd=$(pwd)"; command -v python3.13; python3.13 -V; } >> "$M" 2>&1
    echo "stage=swig-import $(date +%T)" >> "$M"
    timeout 120 python3.13 -u -X importtime -c "from dtOOPythonSWIG import baseContainer; print(\"SWIG-OK\", flush=True)" >> "$M" 2>&1
    echo "stage=swig-rc=$? $(date +%T)" >> "$M"
    echo "stage=generate-machines $(date +%T)" >> "$M"
    timeout 600 python3.13 -u experimentell/hex3d_algohex/scripts/generate_machines.py --export --only machine_0001 >> "$M" 2>&1
    echo "stage=gen-rc=$? $(date +%T)" >> "$M"
    { echo "stage=machine-0001-ls"; ls -la data/dataset/sobol/machine_0001/; echo "parts=$(ls data/dataset/sobol/machine_*/mesh.msh.part 2>/dev/null | wc -l)"; echo "stage=residues"; ls -la /tmp; find /tmp /root -maxdepth 4 -name "dtoo_build*" 2>/dev/null; find /tmp -maxdepth 2 -name "design_*" 2>/dev/null; } >> "$M" 2>&1
    echo "stage=done $(date +%T)" >> "$M"
    exec true
' < /dev/null
R=$?
echo "[bisect4] step rc=$R"

# ------------------------------------------------------------------ harvest
echo "[bisect4] ================ probe_markers.txt ================"
if [ -f "$REPO/probe_markers.txt" ]; then
    cat "$REPO/probe_markers.txt"
else
    echo "(no marker file -- the container never wrote inside /repo)"
fi
PF="$REPO/$DATASET"
OK=0; FAIL=0; PART=0
for d in "$PF"/machine_00*; do
    [ -d "$d" ] || continue
    if [ -f "$d/mesh.msh" ]; then OK=$(( OK + 1 ))
    elif [ -f "$d/mesh.msh.part" ]; then PART=$(( PART + 1 ))
    elif [ -f "$d/export_error.txt" ]; then FAIL=$(( FAIL + 1 )); fi
done
echo "[bisect4] machines: $OK mesh.msh, $PART mesh.msh.part, $FAIL export_error.txt"
echo "[bisect4] log: $LOGDIR/dtoo_bisect4_${SLURM_JOB_ID:-manual}.log"
echo "[bisect4] keep this log: git add experimentell/hex3d_algohex/logs && git commit -m 'cluster: dtoo bisect4 log'"
wait
