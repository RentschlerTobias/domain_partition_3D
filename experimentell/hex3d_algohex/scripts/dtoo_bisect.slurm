#!/bin/bash
# dtOO container start bisect. One sbatch job, stepwise, timeout-wrapped.
#
#   sbatch experimentell/hex3d_algohex/scripts/dtoo_bisect.slurm
#
# WHY. Job 6868031 passed every guard, printed 'container dtOO already
# unpacked, skipping', and then sat silent for the full 30-minute slot with
# ZERO container output -- no '[sampler] exporting', nothing. The same was
# true of the interactive probe. Something between 'srun enroot start' and
# python's first print blocks, and guessing further from the outside is
# wasteful. This job isolates it in one run.
#
# Steps, each wrapped in timeout so a hang cannot eat the slot:
#   1  enroot start --root --rc <rc> dtOO bash -c 'echo'   (no mounts)
#   1b same but WITHOUT --rc            -- only if 1 fails
#   0b fresh enroot create of a probe container, then start it
#                                       -- only if 1 fails (is dtOO broken?)
#   2  + the two mounts, both env files sourced, echo SOURCED-OK
#   3  + python3.13 -V and one trivial import
#   4  + adapter import and design_bounds()      (the import seam)
#   5  + a real bounded export window (10 min), then the harvest counts
#
# Read the log top to bottom: the first step that times out (rc=124) or fails
# (rc!=0) is the answer. Every step names itself and stamps the time.
#
# Env facts, same as batch_dtoo_export.slurm: container lives in the
# eigenfrequencies workspace store, --rc is a HOST path, mounts mount the
# repo at /repo and the eigenfrequencies checkout at /ef.
# ---------------------------------------------------------------------------
#SBATCH --job-name=dtoo_bisect
#SBATCH --output=experimentell/hex3d_algohex/logs/dtoo_bisect_%j.out
#SBATCH --partition=dev_cpu_il
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

set -u

WS="${WS:-$(ws_find hex3d 2>/dev/null)}"
REPO="${REPO:-${SLURM_SUBMIT_DIR:-$WS/domain_partition_3D}}"
EIG_WS="${EIG_WS:-$(ws_find eigenfreq 2>/dev/null)}"
SQSH="${SQSH:-$EIG_WS/enroot-images/dtOO.sqsh}"
DATASET="${DATASET:-data/dataset/sobol}"
COUNT="${COUNT:-64}"

module load "${PYTHON_MODULE:-devel/python/3.13.3-gnu-14.2}" 2>/dev/null || true
export REPO
export ENROOT_DATA_PATH="${ENROOT_DATA_PATH:-$EIG_WS/enroot-data}"
export ENROOT_TEMP_PATH="${ENROOT_TEMP_PATH:-/tmp/${USER:-$(id -un)}-enroot}"
export ENROOT_RUNTIME_PATH="${ENROOT_RUNTIME_PATH:-/tmp/${USER:-$(id -un)}-enroot-run}"
mkdir -p "$ENROOT_TEMP_PATH" "$ENROOT_RUNTIME_PATH" "$ENROOT_DATA_PATH"

cd "$REPO" || exit 1
LOGDIR="$REPO/experimentell/hex3d_algohex/logs"
mkdir -p "$LOGDIR"
LOG="$LOGDIR/dtoo_bisect_${SLURM_JOB_ID:-manual}.log"
exec > >(tee "$LOG") 2>&1

EF_ROOT=""
for c in "$EIG_WS/eigenfrequencies" "$HOME/eigen/eigenfrequencies" "$HOME/eigenfrequencies"; do
    [ -f "$c/src/eigenfrequencies/adapters/dtoo/adapter.py" ] && { EF_ROOT="$c"; break; }
done
RC_SCRIPT="$EF_ROOT/cluster/enroot_rc.sh"
CONTAINER="${CONTAINER:-dtOO}"

echo "[bisect] host $(hostname), job ${SLURM_JOB_ID:-none}, $(date -Is)"
echo "[bisect] repo $REPO"
echo "[bisect] ef   $EF_ROOT"
echo "[bisect] img  $SQSH"
echo "[bisect] data $ENROOT_DATA_PATH"
echo "[bisect] $(enroot version 2>&1 | head -1)"
du -sh "$ENROOT_DATA_PATH/$CONTAINER" 2>/dev/null || echo "[bisect] no container dir $ENROOT_DATA_PATH/$CONTAINER"

# --- shared helper: run one container command under a timeout ----------------
# rc is echoed by the caller: 0 ok, 124 timed out, other = failure.
run_step() {  # run_step <timeout_s> <label> <container> [--] <args...>
    local t="$1" label="$2"; shift 2
    echo "[bisect] --- $label --- $(date +%T)"
    timeout "$t" srun -n1 -N1 enroot start --root "$@" < /dev/null
    echo "[bisect] --- $label rc=$? --- $(date +%T)"
}

# --- step 1: bare start, --rc, no mounts ------------------------------------
run_step 180 "step1 start+rc" --rc "$RC_SCRIPT" "$CONTAINER" \
    bash -c 'exec echo STEP1-OK'
rc1=$?

# --- step 1b: same without --rc (only if step 1 failed) ---------------------
if [ "$rc1" -ne 0 ]; then
    run_step 180 "step1b start no-rc" "$CONTAINER" \
        bash -c 'exec echo STEP1B-OK'
fi

# --- step 0b: is dtOO itself broken? create a fresh probe container ---------
# Unpacking 2.11 GB from Lustre takes minutes; 600 s covers it. Only on
# failure of step 1, because a healthy dtOO makes this pointless.
if [ "$rc1" -ne 0 ]; then
    probe="dtOO_bisect.$SLURM_JOB_ID"
    echo "[bisect] --- step0b fresh create $probe --- $(date +%T)"
    timeout 600 enroot create --name "$probe" "$SQSH"
    echo "[bisect] --- step0b create rc=$? --- $(date +%T)"
    if enroot list 2>/dev/null | grep -qxF "$probe"; then
        run_step 180 "step0b start fresh" "$probe" bash -c 'exec echo FRESH-OK'
        enroot remove -f "$probe" 2>/dev/null || true
    fi
fi

# --- step 2: the real invocation minus python -------------------------------
# Both env files BEFORE python3.13, the container trap this repo documented.
run_step 240 "step2 mounts+source" \
    --rc "$RC_SCRIPT" --mount "$REPO:/repo" --mount "$EF_ROOT:/ef" \
    "$CONTAINER" bash -c '
        source /usr/lib/openfoam/openfoam2606/etc/bashrc
        source /dtOO-install/bin/env.sh
        exec echo SOURCED-OK'

# --- step 3: python itself ---------------------------------------------------
run_step 240 "step3 python" \
    --rc "$RC_SCRIPT" --mount "$REPO:/repo" --mount "$EF_ROOT:/ef" \
    "$CONTAINER" bash -c '
        source /usr/lib/openfoam/openfoam2606/etc/bashrc
        source /dtOO-install/bin/env.sh
        command -v python3.13
        python3.13 -V
        PYTHONUNBUFFERED=1 python3.13 -c "import sys; print(sys.executable); import yaml; print(\"yaml ok\")"
        exec echo PY-OK'

# --- step 4: the import seam -------------------------------------------------
# generate_machines.py's first print comes AFTER this import, which is why a
# hang here looks exactly like a dead job.
run_step 300 "step4 adapter import" \
    --rc "$RC_SCRIPT" --mount "$REPO:/repo" --mount "$EF_ROOT:/ef" \
    "$CONTAINER" bash -c '
        source /usr/lib/openfoam/openfoam2606/etc/bashrc
        source /dtOO-install/bin/env.sh
        export EIGENFREQUENCIES_ROOT=/ef
        PYTHONUNBUFFERED=1 python3.13 -u -c "
import sys, time
t=time.time(); sys.path.insert(0, \"/ef/src\")
from eigenfrequencies.adapters.dtoo.adapter import DtooAdapter
print(\"adapter import\", round(time.time()-t,1), \"s\", flush=True)
t=time.time(); b=DtooAdapter(\"/ef/adapters/machines/tistos.yaml\").design_bounds()
print(\"design_bounds\", len(b), \"params\", round(time.time()-t,1), \"s\", flush=True)
"
        exec echo ADAPTER-OK'

# --- step 5: bounded export window ------------------------------------------
# 600 s of real exporting. Even two finished machines prove the chain; the
# full 64 is a normal batch run after this.
echo "[bisect] --- step5 bounded export 600s --- $(date +%T)"
timeout 600 srun -n1 -N1 enroot start --root \
    --rc "$RC_SCRIPT" --mount "$REPO:/repo" --mount "$EF_ROOT:/ef" \
    "$CONTAINER" bash -c '
        source /usr/lib/openfoam/openfoam2606/etc/bashrc
        source /dtOO-install/bin/env.sh
        export EIGENFREQUENCIES_ROOT=/ef
        export PYTHONUNBUFFERED=1
        cd /repo
        python3.13 experimentell/hex3d_algohex/scripts/generate_machines.py \
            --export --count '"$COUNT"'' < /dev/null
echo "[bisect] --- step5 rc=$? --- $(date +%T)"

# --- harvest -----------------------------------------------------------------
PF="$REPO/$DATASET"
OK=0; MISS=0; FAIL=0
for d in "$PF"/machine_00*; do
    [ -d "$d" ] || continue
    if [ -f "$d/mesh.msh" ]; then OK=$(( OK + 1 ))
    elif [ -f "$d/export_error.txt" ]; then FAIL=$(( FAIL + 1 ))
    else MISS=$(( MISS + 1 )); fi
done
echo
echo "[bisect] ================ bisect summary ================"
echo "[bisect] machines        $OK exported, $FAIL failed, $MISS not attempted"
[ "$FAIL" -gt 0 ] && ls "$PF"/*/export_error.txt 2>/dev/null | head -3
echo "[bisect] log file:       $LOG"
echo "[bisect] keep this log:  git add experimentell/hex3d_algohex/logs && git commit -m 'cluster: dtOO bisect job ${SLURM_JOB_ID:-manual} log'"
echo "[bisect] ================================================"
wait
