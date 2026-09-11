#!/bin/bash
# dtoo_bisect2.slurm — split "mounts" from "sourcing" for the silent hang.
#
# v1 (job 6869278): step1 = `enroot start --root --rc dtOO bash -c 'exec echo ..'`
# with NO mounts -> OK in 4 s. Then step2 = same + both mounts + sourcing both
# env files -> rc=124 after 240 s with ZERO container output, and step3 (same
# preamble + python) inherited the hang.
#
# v1 never isolated WHAT hangs. The eigenfrequencies repo documents two traps
# that can look like this: a shell loop when cd precedes the OpenFOAM source
# (probe_openfoam_cd_order.sh) and mount/argv path mismatches under the
# /home -> /pfs/data6 symlink (commit 135ec9d). Our payloads do not cd before
# sourcing, so v2 measures the remaining variables directly:
#
#   step0  host checks: du rc, ownership, unreadable entries, stale mounts
#   step1  start + --rc, echo                        sanity (v1: 4 s)
#   step2  start + --rc + BOTH mounts, echo          the discriminator
#   step3  one mount at a time / without :rw         only if step2 hung
#   step4  fresh container create + fresh start      only if step2 hung
#   step5  NO mounts, source both env files, echo    sourcing alone
#   step6  BOTH mounts, source both env files, echo  v1 step2 repeat
#   step7  bounded export smoke                      only if 5 and 6 pass
#
# Each step is timeout-bounded and echoes its rc + timestamp; the FIRST rc=124
# (or kill) names the culprit. Output is tee'd to logs/dtoo_bisect2_<jobid>.log,
# committed like every other cluster log.
#SBATCH --job-name=dtoo_bisect2
#SBATCH --output=dtoo_bisect2_%j.out
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
PYTHON_MODULE="${PYTHON_MODULE:-devel/python/3.13.3-gnu-14.2}"
COUNT="${COUNT:-64}"
DATASET="${DATASET:-data/dataset/sobol}"
CONTAINER="${CONTAINER:-dtOO}"

module load "$PYTHON_MODULE" 2>/dev/null || true
export PATH="$REPO/.venv/bin:$PATH"
export ENROOT_DATA_PATH="${ENROOT_DATA_PATH:-$EIG_WS/enroot-data}"
ENROOT_TEMP_PATH="${ENROOT_TEMP_PATH:-/tmp/${USER:-$(id -un)}-enroot}"
ENROOT_RUNTIME_PATH="${ENROOT_RUNTIME_PATH:-/tmp/${USER:-$(id -un)}-enroot-run}"
mkdir -p "$ENROOT_TEMP_PATH" "$ENROOT_RUNTIME_PATH"

cd "$REPO" || exit 1
LOGDIR="$REPO/experimentell/hex3d_algohex/logs"
mkdir -p "$LOGDIR"
LOG="$LOGDIR/dtoo_bisect2_${SLURM_JOB_ID:-manual}.log"
exec > >(tee "$LOG") 2>&1

echo "[bisect2] host $(hostname), job ${SLURM_JOB_ID:-none}, $(date -Is)"
echo "[bisect2] repo $REPO"
echo "[bisect2] data $ENROOT_DATA_PATH"

SQSH="$EIG_WS/enroot-images/dtOO.sqsh"
[ -f "$SQSH" ] || { echo "[bisect2] no dtOO image: $SQSH"; exit 1; }
[ -d "$REPO/$DATASET" ] || { echo "[bisect2] no dataset root: $REPO/$DATASET"; exit 1; }
EF=""
for c in "$EIG_WS/eigenfrequencies" "$HOME/eigen/eigenfrequencies" "$HOME/eigenfrequencies"; do
    [ -f "$c/src/eigenfrequencies/adapters/dtoo/adapter.py" ] && { EF="$c"; break; }
done
[ -n "$EF" ] || { echo "[bisect2] no eigenfrequencies checkout"; ls -la "$HOME" | grep -i ein || true; exit 1; }
RC="$EF/cluster/enroot_rc.sh"
[ -f "$RC" ] || { echo "[bisect2] no rc script: $RC"; exit 1; }

echo "[bisect2] ef   $EF"
echo "[bisect2] img  $SQSH"
echo "[bisect2] enroot $(enroot version 2>/dev/null || true)"

run_step() {  # run_step <timeout_s> <label> <enroot start args...>
    local t="$1" label="$2"; shift 2
    echo "[bisect2] --- $label --- $(date -Is)"
    timeout "$t" srun -n1 -N1 enroot start --root "$@" < /dev/null
    local rc=$?
    echo "[bisect2] --- $label rc=$rc --- $(date -Is)"
    return $rc
}

# ---------------------------------------------------------------- step 0 host
echo "[bisect2] --- step0 host checks --- $(date -Is)"
du -sh "$ENROOT_DATA_PATH/$CONTAINER" 2>&1
echo "[bisect2] step0 du rc=$?"
stat -c 'owner=%U group=%G mode=%a' "$ENROOT_DATA_PATH/$CONTAINER" 2>&1 || true
echo "[bisect2] unreadable entries (maxdepth 4): $(find "$ENROOT_DATA_PATH/$CONTAINER" -maxdepth 4 ! -readable 2>/dev/null | wc -l)"
findmnt -o TARGET,SOURCE,FSTYPE 2>/dev/null | grep -Ei 'enroot|/repo|/ef' || echo "[bisect2] no stale enroot/repo/ef mounts on host"
enroot list 2>&1 | sed 's/^/[bisect2] enroot list: /'

# ---------------------------------------------------------------- step 1 sanity
run_step 60 "step1 start+rc" --rc "$RC" "$CONTAINER" \
    bash -c 'date -Is; exec echo STEP1-OK'
r1=$?

# ---------------------------------------------------------------- step 2 mounts
run_step 90 "step2 mounts-only" --rc "$RC" \
    --mount "$REPO:/repo:rw" --mount "$EF:/ef:rw" "$CONTAINER" \
    bash -c 'date -Is; exec echo MOUNT-OK'
r2=$?

if [ $r2 -ne 0 ]; then
    run_step 60 "step3a repo-mount-only" --rc "$RC" --mount "$REPO:/repo:rw" "$CONTAINER" \
        bash -c 'exec echo REPO-MOUNT-OK'
    run_step 60 "step3b ef-mount-only" --rc "$RC" --mount "$EF:/ef:rw" "$CONTAINER" \
        bash -c 'exec echo EF-MOUNT-OK'
    run_step 60 "step3c repo-mount-hostpath-no-rw" --rc "$RC" --mount "$REPO:$REPO" "$CONTAINER" \
        bash -c 'exec echo REPO-HOSTPATH-OK'

    FRESH="dtOO_fresh.${SLURM_JOB_ID:-$$}"
    echo "[bisect2] --- step4 fresh container create (720 s) name=$FRESH --- $(date -Is)"
    timeout 720 enroot create --name "$FRESH" "$SQSH"
    rc4=$?
    echo "[bisect2] --- step4 create rc=$rc4 --- $(date -Is)"
    if [ $rc4 -eq 0 ]; then
        run_step 90 "step4b fresh start+both-mounts" --rc "$RC" \
            --mount "$REPO:/repo:rw" --mount "$EF:/ef:rw" "$FRESH" \
            bash -c 'exec echo FRESH-MOUNT-OK'
    fi
    enroot remove -f "$FRESH" 2>/dev/null || true
else
    echo "[bisect2] mounts are fine; skipping step3/4"
fi

# ------------------------------------------------------------ step 5/6 sourcing
if [ $r2 -eq 0 ]; then
    run_step 90 "step5 source-no-mounts" --rc "$RC" "$CONTAINER" \
        bash -c 'source /usr/lib/openfoam/openfoam2606/etc/bashrc; source /dtOO-install/bin/env.sh; exec echo SOURCED-NOMOUNT-OK'
    r5=$?
    run_step 90 "step6 source+both-mounts" --rc "$RC" \
        --mount "$REPO:/repo:rw" --mount "$EF:/ef:rw" "$CONTAINER" \
        bash -c 'source /usr/lib/openfoam/openfoam2606/etc/bashrc; source /dtOO-install/bin/env.sh; exec echo SOURCED-OK'
    r6=$?

    if [ $r5 -eq 0 ] && [ $r6 -eq 0 ]; then
        echo "[bisect2] --- step7 export smoke (600 s) --- $(date -Is)"
        timeout 600 srun -n1 -N1 enroot start --root --rc "$RC" \
            --mount "$REPO:/repo:rw" --mount "$EF:/ef:rw" "$CONTAINER" \
            bash -c 'source /usr/lib/openfoam/openfoam2606/etc/bashrc; source /dtOO-install/bin/env.sh; export EIGENFREQUENCIES_ROOT=/ef; export PYTHONUNBUFFERED=1; cd /repo; python3.13 experimentell/hex3d_algohex/scripts/generate_machines.py --export --count '"$COUNT"'' < /dev/null
        r7=$?
        echo "[bisect2] --- step7 rc=$r7 --- $(date -Is)"
        echo "[bisect2] mesh.msh after smoke: $(find "$REPO/$DATASET" -name mesh.msh 2>/dev/null | wc -l)"
    else
        echo "[bisect2] sourcing failed (r5=$r5 r6=$r6); skipping export smoke"
    fi
fi

echo "[bisect2] ================ summary ================"
echo "[bisect2] r1=${r1:-n/a} start | r2=${r2:-n/a} mounts | rc4=${rc4:-n/a} fresh-create | r5=${r5:-n/a} source | r6=${r6:-n/a} source+mounts | r7=${r7:-n/a} export"
echo "[bisect2] first rc=124 (or kill) names the culprit; interpret top-down."
echo "[bisect2] log $LOG"
echo "[bisect2] commit it: git add experimentell/hex3d_algohex/logs && git commit -m 'cluster: dtoo bisect2 job ${SLURM_JOB_ID:-manual} log'"
wait
