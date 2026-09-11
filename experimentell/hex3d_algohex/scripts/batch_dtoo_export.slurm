#!/bin/bash
# T14 pass 2 of 2: export mesh.msh per machine, inside the dtOO container.
#
#   sbatch --test-only experimentell/hex3d_algohex/scripts/batch_dtoo_export.slurm
#   COUNT=64 sbatch      experimentell/hex3d_algohex/scripts/batch_dtoo_export.slurm
#
# The draw already happened locally (generate_machines.py --preview): every
# machine_00NN/ has its params.json in the repo, and sampler.json carries the
# state. This job ONLY calls adapter.export_mesh per pending machine, which is
# why it can run inside the container without scipy.
#
# IDEMPOTENT, needed: a dtOO export has no measured cost yet, so no walltime
# is budgeted right. If the job runs out of the 30-minute dev_cpu_il cap, the
# finished mesh.msh files make a resubmit a no-op for them; export failures
# are recorded in machine_00NN/export_error.txt and count toward the yield
# number T14 wants, so resubmitting is cheap and safe in both directions.
#
# Requires (one-time, cluster/eigenfrequencies/cluster/enroot_dtoo_import.md):
#   ws eigenfreq: enroot-images/dtOO.sqsh   imported via enroot import, zstd opts
#   home: ~/eigen/eigenfrequencies          the sibling repo checkout
#   this repo: data/dataset/sobol/          from the local --preview draw
#
# Both env files BEFORE python3.13 (or: libPstream.so / libTKFeat.so.7.9 not
# found) and NO set -u inside the container (OpenFOAM's bashrc reads unset
# variables) -- carried over verbatim from submit_dtoo_enroot_smoke.sh.
#
# START BY NAME, NOT BY --root <sqsh>. `enroot start --root <sqsh>` mounts the
# squashfs through squashfuse, which decompresses every read in one user-space
# process. Measured in the eigenfrequencies repo: dtOO's CreateStates takes 8 s
# against an unpacked image and had not finished after 20 minutes under
# squashfuse (37 % CPU the whole time -- not a hang, a crawl). That is what
# job 6865391's empty 30-minute TIMEOUT was. This job unpacks the image ONCE
# into the shared "$EIG_WS/enroot-data" store and then starts the container
# "dtOO" by name; the next job, and every resume resubmission, skips the
# unpack outright. The pattern (including the race-safe install) is copied
# from eigenfrequencies cluster/submit_hydroflow_opt.sh -- their measurement,
# their fix; our log of it is docs/cluster-enroot-findings.md.
#
# --rc /ef/cluster/enroot_rc.sh: without it the dtOO image's own command
# script sources OpenFOAM, whose etc/config.sh/functions runs `eval` on the
# argv. That destroys the quoting of `bash -c '...'` and runs the command
# TWICE (their commit fae9a3d); the first, broken execution then decides the
# exit code. The rc script is one line, `exec "$@"`. Its path must resolve
# INSIDE the container, which is why it lives in the mounted eigenfrequencies
# checkout. read eigenfrequencies/cluster/enroot_rc.sh before changing this.
#SBATCH --job-name=hex3d_dtoo_export
#SBATCH --output=hex3d_dtoo_export_%j.out
#SBATCH --partition=dev_cpu_il
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

set -u

# REPO: trust the directory sbatch was invoked from (SLURM_SUBMIT_DIR) over a
# $WS that leaked in from an interactive shell that sourced another repo's
# cluster_env.sh -- that leak sent job 6865252 to look for the dataset inside
# ws eigenfreq. WS then only matters for enroot scratch, not for finding us.
REPO="${REPO:-${SLURM_SUBMIT_DIR:-$WS/domain_partition_3D}}"
# Logs survive the job and git: one file per job under the repo's logs/
# directory, then commit them (mirrors the eigenfrequencies convention, see
# experimentell/hex3d_algohex/logs/README.md). Everything printed below --
# guards included -- lands in that file, because the guards are exactly the
# evidence a failed job leaves behind.
LOGDIR="$REPO/experimentell/hex3d_algohex/logs"
mkdir -p "$LOGDIR"
LOG="$LOGDIR/batch_dtoo_export_${SLURM_JOB_ID:-manual}.log"
exec > >(tee "$LOG") 2>&1
WS="${WS:-$(ws_find hex3d 2>/dev/null)}"
# The eigenfrequencies checkout lives under the eigenfreq workspace (home
# holds only symlinks into it); resolving it via ws_find keeps the ws_ naming
# consistent with everything else on the cluster. Home paths stay as fallback
# candidates. Override with EIGENFREQUENCIES_ROOT once pinned.
EIG_WS="${EIG_WS:-$(ws_find eigenfreq 2>/dev/null)}"
EF_CANDS="${EIGENFREQUENCIES_ROOT:-}"
EF_CANDS="${EF_CANDS:+$EF_CANDS }${EIG_WS:+$EIG_WS/eigenfrequencies}"
EF_CANDS="$EF_CANDS $HOME/eigen/eigenfrequencies $HOME/eigenfrequencies"
EF_ROOT=""
for c in $EF_CANDS; do
    [ -d "$c/src" ] && { EF_ROOT="$c"; break; }
done
if [ -z "$EF_ROOT" ]; then
    echo "[dtoo] no eigenfrequencies checkout found. Probed:"
    for c in $EF_CANDS; do echo "[dtoo]   $c/src"; done
    echo "[dtoo] locate it on a login node, e.g.:  ls -d /home/st/*/${USER}/eigen/eigenfrequencies \$HOME/eigen*"
    echo "[dtoo] then resubmit with:              EIGENFREQUENCIES_ROOT=<that path> sbatch $0"
    exit 1
 fi
# The dtOO image lives in the EIGENFREQUENCIES workspace (ws eigenfreq), where
# the existing smoke-test chain already imported and uses it -- one image, one
# location, started by both pipelines. Do NOT copy it into this workspace.
# (EIG_WS is resolved above, together with the checkout candidates.)
SQSH="${SQSH:-$EIG_WS/enroot-images/dtOO.sqsh}"
COUNT="${COUNT:-64}"
DATASET="data/dataset/sobol"

module load devel/python/3.13.3-gnu-14.2 2>/dev/null || true
export ENROOT_TEMP_PATH="${ENROOT_TEMP_PATH:-/tmp/${USER:-$(id -un)}-enroot}"
export ENROOT_RUNTIME_PATH="${ENROOT_RUNTIME_PATH:-/tmp/${USER:-$(id -un)}-enroot-run}"
mkdir -p "$ENROOT_TEMP_PATH" "$ENROOT_RUNTIME_PATH"

[ -d "$REPO/$DATASET" ] || { echo "[dtoo] no dataset root: $REPO/$DATASET -- run the --preview draw first"; exit 1; }
[ -f "$SQSH" ] || { echo "[dtoo] no dtOO image: $SQSH (ws eigenfreq must hold this workspace's imported image per enroot_dtoo_import.md)"; exit 1; }
[ -d "$EF_ROOT/src" ] || { echo "[dtoo] no eigenfrequencies checkout: $EF_ROOT"; exit 1; }

cpus=$(nproc)
echo "[dtoo] host $(hostname), job ${SLURM_JOB_ID:-none}, $cpus cores"
echo "[dtoo] repo $REPO"
echo "[dtoo] image $SQSH"

# ------------------------------------------------- container: unpack once ever
# One persistent store for both pipelines, beside the image it came from, so a
# resubmission -- or the eigenfrequencies pipeline -- never unpacks twice.
export ENROOT_DATA_PATH="${ENROOT_DATA_PATH:-$EIG_WS/enroot-data}"
mkdir -p "$ENROOT_DATA_PATH"
CONTAINER="${CONTAINER:-dtOO}"
if ! enroot list 2>/dev/null | grep -qxF "$CONTAINER"; then
    # Race-safe install, copied from submit_hydroflow_opt.sh: unpack under a
    # private name, move into place atomically. The loser of a race deletes
    # its copy and uses the winner's. A container that is present but whose
    # name enroot does not list is a broken install, not a slow one.
    tmp_name="$CONTAINER.${SLURM_JOB_ID:-$$}"
    echo "[dtoo] unpacking $SQSH -> $ENROOT_DATA_PATH/$CONTAINER (one time, minutes from Lustre)"
    enroot create --name "$tmp_name" "$SQSH" \
        || { echo "[dtoo] enroot create failed: $SQSH"; exit 1; }
    if mv -T "$ENROOT_DATA_PATH/$tmp_name" "$ENROOT_DATA_PATH/$CONTAINER" 2>/dev/null; then
        echo "[dtoo] container $CONTAINER installed"
    else
        rm -rf "$ENROOT_DATA_PATH/$tmp_name"
        enroot list 2>/dev/null | grep -qxF "$CONTAINER" \
            || { echo "[dtoo] container name metadata mismatch: $ENROOT_DATA_PATH/$CONTAINER"; exit 1; }
    fi
else
    echo "[dtoo] container $CONTAINER already unpacked, skipping"
fi

t0=$(date +%s)
# --rc: OpenFOAM's eval-of-argv double-execution trap, see the header.
srun -n1 -N1 enroot start --root \
    --rc /ef/cluster/enroot_rc.sh \
    --mount "$REPO:/repo:rw" \
    --mount "$EF_ROOT:/ef:rw" \
    "$CONTAINER" bash -c '
    source /usr/lib/openfoam/openfoam2606/etc/bashrc
    source /dtOO-install/bin/env.sh
    export EIGENFREQUENCIES_ROOT=/ef
    # unbuffered: otherwise 30 min of export shows an empty .out file, which
    # is what job 6865391 did -- debugging a TIMEOUT from zero evidence.
    export PYTHONUNBUFFERED=1
    cd /repo
    python3.13 experimentell/hex3d_algohex/scripts/generate_machines.py --export --count '"$COUNT"'
'
DT=$(( $(date +%s) - t0 ))

# ------------------------------------------------------------------ harvest
PF="$REPO/$DATASET"
OK=0; MISS=0; FAIL=0
for d in "$PF"/machine_00*; do
    [ -d "$d" ] || continue
    if [ -f "$d/mesh.msh" ]; then OK=$(( OK + 1 ))
    elif [ -f "$d/export_error.txt" ]; then FAIL=$(( FAIL + 1 ))
    else MISS=$(( MISS + 1 )); fi
done
echo
echo "[dtoo] ================ yield summary ================"
echo "[dtoo] wall            ${DT}s"
echo "[dtoo] machines        $OK exported, $FAIL failed, $MISS not attempted"
[ "$OK" -gt 0 ] && echo "[dtoo] mesh.msh        $(du -sh "$PF" | cut -f1) dataset total"
echo "[dtoo] yield           $OK / $(( OK + FAIL )) attempted; with MISS >= 1, resubmit to finish"
echo "[dtoo] feed these into batch_samples.slurm:"
echo "[dtoo]   NS=\"2000 8000\" CANDS=\"$DATASET/machine_00*\" \\"
echo "[dtoo]       sbatch experimentell/hex3d_algohex/scripts/batch_samples.slurm"
echo "[dtoo] ============================================"
[ "$MISS" -gt 0 ] && echo "[dtoo] $MISS machines have neither mesh.msh nor export_error.txt: job may have been cut before reaching them."

# Persist the evidence: the log is the run's record and belongs in git, like
# eigenfrequencies' cluster/logs/<jobid>/. The yield number T14 needs is in it.
echo "[dtoo] keep this log:  git add experimentell/hex3d_algohex/logs && git commit -m 'cluster: dtOO export job ${SLURM_JOB_ID:-manual} log'"
echo "[dtoo] log file:       $LOG"
# let the tee process drain before the shell exits, or the tail is lost
wait
