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
#SBATCH --job-name=hex3d_dtoo_export
#SBATCH --output=hex3d_dtoo_export_%j.out
#SBATCH --partition=dev_cpu_il
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

set -u

WS="${WS:-$(ws_find hex3d 2>/dev/null)}"
REPO="${REPO:-$WS/domain_partition_3D}"
# The eigenfrequencies checkout: their own scripts disagree on the location
# (cluster_env.sh defaults to $HOME/eigenfrequencies, the import docs show
# $HOME/eigen/eigenfrequencies), so probe instead of guessing. Override with
# EIGENFREQUENCIES_ROOT once the real path is known.
EF_CANDS="${EIGENFREQUENCIES_ROOT:-}"
EF_CANDS="${EF_CANDS:+$EF_CANDS }/home/st/st_us-042020/${USER}/eigen/eigenfrequencies"
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
EIG_WS="${EIG_WS:-$(ws_find eigenfreq 2>/dev/null)}"
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
echo "[dtoo] container $SQSH"

# Enroot writes unpack-state beside the image; keep it out of the checkout.
export ENROOT_DATA_PATH="${ENROOT_DATA_PATH:-$WS/enroot-data}"

t0=$(date +%s)
srun -n1 -N1 enroot start --root \
    --mount "$REPO:/repo:rw" \
    --mount "$EF_ROOT:/ef:rw" \
    "$SQSH" bash -c '
    source /usr/lib/openfoam/openfoam2606/etc/bashrc
    source /dtOO-install/bin/env.sh
    export EIGENFREQUENCIES_ROOT=/ef
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
