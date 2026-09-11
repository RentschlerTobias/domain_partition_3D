#!/usr/bin/env bash
# One-machine dtOO export probe, foreground, live output. For diagnosing the
# empty 30-min batch runs (6865391): shows container import, env.sh sourcing
# and the first export_mesh call with PYTHONUNBUFFERED.
#   bash experimentell/hex3d_algohex/scripts/dtoo_probe_one.sh [count]
set -u
CW=$(cd "$(dirname "$0")" && pwd)
REPO="${REPO:-$(cd "$CW/../../.." && pwd)}"
EIG_WS="${EIG_WS:-$(ws_find eigenfreq 2>/dev/null)}"
EF_ROOT="${EF_ROOT:-}"
if [ -z "$EF_ROOT" ]; then
    for c in "$EIG_WS/eigenfrequencies" "$HOME/eigen/eigenfrequencies" "$HOME/eigenfrequencies"; do
        [ -d "$c/src" ] && { EF_ROOT="$c"; break; }
    done
fi
[ -d "$EF_ROOT/src" ] || { echo "no checkout, probed: $EIG_WS/eigenfrequencies \$HOME/eigen/eigenfrequencies \$HOME/eigenfrequencies"; ls -la "$HOME" | grep -i ein; exit 1; }
SQSH="${SQSH:-$EIG_WS/enroot-images/dtOO.sqsh}"
COUNT="${1:-1}"
export ENROOT_TEMP_PATH="${ENROOT_TEMP_PATH:-/tmp/$USER-enroot}"
export ENROOT_RUNTIME_PATH="${ENROOT_RUNTIME_PATH:-/tmp/$USER-enroot-run}"
mkdir -p "$ENROOT_TEMP_PATH" "$ENROOT_RUNTIME_PATH"
[ -f "$SQSH" ] || { echo "no image: $SQSH"; exit 1; }
[ -d "$EF_ROOT/src" ] || { echo "no checkout: $EF_ROOT"; exit 1; }
# --rc takes a HOST path (enroot bind-mounts the file); /ef/... is the
# container-side mapping and does NOT exist here (job 6867987).
RC_SCRIPT="$EF_ROOT/cluster/enroot_rc.sh"
[ -f "$RC_SCRIPT" ] || { echo "no --rc command script: $RC_SCRIPT"; exit 1; }

# Same container discipline as batch_dtoo_export.slurm: shared persistent
# store, start by name, --rc against OpenFOAM's eval trap. Read the header of
# batch_dtoo_export.slurm for the measurements behind this.
export ENROOT_DATA_PATH="${ENROOT_DATA_PATH:-$EIG_WS/enroot-data}"
mkdir -p "$ENROOT_DATA_PATH"
CONTAINER="${CONTAINER:-dtOO}"
if ! enroot list 2>/dev/null | grep -qxF "$CONTAINER"; then
    tmp_name="$CONTAINER.${SLURM_JOB_ID:-$$}"
    echo "[probe] unpacking $SQSH -> $ENROOT_DATA_PATH/$CONTAINER (one time)"
    enroot create --name "$tmp_name" "$SQSH" && {
        mv -T "$ENROOT_DATA_PATH/$tmp_name" "$ENROOT_DATA_PATH/$CONTAINER" 2>/dev/null \
            || rm -rf "$ENROOT_DATA_PATH/$tmp_name"
    }
fi

ENROOT_CMD="enroot start --root --rc $RC_SCRIPT \
    --mount \"$REPO:/repo:rw\" --mount \"$EF_ROOT:/ef:rw\" \
    \"$CONTAINER\" bash -c \"
source /usr/lib/openfoam/openfoam2606/etc/bashrc
source /dtOO-install/bin/env.sh
export EIGENFREQUENCIES_ROOT=/ef PYTHONUNBUFFERED=1
cd /repo
python3.13 experimentell/hex3d_algohex/scripts/generate_machines.py --export --count $COUNT
\""

if [ -n "${SLURM_JOB_ID:-}" ]; then
    # already inside a salloc/sbatch allocation: run the container directly,
    # a second srun -p here would request a NEW allocation from the login
    # queue and sit there silently (measured on 2026-09-11)
    bash -c "$ENROOT_CMD"
else
    srun -p dev_cpu_il -t 10 --cpus-per-task=8 --mem=16G -n1 -N1 bash -c "$ENROOT_CMD"
fi
