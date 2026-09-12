#!/usr/bin/env bash
# Container-side export loop for batch_dtoo_export.slurm.
#
# One fresh python3.13 process per machine: the dtOO adapter keeps SWIG/Gmsh
# state inside the process, and the sibling repo's turbine_runner/dtoo_cfd_build.py
# documents that repeated SWIG work in one process risks segfaults
# (labeledVectorHandling dt__mustCast). A per-machine process isolates that,
# survives single-machine failures, and gives one log line per machine.
#
# Runs INSIDE the dtOO container with /repo mounted; the batch driver execs
# into this script (exec also disarms OpenFOAM's eval-of-argv double-execution).
#
# Usage (in container): bash /repo/experimentell/hex3d_algohex/scripts/export_all_machines.sh [dataset_dir]
set -u

DS="${1:-data/dataset/sobol}"
cd /repo || exit 1

# Everything is also appended to a per-job progress file: stdout through the
# srun step dies with a SIGKILLed step (see the batch driver's history), while
# a file under /repo survives and can be committed from the login node later.
PROG="/repo/experimentell/hex3d_algohex/logs/export_progress_${SLURM_JOB_ID:-manual}.log"
mkdir -p "$(dirname "$PROG")"
exec > >(tee -a "$PROG") 2>&1
echo "[dtoo] progress log: $PROG"

# dtOO's case tree is baked into the container at /dtOO (every cluster config
# in the sibling repo uses /dtOO/build/test/tistos); the machine yaml's
# ~/dtOO default resolves to /root/dtOO in the container, which does not
# exist and raises FileNotFoundError. export.py documents DTOO_CASE_DIR as
# the override, so use it to point at a WRITABLE copy.
#
# The copy is required: dtOO chdirs into the case directory and opens
# machineSave.xml ReadWrite, but the case tree inside the image sits on the
# enroot rootfs, which is mounted read-only ("Failed to open fileName =
# machineSave.xml", dtXmlParser::checkFile). The sibling repo fixes this the
# same way (_stage_case_dir): stage a fresh case per machine so a leftover
# copy cannot carry the previous machine's written-back state, and
# dereference symlinks with cp -L since they point into the read-only image.
CASE_SRC="${DTOO_CASE_DIR:-/dtOO/build/test/tistos}"
STAGE_ROOT="/tmp/dtoo-case-${SLURM_JOB_ID:-manual}"
echo "[dtoo] case source: $CASE_SRC"

RC_ALL=0
for d in "$DS"/machine_00*; do
    [ -d "$d" ] || continue
    name=$(basename "$d")
    if [ -f "$d/mesh.msh" ]; then
        echo "[dtoo] $name mesh.msh present, skipping"
        continue
    fi
    stage="$STAGE_ROOT/$name"
    rm -rf "$stage"
    mkdir -p "$stage"
    if ! cp -rL "$CASE_SRC/." "$stage/"; then
        echo "[dtoo] $name case stage failed (source: $CASE_SRC)"
        RC_ALL=1
        continue
    fi
    echo "[dtoo] exporting $name $(date -Is)"
    DTOO_CASE_DIR="$stage" python3.13 -u experimentell/hex3d_algohex/scripts/generate_machines.py --export --only "$name"
    rc=$?
    if [ "$rc" -ne 0 ]; then
        echo "[dtoo] $name FAILED rc=$rc"
        RC_ALL=1
    fi
    if [ ! -f "$d/mesh.msh" ]; then
        echo "[dtoo] $name produced no mesh (see export_error.txt)"
        RC_ALL=1
    fi
done
rm -rf "$STAGE_ROOT"
echo "[dtoo] pass complete (rc=$RC_ALL)"
exit "$RC_ALL"
