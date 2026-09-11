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

RC_ALL=0
for d in "$DS"/machine_00*; do
    [ -d "$d" ] || continue
    name=$(basename "$d")
    if [ -f "$d/mesh.msh" ]; then
        echo "[dtoo] $name mesh.msh present, skipping"
        continue
    fi
    echo "[dtoo] exporting $name"
    python3.13 -u experimentell/hex3d_algohex/scripts/generate_machines.py --export --only "$name"
    rc=$?
    if [ "$rc" -ne 0 ]; then
        echo "[dtoo] $name FAILED rc=$rc"
        RC_ALL=1
    fi
done
echo "[dtoo] pass complete (rc=$RC_ALL)"
exit "$RC_ALL"
