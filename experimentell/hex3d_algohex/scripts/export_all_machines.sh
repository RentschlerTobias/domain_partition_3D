#!/usr/bin/env bash
# Container-side export loop for batch_dtoo_export.slurm.
#
# One fresh python3.13 process per machine: the dtOO adapter keeps SWIG/Gmsh
# state inside the process, and the sibling repo's turbine_runner/dtoo_cfd_build.py
# documents that repeated SWIG work in one process risks segfaults
# (labeledVectorHandling dt__mustCast). A per-machine process isolates that,
# survives single-machine failures, and gives one log line per machine.
# The same isolation makes the loop parallelizable: no two processes share a
# case copy, a machine directory or an output path, so xargs -P is race-free.
#
# Runs INSIDE the dtOO container with /repo mounted; the batch driver execs
# into this script (exec also disarms OpenFOAM's eval-of-argv double-execution).
#
# Usage (in container): bash /repo/experimentell/hex3d_algohex/scripts/export_all_machines.sh [dataset_dir] [run_id] [parallel]
#
# run_id: enroot does not propagate SLURM_JOB_ID into the container, so the
# driver injects its own host-side job id as the second argument; without it
# every run would append to the same *_manual.* progress and stage paths.
#
# parallel: number of concurrent exports (xargs -P). Each export is one
# python3.13 process with its own /tmp case stage and ~1 GB RAM; the batch
# driver passes a default sized for its 64-core / 128 GB request. The PAR
# environment variable overrides the third argument for quick experiments.
set -u

DS="${1:-data/dataset/sobol}"
RUN="${2:-manual}"
PAR="${3:-${PAR:-4}}"
cd /repo || exit 1

# Everything is also appended to a per-job progress file: stdout through the
# srun step dies with a SIGKILLed step (see the batch driver's history), while
# a file under /repo survives and can be committed from the login node later.
PROG="/repo/experimentell/hex3d_algohex/logs/export_progress_${RUN}.log"
mkdir -p "$(dirname "$PROG")"
exec > >(tee -a "$PROG") 2>&1
echo "[dtoo] progress log: $PROG (parallel=$PAR)"

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
STAGE_ROOT="/tmp/dtoo-case-${RUN}"
mkdir -p "$STAGE_ROOT"
echo "[dtoo] case source: $CASE_SRC"

# The worker runs in a child bash spawned by xargs, so everything it touches
# travels through the environment (exported function + exported variables).
export DS CASE_SRC STAGE_ROOT

export_one() {
    # The locals split across two statements on purpose: bash expands all
    # arguments of one `local` before assigning any of them, so
    # `local name="$1" d="$DS/$name"` would expand $name while still unset.
    local name="$1"
    local d="$DS/$name" stage="$STAGE_ROOT/$name" rc
    if [ -f "$d/mesh.msh" ]; then
        echo "[dtoo] $name mesh.msh present, skipping"
        return 0
    fi
    rm -rf "$stage"
    mkdir -p "$stage"
    if ! cp -rL "$CASE_SRC/." "$stage/"; then
        echo "[dtoo] $name case stage failed (source: $CASE_SRC)"
        return 1
    fi
    echo "[dtoo] exporting $name $(date -Is)"
    DTOO_CASE_DIR="$stage" python3.13 -u experimentell/hex3d_algohex/scripts/generate_machines.py --export --only "$name"
    rc=$?
    if [ "$rc" -ne 0 ]; then
        echo "[dtoo] $name FAILED rc=$rc"
        return 1
    fi
    if [ ! -f "$d/mesh.msh" ]; then
        echo "[dtoo] $name produced no mesh (see export_error.txt)"
        return 1
    fi
    return 0
}
export -f export_one

todo=()
for d in "$DS"/machine_00*; do
    [ -d "$d" ] || continue
    name=$(basename "$d")
    if [ -f "$d/mesh.msh" ]; then
        echo "[dtoo] $name mesh.msh present, skipping"
        continue
    fi
    todo+=("$name")
done
echo "[dtoo] pending: ${#todo[@]} machines, parallel=$PAR"

# xargs exit codes: 0 all child runs succeeded, 123 at least one exited 1-125.
XRC=0
if [ "${#todo[@]}" -gt 0 ]; then
    printf '%s\n' "${todo[@]}" | xargs -r -P "$PAR" -n 1 bash -c 'export_one "$@"' _
    XRC=$?
fi
rm -rf "$STAGE_ROOT"
echo "[dtoo] pass complete (xargs rc=$XRC; 0 = all ok, 123 = at least one failed)"
if [ "$XRC" -eq 0 ]; then exit 0; else exit 1; fi
