#!/usr/bin/env bash
# Container-side marker trace for the dtOO start diagnosis (bisect5).
#
# Runs INSIDE the dtOO container and writes stage-by-stage evidence to a
# file under the mounted repo. Reason: when the outer srun step is SIGKILLed
# at the job time limit, its stdout is lost (jobs 6865391, 6870688), while a
# file survives and stays inspectable from the host afterwards.
#
# Usage (inside the container): bash probe_markers.sh [machine_name]
# Default marker file: /repo/experimentell/hex3d_algohex/logs/probe_markers5.txt
# Override with MARKER_FILE=<path>.
set -u

M=${MARKER_FILE:-/repo/experimentell/hex3d_algohex/logs/probe_markers5.txt}
MACHINE=${1:-machine_0001}
DS=data/dataset/sobol

mark() { echo "$(date -Is 2>/dev/null || echo no-date) $*" >> "$M"; }

mark "container alive: $(uname -n), python3.13=$(command -v python3.13 || echo MISSING)"
cd /repo || { mark "FAIL: cd /repo"; exit 3; }
mark "cd /repo OK"
mark "gen script present: $( [ -f experimentell/hex3d_algohex/scripts/generate_machines.py ] && echo yes || echo NO )"

timeout 480 python3.13 -u experimentell/hex3d_algohex/scripts/generate_machines.py --export --only "$MACHINE" >> "$M" 2>&1
rc=$?
mark "generate_machines --only $MACHINE rc=$rc"

if [ -f "$DS/$MACHINE/mesh.msh" ]; then
    mark "mesh.msh present: $(ls -l "$DS/$MACHINE/mesh.msh")"
elif [ -f "$DS/$MACHINE/mesh.msh.part" ]; then
    mark "mesh.msh.part present (killed mid-write or export still running)"
else
    mark "no mesh file written"
fi

exit "$rc"
