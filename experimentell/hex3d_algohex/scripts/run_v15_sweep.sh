#!/bin/bash
# The four AlgoHex runs of the v15-v18 iteration, STRICTLY sequential:
# ~6 GB peak on a 7.7 GB box, two at once OOM (that killed v7 and v8).
# Named so no waiting shell's command line matches "run_algohex" (HANDOFF.md).
set -u
cd /root/repos/duty/quadmesh/domain_partition_3D
PY=/root/repos/duty/quadmesh/.venv/bin/python
R=experimentell/hex3d_algohex/run_algohex.py

run () {   # tag  input  [extra flags]
  local tag=$1 in=$2; shift 2
  echo "=== $(date -Is) start $tag ($in $*) ==="
  $PY $R --tag "$tag" --in-vtk "data/T1_9/$in" -- -n 60000 "$@"
  echo "=== $(date -Is) end $tag rc=$? ==="
}

run v15 T1_9_tet_v9.vtk                        # merge all + remesh
run v16 T1_9_tet_v10.vtk                       # merge hub only + remesh
run v17 T1_9_tet_v11.vtk                       # ring kink 70 + remesh
run v18 T1_9_tet_v9.vtk --full-constraints     # merge all + remesh, constrained
echo "=== all four done ==="
