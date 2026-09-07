#!/bin/bash
# Build the three AlgoHex inputs of the v15-v18 iteration, sequentially.
# Named so no waiting shell's command line matches "tet_prep" (HANDOFF.md).
set -u
cd /root/repos/duty/quadmesh/domain_partition_3D
PY=/root/repos/duty/quadmesh/.venv/bin/python
T=experimentell/hex3d_algohex/tet_prep_v5.py
E=.omo/evidence/hex3d-v15-iteration

echo "=== v9: remesh 0.035 + merge all (v15 input) ==="
$PY $T --remesh-ogrid 0.035 --merge-interfaces all --out T1_9_tet_v9.vtk

echo "=== v10: remesh 0.035 + merge hub only (v16 input) ==="
$PY $T --remesh-ogrid 0.035 --merge-interfaces hub --out T1_9_tet_v10.vtk

echo "=== v11: remesh 0.035 + ring kink 70 (v17 input) ==="
$PY $T --remesh-ogrid 0.035 --ring-kink 70 --out T1_9_tet_v11.vtk

echo "=== done ==="
