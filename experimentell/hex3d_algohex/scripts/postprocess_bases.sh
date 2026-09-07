#!/bin/bash
# Post-process one or more AlgoHex runs into block structures, sequentially.
#   postprocess_bases.sh v14u:v14 v13u:v13 v11u:v11 v15 v16 v17 v18
# "out:src" writes T1_9_blocks_<out>.vtk from T1_9_hex_<src>.ovm, so a
# re-postprocessed variant never overwrites the deliverable it is compared
# against. Labels always come from the UNMERGED T1_9_tet_v5.vtk: the labeller
# is independent of which input the run used, and reattach needs hub, shroud
# and O-grid interface told apart.
set -u
cd /root/repos/duty/quadmesh/domain_partition_3D
PY=/root/repos/duty/quadmesh/.venv/bin/python
C=experimentell/hex3d_algohex/clean_blocks.py
ROUNDS=${UNTANGLE_ROUNDS:-6}

for spec in "$@"; do
  out=${spec%%:*}; src=${spec##*:}
  echo "=== $(date -Is) postprocess $src -> T1_9_blocks_$out.vtk (untangle-rounds $ROUNDS) ==="
  $PY $C output/hex3d_algohex/T1_9_hex_$src.ovm \
        --input-vtk data/T1_9/T1_9_tet_v5.vtk \
        --collapse-rounds 5 --untangle --untangle-rounds "$ROUNDS" \
        --out output/hex3d_algohex/deliverable/T1_9_blocks_$out.vtk
  echo "=== $(date -Is) done $out rc=$? ==="
done
echo "=== all postprocessing done ==="
