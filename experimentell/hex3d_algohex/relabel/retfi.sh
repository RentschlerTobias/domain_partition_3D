#!/bin/bash
# Re-run the TFI refill for already relabelled samples that were exported without it.
# retfi.sh <name> <core>
n=$1; c=$2
H=$(cd "$(dirname "$0")/.." && pwd)                 # experimentell/hex3d_algohex
W=${RELABEL_WORK:?set RELABEL_WORK (output: out/, locks/, lane logs)}
SRC=${DATASET_BATCH:-/opt/stack/meshtron/data/hex3d_algohex/batch}
PY=${PY:-/opt/venv-sci/bin/python}
D=$W/out/$n; S=$SRC/$n; E=$H
export OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
grep -q "refilled mesh" $D/tfi.log 2>/dev/null && exit 0
taskset -c $c $PY $H/relabel/npz_surface_vtk.py $S/sample.npz $D/surface.vtk
cd $E && taskset -c $c $PY tfi.py $D/blocks.vtk --input-vtk $D/surface.vtk --target-h 0.05 \
  --require-lattices --solve-divisions --apply-divisions --out $D/tfi.vtk > $D/tfi.log 2>&1
rm -f $D/surface.vtk
