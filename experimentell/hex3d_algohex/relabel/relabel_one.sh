#!/bin/bash
# Relabel one dataset sample with the beam sheet collapse, from its stored blocks.vtk
# (no AlgoHex): relabel_one.sh <sample_dir_name> <cores> <workers>
# Writes $RELABEL_WORK/out/<name>/sample.npz (+ logs); the dataset itself is not touched.
set -u
n=$1; C=$2; W=$3
H=$(cd "$(dirname "$0")/.." && pwd)                 # experimentell/hex3d_algohex
W=${RELABEL_WORK:?set RELABEL_WORK (output: out/, locks/, lane logs)}
SRC=${DATASET_BATCH:-/opt/stack/meshtron/data/hex3d_algohex/batch}
PY=${PY:-/opt/venv-sci/bin/python}
S=$SRC/$n; D=$W/out/$n; mkdir -p $D $W/locks
E=$H
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
[ -f $D/sample.npz ] && { echo "skip $n"; exit 0; }
# claim the sample (atomic): lanes may overlap since 2026-10-05 19:40 (6 lanes)
mkdir $W/locks/$n 2>/dev/null || { echo "skip(locked) $n"; exit 0; }
N=${n##*_n}
taskset -c $C $PY $H/relabel/npz_surface_vtk.py $S/sample.npz $D/surface.vtk
$PY -c "import numpy as np,sys; open('$D/params.json','w').write(str(np.load('$S/sample.npz',allow_pickle=True)['params']))"
t0=$(date +%s)
cd $E
taskset -c $C $PY beam_collapse.py $S/blocks.vtk --input-vtk $D/surface.vtk --out $D/blocks.vtk \
  --width 16 --depth 10 --guard struct --workers $W > $D/beam.log 2>&1 || { echo "FAIL beam $n"; exit 1; }
taskset -c $C $PY tfi.py $D/blocks.vtk --input-vtk $D/surface.vtk --target-h 0.05 --require-lattices \
  --solve-divisions --apply-divisions --out $D/tfi.vtk > $D/tfi.log 2>&1 || { echo "FAIL tfi $n"; exit 2; }
taskset -c $C $PY export_sample.py $D/blocks.vtk --tet $D/surface.vtk --target-h 0.05 --n $N \
  --params $D/params.json --out $D/sample.npz --check > $D/export.log 2>&1 || { echo "FAIL export $n"; exit 3; }
rm -f $D/surface.vtk $D/blocks.msh $D/blocks_nfaces.vtk $D/blocks_quality.vtk
echo "OK $n $(( $(date +%s) - t0 ))s $(grep -oE "\[beam\] start \{'cells': [0-9]+, 'blocks': [0-9]+" $D/beam.log | grep -oE '[0-9]+$') -> $(grep -E '^\[beam\] result' $D/beam.log | grep -oE "'blocks': [0-9]+, 'cuboids': [0-9]+")"
