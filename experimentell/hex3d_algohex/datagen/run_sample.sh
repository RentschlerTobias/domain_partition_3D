#!/bin/bash
# One machine: dtOO -> tet_prep_v5 -> AlgoHex (n=2000) -> beam_collapse -> tfi -> export.
#   run_sample.sh <run_dir>   (run_dir holds params.json)
# Pin the whole call with taskset to cores of its own: AlgoHex ignores the thread env vars,
# so samples sharing cores oversubscribe them (see README.md). Beam search: BEAM_WORKERS.
# The dtOO case/ (with the exact blade geometry) is kept for the later blade insertion.
set -u
R=$(readlink -f "$1"); name=$(basename "$R")
E=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd); PY=${PY:-/opt/venv-sci/bin/python}
REMESH=${REMESH_OGRID:-}; BW=${BEAM_WORKERS:-2}; N=${ALGOHEX_N:-2000}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OSLO_LOCK_PATH=/tmp FOAM_SIGFPE=0
export PYTHONDONTWRITEBYTECODE=1
[ -f "$R/sample.npz" ] && { echo "skip $name"; exit 0; }
st() { local n=$1; shift; /usr/bin/time -f "[stage] $n wall=%es maxrss=%MkB exit=%x" "$@" > "$R/$n.log" 2>&1
       local rc=$?; echo "[$name] $(tail -1 "$R/$n.log")"; [ $rc -eq 0 ] || { echo "[$name] FAIL $n"; exit 1; }; }
rm -rf "$R/case"; cp -r /dtOO/demo/tistos "$R/case"; cp "$E/datagen/gen_mesh.py" "$R/case/"
cd "$R/case" && st dtoo /opt/venv-dtoo/bin/python gen_mesh.py ../params.json
cd "$E"
st tet_prep $PY tet_prep_v5.py --msh "$R/case/mesh.msh" --out-dir "$R" --out tet.vtk ${REMESH:+--remesh-ogrid $REMESH}
st algohex $PY run_algohex.py --backend native --prefix hex --tag x --in-vtk "$R/tet.vtk" --out-dir "$R" -- -n $N
st beam $PY beam_collapse.py "$R/hex_hex_x.ovm" --input-vtk "$R/tet.vtk" --out "$R/blocks.vtk" \
    --width 16 --depth 10 --guard struct --workers $BW
st tfi $PY tfi.py "$R/blocks.vtk" --input-vtk "$R/tet.vtk" --target-h 0.05 --require-lattices \
    --solve-divisions --apply-divisions --out "$R/tfi.vtk"
st export $PY export_sample.py "$R/blocks.vtk" --tet "$R/tet.vtk" --target-h 0.05 --n $N \
    --params "$R/params.json" --out "$R/sample.npz" --check
echo "[$name] OK $(grep -E '^\[beam\] result' "$R/beam.log" | grep -oE "'blocks': [0-9]+, 'cuboids': [0-9]+")"
