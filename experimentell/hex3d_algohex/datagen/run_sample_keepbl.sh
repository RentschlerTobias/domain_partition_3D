#!/bin/bash
# Variant of run_sample.sh with the hub/shroud boundary layer KEPT in the AlgoHex domain (only the blade O-grid is
# cut): tet_prep_v5 --keep-prisms --remesh-all 0, AlgoHex, beam collapse with a relaxed Hausdorff guard (0.1: with
# 0.05 the collapse stalls at 30 blocks on tistos_s1_0131, with 0.1 it reaches 12 blocks isomorphic to the old
# topology), tfi, export, curve_postprocess.
#   run_sample_keepbl.sh <run_dir>   (run_dir holds params.json; case/mesh.msh is reused when present)
set -u
R=$(readlink -f "$1"); name=$(basename "$R")
E=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd); PY=${PY:-/opt/venv-sci/bin/python}
BW=${BEAM_WORKERS:-4}; N=${ALGOHEX_N:-2000}; HTOL=${BEAM_HAUSDORFF:-0.1}; BWID=${BEAM_WIDTH:-32}; BDEP=${BEAM_DEPTH:-20}
TO_DTOO=${TO_DTOO:-1800} TO_TET_PREP=${TO_TET_PREP:-900} TO_ALGOHEX=${TO_ALGOHEX:-10800}
TO_BEAM=${TO_BEAM:-7200} TO_TFI=${TO_TFI:-600} TO_EXPORT=${TO_EXPORT:-600}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OSLO_LOCK_PATH=/tmp FOAM_SIGFPE=0
export PYTHONDONTWRITEBYTECODE=1
[ -f "$R/sample.npz" ] && { echo "skip $name"; exit 0; }
st() { local n=$1 t=$2; shift 2
       /usr/bin/time -f "[stage] $n wall=%es maxrss=%MkB exit=%x" timeout -k 60 $t "$@" > "$R/$n.log" 2>&1
       local rc=$?; echo "[$name] $(tail -1 "$R/$n.log")"
       [ $rc -eq 124 ] && echo "[$name] TIMEOUT after ${t}s"
       [ $rc -eq 0 ] || { echo "[$name] FAIL $n"; exit 1; }; }
if [ ! -f "$R/case/mesh.msh" ]; then
    rm -rf "$R/case"; cp -r /dtOO/demo/tistos "$R/case"; cp "$E/datagen/gen_mesh.py" "$R/case/"
    cd "$R/case" && st dtoo $TO_DTOO /opt/venv-dtoo/bin/python gen_mesh.py ../params.json
fi
cd "$E"
st tet_prep $TO_TET_PREP $PY tet_prep_v5.py --msh "$R/case/mesh.msh" --out-dir "$R" --out tet.vtk --keep-prisms --remesh-all 0
st algohex $TO_ALGOHEX $PY run_algohex.py --backend native --prefix hex --tag x --in-vtk "$R/tet.vtk" \
    --out-dir "$R" -- -n $N
st beam $TO_BEAM $PY beam_collapse.py "$R/hex_hex_x.ovm" --input-vtk "$R/tet.vtk" --out "$R/blocks.vtk" \
    --width $BWID --depth $BDEP --guard struct --hausdorff-tol $HTOL --workers $BW
st tfi $TO_TFI $PY tfi.py "$R/blocks.vtk" --input-vtk "$R/tet.vtk" --target-h 0.05 --require-lattices \
    --solve-divisions --apply-divisions --out "$R/tfi.vtk"
st export $TO_EXPORT $PY export_sample.py "$R/blocks.vtk" --tet "$R/tet.vtk" --target-h 0.05 --n $N \
    --params "$R/params.json" --out "$R/sample.npz" --check
st curves 600 $PY curve_postprocess.py "$R/sample.npz" --out-dir "$R"
echo "[$name] OK $(grep -E '^\[beam\] result' "$R/beam.log" | grep -oE "'blocks': [0-9]+, 'cuboids': [0-9]+")"
