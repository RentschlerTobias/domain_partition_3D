#!/bin/bash
# One canada pilot sample: dtOO (runner + guide vane) -> per passage: envelope cut -> AlgoHex n=2000
# -> beam collapse -> TFI. Runner: reduced mode (delta 1.25 x O-grid); guide vane: fluid mode
# (delta = 1/2 clearance, size 0.08). Pin the whole call with taskset to cores of its own.
# usage: run_sample.sh <sample_dir> <beam_workers>
R=$(readlink -f "$1"); BW=${2:-8}; name=$(basename "$R")
E=/opt/stack/domain_partition_3D/experimentell/hex3d_algohex; PY=/opt/venv-sci/bin/python
export PYTHONDONTWRITEBYTECODE=1
st() { local n=$1 t=$2; shift 2
       /usr/bin/time -f "[stage] $n wall=%es maxrss=%MkB exit=%x" timeout -k 60 $t "$@" > "$R/$n.log" 2>&1
       local rc=$?; echo "[$name] $(tail -1 "$R/$n.log")"
       [ $rc -eq 124 ] && echo "[$name] TIMEOUT $n after ${t}s"
       return $rc; }
rm -rf "$R/case"; cp -rL /work/canada/case "$R/case"; rm -f "$R/case/"*.log
cp /work/canada/pilot/gen_mesh_multi.py "$R/case/"
cd "$R/case" && st dtoo 2400 /opt/venv-dtoo/bin/python gen_mesh_multi.py ../params.json ru_gridGmsh:../ru.msh gv_gridGmsh:../gv.msh \
  || { echo "[$name] FAIL dtoo"; exit 1; }
cd "$E"
for comp in ru gv; do
  D="$R/$comp"; mkdir -p "$D"
  if [ $comp = ru ]; then MODE="--size-max 0.12 --hausd 0.005"; else MODE="--domain fluid --hausd 0.005 --ridges"; fi
  st ${comp}_envelope 1800 $PY envelope_cut.py "$R/$comp.msh" "$D" $MODE || { echo "[$name] FAIL ${comp}_envelope"; continue; }
  st ${comp}_algohex 7200 $PY run_algohex.py --backend native --prefix hex --tag x --in-vtk "$D/tet.vtk" --out-dir "$D" -- -n 2000 \
    || { echo "[$name] FAIL ${comp}_algohex"; continue; }
  [ -f "$D/hex_hex_x.ovm" ] || { echo "[$name] FAIL ${comp}_algohex (no ovm)"; continue; }
  st ${comp}_beam 3600 $PY beam_collapse.py "$D/hex_hex_x.ovm" --input-vtk "$D/tet.vtk" --out "$D/blocks.vtk" \
      --width 16 --depth 10 --guard struct --workers $BW || { echo "[$name] FAIL ${comp}_beam"; continue; }
  st ${comp}_tfi 900 $PY tfi.py "$D/blocks.vtk" --input-vtk "$D/tet.vtk" --target-h 0.05 --require-lattices \
      --solve-divisions --apply-divisions --out "$D/tfi.vtk" || echo "[$name] FAIL ${comp}_tfi"
  echo "[$name] ${comp} $(grep -E '^\[beam\] result' "$D/beam.log" | grep -oE "'blocks': [0-9]+, 'cuboids': [0-9]+")"
done
echo "[$name] DONE"
