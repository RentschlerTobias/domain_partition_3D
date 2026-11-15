#!/bin/bash
# run_ah.sh <dir> <n> <cpu list> <threads> [timeout s] [extra HexMeshing args]
# patched AlgoHex (QGP3D fix) on <dir>/tet.vtk; writes hex_n<n>.ovm, metrics_n<n>.json, hexmeshing_n<n>.log,
# done_n<n> (exit code) and ah_n<n>.pid (PID of the timeout wrapper)
d=$1; n=$2; cpus=$3; th=$4; to=${5:-25200}; shift 5
cd "$d" || exit 1
echo $$ > ah_n${n}.pid
rm -f done_n${n}
start=$(date +%s)
OMP_NUM_THREADS=$th timeout $to taskset -c $cpus /work/datagen_keepbl/algohex_fix/HexMeshing \
  -i tet.vtk -o hex_n${n}.ovm -j metrics_n${n}.json -n $n "$@" > hexmeshing_n${n}.log 2>&1
rc=$?
echo "$rc $(( $(date +%s) - start ))" > done_n${n}
