#!/bin/bash
# One data-generation lane: run_lane.sh <cores e.g. 0,1>. Walks $DATAGEN_WORK/runs/* in order and
# claims each sample atomically (mkdir locks/<name>); skips finished samples. Each lane owns its
# cores, so AlgoHex (which ignores thread env vars) never shares cores with another sample.
# Stop: touch STOP (the current sample finishes).
H=$(cd "$(dirname "$(readlink -f "$0")")" && pwd); cd ${DATAGEN_WORK:?set DATAGEN_WORK}; c=$1
mkdir -p locks
for d in runs/*/; do
  d=${d%/}
  [ -f STOP ] && break
  n=$(basename $d); [ -f $d/sample.npz ] && continue
  mkdir locks/$n 2>/dev/null || continue
  taskset -c $c $H/run_sample.sh $d > $d/run.log 2>&1; tail -1 $d/run.log
done
touch lane_$c.done
