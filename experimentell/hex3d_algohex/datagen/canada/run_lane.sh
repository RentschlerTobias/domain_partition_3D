#!/bin/bash
# one lane: claim samples via mkdir locks, run them pinned to the given cores; stop with `touch STOP`
CORES=$1; P=/work/canada/pilot; nc=$(echo $CORES | tr ',' '\n' | awk -F- '{s+= ($2==""?1:$2-$1+1)} END{print s}')
for d in $P/runs/canada_p_*; do
  [ -e $P/STOP ] && break
  n=$(basename $d); mkdir $P/locks/$n 2>/dev/null || continue
  taskset -c $CORES bash $P/run_sample.sh $d $nc >> $P/lanes.log 2>&1
done
touch $P/lane_$CORES.done
