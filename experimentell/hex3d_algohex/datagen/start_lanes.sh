#!/bin/bash
# Start one lane per core pair: start_lanes.sh  (LANES default "0,1 2,3 4,5"; log lanes.log)
#   DATAGEN_WORK=<dir with runs/<name>/params.json> LANES="0,1 2,3 4,5" start_lanes.sh
# Stop: touch $DATAGEN_WORK/STOP (running samples finish). Resume: rm STOP and start again.
H=$(cd "$(dirname "$(readlink -f "$0")")" && pwd); cd ${DATAGEN_WORK:?set DATAGEN_WORK}
rm -f lane_*.done
for c in ${LANES:-0,1 2,3 4,5}; do "$H/run_lane.sh" $c >> lanes.log 2>&1 & done
wait; echo "[$(date +%F' '%T)] all lanes done" >> lanes.log
