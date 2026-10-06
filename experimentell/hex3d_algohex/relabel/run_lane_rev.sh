#!/bin/bash
# Extra relabel lane working all_n2000.txt from the END: run_lane_rev.sh <id> <cores>
id=$1; C=$2; D=$(cd "$(dirname "$0")" && pwd); cd ${RELABEL_WORK:?set RELABEL_WORK}
tac all_n2000.txt | awk -v k=$id 'NR % 2 == k' | while read n; do
  timeout 3600 $D/relabel_one.sh $n $C 2 || echo "FAIL/TIMEOUT $n rc=$?"
done
touch lane_rev$id.done
