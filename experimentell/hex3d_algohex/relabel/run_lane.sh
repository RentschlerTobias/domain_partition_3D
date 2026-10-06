#!/bin/bash
# One relabel lane: run_lane.sh <lane 0..3>; takes every 4th sample of $RELABEL_WORK/all_n2000.txt
# (e.g. ls $DATASET_BATCH | grep _n2000 > $RELABEL_WORK/all_n2000.txt),
# pinned to cores 2k,2k+1, 1 h budget per sample; writes lane<k>.done at the end.
k=$1; D=$(cd "$(dirname "$0")" && pwd); cd ${RELABEL_WORK:?set RELABEL_WORK}
awk -v k=$k 'NR % 4 == k' all_n2000.txt | while read n; do
  timeout 3600 $D/relabel_one.sh $n $((2*k)),$((2*k+1)) 2 || echo "FAIL/TIMEOUT $n rc=$?"
done
touch lane$k.done
