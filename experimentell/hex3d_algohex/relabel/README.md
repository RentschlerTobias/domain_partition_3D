# Beam relabel of the n2000 dataset (2026-10-05)

Environment: `RELABEL_WORK` (output dir, required), `DATASET_BATCH` (default
/opt/stack/meshtron/data/hex3d_algohex/batch), `MESHTRON_ROOT` (for audit.py), `PY`.

Applies the beam sheet collapse (`beam_collapse.py`, width 16 / depth 10, struct guard) to
all 369 `machine_*_n2000` samples of `/opt/stack/meshtron/data/hex3d_algohex/batch`, starting
from their stored `blocks.vtk` (no AlgoHex re-run). The dataset itself is not touched.
n8000 samples were dropped (same geometries, ~1 h per sample); n1000/n1500 are not used either:
different AlgoHex n only change the initial topology, after the collapse they coincide.

| file | purpose |
|---|---|
| `$RELABEL_WORK/all_n2000.txt` | the 369 sample names (`ls $DATASET_BATCH \| grep _n2000`) |
| `npz_surface_vtk.py` | labelled surface of a sample.npz -> AlgoHex input VTK |
| `relabel_one.sh <name> <cores> <workers>` | one sample -> `$RELABEL_WORK/out/<name>/` (blocks.vtk, sample.npz, tfi.vtk, logs); atomic claim in `locks/` |
| `run_lane.sh <k>`, `run_lane_rev.sh <k> <cores>` | lanes (forward every 4th sample on cores 2k,2k+1 / backward); write `lane*.done` |
| `retfi.sh <name> <core>` | TFI refill for samples exported without it |
| `audit.py` | -> `audit.csv`: ok flag, reason codes, blocks before/after, topo code |

Result (`audit.csv`): 369/369 relabelled, 351 ok, 18 rejected (beam:inverted 10,
beam:1_non_cuboid 7, beam:excess_faces 7, topo:ValueError 4, blocks:mixed_orientation 3).
Blocks before: 22 x223, 12 x70, 16 x55, 75 x10, rest 36-66. After (ok): 12 x338, 10 x5,
42 x4, 16 x3, 36 x1. Canonical topology `12bl/49a0142abb` (rows 3-2-3-3-1): 337/351.
The 42/36/16-block samples are unfinished collapses, the 10-block ones over-collapsed;
training uses only `49a0142abb` (meshtron reports/blockgen_v2_canonical_best_case.md).
