# New tistos samples: dtOO -> AlgoHex -> beam collapse (2026-10-06)

Generates new training samples for the block generator from fresh design vectors.
One AlgoHex n per geometry (n=2000): different n only change the initial AlgoHex
topology; after the beam collapse they end in the same canonical structure.
The samples contain only the core blocks; the blade is inserted after generation
(docs/decisions/2026-10-06-blade-inserted-after-generation.md). The dtOO `case/`
directory (exact blade geometry) is kept per sample for that step.

| file | purpose |
|---|---|
| `draw_new.py <root> <n> <seed> <prefix>` | scrambled Sobol in the 30-D tistos box (bounds from /dtOO/demo/tistos/build.py) -> `<root>/<prefix>_<i>/params.json`; prints the distance to the existing dataset geometries (`DATASET_BATCH`) |
| `gen_mesh.py` | dtOO build of one design vector -> case/mesh.msh (copied into the case by run_sample.sh) |
| `run_sample.sh <run_dir>` | dtOO -> tet_prep_v5 -> AlgoHex -> beam_collapse -> tfi -> export_sample (sample.npz); skips finished runs |
| `run_lane.sh <cores>` | one lane over `$DATAGEN_WORK/runs/*`, atomic claims in `locks/`, stops on `STOP` |
| `start_lanes.sh` | one lane per core pair (`LANES`, default `"0,1 2,3 4,5"`), log `lanes.log` |

    python draw_new.py $DATAGEN_WORK/runs 400 20261006 tistos_s1
    DATAGEN_WORK=<dir> LANES="0,1 2,3 4,5" ./start_lanes.sh      # stop: touch $DATAGEN_WORK/STOP

The audit (../relabel/audit.py) and the dataset builder of meshtron
(scripts/blockgen_v2/build_dataset.py) read the same sample.npz layout.

## Measured (courier, 2 own cores per sample)
dtOO ~110 s, tet_prep ~5 s, AlgoHex 520-1450 s (mean ~16 min, 2.0-2.5 GB RAM), beam 190-290 s,
tfi + export ~2 s -> ~20 min per sample, ~9 samples/h on 6 cores, ~3 per hour per 2-core lane.
First batch (seed 20261006, prefix tistos_s1, 400 drawn, interrupted after 11): all 11 finished
samples 12 blocks / 12 cuboids; 1 AlgoHex segfault (exit -11 after parametrisation, ~10 % expected).
The new design vectors are as far from the 378 existing geometries (median 1.53 box units) as
those are from each other (1.55).

## Pitfalls
- AlgoHex ignores OMP/OpenBLAS thread variables. Five samples sharing six cores (one xargs pool
  pinned to 0-5) gave load ~230, AlgoHex 41-91 min and beam 4-39 min. Give every sample its own
  cores (one lane per core pair).
- `gen_mesh.py` must be copied into the dtOO case before the build; run_sample.sh does that.
