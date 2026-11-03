# Porting the block-generation pipeline to a cluster

Status: 2026-11-03. Branch `feat/keep-bl` (tistos keep-BL work) and `feat/canada`
(canada case) of `domain_partition_3D`.

## What has to move

| Piece | Where it is today | Notes |
|---|---|---|
| Code | this repo | `experimentell/hex3d_algohex/` is the working code, `dp3d/` the library |
| Python env | `/opt/venv-sci` (Python 3.12.3) | `pip install -r requirements.txt`; verified numpy 2.4.4, meshio 5.3.5, gmsh 4.15.2 |
| dtOO env | `/opt/venv-dtoo` | only needed for `datagen/gen_mesh.py` (case/mesh.msh) |
| gmsh | `/usr/local/bin/gmsh` | 4.15.2 |
| AlgoHex | `/opt/algohex` (binary), `/opt/algohex-src` | 2.0-2.5 GB RAM per sample, CPU only |
| mmg | `/work/tools/mmg-install/bin/{mmg2d_O3,mmg3d_O3,mmgs_O3}` | optional remeshing |
| Data | `/work/datagen/...` | one 16 GB archive, see `docs/DATA.md` |

The Apptainer image itself is **not** in the repo; rebuild it from the list above or
reuse the image the data was produced with.

## Paths and environment

The pipeline code takes explicit paths on the command line. The *defaults* still point
at the old sandbox layout, so on a new machine either pass paths explicitly or export:

| Variable | Default | Meaning |
|---|---|---|
| `DP3D_DATA_ROOT` | `/work` | root of the working data |
| `DATAGEN_WORK` | `$DP3D_DATA_ROOT/datagen` | lane directory containing `runs/` |
| `DP3D_KEEPBL_ROOT` | `$DP3D_DATA_ROOT/datagen_keepbl` | keep-BL case directory |

The datagen lane scripts (`experimentell/hex3d_algohex/datagen/run_lane.sh`,
`start_lanes.sh`, `run_sample.sh`, `run_sample_keepbl.sh`) already read
`DATAGEN_WORK`; they contain no absolute paths.

Scripts that still carry a `/work/...` default (pass the path explicitly on the cluster):

- `blade_check.py`, `audit_tfi_blade.py`, `grading_check.py`, `viz_inverted_cells.py`
- `viz_merge_curves.py`, `viz_singular.py`, `export_features_vtk.py`, `p1_export_vtk.py`
- `canada_*.py`, `calib_*.py`, `collapse_probe.py`, `curve_merge.py`, `sheet_beam.py`
- `datagen/check_keepbl_sample.py`, `datagen/draw_new.py`

Find them again with:

    git grep -n '/work/' -- '*.py' '*.sh'

## Data

Copy the archive in `docs/DATA.md` (16 GB, SHA256 recorded), extract it so that
`runs/` lands under `$DATAGEN_WORK`, and point `DATAGEN_WORK` at it. Nothing else is
required to re-run the keep-BL lane; the keep-BL case itself
(`datagen_keepbl/runs/tistos_s1_0131`, 388 MB) is not in the archive and has to be
copied separately.

## Smoke test on the new machine

    export DP3D_DATA_ROOT=/path/to/work
    export DATAGEN_WORK=$DP3D_DATA_ROOT/datagen
    cd experimentell/hex3d_algohex
    python blade_check.py s1 --run $DP3D_DATA_ROOT/datagen/runs/tistos_s1_0131
    python grading_check.py --run-dir $DP3D_DATA_ROOT/datagen_keepbl/runs/tistos_s1_0131
    bash -n datagen/run_sample_keepbl.sh

`blade_check.py s1` must report 16 blocks, lattices present, 0 cells with signed
volume <= 0 and 105/105 shared interface vertices. If it does, the environment is
complete.

## Cluster specifics

- Everything is CPU-only; AlgoHex dominates the runtime (~20 min per design on
  2 cores, ~9 samples/h on 6 cores).
- The tistos case needs `dtOO` only for fresh geometries; existing runs already carry
  `case/mesh.msh`.
- `/work` here is container-local and ephemeral; on the cluster put the data on the
  shared filesystem and keep `DP3D_DATA_ROOT` pointing there.
