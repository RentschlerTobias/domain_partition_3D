# `blocks_core.tar.zst` — the structures that cannot be recomputed cheaply

`output/` and `data/T1_9/T1_9_tet*.vtk` are gitignored, and this box has no
off-box store: no rclone, restic or borg, no git-lfs, no configured SSH host.
GitHub is the only one, so the subset that is both irreplaceable and small
lives here as a tracked archive — 22.4 MiB of ASCII VTK compressed to 5.5 MiB
by `zstd -19`.

```
sha256  94666a40037552276be307853410c20d997d6bbe75aa9d1262e3bab1589b09c4
```

## Restore

Paths inside the archive are repo-relative, so extracting from the repo root
puts every file back where the code looks for it:

```bash
tar -I zstd -xf experimentell/hex3d_algohex/fixtures/blocks_core.tar.zst
```

Then the tests run — with the project interpreter, not the host `python`, which
has no `meshio`:

```bash
PY=/root/repos/duty/quadmesh/.venv/bin/python
$PY -m pytest experimentell/hex3d_algohex/tests/ -q
```

## Rebuild

```bash
D=output/hex3d_algohex/deliverable
tar cf - \
  $D/T1_9_blocks_n2000.vtk $D/T1_9_blocks_n8000.vtk \
  $D/T1_9_blocks_v11.vtk $D/T1_9_blocks_v11m.vtk $D/T1_9_blocks_v16m.vtk \
  $D/T1_9_blocks_*.divisions.json \
  data/T1_9/T1_9_tet_v5.vtk data/T1_9/T1_9_tet_v11.vtk \
  | zstd -19 -T2 -o experimentell/hex3d_algohex/fixtures/blocks_core.tar.zst
```

## What is in it, and who needs it

| file | size | needed by |
|---|---|---|
| `T1_9_blocks_n2000.vtk` | 0.15 MB | T6, the `-n` measurement |
| `T1_9_blocks_n8000.vtk` | 0.56 MB | T6 |
| `T1_9_blocks_v11.vtk` | 4.59 MB | T6 (as the n=60000 structure), T7, and **both tests** — `tests/test_divisions.py` and `tests/test_refill.py` hardcode this path |
| `T1_9_blocks_v11m.vtk` | 6.37 MB | T12, the 12-block question |
| `T1_9_blocks_v16m.vtk` | 6.91 MB | the basis comparison in `basis_report.py` |
| `T1_9_blocks_*.divisions.json` | 0.01 MB | the solved conforming division counts, minutes of MILP each |
| `T1_9_tet_v5.vtk` | 2.24 MB | the AlgoHex input behind v11 (`RUNS.md`), and the input T3 runs through enroot |
| `T1_9_tet_v11.vtk` | 2.62 MB | the input behind v16/v17 |

The tet inputs are reproducible via `tet_prep_v5.py`, but they are the anchor
for every measurement in `RUNS.md`, and 4.9 MB is not worth the risk of a
re-meshing that shifts a vertex.

## What is deliberately NOT in it

Everything else under `output/`. Losing it costs compute, not information:
about 12 min of AlgoHex per structure plus the sheet collapse, which is 1 / 3 /
45 min at n=2000 / 8000 / 60000. That is cluster work in any case, and the
243 MB of `*_hex*.ovm` still on disk skip the AlgoHex half of it.

Two things this archive does not protect, by design: the ~1 MB per-sample
neutral sample files of T9 (they do not exist yet, and once they do they belong
to the dataset, not to a fixture), and the generative samples in `output/.../gen/`.
Add the latter here only if a measurement starts depending on a specific one.
