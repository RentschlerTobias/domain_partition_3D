# Handoff — hex3d_algohex, state and how to continue

Written for someone (or some session) starting with no context. Read this
first, then `README.md`. `PROGRESS.md` is the chronological log — long, and
only worth opening for a specific question.

---

## What this branch does

Generates a block-structured hexahedral decomposition of the T1_9 runner
passage via a genuine 3D frame field, as an alternative to extruding the 2D
surface partitions in `dp3d/`.

```
source MSH ──> tet_prep_v5.py ──> AlgoHex (frame field → IGM → HexEx) ──> hex mesh
           ──> clean_blocks.py (base complex, cleanup, sheet collapse) ──> blocks
           ──> reattach.py (blade O-grid + hub/shroud layer) ──> full domain
           ──> tfi.py (transfinite fill) ──> [not finished]
```

## Where it stands

**Recommended basis: run v11.** After post-processing:

| | value |
|---|---|
| blocks | **16**, every one with exactly 6 faces |
| cuboids | 14 (2 fail the cube-adjacency test — see "the pinch") |
| inverted cells | **0** |
| min / mean scaled Jacobian | 0.152 / 0.979 |
| `HexBlockValidator` | **VALID** |
| Hausdorff to the input surface | 0.0341 |

Full domain after re-attachment: **53 blocks, 165 356 cells, 0 inverted**,
first wall cell 0.000890 (source mesh reference: 0.000889).

Deliverables in `output/hex3d_algohex/deliverable/`, per run tag:
`T1_9_blocks_<tag>.{vtk,msh}` plus `_nfaces`, `_quality`, `_edges`, and for
v11 also `_full` and `_full_part`. The `.msh` carries the block edges as 1D
line elements (physical tag = curve id), the hexes carry `block_id`.

## Reproduce

```bash
PY=/root/repos/duty/quadmesh/.venv/bin/python
cd /root/repos/duty/quadmesh/domain_partition_3D
systemctl start docker            # AlgoHex runs in a container

$PY experimentell/hex3d_algohex/tet_prep_v5.py                 # v11 input
$PY experimentell/hex3d_algohex/run_algohex.py --tag v11 \
      --in-vtk data/T1_9/T1_9_tet_v5.vtk -- -n 60000           # ~13 min
$PY experimentell/hex3d_algohex/clean_blocks.py \
      output/hex3d_algohex/T1_9_hex_v11.ovm \
      --collapse-rounds 5 --untangle \
      --out output/hex3d_algohex/deliverable/T1_9_blocks_v11.vtk   # ~45 min
$PY experimentell/hex3d_algohex/reattach.py \
      output/hex3d_algohex/deliverable/T1_9_blocks_v11.vtk --layers 17 \
      --out output/hex3d_algohex/deliverable/T1_9_blocks_v11_full.vtk
$PY experimentell/hex3d_algohex/showcase.py --steps all --tag v11
$PY experimentell/hex3d_algohex/mesh_quality.py
```

**AlgoHex runs must be sequential** — it peaks around 6 GB and this box has
7.7. Two concurrent runs OOM (that is what killed v7 and v8).

## The one open defect: the pinch

Two of the 16 blocks touch each other at a **single shared edge** between two
of their opposite faces — a point contact, at r ≈ 0.58, z ≈ 1.64/1.73, on the
hub-side ring where the O-grid interface meets the boundary-layer interface.

It is **caused by the artificial feature rings** between the cut surfaces.
That is established, not guessed: removing them (run v14) removes the pinch
entirely — 217 of 217 blocks become cuboids. But it costs 96 inverted cells
instead of 21, 217 blocks instead of 16, and 41 blocks under 10 cells. Not
worth it. Two other attempts did **not** move the pinch: `--full-constraints`
(v12) and re-meshing the O-grid interface isotropically (v13).

Whether the pinch actually harms TFI is **still unanswered**. That question
resolves itself when the target-resolution step runs.

## Next steps, in order

1. **Target resolution.** `tfi.direction_classes` already gives the degrees of
   freedom — 11 classes over 48 block axes, validated. What is missing:
   prescribe a cell size h, solve the 11 integer counts, re-fill the blocks by
   TFI at the new counts. This is the step that turns the tooling into a
   generator; everything before it only rebuilds what AlgoHex produced.
2. **Re-measure `FRAMEFIELD_PLAN.md` step 0** afterwards — its numbers are
   from the mislabelled v9 and it carries a banner saying so. It may have
   become unnecessary; TFI defines what "good enough" means.
3. **Dataset generation** — explicitly last, by the user's instruction. See
   `ANALYSIS_PLAN.md` and the TODO at the end of `PROGRESS.md`.

## Modules

| file | role |
|---|---|
| `tet_prep_v5.py` | builds the AlgoHex input. `--remesh-ogrid H`, `--merge-interfaces`, `--keep-prisms` |
| `run_algohex.py` | Docker wrapper, `--tag` / `--in-vtk` |
| `base_complex.py` | singular edges, sheets, block partition — **purely topological** |
| `clean_blocks.py` | the big one: exact surface labels, cavity refill, block merge/split, mesh-level sheet collapse, untangling, `HexBlockValidator` |
| `reattach.py` | puts the blade O-grid and the hub/shroud layer back |
| `tfi.py` | block lattices, 3D Gordon-Hall map, conforming direction classes |
| `mesh_quality.py` | CFD metrics with OpenFOAM limits, as VTK + PNG + HTML |
| `showcase.py` | the pipeline in 11 steps |

## Things that will bite you

**Run tags are not input file numbers.** `vN` and `T1_9_tet_vN.vtk` are
independent numberings and they collide. The table in `RUNS.md` is
authoritative.

**v9's input is no longer regenerable.** A quad-diagonal bug in
`tet_prep_v5` was fixed after v9 ran; `T1_9_tet_v5.vtk` now comes out
correctly labelled. The file v9 used is kept as `T1_9_tet_v5_diagbug.vtk`.
Everything v9-based in older documentation is superseded by v11.

**More feature constraints beat fewer, five times over.** v3 diverged with
437 edges, v7/v8 OOM-died with 166, v14 produced 96 inverted cells with 310.
v5 (630) and v11 (534) are the good ones. Do not "clean up" the feature graph
without measuring.

**Block extraction and sheet collapse are topological.** They never read
coordinates, so inverted cells do not affect them — but inverted cells make
the collapse *guard* useless, which is why `untangle` now runs before the
collapse as well.

**Never classify a surface by a coordinate threshold.** It cuts across the
triangulation. Transfer labels from the input mesh by nearest face; the
machinery is `clean_blocks.SurfaceLabeller`.

**Beware counts that look like progress.** This branch has been fooled at
least six times by a metric moving the right way while the structure got
worse: cuboid *share* rising as the denominator shrank; a soft-min Jacobian
lifting the worst cell while inverting three neighbours; a kink *median*
improving while the tails worsened; non-orthogonality violations rising 4 → 56
purely because 17 layers count the same four bad columns 17 times. Prefer
counts and tail quantiles over means and extrema, and **measure the location
of a defect before proposing a cause** — several proposed fixes in the log
targeted things that were not there.

**Background jobs: never wait with `while pgrep -f <pattern>`.** The pattern
matches the waiting shell's own command line and loops forever. It cost two
silent 20-minute stalls. Launch through a script file whose name does not
appear in the waiting command.

## Documents

| file | content |
|---|---|
| `README.md` | overview, run/input mapping, "What was learned" 1-11 |
| `RUNS.md` | every AlgoHex run v1-v11 with input, runtime, outcome |
| `MESH_QUALITY.md` | the CFD criteria, measured on our own mesh |
| `TFI_RESEARCH.md` | transfinite interpolation literature + implementation plan |
| `SHOWCASE.md` | the 11-step visual analysis |
| `FRAMEFIELD_PLAN.md` | proposed frame-field stage — **numbers superseded**, see its banner |
| `POSTPROCESSING_PLAN.md`, `ANALYSIS_PLAN.md`, `PLAN.md` | earlier plans, done |
| `PROGRESS.md` | full chronological log including every failed attempt |
