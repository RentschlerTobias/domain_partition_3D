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
           ──> tfi.py (prescribe h, solve counts, refill) ──> hex mesh
```

**`tfi.py` is a generator now.** Prescribe a cell size and it solves the
conforming division counts and refills every block:

```bash
$PY experimentell/hex3d_algohex/tfi.py \
      output/hex3d_algohex/deliverable/T1_9_blocks_v11.vtk \
      --target-h 0.05 --apply-divisions --out .../T1_9_blocks_v11_h0.05.vtk
$PY experimentell/hex3d_algohex/reattach.py .../T1_9_blocks_v11_h0.05.vtk \
      --layers 17 --out .../T1_9_blocks_v11_h0.05_full.vtk
```

Measured on v11: core 54 460 → 49 550 cells at h = 0.05, watertight, 0
inverted, min scaled Jacobian 0.1357; assembled 161 619 cells in 53 blocks, 0
inverted. Against the on-disk 17-layer v11 full mesh (165 356 cells) the CFD
metrics tie on five of seven — non-orthogonality 56 violations against 56 —
and pick up two each on aspect ratio (2 → 4) and face flatness (3 → 5), on a
mesh 2 % coarser. Refilling at the counts a block already has returns every
boundary vertex to within 2.2e-15, which is the fold gate `TFI_RESEARCH.md`
asks for.

## Where it stands

**Reference structure: `T1_9_blocks_v11m6` — the v11 basis merged to 6 blocks.**
Decided 2026-09-08 together with what the structure is *for*: it is training
data for a transformer, filled by TFI in post-processing and then run as a CFD
mesh, which makes **block count the primary criterion**. The reasoning, the
alternatives and the measurements are in
`docs/decisions/2026-09-08-hex3d-block-structure-objective.md`.

| structure | blocks | smallest | pinch | inverted | min sJ | Hausdorff | validator |
|---|---|---|---|---|---|---|---|
| v11 (raw basis) | 16 | 224 | 2 | 0 | 0.1524 | 0.03405 | VALID |
| v11m | 12 | 476 | 0 | 0 | 0.1524 | 0.03405 | VALID |
| **v11m6 (reference)** | **6** | 448 | 0 | 0 | 0.1524 | 0.03405 | VALID |
| v16m | 26 | 51 | 0 | 0 | 0.0159 | 0.01377 | VALID |

Merging changes no geometry, so every quality column is inherited unchanged
from the basis; only the topology moves. Every block is a valid TFI lattice.
`v16m` is kept as the low-boundary-error alternative — it is the only thing
v16 still wins, and decision B below is aimed at that gap.

Regenerate with:

```bash
$PY experimentell/hex3d_algohex/merge_ilp.py \
      output/hex3d_algohex/deliverable/T1_9_blocks_v11.vtk --seed all \
      --apply output/hex3d_algohex/deliverable/T1_9_blocks_v11m6.vtk
```

**The v11 basis itself** is unchanged and still the input to all of this:

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

## The pinch, and the basis question it opens

Two of v11's 16 blocks touch at a **single shared edge** between two of their
opposite faces — a point contact at r ≈ 0.58, z ≈ 1.64/1.73, where the O-grid
interface meets the boundary-layer interface on the hub side. `clean_blocks.py
--detect-pinch BLOCKS_VTK` reports it with its location and exits 1; it is a
defect of the face partition, not of the mesh (the cells around every edge of
those blocks form a single fan, checked).

**The pinch is now removable at its cause, and cheaply.** It is created by the
artificial rings between the cut surfaces, but the two rings are separable and
only the root one matters:

| | v11 | v14 | **v16** |
|---|---|---|---|
| rings | both kept | both deleted | root ring deleted |
| feature edges | 534 | 310 | **422** |
| pinch | 2 blocks | none | **none** |
| raw inverted | 21 | 96 | **32** |
| blocks after collapse | **16** | 217 | 103 |
| blocks < 10 cells | **0** | 50 | 9 |
| Hausdorff | 0.0341 | 0.0233 | **0.0138** |
| validator | VALID | INVALID | **VALID** |

`tet_prep_v5.py --merge-interfaces hub --remesh-ogrid 0.035` builds that input.
Dropping only the ring edges whose dihedral kink is below 70° (`--ring-kink
70`, run v17) also removes the pinch and keeps 499 edges, but leaves 3
inverted cells and does not validate.

**Which basis to build on is a real choice, not a formality.** v16 is the only
valid pinch-free basis and halves the boundary error; v11 has six times fewer
blocks, none of them tiny, and a much better worst cell (0.152 against 0.016).
For a hand-editable topology v11 still wins; for boundary fidelity and a clean
complex, v16 does. Full table: `output/hex3d_algohex/basis_scorecard.md`,
regenerate with `basis_report.py --map scripts/bases.json`.

One hypothesis was **refuted** on the way: v14's collapse does not stall
because its inverted cells disarm the guard. Iterating the untangle takes it
from 17 to 13 inverted and its worst cell from −0.81 to −0.30, and it still
stops at exactly 217 blocks. After two rounds no sheet strictly improves
(excess faces, block count) — the stall is structural.

## Where the boundary error actually is

Worth knowing before anyone tries to improve it: the Hausdorff distance of
0.034 is **one surface**, and it is neither a block-size nor a chordal effect.

| surface | faces | max | p99 |
|---|---|---|---|
| inlet / outlet | 2044 | 0.00000 | 0.00000 |
| periodic_A / _B | 3080 | 0.0017 | 0.0012 |
| bl_interface_shroud | 1945 | 0.0033 | 0.0028 |
| bl_interface_hub | 1943 | 0.0084 | 0.0060 |
| **ogrid_interface** | 1430 | **0.03405** | **0.02310** |

All 200 worst faces lie on `ogrid_interface`. Block count does not move the
number at all (v11 and v11m agree to five decimals — merging changes no
geometry), and refilling 10 % coarser does not either. Boundary vertices are
as far off the surface as face interiors (mean 0.00116 against 0.00135), so
the error is **not chordal**: spline or higher-order block edges would address
a term of about 0.0002.

`ogrid_interface` is not real geometry — it is the cut face towards the blade
O-grid that `reattach.py` glues back on. Measured against that O-grid block's
own surface, the core sits median **0.00214**, p95 0.01566, max 0.02734 away,
at a local cell edge of 0.042. Note that `reattach.interface_gap` reports this
as median 0.0218 because it compares CENTROIDS of differently sized quads; it
overstates the gap by an order of magnitude and should be replaced by a
point-to-triangle measure.

## Next steps, in order

1. **Implement the projection (decision B/G2).** Pull the `ogrid_interface`
   face grids onto the O-grid block's surface *inside* `tfi.refill_block`,
   while it resamples them — the Gordon-Hall fill then absorbs the motion. It
   must be applied per VERTEX at complex level, not per block face, or the
   coordinate weld in `refill_complex` breaks; `check_watertight` catches that
   and refuses to write, so it cannot pass silently. This closes both the
   boundary error and the assembly seam in one step.
2. **Finish the v18 comparison.** The raw run is done and interesting
   (`--full-constraints` halves the raw block count on a 310-edge input,
   242 → 120); the block structure is not measured. `postprocess_bases.sh v18`.
   Needs a full `clean_blocks` run — see the CPU note below.
3. **Plateau moves in the sheet collapse.** The collapse still accepts only
   strict improvement in (excess, blocks) and therefore stops in a local
   optimum every time. Allowing an equal-cost round with a tabu list is the
   cheap test of whether a locally-neutral collapse unlocks a later one;
   `FRAMEFIELD_PLAN.md` §2.2 names the ILP version (Duan 2023) as the real fix.
4. **Wall-normal grading in the core, if it is ever wanted.** `tfi.py` has
   `clustered_fractions` but the complex refill does not use it: a one-sided
   distribution is not invariant under the mirror relating two blocks' views
   of a shared face, so clustering a direction class would tear the seam.
   Today the grading comes from `reattach.py`, which is where the first cell
   height is set, and the core's outer faces are interfaces, not walls.
5. **Dataset generation** — explicitly last, by the user's instruction. See
   `ANALYSIS_PLAN.md` and the TODO at the end of `PROGRESS.md`.

**Compute note.** This box is a 2-vCPU VPS with a fair-use CPU limit on
SUSTAINED load. Twelve hours of two parallel `clean_blocks` runs got it
throttled to ~10 % of its own cores (steal 90 %); ~90 minutes of sequential
AlgoHex runs did not. Run one compute job at a time.

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

**A relative guard is only as good as what it points at.** The sheet-collapse
guard rejects anything that lowers the worst cell *relative to now*, so
improving the mesh first raises the bar and blocks the collapse. That silently
cost v11 six blocks (16 → 22) once the pre-collapse untangle was added, and it
cost them again when a repair was tried between rounds. The bar is now fixed
at the raw AlgoHex quality for the whole collapse. Any future "let us clean
this up first" idea should be checked against that pattern.

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
