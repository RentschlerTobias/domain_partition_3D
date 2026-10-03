# Handoff — hex3d_algohex, state and how to continue

Written for someone (or some session) starting with no context. Read this
first, then `README.md` (the working reference: pipeline, commands, modules,
what was learned). `PROGRESS.md` is the chronological log — long, and only
worth opening for a specific question.

---

## What this does

Generates the minimal hexahedral block structure of a turbomachinery flow
passage — from a mesh of our design system, design tool Object-Oriented
(dtOO) — via a genuine 3D frame field, and exports it as one `sample.npz` per
machine — the training data of the block-structure transformer in `meshtron`.

```
dtOO mesh.msh ─ tet_prep_v5 ─▶ AlgoHex ─▶ fine hex mesh ─ base_complex ─▶ raw partition
  ─ sheet collapse ─▶ minimal blocks ─ block_edges / tfi ─▶ refill ─ export_sample ─▶ sample.npz
                      (reattach / ogrid_extrude: blade O-grid + boundary layer back on, full domain)
```

## Where it stands (2026-10-03)

- **Dataset generation runs.** Machines are Sobol samples of the 30 dtOO
  design parameters (`scripts/generate_machines.py`); one production job
  (`scripts/batch_generate.slurm`) exports their meshes with dtOO and samples
  them on bwUniCluster with enroot, refilling its queue until the walltime
  drain. The batch dataset holds **739 samples**.
- **The labels are noisy, and the cause is known.** The singularity graph and
  the raw partition are identical across runs of one geometry; the greedy
  sheet collapse is not — machine_0004 ended at 12, 22 or 75 blocks. Over the
  dataset: 29 combinatorial / 113 labelled topologies over 739 samples, and
  the design parameters carry no signal about the topology.
- **The fix is measured.** `beam_collapse.py` (beam search over collapse
  orders, `struct` guard) ends 17 of 17 runs at 12 blocks / 12 cuboids, 0
  inverted, VALID, one shared topology —
  [`2026-09-28-beam-collapse-relabelling.md`](../../docs/decisions/2026-09-28-beam-collapse-relabelling.md).
  It is **not yet wired** into `scripts/sample_one.sh`.
- **Blade O-grid.** `ogrid_extrude.py` extrudes a *conforming* blade O-grid
  from the core's block faces (26 blocks on the smoke sample, 0 inverted after
  a joint TFI); `reattach.py` still reuses the dtOO O-grid verbatim.
- **Rejected on measurement:** the block-level cut-set search as a collapse
  substitute (25 blocks at best —
  [decision](../../docs/decisions/2026-09-30-cutset-search-not-a-collapse-substitute.md));
  leading/trailing-edge feature curves as extra AlgoHex constraints
  ([decision](../../docs/decisions/2026-10-03-le-te-feature-edges-rejected.md)).

## Next steps, in order

1. **Relabel with the beam search.** Make the batch scripts call
   `beam_collapse.patch` instead of `--collapse-rounds 5`, and re-run the
   collapse for all 739 samples from the saved `blocks.vtk` + `sample.npz`
   (no AlgoHex rerun needed), written beside the original as
   `sample_beam.npz`, never overwriting. Estimate ~18–20 h on 6–8 cores. The
   open question it answers: does the whole Sobol set land on one canonical
   topology?
2. **Track the quality cost of starting from the greedy intermediate.** On one
   run the min scaled Jacobian was 0.03 from `blocks.vtk` against 0.33 from
   the raw AlgoHex mesh; block corners and edges are unaffected in principle.
3. **Put the blade back for the generated meshes too.** meshtron's inference
   chain refills the core passage only; `ogrid_extrude.py` / `reattach.py`
   live here and are not called from it.
4. **The watertight check** (`tfi.check_watertight`) flags faces used by
   three or more cells and returns the boundary-face count, but both callers
   discard that count; a seam that comes apart shows up only there, as extra
   boundary faces. Compare it against the boundary before the refill (see the
   failed weld in the T1_9 section below).
5. **Block-level T-junctions** (17 % of the corpus, measured in meshtron):
   two blocks touching across only part of a side, which the four-corner face
   format cannot express. Check whether the beam-relabelled structures still
   carry them.
6. Older items, still open: leave the O-grid projection
   (`tfi.project_ogrid_interface`) off; wall-normal grading in the core is
   not used (the grading comes from `reattach.py`).

## Reproduce

One sample, end to end (the unit the batch jobs run):

```bash
experimentell/hex3d_algohex/scripts/sample_one.sh <tet.vtk> <n> <outdir> [params.json]
```

The beam collapse on top of it:

```bash
python experimentell/hex3d_algohex/beam_collapse.py <outdir>/hex_hex_<name>.ovm \
    --input-vtk <tet.vtk> --out <outdir>/blocks_beam.vtk \
    --width 16 --depth 6 --guard struct
```

The T1_9 v11 reference (paths as on the original VPS):

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

## T1_9 reference state (measured up to 2026-09-10)

This is the single-geometry state the dataset route grew out of; its
measurements still hold. **What the structure is for** was settled on 2026-09-08: it is training data
for a transformer, filled by TFI in post-processing and then run as a CFD
mesh, which makes **block count the primary criterion**. Reasoning and
measurements: `docs/decisions/2026-09-08-hex3d-block-structure-objective.md`.

**The reference is v11, unmerged.** Merging the blocks was tried down to 6 and
measured: it is free in every static metric and it wrecks the refill.

| structure | blocks | rule | cells at h=0.05 | inverted | boundary p95 | boundary max |
|---|---|---|---|---|---|---|
| **v11** | **16** | max | 49 550 | **0** | **0.0060** | **0.0343** |
| v11m | 12 | max | 76 025 | 3 | 0.1845 | 0.3710 |
| v11m | 12 | median | 47 256 | 0 | 0.1739 | 0.3743 |
| v11m6 | 6 | max | 249 100 | 2 | 0.2162 | 0.4466 |

> **That table is real but its explanation was wrong twice, and the second
> correction arrived last.** Read this before acting on it.
>
> The boundary numbers for the merged structures are NOT geometric error. They
> come from a **failed weld**. Counted on the refilled meshes:
>
> | mesh | cells | boundary faces | boundary vertices |
> |---|---|---|---|
> | v11 refill | 49 550 | 9 864 | 9 864 |
> | v11m refill | 47 256 | **12 482** | 11 888 |
>
> v11m has 2 618 MORE boundary faces than v11 while having fewer cells. A
> coarser mesh has less boundary, not more, so those faces are interior faces
> whose two sides failed to weld. Their vertices sit inside the domain, which
> is why they measure a tenth of a radius from the surface. v11's worst 200
> boundary vertices all lie on `ogrid_interface`, the one genuinely bad
> surface; v11m's spread over three unrelated surfaces, which is what seams
> look like.
>
> Measured directly, the merged structure's boundary faces are **fine**:
> resampling them at the counts the MILP actually solves gives
>
> | structure | linear p95 | linear max |
> |---|---|---|
> | v11, 16 blocks | 0.00625 | 0.03343 |
> | v11m, 12 blocks | **0.00493** | 0.03343 |
>
> better than v11. So merging does not wreck the refill; a bug in it does.
> `tfi.check_watertight` did not catch this because it tests for faces used by
> THREE or more cells — a seam that comes apart produces extra faces used by
> ONE, and the check has no boundary-count balance. That is the first thing to
> fix.
>
> The class-spread analysis below stands as a measurement (a class does carry
> one count, and merging does raise the spread to 16.7x) but it has NOT been
> shown to be what produced these numbers.

A **direction class carries one cell count for all its axes**. Merging unions
axes into fewer classes, and a class holding axes of 0.06 and 1.01 has no good
count: 1 leaves the long axis 20x too coarse, 20 leaves the short one 16x too
fine. Measured on cand_002:

| structure | classes | worst axis-length spread in a class |
|---|---|---|
| 22 blocks | 12 | **1.3x** |
| 7 blocks | 5 | 16.7x |

and not one of the 34 possible single-pair merges stays under 2x (best 2.1x,
worst 12.5x). `merge_ilp.class_spread` computes this in milliseconds.

`v11m` and `v11m6` stay on disk as the smallest valid block structures. If a
compact topology is what a model should learn, merge for the learning target
and refill the unmerged structure — they are different artifacts, and whether
the merged one can also be refilled is **open**, pending the weld fix.

**The pinch does not harm TFI.** Open since the first handoff, now measured:
v11 refills with 0 inverted cells and an unchanged boundary while carrying it.

**The structured mesh beats the source mesh**, measured with one
implementation (`mesh_quality.mixed_metrics`) over both:

| metric (OpenFOAM limit) | source hybrid, 257 219 cells | structured, 161 619 hexes |
|---|---|---|
| non-orthogonality > 65° | 5105 | **56** |
| skewness > 4 | 4 | **0** |
| skewness p95 | 0.472 | **0.048** |
| face weight < 0.05 | 34 | **0** |
| aspect ratio p95 | **49.1** | 65.8 |

91x fewer non-orthogonality violations with 37 % fewer cells. The aspect ratio
is the one metric that is worse, by design: 17 boundary-layer layers at a
first cell height of 8.9e-4.

| structure | blocks | smallest | pinch | inverted | min scaled Jacobian | Hausdorff | validator |
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
| `README.md` | the working reference: pipeline, commands, modules, what was learned |
| `DATASET_PIPELINE.md` | the dataset plan, task by task, with a resume block |
| `DATA_GENERATION.md` | how the geometries are generated (dtOO), first candidates |
| `RUNS.md` | every T1_9 AlgoHex run with input, runtime, outcome |
| `MESH_QUALITY.md` | the CFD criteria, measured on our own mesh |
| `PROGRESS.md`, `SESSION_*.md` | chronological logs including every failed attempt |
| `../../docs/decisions/` | every decision with its measurements and rejected options |
| `../../docs/cluster-enroot-findings.md` | every cluster-side trap |
