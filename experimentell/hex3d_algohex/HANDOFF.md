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

---

## P3a — the canada `-n` sweep (2026-10-30)

Question: does the target cell count `-n` decide whether the hex complex can be
turned into a valid block structure? Four replays of one parametrisation state
(`/work/canada/dt/cad1/state/param_state.ovm`, 352949 input tets) at
n = 2000 / 8000 / 20000 / 100000, each run exactly as

```
/work/tools/algohex-size/build/Build/bin/HexMeshing --load-param-state <state> \
  -o /work/canada/sweep/n<N>/hex_n<N>.ovm \
  -j /work/canada/sweep/n<N>/metrics_n<N>.json -n <N>
```

and ranked by bad valency, then invalid tets. Collector: `canada_sweep.py`.

| n | valid_volume | bad valency | invalid tets | hexes | rank |
|---|---|---|---|---|---|
| 20000 | 0.998956 | 4 | 7157 | 19915 | 1 |
| 100000 | 0.999818 | 6 | 2551 | 105837 | 2 |
| 8000 | 0.998049 | 13 | 8960 | 7888 | 3 |
| 2000 | 0.996418 | 106 | 10930 | 1823 | 4 |

`bad valency` is `Parametrization.n_invalid_valencies`; its seamless variant is
0 at every n, i.e. the parametrization stage repairs every invalid valency it
finds. `invalid tets` is `n_invalid_param_tets`; its seamless variant is 476 /
440 / 459 / 467, i.e. the seamless repair fixes all but a few hundred of them at
every n. The log's `remaining defects: #invalid_tets = 596 and
#invalid_edge_valencies = 0` is the same for all four n — it describes the
input tet mesh restored from the state, not the run. The n = 20000 replay
reproduces the prior reference (`state/metrics_n20k.json`: valid_volume
0.9989558, 4 bad valencies, 7156 invalid tets; `state/hex_n20k.ovm` has 19915
polyhedra against the "~19785 hexes" of the brief).

**The beam stage, run with the flags of `datagen/run_sample_keepbl.sh`**
(`--width 32 --depth 20 --guard struct --hausdorff-tol 0.15 --workers 4`):

| lane | role | result |
|---|---|---|
| n = 20000 | best 1 | **REJECT** — raw hex mesh off the input surface: Hausdorff 1.2234 > 0.15 (164 blocks, 19922 cells) |
| n = 8000 | diagnostic 3rd | **REJECT** — Hausdorff 1.2664 > 0.15 (125 blocks, 7888 cells) |
| n = 100000 | best 2 | gate passed (Hausdorff 0.1005); base 330 blocks / 105837 cells, **10 inverted cells, min scaled Jacobian -0.9906**, 634 singularities, 35 tiny blocks, 86 excess; depth 1 gave 47 admissible collapses, best 210 blocks, in 1095 s; still running when the report was written |

`n = 8000` was run as a third, clearly labelled diagnostic lane because the
n = 100000 beam cannot finish (below); it is the only other lane that produces
a recorded beam row.

**Which n removes inverted cells? None.** The inverted cells are not a
parametrization defect and raising n does not remove them — it makes them
worse. After `fill_cavities` + `untangle` the raw mesh still holds

| n | inverted cells before untangle | after untangle | min scaled Jacobian |
|---|---|---|---|
| 8000 | 52 | 6 | -0.3983 |
| 20000 | 61 | 5 | -0.1498 |
| 100000 | 148 | 10 | -0.9906 |

`untangle` moves vertices and reduces the count (61 → 5 at n = 20000, 148 → 10
at n = 100000) but always leaves a residue, and the residue is what the
collapse guard measures. The log lines are reproducible: the prior canada runs
(`/work/canada/dt/cad1/beam_n20k.log`, `beam_n8k.log`,
`state/beam_n100k.log`) show the same 5 inverted / min -0.1498, the same
Hausdorff 1.2234 / 1.2662 and the same 10 inverted / min -0.9906.

**Does beam reach a small valid structure? No — two different blockers.**

1. **n ≤ 20000: the raw hex mesh does not hug the input surface.** The gate
   `max_base_hausdorff = 0.15` fires before the search starts, so there is no
   collapse, no `blocks.vtk`, no TFI and no validator run. The blocker is a
   *surface distance*, not a cell count: max Hausdorff 1.2234 (n = 20000) and
   1.2664 (n = 8000) against 0.15, i.e. eight times the gate. Nothing in the
   `-n` direction fixes it, because n = 20000 is already the best of the four
   on valency and invalid tets and it still fails by 1.07.
2. **n = 100000: the gate passes and the guard cannot.** The search starts
   (330 → 210 blocks at depth 1) but every admissible candidate keeps
   min_sj ≈ -0.99, because the guard `min sJ >= 0.05` is measured against the
   same 10 inverted cells (scaled Jacobian -0.9906, a fully flipped hexahedron).
   The prior completed run of this very mesh (`state/beam_n100k.log`,
   `guard=abs`, width 6, depth 6) ended `[beam] depth 1: no admissible
   collapse, stop` and HexBlockValidator reported **INVALID** with one hard
   failure, `10 inverted cells (scaled Jacobian <= 0)`: 326 blocks, cuboid_fraction
   0.871, 35 tiny blocks, block face counts {6: 284, 7: 21, 8: 9, 9: 2,
   10: 10}, min scaled Jacobian -0.9906, boundary Hausdorff to input surface 0.1005.
   So at the one n where the beam
   gets to run, the answer is still a structure with 10 inverted cells.
   It also cannot be finished here: depth 1 took 1095 s on 105837 cells and
   depth 2 must repeat that for each of the 32 beam entries (~2.4 h per level,
   20 levels).

**Net answer for the pipeline:** the canada domain is blocked before the block
structure, at two places that are independent of `-n`. At n ≤ 20000 the raw
hex mesh is off the input surface by ~1.2 (gate 0.15); at n = 100000 it is on
the surface (0.10) but 10 inverted cells (scaled Jacobian -0.9906) keep the collapse
guard at min_sj >= 0.05 unreachable. Both are HexEx/extraction-stage defects,
not parametrization ones — which is consistent with the metrics: every invalid
*valency* is repaired (seamless 0) while inverted *cells* survive and grow
with n.

**Reproduce / re-check:**

```
cd experimentell/hex3d_algohex
/opt/venv-sci/bin/python canada_sweep.py run /work/canada/sweep/report.json \
    --no-replay --beam-ns 20000 8000 100000   # collect, rank, beam the lanes
/opt/venv-sci/bin/python canada_sweep.py --check /work/canada/sweep/report.json
```

`canada_sweep.py` never deletes or overwrites a canada run: a replay is refused
while `hex_n<N>.ovm` exists (unless `--force`), a beam lane is skipped when its
log already holds a result, and logs are appended to. Two artefacts of the
session are left in place rather than cleaned up:
`/work/canada/sweep/n32000` is a stray lane from a first command whose loop
computed `2000*4**i` (32000 instead of 20000 — it completed and is kept), and
`n2000/hexmeshing_n2000.log` holds two attempts after an accidental relaunch
(valid_volume 0.99479 → 0.99642; the ADMM is not bit-reproducible). The sweep
process that launched the n=100000 beam lane is still running: when it ends it
rewrites `report.json` from its own two-lane selection (n=20000 + n=100000),
which also satisfies the `--check` contract, and the n=8000 diagnostic lane
stays on disk in `beam_n8000/`.

## P3b — the 10 inverted cells of canada n = 100000 (2026-10-30)

**Answer: 9 of the 10 are a coordinate defect and are repaired locally; the
10th is topologically locked and cannot be fixed by moving vertices.** The mesh
goes from 10 inverted cells (min sJ −0.9906) to **1** (min sJ −0.0438), with the
Hausdorff distance to the input surface unchanged at 0.10053 (gate 0.15).

**Tool.** `canada_repair_probe.py` (standalone; `clean_blocks.py` defaults are
untouched). `run` replays the beam's own preprocessing, localizes the cells,
repairs only their neighbourhood and writes `/work/canada/sweep/repair_probe.json`
plus `repair_probe_stages.npz` (both mesh stages); `--check` re-derives every
count, minimum and centroid from the saved geometry instead of trusting the
report. The running width-32 lane `beam_n100000` was not touched.

**Where they are.** The 10 cells form five clusters, four of them single cells:

| cell | centroid | surface label | singular edges | volume | min sJ | max corner dist to surface |
|---|---|---|---|---|---|---|
| 85   | (5.429, 1.318, 6.380)   | 5  | 1 | 0.00182 | -0.4045 | 0.162 |
| 90   | (2.307, 0.963, −3.944)  | 10 | 2 | 0.00232 | -0.1262 | 0.115 |
| 102  | (5.257, 1.109, 9.941) | 12 | 1 | 0.00105 | -0.2334 | 0.124 |
| 117  | (0.589, 0.512, 0.242) | 2  | 1 | 0.00165 | -0.1624 | 0.094 |
| 276  | (2.208, 0.908, 2.442)   | 3  | 0 | 0.00160 | **−0.9906** | 0.203 |
| 852 | (2.209, 0.889, 2.538) | 3 | 1 | 0.00078 | -0.1737 | 0.285 |
| 2054 | (2.204, 0.860, 2.643) | — | 1 | 0.00041 | -0.1265 | 0.285 |
| 15707 | (1.563, 0.727, 0.238) | 2  | 0 | 0.00083 | -0.1191 | 0.091 |
| 17799 | (1.506, 0.666, 0.243) | 2  | 0 | 0.00050 | −0.0846 | 0.090 |
| 17800 | (1.618, 0.660, 0.233) | 2  | 0 | 0.00082 | -0.1584 | 0.093 |

Four lone cells, each on a *different* boundary surface (labels 5, 10, 12, 2)
with 1–2 singular edges; a three-cell chain at (2.21, 0.89, 2.44–2.64) in a
valence-6 singularity cluster, which holds the fully flipped −0.9906 cell; and
three cells stacked on the z ≈ 0.2 plane (label 2) with no singular edge at all.
No cell touches two surface labels, i.e. none sits on a feature curve — the
labelling/valency machinery is not what fails here, the extraction is.

**What fixed them.** A local one-sided barrier `sum max(0, eps − J)²` on the
scaled Jacobians of the inverted cells and their neighbourhood — the objective
of `clean_blocks._smooth`, with the barrier grid and the ring widening driven
from the probe — accepted only when the *global* inverted count drops. Two
sweeps did it: (rings 2, eps 0.15) 10 → 2, then (rings 2, eps 0.05) 2 → 1
(~30 s of local work). Every wider or differently weighted sweep of the grid
(rings 2/3, eps 0.02–0.15) is worse and keeps 2–10 inverted. Nothing outside
the neighbourhoods moved: cell count unchanged (105837), Hausdorff unchanged.

**The cell that cannot be fixed.** `852` at (2.209, 0.889, 2.538):
volume 0.00126, edges 0.0519–0.4409 (aspect 8.5), all eight vertices interior,
one singular edge. It is one half of a two-cell doublet: at the start of the
repair the shared quad was non-planar by 0.046, i.e. 31% of a cell edge. The
other half (cell `2054` at (2.204, 0.860, 2.643)) is repaired — it
now sits at sJ 0.004 — so the pair was genuinely folded, not merely displaced,
and the surviving half is the one that cannot be straightened. Measured floor: a local
barrier with 2, 3, 4 and 6 rings of free vertices (12 → 965 vertices, 32 →
1015 cells), eps 0.02–0.15, lifts its own min sJ from -0.1737 to at best
−0.0368 and never clears 1 inverted; every setting that pushes a single cell
further inverts neighbours instead (up to 10). All eight vertices are interior,
so the input surface is not the constraint. The corner tetrahedra of this cell
are ordered opposite to its neighbours; making them positive means pushing
those neighbours through zero volume, i.e. the defect is a parametrization
fold, not a coordinate error. A topological op (collapsing the sheets through
the doublet, which is what the beam does elsewhere) can remove it; vertex
motion cannot.

**Bounded beam.** Not run: step 3 is conditional on 0 inverted, and with the
worst cell at −0.0438 no collapse can meet `--guard abs` (`min sJ >= 0.05`).
The report records `beam: SKIPPED` with that reason. Note the earlier completed
width-6/depth-6 lane on the unrepaired mesh also stopped at
`depth 1: no admissible collapse` — this is the same guard, and the gap it has
to close is now 0.09 in sJ (from −0.0438 to 0.05) instead of 1.04, but it is
still a wall for vertex motion alone.

**Consequence for the sweep.** The `-n` knob is exhausted: n = 100000 passes
the surface gate and is now *one* folded hex away from a clean mesh, but that
hex is not a coordinate defect. The next lever is upstream — the HexEx param
state (`cad1/state/param_state.ovm`) or a topological removal of the
doublet — not another replay at a different n.

**Reproduce / re-check:**

```
cd experimentell/hex3d_algohex
/opt/venv-sci/bin/python canada_repair_probe.py run        # ~27 min
/opt/venv-sci/bin/python canada_repair_probe.py --check /work/canada/sweep/repair_probe.json
```

`repair_probe.json` holds every cell (centroid, bbox, edge lengths, volume,
min sJ, surface distance, labels, singular edges), every repair sweep with its
inverted count, the three wide probes that establish the floor, and the beam
record; `repair_probe_stages.npz` holds both mesh stages so `--check` verifies
the counts and centroids against geometry.

---

## HexMe calibration sweep (P4b)

`CALIBRATION.md` is the hexme calibration matrix (6 shapes × 3 feature
settings at n=8000): one table per geometry, the status vocabulary
(works/diverged/segfault/timeout/running), the AlgoHex failure modes verbatim
and a synthesis. It is built and checked from `/work/calib/report_hexme.json`
by `calib_report.py`:

```
/opt/venv-sci/bin/python calib_report.py --build /work/calib/hexme_sweep \
    --out /work/calib/report_hexme.json
/opt/venv-sci/bin/python calib_report.py --check /work/calib/report_hexme.json --rows 18
```

What it settles for the next sweep: (1) `-d` alone is a no-op on tagged inputs
— A and B are the same mode; only `--force-feature-threshold` (C) changes the
feature network; (2) the prescribed cell colors collapse the torus
(valid_volume 0.699, 0 hexes) while C meshes it, so C is the default for
third-party/hexme meshes; (3) `-n` is a target, not a budget, in C mode
(+41 % on n04u at -n 8000); (4) `--keep-prisms` is a no-op on every hexme shape
(no prism cells, no tet_prep step in the path) and only matters in this keep-BL
chain. `n=8000` is the working point (the n=2000 smoke was too coarse).

---

## Charité LV sweep (12 runs, n=2000)

See the "Charité LV sweep (n=2000)" section of `CALIBRATION.md` for the
per-case tables (Case09_pre/post, Case199_pre/post): 4 cases × 3 feature
settings, 6 works / 6 diverged, all failures parametrization-side
(exit 241/245/250). Built and checked from `/work/calib/report_charite.json`
by `calib_report.py` (`--check --rows 12`).

---

## T19 — the inverted cells of canada n = 100000 as VTK (2026-10-30)

`canada_export_vtk.py` writes the two saved stages of the n = 100000 run
(`/work/canada/sweep/repair_probe_stages.npz`, `/work/canada/sweep/repair_probe.json`)
out as legacy VTK into `/work/canada/sweep/inspect/`, so the defect that
`canada_repair_probe.py` could only report as numbers becomes something you
can look at. sJ here is the minimum corner scaled jacobian
(`clean_blocks.scaled_jacobians`; `ovm_io.scaled_jacobian` is the same thing
for one hex): a hex can keep a positive signed volume and still have inverted
corners, so this — not the volume — is the validity test.

| file | what it is | cells | `inverted == 1` |
| --- | --- | ---: | ---: |
| `hex_n100000.vtk` | the pre-repair mesh — what the beam measures, i.e. the `fill_cavities` + `untangle` stage of `hex_n100000.ovm` | 105837 | 10 |
| `inverted_cells.vtk` | those 10 hexes plus their 1-ring neighbourhood | 140 | 10 |
| `repaired_n100000.vtk` | the post-repair mesh | 105837 | 1 |

All three are ASCII legacy `UNSTRUCTURED_GRID` files (VTK 2.0 header,
`VTK_HEXAHEDRON`, one `SCALARS` block per array). The two full meshes keep the
mesh's own point array and cell order, so the cell id in ParaView's spreadsheet
view is the cell id in `repair_probe.json`; `inverted_cells.vtk` keeps the same
point array and the same cell ids and only drops the 105697 cells that
are not in the 1-ring, so the defect is findable without the 105837 cells.
`--check` re-measures the sJ of both stages from the npz, re-reads the three
files, and confirms that the stage's connectivity is the `.ovm`'s — the `.ovm`
and the stage are the same 105837 hexes, the same vertex sets; only the points
moved (first `untangle`, then the repair).

Cell data:

* `scaled_jacobian` (float) — the cell's sJ, negative exactly where the cell is inverted;
* `inverted` (int) — 1 for an inverted cell, 0 else;
* `distance_to_inverted` (int, `inverted_cells.vtk` only) — 0 for an inverted
  cell, 1 for a cell that shares a vertex with one.

To see it in ParaView:

1. Open the file (`File → Open`; an ASCII legacy file is picked up by the
   *legacy VTK reader* from its header).
2. Add a `Threshold` filter (`Filters → Common → Threshold`), set `Scalars` to
   `inverted` with `Minimum = Maximum = 1`, Apply. What is left are the
   10 inverted cells — in `inverted_cells.vtk` use
   `distance_to_inverted` from 0 to 1 to get the 140-cell patch instead
   (130 neighbours, which is the part a repair has to move).
3. Color by `scaled_jacobian` (diverging map, no "Rescale to Data Range") to
   see which of the 10 are collapsed (sJ near -1) and which are
   only slightly folded.
4. Click a cell and read `Cell Id` in the `Information` tab (or open the
   *Spreadsheet View*) to get the id to feed back into the probe.

The 10 inverted cells of the pre-repair stage (centroid in mesh
coordinates, sJ = minimum corner scaled jacobian; negative means inverted):

| cell | x | y | z | sJ |
| --- | ---: | ---: | ---: | ---: |
| 85 | 5.428614 | 1.317982 | 6.380345 | -0.404510 |
| 90 | 2.307035 | 0.963002 | -3.944410 | -0.126242 |
| 102 | 5.257099 | 1.109242 | 9.941260 | -0.233411 |
| 117 | 0.588783 | 0.511962 | 0.242499 | -0.162434 |
| 276 | 2.207800 | 0.907584 | 2.442037 | -0.990573 |
| 852 | 2.209148 | 0.888689 | 2.537983 | -0.173731 |
| 2054 | 2.204447 | 0.859904 | 2.643421 | -0.126488 |
| 15707 | 1.563037 | 0.727327 | 0.237509 | -0.119079 |
| 17799 | 1.505772 | 0.666320 | 0.243490 | -0.084556 |
| 17800 | 1.617603 | 0.660174 | 0.232996 | -0.158389 |

After the repair 9 of the
10 are clean and one remains: cell 852 (centroid
(2.208712, 0.894427, 2.543348),
sJ -0.043762) — the folded cell of a two-cell doublet in a
valence-6 singularity cluster, which the beam reports as topologically locked
(`repair.remaining[0].reason`). That is the cell to look at in
`repaired_n100000.vtk` with the same threshold.

Reproduce / verify:

```
# (re)write the three files (~3 s: reads the npz, no AlgoHex, no beam)
/opt/venv-sci/bin/python experimentell/hex3d_algohex/canada_export_vtk.py run
# verify them (~30 s): the three files exist, cell counts match the stages and
# the report, the sJ arrays are the stages' sJ, 10 / 1 cells are inverted, the
# flagged ids match repair_probe.json, the 1-ring is complete, and the stage's
# connectivity is the .ovm's (--no-ovm skips that last one)
/opt/venv-sci/bin/python experimentell/hex3d_algohex/canada_export_vtk.py --check /work/canada/sweep/inspect
```

---

## T33 — keep-BL tistos: feature curves as VTK (2026-10-31)

`export_features_vtk.py` reads the keep-BL tistos tet input
`/work/datagen_keepbl/runs/tistos_s1_0131/tet.vtk` (legacy ASCII VTK,
`CELL_DATA 'color'` with nonzero = feature element), counts the feature graph
and exports the feature **curves** as explicit VTK_LINE cells to
`/work/datagen_keepbl/inspect/keepbl_features.vtk` (all source points copied,
so the lines keep the source vertex indices; feature points as vertex cells;
per-line `CELL_DATA 'feature_curve_id'`, -1 on vertices).

Counts, source → export (must be identical):

| | feature edges | curves | curve sizes / topology |
|---|---|---|---|
| keep-BL source | 580 | 3 | 462 open, 72 closed, 46 closed |
| keep-BL export | 580 | 3 | same ids 1/2/3 |

`--check` re-parses the export and asserts **exact edge-set equality** plus the
same `feature_curve_id` partition (verified to also fail on a tampered export,
exit 1), exits 0 only if equal. Feature points: 8, feature tris: 26728,
surface tags on tris [1, 2, 3, 4, 7, 8, 9]. For context, the stripped
reference `/work/datagen/runs/tistos_s1_0131/tet.vtk`: 534 edges, 3 curves
(310 open, 112 + 112 closed), tags [1..7]; the difference is the extra
keep-BL surfaces (tags 8, 9; the mesh header says T1_9) — not an export
artefact.

```
# export + report (~4 s)
/opt/venv-sci/bin/python experimentell/hex3d_algohex/export_features_vtk.py
# verify export vs source (exact)
/opt/venv-sci/bin/python experimentell/hex3d_algohex/export_features_vtk.py --check /work/datagen_keepbl/inspect/keepbl_features.vtk
```

---

## T40 — canada n = 100000, struct-guard beam + TFI attempt (2026-11-02)

**Answer: the struct guard buys a *structurally* valid subdivision, but TFI
cannot refill it, so the lane stops here.** The depth-1 beam result is a real
block complex (no block touches two opposite or two adjacent faces, the face
count is 6 for 180 of 210 blocks), the hexahedral validator still calls it
INVALID because 9 cells are inverted, and the TFI refill **fails** — twice over.

**Tool.** `canada_struct_beam.py` (standalone, in this directory).
`run --depth 1` runs the beam, the TFI attempt and the export, and writes
`/work/canada/sweep/struct_beam/report.json`;
`--check REPORT.json` re-derives every number from the logs and the files on
disk and exits 0 only if they agree. The lane was started by T20; T40 finished
it (the report carries `"continued_by": "T40"`).

### The beam (depth 1)

Guard `struct`: inverted cells allowed (≤ 148), min sJ ≥ 0.05, Hausdorff ≤ 0.15,
width 4, depth 1. The base complex is 330 blocks / 105 837 cells with 10
inverted cells; after the untangle the interior-only pass fails and 31 boundary
vertices slide on the input surface.

| | cells | blocks | cuboids | excess | inverted | min sJ | Hausdorff | wall |
|---|---|---|---|---|---|---|---|---|
| base | 105 837 | 330 | 287 | 86 | 10 | −0.9906 | 0.10053 | — |
| depth 1 (best of 4 states) | 103 824 | **210** | 180 | 66 | **9** | −0.9906 | 0.10053 | 1174 s |

47 of the states at depth 1 are admissible. The four states are
210 bl/ex 66/sJ −0.99, 246 bl/ex 69, 186 bl/ex 74, 275 bl/ex 78 — so **210
blocks is the smallest** the guard allows at depth 1, and the excess (66
non-cuboid block faces) is what the guard buys.

**Is the structure a valid subdivision?** Yes, structurally: the validator's
only failure is `9 inverted cells (scaled Jacobian <= 0)`, and the report records
`beam.structurally_valid = true` (the inverted cells are a cell-quality
failure, not a structural one). Face census:
180 blocks with 6 faces, 20 with 7, 4 with 12, 4 with 10, 1 with 11, 1 with 5;
11 tiny blocks, min 2 cells per block; mean sJ 0.9489. Blocks touching one
surface on two opposite faces: 0; on two adjacent faces: 0.

### The TFI attempt — two failures

`TFI cannot build the lattices`, and it fails **before** the `--require-lattices`
gate, so the run never reaches it:

```
[tfi] 9 blocks are not lattices and will be passed through unchanged:
      [8, 15, 25, 26, 83, 86, 88, 92, 96]
[tfi] conformity of the existing divisions: 1 classes disagree
AssertionError: class 0 below its lower bound      (tfi.py:492, solve_block_divisions)
```

Exit code 1, no `tfi_blocks.vtk` written. Reproduced with `--solve-divisions`
alone (without `--apply-divisions`), so it is the **division solve** that dies,
not the refill. The root cause is visible in the line above the crash: the
beam's own divisions are **not conforming** — one direction class disagrees —
so `solve_block_divisions` cannot find counts that satisfy every block's lower
bound. The 9 blocks that are not lattices are a second, independent problem.

The plain TFI refill (no `--solve-divisions`, at the existing divisions)
answers the other half of the question — for the 201 blocks that *do* carry a
lattice, does the refill hold up? **No, it makes the mesh worse:**

| metric | before | after |
|---|---|---|
| scaled Jacobian (min) | −0.9906 (9 bad) | −0.9906 (**102 bad**) |
| non-orthogonality (deg) | 81.14 (30 bad) | 81.78 (**138 bad**) |
| skewness | 0.869 (0 bad) | 1.199 (0 bad) |
| aspect ratio | 20.34 (0 bad) | 13.89 (0 bad) |
| vol ratio | 0.118 (0 bad) | 0.003 (**5 bad**) |
| face weight | 0.0105 (2 bad) | 0.0005 (**86 bad**) |
| face flatness | 0.0000 (257 bad) | 0.0000 (345 bad) |

So the honest summary is: **the struct guard preserves the block topology but
destroys the conformality that TFI needs.** The inverted cells are the visible
symptom; the non-conforming divisions are the cause.

### What this means for the pipeline

- The **export is skipped**: with no lattices there is nothing for
  `export_sample.py --check` to round-trip, so the report records
  `export.status = "skipped"`, reason `TFI cannot build the lattices`. No
  `sample_100000.npz` is produced from this lane.
- Depth 2/3 are out of reach for the same reason and for the clock: the
  measured per-state cost is 1174 s, so depth 2 ≈ 5870 s (fits the timeout)
  and depth 3 ≈ 10 566 s (does not). Deeper is not the blocker — **conformity
  is**.
- The fix is on the beam side, not the TFI side: the collapse has to leave a
  **conforming** division per direction class (every block of a class sharing
  the same counts along its class axis), and the 9 blocks whose shape is not a
  lattice have to be reshaped before the refill, not passed through.

### Files and commands

```
/work/canada/sweep/struct_beam/
  beam_d1.log        the depth-1 beam run (the reused log)
  beam_d1.pid        its pid (the process has exited)
  blocks.vtk         the 210-block result (103 824 cells)
  blocks_quality.vtk per-cell scaled Jacobian, for the inverted cells
  blocks_edges.vtk   the block edge graph
  tfi.log            the --solve-divisions run: the AssertionError
  tfi_rebuild.log    the plain refill at the existing divisions: the metric table
  report.json        the report
```

```
# run the lane (reuses beam_d1.log, tfi.log and tfi_rebuild.log if present)
/opt/venv-sci/bin/python experimentell/hex3d_algohex/canada_struct_beam.py run --depth 1
# verify the report against the logs and the files on disk
/opt/venv-sci/bin/python experimentell/hex3d_algohex/canada_struct_beam.py --check /work/canada/sweep/struct_beam/report.json
```

---

## Documents

| file | content |
|---|---|
| `README.md` | the working reference: pipeline, commands, modules, what was learned |
| `DATASET_PIPELINE.md` | the dataset plan, task by task, with a resume block |
| `DATA_GENERATION.md` | how the geometries are generated (dtOO), first candidates |
| `RUNS.md` | every T1_9 AlgoHex run with input, runtime, outcome |
| `MESH_QUALITY.md` | the CFD criteria, measured on our own mesh |
| `CALIBRATION.md` | the hexme P4b calibration sweep: 6 shapes × 3 feature settings at n=8000, per-geometry status tables + synthesis; and the Charité LV sweep: 4 cases × 3 settings at n=2000 |
| `PROGRESS.md`, `SESSION_*.md` | chronological logs including every failed attempt |
| `../../docs/decisions/` | every decision with its measurements and rejected options |
| `../../docs/cluster-enroot-findings.md` | every cluster-side trap |

## T43/T48 - trailing-edge ring raggedness

**Source:** `/work/canada/sweep/inspect/hex_n100000.vtk` (105837 hexes, 116806 points). Axis x; trailing-edge region x > 9.0 (first 0.5-bin with z-span < 25% of inlet 13.900095, te_zspan_at_te0=3.1927825). Tip vertex (11.0,0.0,10.0). Rings extracted as outer-skin vertex chains in 36 polar bins around local channel centre, closed per cross-section (`rings.json` meta).

**PNGs:** `/work/canada/sweep/inspect/trailing_edge_rings.png`, `/work/canada/sweep/inspect/trailing_edge_jacobian.png`, `/work/canada/sweep/inspect/trailing_edge_rings_control.png`.

**TE rings (x > 9.0), n=25:** per-ring `station_x`, vertex counts 6-54, full 360° coverage. Geometric raggedness (wobble/spikes):
- `dev_axial_max` up to 0.2925 (median across TE rings ~0.088), `dev_axial_median` up to ~0.070; `dev_spline_max` up to 0.3211 (median ~0.162), `dev_spline_median` up to ~0.101.
- `kink_max_deg` 180° on most rings (pinched/pointing joints), `kink_p90_deg` ranges 85-136° (median ~92°).
- Ring spacing: `gap_mean` 0.004-0.327 (median ~0.027), `gap_jump`/`gap_station_jump` show large jumps (e.g. 0.154, 0.375) between consecutive rings.
- Normalised by spacing, rings become very "stretched": `dev_axial_over_gap` reaches 175x in extreme TE cases (median across TE rings large, often >>1), indicating the axial wobble is comparable to/larger than ring spacing in many places.

**Mid-channel control rings (x~1.56-3.03), n=380:** similar geometric wobble relative form but finer spacing.
- `dev_axial_max` up to ~0.387, `dev_spline_max` up to ~0.241; `kink_p90_deg` can be very small (down to ~0.8-1.4° in tight clusters) with `kink_max_deg` still 180° in places.
- `gap_mean` often 0.001-0.206, `gap_jump` small in tight packs but `dev_axial_over_gap` explodes to 10-136x in the tightest spacing regions (e.g. stations ~1.566 give dev_axial_over_gap 50-175x). So normalised raggedness is also severe in tight control regions.

**Cell quality / relation to inversion:** no inverted cells in TE region (x>9.0). `sJ_min` per ring >= 0.8459 in TE (mesh mean sJ_high), `n_inverted` 0 for all listed rings. Correlation between ring wobble and worst sJ in ring is negative (-0.65 over TE) - more ragged rings touch slightly lower-sJ cells but never inverted; raggedness is geometric (post-processed vertex columns) not due to bad/inverted cells.

**Numbers (representative):** TE ring with largest axial deviation `station_x`~9.56635 has `dev_axial_max` 0.2925, `kink_p90_deg` 85.4°, `gap_mean` 0.0432, `dev_axial_over_gap` 6.77. A tight TE ring `station_x`~1.565875 has `dev_axial_max` 0.3275, `dev_spline_max` 0.175, `kink_p90_deg` 88.25°, `gap_mean` 0.0024, `dev_axial_over_gap` 136.2. The 10 inverted cells of the mesh sit upstream (centroids down to z~-3.94, x>=0.59) far from x>9.0; raggedness is orthogonal to the inversion issue.
