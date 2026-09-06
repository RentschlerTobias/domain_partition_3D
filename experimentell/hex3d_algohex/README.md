# hex3d_algohex — 3D frame-field hexahedral blocks for T1_9

Branch `feat/algohex-3d-frame-field`.

Generates geometry-aligned hexahedral **blocks** for the T1_9 runner passage
via a genuine 3D octahedral frame field, as an alternative to extruding the
2D surface partitions of `dp3d/`. Pipeline:

```
tet mesh + feature tags  ->  AlgoHex (frame field -> IGM -> HexEx)  ->  hex mesh
                         ->  base complex (sheets)  ->  coarse blocks  ->  TFI
```

The motivation for going 3D: `experimentell/3d_extrapolation/hexa_interpolation.py:4-8`
records that the hub and shroud cross fields are *topologically different*
(hub four idx=−1 singularities, shroud four idx=+1 — a real effect of blade
twist). A per-surface 2D field with a ruled lift cannot represent that, which
is why the existing 3D lift needs a morph hack. A 3D field represents it
natively.

**Status**: Stages 0–6 done. Best result: **v11 → 16 blocks, all 6-faced, 0
inverted cells, validator VALID** (`clean_blocks.py --collapse-rounds 5
--untangle` on `T1_9_hex_v11.ovm`). Older text below still describes the v9
path; see PROGRESS.md "Runs v10 and v11".

Stages 0–5 done. A valid hex mesh, a block decomposition and its
postprocessing (`clean_blocks.py`) exist. Two usable endpoints: **81 blocks /
73 cuboids (90 %)** with the mesh nearly untouched, or **42 blocks / 36
cuboids (86 %)** after a sheet collapse and untangling — half the blocks and
**0 inverted cells, `HexBlockValidator` VALID**. TFI is next. See "What was
learned", points 6–10.

---

## Quick start

```bash
PY=/root/repos/duty/quadmesh/.venv/bin/python

# 1. build the AlgoHex input (reduced domain, the best-performing variant)
$PY experimentell/hex3d_algohex/tet_prep_v5.py            # -> data/T1_9/T1_9_tet_v5.vtk

# 2. run AlgoHex  (~45 min, needs the machine to itself, see "Memory" below)
$PY experimentell/hex3d_algohex/run_algohex.py \
      --tag v9 --in-vtk data/T1_9/T1_9_tet_v5.vtk -- -n 60000

# 3. blocks from the result  (raw base complex; see "What was learned" 6 --
#    its own surface classifier gives 48 %, not the 74 % on record)
$PY experimentell/hex3d_algohex/block_faces.py output/hex3d_algohex/T1_9_hex_v9.ovm

# 4. block-structure postprocessing: report, cleanup, validation, export
$PY experimentell/hex3d_algohex/clean_blocks.py

#    optional: mesh-level sheet collapse (Gao), ~15 min per round, halves the
#    block count on v9; changes the hex mesh, so it is off by default
$PY experimentell/hex3d_algohex/clean_blocks.py --collapse-rounds 5 --untangle \
      --out output/hex3d_algohex/deliverable/T1_9_blocks_v9_gao.vtk
$PY experimentell/hex3d_algohex/clean_blocks.py \
      output/hex3d_algohex/cylinder_hex.ovm --no-input-vtk --out /tmp/cyl.vtk

# 5. re-attach the parts cut out of the AlgoHex domain -> full-domain blocks
$PY experimentell/hex3d_algohex/reattach.py

# inspection files for ParaView
$PY experimentell/hex3d_algohex/export_vtk.py
$PY experimentell/hex3d_algohex/plot_stages.py
```

---

## The two usable results

| | **v5** (full domain) | **v9** (reduced domain) |
|---|---|---|
| input | `data/T1_9/T1_9_tet_v2.vtk` | `data/T1_9/T1_9_tet_v5.vtk` |
| tets in | 71 415 | 42 218 |
| hex cells out | 60 612 | 56 661 |
| coverage | 100 % of 6.5253 | 99.8 % of 5.4606 |
| **inverted cells** (scaled Jacobian ≤ 0) | 96 | **2** |
| worst scaled Jacobian | −0.9956 | **−0.0607** |
| mean scaled Jacobian | 0.9651 | 0.9664 |
| non-manifold boundary edges | 0 | 0 |
| runtime | 1 h 17 | 45 min |
| **blocks** | 307 | **82** |
| cuboids | 290 (94 %) | 61 (74 %, 85 % of cells) |
| after `clean_blocks.py` | — | **81 blocks, 73 cuboids (90 %)** |
| after `--collapse-rounds 5 --untangle` | — | **42 blocks, 36 cuboids (86 %), 0 inverted, VALID** |

v9 is the recommended base: 48× fewer inverted cells and a far coarser block
structure, which is what TFI wants. v5 has the higher cuboid share.

v9's domain excludes both boundary layers. They are re-attached by
`reattach.py` (stage 6): full domain 103 157 cells in **73 blocks**.

---

## Run ↔ input-file mapping

Two independent numberings collide here — a run tag `vN` is **not** the same
as an input file `T1_9_tet_vN.vtk`. Explicitly:

| run | input file | domain | feature edges | `-n` | outcome |
|---|---|---|---|---|---|
| v1 | `T1_9_tet.vtk` (classifySurfaces tags) | full | 853 (427 flat) | 10000 | 95.4 % coverage |
| v2 | `T1_9_tet.vtk` (analytic tags) | full | 496 | 10000 | SIGSEGV — dangling feature curves |
| v3 | `T1_9_tet.vtk` (analytic, cleaned) | full | 437 | 10000 | IGM diverged (energy 3.3e48) |
| v4 | `T1_9_tet_v1tags.vtk` | full | 853 | 60000 | 97.9 % coverage, 124 inverted |
| **v5** | `T1_9_tet_v2.vtk` | full, original boundary | 630 | 60000 | **valid IGM, 100 % coverage** |
| v6 | `T1_9_tet_v3.vtk` | minus blade O-grid | 646 | 60000 | aborted (stopped for memory) |
| v7 | `T1_9_tet_v4.vtk` | minus O-grid + BL | 166 | 60000 | OOM in quantization |
| v8 | `T1_9_tet_v4.vtk` | same as v7 | 166 | 15000 | OOM — `-n` is not the cause |
| **v9** | `T1_9_tet_v5_diagbug.vtk` | reduced, exact labels | 550 | 60000 | **completed, 2 inverted cells** |
| **v11** | `T1_9_tet_v5.vtk` | same domain, **labels corrected** | 534 | 60000 | **13 min, 21 inverted, no cavities — best block structure** |
| v10 | `T1_9_tet_v6.vtk` | only the blade layer cut | 854 | 60000 | 1 h 29, 19 inverted, 198 blocks |

Logs and metrics per run: `output/hex3d_algohex/hexmeshing_<tag>.log`,
`T1_9_hex_metrics_<tag>.json`.

---

## Modules

| file | role |
|---|---|
| `tet_prep.py` | Stage 1 base: topological boundary extraction from the hybrid MSH, analytic surface classification, feature graph, AlgoHex VTK writer. Shared helpers used by all later variants. |
| `tet_prep_v2.py` | keeps the **original** boundary triangulation and meshes only the interior (full domain). Input for v5. |
| `tet_prep_v3.py` | reduced domain: blade O-grid removed. Input for v6. |
| `tet_prep_v5.py` | reduced domain with **exact** labels matched against the MSH's own tagged 2D elements. Input for v9. **Current best.** |
| `run_algohex.py` | Docker wrapper for `HexMeshing`, `--tag`/`--in-vtk`, checkpoints the seamless map |
| `ovm_io.py` | OpenVolumeMesh reader, hex VTK/MSH writers, `scaled_jacobian`, `_hex_volume` |
| `base_complex.py` | singular edges, sheet propagation, block partition |
| `block_faces.py` | per-sheet labelling, block faces, cuboid test |
| `clean_blocks.py` | Stage 5: block-structure postprocessing. Exact surface labels by nearest-face lookup; cavity detection and refill; block merge/split on the cut set; optional mesh-level sheet collapse (Gao); `HexBlockValidator`; before/after report |
| `reattach.py` | Stage 6: re-attaches the removed parts — blade O-grid reused verbatim (5 ready-made cuboid blocks), hub/shroud boundary layer regenerated by radial extrusion; block edges per part |
| `showcase.py` | the pipeline in 11 steps, each as VTK (authoritative) + PNG (pyvista) + interactive HTML (plotly) — see `SHOWCASE.md` |
| `export_vtk.py`, `plot_stages.py` | ParaView exports and figures (`plot_stages.py` is superseded by `showcase.py`) |

`tet_prep_v4.py` does not exist — the numbering skips it because
`T1_9_tet_v4.vtk` is produced by an ad-hoc variant of `tet_prep_v3`.

---

## Running AlgoHex

The build lives in a **named Docker volume**, not in the image (the compile
was OOM-killed repeatedly and had to be made resumable):

```bash
docker run --rm --network=host \
  -v $PWD:/work -v algohex-build-cache:/app/build algohex-configured \
  /app/build/Build/bin/HexMeshing -i ... -o ...
```

A self-contained image `algohex:portable` also exists with the binaries baked
in — that is the one to export to a cluster:

```bash
docker save algohex:portable -o algohex.tar
apptainer build algohex.sif docker-archive://algohex.tar     # bwUniCluster
```

**Memory**: this machine has 7.7 GiB. AlgoHex peaks hard during quantization;
two concurrent runs OOM. Run them **sequentially**. `--network=host` is
required — the default bridge network has no DNS here.

---

## A naming note that matters

The surfaces `bl_interface_hub`, `bl_interface_shroud` and `ogrid_interface`
are **interfaces to cells cut out of the domain**, not walls. In the reduced
domain the blade wall and the hub/shroud walls are not present at all — they
belong to the removed O-grid and prism layer, which are re-attached later by
`reattach.py`.

They used to be called `shell_hub`, `shell_shroud` and `shell_blade`, and
that cost three wrong diagnoses in a row: "the block structure at the blade"
repeatedly meant "blocks touching the O-grid cut face", which is a different
thing one cell layer away. Older entries in `PROGRESS.md` still use the old
names; they are a chronological log and were deliberately not rewritten.

---

## What was learned (the non-obvious parts)

**1. Feature constraints: more is better, counter-intuitively.**
Four runs point the same way. Fewer, "cleaner" feature curves consistently
produced a *worse* integer-grid map:

| run | feature edges | surfaces | quantization |
|---|---|---|---|
| v1 | 853 (427 geometrically flat) | 21 | completed |
| v5 | 630 | 7 | 0 invalid tets |
| v3 | 437 (clean) | 7 | IGM diverged |
| v7/v8 | 166 | 3 | OOM |

The extra curves evidently act as regularisation, pinning the octahedral
field into a simpler, hex-meshable singularity graph. Removing them gives a
smoother field that is *not* hex-meshable. The v7/v8 OOM was therefore not a
memory-size problem: the same domain with proper labels (v9, 550 edges)
completes in 45 min.

**2. IGM validity does not predict mesh quality.**
v9's parametrization is the worse one by AlgoHex's own numbers (9030 invalid
tets vs v5's 0, HexEx reporting `Invalid Input`) yet its extracted mesh has
2 inverted cells against v5's 96. Judge the extracted hexes, never
`valid_volume`.

**3. Never classify a surface by a raw coordinate threshold.**
It cuts across the triangulation and produces zigzag feature curves. This
happened twice: a radius split in `tet_prep_v5` (p95 kink 174.8°, 318 of 679
curve nodes above 30°) and again in the v9 block classifier (cuboid share
down to 44 %). Transferring exact labels from the input mesh by nearest-face
lookup fixes both (p95 8.2°; 74 %).

**4. The source MSH's tagged 2D elements are closed but not 2-manifold.**
724 edges carry three faces, because the hex O-grid's *interior* block
interfaces are tagged as 2D elements too. Dropping the 10 interior geom ids
(`10,11,12,16,17,21,22,26,27,31`) leaves a closed 2-manifold that preserves
the original quads — the outer ids are
`1,2,3,4,5,6,7,8,9,13,14,15,18,19,20,23,24,25,28,29,30`.

**5. A hex volume bug that invalidated an early conclusion.**
`_hex_volume` returned exactly 1/3 of the true volume (two of five tets had
swapped vertices). That made a 95 %-complete v1 mesh look like "31.4 %
coverage, full of holes" and drove a whole wrong diagnosis. The same swapped
tets are still present in
`experimentell/3d_extrapolation/hexa_interpolation.py:224-237`
(`hexa_cell_volumes`) — harmless there since it only feeds statistics, but
**unfixed**. The tell was that AlgoHex's own cylinder demo also measured
"32 %".

**6. `block_faces.physical_of` does not reproduce the documented baseline.**
It still classifies by analytic coordinate thresholds and gives 39/82 cuboids
(48 %). The nearest-face label transfer that produced the recorded 61/82
(74 %) was never committed; it now lives in `clean_blocks.SurfaceLabeller`
and reproduces the documented numbers exactly. Point 3 applies to the
committed classifier too, not only to the abandoned radius split.

**7. The v9 hex mesh is not solid — two internal cavities (now repaired).**
Its boundary has three connected components: the outer surface and two closed
surfaces of 28 and 22 quads, enclosing 0.00170 and 0.00145. That is the
missing 0.2 % of "99.8 % coverage". Both are structured boxes (8 valence-3
corners, everything else valence 4) — 8 and 5 missing cells — and all the
vertices needed to refill them exist. They are the reason 8 blocks (7046
cells, among them the 3395- and 1530-cell ones) are not cuboids. Refilled by
`clean_blocks.fill_cavities` via corner peeling — 14 cells, min scaled
Jacobian 0.848, volume matching the measured cavity to 4 decimals — which
takes the structure to 73/81 cuboids and 98 % of cells. (The second cavity is
1x2x3, not the 1x1x5 first guessed: quad and vertex counts fit both.)

**8. Merging blocks cannot simplify this base complex.**
Sliver-sheet collapse and tiny-block absorption were both rejected in every
single case, always at exactly +3 excess faces. The sheets are bounded by
singular edges, so dropping one merges two blocks without merging their side
neighbours and all four side faces of the union stay split. This is
structural, not a threshold. Real simplification needs Gao et al.'s
mesh-level sheet collapse, which removes a layer of hexes.

Related trap: scoring by the *number* of non-cuboid blocks instead of their
severity looks like a big win (82 → 65 blocks, 78 % → 91 % cuboids) while the
absolute cuboid count falls and 17- and 27-faced blocks appear. The share
only rises because the denominator shrinks.

**9. Sheet collapse is safe here, and one sheet does almost all the work.**
All 145 mesh sheets (edge parallel classes) were screened: none is
self-intersecting and **none adds an inverted cell**. The expected risk —
welding vertices wrecks element quality — did not materialise. But only 1 of
145 improves the block structure, and it does so dramatically: 82 → 42
blocks, excess faces 24 → 10, tiny blocks 8 → 2, inverted cells 2 → 1, and
the boundary moves *closer* to the input surface (0.185 → 0.100) because the
quad holding the old maximum was in the removed layer. A second round finds
nothing. Cuboid share drops 90 % → 86 %, but on half as many blocks.

**11. A quad-diagonal bug mislabelled the v9 input.**
`orient_and_triangulate` reverses an inward-pointing quad and only then
splits it on the 0-2 diagonal — for a reversed quad that is the *other*
diagonal, so its triangles never matched the lookup tables, which registered
only one diagonal. Present in both `tagged_2d_lookup` and
`removed_face_kind`. Effect on the v9 input: `bl_interface_hub` swallowed the
shroud side (10 346 triangles spanning r 0.577–1.804, `bl_interface_shroud` empty)
instead of 2178 + 7912. The corrected counts match the MSH's own wall
triangulations exactly (hub 2178, shroud 7912), which is the cross-check that
settles it — a prism layer's inner interface carries the wall's own
triangulation.

Two earlier entries here were consequences of this bug, not facts about the
geometry: `bl_interface_hub` being one connected shell wrapping both sides, and the
`bl_interface_hub|ogrid_interface` staircases of point 10 below — the feature curve
AlgoHex was given was itself wrong. A rerun on the corrected input is the
cheapest open experiment on this branch.

**10. The last inverted cell needs the boundary to move — but not to deform.**
Smoothing that only moves *interior* vertices cannot untangle it: 4 of its 8
vertices are on the domain boundary and the inversion is already in its
boundary quads, so it stalls at −0.0451. Letting those 4 slide *along* the
input surface untangles it (min scaled Jacobian **+0.0143**, mesh VALID) with
a measured boundary drift of 0.000000 and unchanged Hausdorff — the surface
constraint is what makes the extra freedom safe.

And a metric trap that recurred three times in this stage: maximising the
soft *minimum* scaled Jacobian lifts the worst cell (−0.0582 → −0.0453) while
pushing three neighbours below zero, turning 1 inverted cell into 4. When the
quantity to minimise is a *count*, an extremum or an average will happily
report progress while the structure degrades. Use a one-sided barrier and
rank by the count.

---

## Reduced-domain strategy

Rather than making the frame field resolve the thin boundary layers, they are
cut out of the AlgoHex domain and re-attached afterwards as their own blocks.
Conformity is not an issue: TFI regenerates each block's interior from its
boundary curves, so the parts only have to agree at **block** level, where
the conforming-division MILP (`dp3d/tmesh.py:743`) handles it.

Cell census of the source MSH (measured), which settles what can be reused:

| type | count | median wall distance | role |
|---|---|---|---|
| prism | 121 080 | 0.0120 | boundary layer on hub + shroud |
| hex | 44 800 | 0.047 | O-grid around the blade |
| tet | 89 547 | 0.274 | free core |
| pyramid | 1 792 | 0.286 | O-grid → tet transition |

Faces lying exactly on the cylinders: hub 2178 triangles + 1120 quads,
shroud 7912 triangles + 1120 quads. The quads are exactly the blade O-grid
footprint. A further 2699 (hub) / 2996 (shroud) quads sit near the wall but
on *other* surfaces — they are the quad **side faces** of the prism layer,
which is why the near-wall zone looks fully quad-meshed while the cylinder
surface itself is triangulated.

**Consequence**: the blade O-grid is genuine hexahedra and can be reused as
blocks directly. The hub/shroud boundary layer is triangular prisms, which
are *not* hexahedral blocks, so it must be regenerated — either by tanh wall
clustering inside the AlgoHex blocks (`dp3d/tmesh.py:811 edge_fractions`) or
by extruding from the AlgoHex surface.

---

## Deliverables

`output/hex3d_algohex/deliverable/` (all verified readable with `meshio`):

| file | content |
|---|---|
| `T1_9_hexmesh_v9.{vtk,msh}` | the v9 hex mesh, 56 661 cells |
| `T1_9_hexmesh_v9_quality.vtk` | `scaled_jacobian_x1000` — threshold < 0 shows the 2 bad cells |
| `T1_9_blocks_v9.vtk` | `block_id`, 82 blocks |
| `T1_9_blocks_v9_nfaces.vtk` | `n_block_faces` — threshold ≠ 6 shows the 21 non-cuboids |
| `T1_9_v9_input_*.vtk` | the AlgoHex input: surface (`surface_id`), tets, feature graph |
| `T1_9_hexmesh_v5.*`, `T1_9_blocks_v5.*` | same for the full-domain run |
| `T1_9_blocks_v9_clean.{vtk,msh}` | postprocessed blocks, `block_id`, 81 blocks, cavities filled |
| `T1_9_blocks_v9_clean_nfaces.vtk` | `n_block_faces` — threshold ≠ 6 shows the 8 non-cuboids |
| `T1_9_blocks_v9_clean_edges.vtk` | block-edge wireframe, `block_edge_id`, 423 curves |
| `T1_9_blocks_v9_gao.{vtk,msh}` | same after sheet collapse + untangling: 42 blocks, 54 360 cells, 0 inverted |
| `T1_9_blocks_v9_gao_edges.vtk` | block-edge wireframe, 234 curves |
| `T1_9_blocks_v9_full.{vtk,msh}` | **full domain**: core + blade O-grid + boundary layer, 103 157 cells, **73 blocks**, 5 inverted |
| `T1_9_blocks_v9_full_part.vtk` | `part` — 0 = AlgoHex core, 1 = O-grid, 2 = boundary layer |
| `T1_9_blocks_v9_full_edges.vtk` | block-edge wireframe of the full domain, 678 curves |

The `.msh` files carry the block edges as **1D line elements** next to the
hexes, physical tag = curve id (1-based); the hexes keep `block_id` as their
physical tag. A solid hex mesh shows nothing of the block structure, so the
wireframe is what makes it visible — and it is the entity the conforming
division MILP will tag later. To regenerate it for an existing deliverable
without rerunning the pipeline:

```bash
$PY experimentell/hex3d_algohex/clean_blocks.py --edges-only \
      output/hex3d_algohex/deliverable/T1_9_blocks_v9_gao.vtk
```
| `T1_9_walls_tri_vs_quad.vtk` | `is_quad` on hub+shroud: 1 = triangle, 2 = quad |

Blocks are **volumetric**, not just a surface partition: 50 540 of 60 612
cells are purely interior, and all 27 323 sheet faces lie inside the volume.
ParaView renders only the outer hull of an unstructured grid — use `Clip` or
`Threshold` on `block_id` to see it.

---

## Documents

- `PLAN.md` — original stage plan and background
- `PROGRESS.md` — full chronological log, including failed attempts and why
- `POSTPROCESSING_PLAN.md` — block-structure cleanup (done, stages 5–5c)
- `RUNS.md` — every AlgoHex run: input, domain, runtime, outcome, and the
  mesh/block metrics side by side
- `SHOWCASE.md` — the 11-step visual analysis and what each metric caught
- `ANALYSIS_PLAN.md` — the proposal `SHOWCASE.md` was built from
- `FRAMEFIELD_PLAN.md` — proposed next stage: reducing the block count by
  manipulating the frame field. Includes the singular-graph census, the
  AlgoHex flags that control singular-graph optimisation, and a literature
  review. **Not started.**
