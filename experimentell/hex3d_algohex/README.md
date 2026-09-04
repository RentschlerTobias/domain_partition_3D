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

**Status**: Stages 0–4 done. A valid hex mesh and a block decomposition
exist. Postprocessing of the block structure is the next step — see
`POSTPROCESSING_PLAN.md`. TFI is after that.

---

## Quick start

```bash
PY=/root/repos/duty/quadmesh/.venv/bin/python

# 1. build the AlgoHex input (reduced domain, the best-performing variant)
$PY experimentell/hex3d_algohex/tet_prep_v5.py            # -> data/T1_9/T1_9_tet_v5.vtk

# 2. run AlgoHex  (~45 min, needs the machine to itself, see "Memory" below)
$PY experimentell/hex3d_algohex/run_algohex.py \
      --tag v9 --in-vtk data/T1_9/T1_9_tet_v5.vtk -- -n 60000

# 3. blocks from the result
$PY experimentell/hex3d_algohex/block_faces.py output/hex3d_algohex/T1_9_hex_v9.ovm

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

v9 is the recommended base: 48× fewer inverted cells and a far coarser block
structure, which is what TFI wants. v5 has the higher cuboid share.

v9's domain excludes both boundary layers, which are meant to be re-attached
as their own blocks afterwards (see "Reduced-domain strategy").

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
| **v9** | `T1_9_tet_v5.vtk` | reduced, exact labels | 550 | 60000 | **completed, 2 inverted cells** |

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
| `export_vtk.py`, `plot_stages.py` | ParaView exports and figures |

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
| `T1_9_walls_tri_vs_quad.vtk` | `is_quad` on hub+shroud: 1 = triangle, 2 = quad |

Blocks are **volumetric**, not just a surface partition: 50 540 of 60 612
cells are purely interior, and all 27 323 sheet faces lie inside the volume.
ParaView renders only the outer hull of an unstructured grid — use `Clip` or
`Threshold` on `block_id` to see it.

---

## Documents

- `PLAN.md` — original stage plan and background
- `PROGRESS.md` — full chronological log, including failed attempts and why
- `POSTPROCESSING_PLAN.md` — the next step, block-structure cleanup
