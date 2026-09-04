# Postprocessing of the 3D block structure (T1_9, AlgoHex v9)

## Context

Stage 2 and the block extraction are done: run **v9** (reduced domain, blade
O-grid and hub/shroud boundary layer cut out) gives a hex mesh of 56 661
cells at 99.8 % coverage with only **2** inverted cells, and its base complex
yields **82 blocks, 61 of them cuboids (74 %)**.

That raw structure is not yet usable. Measured on v9:

- **19 sheets**, sizes 2136 … 154, 54, **27** faces — the small ones are
  slivers that separate almost nothing
- **8 blocks under 10 cells** (smallest 3), 36 under 100
- **21 non-cuboid blocks**: 9 touch the *same* surface on two separate
  patches (e.g. `shell_hub` twice), 12 simply have the wrong face count
  (5, 7, 9, up to 12 faces)

This is structurally the same situation the 2D / unwrapped-surface pipeline
already solves: a raw separatrix network is never directly usable either, and
`dp3d` spends most of its code on cleaning it. The 3D block structure needs
the analogous pass. **TFI is explicitly secondary** — it only becomes
meaningful once the block structure is sound.

**Goal**: measurably better than the current state (82 / 61 cuboids, 8 tiny
blocks). Not a hard all-or-nothing criterion.

## What already exists (to mirror, and partly to reuse)

The surface pipeline's postprocessing, which this mirrors one dimension up:

| operation | location | effect |
|---|---|---|
| singularity pair annihilation (Kowalski 2015) | `dp3d/clean_separatrix.py:131` | cancels +1/−1 pairs |
| `_drop_degenerate_corner_seps(min_len=0.1)` | `dp3d/partition_surface.py:769` | drops stub separatrices below a length |
| `_snap_separatrix_endpoints(radius=0.045, bnd_radius=0.05)` | `dp3d/partition_surface.py:922` | snaps ends onto singularity/corner/boundary and **splits** the boundary there so the T-junction becomes a shared graph node |
| `_close_helical_streamlines` | `dp3d/partition_surface.py:604` | closes near-closing curves |
| `collapse_seam_wedges(gap, wall_tol)` | `dp3d/tmesh.py:315` | collapses thin wedges |
| `straighten_sing_connectors(max_len=0.1)` | `dp3d/tmesh.py:193` | straightens short singularity-to-singularity connectors |
| `fix_start_kinks(angle_deg=15, frac=0.25)` | `dp3d/tmesh.py:120` | removes kinks near a curve start |
| `QuadPartitionValidator` | `dp3d/field/quad_partition_validator.py:21` | hard checks (Euler characteristic, valence, inversions, angles, aspect ratio, boundary match) + soft quality score + human-readable diagnostics |

Our 3D side so far: `experimentell/hex3d_algohex/base_complex.py`
(`singular_edges`, `sheet_faces`, `blocks_from_cut`) and `block_faces.py`
(`label_sheets`, `block_faces`, `physical_of`), plus `ovm_io.scaled_jacobian`
and `ovm_io._hex_volume` for the quality gates.

## Design: the 3D analogues

The 2D network is *separatrix curves meeting at singularities*; the 3D
network is *sheets meeting at singular arcs*. Every cleaning operation maps
one dimension up:

| 2D (existing) | 3D (to build) |
|---|---|
| drop stub separatrix (`min_len`) | **collapse a sliver sheet** (few faces / thin) |
| collapse seam wedge (thin quad) | **collapse a chord** (thin block layer) |
| snap endpoint to nearby target (`radius`) | **merge nearly-coincident block corners / short block edges** |
| singularity pair annihilation | **cancel a short singular arc** between two nodes |
| straighten connectors | straighten block edges, smooth block faces |
| `QuadPartitionValidator` | `HexBlockValidator` |

The collapse operations are exactly those of Gao et al. 2017 (*Robust
structure simplification for hex re-meshing*), already cited in
`base_complex.py`. Implemented in Python here rather than pulling in the
reference C++ tool — decided against a second C++ build after the AlgoHex
experience; the tool stays available as a fallback if the quality guarantees
turn out to matter.

## Implementation

New module `experimentell/hex3d_algohex/clean_blocks.py`, deliberately named
and structured after `clean_separatrix.py`, with a thin `postprocess()`
driver mirroring `StreamlinePostProcessor`.

**Step 1 — diagnose (`report_structure`)**
Sheet sizes, block sizes, face count per block, which blocks repeat a
surface label, inverted cells, per-block scaled Jacobian. Prints the same
"before" table used as the baseline to beat. Reuses `label_sheets` and
`block_faces`.

**Step 2 — collapse sliver sheets (`collapse_small_sheets`)**
The 3D analogue of `_drop_degenerate_corner_seps`. A sheet whose face count
is below a threshold (start at ~5 % of the median, i.e. the 27/54/154-face
ones) separates almost nothing: remove it from the cut set, which merges the
blocks on either side. Iterate until no sheet is below threshold.
*Guard*: never remove a sheet whose two sides carry different physical
surfaces — that would merge a block across e.g. the blade shell.

**Step 3 — absorb degenerate blocks (`absorb_tiny_blocks`)**
Blocks below a cell threshold (the 3-, 7-, 9-cell ones) are merged into the
neighbour they share the largest interface with, by dropping that shared
sheet patch. Analogue of the stub filter.

**Step 4 — repair non-cuboid blocks (`fix_non_cuboid`)**
Two distinct causes were measured and need different treatment:
- *duplicate surface label* (9 blocks): the block legitimately touches one
  surface on two opposite sides. Only a counting problem if the two patches
  are genuinely opposite — accept those as cuboids, and only split when the
  patches are adjacent. Fix the classification first, before changing any
  geometry.
- *wrong face count* (12 blocks): split along the sheet that already cuts
  furthest into the block, or absorb if small.

**Step 5 — validate (`HexBlockValidator`)**
Modelled directly on `QuadPartitionValidator` (same three-part shape:
`is_valid()` hard checks, `quality_score()` soft metrics, `diagnostics()`
human-readable):
- hard: every block a topological cuboid (6 faces, each on ONE surface);
  no inverted cells (`ovm_io.scaled_jacobian`); block complex covers the
  whole mesh
- soft: block count, cells per block, min/mean scaled Jacobian, block aspect
  ratio, Hausdorff distance of the block boundary to the input surface
  (`data/T1_9/T1_9_tet_v5.vtk`) so simplification cannot silently deform the
  geometry

**Step 6 — export**
`T1_9_blocks_v9_clean.vtk` / `.msh` with `block_id`, plus a before/after
table. Reuses `ovm_io.write_hex_vtk` / `write_hex_msh` and
`export_vtk.write_vtk`.

### Guard against a mistake already made twice

Never classify a surface by a raw coordinate threshold. It cut across the
triangulation and produced zigzag feature curves in `tet_prep_v5` (p95 kink
174.8°, 318 of 679 nodes above 30°) and dropped the cuboid rate to 44 % in
the v9 block classifier. Always transfer exact labels from the input mesh by
nearest-face lookup, as the corrected v9 run does (74 %).

## Verification

1. `report_structure` on v9 before and after — the "before" numbers are the
   baseline: 82 blocks / 61 cuboids (74 %) / 8 blocks < 10 cells / 2 inverted
   cells. Success = cuboid share up and tiny blocks down, with **0 new**
   inverted cells.
2. `HexBlockValidator.diagnostics()` must not report a block whose face lies
   on two different physical surfaces (that would break TFI later).
3. Hausdorff distance of the simplified block boundary to
   `T1_9_v9_input_surface.vtk` below the tolerance printed in the report —
   catches over-aggressive collapsing.
4. Sanity run on `output/hex3d_algohex/cylinder_hex.ovm`: its base complex is
   the 5-block O-grid, which is already minimal, so postprocessing must leave
   it at 5 blocks and change nothing. A regression here means the collapse
   criteria are too aggressive.
5. Visual: `T1_9_blocks_v9_clean.vtk` coloured by `block_id`, plus the
   `n_block_faces` field to confirm the non-cuboids are gone.

## Afterwards (not this plan)

Assembly with the cut-out boundary layers (blade O-grid blocks, regenerated
hub/shroud layer) and then the conforming-division MILP
(`dp3d/tmesh.py:743 solve_edge_divisions`) plus tanh clustering
(`dp3d/tmesh.py:811 edge_fractions`) for the 3D TFI fill.

---

## Baseline to beat (measured on v9, 2026-09-04)

Reproduce with `block_faces.py output/hex3d_algohex/T1_9_hex_v9.ovm`:

```
19 sheets, 82 blocks, 118 singular edges
sheet sizes (faces): 2136 2010 1758 1323 1056 837 783 702 675 533
                     440 440 374 324 297 270 154 54 27
block sizes:  7650 4250 4250 4080 3400 3395 3150 1750 ... 9 9 7 7 3 3
              8 blocks < 10 cells, 36 blocks < 100 cells
faces per block: {5: 4, 6: 61, 7: 8, 8: 1, 9: 6, 10: 1, 12: 1}
  -> cuboids 61/82 (74 %), covering 48293/56661 cells (85 %)
21 non-cuboid, of which 9 repeat a surface label:
  (block, cells, faces, repeated labels)
  (5, 3395, 9, {shell_hub: 2, shell_blade: 2})
  (19, 1530, 9, {shell_hub: 2, shell_blade: 2})
  (18, 1015, 9, {shell_hub: 2, shell_blade: 2})
  (33,  170, 12, {shell_hub: 4, shell_blade: 3})
  (8,   141, 10, {shell_hub: 2, shell_blade: 2, sheet 13: 2})
mesh quality: 2 cells with scaled Jacobian <= 0, min -0.0607, mean 0.9664
              0 non-manifold faces
```

Success = cuboid share up, tiny blocks down, **0 new inverted cells**.

## Note on the run/file numbering

A run tag `vN` is not the same as an input file `T1_9_tet_vN.vtk`. v9 uses
`T1_9_tet_v5.vtk`. See the mapping table in `README.md` before touching any
input file.

## Environment

```bash
PY=/root/repos/duty/quadmesh/.venv/bin/python
# AlgoHex, if a rerun is needed (sequential only, 7.7 GiB machine):
docker run --rm --network=host -v $PWD:/work \
  -v algohex-build-cache:/app/build algohex-configured \
  /app/build/Build/bin/HexMeshing ...
```
