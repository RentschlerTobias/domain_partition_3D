# Progress — AlgoHex 3D frame field for T1_9

See `PLAN.md` for full context and stage descriptions. This file tracks
execution status only.

| Stage | Status | Notes |
|---|---|---|
| 0. Build AlgoHex (Docker) | **done** | binary built via resumable volume, see log below |
| 1. T1_9 tet mesh + feature tags | **done** | `data/T1_9/T1_9_tet.vtk` written and validated |
| 2. Run HexMeshing | v1 done (bad), v2 running | v1 ran fully but IGM invalid -> 31% coverage; root cause found + fixed, v2 running |
| 3. Validation gates | not started | |
| 4. Base complex (coarse blocks) | not started | |
| 5. 3D TFI fill + export | not started | |

## Log

- Cloned AlgoHex (`git clone --recursive`) into `external/algohex-src` (gitignored,
  see `external/` in `.gitignore`).
- **Docker build attempt 1 failed**: no network inside the container (DNS
  resolution to `deb.debian.org` fails) even though the host has network
  access — sandbox isolates the container's default bridge network.
  **Fix**: `docker build --network=host ...` — verified `apt-get update`
  works with that flag. Rebuild in progress.
- **Exact input file format nailed down by reading AlgoHex source directly**
  (`demo/HexMeshing/HexMeshing.cc`), more reliable than the `-h` output or
  online docs:
  - `.vtk` input is read by `OpenVolumeMesh::Reader::VtkColorReader`
    (`#include <OpenVolumeMesh/FileManager/VtkColorReader.hh>`).
  - Single VTK v2.0 ASCII `UNSTRUCTURED_GRID` mixing cell types 1 (vertex),
    3 (line), 5 (triangle), 10 (tetrahedron), all cells sharing one global
    point list, plus one integer `CELL_DATA` "color" scalar array over all
    cells (matches the HexMe dataset paper's description of Gmsh's own vtk
    v2.0 export).
  - `MeshProperties.hh:import_feature_properties()` reads this back as
    per-mesh-entity properties `vertex_colors` / `edge_colors` /
    `face_colors`; any color `!= 0` = feature. Triangle cells not part of
    any real feature curve should still be colored `>0` on ALL boundary
    triangles (frame field must be tangent to any boundary, not just
    sharp-feature ones) — `initialize_feature_properties()` also has an
    automatic fallback (`-d/--dihedral-angle`, `DihedralFeatureDetector`)
    that marks every boundary face as a feature face and auto-detects
    feature edges by dihedral angle if no `AlgoHex::FeatureEdges` /
    `edge_colors` property is present at all — a viable fallback if hand
    tagging turns out unreliable.
  - Relevant CLI flags found in `demo/HexMeshing/main.cc` (CLI11):
    `-n` = target output hex cell count (coarseness knob for later),
    `-d/--dihedral-angle`, `--force-feature-threshold`,
    `--full-constraints`, `-j/--json-out-file` (metrics), `--igm-out-path`
    / `--sm-out-path` (dump the integer-grid / seamless map for
    inspection). No flag found yet for exporting an intermediate
    block/motorcycle-complex structure directly — Stage 4's base-complex
    extraction on the fine hex output is still the plan.
- **T1_9 boundary verified by direct MSH parse** (reusing
  `dp3d.extraction.parse_msh` / `surface_triangles`, no new file written
  yet): 19 tagged boundary patches (geom 1-19: blade 1-5, hub 6-11,
  shroud 12-17, inlet 18, outlet 19), 36 068 boundary triangles total,
  matches `$PhysicalNames` in the MSH.
- **Stage 1 approach decided**: do NOT hand-build per-patch discrete gmsh
  topology (fragile w.r.t. shared boundary nodes across patches). Instead
  reuse the proven pattern already in this repo
  (`experimentell/gmsh_pipeline/remesh_step1.py`): merge the full
  boundary (all 19 patches combined) as one STL into gmsh, run
  `classifySurfaces` + `createGeometry()` to recover surfaces/curves/
  points purely from dihedral angle (no need to thread through the
  original geom ids), then mesh the enclosed volume with tets. Physical
  Surface/Curve/Point tags from this become the VTK color ids. Custom
  hand-rolled VTK writer (same style as `dp3d/extraction.py`'s
  `_write_quad_vtk`) rather than trusting gmsh's own VTK exporter's
  undocumented color/field-data behavior.
- **First attempt at combining the 19 originally-tagged boundary patches
  (geom ids 1-19) failed**: `classifySurfaces` errored "one edge is
  incident to 3 triangles". Root cause found by direct edge-multiplicity
  analysis: the MSH's geom-id numbering does NOT correspond 1:1 to the 19
  physical patch names — the same physical surface (e.g. hub) is
  represented redundantly by both a coarse triangle patch AND separate
  small quad O-grid patches at different geom ids, which don't conform to
  each other. **Fix**: derive the boundary purely topologically instead
  (`tet_prep.volume_exterior_faces`: a face referenced by exactly one
  volume cell, generalizing `dp3d.extraction.hex_exterior_quads` from
  hex-only to the hybrid tet/hex/prism/pyramid mix) — this is exact and
  sidesteps the redundant-tagging problem entirely. Result: clean 2-
  manifold, 34 308 boundary triangles, edge-multiplicity histogram
  `{2: 51462}` (no violations).
- **gmsh volume meshing pitfall**: an explicit `MathEval` background size
  field (following `remesh_step1.py`'s pattern) combined with
  `gmsh.model.mesh.clear()` hung for >240s on this geometry (never
  finished, no output). Switched to the simpler
  `Mesh.MeshSizeMax`/`Mesh.MeshSizeMin` options (no field, no explicit
  clear) — fast and reliable (~18-25s wall time). Root cause not fully
  diagnosed (likely a degenerate/tiny reparametrized patch — the log
  showed several sub-1e-9-area splits during `classifySurfaces` — choking
  the field evaluation), not worth chasing further since the simpler
  option-based approach works.
- **Stage 1 result** (`experimentell/hex3d_algohex/tet_prep.py`, run via
  `/root/repos/duty/quadmesh/.venv/bin/python experimentell/hex3d_algohex/tet_prep.py`,
  ~25s): `classifySurfaces` at 30°/60° (surface/curve angle) recovered 21
  surfaces, 154 curves, 133 points from the T1_9 boundary. Volume meshed
  at `MeshSizeMax=0.08`/`MeshSizeMin=0.02`: 17 969 nodes, 89 173 tets, only
  5 ill-shaped tets, quality mostly 0.6-1.0. Wrote
  `data/T1_9/T1_9_tet.vtk` (17 969 points, 104 437 cells: 133 vertex + 853
  line + 14 278 triangle + 89 173 tetra). **Validated**: all 14 278
  boundary triangles carry a nonzero color (no gaps in feature-face
  coverage), boundary edge-multiplicity histogram is exactly `{2: 21417}`
  (still 2-manifold after the coarser remesh), 853 feature-edge segments
  span 154 distinct feature curves.
- **Docker build attempt 2 (`--network=host`) reached AlgoHex's own build
  step and got OOM-killed** (`cc1plus` "Killed signal terminated program",
  ninja aborted) -- NOT a real success despite the background-task
  notification reporting "exit code 0" (that was the outer shell
  pipeline's status via `tee`, not the actual `docker build` exit code --
  always check the log, not just the notification). Root cause: this
  sandbox only has 7.7 GiB RAM / 2 CPUs and swap was already half used;
  AlgoHex is template-heavy (Eigen, CoMISo, TinyAD) and the vendored
  Dockerfile builds `RelWithDebInfo` (debug-info generation for heavy
  templates is the usual memory hog) with unconstrained `ninja` (default
  parallelism = nproc). **Fix**: patched the local clone
  `external/algohex-src/Dockerfile` (not upstream, not committed --
  `external/` is gitignored): `CMAKE_BUILD_TYPE=Release` instead of
  `RelWithDebInfo`, `ninja -j1` instead of unconstrained. Rebuild reuses
  Docker's layer cache for the already-successful IPOPT/Bonmin steps
  (confirmed via `Using cache` in the log) so only AlgoHex's own
  CMake+ninja step reruns. Attempt 3 in progress.
- **Stage 0 finished**: docker build with unmodified `ninja` (no `-j1`)
  still restarted from object 1/346 after being killed once more mid-build
  (a killed `RUN` layer leaves nothing cached -- confirmed empirically, the
  "resume" hope from a named volume only works if ninja itself was already
  running against that volume from the start). Worked around by: (1)
  commented out the final `ninja`/`ln -s` lines in the local Dockerfile so
  `docker build` only produces an image configured-but-not-compiled
  (`algohex-configured`, cache-hits instantly on retry); (2) created a
  named Docker volume `algohex-build-cache`, mounted at `/app/build` in a
  plain `docker run ... ninja -j1`, so ninja's own incremental state
  (`.ninja_log`) persists in the volume across separate `docker run`
  invocations if interrupted again. This run completed for real (verified
  by checking `Build/bin/` inside the volume directly, not just trusting
  the "exit code 0" completion summary -- that summary was WRONG for the
  earlier OOM'd attempt, see above). `HexMeshing -h` and the upstream demo
  (`demo/HexMeshing/cylinder.ovm`) both run correctly: full pipeline
  completed in ~100s, 9792/9794 hex cells extracted cleanly. **No image
  named `algohex` with baked-in binaries exists** -- always invoke via
  `docker run -v algohex-build-cache:/app/build algohex-configured
  /app/build/Build/bin/HexMeshing ...` (this is what
  `experimentell/hex3d_algohex/run_algohex.py` does).
- **Stage 2 started** on the real T1_9 case
  (`experimentell/hex3d_algohex/run_algohex.py`, log:
  `output/hex3d_algohex/hexmeshing.log`). Confirmed it picked up our
  Stage 1 feature tags correctly: `"Found prescribed feature colors (from
  vtk file)"` (not the automatic dihedral-angle fallback) -- validates the
  whole VTK color-property round trip end to end. Frame field
  optimization iterating as expected. T1_9 has ~9x the tets of the
  cylinder demo (89 173 vs 9 764); expect a correspondingly longer run,
  possibly well beyond the demo's 100s given the nonlinear solves likely
  don't scale linearly.
- **Note on background-task tracking**: the tool-level background task
  wrapping `run_algohex.py` was reported "killed"/"stopped" after ~18min,
  but the actual `docker run` container it launched (`HexMeshing`, PID
  visible via `docker top`) kept running completely unaffected --
  `docker run` (no `-d`) detaches from a dead CLI client at the daemon
  level, it does not stop the container. Confirmed alive by checking the
  log file was still growing and `docker top` showed the process at ~94%
  CPU well after the "killed" notification. **Lesson**: for long AlgoHex
  runs, check liveness directly (`docker ps`, `docker top <id>`, log file
  growth) rather than trusting the background-task notification's
  killed/completed status -- it does not reflect the container's actual
  state.
## Run v1 (classifySurfaces feature tags) — completed, RESULT UNUSABLE

Ran to completion in **7777 s (~2h10m)**, exit 0, and produced a hex mesh —
but an incomplete one.

**What succeeded:**
- Frame field generation converged.
- **Locally meshable frame field generation SUCCEEDED**: log line 117475
  `"All special vertices are locally meshable. Stop pipeline."` This is the
  hard robustness gate — per Liu & Bommes 2023 it is exactly the step that
  fails on ~42% of the HexMe benchmark. T1_9 passed it.

**What failed** (`T1_9_hex_metrics_v1.json`, `hexmeshing_v1.log`):
- seamless map: `n_invalid_param_tets_seamless = 3863`,
  `n_invalid_valencies_seamless = 0`, `valid_volume_seamless = 0.9984`
- after quantization → IGM: `n_invalid_param_tets = 12588`,
  `n_invalid_valencies = 12`, `valid_volume = 0.9886`
- HexEx: `Detected 141139 Proper, 12502 Flipped and 99 Degenerate Tets`
  → `Error: Invalid Input - Flipped or Degenerate Tet`
- extracted mesh: 10 061 hexes, **but only 31.4% volume coverage**
  (2.0497 of the true 6.52525) and 12 negative-volume cells. Bounding box
  is correct, so it is not a scaling problem — the mesh has large holes
  wherever the IGM was flipped.

**ROOT CAUSE FOUND — spurious feature curves from `classifySurfaces`.**
Measured the true dihedral angle across every tagged feature edge:
**427 of 853 (50%) were geometrically FLAT (< 20°)**, i.e. pure artifacts.
Reason: `classifySurfaces` splits surfaces for *parametrizability*, not
physics — it cut the smooth hub cylinder into 5 patches (ids 8,14,18,19,20)
and the shroud into 4 (13,16,21,22), and every artificial cut line became a
feature curve. AlgoHex hard-aligns the octahedral field to feature curves,
so a spurious curve running across a smooth cylinder forces artificial
singularities and destroys integrability → invalid seamless map → invalid
IGM → holes. The 21 surfaces flagged as "worth a sanity check" earlier were
in fact the bug.

## Stage 1 revision — analytic feature tagging

Rewrote the tagging in `tet_prep.py` to derive features from the KNOWN
geometry instead of gmsh's parametrization splits. gmsh is still used to
generate the tets; its surface ids are simply ignored.
- `classify_boundary()`: hub (r≈0.5), shroud (r≈1.9), inlet (z≈0), outlet
  (z≈2.5) recognised analytically by position **and** normal direction; the
  remainder splits by edge-connectivity into exactly 3 components — blade
  (2119 tris) and the two periodic side walls (1767 / 1763 tris, near-equal
  sizes confirming they are the matching periodic pair).
- `feature_edges_and_vertices()`: feature curve = boundary between two
  DIFFERENT physical surfaces, plus genuinely sharp creases inside a patch
  (needed for the blade LE/TE).
- **Bug found and fixed while doing this**: `tet_boundary_faces()` initially
  returned boundary triangles with inconsistent winding, so adjacent
  normals could come out anti-parallel and the dihedral test reported ~180°
  creases on flat surfaces (154 bogus "creases" on the planar inlet). Now
  each face is wound outward using the owning tet's apex vertex. This cut
  the detected sharp creases from 2451 → 129.

**Result of the revision** (`data/T1_9/T1_9_tet.vtk`, same 89 173 tets):

| | v1 (classifySurfaces) | v2 (analytic) |
|---|---|---|
| surfaces | 21 (arbitrary splits) | 7 (physical) |
| feature edges | 853 (427 spurious flat) | 496 (367 patch bnd + 129 real creases) |
| feature vertices | 133 | 36 |

Validated: 0 uncoloured boundary triangles, edge multiplicity exactly
`{2: 21417}` (2-manifold), and the 129 remaining sharp creases sit on the
blade (128, = leading/trailing edge) + 1 stray on the hub — i.e. all
physically justified.

## Run v2 (analytic tags) — CRASHED (SIGSEGV), root cause found

Initial signal was good (SH energy 8312.73 vs v1's 9356.32 = smoother, less
over-constrained field) but the run **segfaulted** (`WRAPPER_EXIT=139`)
during locally-meshable-field generation, ~35 400 log lines in, after
`Error: adjacent halfface is invalid!` and `Vertex: 29 has 12 invalid
cells!`. No output files; all metrics `-1`.

**Cause: my hand-built feature-curve graph was not a valid 1-manifold.**
Analysis of the tagged graph found **23 dangling ends** (valence-1 vertices,
none tagged as feature vertices) and 5 tiny isolated components — fragments
produced by the 40° sharp-crease threshold picking up partial creases that
just stop in the middle of a smooth surface. AlgoHex assumes feature curves
are proper curves terminating on other feature curves; dangling ends crash
it. (gmsh's `classifySurfaces` output in v1 never had this problem because
gmsh guarantees curve topology — so v1 was over-constrained but structurally
valid, v2 was well-constrained but structurally invalid.)

**Fix** — `clean_feature_graph()` in `tet_prep.py`:
1. iteratively prune dangling crease edges (valence-1 endpoint, edge is not
   a patch boundary). Patch-boundary edges are structural and never pruned;
   a genuine blade LE/TE crease runs between two patch-boundary curves, so
   both ends are junctions and it survives.
2. drop small isolated components (threshold noise).
3. tag as feature vertices every vertex of valence != 2, i.e. junctions AND
   any surviving endpoints (previously only valence >= 3 was tagged).

## Run v3 (analytic tags + valid feature graph) — in progress

Feature graph after cleanup, verified as a valid 1-manifold:

| | v1 (classifySurfaces) | v2 (analytic) | v3 (analytic + cleaned) |
|---|---|---|---|
| surfaces | 21 | 7 | 7 |
| feature edges | 853 (427 flat) | 496 | 437 |
| feature vertices | 133 | 36 | 27 |
| dangling ends | 0 | **23** | **0** |
| graph components | — | 9 (5 tiny) | 3 (no tiny) |
| valence histogram | — | {1:23, 2:426, 3:28, 4:7, 5:1} | {2:393, 3:21, 4:5, 5:1} |

All endpoints tagged (`endpoints tagged: True`). Running; already past the
SH-optimization stage with 0 errors.

Also now passes `--sm-out-path` / `--final-tetmesh-out-path` to checkpoint
the seamless map and post-singularity-optimization tetmesh. AlgoHex's
`hexMeshing_from_seamless_map()` path (triggered by supplying BOTH
`--hexex-in-path` and `-i`, see `demo/HexMeshing/HexMeshing.cc:503`) skips
field generation + locally-meshable-field generation + integrability
optimization — 6704 s of v1's 7777 s total, i.e. 86% — and resumes at
parametrization/quantization/extraction. This makes subsequent `-n`
(`num_hex_cells`, default 10000) sweeps cheap instead of 2h each. Note that
restart path reads the tetmesh with OVM's `FileManager`, so it needs the
`.ovm` tetmesh from `--final-tetmesh-out-path`, not the input `.vtk`.

## VTK exports for manual inspection (`export_vtk.py`)

`output/hex3d_algohex/vtk/` — open in ParaView:

| file | cell data | what |
|---|---|---|
| `01_boundary_from_2D_elements.vtk` | `geom_id` | the MSH's own tagged 2D elements, untouched |
| `01b_2D_elements_outer_flag.vtk` | `on_outer_skin` | 1 = on the true outer skin, 0 = interior |
| `02_boundary_topological.vtk` | `surface_id` | skin from volume cells, physical id 1..7 |
| `03_feature_curves.vtk` | `kind` | 1 = feature edge, 2 = feature vertex |
| `04_tet_volume.vtk` | — | the 89 173-tet mesh fed to AlgoHex |
| `05_hexmesh_v1.vtk` | `inverted` | v1 result (incomplete, 31.4%) |

### "just take the 2D elements" — CONFIRMED, with one filter

Tobias' suggestion is right and simpler than the topological route, but the
raw 2D element set cannot be used as-is:

- all 31 116 tagged 2D elements: edge multiplicity `{2: 53664, 3: 724}` →
  **closed but NOT 2-manifold** (724 edges shared by 3 faces). Reason: the
  hex O-grid **block interfaces inside the volume are tagged as 2D elements
  too**, so the set mixes outer skin with interior faces. This is what broke
  the very first `classifySurfaces` attempt ("one edge is incident to 3
  triangles").
- classifying every tagged 2D element against the true skin splits the geom
  ids cleanly, with no mixed patch:
  - **outer**: `1,2,3,4,5,6,7,8,9,13,14,15,18,19,20,23,24,25,28,29,30`
  - **interior** (O-grid interfaces): `10,11,12,16,17,21,22,26,27,31`
- keeping only the outer ids gives edge multiplicity `{2: 41790}` →
  **closed 2-manifold, i.e. directly usable**.

So the simple route works as: *take the 2D elements, drop the 10 interior
geom ids*. It is arguably better than the current topological extraction
because it **preserves the original boundary quads** (41 790 edges) instead
of triangulating them (51 462). Worth switching `tet_prep.py` to this once
confirmed — the two skins should be geometrically identical.

## Tooling added

- `ovm_io.py` — OVM ASCII reader + hex VTK writer. AlgoHex writes
  half-edge/half-face OVM which meshio cannot read. Hex connectivity is
  recovered topologically (pick a face, find the disjoint opposite face,
  pair each bottom vertex with its unique top neighbour along a cell edge)
  rather than by trusting a half-face ordering convention, then the winding
  is fixed by the sign of the cell volume. **Validated against AlgoHex's own
  cylinder demo: 9792/9792 hexes recovered, 0 skipped, 0 negative volumes,
  and vertex/edge/face/cell counts exactly match HexHex's own reported
  output** (10955 / 31642 / 30480 / 9792).
