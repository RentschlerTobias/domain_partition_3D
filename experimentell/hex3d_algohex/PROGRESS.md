# Progress — AlgoHex 3D frame field for T1_9

See `PLAN.md` for full context and stage descriptions. This file tracks
execution status only.

**See `README.md` for the branch overview, the run/input-file mapping and
the key findings. This file is the chronological log, including the failed
attempts and why they failed.**

| Stage | Status | Notes |
|---|---|---|
| 0. Build AlgoHex (Docker) | **done** | resumable volume build; `algohex:portable` image also exists |
| 1. Tet mesh + feature tags | **done** | best variant `tet_prep_v5.py` -> `data/T1_9/T1_9_tet_v5.vtk` |
| 2. Run HexMeshing | **done** | v5 (full domain) and v9 (reduced) both usable; v9 has 2 inverted cells vs v5's 96 |
| 3. Validation gates | **done** | scaled Jacobian, manifoldness, coverage, boundary conformity |
| 4. Base complex (coarse blocks) | **done** | v9: 19 sheets -> 82 blocks, 61 cuboids (74%) |
| 5. Block postprocessing | **next** | see `POSTPROCESSING_PLAN.md` |
| 6. Assembly + 3D TFI fill | not started | re-attach blade O-grid, regenerate hub/shroud BL, MILP + TFI |

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

## Run "fast" (`--without-integrable-field-optimization`) — FAILED, but informative

Launched in parallel with v3 on the second core to get a same-day answer,
using the corrected v3 feature tags but skipping the 86%-of-runtime
integrability phase. Result: **SIGSEGV after 3326 s**, 166 404 `ERROR :(`
lines during HexEx extraction, **no output mesh at all**.

Parametrization metrics vs v1 (`T1_9_hex_metrics_fast.json`):

| | v1 (bad tags, WITH integrability) | fast (good tags, WITHOUT) |
|---|---|---|
| invalid tets (seamless) | 3 863 | **16 882** |
| invalid valencies (seamless) | 0 | **137** |
| invalid tets (IGM) | 12 588 | **37 538** |
| invalid valencies (IGM) | 12 | **331** |
| valid_volume (IGM) | 0.9886 | 0.9625 |
| parametric_volume | 10 190 | **-65 898** (negative!) |

**Conclusions:**
1. `--without-integrable-field-optimization` is **not a usable shortcut** on
   this geometry. It is not merely slower-but-worse, it is catastrophic:
   negative parametric volume and a crash in extraction.
2. **This corrects an over-attribution made earlier.** The spurious feature
   curves were a real, measured defect (427 of 853 edges geometrically
   flat) and fixing them is justified — but this run shows the
   *integrability optimization* dominates IGM validity far more than tag
   quality does: with good tags but no integrability phase the seamless map
   is 4x worse than v1's (16 882 vs 3 863 invalid tets). So the earlier
   claim that bad feature tags were "the root cause" of v1's 31.4% coverage
   was too strong; they were *a* contributing defect, not the whole story.
3. Therefore the fast probe did **not** answer the open question (does the
   tag fix repair extraction?) — it removed the phase that matters most.
   Only v3 can answer it.

## Run v3 (analytic tags, full pipeline) — COMPLETED THE HARD PHASES, IGM BLEW UP

Ran 11 h, cleared local meshability, completed all 9 integrability rounds,
wrote the seamless-map checkpoint (132 MB) — then the **IGM diverged
numerically** and HexEx refused the input.

| | v1 (classifySurfaces tags) | **v3 (analytic tags)** | fast (no integrability) |
|---|---|---|---|
| invalid tets (seamless) | 3 863 | 17 878 | 16 882 |
| **invalid valencies (seamless)** | **0** | **10** | 137 |
| invalid tets (IGM) | 12 588 | 103 867 | 37 538 |
| invalid valencies (IGM) | 12 | 870 | 331 |
| valid_volume (IGM) | 0.9886 | **0.4606** | 0.9625 |
| final_energy (IGM) | 91.0 | **3.26e+48** | 45.6 |
| parametric_volume | 10 190 | **6.08e+62** | -65 898 |
| HexEx | 12 502 flipped | 103 647 flipped | crash |

**The analytic-tag hypothesis is DISPROVEN. v1's "bad" tags are the best of
the three.** Stated plainly because two earlier entries in this file claimed
the opposite: the spurious feature curves were a genuine, measured defect
(427/853 edges flat) but removing them made the *result* worse, not better.

**The decisive metric is `n_invalid_valencies_seamless`** — invalid vertex
valencies in the seamless map, i.e. singularity configurations that are not
hex-meshable. v1: **0**. v3: 10. fast: 137. This gates everything downstream:
`hexMeshing_from_seamless_map()` (and the equivalent code in the main path)
only runs `parametrize_robust_quantization()` when the seamless map comes
back fully valid. With 10 invalid valencies v3's quantization ran off the
rails (energy 1e48, parametric volume 1e62).

**Why the extra feature curves helped**, most plausible reading: they act as
regularisation. Pinning the octahedral field along many extra curves
constrains it into a simpler, hex-meshable singularity graph. Removing them
gave the field more freedom and it relaxed into a configuration that is
smoother but *not* hex-meshable. Fewer, cleaner constraints are not
automatically better for integer-grid maps.

v3's extraction was killed after hours stuck at `Processing edge 0 of
238061` — with 103 647 flipped tets it could not produce anything usable.
Outputs archived as `*_v3.*`.

## Run v4 (v1 tags + `-n 60000`) — in progress

Back to the tag set that produced a hex-meshable seamless map, changing the
one lever that plausibly explains v1's *extraction* failure: `-n`
(`num_hex_cells`, default **10000**). For a thin, twisted blade passage a
10k-cell integer grid is coarse — thin regions quantize to zero thickness
and flip, which matches v1's symptom (valid seamless map, only 12 invalid
IGM valencies, yet 31.4% coverage with holes). Raising to 60000 gives the
quantization enough resolution to keep thin regions at >= 1 layer.
Input: `data/T1_9/T1_9_tet_v1tags.vtk`. Checkpointing enabled.

## CORRECTION — the "31.4% coverage" figure was a bug in MY volume function

**Every earlier statement in this file about coverage / holes / inverted
cells is wrong and is superseded by this section.**

`ovm_io._hex_volume` returned exactly **1/3** of the true volume: unit cube
→ 0.3333, 2x3x4 box → 8.0 instead of 24.0. Cause: two of the five tets in
the decomposition had their last two vertices swapped — `(0,5,4,7)` instead
of `(0,5,7,4)` and `(2,7,6,5)` instead of `(2,7,5,6)` — so they contributed
with the wrong sign. The tell was there all along and I missed it: the
*cylinder demo*, a known-good AlgoHex result, also came out at "32%".

This bug is **inherited from `experimentell/3d_extrapolation/hexa_interpolation.py`
(`hexa_cell_volumes`, lines 224-237), which has the same swapped tets.**
There it is largely harmless (the value only feeds volume statistics; cell
validity is judged by `hexa_corner_jacobians`), but worth fixing — flagged,
not yet changed, since it is outside this branch's scope.

### Corrected results (domain volume 6.5253)

| | cells | coverage | inverted (scaled Jacobian <= 0) |
|---|---|---|---|
| cylinder demo (reference) | 9 792 | 98.5 % | — |
| **v1** (`-n` 10000, classifySurfaces tags) | 10 061 | **95.4 %** | — |
| **v4** (`-n` 60000, same tags) | 60 100 | **97.9 %** | **124 (0.21 %)** |

The residual ~2 % is largely the piecewise-flat hex boundary approximating
curved cylinders — the perfect cylinder demo shows the same 1.5 % deficit.

**So Stage 2 did not fail. v1 was already a 95 %-complete hex mesh.** The
"holes" narrative, and the plot `step04_hexmesh_v1.png` captioned
"INCOMPLETE ... 31.4%", were artifacts of the volume bug. What HexEx reports
as `Error: Invalid Input - Flipped or Degenerate Tet` refers to its *input
IGM*; it sanitizes and still extracts a usable mesh.

### v4 validity audit (the honest numbers)

- boundary: 10 764 faces, **0 non-manifold faces**; boundary edge
  multiplicity `{2: 21522, 4: 3}` → 3 defective edges
- boundary vertices land on the domain: 1955 on hub r=0.5, 1999 on shroud
  r=1.9, 1386 on inlet, 935 on outlet, 4488 on blade/periodic
- **scaled Jacobian: mean 0.9314, min -0.9956; 124 cells <= 0**, 201 <= 0.1
- note: signed *volume* said "0 inverted" for these same cells — a sheared
  hex can have positive volume with inverted corners, so scaled Jacobian
  (corner Jacobians) is the correct validity test, exactly as the docstring
  of `hexa_corner_jacobians` in `hexa_interpolation.py` already states.

Verdict: v4 is a good, near-complete, mostly high-quality hex mesh, but
**not yet valid** — 124 inverted cells and 3 non-manifold boundary edges
must be repaired before TFI / CFD use.

### What this means for the earlier conclusions

- The v3 comparison **still stands**: it is based on AlgoHex's own JSON
  metrics (`n_invalid_valencies_seamless` 0 vs 10, IGM energy 91 vs 3.3e48),
  not on my volume function. Analytic tagging still made the IGM worse.
- But the framing "v1 failed, so the tags must be at fault" was built on the
  bad coverage number. v1 had in fact essentially succeeded, so there was
  never a failure to explain.

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

## Stage 4 status: partition works, structured block extraction does NOT yet

**Works and is validated** (`base_complex.py`): partitioning the hex mesh
along base-complex sheets. Cylinder demo -> exactly 5 components, matching
the textbook O-grid (core + 4 quadrants), confirmed geometrically:
core r 0.68-7.37, three quadrants spanning ~82-87 deg in theta and one
crossing the +-180 deg branch cut. T1_9 v5 -> 307 components.

**Does NOT work yet** (`block_structure.py`): turning a component into a
structured (i,j,k) block, which is what TFI actually needs. The grid walk
still reports coordinate conflicts on the cylinder core, i.e. different
paths to the same cell disagree, so no consistent (i,j,k) is produced.
Current state: 0 of 5 cuboids on the cylinder -- impossible for an O-grid,
so the walk is still wrong, not the mesh.

Two real bugs were found and fixed along the way (both mine):
- `AX_FACES` and `AX_EDGES` were indexed inconsistently: local axis a is the
  NORMAL of face pair `AX_FACES[a]`, so the parallel edges of axis a are the
  ones running along that normal. The tables were offset, so a face's own
  edges were classified as its normal axis and frame propagation collapsed
  after 3 cells. After the fix the walk reaches all 4896 cells of the block.
- the sign of the frame transfer was copied verbatim instead of being
  derived from the edge orientation; now `ss * oa * ob`.

Two earlier attempts at a quick cuboid test were also wrong and are recorded
so they are not repeated:
- counting edge parallel classes and expecting 3: a cuboid actually has
  `i+j+k` classes, one per layer, so 49 classes for the cylinder core is
  correct, not a defect.
- counting boundary patches with a 45 deg coplanarity test: unreliable on
  curved hub/shroud walls.

**So: the answer to "is extracting the hexa blocks the next step" is yes,
and it is genuinely unfinished work, not a formatting step.** What remains:
1. fix the grid walk so each component gets consistent (i,j,k) (or detect
   components that are not cuboids and split them);
2. split blocks whose boundary spans several physical surfaces, so every
   block face lies on exactly one of hub/shroud/blade/inlet/outlet/periodic;
3. only then the existing conforming-division MILP (`dp3d/tmesh.py:743`)
   and tanh clustering (`edge_fractions`, `dp3d/tmesh.py:811`) can be reused
   for the 3D TFI fill.

Also still open on the mesh itself: 96 inverted cells (scaled Jacobian <= 0)
in v5 despite a fully valid IGM.

## Reduced-domain runs (v6, v7) — Tobias' cut-out strategy

Idea: do not make the frame field resolve the boundary layers at all. Cut
them out of the AlgoHex domain, EXTRACT blocks from the result, and assemble
them with the boundary layers' own block structure afterwards. Conformity is
a non-issue because TFI regenerates each block's interior from its boundary
curves, so the parts only have to agree at BLOCK level, where the
conforming-division MILP (dp3d/tmesh.py:743) handles it. My earlier
objection about a non-conforming cell-level interface was reasoning at the
wrong level and is withdrawn.

Domains (all verified closed 2-manifold):

| domain | cells | boundary tris | tets fed to AlgoHex |
|---|---|---|---|
| full | 257 219 | 34 308 | 71 415 (v5) |
| minus blade O-grid (v6) | 212 419 | 29 828 | 62 874 |
| minus O-grid AND prism BL (v7) | 91 339 | 18 548 | 42 218 |

The v7 domain detaches from the walls entirely (r 0.577-1.804 instead of
0.498-1.900), so no wall, blade edge or hub/shroud intersection curve is
left. Its feature graph is 166 edges with valence uniformly 2 — only closed
loops, no branches, no dangling ends: by far the cleanest input so far.

Convergence to "All special vertices are locally meshable":
v1/v4 ~2 h (42 repair iterations), v5 ~15 min (32), v6 ~11 min, v7 ~8 min.
Removing the boundary layers helps markedly.

**v7 was OOM-killed** (`exit=137`, kernel: `Out of memory: Killed process
HexMeshing total-vm:5841700kB`) at the quantization stage after 27 min --
my fault for running v6 and v7 concurrently on a 7.7 GiB box. AlgoHex peaks
hard during quantization; these runs must be sequential. v6 was stopped
(round 3 of 9 after 1 h) and v7 restarted alone.

### Cell types in the source MSH (measured, settles what can be reused)

| type | count | median wall distance | role |
|---|---|---|---|
| prism | 121 080 | 0.0120 | boundary layer on hub+shroud |
| hex | 44 800 | 0.047 | O-grid around the blade |
| tet | 89 547 | 0.274 | free core |
| pyramid | 1 792 | 0.286 | O-grid to tet transition |

Prisms have height/sqrt(base area) = 0.20, i.e. flat stacked layers.

Faces lying exactly ON the cylinders (all nodes at r=const, radial normal):
hub 2178 triangles + 1120 quads, shroud 7912 triangles + 1120 quads. The
quads are exactly the blade O-grid footprint (geom 7/13/18/23/28 and
8/14/19/24/29). Additionally 2699 (hub) and 2996 (shroud) quads sit close to
the wall but on OTHER surfaces, with zero triangles among them -- these are
the quad SIDE faces of the prism layer on inlet/outlet/periodic. That is why
the near-wall zone looks fully quad-meshed while the cylinder surface itself
is triangulated.

Consequence for assembly: the blade O-grid is genuine hexahedra and can be
reused as blocks directly. The hub/shroud boundary layer is triangular
prisms, which are NOT hexahedral blocks, so it has to be regenerated --
either by tanh wall clustering inside the AlgoHex blocks (`edge_fractions`,
dp3d/tmesh.py:811) or by extruding from the AlgoHex surface.

## v9 — reduced domain with exact labels: the memory blowup was a LABELLING problem

Same reduced domain as v7/v8 (no blade O-grid, no prism boundary layer,
42 218 tets), only the boundary labelling changed: faces matched against the
MSH's own tagged 2D elements (geom 3 inlet, 4 outlet, 5/6 periodic), the
newly exposed shell split by which cell type was removed behind it.
550 feature edges over 6 surfaces instead of 166 over 3.

**It completed** (exit 0, 45 min) where v7 and v8 both OOM-died at the
identical point. So the OOM was not a memory-size problem and not driven by
`-n` (which changed nothing between v7 and v8) -- it came from an
under-constrained field. This is the fourth time in this project that fewer
feature constraints produced a worse result.

### v9 vs v5 (the two usable results)

| | v5 (full domain) | v9 (reduced) |
|---|---|---|
| cells | 60 612 | 56 661 |
| coverage | 100 % of 6.5253 | 99.8 % of 5.4606 |
| **inverted cells (scaled Jacobian <= 0)** | **96** | **2** |
| worst scaled Jacobian | -0.9956 | **-0.0607** |
| mean scaled Jacobian | 0.9651 | 0.9664 |
| non-manifold faces | 0 | 0 |
| runtime | 1 h 17 | 45 min |
| blocks | 307 | **82** |
| cuboids | 290 (94 %) | 61 (74 %), 85 % of cells |

v9's mesh is markedly cleaner (48x fewer inverted cells, and the two that
remain are barely inverted) and its block structure is much coarser, which
suits TFI better. v5 has the higher cuboid rate.

**IGM validity does not predict mesh quality.** v9's parametrization was the
worse one by AlgoHex's own numbers (9030 invalid tets vs v5's 0, HexEx
reporting `Invalid Input`), yet its extracted mesh has 2 inverted cells
against v5's 96. Judge the extracted hexes, not `valid_volume`.

### Repeated own-goal worth recording

Splitting a surface by a coordinate threshold cuts across the triangulation
and creates zigzag feature curves. It happened twice: first in tet_prep_v5
(radius threshold, p95 kink 174.8 deg, 318 of 679 nodes above 30 deg; fixed
by splitting on removed-cell type -> p95 8.2 deg), then again in the v9
block-face classifier (r < 0.95 / r > 1.5), which dropped the cuboid rate to
44 %. Transferring the exact labels from the input tet mesh by nearest-face
lookup restored it to 74 %. Never classify by raw coordinate thresholds when
an exact labelling is available.
