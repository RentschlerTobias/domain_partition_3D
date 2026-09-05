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

---

## Stage 5: block-structure postprocessing (`clean_blocks.py`)

Executing `POSTPROCESSING_PLAN.md`. Three of its assumptions did not survive
contact with the data; all numbers below are measured on v9.

### The documented baseline was not reproducible from the repo

`block_faces.py output/hex3d_algohex/T1_9_hex_v9.ovm` gives **39/82 cuboids
(48 %)**, not the 61/82 (74 %) recorded as the baseline. The committed
`block_faces.physical_of` still classifies by analytic coordinate thresholds
(`abs(r - R_HUB) < 0.03`, `theta < -9`); the nearest-face label transfer that
produced 74 % was never committed. Rebuilt in `clean_blocks.SurfaceLabeller`
(exact point-to-triangle distance over the 8 nearest input faces), which
reproduces the documented figures exactly: `{5: 4, 6: 61, 7: 8, 8: 1, 9: 6,
10: 1, 12: 1}`, 61/82, 48293/56661 cells. Baseline confirmed.

### 1. The duplicate surface labels were indeed a counting error

All 8 blocks that touch one surface twice do so on **opposite** faces, 0 on
adjacent faces. They are legitimate cuboids in that respect, and collapsing
or splitting them would have damaged sound structure. Cause: `shell_hub` is a
single connected shell spanning r = 0.577 … 1.804 -- hub side and shroud side
are one surface, because `tet_prep_v5`'s connected-component split found only
one component (which is also why v9 has 6 surfaces, not 7). A passage block
touching it top and bottom is a cuboid.

A second, independent counting error: interior block faces were labelled by
*sheet id*, so one interface covered by two sheets read as two faces.
Labelling by the neighbouring block instead is the correct base-complex
definition. Classification alone, no geometry touched: **61 -> 64 cuboids**.

### 2. The v9 hex mesh has two internal cavities (previously unrecorded)

The boundary has 3 connected components, not 1: the outer surface (10 822
quads) plus two closed surfaces of 28 and 22 quads enclosing 0.00170 and
0.00145 volume. That is where "99.8 % coverage" goes. Both are structured
boxes -- exactly 8 valence-3 corners, all other vertices valence 4, and for
the 28-quad one 28 of 56 edges fold by ~90 deg -- so they are a 1x2x4 grid of
8 missing cells and (corrected below when the fill actually ran: 1x2x3, not
1x1x5 -- both fit 22 quads and 24 vertices, the quad counts alone do not
decide) 6 missing cells, and every vertex needed to refill them already
exists.

They cause **8 of the 18 non-cuboid blocks** (7046 cells, including the
3395-, 1530- and 1015-cell blocks): the cavity surface chips 1-6 extra faces
off each. Refilling would take the structure to 73/81 cuboids. Not done:
refilling adds cells to the hex mesh, and nothing else in this module changes
the mesh. Reconstructing the lattice is the remaining work -- the
parallelogram rule stalls after 8 vertices because adjacent quads share only
two, and opposite-face matching fails because the interior thickness edges of
a 1-thick box are not on the cavity surface.

### 3. Merge-based collapse does not work on this structure

The plan expected steps 2 and 3 to do the work. Measured:

| operation | result |
|---|---|
| `collapse_small_sheets` | 0 of 3 sliver sheets accepted, each at +3 excess faces |
| `absorb_tiny_blocks` | 0 of 8 accepted; all 40 candidate pairs cost +3 |
| `merge_non_cuboids` | 2 accepted |
| `split_non_cuboid` | 1 accepted |

The +3 is structural, not a threshold. These sheets are bounded by singular
edges, so dropping one merges two blocks but not their side neighbours, and
each of the four side faces of the union stays split in two. Getting past it
requires Gao et al.'s *mesh-level* sheet collapse, which removes a layer of
hexes and merges the side neighbours in the same step -- a geometry change,
deliberately out of scope here.

A first version of the acceptance rule counted non-cuboid *blocks* instead of
excess faces. It looked like a large win -- 82 -> 65 blocks, cuboid share
78 % -> 91 % -- but the absolute cuboid count *fell* (64 -> 59) and it had
built blocks with 17 and 27 faces, because merging a 9-faced and a 7-faced
block into a 12-faced one reduces the count 2 -> 1. The share only rose
because the denominator shrank. Summing excess faces is monotone and does not
have this failure mode.

The plan's sliver threshold is also self-contradictory: "~5 % of the median,
i.e. the 27/54/154-face ones" -- the median sheet is 533 faces, so 5 % is 26
and selects none of them. `sheet_frac` defaults to 0.30.

### Result

| | baseline | after |
|---|---|---|
| blocks | 82 | 81 |
| cuboids | 61 (74 %) | **65 (80 %)** |
| cells in cuboids | 48 293 (85 %) | 48 533 (86 %) |
| blocks < 10 cells | 8 | 8 |
| inverted cells | 2 | 2 (**0 new**) |
| adjacent same-surface faces | -- | 0 |

Of the +4 cuboids, +3 come from fixing the classification and +1 from the one
accepted split. Tiny blocks are unchanged: every candidate absorption made
the structure worse by the cuboid criterion, and the plan's two goals
("cuboid share up", "tiny blocks down") turn out to conflict here -- a 3-cell
block that is a proper cuboid is fine for TFI, just small.

Since the whole cleanup only edits the cut set, no vertex moves and no cell
is added or removed: "0 new inverted cells" and "block boundary stays on the
input surface" (Hausdorff 0.325 before and after, unchanged) hold by
construction, not by tolerance. They are not evidence that the cleanup is
good.

Regression on `cylinder_hex.ovm` (plan criterion 4): 5 blocks in, 5 blocks
out, 100 % cuboids, validator VALID, nothing changed.

```bash
PY=/root/repos/duty/quadmesh/.venv/bin/python
$PY experimentell/hex3d_algohex/clean_blocks.py
$PY experimentell/hex3d_algohex/clean_blocks.py \
      output/hex3d_algohex/cylinder_hex.ovm --no-input-vtk --out /tmp/cyl.vtk
```

---

## Stage 5b: the two mesh-level operations

Both blockers from stage 5 were then removed, on request. Both change the hex
mesh, which stage 5a had ruled out; rereading the plan, that restriction was
too tight -- the plan's own validator checks a Hausdorff distance "so
simplification cannot silently deform the geometry" and offers the Gao
reference tool "if the quality guarantees turn out to matter", neither of
which makes sense for a cut-set-only cleanup.

### Cavity refill (`fill_cavities`)

Filled by **corner peeling**: at a corner of the missing region exactly three
cavity faces meet at one vertex; they fix seven of the cell's eight vertices,
and the eighth is the only vertex joined to at least two of the three face
diagonals. Emit the cell, XOR its six faces into the surface, repeat. The
surface empties exactly when the region is filled.

Two earlier attempts failed and are worth recording: propagating lattice
coordinates by the parallelogram rule stalls after 8 of 30 vertices, because
two adjacent quads share only two; and matching opposite quads finds nothing,
because the through-thickness edges in the middle of a one-cell-thick slab
belong only to missing cells and so lie on neither the cavity surface nor the
mesh.

Result: **14 cells** (8 + 6), min scaled Jacobian **0.848**, volume 0.003145
against the 0.003141 measured as enclosed, face multiplicity only 1 and 2,
boundary components 3 -> 1. The second cavity is a **1x2x3** box, not the
1x1x5 guessed in stage 5 -- both fit 22 quads and 24 vertices, so the counts
alone never decided it; the fold-edge count (20 flat) does.

Block structure: **72 cuboids of 82 (88 %)** before any cleanup, up from 64,
and 98 % of cells in cuboids instead of 86 %.

### Sheet collapse (`collapse_mesh_sheets`, Gao et al. 2017)

145 sheets (edge parallel classes), sizes 316 … 2866 cells, none fully
interior. Screened all 145: **none is self-intersecting and none adds an
inverted cell** -- every collapse keeps min scaled Jacobian at -0.0607, the
two pre-existing bad cells. So the hard constraint was never the binding one
here, contrary to the expectation that welding vertices would be the risk.

Greedy, accepting only strict improvement in (excess faces, block count) with
no new inverted cell and no worse Hausdorff:

| round | sheet | cells | blocks | cuboids | excess | singular edges | tiny | inverted | Hausdorff |
|---|---|---|---|---|---|---|---|---|---|
| start | — | 56 675 | 82 | 72 | 24 | 125 | 8 | 2 | 0.1851 |
| 1 | 7 | 54 360 | **42** | 36 | **10** | 96 | **2** | **1** | **0.1000** |
| 2 | — | no sheet improves, stop |

One collapse halves the block count. It also removes one of the two inverted
cells and moves the boundary *closer* to the input surface, because the quad
carrying the old maximum was in the removed layer. Only 1 of 145 sheets
helped at all; 35 of the 36 sampled in the pre-screen left the structure at
82 blocks.

After the collapse the cut-set steps have nothing left to do (0 merges, 0
splits), which is consistent with stage 5's finding about their reach.

### Where it lands

| | baseline | + classification & cut-set | + cavity fill | + sheet collapse |
|---|---|---|---|---|
| cells | 56 661 | 56 661 | 56 675 | 54 360 |
| blocks | 82 | 81 | 81 | **42** |
| cuboids | 61 (74 %) | 65 (80 %) | **73 (90 %)** | 36 (86 %) |
| cells in cuboids | 85 % | 86 % | 98 % | 98 % |
| blocks < 10 cells | 8 | 8 | 8 | **2** |
| worst block | 12 faces | 13 | 13 | **9** |
| inverted cells | 2 | 2 | 2 | **1** |

The two endpoints are genuinely different trade-offs, not one dominating the
other. The cavity-fill result has the highest cuboid *share* (90 % of 81);
the sheet-collapse result has half as many blocks, a quarter as many tiny
blocks, fewer defective blocks in absolute terms (6 vs 8) and one inverted
cell fewer, at 86 % of 42. For the conforming-division MILP, which scales
with block count, the 42-block structure is the better input.

Cylinder regression holds for both paths: 5 blocks in, 5 out, 100 % cuboids,
and in the collapse path no sheet passes the acceptance test at all.

```bash
PY=/root/repos/duty/quadmesh/.venv/bin/python
$PY experimentell/hex3d_algohex/clean_blocks.py                    # 81 blocks
$PY experimentell/hex3d_algohex/clean_blocks.py --collapse-rounds 5 \
      --out output/hex3d_algohex/deliverable/T1_9_blocks_v9_gao.vtk  # 42, ~30 min
```

---

## Stage 5c: quality repair (Gao's second half) -- the mesh is now VALID

Framing first, because it changes what this step is for: in Gao's pipeline
the smoothing repairs damage done by the collapses. **Our collapses did no
damage** -- none of the 145 sheets added an inverted cell. The only thing
left to repair was the single cell inherited from HexEx.

`untangle` maximises element quality over the cells around the inverted ones
by moving vertices, in two stages.

**Stage 1, interior vertices only.** Attractive because it cannot change
anything else: the domain boundary is untouched to the bit, so the Hausdorff
distance cannot move, and connectivity is untouched, so the block topology
cannot move either. It is not enough. It lifts the worst cell from **-0.0582
to -0.0451** and stops; 4 of that cell's 8 vertices are on the domain
boundary, and the inversion is already present in its boundary quads.

**Stage 2, boundary vertices sliding on the input surface**, penalised by
their squared distance to it. This untangles: **1 inverted cell -> 0, min
scaled Jacobian -0.0451 -> +0.0143**, 8 vertices moved (4 of them on the
boundary), largest displacement 0.027 against a local edge length of ~0.05.
Measured boundary drift **0.000000** and Hausdorff unchanged at 0.09996 -- the
sliding vertices stayed on the surface to 6 decimals, so the domain is not
deformed.

### Two objective bugs, both the same shape as the earlier one

1. The natural objective -- maximise the soft MINIMUM scaled Jacobian --
   is actively harmful. It raised the worst cell to -0.0453 while dragging
   three neighbours below zero: **1 inverted cell became 4**, and the
   headline number ("worst cell improved") said it was working. Replaced by
   a one-sided barrier, sum of max(0, eps - J)^2, which pushes everything
   above eps and leaves good cells alone. Acceptance ranks (number of
   inverted cells, worst value), never the worst value alone.

2. The escalation to stage 2 was keyed on stage 1 *failing*. Stage 1 does
   not fail -- it improves the worst value while leaving the cell inverted --
   so the escalation never fired and the run stopped one step short of the
   repair. Keyed on the inverted-cell count now.

That is three times in this stage that a plausible scalar metric moved the
right way while the structure got worse. The pattern each time: the metric
was an average or an extremum where the thing being minimised is a *count*.

### Final state

```
$PY experimentell/hex3d_algohex/clean_blocks.py --collapse-rounds 5 --untangle
```

| | baseline | + classification & cut-set | + cavity fill | + sheet collapse | + untangle |
|---|---|---|---|---|---|
| cells | 56 661 | 56 661 | 56 675 | 54 360 | 54 360 |
| blocks | 82 | 81 | 81 | **42** | **42** |
| cuboids | 61 (74 %) | 65 (80 %) | 73 (90 %) | 36 (86 %) | 36 (86 %) |
| cells in cuboids | 85 % | 86 % | 98 % | 98 % | 98 % |
| blocks < 10 cells | 8 | 8 | 8 | **2** | **2** |
| worst block | 12 faces | 13 | 13 | **9** | **9** |
| inverted cells | 2 | 2 | 2 | 1 | **0** |
| min scaled Jacobian | −0.0607 | −0.0607 | −0.0607 | −0.0582 | **+0.0143** |
| `HexBlockValidator` | — | INVALID | INVALID | INVALID | **VALID** |

The validator passes for the first time in this project. What remains is 6
non-cuboid blocks (7, 7, 7, 8, 8, 9 faces) and 2 blocks under 10 cells; no
operation implemented here improves them further.

Cylinder regression unaffected: `untangle` finds no inverted cell and returns
immediately, so the 5-block result stands on every path.

---

## Block edges as 1D elements

A solid hex mesh renders as an opaque brick — the block structure is
invisible. The block edges are now written as 2-node line elements into the
`.msh` alongside the hexes (physical tag = curve id, 1-based; the hexes keep
`block_id`), plus a standalone `*_edges.vtk` wireframe.

Definition (`block_edge_curves`): a mesh edge lies on a block edge iff some
block sees it on the border between two of its own faces. Block membership
alone does not find these -- an edge along the corner between two boundary
faces of one block is incident to exactly one block, like every other edge of
that block's surface -- so this needs the face partition and therefore the
surface labels. Curves are split at block corners, so each curve is one edge
of the block topology, which is the entity `dp3d/tmesh.py:743
solve_edge_divisions` will need to assign a division count to.

| | segments | curves | segments per curve (min / median / max) |
|---|---|---|---|
| `T1_9_blocks_v9_clean` (81 blocks) | 3589 | 423 | 1 / 7 / 28 |
| `T1_9_blocks_v9_gao` (42 blocks) | 2628 | 234 | 1 / 10 / 35 |

Verified after writing: every line node is a node of the hex mesh, hex
physical tags cover exactly the expected block count, line tags are
contiguous over the curves.

`--edges-only BLOCKS_VTK` re-derives them from an existing deliverable
(`BlockStructure.from_arrays` rebuilds the face partition from the stored
`block_id`) rather than rerunning the ~30 min pipeline.

---

## Block-edge zigzag: diagnosed, and NOT smoothable

Requested as a direct fix. It is not one, and the measurements say why.

**The zigzag curves are the `shell_hub | shell_blade` label boundary.** Kink
angles over all 234 block-edge curves: median 1.4 deg, p95 11.2, 63 above 30.
Split by curve type:

| type | curves | vertices | median | p95 | > 30 deg |
|---|---|---|---|---|---|
| **hub \| blade** | 14 | 160 | 5.2 | **96.5** | **45** |
| block \| block | — | 576 | 1.9 | 11.0 | 4 |
| surface \| block | — | 1216 | 1.2 | 8.3 | 14 |
| surface \| surface | — | 442 | 0.5 | 3.9 | 0 |

**The hex mesh is not aligned to that feature curve at all**: of all hex
boundary vertices, **zero** lie within 0.010 of the true `shell_hub |
shell_blade` curve of the input triangulation; the staircases run beside it at
a median distance of 0.074 (1–3 cells). So there is no smooth edge chain to
route the label boundary onto, and no nearby geometry to smooth toward.
Smoothing them would only shrink the staircase amplitude while dragging the
surface mesh off the geometry — cosmetics over an alignment defect, and
exactly the class of mistake README "What was learned" 3 warns about.

Of the remaining 18 kinks on curves that *are* smoothable, **16 sit directly
next to a pinned block corner** and 14 are on curves only 3–4 vertices long:
they are short block edges turning at their corner, i.e. block topology, not
roughness. 10 of the 18 are near-perpendicular (60–120 deg), so staircases
too.

Conclusion: all 63 kinks above 30 deg are structural. The other 97 % of the
wireframe is already smooth.

### Two failed smoothing attempts, both caught by their guards

`smooth_block_edges` is implemented, guarded, and is a **no-op on v9** —
which is the correct behaviour, not a missing feature.

1. *Global damped sweep* (alpha = 0.5, all 2234 vertices at once): rejected
   0/30 iterations. The touched cells start at a minimum scaled Jacobian of
   only 0.0986, so one sweep inverts a cell.
2. *Per-vertex with backtracking, accepting "no worse than 80 % of current"*:
   65 626 moves accepted and the result was **worse** — median kink 1.4 -> 0.5
   deg, but p95 11.2 -> 15.0, kinks above 30 deg **63 -> 93**, and the touched
   cells fell from 0.0986 to 0.0047. The relative-to-current floor compounds
   over sweeps. Fixed by measuring the floor against the ORIGINAL mesh and by
   accepting a sweep only on (count > 30 deg, p95, median) lexicographically.
   With that, 0 sweeps are accepted: a Laplacian sweep removes none of the 18.

That is the fifth time in this project that a central-tendency metric moved
the right way while the structure got worse. The deliverable was restored to
the pre-smoothing state and verified byte-identical in its point array.

The defect is a frame-field alignment failure, not a postprocessing gap, and
is carried into `FRAMEFIELD_PLAN.md`.

## Singular graph census (input to `FRAMEFIELD_PLAN.md`)

```
125 interior singular edges     valence 3: 22   valence 5: 103
130 singular vertices, 10 nodes, all degree 1
5 singular arcs, lengths 1.006 … 1.226, all ending on the boundary
0 arcs below length 0.10
```

| # | edges | valence | length | r start → end | z start → end |
|---|---|---|---|---|---|
| 1 | 27 | 5 | 1.226 | 1.80 → 0.59 | 1.55 → 1.58 |
| 2 | 27 | 5 | 1.208 | 1.80 → 0.60 | 1.05 → 1.04 |
| 3 | 27 | 5 | 1.202 | 0.60 → 1.80 | 1.11 → 1.12 |
| 4 | 22 | **3** | 1.028 | 1.00 → 1.06 | 1.49 → 2.50 |
| 5 | 22 | 5 | 1.006 | 0.83 → 0.87 | 1.51 → 2.50 |

The graph is far cleaner than the "too many singularities" hypothesis
assumes: no zipper nodes, no complex nodes, no short arcs, nothing for the
standard simplification operators to act on. Arcs 1–3 are the 3D expression
of the hub/shroud topology difference that motivated going 3D at all and are
probably irreducible. Arcs 4 and 5 are a valence-3 / valence-5 pair running
parallel to the outlet — the 3D analogue of the +1/-1 pair the 2D pipeline
already annihilates, and the one real candidate.

---

## Stage 6: re-attaching the removed parts (`reattach.py`)

Done before the frame-field stage, on request. The two removed parts needed
completely different treatment.

### Blade O-grid — reused verbatim, and it was already a block structure

44 800 hexes in the source MSH (geom regions 2–6), min scaled Jacobian
**0.5217**, 0 inverted. Its base complex has **0 singular edges**, so as a
mesh it is one regular ring — but the five geom regions turn out to be
**exactly five topological cuboids** (16 000 / 6 400 / 16 000 / 3 200 / 3 200
cells, 6 faces each). The block decomposition was already in the file.
Nothing to generate.

### Hub/shroud boundary layer — regenerated by extrusion

The 121 080 cells are triangular PRISMS (all geom region 1), so nothing there
is reusable as hex blocks. Rebuilt by extruding each AlgoHex block face
carrying the `shell_hub` label out to its wall: 29 such faces, 4 134 quads.

Which wall a face belongs to is decided by nearest-face lookup against the
MSH's own tagged wall triangulations, never by a radius threshold —
`shell_hub` is one connected shell spanning r = 0.577 … 1.804 and wrapping
both sides, so a threshold would cut through it. Result: 11 hub, 15 shroud.

Three defects found and handled, each measured:

**1. An indexing bug, caught by the inverted-cell count.** The point offset
per face used `len(newP)`, which counts appended *arrays* (0, 2, 4, …) rather
than points, so every face after the first referenced wrong nodes:
**3 606 of 4 134 cells inverted**, layer thickness ranging to 1.85 (the full
channel width). After the fix: 42 inverted, thickness 0.087 … 0.175.

**2. Nearest-point projection warps oblique patches.** Both walls are
**exact cylinders** — measured on the MSH's own wall triangulations, radius
standard deviation **0.000000** for both (hub 0.50000, shroud 1.90000) — so
the offset is done radially, which is the correct construction here and not
an approximation. Projecting each vertex to its own closest wall point left
42 inverted cells; radial offset leaves **5**. The wall *assignment* is still
by nearest-face lookup.

**3. Three `shell_hub` patches are not boundary-layer faces at all** and are
skipped, with a message: a 46-quad patch at r = 0.591 (0.083 … 0.417 from the
wall), an 82-quad patch at r = 1.718 (0.098 … 0.856), and a 9-quad patch at
**r = 0.977 — mid-passage, 0.48 from the hub and 0.92 from the shroud**.
Extruding that one would sweep a block across half the channel.

The remaining **5 inverted cells are inherited, not created**: their base
quads are slivers in the AlgoHex boundary itself, short edge 0.003 … 0.005
against an extrusion thickness of 0.14 … 0.175, i.e. 30–50 : 1.

### Assembly

| part | cells | blocks | inverted |
|---|---|---|---|
| AlgoHex core | 54 360 | 42 | 0 |
| blade O-grid | 44 800 | 5 | 0 |
| boundary layer | 3 997 | 26 | 5 |
| **total** | **103 157** | **73** | **5** |

Mean scaled Jacobian 0.9108.

**Conformity, measured rather than assumed.** The parts do not share nodes,
by design — the reduced-domain strategy relies on TFI regenerating block
interiors, so the parts only have to agree at block level. `interface_gap()`
checks the assumption: at the `shell_blade` interface the 1 307 AlgoHex quads
sit at median **0.0418**, p95 0.0654, max 0.1144 from the nearest O-grid
boundary quad — about one cell, so the removed and re-attached volumes are
the same region.

### Block edges are computed PER PART

Computing them on the assembled mesh over-segments badly — 1 108 curves
against 234 for the core alone — because every boundary-layer slab is its own
closed surface, so a feature-angle segmentation of the assembly invents patch
borders at every non-conforming interface. Each part has an exact labelling
of its own, so they are computed separately and concatenated:

| part | segments | curves |
|---|---|---|
| core | 2 628 | 234 |
| O-grid | 948 | 40 |
| boundary layer | 2 874 | 404 |
| **total** | **6 450** | **678** |

Written into `T1_9_blocks_v9_full.msh` as 1D line elements (physical tag =
curve id), plus `T1_9_blocks_v9_full_edges.vtk` as a standalone wireframe and
`T1_9_blocks_v9_full_part.vtk` with `part` (0 = core, 1 = O-grid, 2 = layer).

### Consequence for the zigzag and for the frame-field stage

`shell_blade` is no longer a domain boundary now that the O-grid is attached.
The `shell_hub | shell_blade` label boundary — which carried 45 of the 63
block-edge kinks above 30° and which no smoothing could fix — is gone from
the assembled structure. That was the coupling predicted in
`FRAMEFIELD_PLAN.md` §7, and it means step 1 of that plan (feature alignment)
should be re-measured before it is run.

Also relevant to that plan: the two singular arcs the user identified run
from **`shell_blade` to the outlet** (valence 3 and valence 5, 22 edges
each) — verified by looking up the surface label at each arc endpoint. Both
end on the outlet plane; neither touches the inlet. With the O-grid attached,
their `shell_blade` ends are interior.
