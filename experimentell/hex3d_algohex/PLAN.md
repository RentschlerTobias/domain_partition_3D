# 3D octahedral frame field → IGM → libHexEx hex blocks for T1_9

## Context

The current `dp3d` pipeline builds a 2D cross field *per surface* (hub, shroud
separately) in the unwrapped `(s,t)` domain, partitions each into quad
blocks, and — for the 3D lift — has to **morph** the hub block skeleton onto
the shroud (`experimentell/3d_extrapolation/hexa_interpolation.py`). That
morph exists only because the hub and shroud 2D cross fields are
**topologically different** (hub: four idx=−1 singularities; shroud: four
idx=+1 — a real, physical effect of blade twist, not a solver bug — see
`hexa_interpolation.py:4-8`). A per-surface 2D field with a ruled/extruded
lift cannot represent a genuinely 3D singularity graph; the morph is a
workaround, not a generalization to "beliebige Volumen" (arbitrary volumes).

The requested approach is the one used by the pipeline in
[arXiv:2603.12820 (NeurFrame)](https://arxiv.org/html/2603.12820v1) and by
the broader frame-field-hexmeshing literature it builds on: a genuinely 3D
**octahedral (24-symmetry) frame field** on a tet mesh of the volume, an
**integer-grid map (IGM)** aligned to it (CubeCover-style), and
**HexEx/libHexEx** ([Lyon et al. 2016](https://doi.org/10.1145/2897824.2925976))
to robustly extract a hex mesh from the (imperfect) IGM. This is real,
working, cited machinery — not something to hand-roll. The honest caveat,
stated by the papers themselves: smooth octahedral fields frequently have
singularity graphs that do **not** admit a valid IGM
([Liu & Bommes 2023, "Locally Meshable Frame Fields"](https://www.algohex.eu/publications/locally-meshable-frame-fields/)
report the HexMe benchmark success rate at 58%, up from 2% for naive
fields). So this can fail on a genuinely hard geometry — T1_9's twisted
blade passage is exactly the kind of case that stresses this. Budget for
iteration on feature tagging / field constraints, not a first-try success.

**Chosen implementation**: [AlgoHex](https://github.com/cgg-bern/AlgoHex)
(cgg-bern, AGPLv3) bundles all four pipeline stages behind one CLI
(`HexMeshing`): frame field init → locally-meshable frame field → seamless
parametrization + quantization (via
[QGP3D](https://github.com/HendrikBrueckler/QGP3D), which internally builds
a **3D Motorcycle Complex** — i.e. it already computes a block
decomposition as part of quantization) → **libHexEx** extraction. Straight
through, as requested, rather than reimplementing octahedral-field
optimization or CubeCover from scratch.

**Scope for this plan** (per your answers): T1_9 only, single passage (no
periodicity handling — treat the two pitch-side surfaces as ordinary
feature surfaces for now), coarse blocks obtained by extracting the base
complex from AlgoHex's fine hex output (Gao et al. 2015/2017), with a
one-time check of whether QGP3D's internal Motorcycle Complex can be
exported directly (would make the base-complex step largely unnecessary).

**License note** (not a blocker, worth knowing): AlgoHex is AGPLv3,
libHexEx/QGP3D are GPLv3. Fine for internal research use in this repo;
matters only if you ever distribute a binary or run it as a network
service — flagging so it's a conscious choice later, not a surprise.

## Environment check (done during planning)

- Docker is available and the daemon is reachable (`docker info` works).
  AlgoHex ships an official `Dockerfile` (Debian trixie base, builds
  Ipopt 3.14.19 + Bonmin from COIN-OR via `coinbrew`, then AlgoHex via
  CMake+Ninja) — this is the CI-validated build path.
- **No cmake/ninja on the host** (Arch, pacman has neither pre-fetched).
  Building AlgoHex's C++ stack (OpenVolumeMesh, libHexEx, CoMISo, TinyAD,
  QGP3D, Eigen, IPOPT, Bonmin) from source directly on the host is a real
  risk (IPOPT/Bonmin builds are notoriously fragile). **Use the Docker
  image**, not a host build — avoids that risk entirely and matches what
  upstream actually tests.
- venv at `/root/repos/duty/quadmesh/.venv` has numpy/scipy/torch/meshio/
  gmsh/networkx/trimesh already — sufficient for all the Python-side work
  (tet prep, VTK conversion, base-complex extraction, 3D TFI).
- AlgoHex's exact CLI flags (`HexMeshing -h`) and the precise OVM/VTK
  feature-tag property format are **not documented online** in enough
  detail to plan against blindly — Stage 0 includes actually building the
  image and running `-h` / inspecting `demo/HexMeshing/cylinder.ovm` to
  learn this empirically before writing the tet-prep converter.

## T1_9 geometry (verified by parsing the MSH directly)

`data/T1_9/T1_9_ru_gridGmsh.msh` (Gmsh 2.2 ASCII), boundary tagged by
`$PhysicalNames` + per-element geometrical-entity tag (the *second* tag,
which is what `dp3d/extraction.py` already keys off for hub=1/shroud=2):

- blade: 5 patches (geom 1-5), hub: 6 patches (geom 6-11), shroud: 6
  patches (geom 12-17), inlet (geom 18, z=0 plane), outlet (geom 19,
  z=2.5 plane).
- hub cylinder r=0.500, shroud cylinder r=1.900, one blade passage spans
  θ ≈ −68.6°…+83.1° (single-passage sector, not the full 360°).
- Volume is **hybrid**: 89 547 tets, 44 800 hexes (O-grid core),
  121 080 prisms, 1 792 pyramids — this is NOT a pure tet mesh, so it
  cannot be fed to AlgoHex as-is.

## Implementation stages

### Branch
`feat/algohex-3d-frame-field` off `master`.

### Stage 0 — build AlgoHex (Docker)
```
git clone --recursive https://github.com/cgg-bern/AlgoHex.git external/algohex-src   # NOT committed, add to .gitignore
docker build -t algohex external/algohex-src
docker run --rm algohex HexMeshing -h                     # discover real CLI flags
docker run --rm -v ...:/data algohex HexMeshing -i /data/demo/HexMeshing/cylinder.ovm -o /tmp/out.ovm  # sanity run on upstream demo
```
Inspect the demo `.ovm`/feature-tag encoding (OpenVolumeMesh property
names for feature faces/edges/vertices — QGP3D's README shows the
programmatic API marks feature faces/edges/vertices explicitly, so the
file format almost certainly has a matching property block) to nail down
exactly what Stage 1's converter must emit. Also check at this point
whether `HexMeshing` (or a QGP3D debug flag) can dump the intermediate
Motorcycle-Complex block decomposition directly — if so, Stage 4's
base-complex extraction shrinks to "read this file" instead of
reimplementing Gao et al.

### Stage 1 — tet mesh + feature tags for T1_9
New module `experimentell/hex3d_algohex/tet_prep.py`:
- Extend `dp3d/extraction.py`'s MSH parser (reuse `parse_msh`,
  `surface_triangles`) to pull **all** boundary patches (blade, hub,
  shroud, inlet, outlet — geom 1-19), not just hub/shroud, keeping the
  per-patch geom id as a face group.
- Tetrahedralize the enclosed volume from this boundary via gmsh's 3D
  Delaunay (gmsh already a dependency; this matches how the HexMe
  benchmark itself was generated, so it's a well-trodden input shape for
  AlgoHex) — do **not** attempt to split the existing hybrid
  hex/prism/pyramid core into tets; a fresh boundary-conforming Delaunay
  tet mesh gives better-conditioned tets for the frame-field optimizer.
- Tag feature **surfaces**: hub, shroud, inlet, outlet, blade, and the two
  periodic-pitch faces (all boundary-aligned, i.e. field should be
  tangent). Tag feature **curves**: boundary-patch intersections (e.g.
  blade↔hub, blade↔shroud, hub↔inlet, hub↔outlet) as sharp feature edges —
  analogous to how `dp3d/partition_surface.py` already treats the
  hub/shroud seam and blade wall as hard boundaries in 2D.
- Write `data/T1_9/T1_9_tet.vtk` (VTK v2.0 ASCII, matching HexMe's format)
  with feature tags as cell/field data per Stage 0's findings.

### Stage 2 — run HexMeshing
`experimentell/hex3d_algohex/run_algohex.py`: shells out to the Docker
image (`docker run -v <repo>/data/T1_9:/data algohex HexMeshing -i
/data/T1_9_tet.vtk -o /data/T1_9_hex.ovm [flags from Stage 0]`), captures
stdout/stderr and the gate messages AlgoHex itself prints (it already
warns when the output is an incomplete hex mesh). Convert `.ovm` → `.vtu`
for inspection with `meshio` (or the OVM Python bindings if present in the
image) and reuse the existing plot style from
`experimentell/3d_extrapolation/hexa_interpolation.py:plot_3d`.

### Stage 3 — validation gates
Reuse the exact validity check already implemented in
`hexa_interpolation.py:224-265` (`hexa_cell_volumes`,
`hexa_corner_jacobians` — 8-corner-Jacobian trilinear hex validity), applied
to AlgoHex's *arbitrary-topology* output rather than the ruled grid it was
written for (generalize the `(nv,nu,nw,3)` grid-indexing assumption to a
plain hexahedron→8-node array). Gates: 0 sign-flipped corner Jacobians;
all 19 boundary feature surfaces still present on the extracted mesh
boundary (HexMe found feature-surface preservation only 34.6% on average
across methods — this is the gate most likely to fail first and need
iteration on Stage 1's feature tagging or AlgoHex's boundary-alignment
weight flags).

### Stage 4 — coarse blocks (base complex)
New `experimentell/hex3d_algohex/base_complex.py`. Per
[Gao et al. 2015/2017](https://dl.acm.org/doi/10.1145/3130800.3130848):
start from all quad facets incident to a **singular edge** (interior edge
with hex-valence ≠ 4, or boundary edge with quad-valence ≠ 2 / a sharp
feature angle) and iteratively expand through the **opposite** facet
across each **regular** edge until termination — this partitions the hex
mesh's faces into base-complex sheets; the hex groups bounded by sheets
are the coarse blocks. Implement with `networkx` (already a dependency)
over the hex mesh's face-adjacency graph.

Output the block structure in the **same shape** `dp3d/tmesh.py` already
consumes downstream (`result["blocks"]`: cycle, corners, side_chains) so
the existing edge-conforming machinery can be reused as-is rather than
reimplemented for 3D:
- `solve_edge_divisions` (`dp3d/tmesh.py:743`, scipy MILP, opposite-side
  cell-count conformity) — generalizes directly to 3D blocks (12 edges,
  3 independent opposite-pairs per block instead of 2).
- `edge_fractions` (`dp3d/tmesh.py:811`, tanh boundary-layer clustering).

### Stage 5 — 3D TFI fill + export
New `experimentell/hex3d_algohex/tfi3d.py`: trilinear Coons blend per
block (8 corners + 12 edges + 6 faces, straightforward 3D generalization
of `dp3d/tmesh.py:_coons`), reusing the `w`-blend / cylindrical-vs-
Cartesian lessons already learned in `hexa_interpolation.py:199-221`
("interpolation in cylinder coords, not Cartesian, to avoid shear/flip on
twisted blade passages" — likely applies here too if any block spans the
blade twist). Export via the existing `export_vtk` pattern
(`hexa_interpolation.py:268-287`). Gates: 0 inverted TFI cells
(`hexa_corner_jacobians` again), matching grid dimensions across shared
block faces.

## New files
```
experimentell/hex3d_algohex/
  tet_prep.py       # Stage 1: full-boundary extraction + gmsh tet + feature tags
  run_algohex.py    # Stage 2: docker wrapper around HexMeshing
  base_complex.py   # Stage 4: singular-edge sheet extraction -> blocks
  tfi3d.py          # Stage 5: trilinear Coons + export
  README.md         # build instructions, discovered CLI flags, known risks
```
`external/` added to `.gitignore` (AlgoHex source + Docker build context
— large, AGPL, not vendored into this repo).

## Verification
1. `docker run --rm algohex HexMeshing -h` succeeds and the upstream demo
   (`cylinder.ovm`) round-trips to a valid hex mesh — confirms the Docker
   build itself is sound before touching T1_9.
2. Stage 1 output: `T1_9_tet.vtk` boundary triangle count matches the sum
   of the 19 tagged patches' triangle counts from the MSH (no dropped
   boundary); volume closed (watertight boundary check via `trimesh` or
   `meshio`, already available).
3. Stage 2/3: run the full T1_9 case through AlgoHex, report the gate
   numbers (inverted cells, feature-surface preservation) the same way
   `hexa_interpolation.py` already prints `[hexa] GATE ...` lines.
4. Stage 4/5: block count + edge-conformity report (reuse
   `check_edge_conformity`, `dp3d/tmesh.py:1109`, generalized to 3D), plot
   the coarse block wireframe (style of `hexa_interpolation.py:plot_3d`)
   and one interior TFI slice for visual sanity.
5. If Stage 2/3 fails to produce a valid IGM on the full T1_9 passage
   (plausible per the literature's own success-rate numbers), fall back to
   running Stage 0-3 on a simplified sub-volume first (e.g. hub-only thin
   shell without the full blade twist) to isolate whether the failure is
   geometry-general or specific to the twisted passage — report findings
   rather than silently giving up.
