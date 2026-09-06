# Mesh quality criteria for CFD, measured on our own mesh

What a CFD solver cares about, why, and where **our** current mesh stands.
Every number below is measured by `mesh_quality.py`, not quoted from a table.

```bash
PY=/root/repos/duty/quadmesh/.venv/bin/python
$PY experimentell/hex3d_algohex/mesh_quality.py            # full domain
$PY experimentell/hex3d_algohex/mesh_quality.py MESH.vtk --no-render
```

Output in `output/hex3d_algohex/quality/`: one `.vtk` per metric (threshold
it in ParaView), plus a `.png` and an interactive `.html`. In every figure
the mesh is drawn translucent and the cells that **violate** the limit are
drawn opaque in red, so a figure shows *where* a metric fails, not only that
it does.

The limits are OpenFOAM's own defaults from
`etc/caseDicts/meshQualityDict`, because that is the solver this project
feeds. They are guidance, not physics: how much a solver tolerates depends on
the scheme and the flow. Simple incompressible CFD is far less sensitive to a
few bad cells than compressible flow with strong density gradients.

---

## Where we stand

Measured on `T1_9_blocks_v11_full.vtk` — the full domain, 103 148 cells:

| metric | min | median | max | limit | violations |
|---|---|---|---|---|---|
| scaled Jacobian | **0.1150** | 0.9712 | 0.9997 | > 0 | **none** |
| non-orthogonality [°] | 0.35 | 12.32 | **68.41** | ≤ 65 | 4 cells |
| skewness | 0.0017 | 0.0160 | **0.4594** | ≤ 4 | none |
| aspect ratio | 1.02 | 1.96 | **1528** | ≤ 1000 | 2 cells |
| volume ratio | **0.1910** | 0.7271 | 0.9988 | ≥ 0.01 | none |
| face weight | **0.1884** | 0.4237 | 0.4999 | ≥ 0.05 | none |
| face flatness | **0.0000** | 0.9995 | 1.0000 | ≥ 0.8 | 3 cells |

**9 cells of 103 148 — 0.009 % — violate an OpenFOAM default.** Split by part:

| part | cells | min scaled Jac. | max non-ortho | max skew | max aspect |
|---|---|---|---|---|---|
| AlgoHex core | 54 460 | 0.1524 | 65.38° | 0.459 | 1528 |
| blade O-grid | 44 800 | 0.5217 | 58.31° | 0.192 | 177 |
| hub/shroud BL | 3 888 | 0.1150 | 68.41° | 0.245 | 13.8 |

All 2 aspect-ratio and all 3 flatness violations are in the AlgoHex core; the
4 non-orthogonality violations split 2 core / 2 boundary layer. The
re-attached O-grid — the part we did not generate — is the cleanest of the
three on every metric.

---

## The metrics

### Scaled Jacobian — the only hard one

The minimum over the 8 corners of the determinant of the corner frame,
normalised by the three edge lengths. Positive everywhere means the cell is
not turned inside out.

$$J_c = \frac{(\mathbf{u} \times \mathbf{v}) \cdot \mathbf{w}}{|\mathbf{u}||\mathbf{v}||\mathbf{w}|}, \qquad
\text{scaled Jacobian} = \min_{c \in \text{8 corners}} J_c$$

A negative value is fatal: OpenFOAM aborts on a negative cell volume. Note a
sheared cell can have positive *volume* and still an inverted corner, which
is why the corner test is the right one — `hexa_interpolation.py` uses corner
Jacobians for exactly this reason.

**Ours: 0.1150 minimum, nothing inverted.** File
`quality_scaled_jacobian.vtk` (field ×1000).

### Non-orthogonality — the one that matters most for hex meshes

The angle between the vector joining two neighbouring cell centres and the
normal of the face they share. It enters the diffusion term directly: a
non-orthogonal face needs an explicit correction, which is lagged and
therefore both diffusive and destabilising.

$$\theta = \arccos\left(\frac{\mathbf{d} \cdot \mathbf{S}_f}{|\mathbf{d}||\mathbf{S}_f|}\right)$$

OpenFOAM's `maxNonOrtho` default is **65°**; above ~70° you need
`nNonOrthogonalCorrectors` in `fvSolution`, which costs iterations and only
compensates numerically.

For *structured* hex meshes this is the crucial metric — more so than
skewness, which dominates for tetrahedra. That is the whole argument for a
block-structured mesh in the first place.

**Ours: max 68.41°, 4 cells above 65°.** The figure shows them along the
blade O-grid, where the wall-normal direction turns fastest.

### Skewness — the offset between face centre and the centre-line crossing

Where the line joining the two cell centres pierces the shared face, measured
against the face's own centre, normalised by the centre distance. It biases
the interpolation of face values.

$$\text{skew} = \frac{|\mathbf{x}_{\text{pierce}} - \mathbf{x}_f|}{|\mathbf{d}|}$$

OpenFOAM: `maxInternalSkewness 4`, `maxBoundarySkewness 20`.

**Ours: max 0.4594 — an order of magnitude inside the limit.** This is the
expected result for a block-structured mesh and is the metric where such a
mesh is naturally strongest.

### Aspect ratio — cheap where the flow allows it

Longest over shortest cell edge. Ideally 1, but a boundary layer *needs*
large values: resolving a wall gradient means thin cells, and stretching them
along the wall is exactly how you avoid paying for it in cell count. High
aspect ratio is only harmful where the gradient is large in the long
direction too.

`checkMesh` warns above **1000**; boundary layers routinely and legitimately
exceed it.

**Ours: max 1528, 2 cells.** Both in the AlgoHex core, not in a boundary
layer — that is where it is a defect rather than a design choice.

### Volume ratio and face weight — smoothness of the size transition

`vol_ratio` is the smaller over the larger volume of two neighbours
(OpenFOAM `minVolRatio 0.01`); `face_weight` is how centrally the shared face
sits between the two centres (`minFaceWeight 0.05`, ideal 0.5). Together they
measure how abruptly cell size changes. An abrupt jump adds numerical
diffusion; the usual guidance is to keep growth below 15–20 % per cell, and
1.0–1.2 from a wall outward.

**Ours: 0.1910 and 0.1884 minimum, both comfortably inside.** Median face
weight 0.42 against an ideal 0.5.

### Face flatness — a quad in 3D is generally not planar

A hexahedral face has four nodes that need not be coplanar. Splitting it into
two triangles two different ways gives two different normals; the cosine
between them measures the warp. A badly warped face makes the flux through it
ambiguous.

**Ours: 3 cells at flatness 0, all in the AlgoHex core.** Those are the
sliver quads already known from the boundary-layer extrusion (short edge
0.003 against 0.05 neighbours).

---

## What is deliberately *not* in here

**y+.** The wall-normal first-cell height that a turbulence model needs —
y⁺ ≈ 1 for wall-resolved, 30–300 for wall functions. It cannot be evaluated
on a mesh alone: it depends on the flow (velocity, viscosity, wall shear), so
it is a *sizing input* to the mesh, not a property of it. It belongs in the
TFI stage, where `edge_fractions` (tanh clustering, `dp3d/tmesh.py:811`)
places the first cell — and it needs a target flow condition to be chosen.

**Grid convergence.** Ultimately no metric substitutes for solving on
successively finer meshes. The metrics tell you a mesh will not break the
solver; only a grid study tells you the answer is right.

---

## How this connects to TFI

`HexBlockValidator` checks the **block structure** — is each block a
topological cuboid. `mesh_quality.py` checks the **mesh** a solver integrates
over. They are different questions, and the second is what TFI decides: the
block structure fixes the boundary curves, and TFI fills the interior.

That makes this the gate for the TFI stage, and the reason it was built
before it. The numbers above are the *before*: whatever TFI produces has to
be compared against the same seven metrics on the same mesh, and the two
places our current mesh is weakest — non-orthogonality along the blade
O-grid, and 3 sliver cells in the core — are exactly where a transfinite fill
either helps or does not.

Sources: [OpenFOAM meshQualityDict](https://github.com/OpenFOAM/OpenFOAM-7/blob/master/etc/caseDicts/mesh/generation/meshQualityDict) ·
[OpenFOAM mesh quality guide](https://www.openfoam.com/documentation/guides/latest/doc/guide-meshing-snappyhexmesh-meshquality.html) ·
[SimScale, Mesh Quality](https://www.simscale.com/docs/simulation-setup/meshing/mesh-quality/) ·
[Resolved Analytics, Mesh Quality and CFD Solution Accuracy](https://www.resolvedanalytics.com/cfd-in-practice/what-is-effect-of-cfd-mesh-on-solution-accuracy) ·
[wolfdynamics, meshing preliminaries and quality assessment](https://www.wolfdynamics.com/wiki/meshing_preliminaries_and_quality_assessment.pdf)
