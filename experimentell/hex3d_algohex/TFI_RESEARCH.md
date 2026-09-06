# Transfinite interpolation for CFD meshes — research notes

Background for the TFI stage. The question this had to answer: what exactly
do we implement, what is known to go wrong, and what does the state of the
art add on top of the textbook formula.

---

## 1. What we have and what is missing

`dp3d/tmesh.py` already contains a complete TFI pipeline — but a **2D** one:

```python
def _coons(S, N, W, E):
    """S: c0->c1, N: c3->c2 (both len n_u); W: c0->c3, E: c1->c2 (len n_v).
    Returns grid (n_v, n_u, 2)."""
```

Four boundary curves, a 2D grid. `tfi_fill` builds one Coons patch per quad
block. And `experimentell/3d_extrapolation/hexa_interpolation.py` gets its
hexahedra by connecting *two* such 2D grids (hub and shroud) — the ruled lift
whose limitations are the reason this branch exists (README, motivation).

**A genuine 3D transfinite interpolation does not exist in the repo.** What
is directly reusable:

| piece | file | reuse |
|---|---|---|
| tanh wall clustering | `dp3d/tmesh.py:811 edge_fractions` | as is, 1D |
| conforming divisions (MILP) | `dp3d/tmesh.py:743 solve_edge_divisions` | idea yes, constraints no — they encode "opposite sides of a *quad*"; a hex block has 3 families of 4 parallel edges |
| Coons patch | `dp3d/tmesh.py:948 _coons` | needs the 3D analogue |
| Thomas–Middlecoff smoothing | `dp3d/tmesh.py:972 _tm_smooth` | 2D elliptic; for 3D our `clean_blocks.untangle` is the closer starting point |

---

## 2. The formula

Transfinite interpolation is due to **Gordon & Hall (1973)**, extending the
**Coons patch** (1967). The name is theirs: the interpolant matches the
boundary data at a *non-denumerable* set of points, not at finitely many
nodes as classical interpolation does.

The construction is a **Boolean sum of projection operators**. With
univariate projectors $P_1, P_2, P_3$ in the three parameter directions,

$$P_1 \oplus P_2 \oplus P_3 = P_1 + P_2 + P_3 - P_1P_2 - P_1P_3 - P_2P_3 + P_1P_2P_3$$

For a hex block parametrised over $(u,v,w) \in [0,1]^3$, with the six
boundary faces $S$, the twelve edge curves $C$ and the eight corners $X$,
linear blending gives

$$
\mathbf{x}(u,v,w) = \underbrace{\sum_{i=1}^{3} F_i}_{\text{6 faces}}
- \underbrace{\sum_{i<j} E_{ij}}_{\text{12 edges}}
+ \underbrace{V}_{\text{8 corners}}
$$

with, for example,

$$F_1 = (1-u)\,S_{u=0}(v,w) + u\,S_{u=1}(v,w)$$
$$E_{12} = \sum_{a,b \in \{0,1\}} \beta_a(u)\beta_b(v)\, C_{ab}(w), \qquad \beta_0(t)=1-t,\ \beta_1(t)=t$$
$$V = \sum_{a,b,c \in \{0,1\}} \beta_a(u)\beta_b(v)\beta_c(w)\, X_{abc}$$

The subtraction of the edge terms and the re-addition of the corners is not
cosmetic: without it every corner would be counted three times and every edge
twice, and the map would not reproduce its own boundary.

Perronnet's [transfinite interpolation catalogue](https://www.ljll.fr/perronnet/transfini/transfini.html)
gives explicit forms for the square and the **cube** (Coons, Gordon, Hall),
plus the triangle, tetrahedron and pentahedron — the most direct source for
the 3D formula.

---

## 3. What goes wrong, and what to do about it

### Folding

Gordon & Hall define the transfinite element as an **invertible** map. Folding
is exactly the failure of that hypothesis: the Boolean sum's correction terms
overshoot when boundary curvature is strong relative to block size, and
$\det J \le 0$ appears in the interior *even though all six boundary faces are
perfectly valid*.

This is the single most important thing to expect. It is why the TFI stage
needs the quality gate (`mesh_quality.py`) rather than a visual check, and
why an untangling step has to be part of the plan rather than a contingency.

### Parametrisation mismatch

TFI only produces a sensible grid **if the bounding surfaces meet at the
edges and are parametrised in the same direction** — otherwise grid lines
cross. In our case this is not a detail: the block faces come out of the
base complex as *sets of mesh quads*, and turning each into a consistently
oriented $(s,t)$ patch is a real step, not bookkeeping. The conforming
division MILP exists precisely so opposite faces of a block agree on their
counts.

### The standard remedy: TFI first, elliptic smoothing after

The near-universal workflow is algebraic TFI for the initial grid, then an
elliptic (Winslow / Thompson–Thames–Mastin) solve to smooth it. Elliptic
systems suppress boundary singularities and stop gradient discontinuities
propagating inward, which is what buys orthogonality; the cost is solving a
PDE per block.

Orthogonality and spacing are steered through **control functions** — Steger
& Sorenson, or Hilgenstock's variant — constructed on the boundaries and
propagated inward, classically by transfinite Lagrangian interpolation of the
forcing terms themselves.

### Or: improve the algebraic stage instead

**Allen (2008), "Towards automatic structured multiblock mesh generation
using improved transfinite interpolation"** modifies TFI to include
orthogonality and spacing control directly, plus an aspect-ratio-based
smoothing that removes grid crossover — avoiding the elliptic solve
altogether. For a first implementation this is the more attractive route:
one pass, no PDE, and the crossover fix targets exactly the failure mode
above.

---

## 4. What turbomachinery practice does

**NASA's TIGGERC** is the closest published match to our plan: block boundary
points distributed by a **hyperbolic tangent or algebraic** stretching, block
interiors filled by **transfinite interpolation**, producing H-, C-, I- or
O-grid topologies. That is exactly `edge_fractions` + 3D TFI.

On topology, the literature splits blade meshing into **passage-centred**
(periodic boundaries cut through leading and trailing edge) and
**blade-centred** (periodics away from the blade). Passage-centred gives less
skewed cells for highly cambered blades in narrow passages but is more prone
to solver instability, because periodic boundaries then sit in high-gradient
regions. Our structure is effectively blade-centred: an O-grid around the
blade inside an H-type passage, which is the standard combination — the H
grid carries pitchwise periodicity, the O grid resolves the blade boundary
layer.

For the near-wall distribution itself, **Vinokur (1983), "On one-dimensional
stretching functions for finite-difference calculations"** is the definitive
reference for two-sided tanh/sinh stretching with prescribed end spacings.
`dp3d`'s `edge_fractions` already implements a tanh clustering with a bisection
solve for the stretching parameter, so this is a check of our formulation
rather than new work.

**Sauer et al. (2023), DLR**, is already cited in this repo (`dp3d/tmesh.py`
follows its section 2.7 for hanging seam junctions). Their optimisation-based
multiblock generator takes a topological block partition plus block sizes and
distance constraints — and explicitly notes that **initialisation is done by
basic algebraic techniques: orthogonal offsets, equidistant lines, and
transfinite interpolations.** In other words, TFI is the initialiser and the
optimisation is the quality step. That is a useful framing for where to stop:
if TFI plus untangling gets us inside the quality limits, the optimisation
layer is optional.

---

## 5. Implementation plan

Ordered so the risky part is answered first.

**Step 1 — the core map, on one block.** Trilinear Gordon–Hall over the six
faces, with the block's own boundary quads as the face data. Immediately
gated by `mesh_quality.py`: scaled Jacobian > 0 everywhere, and the seven
metrics compared against the pre-TFI mesh. If a single block folds, that is
the answer to the feasibility question and it arrives on day one.

**Step 2 — consistent face parametrisation.** Turn each block face's quad set
into an $(n_s \times n_t)$ ordered grid, with the four blocks sharing an edge
agreeing on its sampling. This is where "the bounding surfaces must be
parametrised the same way" becomes concrete work.

**Step 3 — conforming divisions.** Extend the MILP from quad conformity
(2 pairs of opposite sides) to hex conformity (3 families of 4 parallel
edges), reusing `solve_edge_divisions`' structure.

**Step 4 — wall clustering.** `edge_fractions` on the wall-normal edge family,
with the first cell height as an input. Note that y⁺ cannot be evaluated on
the mesh alone — it needs a target flow condition, so this step needs a
number from outside the geometry.

**Step 5 — untangling / smoothing, only if step 1 says it is needed.**
Either Allen's aspect-ratio smoothing inside the algebraic stage, or an
elliptic pass. `clean_blocks.untangle` already does barrier-based quality
repair with boundary vertices sliding on the input surface, and that
machinery transfers directly.

The gate at every step is `MESH_QUALITY.md`'s seven metrics against the
current mesh — 0.009 % of cells outside an OpenFOAM default. TFI has to beat
that, not merely produce something.

---

Sources: [Gordon & Hall, Transfinite element methods](https://link.springer.com/article/10.1007/BF01436298) ·
[Transfinite interpolation (overview)](https://en.wikipedia.org/wiki/Transfinite_interpolation) ·
[Perronnet, transfinite interpolation catalogue](https://www.ljll.fr/perronnet/transfini/transfini.html) ·
[Allen, improved transfinite interpolation](https://www.researchgate.net/publication/229806466_Automatic_Structured_Multiblock_Mesh_Generation_Using_Robust_Transfinite_Interpolation) ·
[Zhang et al., structured mesh generation with smoothness controls](https://onlinelibrary.wiley.com/doi/10.1002/fld.1150) ·
[Elliptic grid generation with orthogonality and spacing control](https://www.researchgate.net/publication/23582709_Elliptic_grid_generation_with_orthogonality_and_spacing_control_on_an_arbitrary_number_of_boundaries) ·
[Ali & Tucker, Multiblock structured mesh generation for turbomachinery flows](https://www.researchgate.net/publication/273258150_Multiblock_Structured_Mesh_Generation_for_Turbomachinery_Flows) ·
[Sauer et al., optimization based multi-block-structured grid generation](https://elib.dlr.de/195590/1/Sauer_2023_Numerical_Meth_Engineering-An_optimization_based_multi%E2%80%90block%E2%80%90structured_grid_generation_method.pdf) ·
[GridPro, meshing turbine blades](https://blog.gridpro.com/the-art-and-science-of-meshing-turbine-blades/)
