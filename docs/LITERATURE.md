# Literature

Every source this repository actually leans on, in one place, grouped by the
role it plays in the pipeline rather than alphabetically — the question being
answered is usually "what do we already rely on for X?", not "what did someone
publish".

**Status column.** The distinction matters more than the citation:

| | meaning |
|---|---|
| **runs** | the method is executing in our pipeline, ours or as a dependency |
| **reimplemented** | we wrote it from the paper; our code, their method |
| **informed** | read, shaped a decision, no code |
| **rejected** | read, tried or evaluated, measured worse or inapplicable — the reason is recorded |
| **shortlist** | identified as the likely next step, not started |

Where a full bibliographic record was never established in this branch, the
entry says so and points at the repo text that cites it. Please do not fill
those in from memory — look them up.

---

## A. The hex-meshing chain we run

| work | what it does | where | status |
|---|---|---|---|
| **AlgoHex** — [github.com/cgg-bern/AlgoHex](https://github.com/cgg-bern/AlgoHex), AGPLv3, pinned at commit `3519289` | one CLI (`HexMeshing`) over all four stages: frame-field init → locally-meshable field → seamless parametrisation + quantisation → extraction | `run_algohex.py`, `external_patches/Dockerfile.portable` | **runs** |
| **Liu & Bommes 2023, _Locally Meshable Frame Fields_** — [algohex.eu](https://www.algohex.eu/publications/locally-meshable-frame-fields/) | repairs a non-meshable octahedral field into a locally meshable one; lifts IGM success on the HexMe benchmark from 2 % to 58 % | AlgoHex's `--generate-locally-meshable-field` stage | **runs** |
| **Lyon et al. 2016, HexEx / libHexEx** — [doi:10.1145/2897824.2925976](https://doi.org/10.1145/2897824.2925976) | robustly extracts a hex mesh from an imperfect integer-grid map | inside AlgoHex | **runs** |
| **QGP3D** — [github.com/HendrikBrueckler/QGP3D](https://github.com/HendrikBrueckler/QGP3D), GPLv3 | quantised global parametrisation; internally builds a 3D Motorcycle Complex, so it already computes a block decomposition during quantisation | inside AlgoHex | **runs** |
| **NeurFrame** — [arXiv:2603.12820](https://arxiv.org/html/2603.12820v1) | the pipeline shape this branch follows: octahedral field on a tet mesh → integer-grid map → robust extraction | `PLAN.md†` names it as the reference architecture | **informed** |

Note for anyone reading AlgoHex's licence chain: AlgoHex is AGPLv3 and
libHexEx/QGP3D are GPLv3. Fine for internal research in this repo; it
constrains redistribution, not use.

## B. Block extraction — turning a fine hex mesh into a coarse complex

| work | what it does | where | status |
|---|---|---|---|
| **Gao et al. 2015/2017** — [doi:10.1145/3130800.3130848](https://dl.acm.org/doi/10.1145/3130800.3130848) | the base complex: seed from every facet incident to a **singular** edge, expand through the **opposite facet across each regular edge** until termination; the hex groups bounded by the resulting sheets are the coarse blocks | `base_complex.py` (`singular_edges`, `sheet_faces`, `blocks_from_cut`) | **reimplemented** |
| **Brückler, Gupta, Mandad, Campen 2022, _The 3D Motorcycle Complex for Structured Volume Decomposition_** | a principled coarse decomposition, provably finer-grained than the base complex but better-behaved | `FRAMEFIELD_PLAN.md†` §2.3 | **informed** — it produces *more* blocks, not fewer, so it is relevant only if block *quality* becomes the binding problem rather than block count |

The expansion rule in Gao's construction is easy to get wrong in a way that
still runs: stepping to the opposite face **within** a hex walks a 1D chain of
faces and never closes into a separating surface, which left 9384 of 9792
hexes in a single component on the cylinder test. The expansion must go across
the **edges** of a face. The working definition of "opposite facet across edge
g" is: rotate twice around g. See the docstring of `base_complex.sheet_faces`.

## C. Structure simplification — fewer blocks from the same mesh

| work | what it does | where | status |
|---|---|---|---|
| **Gao et al. 2017, _Robust structure simplification for hex re-meshing_** | sheet/chord collapse with quality repair afterwards | `clean_blocks.py` (`mesh_sheets`, `collapse_sheet`, `collapse_mesh_sheets`) | **reimplemented** |
| **Xu et al. 2021, _Singularity structure simplification via weighted ranking_** — [github.com/ohehe/HexMeshSimplification](https://github.com/ohehe/HexMeshSimplification) | ranks collapsible sheets/chords by a weighted function (valence prediction + element quality + width) instead of thickness alone; fixes thickness ranking's early-termination and closed-loop failures | `FRAMEFIELD_PLAN.md†` §2.2 | **shortlist** |
| **Duan et al. 2023, _Singularity structure simplification via integer linear program_** | plans a **global** collapse strategy by ILP instead of greedy ranking; relaxes singularity constraints; handles bad singularities by sheet **inflation**; preserves sharp features | `FRAMEFIELD_PLAN.md†` §2.2, §5 | **shortlist — the most likely fix for our measured early stop** |
| **Duan et al. 2024, _Feature-aware singularity structure optimization_** | adds alignment of singularities to feature lines, plus inflation for high-valence singularities | `FRAMEFIELD_PLAN.md†` §2.2 | **shortlist** |

Two things worth carrying over from our own implementation of Gao's collapse,
because both cost this branch real time:

**What the collapse is.** It does not merge neighbouring hexes into larger
ones. `mesh_sheets` groups edges into sheets as equivalence classes of
*parallel* edges (two edges share a sheet if they are parallel inside some hex,
transitively) — the dual view, and a **different object** from the *separating*
sheets of the base complex, which are bounded by singular edges.
`collapse_sheet` then welds the two endpoints of every edge of one sheet: each
hex of that sheet loses distinct vertices, degenerates and drops out, and the
two sides of the deleted layer become direct neighbours. So a whole layer is
deleted, and the point is not cell reduction but removing singularities so
that the base complex computed *afterwards* finds fewer blocks.

**Repair after, never before.** Our greedy version accepts a sheet only if it
adds no inverted cell, does not push the boundary further from the input
surface, and strictly improves (excess faces, then block count). That guard is
a comparison, so whatever it compares against sets the height of the bar — and
it must stay pinned to the raw mesh as AlgoHex produced it. Measured twice:
a pre-collapse untangle lifts v11 from 21 inverted / min sJ −0.4938 to
6 / −0.3116, and the collapse that used to reach 16 blocks then stops at 22; a
repair between rounds lifted the worst cell from −0.19 to +0.15 after round 2,
after which no sheet could clear the new bar at all. Gao's pipeline repairs
after simplifying for exactly this reason: quality damage from a collapse is
repairable, structure thrown away is not.

## D. Frame fields and singularity control

| work | what it does | status |
|---|---|---|
| **Li et al. 2012, _All-hex meshing using singularity-restricted field_** | tet split operations that force field singularities into hex-meshable configurations | **rejected** — superseded by AlgoHex's own stage |
| **Jiang et al. 2014, _Frame field singularity correction_** | modifies rotational transitions between tets to repair invalid singularities, then re-smooths | **rejected** — superseded |
| **Liu et al. 2018, _Singularity-constrained octahedral fields_** | enumerates valid local configurations, generalises Hopf–Poincaré to octahedral fields, generates a field with a **prescribed** singularity graph | **informed** — the only route to actually dictating the graph; requires solving a large nonlinear mixed-integer algebraic system |
| **Palmer, Bommes, Solomon 2020, _Algebraic representations for volumetric frame fields_** | odeco frames, SDP relaxation | **informed** — alternative field representation, no direct singularity control |

The load-bearing conclusion of this section, recorded in `FRAMEFIELD_PLAN.md†`:
AlgoHex already implements the state of the art for *meshability*. It does not
optimise for *fewness*. That gap is the opening, and section C is where it is
addressed — not here.

Related decision: **D, in `docs/decisions/2026-09-11-hex3d-dataset-pipeline.md`**
rejects the frame field as the source of curved block edges. It is never
written out, it is a per-tet direction rather than a curve, and it is the stage
that produced the mesh in the first place. The fine-mesh polyline is used
instead (`block_edges.py`).

## E. Layout embedding — optimising geometry with the topology held fixed

| work | what it does | status |
|---|---|---|
| **Heuschling, Lim, Kobbelt 2026, _Embedding Optimization of Layouts via Distortion Minimization_** — Eurographics 2026 / Computer Graphics Forum 45(2), RWTH Aachen. Code announced at [github.com/7-AlexH/layout-embedding-optimization](https://github.com/7-AlexH/layout-embedding-optimization) | takes a target surface, a layout **connectivity**, and an initial embedding, and optimises the embedding **geometrically** while preserving connectivity strictly. Repositions layout nodes, embeds arcs as piecewise geodesic curves, inserts extra nodes along arcs for flexibility. Multi-resolution: optimise coarse, prolongate to high-resolution | **informed, added 2026-09-11** |

This is a category the branch had no entry for, and no pipeline stage for
either: our block corners and edges sit wherever AlgoHex's integer-grid map put
them, and nothing ever improves that placement. Note what the paper does *not*
do — it never changes the layout's topology, so it is no help for the fewness
question of section C.

It became interesting on 2026-09-11 because both defects measured in T5/T9 are
embedding defects rather than topology defects:

* a **collapsed block edge** in v11, chord 7.6e-5 against a median of 0.7195 —
  a factor of 10⁴, i.e. two layout nodes practically coincident. Node
  repositioning is exactly the operation for that.
* **winding edges** in cand_001, arc/chord up to 4.506, which drive the cubic
  fit residual to 32.7 %. An arc that wraps 4.5x its chord is a poor embedding
  of that arc, and a geodesic re-embedding would shorten it.

Two honest limits. It is a **surface** method (disk-like patches on a 3D
surface), so it transfers to our quad shell but not to interior block edges,
which have no surface to be geodesic on. And the announced code is not out yet.

## F. Transfinite interpolation and structured grid generation

| work | what it does | where | status |
|---|---|---|---|
| **Coons 1967, the Coons patch** | the bilinearly blended surface patch from four boundary curves | `dp3d/tmesh.py::_coons`, `curved_refill.py::_coons` | **reimplemented** |
| **Gordon & Hall 1973, _Transfinite element methods_** — [doi:10.1007/BF01436298](https://link.springer.com/article/10.1007/BF01436298) | transfinite interpolation: the interpolant matches the prescribed data on a whole continuum (the boundary), not at finitely many points — hence "transfinite". Defines the transfinite element as an **invertible** map | `tfi.py::tfi`, the trilinear boolean sum P1 ⊕ P2 ⊕ P3 | **reimplemented** |
| **Perronnet, transfinite interpolation catalogue** — [ljll.fr/perronnet/transfini](https://www.ljll.fr/perronnet/transfini/transfini.html) | explicit forms for the square and the **cube** (Coons, Gordon, Hall) | `TFI_RESEARCH.md†` — the 3D form was taken from here | **informed** |
| **Allen (2008), _Towards automatic structured multiblock mesh generation using improved transfinite interpolation_** — [ResearchGate](https://www.researchgate.net/publication/229806466_Automatic_Structured_Multiblock_Mesh_Generation_Using_Robust_Transfinite_Interpolation) (full record not established in this branch, see `TFI_RESEARCH.md†:114`) | improved/robust TFI for multiblock generation | `TFI_RESEARCH.md†` | **informed** |
| **Vinokur (1983), _On one-dimensional stretching functions for finite-difference calculations_** (record as cited at `TFI_RESEARCH.md†:141`) | the near-wall point distribution itself | `tfi.py::clustered_fractions` | **informed** |
| **Winslow / Thompson–Thames–Mastin**, elliptic grid generation | elliptic smoothing of a structured grid | `TFI_RESEARCH.md†:102` | **informed** |
| **Thomas–Middlecoff** elliptic smoothing | elliptic smoothing with control functions that preserve boundary spacing into the interior, unlike plain Winslow which relaxes towards uniform | `dp3d/tmesh.py:973` | **reimplemented** (2D route) |
| **Elliptic grid generation with orthogonality and spacing control** — [ResearchGate](https://www.researchgate.net/publication/23582709_Elliptic_grid_generation_with_orthogonality_and_spacing_control_on_an_arbitrary_number_of_boundaries) | orthogonality and spacing control on an arbitrary number of boundaries | `TFI_RESEARCH.md†` | **informed** |
| **Zhang et al., structured mesh generation with smoothness controls** — [doi:10.1002/fld.1150](https://onlinelibrary.wiley.com/doi/10.1002/fld.1150) | smoothness controls | `TFI_RESEARCH.md†` | **informed** |
| **Ali & Tucker, _Multiblock structured mesh generation for turbomachinery flows_** — [ResearchGate](https://www.researchgate.net/publication/273258150_Multiblock_Structured_Mesh_Generation_for_Turbomachinery_Flows) | the turbomachinery-specific multiblock practice | `TFI_RESEARCH.md†` | **informed** |
| **Sauer et al. 2023 (DLR), _An optimization-based multi-block-structured grid generation method_** — [elib.dlr.de](https://elib.dlr.de/195590/1/Sauer_2023_Numerical_Meth_Engineering-An_optimization_based_multi%E2%80%90block%E2%80%90structured_grid_generation_method.pdf) | optimisation-based multiblock generation | cited in `dp3d/tmesh.py`, `TFI_RESEARCH.md†:148` | **informed** |
| **GridPro, _The art and science of meshing turbine blades_** — [blog.gridpro.com](https://blog.gridpro.com/the-art-and-science-of-meshing-turbine-blades/) | industrial practice for exactly our geometry class | `TFI_RESEARCH.md†` | **informed** |
| [Transfinite interpolation, overview](https://en.wikipedia.org/wiki/Transfinite_interpolation) | orientation only | | **informed** |

One implementation invariant that comes from our own measurements rather than
from any of these: the resampling in `tfi.refill_block` is parameterised by
normalised **INDEX**, not arc length. That is what makes neighbouring blocks
weld, and it must not be changed. Arc length is used in the *other* direction,
for subdivision (T8), where dividing a shared edge by arc length of the true
curve is what makes two blocks agree on the division points.

## G. The 2D route — surface partition and streamlines

Cited by `dp3d/`, the 2D/surface predecessor of this branch. Kept here because
the 3D work reuses its vocabulary and, in T5, lifted its cubic fit.

| work | what it does | where | status |
|---|---|---|---|
| **Kowalski 2015** (record as cited in `dp3d/`, e.g. `partition_surface.py:1043` "sec 4", `clean_separatrix.py:263` "sec 3") | singularity pair annihilation, the 3/5 invariant, separatrix blending between two singularities (Eq 28), the singularity-efficiency metric | `dp3d/partition_surface.py`, `dp3d/clean_separatrix.py`, `dp3d/field/quad_partition_validator.py` | **reimplemented** |
| **Xiao et al. 2020** (record as cited at `dp3d/field/streamline_merging.py:204` "Algorithm 2", `partition_surface.py:1094` "Case 2, Fig 8b") | stack-based streamline merging | `dp3d/field/streamline_merging.py` | **reimplemented** |
| **The 2D higher-order edge study**, `meshtron/docs/ho_quad_transformer/06_edge_geometry_study.md` | in-house, not published: compares hermite / quadratic / cubic Bezier edge representations over 330 164 edges. `cubic_bezier` wins both axes at once — 4 tokens per edge and 0.74 % median error — and reports **0 inflections** across the whole dataset | the reference numbers T5 is measured against; `block_edges.py` lifted `_fit_cubic_bezier` from `meshtron/prototype_twostage.py` | **reimplemented in 3D** |

A caution recorded in `block_edges.py::curvature_sign_changes`: the 2D study's
inflection test (sign changes of smoothed curvature) does **not** transfer to
3D block edges, and that was established by measurement, not preference. On
the 85 v11 edges that provably do not cross their chord, |k|·chord² reaches
5.995, while the 4 that genuinely S-curve only reach 2.740 — the noise exceeds
the signal, so no threshold separates them. Curvature differentiates twice and
amplifies the fine mesh's own ~1 %-of-chord vertex jitter by ~1/h². The
replacement integrates instead: how often the polyline crosses its own chord.

## H. Mesh quality — the CFD side

`mesh_quality.py` follows OpenFOAM's `checkMesh` conventions, and its limits
are OpenFOAM's own defaults rather than invented thresholds. Details and the
measured numbers are in `MESH_QUALITY.md`.

* [OpenFOAM `meshQualityDict`](https://github.com/OpenFOAM/OpenFOAM-7/blob/master/etc/caseDicts/mesh/generation/meshQualityDict) — the defaults themselves — **runs**
* [OpenFOAM mesh quality guide](https://www.openfoam.com/documentation/guides/latest/doc/guide-meshing-snappyhexmesh-meshquality.html) — **informed**
* [SimScale, Mesh Quality](https://www.simscale.com/docs/simulation-setup/meshing/mesh-quality/) — **informed**
* [Resolved Analytics, mesh quality and CFD solution accuracy](https://www.resolvedanalytics.com/cfd-in-practice/what-is-effect-of-cfd-mesh-on-solution-accuracy) — **informed**
* [wolfdynamics, meshing preliminaries and quality assessment](https://www.wolfdynamics.com/wiki/meshing_preliminaries_and_quality_assessment.pdf) — **informed**

## J. In-house tooling this pipeline depends on

Not literature, but the same question applies — what do we already rely on?

| piece | where | what it gives us |
|---|---|---|
| **dtOO** | `github.com/ihs-ustutt/dtOO`, case `demo/tistos` | the parametric geometry generator behind every input mesh. All 30 `params.json` keys are its const-values |
| **`DtooAdapter`** | `eigenfrequencies/src/eigenfrequencies/adapters/dtoo/adapter.py` | `export_mesh({label: value}) -> .msh` and `design_bounds() -> {label: (min,max)}`. Imports `dtOOPythonSWIG` directly, so it runs INSIDE the dtOO container |
| **tistos machine config** | `eigenfrequencies/adapters/machines/tistos.yaml` | the 30 design parameters **with bounds** — checked key for key against our `params.json` |
| **dtOO on the cluster** | `eigenfrequencies/cluster/enroot_dtoo_import.md`, `submit_dtoo_enroot_smoke.sh` | the enroot route, plus the two env-sourcing traps |
| **cluster sizing** | `eigenfrequencies/docs/cluster-resource-sizing.md` | measured partitions; read our `docs/cluster-enroot-findings.md` §7 before trusting its recommendation |

## I. Benchmarks

* **HexMe** — the hex-meshing benchmark the AlgoHex papers report against
  (2 % → 58 % IGM success for locally meshable fields; feature-surface
  preservation only 34.6 % on average across methods). Used here as the
  reference for what "hard" means and as the format model for
  `data/T1_9/T1_9_tet.vtk`. See `PLAN.md†:27`, `:130`, `:151`. — **informed**

---

## Where the citations live in the repo

This file is the index; the argument is usually in the document that cites the
work, and those are worth reading before acting on an entry here.

| document | what it argues |
|---|---|
| `experimentell/hex3d_algohex/PLAN.md†` | why AlgoHex and not a hand-rolled field/IGM/extraction chain |
| `experimentell/hex3d_algohex/FRAMEFIELD_PLAN.md†` | the frame-field and simplification survey, §2.1-2.3, and why the opening is *fewness* |
| `experimentell/hex3d_algohex/TFI_RESEARCH.md†` | transfinite interpolation, from Coons and Gordon–Hall to the industrial practice |
| `experimentell/hex3d_algohex/MESH_QUALITY.md` | the OpenFOAM metrics and our measured values |
| `docs/decisions/2026-09-11-hex3d-dataset-pipeline.md` | the ten dataset-pipeline decisions, with the alternatives that were rejected |
| `docs/decisions/2026-09-08-hex3d-block-structure-objective.md` | what the block structure is for, and the Coons/Gordon-Hall boundary-error numbers |

† Planning document removed in the 2026-10-03 docs cleanup; read it with
`git show a583e59:experimentell/hex3d_algohex/<file>`.
