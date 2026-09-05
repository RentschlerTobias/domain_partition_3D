# Reducing block count by manipulating the frame field (T1_9, v9)

Plan for the next stage. **Not started — awaiting go.**

The working hypothesis behind this stage is: *the remaining blocks that could
be avoided are caused by frame-field singularities, so manipulating the field
to reduce singularities will reduce the block count.*

Section 1 tests that hypothesis against measurements before anything is
planned on top of it, because it is only **partly** true and the part that is
false would send the work in the wrong direction.

---

## 1. Measured starting point

### 1.1 The singular graph is already close to minimal

Measured on v9 after the cavity refill (`clean_blocks.BlockStructure`,
`base_complex.singular_edges`):

```
125 interior singular edges      valence 3: 22   valence 5: 103
130 singular vertices, 10 nodes  all of degree 1 (arc endpoints)
5 singular arcs                  lengths 1.006 … 1.226
0 arcs shorter than 0.10         nothing to shrink
5 of 5 arcs end on the domain boundary
```

| # | edges | valence | length | r (start → end) | z (start → end) |
|---|---|---|---|---|---|
| 1 | 27 | 5 | 1.226 | 1.80 → 0.59 | 1.55 → 1.58 |
| 2 | 27 | 5 | 1.208 | 1.80 → 0.60 | 1.05 → 1.04 |
| 3 | 27 | 5 | 1.202 | 0.60 → 1.80 | 1.11 → 1.12 |
| 4 | 22 | **3** | 1.028 | 1.00 → 1.06 | 1.49 → **2.50** |
| 5 | 22 | 5 | 1.006 | 0.83 → 0.87 | 1.51 → **2.50** |

Three observations that constrain everything below:

**(a) There is no clutter to clean up.** Every arc is long, of uniform
valence, and terminates on the boundary at both ends. There are no zipper
nodes, no complex nodes (all 10 nodes have degree 1), no short arcs. The
generic simplification machinery — shrinking complex singular edges,
collapsing zipper nodes, cancelling stubs — has nothing to act on here. This
is important: it is the opposite of the situation those methods are designed
for, and the reason to expect only modest returns from them.

**(b) Arcs 1–3 are very probably irreducible.** They run radially,
shroud → hub, at mid-passage, and are the 3D expression of the effect that
motivated going 3D in the first place: `hexa_interpolation.py:4-8` records
that hub and shroud carry *topologically different* cross fields (hub four
idx = −1, shroud four idx = +1), a real consequence of blade twist. A field
without those arcs cannot represent that difference. Removing them is not a
simplification, it is a misrepresentation.

**(c) Arcs 4 and 5 are the real candidate.** They are a valence-3 / valence-5
pair — a −1 and a +1 — running roughly parallel from z ≈ 1.5 to the outlet at
z = 2.50, at radii 1.00 and 0.83. That is exactly the configuration the 2D
pipeline already annihilates (`dp3d/clean_separatrix.py:131`, Kowalski 2015)
one dimension down. If this pair can be merged or cancelled, the graph drops
from 5 arcs to 3.

### 1.2 But the block count is not dominated by the singular graph

42 blocks from 5 arcs. The blocks come from the *sheets* the arcs emit and
from how those sheets intersect each other and the boundary patches — not
from the arc count directly. Halving the arc count does not halve the block
count, and the plan must not assume it does. Quantifying that relation is
step 0 below.

### 1.3 A second, independent defect points at the same place

The block-edge wireframe has 63 kinks above 30°. Their distribution
(measured):

| curve type | curves | vertices | median kink | p95 | > 30° |
|---|---|---|---|---|---|
| **shell_hub \| shell_blade** | 14 | 160 | 5.2° | **96.5°** | **45** |
| block \| block | — | 576 | 1.9° | 11.0° | 4 |
| surface \| block | — | 1216 | 1.2° | 8.3° | 14 |
| surface \| surface | — | 442 | 0.5° | 3.9° | 0 |

45 of the 63 sit on the 14 curves that separate `shell_hub` from
`shell_blade`. Those are staircases, and they cannot be smoothed:

```
hex boundary vertices within 0.010 of the true shell_hub|shell_blade
feature curve:  0
staircase distance to that curve: median 0.074, max 0.133  (1–3 cells)
```

**The hex mesh is not aligned to that feature curve at all.** There is no
edge chain to route the label boundary onto and no geometry to smooth toward.
This is a field-alignment failure, and it is the same root cause as the block
count: the field is not doing what the feature constraints asked of it. Any
experiment below should be evaluated on this metric too, not only on block
count — it is the more sensitive of the two.

---

## 2. What the literature offers

Three distinct families. They act at different stages and are not
alternatives to one another.

### 2.1 Correct/constrain the field before parametrization

| work | what it does | relevance here |
|---|---|---|
| Li et al. 2012, *All-hex meshing using singularity-restricted field* | tet split operations that force the field's singularities into hex-meshable configurations | superseded by AlgoHex's own stage |
| Jiang et al. 2014, *Frame field singularity correction* | modifies rotational transitions between tets to repair invalid singularities, then re-smooths | superseded |
| **Liu et al. 2018, *Singularity-constrained octahedral fields*** | enumerates valid local configurations, generalises Hopf–Poincaré to octahedral fields, and generates a field with a **prescribed** singularity graph | **the only route to actually dictating the graph.** Requires solving a large nonlinear mixed-integer algebraic system |
| Palmer, Bommes, Solomon 2020, *Algebraic representations for volumetric frame fields* | odeco frames, SDP relaxation | alternative field representation, no direct singularity control |
| **Liu & Bommes 2023, *Locally Meshable Frame Fields*** | turns a non-meshable field into a locally meshable one; lifts IGM success 2 % → 58 % on HexMe | **this is already running in our pipeline** — it is AlgoHex's `--generate-locally-meshable-field` stage |

Key point: AlgoHex already implements the state of the art for *meshability*.
It does not optimise for *fewness*. That gap is the opening.

### 2.2 Simplify the structure after meshing

| work | method | code |
|---|---|---|
| Gao et al. 2017, *Robust structure simplification for hex re-meshing* | sheet/chord collapse with quality repair | — (already reimplemented here, `collapse_mesh_sheets`) |
| Xu et al. 2021, *Singularity structure simplification via weighted ranking* | ranks collapsible sheets/chords by a weighted function (valence prediction + element quality + width) instead of thickness alone; fixes the early-termination and closed-loop failures of thickness ranking | [github.com/ohehe/HexMeshSimplification](https://github.com/ohehe/HexMeshSimplification) |
| Duan et al. 2023, *Singularity structure simplification via integer linear program* | plans a **global** collapse strategy with an ILP instead of greedy ranking; relaxes singularity constraints for further simplification; handles bad singularities by **sheet inflation**; preserves sharp features | — |
| Duan et al. 2024, *Feature-aware singularity structure optimization* | adds **alignment of singularities to feature lines**, plus inflation to fix high-valence singularities | — |

This family is directly relevant because our greedy collapse **stopped after
one sheet**: round 2 found no sheet that improves. A global ILP plan can
accept a locally-neutral collapse that unlocks a later one, which greedy
cannot. This is the single most likely explanation for our early stop.

### 2.3 Better block extraction from the same mesh

Brückler, Gupta, Mandad, Campen 2022, *The 3D Motorcycle Complex for
Structured Volume Decomposition* — a principled coarse block decomposition,
provably finer-grained than the base complex but with better-behaved
structure. Our `base_complex.py` is the naive version. Worth knowing about,
but it produces *more* blocks, not fewer; relevant only if block *quality*
becomes the binding constraint rather than block count.

---

## 3. What AlgoHex already exposes (checked in the source)

`external/algohex-src/demo/HexMeshing/main.cc` and `src/AlgoHex/Args.hh`.
This is the cheapest lever by a wide margin — no new code, only reruns.

**Singular graph optimisation** (defaults from `Args.hh`):

| flag | default | meaning |
|---|---|---|
| `--field-alignment-weight` | 10.0 | weight of aligning singular edges to the field |
| `--ces-weight` | 1.0 | weight of **shrinking complex singular edges** |
| `--tps-weight` | 1.0 | weight of **shrinking edges at zipper nodes** |
| `--rgl-weight` | 1.0 | regularisation |
| `--rp-weight` | 0.1 | weight of separating twisted singular arc pairs |
| `--dpath-weight` | 100.0 | aligning the dual path to the field when fixing zipper nodes |
| `--merge-turning-point` | **true** | allow merging zipper nodes into valence +1 arcs |
| `--first-iters` / `--second-iters` / `--second-inner-iters` | 13 / 15 / 6 | iteration budget of the two optimisation phases |
| `--early-stop-iter` | — | stop once all singular/feature vertices are locally meshable |
| `--push-boundary-singular-circle` | — | push boundary singular circles into the interior |
| `--fix-turning-point-from-shortest-arc` | — | order in which zipper nodes are fixed |

**Feature handling:**

| flag | default | meaning |
|---|---|---|
| `-d, --dihedral-angle` | **70.0** | dihedral angle defining a feature edge |
| `--full-constraints` | **false** | **fully constrain the frames at features** |
| `--force-feature-threshold` | false | force the feature threshold |
| `-p, --penalty` | — | penalty of the normal alignment |

**`--full-constraints=false` is the prime suspect for §1.3.** Our
`shell_hub | shell_blade` interface is not a geometric feature — there is no
dihedral kink there, it is an artificial cut surface — so the field has no
geometric reason to align to it and is only softly asked to. Turning on full
constraints is a one-flag experiment.

**Diagnostics we are not yet using:**

| flag | value |
|---|---|
| `--animation-out-path` | saves the singularity graph development across the locally-meshable-field iterations, as `.ovm` |
| `--locally-non-meshable-out-path` | dumps the configurations that failed the meshability test |
| `--with-local-meshability-test` | runs the test explicitly |
| `--sub-hexmesh-out-path` | saves the boundary layer of the output hex mesh |

**Checkpoints already on disk** make reruns much cheaper than a cold run:
`T1_9_seamless_v9.hexex` (seamless map) and `T1_9_final_tet_v9.ovm`
(post-singularity-optimisation tet mesh). Passing `--hexex-in-path` together
with `-i` skips field generation, locally-meshable-field generation and
integrability optimisation — 86 % of the runtime (6704 s of 7777 s on the
earlier full-domain run). Experiments that change only quantization or
extraction cost minutes; experiments that change the field cost ~45 min each.

---

## 4. Proposed experiments, in order

Each step states what it costs, what it would prove, and **what result would
falsify it** — because in this project four plausible hypotheses have already
turned out wrong (README "What was learned" 1, 2, 5, 8).

### Step 0 — Establish the arc → block relation (no reruns, ~1 h)

Before trying to reduce arcs, measure what they actually cost. For each of
the 5 arcs, count the blocks whose existence depends on the sheets emanating
from it: remove that arc from the singular-edge set, recompute the base
complex, count blocks. Purely local to `clean_blocks.py`.

*Falsifier*: if removing arcs 4+5 changes the block count by only a handful,
the whole hypothesis behind this stage is wrong and steps 2–4 should be
dropped in favour of §2.2 (ILP simplification), which acts on the structure
directly. **This is a genuine possibility given §1.2 and should be checked
first.**

### Step 1 — Feature alignment (2 reruns × ~45 min)

Rerun v9 with `--full-constraints` on, and separately with a lower
`--dihedral-angle`. Evaluate on the §1.3 metric (staircase kinks, and hex
vertices within 0.010 of the `shell_hub|shell_blade` curve) as well as block
count.

*Why first*: it is the cheapest field-level change, it targets a defect we
have measured precisely, and this project has already learned twice that
feature constraints dominate outcomes (README 1: fewer, cleaner features
consistently gave *worse* results).

*Falsifier*: if the hex mesh still has zero vertices near the feature curve,
the constraint is not being applied to that curve at all, and the problem is
in how `tet_prep_v5` emits feature edges, not in the field solver.

### Step 2 — Singular graph optimisation weights (3–5 reruns × ~45 min)

Sweep `--ces-weight`, `--tps-weight`, `--rp-weight` upward and
`--field-alignment-weight` in both directions. Use `--animation-out-path` on
one run to see whether the optimiser ever considers merging arcs 4 and 5.

*Expectation, stated honestly*: **low.** These weights act on complex edges
and zipper nodes, and §1.1(a) shows we have none. The run is cheap and the
animation output is diagnostic, but a null result is the likely outcome and
should be accepted quickly rather than tuned around.

### Step 3 — Target the ±1 arc pair directly (1–2 weeks)

The real candidate from §1.1(c). Two routes:

**3a, post-hoc (preferred).** Cancel the pair at the *mesh* level with the
sheet-inflation/collapse machinery of Duan et al. 2023/2024 rather than at
the field level. We already have `collapse_mesh_sheets`; inflation is the
missing inverse operation, and it is what the ILP paper uses to remove
"unreasonable" singularities. Extends existing code, no C++ build, and it is
testable against the existing 42-block result.

**3b, field-level.** Prescribe the singularity graph with Liu et al. 2018 and
generate a field with arcs 4 and 5 removed. Correct in principle, but it
means a second C++ build chain (the AlgoHex experience argues strongly
against that, POSTPROCESSING_PLAN.md already records the decision) and the
prescribed graph must be hex-meshable, which for this geometry we cannot yet
verify independently.

*Falsifier for 3a*: if the two arcs are not separated by a single sheet, they
cannot be cancelled by one inflation/collapse pair and the operation count
grows quickly.

### Step 4 — Global ILP instead of greedy collapse (1–2 weeks)

Replace the greedy `collapse_mesh_sheets` acceptance rule with the ILP
formulation of Duan et al. 2023. Directly addresses the measured early stop
(round 2: "no sheet improves"), which greedy cannot escape by construction.
Independent of steps 1–3 and can be done in parallel.

Reference implementation for the weighted-ranking predecessor is available
([HexMeshSimplification](https://github.com/ohehe/HexMeshSimplification));
the ILP paper has no public code, so this is a reimplementation against the
paper, with the existing sheet machinery as the substrate.

---

## 5. Success criteria

| metric | now | target |
|---|---|---|
| blocks | 42 | < 30 |
| non-cuboid blocks | 6 | 0 |
| blocks < 10 cells | 2 | 0 |
| kinks > 30° on block edges | 63 | < 10 |
| hex vertices on the `shell_hub\|shell_blade` curve | 0 | > 0 |
| inverted cells | 0 | **0 — hard, unchanged** |
| boundary Hausdorff to input surface | 0.0999 | no worse |

The inverted-cell count and the Hausdorff distance stay hard constraints, as
in the previous stage.

---

## 6. Risks

1. **The hypothesis may be wrong** (§1.2). Step 0 exists to find that out for
   a day's work rather than a month's. Take its answer seriously.
2. **Arcs 1–3 are load-bearing** (§1.1b). Any method that removes them
   produces a field that cannot represent the hub/shroud topology difference,
   and the mesh will be worse, not simpler. Guard every experiment with the
   inverted-cell and Hausdorff checks.
3. **A second C++ build chain.** Explicitly decided against once already.
   Step 3b would reopen that; step 3a is the route that does not.
4. **Metric traps.** Four times in the previous stage a plausible scalar
   metric improved while the structure degraded (cuboid *share* rising as the
   denominator shrank; soft-min Jacobian raising the worst cell while
   inverting three more; kink *median* improving while the tails worsened).
   Every criterion in §5 is a count or a tail quantile for that reason. Do
   not add a mean.

---

## 7. Out of scope here, but coupled

Re-attaching the removed boundary layers (blade O-grid, hub/shroud prism
layer) is tracked separately. It is coupled to §1.3: once the blade O-grid is
re-attached, `shell_blade` stops being a domain boundary, the
`shell_hub | shell_blade` label boundary disappears, and with it the 14
staircase curves carrying 45 of the 63 bad kinks. If that work happens first,
step 1 should be re-measured afterwards — it may become unnecessary.
