# hex3d block structure: what it is optimised for

Grilling session, 2026-09-08. Decisions are appended as they are made.

## Design tree

```
Goal: what is the 3D block structure optimised for, and which basis wins?
│
├── x A  Purpose / objective function            DECIDED
├── ? B  Boundary fidelity: where the error is and what fixes it
├── ? C  Basis: v11m vs v16m
├── ? D  Merge objective: include a TFI degrees-of-freedom term?
└── ? E  Remaining compute: v18 post-processing, plateau moves in the collapse
```

## A. Purpose (decided)

**Decision.** The block structure is training data for a transformer. The
learned structure is then filled by TFI in post-processing, and the result is
run as a CFD mesh. Thesis figures are a by-product, not a driver.

**Options considered.** (1) CFD production mesh first, (2) thesis
presentation, (3) ML dataset, (4) reproducibility/canonical form.

**Why.** Stated directly by the user. It resolves a conflict the scorecard
could not: a structure that a model has to learn wants to be compact and
canonical, which argues for few, regular blocks — while pure CFD accuracy
would argue for whichever basis has the lowest boundary error regardless of
how many blocks it takes.

**Consequences.** Block count and topological regularity are first-class
objectives, not tie-breakers. `basis_report.rank` and the cost vector in
`merge_ilp.solve` both encode a ranking that predates this decision and will
need to be revisited. Cheap to change (a cost vector and a sort key); what is
expensive is re-running post-processing for a different basis (1-2 h at full
speed, 4-9 h while the VPS is throttled).

## B. Boundary fidelity — measurement before deciding

The premise under test was "larger blocks mean a larger Hausdorff error,
because a block is a linear approximation". Measured on the v11 core:

| mesh | cells | blocks | Hausdorff (centroids) | dense 5x5 max | p99 | mean |
|---|---|---|---|---|---|---|
| v11 original | 54460 | 16 | 0.03405 | 0.03543 | 0.01485 | 0.00125 |
| v11m (merged only) | 54460 | 12 | 0.03405 | — | — | — |
| v11 refilled at h=0.05 | 49550 | 16 | 0.03197 | 0.03478 | 0.01494 | 0.00132 |

Block count does not move the error at all: merging changes no geometry, and
v11 and v11m agree to five decimals. Coarsening the mesh by 10 % does not move
it either — max and p99 are unchanged, the mean rises by 5 %.

Splitting the error by where it is sampled, and by surface:

| sample | n | max | p99 | mean |
|---|---|---|---|---|
| boundary vertices | 10442 | 0.03493 | 0.01448 | 0.00116 |
| face centres | 10442 | 0.03405 | 0.01541 | 0.00135 |

| surface | n faces | max | p99 |
|---|---|---|---|
| inlet | 1064 | 0.00000 | 0.00000 |
| outlet | 980 | 0.00000 | 0.00000 |
| periodic_A | 1484 | 0.00170 | 0.00123 |
| periodic_B | 1596 | 0.00139 | 0.00107 |
| bl_interface_hub | 1943 | 0.00840 | 0.00598 |
| bl_interface_shroud | 1945 | 0.00326 | 0.00280 |
| **ogrid_interface** | 1430 | **0.03405** | **0.02310** |

Two conclusions. The error is **not chordal** — vertices are as far off the
surface as face interiors (0.00116 against 0.00135 mean), so higher-order or
spline block edges would address a term that is roughly 0.0002 in size. And it
is **entirely one surface**: all 200 worst faces lie on `ogrid_interface`,
every other surface is at or below 0.0084, and the inlet and outlet are exact.

`ogrid_interface` is not real geometry. It is the cut face towards the blade
O-grid that was removed from the AlgoHex domain and is re-attached later by
`reattach.py`. Its measured mismatch against the O-grid block it must meet is
median 0.0218, p95 0.0500, max 0.0584 (`reattach.py` interface report).

## A (amendment). Block count is the primary criterion

Stated by the user after the measurements above: block count is decisively the
biggest criterion. This settles C in advance — v11m has 12 blocks against
v16m's 75 — and it changes what the merge ILP should optimise: its cost vector
currently minimises block count only among groups seeded by *defective*
blocks, which is a repair objective, not a minimisation one.

## B. Align the O-grid cut face to the O-grid block (decided: F2)

**Decision.** The core's `ogrid_interface` is aligned to the re-attached
O-grid block's boundary, not to the input triangulation.

**Options considered.** F1 project onto the input triangulation, F2 project
onto the O-grid block, F3 rely on input re-meshing (v16-style), F4 leave it
non-conforming.

**Why.** The surface carries 100 % of the boundary error (max 0.03405 against
0.0084 for the next-worst surface, inlet and outlet exactly zero), and it is
not real geometry: it is the cut face towards the removed blade O-grid. There
is no "true" geometry to be faithful to there — what is physically meaningful
is that the core and the re-attached O-grid share the same surface.

**Correction found while measuring.** `reattach.interface_gap` compares
CENTROIDS of differently sized quads and reports median 0.0218 / p95 0.0500 /
max 0.0584. The true point-to-surface distance from core vertices to the
O-grid surface is median **0.00214** / p95 0.01566 / max 0.02734, against a
local cell edge of 0.042. The gap is 5 % of a cell at the median, not 50 %.
The metric overstates it by an order of magnitude and should be replaced by a
point-to-triangle measure.

**Consequences.** The projection is a small motion, so it does not threaten
cell validity at the median; the worst vertices move 0.65 of a cell. Cheap to
reverse (delete the projection step). Open: at which stage it is applied.

## B2. The projection happens inside the refill (decided: G2)

**Decision.** The `ogrid_interface` face grids are pulled onto the O-grid
block's surface *while the TFI refill resamples them*, not as a separate
snapping pass afterwards.

**Options considered.** G1 snap after the refill and let `untangle` absorb it,
G2 fold it into the refill, G3 move the O-grid block instead, G4 cut the input
on the O-grid boundary in `tet_prep` from the start.

**Why.** `refill_block` already rebuilds each block's six faces from sample
points and fills the interior by Gordon-Hall, so a boundary that has moved is
absorbed by construction. The motion is small — median 0.00214 against a local
cell edge of 0.042 — but G1 would have to recover it afterwards with the same
local smoothing whose quality guards misfired twice in this session. G3 was
rejected because the O-grid block is reused verbatim from the source MSH and
is therefore the more trustworthy reference of the two. G4 costs a new AlgoHex
run plus post-processing, which is 5-10 hours while the VPS is throttled.

**Consequence to watch.** The projection must be applied per VERTEX at complex
level, not per block face: a vertex on `ogrid_interface` is also a corner of
adjacent blocks' faces, and projecting it in one block but not the other
breaks the coordinate weld in `refill_complex`. `check_watertight` catches
exactly that and the refill refuses to write, so the failure cannot pass
silently.

**Status.** Decided, not yet implemented — it needs a refill run.

## C + D. Reference structure: v11 merged to 6 blocks (`v11m6`)

**Decision.** `T1_9_blocks_v11m6` — the v11 basis, merged by the ILP with
`--seed all` — is the new reference structure. Block count is the primary
criterion (decision A), and the seed set, not the group size, is what
determines it:

| basis | blocks | smallest | non-cuboid | pinch | inverted | min sJ | Hausdorff | validator | TFI classes |
|---|---|---|---|---|---|---|---|---|---|
| v11 | 16 | 224 | 2 | 2 | 0 | 0.1524 | 0.03405 | VALID | 11 |
| v11m | 12 | 476 | 2 | 0 | 0 | 0.1524 | 0.03405 | VALID | 8 |
| **v11m6** | **6** | 448 | 4 | 0 | 0 | 0.1524 | 0.03405 | VALID | 4 |
| v16 | 103 | 3 | 2 | 0 | 0 | 0.0159 | 0.01377 | VALID | 15 |
| v16m | 26 | 51 | 23 | 0 | 0 | 0.0159 | 0.01377 | VALID | 3 |

Every structure is VALID, every block is a valid TFI lattice, and the mesh
quality columns are identical within each family because merging touches no
geometry.

**Why v11 and not v16.** With block count primary, 6 beats 26. v16's remaining
advantage is boundary fidelity (0.01377 against 0.03405), but decision B
removes most of that argument: the entire v11 error sits on `ogrid_interface`,
which the projection targets directly.

**What it costs.** The 6-block structure is unbalanced — one block holds
26 600 of 54 460 cells — and the whole mesh's resolution is then controlled by
4 direction classes instead of 11. For a transformer that has to learn the
structure this is the right trade; for hand-grading a boundary layer inside
the core it would not be.

**Correction recorded.** "Direction class stuck at one division" was presented
earlier in this session as a hard cap on resolution. It is not:
`solve_block_divisions` re-solves every count from the prescribed h, so the
current 1 is descriptive, not constraining. It flags a thin block (aspect
ratio), nothing more.

## Artifacts

`T1_9_blocks_{v11m,v11m6,v16m}.{vtk,msh}` plus `_edges`, `_nfaces`,
`_quality` VTKs, written through `clean_blocks.write_blocks` so the merged
structures carry the same artifacts as every other deliverable. The `.msh`
holds the block edges as 1D elements with the curve id as physical tag.

## Still open

* **E** — v18 post-processing (raw run done, structure not measured) and the
  plateau-move experiment in the sheet collapse. Both need a full
  `clean_blocks` run and wait for the VPS throttle to lift.
* Implementation of G2.
* What "learnable" requires beyond block count: canonical block numbering, a
  stable ordering across geometries. Not yet sharp enough to decide.

## D (revision). The 6-block structure does not carry the generator

Measured immediately after the decision above, on the refill at h = 0.05. The
lower bound of a direction class is the MAX over its axes, so a class that
bundles a long and a short axis forces the long axis's cell count onto the
short one. With few classes, that is every class.

| structure | classes | bound rule | cells | inverted | min sJ | edge p5 … p95 |
|---|---|---|---|---|---|---|
| v11m6 | 4 | max | 249100 | 2 | −0.6318 | 0.0099 … 0.0525 |
| v11m6 | 4 | median | 30240 | 4 | −0.7470 | 0.0314 … 0.1639 |
| v11m | 8 | max | 76025 | 3 | −0.1121 | 0.0102 … 0.0651 |
| v11m | 8 | median | 43488 | **0** | **+0.1127** | 0.0248 … 0.1002 |
| v11 | 11 | max | 49550 | **0** | **+0.1357** | 0.0310 … 0.0684 |
| v11 | 11 | median | 43008 | **0** | **+0.1484** | 0.0329 … 0.0724 |

Target was h = 0.05 and roughly 50 000 cells.

**The trade-off is direct and was not visible in any static metric.** Merging
costs nothing in mesh quality — the merged structures are byte-identical to
their basis in inverted count, min sJ and Hausdorff — but it destroys the
ability to REGENERATE the mesh at a prescribed resolution, because it removes
the degrees of freedom that resolution control needs. v11m6 overshoots the
target by 5x under the max rule and produces cells spread over 5x in size
under the median rule; both leave inverted cells.

So "block count is the primary criterion" holds only down to the point where
the structure still has enough direction classes to be filled. Between 16 and
12 blocks the generator still works; at 6 it does not.

**Status.** `v11m6` is NOT a usable reference for the TFI stage. The choice is
between v11 (16 blocks, 11 classes, works under both rules, tightest cell-size
spread) and v11m (12 blocks, 8 classes, works under the median rule, cells up
to 2x the requested h). Pending the user's call.
