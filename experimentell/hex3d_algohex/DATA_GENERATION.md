# Dataset generation: does the pipeline generalise?

Branch `data_generation`. The block structures are meant to become training
data for a transformer, so the pipeline has to run on geometries it was never
tuned for. Everything up to this point was measured on T1_9 alone.

`data/random_tistos/` holds 21 parametrically generated runners (`investigated/
cand_NNN/mesh.msh`, ~22 MB each, with the design parameters in `params.json`)
plus 8 loose meshes. The meshes are not versioned — 588 MB — but the
**design parameters are**, 27 KB for all 21, because they are the actual input
and every exported sample carries them in its `params` field.

## Where the 21 candidates come from

All five links exist, so a third party can regenerate the geometry without
copying 429 MB from anyone's disk:

| | |
|---|---|
| generator | **dtOO**, `github.com/ihs-ustutt/dtOO`, case `demo/tistos` |
| the parameters | `demo/tistos/build.py` and `machineSave.xml` — **all 30** keys of a `params.json` appear in both, verified 2026-09-11 |
| the values | `data/random_tistos/investigated/cand_*/params.json`, in this repo |
| runtime | `atismer/dtoo-opensuse:stable`, public on Docker Hub |
| mesh-only wrapper | `scripts/dtoo_mesh_only.py` — runs the one `boundedVolume` of type `map3dTo3dGmsh` that produced `T1_9_ru_gridGmsh.msh`, skipping geometry/decompose/solve |

The mesher is the same one that made T1_9, which is why `tet_prep_v5.classify`
reads the same GEOMETRIC entity ids across both — hub 1, shroud 2, inlet 3,
outlet 4, periodic 5 and 6.

What is NOT in any repo here: the sweep driver, i.e. the loop that varied those
30 parameters, wrote `cand_NNN/params.json` and called dtOO 21 times. That is
a small script over a documented mechanism, not a missing capability — but if
the 21 have to be reproduced exactly, the `params.json` values are the record
to feed back in, one case at a time.

## What had to change: one thing

`tet_prep_v5.py` had the source path as a module constant. It now takes
`--msh` and `--out-dir`. Nothing else.

The candidates declare **two physical groups fewer** than T1_9 — no
`aS_ru_inlet_full_0`, no `aS_ru_outlet_full_0` — which looks like a blocker
and is not. `classify()` reads the GEOMETRIC entity ids, and those are
identical across both: hub 1, shroud 2, inlet 3, outlet 4, periodic 5 and 6,
O-grid volume regions 2-6 with 44 800 hexes in both. Same mesher, same block
layout, different blade parameters.

## First candidate, end to end

| stage | cand_000 | T1_9 (v11) |
|---|---|---|
| **tet_prep** | | |
| surfaces classified | all 7 | all 7 |
| rings | 112 / 112 edges, median kink 83.6° / 91.7° | 112 / 112, 81° / 98° |
| feature edges | 530, valence {2: 518, 3: 8} | 534, {2: 522, 3: 8} |
| tets | 43 012 | 42 218 |
| **AlgoHex** (`-n 60000`) | 10.1 min, IGM valid | 13.1 min, IGM valid |
| cells | 61 393 | 61 546 |
| inverted (raw) | **4** | 21 |
| raw blocks | **31** (30 cuboid, 0 tiny) | 117 (114 cuboid, 16 tiny) |
| **clean_blocks** | 28.6 min | ~45 min |
| blocks | 21, 20 cuboid (95 %) | 16, 14 cuboid (88 %) |
| blocks < 10 cells | 0 | 0 |
| inverted | 0 | 0 |
| min / mean scaled Jacobian | **0.2093** / 0.9836 | 0.1524 / 0.9767 |
| pinch | **none** | 2 |
| validator | VALID | VALID |
| **TFI** at h = 0.05 | 12 classes over 60 axes | 11 over 48 |
| refill | 56 869 cells, watertight | 49 550, watertight |
| inverted after refill | **3** (min sJ −0.0575) | 0 |
| cell size p5 / mean / p95 | 0.036 / 0.045 / 0.064 | 0.031 / 0.048 / 0.068 |

**The new geometry is easier than T1_9**, not harder: 4 raw inverted cells
against 21, 31 raw blocks against 117, no pinch, and a better worst cell after
cleanup. T1_9 appears to be an awkward member of this family rather than a
representative one.

**Two things to fix before this scales.**

The refill leaves **3 inverted cells** where T1_9 leaves none. One block of the
21 is a 7-face non-cuboid with no recoverable lattice, so the passthrough
guard fired for the first time — on T1_9 every block was a lattice and that
path had never executed. The inverted cells are the obvious suspects around
the frozen counts at that block's faces, and that is where to look first.

The `T1_9_` prefix is baked into `run_algohex.py`'s output names, so a run on
another geometry still writes `T1_9_hex_gen000.ovm`. Harmless for one sample,
wrong for a dataset.

## Cost per sample

Measured on an unthrottled box, sequential:

| stage | minutes |
|---|---|
| tet_prep | ~3 |
| AlgoHex | 10 |
| clean_blocks (5 collapse rounds, untangle) | 29 |
| TFI refill | ~1 |
| **total** | **~43** |

So the 21 candidates are roughly **16 hours** end to end, one at a time. The
sheet collapse is two thirds of it and is the only stage worth optimising —
and note that the merge ILP, which runs in seconds, does NOT substitute for it
(it produces structures that cannot be refilled; see the 2026-09-08 decision
log).

Two AlgoHex runs must never run concurrently (~6 GB peak each), and this box
throttles hard under hours of sustained full load — see the CPU note in
`HANDOFF.md`.

## The direction classes are physical, and they recur

A direction class is the set of block axes that must carry the same cell
count, because two blocks sharing a face have to agree on the two directions
spanning it and that constraint propagates. They are not an artefact of the
solver: decomposing each axis into the cylindrical basis (e_r, e_theta, e_z)
shows the classes ARE the flow directions, which is what one would hope for a
mesh whose frame field was aligned to the passage.

On cand_002, the largest class:

    22 axes, one from every block, 29 cells
    radial component 1.00, circumferential 0.04, axial 0.05
    axis lengths 1.202 .. 1.224  ->  spread 1.02x

That is the hub-to-shroud direction, threading the entire mesh as a single
degree of freedom. And the pattern repeats across every sample measured:

| sample | classes | radial class | circumferential | meridional | mixed |
|---|---|---|---|---|---|
| cand_001 | 11 | 21 axes, 28 cells, 1.02x | 3 | 4 | 3 |
| cand_002 | 12 | 22 axes, 29 cells, 1.02x | 4 | 5 | 2 |
| cand_003 | 12 | 22 axes, 28 cells, 1.01x | 4 | 5 | 2 |
| cand_004 | 12 | 22 axes, 29 cells, 1.02x | 4 | 5 | 2 |
| cand_005 | 12 | 22 axes, 28 cells, 1.01x | 4 | 5 | 2 |
| cand_006 | 12 | 22 axes, 28 cells, 1.02x | 4 | 4 | 3 |

**Exactly one radial class in every sample**, always spanning every block,
always 28-29 cells, always within 1.02x in axis length. The circumferential
classes separate as cleanly (component 0.96-0.98). The meridional direction
splits over several classes, which is expected: the passage curves, so the
meridian rotates from axial to radial, and the 2-3 "mixed" classes are that
turn.

Two consequences.

**`solve_block_divisions(h_map=...)` finally has a principled use.** The
argument has existed unused; prescribing h per PHYSICAL direction -- finer
radially for the wall layers, coarser circumferentially -- is now meaningful,
and it works precisely because the axis lengths inside these classes are
uniform to 1-2 %. Note the contrast with merging, where a class spread 16.7x
and no single h could serve it.

**For a learned model the class structure is a stable label, not an index.**
"Radial class, 22 axes, 29 cells" means the same thing in every sample of the
family, where block and class NUMBERING does not. Measured on one geometry
family and six samples; whether it survives a different runner type is
untested.

## TODO: a label resolver instead of hard-coded ids

`tet_prep_v5` hard-codes `GEOM_INLET = 3`, `GEOM_HUB = 1`, `OGRID_GEOM =
{2..6}` and so on. Those ids serve exactly two purposes, and it is worth being
precise about them because everything else follows:

1. **What to cut out.** `OGRID_GEOM` selects the volume regions of the blade
   O-grid, which `reduced_boundary` drops and `reattach` puts back verbatim.
2. **What each boundary surface IS.** `classify` maps the tagged 2D elements
   to seven labels, and `feature_graph` then makes a feature edge wherever two
   triangles carry DIFFERENT labels -- no angle criterion at all. The label
   boundaries *are* the feature curves, which is why merging three labels in
   v14 destroyed 224 constraints and the field with them.

A numeric-range convention ("0-100 is cut out, 100+ is geometry") was
considered and rejected. It fails silently when a case exceeds a range, it
carries no meaning in the file itself, and it cannot express which periodic
surface pairs with which.

**Proposed instead: prefixed physical names**, which dtOO already writes and
which are readable in ParaView:

    keep:core          volume, to be meshed by AlgoHex
    skip:ogrid_blade   volume, already block-structured, re-attached verbatim
    wall:hub           real geometry, gets the boundary layer
    flow:inlet         planar flow boundary
    periodic:a         paired by the suffix

The prefix is what the pipeline reads; the rest is free text. A surface with
no known prefix is an ERROR rather than a silent misclassification.

Resolution order: prefixed names, then a per-case JSON mapping, then the
current T1_9 constants as a fallback, so nothing that runs today breaks. The
per-case file is needed regardless for meshes that carry few names --
canadaLight has 12 physical names for 391 2D groups -- and it is more honest
than a range rule that pretends to know.

Note that the naming schemes already differ between machines: tistos writes
`aS_ru_hub_0`, `aS_ru_shroud_0`, `aS_ru_inlet_full_0`; canadaLight writes
`DT_HUB`, `IN_INLET`, `OUT_OUTLET`, `GVRU_WALL`. Both are semantic, neither is
the other.

## The target cell count `-n` barely costs AlgoHex, and dominates ours

The pipeline extracts ~60 000 hexes, collapses them to ~16-22 blocks, and then
TFI throws the interior away and refills at a prescribed h. The fine
resolution is work for the bin. Measured on T1_9, same input, three values of
`-n`:

| structure | cells | blocks | cuboids | tiny | inverted | min sJ | Hausdorff | TFI classes | pinch |
|---|---|---|---|---|---|---|---|---|---|
| n=2000 | 1872 | **16** | **14** | 0 | 0 | 0.1499 | 0.0448 | **11** | 2 |
| n=8000 | 7028 | 19 | 18 | 0 | 0 | 0.1503 | 0.0312 | 9 | **0** |
| n=60000 (v11) | 54460 | **16** | **14** | 0 | 0 | 0.1524 | 0.0341 | **11** | 2 |

n=2000 reproduces v11's block structure exactly from 30x fewer cells.

**Where the time goes**, per AlgoHex's own stage timings:

| stage | n=2000 | n=8000 | n=60000 |
|---|---|---|---|
| frame field | 0.4 | 0.4 | 0.3 min |
| singularities | 1.9 | 1.7 | 1.6 |
| integrability | 4.7 | 3.7 | 3.6 |
| parametrization | 5.1 | 4.9 | 4.8 |
| extraction (HexEx) | 0.0 | 0.0 | 0.1 |
| **AlgoHex total** | 12.2 | 10.7 | 10.3 min |
| **our clean_blocks** | **~1** | ~3 | **~45 min** |

`-n` does not move AlgoHex: the parametrization computes a map, quantization
scales it to the cell count afterwards, and extraction costs seconds. The
saving is entirely on OUR side of the seam, because the sheet collapse
iterates over cells. Per sample: 43 min -> about 17.

**Two caveats.** The boundary error is worse at n=2000 (0.0448 against 0.0341)
because 1872 cells resolve the surface coarsely -- irrelevant for the block
structure, which TFI re-samples from the input geometry, but not for a mesh
used directly. And the invariance is not strict: n=8000 gives a DIFFERENT
structure, 19 blocks with 9 classes and no pinch. There are several valid
collapse endpoints and the resolution decides which one is reached. So `-n` is
a parameter with an effect on the target data, not a free efficiency knob;
sweep a few values over a few geometries before committing to one.

## Checkpoint reuse does not work

`run_algohex.py` and `FRAMEFIELD_PLAN.md` both state that passing
`--hexex-in-path` with `-i` skips field generation and integrability, "86 % of
the runtime", making `-n` sweeps cheap. It was never tested and it fails
immediately: AlgoHex writes the intermediate tet mesh as BINARY OVM (`OVMB`
magic) via `--final-tetmesh-out-path`, and its own `-i` reader rejects that
file with "The specified file might not be in OpenVolumeMesh format! No vertex
section defined!". Both the original VTK and the saved OVM fail the same way.
The `-n` runs above are therefore cold runs.
