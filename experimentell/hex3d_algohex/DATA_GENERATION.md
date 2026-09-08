# Dataset generation: does the pipeline generalise?

Branch `data_generation`. The block structures are meant to become training
data for a transformer, so the pipeline has to run on geometries it was never
tuned for. Everything up to this point was measured on T1_9 alone.

`data/random_tistos/` holds 21 parametrically generated runners (`investigated/
cand_NNN/mesh.msh`, ~22 MB each, with the design parameters in `params.json`)
plus 8 loose meshes. Not versioned — 588 MB of regenerable input, see
`.gitignore`.

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
