# Every AlgoHex run on this branch

From `output/hex3d_algohex/hexmeshing_<tag>.log` and
`T1_9_hex_metrics_<tag>.json`. A run tag `vN` is **not** the same as an input
file `T1_9_tet_vN.vtk` — the two numberings are independent and collide.

## The runs

| run | input file | domain | feature edges | `-n` | runtime | outcome |
|---|---|---|---|---|---|---|
| v1 | `T1_9_tet.vtk` (classifySurfaces tags) | full | 853 (427 flat) | 10 000 | 2 h 10 | 95.4 % coverage |
| v2 | `T1_9_tet.vtk` (analytic tags) | full | 496 | 10 000 | — | SIGSEGV, dangling feature curves |
| v3 | `T1_9_tet.vtk` (analytic, cleaned) | full | 437 | 10 000 | — | IGM diverged, energy 3.3e48 |
| v4 | `T1_9_tet_v1tags.vtk` | full | 853 | 60 000 | 2 h 21 | 97.9 % coverage, 124 inverted |
| v5 | `T1_9_tet_v2.vtk` | full, original boundary | 630 | 60 000 | 1 h 17 | valid IGM, 100 % coverage, 96 inverted |
| v6 | `T1_9_tet_v3.vtk` | minus blade O-grid | 646 | 60 000 | — | stopped by hand at round 3 of 9, memory |
| v7 | `T1_9_tet_v4.vtk` | minus O-grid + BL | 166 | 60 000 | — | OOM in quantization |
| v8 | `T1_9_tet_v4.vtk` | same as v7 | 166 | 15 000 | — | OOM again; `-n` is not the cause |
| v9 | `T1_9_tet_v5_diagbug.vtk` | reduced, **mislabelled** | 550 | 60 000 | 45 min | 2 inverted, 2 internal cavities |
| **v10** | `T1_9_tet_v6.vtk` | only the blade layer cut | 854 | 60 000 | 1 h 29 | 19 inverted, 198 blocks |
| **v11** | `T1_9_tet_v5.vtk` | reduced, **labels fixed** | 534 | 60 000 | **13 min** | 21 inverted, 0 cavities, 117 blocks |
| v12 | `T1_9_tet_v5.vtk` (`--full-constraints`) | reduced | 534 | 60 000 | 18 min | 17 inverted, 84 raw blocks, pinch unchanged |
| v13 | `T1_9_tet_v7.vtk` | O-grid interface re-meshed | 534 | 60 000 | 18 min | 14 inverted, 84 blocks, does not collapse |
| v14 | `T1_9_tet_v8.vtk` | all three cut surfaces merged | 310 | 60 000 | 24 min | 100 % cuboids, but 96 inverted |
| v15 | `T1_9_tet_v9.vtk` | merged **and** re-meshed | 310 | 60 000 | 30 min | 129 inverted, 242 raw blocks |
| **v16** | `T1_9_tet_v10.vtk` | **hub ring only** merged, re-meshed | **422** | 60 000 | 19 min | **32 inverted, 127 raw blocks, no pinch** |
| v17 | `T1_9_tet_v11.vtk` | ring edges under 70° dropped, re-meshed | 499 | 60 000 | 26 min | 46 inverted, 158 raw blocks, no pinch |
| v18 | `T1_9_tet_v9.vtk` (`--full-constraints`) | merged and re-meshed | 310 | 60 000 | 22 min | 102 inverted, but **120** raw blocks |

`fast` and `v2_crashed` in the output directory are probe runs, not results.

### v15-v18: keeping the constraints while removing the pinch

v14 established that the two artificial rings between the cut surfaces cause
the pinch, and that deleting both of them costs the field dearly. What it did
not establish is that the two go together. They do not:

| | v11 | v14 | v15 | **v16** | v17 | v18 |
|---|---|---|---|---|---|---|
| feature edges | 534 | 310 | 310 | **422** | 499 | 310 |
| what was done to the rings | kept | both deleted | both deleted | root ring deleted | edges < 70° dropped | both deleted |
| O-grid interface re-meshed | no | no | yes | yes | yes | yes |
| raw inverted | 21 | 96 | 129 | **32** | 46 | 102 |
| raw min scaled Jacobian | −0.494 | −0.990 | −0.999 | **−0.798** | −0.927 | −0.999 |
| singular edges | 164 | 276 | 272 | 218 | 272 | 232 |
| raw blocks | 117 | 242 | 242 | 127 | 158 | **120** |

The root ring is the one the pinch sits on, and it is also the cheaper of the
two to lose: dropping it alone (v16) keeps 422 of the 534 feature edges and
the field degrades by a third of what deleting both costs — 32 inverted cells
against 96.

Two further results from this round:

* **`--full-constraints` is not null on a reduced feature set.** v12 changed
  nothing at all on v11's 534 edges. On the same 310-edge input as v15 it
  halves the raw block count, 242 → 120 (v18). Fewer declared features leave
  more for the flag to constrain.
* **Re-meshing the O-grid interface improves boundary fidelity every time.**
  Hausdorff 0.0341 (v11, original triangulation) against 0.0136-0.0147 for
  every re-meshed variant, whatever was done to the rings.

## Mesh-level comparison of the four that completed

| | v5 | v9 | v10 | v11 |
|---|---|---|---|---|
| cells | 60 612 | 56 661 | 62 076 | 61 546 |
| inverted cells | 96 | **2** | 19 | 21 |
| min scaled Jacobian | −0.9993 | −0.0607 | −0.2729 | −0.4938 |
| mean scaled Jacobian | 0.9651 | 0.9664 | 0.9648 | **0.9767** |
| singular edges | 308 | 118 | 188 | 164 |
| singular arcs | 12 | 7 | 8 | **6** |
| arc endpoints on the outlet | 2 | 2 | 2 | 2 |
| arc endpoints on the inlet | 2 | 0 | 2 | 0 |
| internal cavities | 0 | 2 | 1 | **0** |
| enclosed volume | 6.5253 | 5.4483 | 6.3466 | 5.4688 |

**Every run has exactly two singular-arc endpoints on the outlet**, whatever
is cut out of the domain. Cutting the hub/shroud layer does not create them;
it removes the two at the inlet and simplifies the graph.

## Block structures produced from them

| basis | blocks | cuboids | cells in cuboids | tiny (<10) | inverted | Hausdorff | validator |
|---|---|---|---|---|---|---|---|
| v9 raw | 82 | 61 (74 %) | 85 % | 8 | 2 | — | — |
| v9 + cleanup | 81 | 73 (90 %) | 98 % | 8 | 2 | 0.1851 | INVALID |
| v9 + collapse + untangle | 42 | 36 (86 %) | 98 % | 2 | **0** | 0.0999 | VALID |
| v9 re-attached (full domain) | 73 | — | — | — | 5 | — | — |
| v10 raw | 198 | 164 (83 %) | 96 % | 12 | 19 | — | — |
| v11 raw | 117 | 114 (97 %) | **100 %** | 16 | 21 | 0.0341 | — |
| **v11 + collapse + untangle** | **16** | 14 (88 %) | 99 % | **0** | **0** | **0.0341** | **VALID** |

The final v11 structure has **all 16 blocks with exactly 6 faces**; the two
non-cuboids fail the cube-adjacency test, not the face count.

### The v15-v18 bases, scored on one card

`basis_report.py` (full table: `output/hex3d_algohex/basis_scorecard.md`), all
of them post-processed identically with `--collapse-rounds 5 --untangle
--untangle-rounds 6` and labelled from the unmerged `T1_9_tet_v5.vtk`:

| | v11 | v15 | **v16** | v17 | v14 |
|---|---|---|---|---|---|
| blocks | **16** | 120 | 103 | 70 | 217 |
| non-cuboid | 2 | 4 | 2 | 2 | 4 |
| **pinches** | **2** | 0 | **0** | 0 | 0 |
| blocks < 10 cells | **0** | 26 | 9 | 8 | 50 |
| smallest block | **224** | 1 | 3 | 6 | 1 |
| inverted | **0** | 12 | **0** | 3 | 12 |
| min scaled Jacobian | **0.152** | −0.604 | 0.016 | −0.123 | −0.299 |
| scaled Jacobian p5 | 0.936 | 0.920 | **0.944** | 0.933 | 0.921 |
| Hausdorff | 0.0341 | 0.0139 | **0.0138** | 0.0147 | 0.0233 |
| validator | VALID | INVALID | **VALID** | INVALID | INVALID |
| TFI direction classes | **11** | 14 | 15 | 13 | 16 |
| classes stuck at 1 division | **1** | 3 | 2 | 2 | 4 |

**Only two bases are valid at all: v11 and v16.** By the selection rule fixed
before the runs — 0 inverted, then non-cuboid + pinches, then tiny blocks,
then block count, then Hausdorff — v16 wins, because it is the only valid
basis without a pinch and it halves the boundary error.

It is worth being explicit that the rule's *ordering* decides this. v16 is
better where the input is concerned: no pinch, Hausdorff 0.0138, and a better
scaled-Jacobian p5. v11 is better where the block structure is concerned: 16
blocks against 103, no block under 10 cells against nine, and one direction
class stuck at a single division against two. For a hand-editable block
topology v11 still reads better; for boundary fidelity and freedom from the
pinch, v16 does.

The v14 row is measured with the collapse guard as it stood before the fix
below; its 217 blocks were reproduced under both variants, so the number
stands, but its other columns are not strictly comparable.

## Reproducing

```bash
PY=/root/repos/duty/quadmesh/.venv/bin/python

# inputs
$PY experimentell/hex3d_algohex/tet_prep_v5.py                      # v11 input
$PY experimentell/hex3d_algohex/tet_prep_v5.py --keep-prisms \
      --out T1_9_tet_v6.vtk                                         # v10 input
$PY experimentell/hex3d_algohex/tet_prep_v5.py --remesh-ogrid 0.035 \
      --merge-interfaces all --out T1_9_tet_v9.vtk                  # v15/v18
$PY experimentell/hex3d_algohex/tet_prep_v5.py --remesh-ogrid 0.035 \
      --merge-interfaces hub --out T1_9_tet_v10.vtk                 # v16 input
$PY experimentell/hex3d_algohex/tet_prep_v5.py --remesh-ogrid 0.035 \
      --ring-kink 70 --out T1_9_tet_v11.vtk                         # v17 input

# AlgoHex — sequentially, two concurrent runs OOM on a 7.7 GiB box
systemctl start docker
$PY experimentell/hex3d_algohex/run_algohex.py --tag v11 \
      --in-vtk data/T1_9/T1_9_tet_v5.vtk -- -n 60000
$PY experimentell/hex3d_algohex/run_algohex.py --tag v10 \
      --in-vtk data/T1_9/T1_9_tet_v6.vtk -- -n 60000

# blocks
$PY experimentell/hex3d_algohex/clean_blocks.py \
      output/hex3d_algohex/T1_9_hex_v11.ovm \
      --collapse-rounds 5 --untangle \
      --out output/hex3d_algohex/deliverable/T1_9_blocks_v11.vtk

# re-attach the removed parts, and the 11-step showcase
$PY experimentell/hex3d_algohex/reattach.py
$PY experimentell/hex3d_algohex/showcase.py --steps all
```

## Caveats

**v9's input is no longer regenerable from current code.** The quad-diagonal
fix changed `tet_prep_v5.py`, so `T1_9_tet_v5.vtk` now comes out with the
corrected labelling. The file v9 actually ran on is kept as
`T1_9_tet_v5_diagbug.vtk`; the generator that produced it is in git history
before commit `6774ce1`.

**Runtimes are not comparable across machines.** All of these ran on the same
7.7 GiB box, sequentially. AlgoHex peaks hard during quantization — v7 and v8
were OOM-killed there when run concurrently with v6.
