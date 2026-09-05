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

`fast` and `v2_crashed` in the output directory are probe runs, not results.

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

## Reproducing

```bash
PY=/root/repos/duty/quadmesh/.venv/bin/python

# inputs
$PY experimentell/hex3d_algohex/tet_prep_v5.py                      # v11 input
$PY experimentell/hex3d_algohex/tet_prep_v5.py --keep-prisms \
      --out T1_9_tet_v6.vtk                                         # v10 input

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
