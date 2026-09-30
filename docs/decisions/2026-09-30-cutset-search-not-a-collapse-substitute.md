# Cut-set sheet search is not a substitute for the cell-level collapse

Follow-up to 2026-09-28-beam-collapse-relabelling. Measured 2026-09-30
(courier sandbox). Question: can the block reduction be done purely on the
SEPARATRIX CUT SET — drop whole labeled sheets (`block_faces.label_sheets`),
merge regions by dual adjacency — instead of running the cell-level
`beam_collapse.py` mesh edits?

## Setup

`experimentell/hex3d_algohex/sheet_beam.py` (also in /work/analysis): a
union-find partition over hexahedra merged across every face whose sheet is
not active, driven by greedy and beam (width 16) subset search over the
labeled sheets. Same census / non-cuboid / excess bookkeeping as the
pipeline's `_structure_stats`. Inputs: raw AlgoHex `.ovm` from the perturb
runs plus smoke, and two cluster samples from `blocks.vtk`.

## Aggregate results (`/work/analysis/sheet_beam_results.csv`, 15 runs x 2 methods + 3 recorded baselines)

| group | samples | raw blocks | greedy final | beam final |
|---|---|---|---|---|
| base_a / base_b (n24 sheets) | 2 | 84 (base_a), 84 (base_b) | 84 (no progress), 22 | 25 / 25 |
| r01_x, r03_x, r10_x (n24 sheets) | 11 | 84 (one outlier 60) | 22 (one stuck at 84/60) | 25 (one outlier 16) |
| smoke (18 sheets) | 1 | 22 | 22 | 22 |
| machine_0004_n2000 (blocks) | 1 | 22 | 22 | 22 |
| machine_0034_n2000 (blocks) | 1 | 12 | 12 | 12 |

Typical values over the 13 n24-sheet runs: raw 84 blocks; greedy reaches 22
(2 runs stuck at 84 with 0 drops); beam caps at 25 almost everywhere (max
25 = 84 - many-sheet state; the single 16 was on the anomalous r10_1 whose
raw complex is only 60 blocks). The 12-block canonical topology is NEVER
produced on a fresh run: sheet-level union produces at best 25 or 22.

## Key measurements (base_a)

* 24 separatrix sheets total; cutting ALL of them gives 84 blocks.
* The canonical 12-block topology (combinatorial WL hash `3c7430d6`, as
  reached by `beam_collapse.py` width 16 depth 6 guard struct on the
  same fine mesh) is NOT reachable in cut-set space.
* Greedy plateaus one step early: for base_a no single first drop reduces
  the block count, so greedy exits immediately with 84 blocks -- the
  partition only improves when MULTIPLE sheets are dropped together
  (interactions between sheets), which is exactly what beam explores;
  beam reaches 25. As a "7 drops -> 12 blocks" narrative, this was a
  target-seeking artifact: with only drops as the move set, the
  reached structures contain only 6 of the 12 cuboids
  (`clean_blocks.cuboid_status` census), not a 12-block solution.
* The mesh-quality guard machinery of the cell-level collapse is not
  replaceable either: weld semantics of the cell collapse (vertex merges
  that merge FACES beyond pure dual adjacency) exceed what "union all
  regions on one side of a sheet" can express. Cut-set merging cannot
  emulate the parallel-edge class het of `collapse_mesh_sheets`.

## Decision

1. `beam_collapse.py` stays: it is the only known path to the canonical
   12-block structure (12/12 cuboids, excess 0, min sJ 0.407,
   Hausdorff 0.0484 on the R1 base_a ring-densify run, 0 inverted).
2. Cut-set search is demoted to a cheap AUDIT / FILTER tool:
   `sheet_beam.py` runs in ~10 s per sample and tells whether the raw
   complex is stuck (e.g. greedy 84 / beam 84 = hopeless input, drop it
   before paying the 25-50 min AlgoHex rerun) -- it is not a relabeller.
3. Do not build a cut-set equivalent of the collapse for dataset
   relabelling; the 2026-09-28 relabelling plan was defined on beam_collapse
   semantics and remains there.

## Consequences / follow-ups

* Any claim of the form "dropping sheet subset S gives topology T" must be
  produced by a cell-collapse trace, not a cut-set union-find.
