# Beam-collapse relabelling: the greedy cleanup was the label noise

Grilling session, 2026-09-26 (courier sandbox), decision recorded 2026-09-28.

## Finding: topology is not a function of the geometry

The AlgoHex -> `clean_blocks` chain produces 12, 22, or 75 blocks for the
SAME geometry (machine_0004, three runs) although the singularity graph and
the raw base complex (84 blocks) are identical in all three runs. Two runs
of the identical dataset geometry have tet meshes differing by exactly 2
tets (gmsh is not deterministic) and end at 12 vs 75 blocks; the cluster
batch run of the same geometry records 22.

All of the difference is introduced by the greedy sheet collapse in
`clean_blocks.collapse_mesh_sheets`:

* it stops at the first local optimum (single-step moves only),
* its Hausdorff ceiling is relative to the raw mesh, so the ceiling value is
  an accident of the tet meshing (the decisive 84->22 sheet at base_b lay
  0.0019 above the ceiling),
* it treats an inverted intermediate cell as a blocker, although the TFI
  refills every block from scratch later.

In the batch dataset this shows as label noise: 29 combinatorial
topologies / 113 labeled topologies over 739 samples, the same geometry
resolving differently at -n 2000 vs -n 8000 in ~65 % of 361 pairs
(only 34.6 % of pairs share the labeled topology), and design parameters
carrying essentially no signal (kNN leave-one-out from the 30 parameters,
0.50-0.58 accuracy, is at or below the majority-class baseline of
0.51-0.59).

## Decision

1. The sheet collapse is optimised globally instead of greedily.
   `beam_collapse.py` does a beam search (width 16, depth 10) over collapse
   orders and keeps the best state reached anywhere in the tree, guarded by
   the `struct` guard: structure only (excess faces, block count, flat
   absolute Hausdorff tolerance 0.05); quality is a tie-break only.
2. The greedy complex is sufficient for the search input — AlgoHex fine
   meshes need not be kept or re-run. The saved `blocks.vtk` from the batch
   dataset plus `sample.npz` surface points suffice (verified: the cluster
   sample machine_0004_n2000 collapses 22 -> 12 from that input alone).
3. Dataset relabelling (planned): re-run the collapse for all 739 samples
   with the beam search and store the result beside the original
   (`sample_beam.npz`), never overwriting. Cost estimate ~18-20 h on at
   most 6-8 cores; whether the whole Sobol set lands on one canonical
   topology is the question this answers.
4. Blade insertion into the block structure is a deterministic
   post-step, not something the transformer generates: `ogrid_extrude.py`
   extrudes every core block face on `ogrid_interface` to the blade wall
   along the dtOO O-grid's own wall-normal grid lines, giving a CONFORMING
   blade O-grid by construction (26 blocks on the smoke sample, 0
   inverted cells after a joint TFI). The hub/shroud boundary layer stays
   a pure extrusion in `reattach` and is not learned either.

## Measured consequences

* 17 of 17 runs (base repeats, 3 perturbation levels of machine_0004,
  smoke sample, 2 from-blocks starts incl. the real cluster sample) end at
  12 blocks / 12 cuboids, 0 inverted, validator VALID, min scaled Jacobian
  0.03-0.53, Hausdorff <= 0.049, and all share one combinatorial topology.
* Same collapse from the greedy intermediate (`blocks.vtk`) can lose a
  little cell quality against a start from the raw AlgoHex mesh (min sJ
  0.03 vs 0.33 on one run) — block corners/edges are unaffected in
  principle, but the effect is tracked, not hidden.
* `reattach.py` additionally fixes the per-block layer fractions bug: layer
  spacing now depends on the vertex alone (cached), so blocks sharing a
  base edge agree on all layers, and the assembly welds the boundary layer
  that was previously internally disconnected (44 blocks with no shared
  vertices). Part offsets are remapped after the weld.

## Consequences / follow-ups

* `sample_one.sh` and the batch scripts should call `beam_collapse.patch`
  instead of `--collapse-rounds 5` for every future sample.
* The singularity graph is a stable learning target; block count is not.
  A GNN should predict the graph (or an integrable field), not a smooth
  frame field.
* Reported to the stack: `conditioning.py` uses `BLADE_LABEL = 5`, which is
  the hub boundary layer; the O-ring is label 7.
