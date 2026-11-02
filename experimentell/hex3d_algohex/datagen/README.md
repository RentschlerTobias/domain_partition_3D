# New tistos samples: dtOO -> AlgoHex -> beam collapse (2026-10-06)

Generates new training samples for the block generator from fresh design vectors.
One AlgoHex n per geometry (n=2000): different n only change the initial AlgoHex
topology; after the beam collapse they end in the same canonical structure.
The samples contain only the core blocks; the blade is inserted after generation
(docs/decisions/2026-10-06-blade-inserted-after-generation.md). The dtOO `case/`
directory (exact blade geometry) is kept per sample for that step.

| file | purpose |
|---|---|
| `draw_new.py <root> <n> <seed> <prefix>` | scrambled Sobol in the 30-D tistos box (bounds from /dtOO/demo/tistos/build.py) -> `<root>/<prefix>_<i>/params.json`; prints the distance to the existing dataset geometries (`DATASET_BATCH`) |
| `gen_mesh.py` | dtOO build of one design vector -> case/mesh.msh (copied into the case by run_sample.sh) |
| `run_sample.sh <run_dir>` | dtOO -> tet_prep_v5 -> AlgoHex -> beam_collapse -> tfi -> export_sample (sample.npz); skips finished runs |
| `run_lane.sh <cores>` | one lane over `$DATAGEN_WORK/runs/*`, atomic claims in `locks/`, stops on `STOP` |
| `start_lanes.sh` | one lane per core pair (`LANES`, default `"0,1 2,3 4,5"`), log `lanes.log` |

    python draw_new.py $DATAGEN_WORK/runs 400 20261006 tistos_s1
    DATAGEN_WORK=<dir> LANES="0,1 2,3 4,5" ./start_lanes.sh      # stop: touch $DATAGEN_WORK/STOP

The audit (../relabel/audit.py) and the dataset builder of meshtron
(scripts/blockgen_v2/build_dataset.py) read the same sample.npz layout.

## Measured (courier, 2 own cores per sample)
dtOO ~110 s, tet_prep ~5 s, AlgoHex 520-1450 s (mean ~16 min, 2.0-2.5 GB RAM), beam 190-290 s,
tfi + export ~2 s -> ~20 min per sample, ~9 samples/h on 6 cores, ~3 per hour per 2-core lane.
First batch (seed 20261006, prefix tistos_s1, 400 drawn, interrupted after 11): all 11 finished
samples 12 blocks / 12 cuboids; 1 AlgoHex segfault (exit -11 after parametrisation, ~10 % expected).
The new design vectors are as far from the 378 existing geometries (median 1.53 box units) as
those are from each other (1.55).

## Kept boundary layer: run_sample_keepbl.sh + check_keepbl_sample.py (2026-10-30)

`run_sample_keepbl.sh <run_dir>` is the keep-BL variant of the pipeline
(tet_prep_v5 `--keep-prisms --remesh-all 0`: only the blade O-grid cut; hub at
r = 0.5 and shroud at r = 1.9 are real walls, their prism layers stay in the
AlgoHex domain; all surfaces and seams remeshed). Reference input for
tistos_s1_0131: /opt/stack/blade_blocking_20261021/keep_bl/tistos_s1_0131
(read-only; holds only tet.vtk - if the run dir has no params.json, copy it
from /work/datagen/runs/<name> and let the script build the case). Beam
collapse runs guard struct width 32 depth 20 hausdorff-tol 0.1 with 4 workers;
block count is checked against the 12-block S1 reference (10-16 acceptable).

Run:  taskset -c 0,1 bash experimentell/hex3d_algohex/datagen/run_sample_keepbl.sh /work/datagen_keepbl/runs/tistos_s1_0131
Check: /opt/venv-sci/bin/python experimentell/hex3d_algohex/datagen/check_keepbl_sample.py /work/datagen_keepbl/runs/tistos_s1_0131

Reference (2026-11-02): the P1 keep-BL reference of tistos_s1_0131 is built on
the **16-core AlgoHex replay** `hex_hex_x_16c.ovm` (2060 hexes, converged),
not on the diverged `hex_hex_x.ovm` run. The stages, in order:

    beam_collapse.py <run>/hex_hex_x_16c.ovm --input-vtk <run>/tet.vtk \
        --out <run>/blocks.vtk --width 32 --depth 20 --guard struct \
        --hausdorff-tol 0.1 --workers 4
    tfi.py <run>/blocks.vtk --input-vtk <run>/tet.vtk --target-h 0.05 \
        --require-lattices --solve-divisions --apply-divisions --out <run>/tfi.vtk
    export_sample.py <run>/blocks.vtk --tet <run>/tet.vtk --target-h 0.05 \
        --n 2000 --params <run>/params.json --out <run>/sample.npz --check
    curve_postprocess.py <run>/sample.npz --out-dir <run>

The 16-core mesh has a different hex connectivity than `hex_hex_x.ovm` (same
2654 point positions, different cell grouping) but collapses to the same
structure: 12 blocks / 12 cuboids, 582 cells, HexBlockValidator VALID, 0
inverted, min scaled Jacobian 0.3814, boundary Hausdorff 0.0590. TFI refills
582 -> 77372 cells at h = 0.05 (0 inverted), export_sample round-trips exact
(12 blocks, 40 vertices, 168 directed edges), curve_postprocess reports 38
feature arcs, complete, G1 fit dev/L max 0.066.

The checker verifies the five P1 acceptance points on the finished run and
writes `<run_dir>/p1_report.json` (exit 0 iff all pass). Stage wall time /
maxrss / exit code come from each `<stage>.log` `[stage]` line. The report
also names the source hex mesh (`source_hex_mesh`, `hex_hex_x_16c.ovm` for the
reference run).

1. HexBlockValidator VALID, 0 inverted cells: parsed from the validator section
   of beam.log (`min_scaled_jacobian > 0`, `scaled Jacobian: N cells <= 0` with N = 0)
   and recomputed from blocks.vtk with ovm_io.scaled_jacobian.
2. Kept prism layers in the block structure: labels 8 (hub) / 9 (shroud) are
   the wall labels of a --keep-prisms tet (tet.vtk); per surface, blocks
   touching = blocks.vtk support cells whose boundary quad's nearest input
   triangle (within 0.06 of it) carries that label (the same nearest-face rule
   clean_blocks' labeller applies); kept prisms = case/mesh.msh prisms whose
   wall tri face sits in the r <= 0.5+3e-3 (hub) / r >= 1.9-3e-3 (shroud) band;
   coverage = fraction of those prisms whose centroid lies inside a support
   cell of a touching block. Pass requires the wall labels present, kept
   prisms > 0 and at least one touching block containing BL cells.
3. O-grid outline coverage after curve_postprocess.py: min seam_coverage from
   the post_quality of `<name>_curves.npz` >= 0.99 and feature_complete True
   (the seam-repair reference measured 97.5 % -> 99.99 %).
4. Block count 10-16 (from beam.log validator `blocks: N` and sample.npz).
5. sample.npz carries the standard export_sample schema (vertices, blocks,
   quad_faces, edges, edge_ctrl, edge_polyline(+offset), dir_class(+count),
   params, quality, provenance; surface_points/tris/tri_label with --tet) and
   `<name>_curves.npz` exists.

The shroud coverage number is small (0.034 for tistos_s1_0131) and that is
geometry, not a lost layer: the 12-block structure is coarse, its outer faces
are planar quads inscribed in the shroud cylinder (r = 1.9) instead of
following it, so the msh prism centroids (a 3e-3-thick shell at
r = 1.897-1.8996) fall a few 1e-3 outside the flat faces. Measured: every one
of the 20543 "outside" shroud prisms lies within 0.0033 of the nearest support
cell's outer radius, and a tolerance sweep (tol = 1e-9 -> 1e-2) only moves the
coverage from 0.036 to 0.102, so the deficit is the flat-face offset and not a
missing layer. The shroud *is* covered: all 97 support boundary quads of the
outer surface reach r = 1.9 and 11 of the 12 blocks own shroud cells. The hub
coverage (0.999) is the honest one - the hub faces are cut by the O-grid and
stay on the cylinder.

## Pitfalls
- AlgoHex ignores OMP/OpenBLAS thread variables. Five samples sharing six cores (one xargs pool
  pinned to 0-5) gave load ~230, AlgoHex 41-91 min and beam 4-39 min. Give every sample its own
  cores (one lane per core pair).
- `gen_mesh.py` must be copied into the dtOO case before the build; run_sample.sh does that.
- dtOO and AlgoHex can hang forever (papaya 2026-10-12: 3 dtOO builds stuck 7-19 h, 4 AlgoHex
  runs 2-10 h, normal max ~7 min / ~70 min). Every stage now runs under `timeout`
  (TO_DTOO 1800 s, TO_ALGOHEX 10800 s, TO_BEAM 3600 s, others 600 s; override per env);
  a timeout logs `TIMEOUT after <t>s` before the `FAIL <stage>` line.
- P1 reference (tistos_s1_0131, 2026-11-02): built on hex_hex_x_16c.ovm; see the kept-BL section above.
