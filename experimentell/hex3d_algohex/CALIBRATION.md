# HexMe calibration sweep — P4b: 6 shapes × 3 settings

**What this is.** The hexme part of the calibration matrix (P4a harness,
`calib_sweep.py`): 6 hexme shapes spanning the families × 3 HexMeshing
settings = 18 runs at `-n 8000`, to find which feature-edge mode and which `n`
range the later sweeps should use. Runs and raw outputs live in
`/work/calib/hexme_sweep/<geometry>_<setting>_8000/`, the machine-readable
report in `/work/calib/report_hexme.json` (built and checked by
`calib_report.py`).

**Matrix.** Shapes (families: `s` = simple/structured, `i` = industrial,
`n` = n-simplex/transition; suffix `c`/`u` = mesh variant):
`s01c_cube`, `s16u_torus`, `i29c_bracket`, `i31u_dlr_f6`,
`n04u_transition_prism`, `i26c_gc_turbo_vol1`. Settings:

| setting | HexMeshing flags | feature edges come from |
|---|---|---|
| A | *(none)* | the input's prescribed colors (`CellEntityIds` on the vertex/line/triangle tag cells → `vertex_colors`/`edge_colors`/`face_colors`) |
| B | `-d 70` | same colors — `-d` only feeds the dihedral fallback, which runs when the input has **no** colors at all |
| C | `-d 70 --force-feature-threshold` | recomputed from the 70° dihedral threshold, colors ignored (`--force-feature-threshold` is boolean in this build; the angle is `-d`) |

**Method.** Each run: `calib_sweep.py run --geometry … --setting … --n 8000
--out /work/calib/hexme_sweep --cores <pair> --timeout 14400`, i.e.
`run_algohex.py --backend native` → `HexMeshing -i … -o hex_hex.ovm -j
hex_hex_metrics.json` under `taskset -c <pair>`, `OMP_NUM_THREADS=2`, one run
per core pair (18 pairs from the sanctioned list; cores 22,23 and 48-71 were
left to the sibling tasks). `n=8000` because the earlier smoke at `n=2000` was
too coarse to separate settings; `--timeout 14400` because `n=8000` on the
large inputs can exceed the default 7200. Per row: status, `valid_volume`
(HexMeshing's own `valid_volume_seamless`), hex count, inverted cells
(min corner scaled Jacobian ≤ 0 over the recovered hexes — the `showcase.py` step06
metric; note `ovm_io` reorients a uniformly inside-out cell, so that case is
not counted), `n`, setting, wall time.

**Status legend** (`calib_report.py::classify`, deterministic):

| status | meaning |
|---|---|
| `works` | exit 0, hex mesh written (`hex_cells > 0`), `valid_volume ≥ 0.9` (the smoke bar) |
| `diverged` | everything else: parametrization failed/covered < 90 %, solver abort, process died on a non-SIGSEGV signal, or no output |
| `segfault` | HexMeshing died on SIGSEGV (exit −11 / 139) |
| `timeout` | killed by the harness after 14 400 s (`stopped: "timeout"`) |
| `running` | no `metrics.json` yet — the run was still executing when this report was cut; the row carries its exact command |

The **failure mode** column carries the AlgoHex failure verbatim: exit code +
last non-empty line of `hexmeshing.log`. The **note** column adds HexMeshing's
own parametrization diagnostic (`hex_hex_metrics.json`: seamless parametric
volume and invalid tet count).

## Run tables (one per geometry)

### s01c_cube

| setting | flags | n | status | valid_volume | hexes | inverted | wall_s | failure mode (exit + last log line) | note |
|---|---|---|---|---|---|---|---|---|---|
| A | default | 8000 | works | 1.000 | 8000 | 0 | 8.8 | — | parametric_volume 8e+03, n_invalid_param_tets 0 |
| B | -d 70 | 8000 | works | 1.000 | 8000 | 0 | 8.9 | — | parametric_volume 8e+03, n_invalid_param_tets 0 |
| C | -d 70 --force-feature-threshold | 8000 | works | 1.000 | 8000 | 0 | 8.8 | — | parametric_volume 8e+03, n_invalid_param_tets 0 |

### s16u_torus

| setting | flags | n | status | valid_volume | hexes | inverted | wall_s | failure mode (exit + last log line) | note |
|---|---|---|---|---|---|---|---|---|---|
| A | default | 8000 | diverged | 0.699 | 0 | — | 3605.6 | exit 0: Closed: Fri Oct 30 20:43:17 2026 | parametric_volume -1.26e-10, n_invalid_param_tets 41134 |
| B | -d 70 | 8000 | diverged | 0.699 | 0 | — | 3600.3 | exit 0: Closed: Fri Oct 30 20:43:12 2026 | parametric_volume -1.26e-10, n_invalid_param_tets 41134 |
| C | -d 70 --force-feature-threshold | 8000 | works | 1.000 | 8742 | 132 | 3491.6 | — | parametric_volume 8e+03, n_invalid_param_tets 17 |

### i29c_bracket

| setting | flags | n | status | valid_volume | hexes | inverted | wall_s | failure mode (exit + last log line) | note |
|---|---|---|---|---|---|---|---|---|---|
| A | default | 8000 | running | — | — | — | &gt;7676 (running) | running: taskset -c 12,13 /opt/venv-sci/bin/python /opt/stack/dp3d_keep_bl/.micromanager/worktrees/job-20261030-133532-3225--T4/experimentell/hex3d_algohex/run_algohex.py --backend native --prefix hex --in-vtk /work/calib/hexme/meshes/i29c_bracket.vtk --out-dir /work/calib/hexme_sweep/i29c_bracket_A_8000 -- -n 8000 | — |
| B | -d 70 | 8000 | running | — | — | — | &gt;7676 (running) | running: taskset -c 14,15 /opt/venv-sci/bin/python /opt/stack/dp3d_keep_bl/.micromanager/worktrees/job-20261030-133532-3225--T4/experimentell/hex3d_algohex/run_algohex.py --backend native --prefix hex --in-vtk /work/calib/hexme/meshes/i29c_bracket.vtk --out-dir /work/calib/hexme_sweep/i29c_bracket_B_8000 -- -d 70 -n 8000 | — |
| C | -d 70 --force-feature-threshold | 8000 | running | — | — | — | &gt;7676 (running) | running: taskset -c 16,17 /opt/venv-sci/bin/python /opt/stack/dp3d_keep_bl/.micromanager/worktrees/job-20261030-133532-3225--T4/experimentell/hex3d_algohex/run_algohex.py --backend native --prefix hex --in-vtk /work/calib/hexme/meshes/i29c_bracket.vtk --out-dir /work/calib/hexme_sweep/i29c_bracket_C_8000 -- -d 70 --force-feature-threshold -n 8000 | — |

### i31u_dlr_f6

| setting | flags | n | status | valid_volume | hexes | inverted | wall_s | failure mode (exit + last log line) | note |
|---|---|---|---|---|---|---|---|---|---|
| A | default | 8000 | running | — | — | — | &gt;7677 (running) | running: taskset -c 18,19 /opt/venv-sci/bin/python /opt/stack/dp3d_keep_bl/.micromanager/worktrees/job-20261030-133532-3225--T4/experimentell/hex3d_algohex/run_algohex.py --backend native --prefix hex --in-vtk /work/calib/hexme/meshes/i31u_dlr_f6.vtk --out-dir /work/calib/hexme_sweep/i31u_dlr_f6_A_8000 -- -n 8000 | — |
| B | -d 70 | 8000 | running | — | — | — | &gt;7677 (running) | running: taskset -c 20,21 /opt/venv-sci/bin/python /opt/stack/dp3d_keep_bl/.micromanager/worktrees/job-20261030-133532-3225--T4/experimentell/hex3d_algohex/run_algohex.py --backend native --prefix hex --in-vtk /work/calib/hexme/meshes/i31u_dlr_f6.vtk --out-dir /work/calib/hexme_sweep/i31u_dlr_f6_B_8000 -- -d 70 -n 8000 | — |
| C | -d 70 --force-feature-threshold | 8000 | running | — | — | — | &gt;7677 (running) | running: taskset -c 24,25 /opt/venv-sci/bin/python /opt/stack/dp3d_keep_bl/.micromanager/worktrees/job-20261030-133532-3225--T4/experimentell/hex3d_algohex/run_algohex.py --backend native --prefix hex --in-vtk /work/calib/hexme/meshes/i31u_dlr_f6.vtk --out-dir /work/calib/hexme_sweep/i31u_dlr_f6_C_8000 -- -d 70 --force-feature-threshold -n 8000 | — |

### n04u_transition_prism

| setting | flags | n | status | valid_volume | hexes | inverted | wall_s | failure mode (exit + last log line) | note |
|---|---|---|---|---|---|---|---|---|---|
| A | default | 8000 | works | 1.000 | 8100 | 102 | 1381.9 | — | parametric_volume 8e+03, n_invalid_param_tets 0 |
| B | -d 70 | 8000 | works | 1.000 | 8100 | 102 | 1372.8 | — | parametric_volume 8e+03, n_invalid_param_tets 0 |
| C | -d 70 --force-feature-threshold | 8000 | works | 1.000 | 11434 | 104 | 1457.3 | — | parametric_volume 8e+03, n_invalid_param_tets 0 |

### i26c_gc_turbo_vol1

| setting | flags | n | status | valid_volume | hexes | inverted | wall_s | failure mode (exit + last log line) | note |
|---|---|---|---|---|---|---|---|---|---|
| A | default | 8000 | running | — | — | — | &gt;7676 (running) | running: taskset -c 32,33 /opt/venv-sci/bin/python /opt/stack/dp3d_keep_bl/.micromanager/worktrees/job-20261030-133532-3225--T4/experimentell/hex3d_algohex/run_algohex.py --backend native --prefix hex --in-vtk /work/calib/hexme/meshes/i26c_gc_turbo_vol1.vtk --out-dir /work/calib/hexme_sweep/i26c_gc_turbo_vol1_A_8000 -- -n 8000 | — |
| B | -d 70 | 8000 | running | — | — | — | &gt;7676 (running) | running: taskset -c 34,35 /opt/venv-sci/bin/python /opt/stack/dp3d_keep_bl/.micromanager/worktrees/job-20261030-133532-3225--T4/experimentell/hex3d_algohex/run_algohex.py --backend native --prefix hex --in-vtk /work/calib/hexme/meshes/i26c_gc_turbo_vol1.vtk --out-dir /work/calib/hexme_sweep/i26c_gc_turbo_vol1_B_8000 -- -d 70 -n 8000 | — |
| C | -d 70 --force-feature-threshold | 8000 | running | — | — | — | &gt;7676 (running) | running: taskset -c 36,37 /opt/venv-sci/bin/python /opt/stack/dp3d_keep_bl/.micromanager/worktrees/job-20261030-133532-3225--T4/experimentell/hex3d_algohex/run_algohex.py --backend native --prefix hex --in-vtk /work/calib/hexme/meshes/i26c_gc_turbo_vol1.vtk --out-dir /work/calib/hexme_sweep/i26c_gc_turbo_vol1_C_8000 -- -d 70 --force-feature-threshold -n 8000 | — |


## Synthesis

**Feature mode: the prescribed cell colors are not reliable; the forced 70°
dihedral threshold is.**

- A (default) and B (`-d 70`) are the *same* mode on every input here: all six
  hexme VTKs carry their own feature tags (`CellEntityIds` on the vertex/line/
  triangle cells), so HexMeshing logs "Found prescribed feature colors (from
  vtk file)" and the `-d` angle is never used — the dihedral detector is only
  the fallback for inputs with no tags at all. Every completed A row and its B
  twin agree exactly (same valid_volume, hex count and inverted count); the
  sweep cannot separate A from B because there is nothing to separate.
- C (`-d 70 --force-feature-threshold`) ignores the tags instead: "No
  prescribed feature tags found. Set boundary faces to feature faces and edges
  w.r.t. dihedral angle threshold 70 to feature edges." That is a genuinely
  different feature network, not a threshold tweak of A.
- On `s01c_cube` both modes work (8000 hexes, valid 1.000, 0 inverted). On
  `n04u_transition_prism` both work too, but C at the same `-n 8000` delivers
  11434 hexes against 8100 for A/B — forcing the threshold adds feature
  constraints that refine ~40 % past the target.
- On `s16u_torus` the prescribed colors **fail**: A and B exit 0 but the
  parametrization collapses (valid_volume 0.699, seamless parametric volume
  ≈ 0, 41134 invalid tets, 0 hexes), while C meshes it (valid 1.000,
  8742 hexes, 132 inverted). The torus's sparse tagged-cell set
  (158 tagged cells) is the wrong constraint network for AlgoHex; the 70°
  dihedral detector recovers the real creases.
- Conclusion: cell colors win nothing and lose the torus; the forced threshold
  meshed every shape it finished. For third-party/hexme-style meshes default
  to C; use A only on our own `tet_prep_v5` output, whose tags are curated.

**n ranges.**

- `n=8000` is the working point for this calibration: it hit 8000 exactly on
  the cube, 8100 on the transition prism (A/B) and 8742 on the torus
  (C). The `n=2000` smoke under-resolved the same shapes (the cube gave
  2197 = 13³), which is why this sweep is `n=8000`.
- In C mode `-n` is a *target*, not a budget: the extra feature constraints
  overshoot it (+41 % on n04u: 11434 vs 8100). Expect a higher delivered count
  on feature-dense shapes, or lower `-n` when a budget matters.
- The three large shapes (`i29c_bracket`, `i31u_dlr_f6`,
  `i26c_gc_turbo_vol1`) had no `metrics.json` when this table was cut
  (≈2 h in); their rows are `running` with the exact command. `i31u_dlr_f6`
  (19.1 M input tets) is the timeout risk at the 14 400 s cap.

**`--keep-prisms` is a no-op on this whole sweep.** `--keep-prisms` is
`tet_prep_v5.py`'s switch to keep the prismatic boundary layer in the T1_9
blade chain. It cannot matter here: the hexme inputs are pre-meshed mixed VTKs
that feed HexMeshing directly (no `tet_prep*` step in the hexme path at all),
and none of the six carries a single prism cell — every input has only vertex
(1), line (3), triangle (5) and tet (10) cells:

| input | verts | lines | tris | tets | prisms |
|---|---|---|---|---|---|
| s01c_cube | 8 | 60 | 412 | 945 | 0 |
| s16u_torus | 1 | 190 | 10 756 | 43 606 | 0 |
| i29c_bracket | 585 | 8 223 | 114 302 | 207 558 | 0 |
| i31u_dlr_f6 | 227 | 6 741 | 533 284 | 19 118 449 | 0 |
| n04u_transition_prism | 6 | 184 | 3 432 | 7 604 | 0 |
| i26c_gc_turbo_vol1 | 441 | 9 181 | 109 604 | 390 673 | 0 |

(`n04u_transition_prism` despite its name carries no prism cells either —
"transition_prism" is the shape of the domain, not a cell type.) So the flag is
a no-op for all 18 rows; where it *does* act is the keep-BL T1_9 chain, which
is calibrated separately.

**Reproduce / re-check:**

```
cd experimentell/hex3d_algohex
bash /work/calib/hexme_sweep/run_all.sh        # 18 runs, 2 cores/pair, timeout 14400
/opt/venv-sci/bin/python calib_report.py --build /work/calib/hexme_sweep \
    --out /work/calib/report_hexme.json        # + print the tables above
/opt/venv-sci/bin/python calib_report.py --check /work/calib/report_hexme.json --rows 18
```

**Where this is used.** HANDOFF.md carries the handoff: the sweep matrix, the
status vocabulary and the C-mode recommendation are the inputs for the charite
and keep-BL calibration sweeps and for the `n` choice in the dataset runs
(`DATASET_PIPELINE.md`).

_Report built 2026-10-30 20:51 UTC from `/work/calib/report_hexme.json` (18 rows)._
