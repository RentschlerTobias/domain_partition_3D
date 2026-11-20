# hexdp task list with progress

Status: `[ ]` open, `[~]` in progress / partly done, `[x]` done, `[?]` waiting for the owner's decision.
Update the status, the date and a one-line note when something changes. Context and rules: `docs/HANDOFF_2026-11-19.md`.
Never start long computations without the owner's go (see the rules there).

## A. Running / to evaluate (state 2026-11-19 evening)
- [~] A1 HexMe agent test A (10 cases, solid + fluid, free method choice) - **evaluated 2026-11-19** (9 of 10;
      the test was aborted, i11b_s1 fluid never started). Report: `/opt/stack/cases/hexdp_agent/HEXME_TEST_A.md`.
      Result 1/9 valid (s01b cube solid only, 1 block). Main findings: only i11b_s1 ran with `--relative`; 3 of the
      5 "raw mesh off surface" rejects would pass `max_raw_rel 0.05` (s16b fluid, i14b solid, i14b fluid -> rerun
      needed); 3 AlgoHex IGM failures (invalid seamless map), not reproducible with an identical tet mesh; open
      cavities (52-417 quads) and wrong boundary topology (i11b_s1, genus 3: Euler characteristic); polycube broken
      (`libgco.so` not on the library path); beam post-check crashes with a ValueError instead of reporting INVALID;
      `size_max` / `remesh_all` silently inert on a VTK input. Open: fall 10 (i11b_s1 fluid) - needs the owner's go.
- [x] A2 canada runner p_0018 (2026-11-19): VALID 18 blocks, min sJ 0.71 with --max-raw-hausdorff 0.10 (relative)
      (`/work/scratch/canada_snap/beam2.log`). Snap discarded: Hausdorff 0.203 (5.7 % of diag) -> 0.216 at step 0.75.
      Suspected: nearest-face labels wrong when the mesh is that far off -> see B9. Next: rerun ru 0016, 0017
      (needs go; ~1 h each on 2 cores).
- [x] A3 radMeridional variation test (2026-11-19, runner + cone only, guide vane fixed, lengths +-L, rest +-L/2):
      +-10 % 4/4, +-25 % 3/4, +-50 % 2/4. Failures: v25_0, v50_1 dtOO `!md.converged()` in the meridional contour after
      2 s (strong changes of l_hub*/l_shroud*); v50_2 hit the 1 h timeout (hang). ~21 min per machine on 2 cores.
      Factors per sample: `/work/scratch/radmer/vartest/*.factors.json`. Suggestion for C3: +-25 % with a 40 min
      timeout and filtering.
- [ ] A4 Torus fluid (HexMe s16b) beam with relative tolerance + snap: rerun with 4 workers, 3 h limit
      (`/work/scratch/snap_torus`, earlier run hit the 1 h limit in the search).
- [~] A5 polycube batch, the 16 smallest genus-0 HexMe models (2026-11-19, `/work/scratch/polycube_batch/`):
      9 of 16 produced blocks, 5 of those valid (s06, s08, n03, s07, n06; 1-9 blocks), 4 invalid with a boundary
      defect: the block boundary is no closed 2-manifold, chi = 3 (odd) or chi = 6 in 3 components where the
      geometry has chi = 2, 1 - the same defect as i11b_s1 under AlgoHex, plus 2-8 `t_junction` warnings. 5 models
      (n08, n09, n14, n04, n13) got `CELLS 0` from `polycube_withHexEx` (no hex mesh at all) and the beam crashed on
      `np.vstack([])` (beam_collapse.py:377) instead of reporting it; 2 (n15, s14) hit the absolute Hausdorff gate.
      Genus of all 57 solid models on record: 28 are genus 0, 29 have 1-11 handles (i01 g=9, i07 g=11).
      After the B12 fix (`f10774b`) rerun in `/work/scratch/polycube_batch2/`: **9 of 9 models that produce blocks are
      valid** (s10 4->7, n01 3->5, n07 18->9, n02 7->7, the rest unchanged), all `check(geo)` clean. Still 7 failed
      for other reasons (5 `CELLS 0`, 2 Hausdorff gate). Caveat: two runs of the same model differ before the fix
      already (n07: 146 vs 144 fine cells, min sJ 0.014 vs 0.333), so polycube/HexEx is not deterministic - check
      that before a dataset relies on it (cause and fix options: B15). n07 keeps min sJ 2.7e-6: valid, but nearly
      degenerate.
      Rerun 2026-11-20 with the B15 patch (best archive individual instead of the worst), `/work/scratch/polycube_batch5/`:
      **11 of 16 valid** (n01b 1, n02b 5, n03b 3, n06b 1, n07b 9, n08b 1, n13b 1, s06b 1, s07b 7, s08b 1, s10b 7
      blocks), 3 still `CELLS 0` (n04b, n09b, n14b), 2 invalid (`s14b` inverted behind the sphere cutout, `n15b`
      min sJ 2.2e-6 -> one cell at -7e-4). Quality clearly better: n07b 0.64, n03b 0.86, n01b 0.86, s08b 0.81
      (n07b was 2.7e-6). n08b (4e-7) and n13b (0.056) are valid but useless - a min-sJ floor is still needed.
      Those numbers were single random draws (8-thread labeling). Deterministic rerun 2026-11-20 (`label_threads: 1`,
      `min_sj: 0.05`), `/work/scratch/polycube_batch6/`: **8 of 16 valid** (n01b 1, n02b 5, n03b 3, n06b 1, n07b 2,
      s06b 1, s07b 7, s10b 7 blocks; min sJ 0.08-0.996), 4 `CELLS 0` (n04b, n08b, n13b, n14b), 3 raw-mesh-off-surface
      (n09b, n15b, s08b, exit 3), 1 inverted (s14b). So the honest single-run yield is 8/16, and the 11/16 above was
      partly luck of the thread schedule - the deterministic labeling of n08b uses 2 directions and yields no hexes.
      Seeded multi-start (2026-11-20, `label_seeds`, default 8): `evolabel_cli` takes a seed and prints one
      `FINAL_SCORE` line; the strategy runs seeds 1..K single-threaded, keeps the best score and restores that
      labeling. Deterministic, ~2 s per seed. Effect on the 8 failures (`/work/scratch/polycube_batch7/`): s08b
      recovered (1 block, min sJ 0.78), n15b now fails the sJ floor (9.2e-6) instead of inverting, s14b unchanged,
      and the four `CELLS 0` are unchanged - on n08b all 8 seeds end with `invalid 6, flip 55-84`, so the labeling
      is never valid and those small sharp models are simply out of polycube's reach (route them to algohex).
      Deterministic yield: 9 of 16.

## B. Methods (next work, owner said: relative Hausdorff + snap first [done], then canada pilot [done], C4Hex ok)
- [x] B1 relative Hausdorff tolerances (c33c13b); now the CLI default (`--absolute` keeps the old numbers 0.05/0.15,
      relative defaults 0.02/0.05 of the bbox diagonal; `tolerances()` + `test_beam_tolerances.py`). Effect on the
      two gate rejects of A5: n15b icosahedron now valid (6 blocks, min sJ 0.21, 34 s); s14b passes the gate and then
      fails with `inverted_cells` (scaled Jacobian 0 in blocks 0/2/3) - a real defect, not the gate.
- [x] B2 snap boundary to surface, stepwise (c33c13b, f3be635)
- [x] B3 blocks around holes / non-cuboid blocks cut automatically (ab959b2, c33c13b)
- [~] B4 C4HexMeshing as alternative quantization/extraction after AlgoHex: built (`/opt/stack/ext/C4HexMeshing`).
      Stage written: `src/hexdp/strategies/algohex/extract_c4hex.py` (runs `c4hex_cli` on a `.hexex`, finds
      `<prefix>_<n>_hex.ovm`, reports failures with exit code and log tail; 5 tests incl. a real run on
      `minimal.hexex`: 10 648 hexes, min sJ 1.0, readable by `hexdp.io.ovm`). AlgoHex writes the map with
      `--checkpoints` (new param `checkpoints` of the algohex strategy). Measured on s16b_torus solid (2026-11-20):
      HexMeshing produced the seamless map (37 MB) but HexEx still found 22185 flipped tets; `c4hex_cli` then stops
      after 2 s with "Reading seamless map failed with error code 5" = `INVALID_CHART` (Reader.hpp: "parametrically
      inverted or degenerate tet"), also with `--force-sanitization`. So C4Hex does not rescue a folded
      *parametrisation*: it replaces HexEx's IGM, not AlgoHex's field. It helps only for maps that are valid but
      quantise badly. Not wired into the strategy yet (`extractor: algohex|c4hex` still open); the 3 AlgoHex gate
      cases stay blocked unless the field improves. The C4Hex repo's own `fandisk_algohex.hexex` sample fails too
      (exit 22 after 17 min), so the bundled samples are not all solvable either.
- [ ] B5 algohex strategy: option `tet_prep: v5|named` (generic tet_prep_named for any named .msh, e.g. radMeridional,
      canada) instead of the tistos-specific v5.
- [ ] B6 mesh-based sweep detection (extrusion / revolution parts from tet meshes; many HexMe parts)
- [ ] B7 polycube with through holes (cuts before labeling; HexDom BSD-3 is a candidate)
- [ ] B8 body O-grid for external-flow cases (body in box), like the blade O-grid
- [ ] B9 snap: (a) when the raw mesh is far off (> ~3 % of the diagonal), the nearest-face surface labels of the
      boundary quads can be wrong -> label by feature-curve-aware matching before projecting (canada ru p_0018:
      snap made it worse); (b) if steps still add inversions, smoothing (Winslow / Mesquite) after projection.
      Also consider raising algohex max_raw_rel 0.05 -> 0.10 (p_0018 valid with 0.10).
- [ ] B10 runner (structural) robustness: 1452/2000 valid; failures in hub layout QP, near-zero sJ
- [ ] B11 save AlgoHex frame fields in future datagen (`--cell-quaternion-out-path` etc.)
- [x] B12 the block boundary must stay a closed 2-manifold (fixed 2026-11-19, `clean_blocks.resolve_hanging_faces`):
      `cut_non_cuboids` can cut a face sheet that ends on a valence-3 singular edge of the polycube base complex, so
      the cut stops inside a block face (the neighbour shares 9 of 18 fine quads, a 1:2 hanging face) and in the
      coarse model the two blocks touch only along an edge -> a pinch, chi 3, t_junction warnings. Not the beam
      search and not the snap (all flags give the identical result). Fix: label the fine quads of every block face
      with their neighbour, cut the block on the fine index planes where the neighbour changes, iterate to a
      fixpoint (cascade), guards `max_blocks`/`min_cells` roll the whole pass back. s10b 4->7 blocks chi 3->2,
      n01b 3->5 chi 6/3->2/1, n07b 18->26, all `check(geo)` clean; n02b rolls back (dense 1:n pattern) and stays
      invalid. 143 unit tests pass. Cost: the block count goes up (4->7, 3->5, 18->26), which the dataset metric
      must accept. Open: n02b; `tistos_0001` had `nonmanifold_face` before the fix - rerun to see if it is the same
      cause.
- [x] B13 (2026-11-19): an empty hex mesh is an error with a message. `hexahedra()` in beam_collapse prints
      "no hexahedra" and exits 2 instead of dying in `np.vstack([])`; the polycube strategy raises a RuntimeError
      that names the labeling, e.g. "2 of 6 directions: [2, 4], invalid labels [0]", and points to algohex. Not a
      `hex_scale` problem (2/4/8 all give 0 cells): the evocube labeling itself is degenerate and differs between
      runs ([1,2,4] / [3,4] / [2,4]), which is more evidence for the non-determinism note in A5. Affects 5 models.
- [~] B14 curved block edges (goal of the dataset: 8 corners + curved edges). Finding: the old tistos pipeline
      never fitted edges to feature curves - `edge_ctrl` (2 points per edge) are the absolute control points of a
      cubic Bezier least-squares fitted to the raw fine-node polyline with fixed endpoints (54 cubic / 6 quadratic /
      24 chord in sample_16c), `edge_polyline` is that raw chain (2-7 nodes at n=2000, not resampled, not projected).
      Faces/volume come from `tfi.py` (cubic face resampling + surface projection + TFI + per-block fallback to
      linear) and never read `edge_ctrl`; real feature curves exist only in the optional `curve_postprocess.py`
      (block edge between two different patches -> 65 arclength points on the seam, corners may move) and
      `blade_insert.py` (only the HEX20 mid-node). The old "curved faces + TFI" attempt failed: 120 and 185 inverted
      cells. Implemented 2026-11-20: `dataset.block_edges(bs, k, labeller)` resamples every boundary block edge along
      its chord and projects the interior samples onto the input surface (endpoints fixed, a move beyond half the
      chord is refused); the exporter builds the labeller from `meta["tet_vtk"]` and records
      `edge_off_surface_max` (inner points) and `corner_off_surface_max` (the block corners) over boundary edges
      only. A projection is refused when a point would move more than half the chord or when the projected interior
      points land on different patches - that is a crease or a rim, where the nearest point zigzags between two
      faces and the straight chord is the better answer until seam curves exist. Measured (inner off-surface, per
      boundary edge): s10b 47 of 54 edges improved, median 0.179 -> 0.000 (worst 0.562 -> 0.181); n01b 8 of 12,
      median 0.041 -> 0.000; n07b 2 of 20 (its block boundary sits 0.48 off the surface and the edges are creases);
      s07b 6 of 50 (median already 0). The old numbers are the linear curves and they explain the "square rim":
      s10b's rim edges have 2 nodes, so interpolation alone can never bend them. The refined export curve on n01b
      has chord 0.324 / arclength 0.439, i.e. a real arc.
      Open: the s10b residual is the *corners* (0.09-0.18 off the surface, the snap was discarded - B9), not the
      edges; the fine mesh and the sJ/Hausdorff metrics stay linear (only the exported target is curved); and the
      tistos-compatible 2-control-point Bezier form is not written yet. Next for feature edges: an edge whose two
      boundary quads lie on different patches is a crease or a rim (s10b's cylinder rim), and nearest-point
      projection may zigzag between the two faces - such edges need the seam curve (patch-pair matching, as
      `curve_postprocess.py` does it) instead of the projection. Untested so far.
- [~] B15 polycube labeling: non-determinism and the worst-individual bug (2026-11-20). Two separate causes, both in
      `evolabel_cli`, found by reading `/opt/stack/ext/evocube`:
      (a) *Non-determinism.* `srand(1)` is fixed, but the 100 mutations of a generation run in an OpenMP parallel
      for (`app/evolabel_cli.cpp:131`) and all draw from the global `rand()` (`include/archive.h:116`,
      `evolabel_cli.cpp:75,142`, ~20 places in `include/labeling_individual.h`), so the scheduling decides which
      individual gets which number. Measured on n08b: two runs with `OMP_NUM_THREADS=1` give an identical
      `labeling.txt` (md5 6b5c8efa), with 8 threads the second differs (04e14632). No seed/thread option exists.
      (b) *Quality.* The archive is sorted by descending score (front = worst, `bestScore()` = `back()`, the
      roulette selection favours high indices), yet `app/evolabel_cli.cpp:241` (and upstream `app/evolabel.cpp:244`)
      takes `archive.getIndiv(0)` - the worst stored individual - as the final labeling. Patched locally to
      `bestIndividual()`; the build cache of `/opt/stack/ext/evocube/build` still points to the pre-refactor path
      `/work/ext/evocube`, so a rebuild needs `ln -s /opt/stack/ext/evocube /work/ext/evocube`.
      Measured after the patch on the five `CELLS 0` models: n08b and n13b now produce 1 block each (n08b min sJ
[1208 characters, dropped: the file on disk has it]

- [x] B16 (2026-11-20): the structural check is not a quality check. `beam_collapse --min-sj F` (default 0 = off,
      relative tolerances untouched) rejects a structure whose worst fine cell is below F: issue `low_jacobian`,
      exit 4, so `report.json` says invalid like any other defect. Both strategies carry `min_sj` (default 0.05).
      Reason: `valid` alone accepted n07 at min sJ 2.7e-6 and n08b at 4e-7. Tests: `test_quality_gate.py`.

## C. Data / machines
- [x] C1 HexMe b models split into solid / fluid (57 free models), agent data dir
- [x] C2 radMeridional export stage (9ea214b) + viz (`/opt/stack/viz/radmeridional/`)
- [?] C3 radMeridional dataset: choose parameter ranges (owner: "test it"; results in A3), include guide vane and
      blade counts? Then datagen (needs go). Geometry is expensive (~20 min / machine).
- [?] C4 canada datagen: shrink the box (gv_M half range), filter dtOO failures (timeout 900 s), gv works; ru needs
      A2/B4. Needs go.
- [x] C5 canada pilot evaluated (`cases/canada/pilot/PILOT_EVAL_2026-11-19.md`)
- [ ] C6 HexMe datagen for both cases (solid via polycube/algohex, fluid via algohex) after A1 / B4 (needs go)
- [ ] C7 tistos guide vane / diffuser (a colleague of the owner knows where they are)
- [ ] C8 more data sources later (MAMBO, Thingi10K licences not checked) - owner: "later"

## D. Housekeeping
- [ ] D1 `/work/hexdp_agent`, `/work/hexdp_agent_refs` unused now (MCP points to NFS) - delete when the owner says so
- [ ] D2 push (owner pushes): not pushed yet ab959b2, c33c13b, f3be635 + the handoff commit (`git log origin/main..`)
- [ ] D3 SciPy vs SLEPc eigensolver comparison for the FEniCS modal analysis (offered, not requested)
- [ ] D4 remaining datagen samples 0993 / 0309 (keep-BL)
- [x] D5 fixes found by A1 (2026-11-19): evocube library path is derived from the build (`_lib_dirs`; prep, label
      and hexex run with `LD_LIBRARY_PATH`), beam post-check reports `not_structured` as INVALID (exit 4) instead
      of a `ValueError`, `size_max`/`remesh_all` warn when tet_prep is skipped (tagged VTK input), `check_environment`
      honours `$HEXMESHING_BIN` and names the PATH fallback, `report.json` of a failed run keeps the stage timings.
      Tests: 134 passed, 3 failed (`test_mcp_server`, identical on HEAD: hexdp is not pip-installed in `/opt/venv-sci`,
      so the entry points are missing). Evocube library path proven by a smoke test on HexMe `s01b_cube` solid
      (empty `LD_LIBRARY_PATH`): 48.6 s, 1 block, valid, min sJ 0.90 - polycube works now. Open: set
      `HEXMESHING_BIN` in the MCP config (config, not code); an install needs the owner's ok.
