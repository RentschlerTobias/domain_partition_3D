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

## B. Methods (next work, owner said: relative Hausdorff + snap first [done], then canada pilot [done], C4Hex ok)
- [x] B1 relative Hausdorff tolerances (c33c13b)
- [x] B2 snap boundary to surface, stepwise (c33c13b, f3be635)
- [x] B3 blocks around holes / non-cuboid blocks cut automatically (ab959b2, c33c13b)
- [~] B4 C4HexMeshing as alternative quantization/extraction after AlgoHex: built (`/opt/stack/ext/C4HexMeshing`).
      Next: run AlgoHex with `--sm-out-path`, feed `c4hex_cli`, convert its hex mesh to .ovm/.vtk for beam; test on a
      case where AlgoHex extraction failed (HexMe s16b_torus solid v2: 22k flipped tets in HexHex). Then wire it in
      as an option of the algohex strategy (param `extractor: algohex|c4hex`) + docs/backends.md.
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
- [ ] D5 fixes found by A1 (each minutes): evocube library path (`libgco.so` in `ext/evocube/build/lib/gco`),
      patched AlgoHex as the default in `check_environment`/MCP, warn or refuse `size_max`/`remesh_all` when the
      input is already a tagged VTK tet mesh, beam post-check: report INVALID instead of the `structured_grids`
      ValueError, `report.json` for failed runs with timings. Details: `cases/hexdp_agent/HEXME_TEST_A.md` §4.
