# hexdp task list with progress

Status: `[ ]` open, `[~]` in progress / partly done, `[x]` done, `[?]` waiting for the owner's decision.
Update the status, the date and a one-line note when something changes. Context and rules: `docs/HANDOFF_2026-11-19.md`.
Never start long computations without the owner's go (see the rules there).

## A. Running / to evaluate (state 2026-11-19 evening)
- [~] A1 HexMe agent test A (10 cases, solid + fluid, free method choice). Evaluate per case: strategies tried,
      blocks, valid, min sJ, time, failure reason; what methods were missing; MCP usability issues.
      Job dirs `/opt/stack/cases/hexdp_agent/runs/`. Write the result to `cases/hexdp_agent/HEXME_TEST_A.md`.
- [~] A2 canada runner p_0018 with stepwise snap + search (`/work/scratch/canada_snap/beam2.log`): did the snap
      keep (Hausdorff, inverted), how many blocks, valid? If good: rerun the other rejected ru samples 0016, 0017.
- [~] A3 radMeridional variation test (`/work/scratch/radmer/vartest/results.txt`): finish, tabulate success vs
      variation width and the failing parameters (factors in `*.factors.json`).
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
- [ ] B9 snap: if the stepwise snap still adds inversions, try smoothing (Winslow / Mesquite) after projection
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
