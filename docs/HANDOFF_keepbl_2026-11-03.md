# HANDOFF — keep-BL / blade O-grid ring (2026-11-03)

For a session with **zero context**. Read top to bottom. Every number below was re-measured
from the files named in §5 on 2026-11-03, not copied from a chat. Where a number comes from a
commit message rather than a file, that is said.

Job `job-20261030-133532-3225`, integration branch **`jev/job-20261030-133532-3225`**
(`d333fc9` at the time of writing). Longer, older history:
`experimentell/hex3d_algohex/HANDOFF.md` and `.../STATE_2026-11-02.md`.

---

## 1. Goal and the user's hard requirements

We are rebuilding the **keep-BL** lane of `experimentell/hex3d_algohex`: the sheet collapse
keeps the O-grid and both boundary layers, and a **blade is inserted into the block structure
BEFORE TFI**. The reference run is `tistos_s1_0131`.

Non-negotiable, verbatim from the user:

| # | Requirement |
|---|---|
| R1 | Blade inserted **BEFORE** TFI (the TFI refill must see the blade ring, not be patched afterwards) |
| R2 | **No re-blocking of untouched blocks** — only the blocks whose faces the ring replaces may move |
| R3 | **No projection after TFI** — `tfi.project_ogrid_interface` / `--project-ogrid` stays OFF |
| R4 | **Keep 12 blocks** — do **not** split the r-band; the ring (blocks 12–15) must merge into the 12-block core |
| R5 | **Block edges must be curves, not chords**, and the curvature must be applied **before** TFI |
| R6 | **Do not re-run AlgoHex** — reuse `/work/datagen_keepbl/runs/tistos_s1_0131/hex_hex_x_16c.ovm` |
| R7 | **keep-BL ≠ prism identity.** Acceptance is "the AlgoHex input covers min r … max r" (radial coverage). Prism identity is informational only |
| R8 | Validation by **signed cell volume** (divergence theorem), not only sJ — a hex can have positive volume and inverted corners |
| R9 | **Re-run the existing ~1000 designs** (`/work/datagen/runs/tistos_s1_0000 … tistos_s1_0999`, exactly 1000 dirs, each with `params.json`). **No new Sobol samples** |

Blade wall-normal class = **class 10**, 4 axes (one per ring block), pinned to **10 layers**
(`--pin-class 10=10`).

---

## 2. State

Commits on `jev/job-20261030-133532-3225`, newest first (the tip was `d333fc9` at the
time of writing — check `git log` for anything newer):

| commit | task | what it fixed / produced |
|---|---|---|
| `6994ca5` | T88 | The handoff doc for a context-free session — this file |
| `47dd28d` | T90 | Decision record: keep 12 blocks and fix the misleading p1 check → `docs/decisions/2026-11-03-keep12-fix-p1-check.md` |
| `7c4a7b1` | T85 | Flip O-grid grading + rebuild ring + TFI. Inverted wall-normal grading in all three extrusion variants in `ogrid_extrude.py` (profile anchored at the blade wall). Rebuilt ring → `blocks_blade_graded.vtk` (16 blocks, 0 cells with signed volume ≤ 0, 105/105 shared interface vertices). Ran the pinned TFI → `tfi_blade_graded.vtk` + `blade_s2c_graded.json`. **Not promoted**, canonical files untouched |
| `528cc6e` | T83 | Curved block edges/faces before TFI: `--curved-edges {off,coons}` in `tfi.py`, boundary-face detection by surface labels (1,2,3,4,7,8,9), Coons-patch integration. Wrote `curved_edges_check.json` |
| `0312cf0` | T84 | Diagnose the TFI ring grading. Added `grading_check.py` + `blade/grading_check.json`. `tfi.py` untouched — it reproduces the input lattice and cannot repair an inverted input |
| `eb5d4ec` | T79 | Detached ring variants + TFI: OgridMap 257.9 s → 3.3 s (quadratic outer-quad scan → face-index lookup), `n` rebound bug in both variant extruders, const 1609 inverted / min sJ −0.9860, smooth 805 / −0.9953 |
| `8937f94` | T81 | Audit `tfi_blade.vtk` + make inversions visible: `audit_tfi_blade.py`, `tfi_blade_audit.{json,vtk,png}` |
| `418acf3` | T75 | Whole-curve singular merge test (canada) |
| `004dcf2` | T77 | Ring: remove extrusion twist (const vs smoothed) |
| `d43beef` | T78 | Fix misleading p1 boundary-layers check |
| `b771d00` | T70 | Blade S2c: pin wall-normal class + fix ring cell volumes |
| `f1e4df1` | T74 | Curve-level singular merge analysis + top-10 visual (`viz_singular.py`) |
| `dd3f68f` | T68 | Canada: finish singular-pair probe (vectorised) |
| `6d2bcfc` | T57 | Blade S2: TFI over passage + blade |
| `74649c1` | T66 | Canada: singular-pair merge with surface-preserving boundary |
| `0b02e8b` | T56 | Blade S1: extrude O-grid ring on P1 blocks |
| `1527e40` | T64 | Canada: local singular-pair merge probe |

**Nothing is promoted.** `blade/blocks_blade.vtk` is byte-identical to `blade/blocks_blade_v1.vtk`
(md5 `e90480dd22e5c413bf7385775db1b881`) — it is still the **original, wrongly graded** ring.
The fixed ring exists only as `blade/blocks_blade_graded.vtk`.

---

## 3. Key findings (with numbers)

### 3.1 Ring grading was inverted — real bug, but NOT the cause

`grading_check.py` verdict on the original ring: **INVERTED (YES)**.
`k=0` is the CORE / O-ring side (exactly on the 105 core↔ring shared interface vertices, distance
0.0); `k=9` is the blade wall. The thin layers landed on the **core** side.

| ring block | gap at blade side | gap at core side | blade/core ratio |
|---|---:|---:|---:|
| 12 | 0.0008484803816348348 | 0.01319324993140425 | **15.549** |
| 13 | 0.0008481359515794814 | 0.01316900406467573 | 15.549 |
| 14 | 0.0008483774554723526 | 0.01317996650650741 | 15.549 |
| 15 | 0.0008487258798610355 | 0.01318681826586806 | 15.549 |

All **20 worst cells sit at k=0** (`worst_cell_k_histogram = {"0": 20}`, `worst_cells_all_at_k0`).
Root cause: `ogrid_extrude.py` reversed the O-grid chain and applied the resulting profile
backwards, so the fine layers landed on the core. Three flip sites, now fixed by T85:
`extrude_ogrid` (line 520), `extrude_ogrid_const` (line 594), `extrude_ogrid_smooth` (line 684)
— the chain is `C = omap.chain(int(v))` and the reversal `[::-1]` is gone from all three.

**T85 fixed it** (`blade_s2c_graded.json`): thin side is now the blade, blade/core ratio
**0.0643** instead of 15.55.

**But the TFI got slightly worse**, so grading is *not* the cause of the ring inversions:

| metric | before (T70, `tfi_blade.vtk`) | after T85 (`tfi_blade_graded.vtk`) |
|---|---:|---:|
| cells sJ < 0 | 215 | **225** |
| cells sJ < 0.3 | 429 | **443** |
| cells signed volume ≤ 0 | 141 | **150** |
| min sJ | −0.9936709408917652 | **−0.9869441191815008** |

(before-numbers from `blade/tfi_blade_audit.json`, cross-checked and CONFIRMED there;
after-numbers from commit `7c4a7b1` — there is no audit JSON for the graded TFI on disk.)

**Conclusion: the cause is the ring's in-plane wall-normal twist** (T77/T79 measured const
1609 / smooth 805 inverted on the pre-fix variants; the twist, not the grading, is what folds
the cells). T85 therefore did not promote anything.

### 3.2 Curve merge — NO-GO

`/work/canada/sweep/pair_probe3/curve_analysis.json` (canada n=100000,
`/work/canada/sweep/n100000/hex_n100000.ovm`): 15 open singular curves (lengths 15.2177 … 0.42,
edge counts 111 … 3), 14 parallel curve pairs within tolerance. Of the **373** pair candidates,
**0** share a curve (`pairs_belonging_to_same_two_curves = 0`) — no two singular pairs can be
merged by collapsing one shared curve.

`curve_merge_test.json`: 3 whole-curve merge attempts, **`any_pass = false`**, `wall = boundary`.
The best (pair 1-5, tangent deviation 2.49°, min_dist 0.0513) failed with
`boundary surface not preserved: boundary vertex 61 deviates 0.891635 from original surface`
(`hausdorff_max 0.891635`, `boundary_vertex_count_change +322`, 94.6 s).
Verdict: *whole-curve merge not possible on canada; primary wall = boundary.*

### 3.3 Curved edges — conformity is the blocker

`curved_edges_check.json` + commit `528cc6e`. Baseline linear TFI: 77372 cells, **0 inverted**,
Hausdorff p95 **0.02212**.

| variant | inverted | Hausdorff p95 |
|---|---:|---:|
| linear baseline | 0 | 0.02212 |
| full Coons on boundary faces | 120 | — |
| Coons + linear edges | **0** | 0.0299 (no improvement) |
| surface-projected interior | **185** | **0.00196** (10× better) |

Best fidelity (p95 0.00196, hausdorff_mean 0.000517) comes with 185 inversions
(min sJ −0.9911). The 0-inverted option gives no fidelity gain. **Edge conformity is the
blocker**, not the interior: an edge-aware Coons that keeps the four edges of a boundary face
conformal, or post-refill untangling, is the open route.

### 3.4 P1 check fix (T78)

`/work/datagen_keepbl/runs/tistos_s1_0131/p1_report.json` (also checked from
`/work/datagen_keepbl/inspect/p1/`, which holds only the three VTKs) — `pass: true`,
`source_hex_mesh: hex_hex_x_16c.ovm`, 12 quality blocks.

* `bl_prism_identity` is **INFORMATIONAL ONLY** (`informational: true`, `gate: false`):
  hub coverage 0.9987, shroud coverage **0.0336** (715 of 21258 prisms inside,
  20543 outside). keep-BL does not require prism identity — TFI rebuilds the 0.027–0.035
  outer layer, so a low shroud value is expected geometry, not a lost layer.
* **`2_boundary_layers_kept`** is the honest gated keep-BL check: from block vertices,
  min r 0.489804 ≤ r_hub+3e-3, max r 1.900002 ≥ r_shroud−3e-3, and both BL bands
  [0.5, 0.5265] (hub) and [1.872, 1.9] (shroud) inside [min r, max r] → **pass**.
* `3_ogrid_outline_coverage`: min seam coverage 0.99988 ≥ 0.99. `4_block_count`: 12 blocks
  in [10, 16] → pass.
* Stage walls in the same file: dtoo 500 s, tet_prep 27 s, beam 918 s, tfi 4.6 s, export 1.2 s,
  curves 1.35 s — but **algohex exit 124 (timeout at 10800 s)**, recovered by the 16-core run.

### 3.5 2D revolve background blocks — NO-GO

T72 (tistos) and T80 (canada) studied revolved 2D r-z background blocks. Both are no-go.
**No commit, no artifact** — `git log --all --grep "T72\|T80"` is empty; the finding exists only
in the task titles. Do not spend time looking for a report.

### 3.6 The P1 reference that works (for contrast)

`/work/datagen_keepbl/runs/tistos_s1_0131` built on `hex_hex_x_16c.ovm` (16 cores, same 2654
points as `hex_hex_x.ovm` but different cell grouping): beam `--guard struct` → **12 blocks /
12 cuboids, 582 cells, 0 inverted, min sJ 0.3814, boundary Hausdorff 0.0590, VALID**; TFI
582 → **77372 cells at h=0.05, 0 inverted, min sJ 0.3814, mean 0.9533**; `sample.npz` 40
vertices / 12 blocks / 168 directed edges; `curve_postprocess` 38 feature arcs, G1 fit dev/L
max 0.066. Exports for ParaView: `/work/datagen_keepbl/inspect/p1/{blocks,tfi,overview}.vtk`.

---

## 4. Open tasks and the exact next step

| task | status | exact next step |
|---|---|---|
| **T69** | blocked | Merge the ring 16 → 12 blocks (blocks 12–15 into the core, no r-band split, no re-blocking of untouched blocks) **and then** TFI. Blocked on a valid ring: the ring is still inverted. Do not attempt the merge on `blocks_blade.vtk`; wait for a ring with 0 inverted cells after TFI |
| **T58** | open | Blade S3: `blade_check.py s3` currently prints *"not implemented yet (no TFI on the ring)"*. Implement it: surface labels + `sample.npz` + curves. The 12-block core already exports (`sample.npz` / `sample_16c.npz`, 40 vertices, 168 edges, `surface_points` 17938 / `surface_tris` 26728 / `surface_tri_label`) and `curve_postprocess.py` produces `sample_16c_curves.npz` (post_* arrays, 24 pairs). The ring variant needs the same from `blocks_blade_graded.vtk` |
| **T59** | open | S4 lane wiring + HANDOFF section. `datagen/run_sample_keepbl.sh` is the keep-BL lane (dtOO → tet_prep → AlgoHex → `beam_collapse --guard struct` on `hex_hex_x_16c.ovm` → TFI → `export_sample --check` → `curve_postprocess`). It is **not** wired into `datagen/run_lane.sh` (no `RUN_SAMPLE` hook, `grep -c RUN_SAMPLE datagen/run_lane.sh` = 0). Add the hook, then add a keep-BL section to `experimentell/hex3d_algohex/HANDOFF.md` |
| **T7** | re-scoped | Re-run the existing **1000** keep-BL designs in `/work/datagen/runs/tistos_s1_0000 … tistos_s1_0999` (each has `params.json`; `tistos_s1_0000` has `tet.vtk` + `params.json` + logs, `tistos_s1_0999` only `params.json` — the old lane is partly finished). **No new Sobol samples.** `run_sample.sh` skips a dir that already has `sample.npz`, so re-running overwrites/extends only what is missing |
| **T6** | open | Lane harness: `datagen/run_lane.sh` / `start_lanes.sh` / `run_sample.sh`. Add the RUN_SAMPLE hook that picks `run_sample_keepbl.sh` vs `run_sample.sh` |
| **T8** | open | Acceptance review of the whole job (A1–A89) against §1's requirements |
| **T9** | open | Docs consolidation: fold `STATE_2026-11-02.md` and this file into `experimentell/hex3d_algohex/HANDOFF.md` |

### The two viz tasks of 2026-11-03 (T87 merged, T86 in flight)

| task | question | where the output lands | tool |
|---|---|---|---|
| **T86** | singular merge curves: single-pair + whole-curve | `/work/canada/sweep/pair_probe3/viz/merge_curves.vtp` | `experimentell/hex3d_algohex/viz_merge_curves.py` (does not exist yet; branch `micromanager/job-20261030-133532-3225--T86`) |
| **T87** | tistos inverted cells + surface + linear/curved blocks | `/work/datagen_keepbl/runs/tistos_s1_0131/blade/viz/` — written 2026-11-03 10:52: `tistos_inverted.vtm` + `tistos_inverted.png` + `tistos_inverted.json` (215 sJ<0 cells, 429 below 0.3, 141 volume<=0, min sJ -0.99367 at cell 87410 block 12) | `experimentell/hex3d_algohex/viz_inverted_cells.py` (in the repo since `c38984e`; branch `micromanager/job-20261030-133532-3225--T87`) |

T87 has landed as `c38984e`; T86 is still in flight (its `git log --oneline <branch>
--not jev/job-…` is empty). The predecessor `viz_singular.py` (T74) already writes
`/work/canada/sweep/pair_probe3/viz/{singular_curves_overview,top10_pairs}.png` — use it as
the style reference. **T89** verifies the p1-vs-r-band decision.

---

## 5. File map

`R = /work/datagen_keepbl/runs/tistos_s1_0131` (the reference run; **not** a git repo)

### Run root `R`

| path | what |
|---|---|
| `hex_hex_x_16c.ovm` | **the AlgoHex output to reuse** (R6). Do not re-run AlgoHex |
| `tet.vtk` | keep-BL tet input, feature `color` (580 feature edges, 3 curves 462/72/46, surface tags 1,2,3,4,7,8,9) |
| `case/mesh.msh` | dtOO mesh; hub r 0.5, shroud r 1.9; prism source for the p1 identity check |
| `blocks.vtk` | the **12-block** core (582 cells, min sJ 0.3814) — the input to the ring build |
| `blocks.divisions.json` | core divisions |
| `sample.npz`, `sample_16c.npz` | the exported sample (576539 B and 576542 B, **different md5** — `sample_16c.npz` is the 16-core export): `vertices` 40, `blocks` 12×8, `quad_faces` 40, `edges` 168, `edge_ctrl` 168×2×3, `edge_polyline` 812, `dir_class` 168 + `dir_class_count` 10, `params`, `quality`, `provenance`, `surface_points` 17938 / `surface_tris` 26728 / `surface_tri_label` |
| `tfi.vtk`, `tfi_16c.vtk`, `tfi_curved_{c,d,e,f,g,h,i,j,k,l,m,n,p}.vtk` | TFI refills of earlier experiments; the `tfi_curved_*` set (13 files) is the T83 curved-edge sweep |
| `sample_16c_curves.npz` | `curve_postprocess` output: `post_*` arrays (vertices 40, edge_ctrl 168×2×3, feature_edge 168, arc_polyline 4940, vertex_role, pairs 24×4, quality) |
| `tistos_s1_0131_curves.npz` | another curve export |
| `beam.log`, `beam16c.log`, `beam_driver.out` | beam collapse logs |
| `curves.log`, `curves16c.log`, `export.log`, `export16c.log` | stage logs |
| `hexmeshing_x.log`, `algohex*.log`, `algohex_*_*.out`, `dtoo.log`, `tet_prep.log` | stage logs |
| `hex_hex_metrics_{x,x_16c,1000,1500,3000}.json` | per-run quality metrics |
| `pipeline.out`, `prev_x_ovm`, `replay16.log`, `tail_driver.out`, `param_state.ovm`, `algohex_replay{,_1000,_1500,_3000}.out`, `algohex_save.out` | driver leftovers |
| `curved_edges_check.json` | the T83 comparison of linear vs curved TFI (baseline p95 0.02212, 77372 cells, 0 inverted) |
| `blocks_{edges,nfaces,quality}.vtk`, `blocks.msh`, `blocks.divisions.json` | the 12-block core's files |
| `blocks_16c.{vtk,msh,divisions.json}`, `blocks_16c_{edges,nfaces,quality}.vtk` | the 16-core collapse variant |
| `blocks_blade_{const,smooth}.{vtk,divisions.json,part.vtk}`, `blocks_blade_graded.{divisions.json,part.vtk}` | the ring variants' companion files |

### `R/blade` (the blade O-grid ring work)

| path | what |
|---|---|
| `blocks_blade.vtk` | **canonical, still the original wrongly graded ring** (md5 `e90480dd…` = `blocks_blade_v1.vtk`) |
| `blocks_blade.vtk.before_swap` | backup, identical to the above |
| `blocks_blade_graded.vtk` | **the T85 fixed ring** (16 blocks, 1482 cells, 0 inverted, 105/105 shared vertices, min sJ 0.18978, mean 0.86714) |
| `blocks_blade_graded_test.vtk` | copy of the original (md5 matches) |
| `blocks_blade_v1.vtk` | the original ring |
| `blocks_blade_part.vtk` | 0 = AlgoHex core, 1 = extruded O-grid |
| `blocks_blade.divisions.json` | ring divisions |
| `blade_s1.json`, `blade_s1_{const,smooth}.json`, `blade_s1_graded.json`, `blade_s1_graded.json.bak` | S1 ring-build reports (16 blocks, 582 core + 900 ring cells, ring ids 12–15, 4 sides, quads/side 30/30/18/12, replaced core blocks 4–7, 225 layers) |
| `tfi_blade.vtk` | TFI of the **original** ring — the current 215-inverted reference. It carries **only** cell data `block_id` (no sJ, no volume), so any quality/inversion work must use `tfi_blade_audit.vtk` |
| `tfi_blade_audit.json` | its audit: 104940 cells, 16 blocks, 215 sJ<0, 429 sJ<0.3, 141 volume≤0, min sJ −0.99367, cross-check CONFIRMED |
| `tfi_blade_audit.vtk`, `tfi_blade_audit.png` | the same as viewable output |
| `tfi_blade_const.vtk`, `tfi_blade_smooth.vtk` | twist variants (const 1609 inverted, smooth 805) |
| `tfi_blade_graded.vtk` | TFI of the T85 ring: 225 sJ<0, 443 sJ<0.3, 150 volume≤0, min sJ −0.98694 |
| `blade_s2c.json`, `blade_s2c_graded.json`, `blade_s2c_{const,smooth}.json`, `blade_variants_s2c.json`, `blade_s2.json` | S2c reports; the graded one carries the old/new grading verdicts and `promote: false` |
| `grading_check.json` | the T84 report — **see the trap in §6**: it was produced while the graded ring was swapped in |
| `ogrid_const.log`, `ogrid_smooth.log`, `tfi_const.log`, `tfi_smooth.log` | ring-build and TFI logs for the const/smooth variants |
| `viz/` | T87's outputs (written 2026-11-03 10:52): `tistos_inverted.vtm`, `tistos_inverted.png`, `tistos_inverted.json` (215 sJ<0, 429 sJ<0.3, 141 volume<=0, min sJ -0.99367 at cell 87410 block 12) |

### Elsewhere

| path | what |
|---|---|
| `/work/datagen_keepbl/inspect/p1/{blocks,tfi,overview}.vtk` | the working P1 12-block export for ParaView (three files, no JSON) |
| `/work/datagen_keepbl/runs/tistos_s1_0131/p1_report.json` | the P1 acceptance report (pass, 12 blocks, source `hex_hex_x_16c.ovm`) |
| `/work/datagen_keepbl/inspect/keepbl_features.vtk` | feature curves of the keep-BL tet input |
| `/work/datagen/runs/tistos_s1_0000 … 0999` | the **1000 existing designs** to re-run (R9) |
| `/work/canada/sweep/pair_probe3/` | `curve_analysis.json`, `curve_merge_test.json`, `pair_probe3.json`, `detect_top10.json`, `viz/` (T74 PNGs, T86 target) |
| `/work/canada/sweep/struct_beam/` | the struct-guard beam report (T61/T40) — the struct-guard mesh is `blocks.vtk` in that directory |
| `/work/canada/sweep/n100000/hex_n100000.ovm` | the canada mesh the curve-merge probes ran on |
| `/work/canada/sweep/inspect/` | canada VTK exports |
| `/work/calib/` | hexme + Charité LV calibration |

### Repo files that matter

| path | what |
|---|---|
| `experimentell/hex3d_algohex/ogrid_extrude.py` | the ring extruder (the grading flip site) |
| `experimentell/hex3d_algohex/grading_check.py` | the grading orientation checker |
| `experimentell/hex3d_algohex/tfi.py` | the TFI refill (`--pin-class`, `--curved-edges`, `--project-ogrid`) |
| `experimentell/hex3d_algohex/blade_check.py` | blade stages s1 / s2 / s2c (s3 not implemented) |
| `experimentell/hex3d_algohex/audit_tfi_blade.py` | the blade TFI audit (sJ + signed volume) |
| `experimentell/hex3d_algohex/curve_merge.py` | singular curve/pair merge probe |
| `experimentell/hex3d_algohex/datagen/run_sample_keepbl.sh` | the keep-BL lane |
| `experimentell/hex3d_algohex/datagen/README.md` | the keep-BL lane reference incl. the five P1 checks |
| `experimentell/hex3d_algohex/HANDOFF.md` | the long handoff (longer history) |
| `experimentell/hex3d_algohex/STATE_2026-11-02.md` | the state digest |
| `docs/decisions/2026-11-03-keep12-fix-p1-check.md` | the T90 decision record: keep 12 blocks, fix the misleading p1 check (prism identity informational, `2_boundary_layers_kept` is the gate) |

---

## 6. Traps that cost hours

| trap | what to do |
|---|---|
| **Foreground commands die at ~120 s** (the shell tool's default timeout is 120000 ms) | detach: `subprocess.Popen(cmd, stdout=log, stderr=log, start_new_session=True)` and poll the log. Never `while pgrep -f <pattern>` — the pattern matches the waiting shell's own command line and loops forever (it cost two silent 20-minute stalls) |
| **The shell mangles `&` and `>` into pipes** | write the command into a `.sh` file and run that file, or use `subprocess.Popen(..., start_new_session=True)` |
| **Piping to `tail` hides stderr** | always append `2>&1` |
| **`finish_task` needs an uncommitted change** and needs the task's check file to exist | create/modify the file first, then call `finish_task` |
| **Worker test symlinks must not be committed** | `output -> /opt/stack/meshtron/data/hex3d_algohex_compressed` is gitignored; check `git status --short` before finishing |
| **`/work` is not a git repo** | everything under `/work` is scratch; the repo's only records of it are commit messages and the JSONs |
| **The main checkout is a linked worktree** | `/opt/stack/dp3d_keep_bl` is on `feat/keep-bl` (`f4075dd`) and is a worktree of `/opt/stack/domain_partition_3D` (which is on `feat/canada`). `git worktree list` shows all task worktrees |
| **`grading_check.py` hardcodes `blade/blocks_blade.vtk`** | it does **not** take a blocks path. The on-disk `grading_check.json` (verdict NOT INVERTED, ratio 0.0643) was produced while the T85 graded ring was swapped into `blocks_blade.vtk`; the file was restored afterwards. Re-running it **today** on the restored original returns **INVERTED (YES), ratio 15.549** (verified 2026-11-03). To audit the graded ring, copy it in (and restore), or point `--run-dir` at a copy of the run |
| **`grading_check.py` reads `blade/tfi_blade_audit.json` for its worst cells** | so a re-run mixes the *current* ring verdict with the *pre-flip* TFI audit's worst cells; the k=0 histogram is about `tfi_blade.vtk`, not about the graded ring |
| **`tfi.py` cannot fix an inverted input** | it reproduces the input lattice's profile. Fix `ogrid_extrude.py`, never `tfi.py` |
| **Class 10 auto-pin** | `tfi.py` pins the wall-normal axis of blocks ≥ 12 to the lattice count; `--pin-class 10=10` is the explicit form (log line: `class 10 (4 axes): pinned 10, free 1 <-- +900 %`) |
| **AlgoHex ignores OMP/OpenBLAS thread vars** | give every sample its own cores (one lane per core pair); five samples on six cores gave load ~230, AlgoHex 41–91 min |
| **dtOO and AlgoHex can hang forever** | every stage runs under `timeout` (TO_DTOO 1800, TO_ALGOHEX 10800, TO_BEAM 3600 — 7200 in the keep-BL lane, TO_TET_PREP 600 — 900, others 600) |

---

## 7. How to resume

```bash
# 1. Where are we?
cd <your worktree>
git log --oneline -15                      # tip was d333fc9 at the time of writing (jev/job-20261030-133532-3225)
git status --short                          # must be empty before you start
ls /work/datagen_keepbl/runs/tistos_s1_0131/blade/
```
Then ask the micromanager agent for `list_tasks` / `job_status` (a subagent may not call them).

```bash
# 2. Rebuild the ring (fixes the grading; ~minutes; writes 16 blocks, 0 inverted)
cd experimentell/hex3d_algohex
/opt/venv-sci/bin/python ogrid_extrude.py /work/datagen_keepbl/runs/tistos_s1_0131/blocks.vtk \
  --msh /work/datagen_keepbl/runs/tistos_s1_0131/case/mesh.msh \
  --input-vtk /work/datagen_keepbl/runs/tistos_s1_0131/tet.vtk \
  --out /work/datagen_keepbl/runs/tistos_s1_0131/blade/blocks_blade_graded.vtk \
  --json /work/datagen_keepbl/runs/tistos_s1_0131/blade/blade_s1_graded.json
```
(`--variant` is `const` / `smooth` / `original`, default `original`; the graded ring above is the
default variant with the T85 profile fix.)

```bash
# 3. Check the grading orientation (writes blade/grading_check.json; reads the canonical
#    blocks_blade.vtk -- swap the graded ring in first if that is what you want to measure)
/opt/venv-sci/bin/python grading_check.py --run-dir /work/datagen_keepbl/runs/tistos_s1_0131
#    ... or a scratch copy:  --run-dir <copy> --out <copy>/blade/gc.json
```

```bash
# 4. TFI, pinned (class 10 = 10 layers, h = 0.05, no projection)
/opt/venv-sci/bin/python tfi.py /work/datagen_keepbl/runs/tistos_s1_0131/blade/blocks_blade_graded.vtk \
  --input-vtk /work/datagen_keepbl/runs/tistos_s1_0131/tet.vtk \
  --target-h 0.05 --pin-class 10=10 \
  --out /work/datagen_keepbl/runs/tistos_s1_0131/blade/tfi_blade_graded.vtk
```
Expected log shape (from `blade/tfi_const.log`): `class 10 (4 axes): pinned 10`, `target h=0.05:
11 classes, counts [6, 8, 8, 10, 13, 16, 20, 24, 27, 27, 30], 1 pinned (optimal)`,
`refill: 1482 -> 104940 cells`. Do **not** add `--project-ogrid` (R3).

```bash
# 5. Audit (signed volume by divergence theorem + min corner sJ; writes json + vtk + png)
/opt/venv-sci/bin/python audit_tfi_blade.py \
  --vtk /work/datagen_keepbl/runs/tistos_s1_0131/blade/tfi_blade_graded.vtk \
  --out-json /work/datagen_keepbl/runs/tistos_s1_0131/blade/tfi_blade_graded_audit.json \
  --out-vtk  /work/datagen_keepbl/runs/tistos_s1_0131/blade/tfi_blade_graded_audit.vtk \
  --out-png  /work/datagen_keepbl/runs/tistos_s1_0131/blade/tfi_blade_graded_audit.png
```
Acceptance for the ring: **0 inverted** (sJ ≤ 0) **and** 0 signed-volume ≤ 0 cells. Today the
graded ring reaches 225 / 150 — it does **not** pass, which is why T69 is blocked.

```bash
# 6. The keep-BL lane, one sample end to end (see datagen/README.md; > 30 min, detach it)
cd experimentell/hex3d_algohex/datagen
./run_sample_keepbl.sh /work/datagen/runs/tistos_s1_0000          # writes p1_report.json
/opt/venv-sci/bin/python check_keepbl_sample.py /work/datagen/runs/tistos_s1_0000
```
`run_sample_keepbl.sh <run_dir>` skips a run that already has `sample.npz`, so on a re-run
delete that file first. `run_lane.sh` (the multi-lane driver) takes `DATAGEN_WORK` and a core
list and still has **no `RUN_SAMPLE` hook** — that is T6/T59.
