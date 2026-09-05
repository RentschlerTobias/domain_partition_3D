# Proposal: step-by-step visual analysis of the 3D pipeline

**Proposal only — awaiting go.** Mirrors the 11-step showcase the 2D/surface
pipeline already has (`dp3d/plotting.py` → `output/plots/step01…step11`,
selection copied into `slides/figures/`).

Deliverables if approved:

1. `plot_stages.py` rewritten into an 11-step showcase with the same
   conventions as `dp3d/plotting.py` (fixed viewpoint, faint background mesh,
   plain + labeled variants, one PNG series per run tag).
2. A VTK/MSH export next to every step, so each figure can be re-inspected at
   full fidelity in ParaView instead of trusting a 2D render.
3. `RUNS.md` — the complete run table, below in draft form.
4. A metrics section that states every number with the command that produced
   it, on real data.

---

## 0. Two constraints worth agreeing on first

**Rendering.** Only `matplotlib` is available in the venv — no `pyvista`, no
`vtk`, no ParaView/`pvpython`. The existing `plot_stages.py` already works
around this with `Poly3DCollection` / `Line3DCollection`, which is fine for
surfaces, curves and wireframes but has no depth sorting, so a dense
volumetric render comes out muddy. Consequences for the plan below:

- steps that show *curves and wireframes* (feature graph, singularity graph,
  block edges) render well and are the backbone of the series;
- steps that show *volumes* are rendered as the boundary surface plus a clip,
  never as a full opaque block;
- every step also writes its VTK, which is the authoritative artifact.

If you want publication-quality 3D, `pip install pyvista` into the venv would
change this substantially (off-screen rendering, proper depth, scalar bars).
That is a change to a shared environment, so I am asking rather than doing.

**The existing `plot_stages.py` is stale.** Its docstring still advertises
`step04_hexmesh_v1.png  Stage 2 v1: extracted hex mesh (31.4% coverage)` —
that 31.4 % is the `_hex_volume` bug (README "What was learned" 5), not a
real coverage. It would be replaced, not extended.

---

## 1. The 11 steps

Each step lists: what it shows, the file it reads, the file it writes, and
the metric it demonstrates with a real measured value.

| # | figure | shows | reads | writes | metric shown |
|---|---|---|---|---|---|
| 01 | `step01_input_surfaces` | the labelled boundary of the AlgoHex input, one colour per surface | `data/T1_9/T1_9_tet_v5.vtk` | `01_input_surfaces.vtk` | 7 surfaces; hub 2178 / shroud 7912 triangles, matching the MSH's own wall triangulations exactly |
| 02 | `step02_feature_graph` | feature curves and feature vertices on that boundary | same | `02_feature_graph.vtk` | 534 feature edges, valence {2: 522, 3: 8}, 8 feature vertices |
| 03 | `step03_domain_variants` | the three domains side by side: full / minus O-grid / minus both | `T1_9_tet_v2/v6/v5.vtk` | — | tets 71 415 / 62 874 / 42 218; feature edges 630 / 854 / 534 |
| 04 | `step04_singular_graph` | **the frame-field singularity graph in 3D**, arcs coloured by valence (3 = red, 5 = blue), endpoints marked | `T1_9_final_tet_v11.ovm` | **`04_singular_graph.vtk`** (lines + `valence`) | v11: 164 singular edges, 6 arcs, 2 endpoints on the outlet, 0 on the inlet |
| 05 | `step05_hexmesh` | the extracted hex mesh: boundary surface + one clip | `T1_9_hex_v11.ovm` | `05_hexmesh.vtk`, `.msh` | 61 546 cells, volume 5.4688 against the domain's 5.4606 |
| 06 | `step06_quality` | scaled Jacobian on the boundary; inverted cells drawn opaque | same | `06_quality.vtk` (`scaled_jacobian`) | min −0.4938, mean 0.9767, 21 inverted (0.03 %) |
| 07 | `step07_sheets` | singular edges of the hex mesh and the separating sheets they emit | same | `07_sheets.vtk` (`sheet_id`) | 21 sheets over 14 202 faces, sizes 2136 … 27 |
| 08 | `step08_blocks_raw` | the raw base complex, one colour per block | same | `08_blocks_raw.vtk` (`block_id`) | 117 blocks, 114 cuboids (97 %), 100 % of cells in cuboids |
| 09 | `step09_postprocessing` | before/after of the cut-set cleanup, non-cuboids highlighted | — | `09_blocks_clean.vtk` (`n_block_faces`) | faces per block {6: 115, 7: 2}; what merging can and cannot do |
| 10 | `step10_sheet_collapse` | the three collapse rounds as a small multiple, 117 → 31 → 22 → 16 | — | `10_collapse_round{1,2,3}.vtk` | blocks 117→31→22→16, tiny 16→0, inverted 21→12→3→2 |
| 11 | `step11_final_blocks` | final structure: block-edge wireframe over a faint mesh, plus the re-attached full domain | `T1_9_blocks_v11.vtk`, `T1_9_blocks_v9_full.vtk` | `11_final_edges.vtk` | 16 blocks, all 6-faced, 0 inverted, VALID; full domain 73 blocks |

**Labeled variants** (`_labeled.png`) for steps 04, 08, 10, 11 — the ones
where an index per entity is what makes the figure readable, following
`showcase_integration_labeled` / `showcase_simplification_labeled` in
`dp3d/plotting.py`.

**New artifact worth calling out:** step 04 writes the singularity graph as a
VTK line mesh with a `valence` cell field. Nothing in the branch exports that
today — it has only ever been printed as a table — and it is the object the
whole `FRAMEFIELD_PLAN.md` stage is about.

---

## 2. Draft of `RUNS.md`

Every AlgoHex run on this branch, from the logs and metrics JSONs.

| run | input | domain | feature edges | runtime | outcome |
|---|---|---|---|---|---|
| v1 | `T1_9_tet.vtk` (classifySurfaces) | full | 853 | 2 h 10 | 95.4 % coverage |
| v2 | `T1_9_tet.vtk` (analytic) | full | 496 | — | SIGSEGV, dangling feature curves |
| v3 | `T1_9_tet.vtk` (analytic, cleaned) | full | 437 | — | IGM diverged (energy 3.3e48) |
| v4 | `T1_9_tet_v1tags.vtk` | full | 853 | 2 h 21 | 97.9 % coverage, 124 inverted |
| v5 | `T1_9_tet_v2.vtk` | full, original boundary | 630 | 1 h 17 | valid IGM, 96 inverted, 307 blocks |
| v6 | `T1_9_tet_v3.vtk` | minus blade O-grid | 646 | — | stopped manually at round 3 of 9 (memory) |
| v7 | `T1_9_tet_v4.vtk` | minus O-grid + BL | 166 | — | OOM in quantization |
| v8 | `T1_9_tet_v4.vtk` | same | 166 | — | OOM; `-n` is not the cause |
| v9 | `T1_9_tet_v5_diagbug.vtk` | reduced, **mislabelled** | 550 | 45 min | 2 inverted, 2 cavities, 82 blocks |
| **v10** | `T1_9_tet_v6.vtk` | only blade layer cut | 854 | 1 h 29 | 19 inverted, 198 blocks, 83 % cuboids |
| **v11** | `T1_9_tet_v5.vtk` | reduced, **labels fixed** | 534 | **13 min** | 21 inverted, 0 cavities, **117 blocks, 97 % cuboids** |

Post-processed results:

| basis | blocks | cuboids | tiny | inverted | Hausdorff | validator |
|---|---|---|---|---|---|---|
| v9 + cleanup | 81 | 73 (90 %) | 8 | 2 | 0.1851 | INVALID |
| v9 + collapse + untangle | 42 | 36 (86 %) | 2 | 0 | 0.0999 | VALID |
| **v11 + collapse + untangle** | **16** | 14 (88 %) | **0** | **0** | **0.0341** | **VALID** |
| v9 full domain (re-attached) | 73 | — | — | 5 | — | — |

---

## 3. The metrics, each with a real example

Not a glossary — every entry names the function, the command, and a measured
value, so the number in a figure can be traced.

| metric | where | what it caught, concretely |
|---|---|---|
| **scaled Jacobian** (min over 8 corners) | `ovm_io.scaled_jacobian`, vectorised `clean_blocks.scaled_jacobians` | v5 has 96 inverted cells against v9's 2 — the reason v9 became the base |
| **cell volume** | `ovm_io._hex_volume` | the corrected version showed a "31.4 % coverage" was really 95 %; the buggy variant is still live in `hexa_interpolation.py` |
| **coverage** = Σ hex volume / domain volume | census script | v11 5.4688 against the domain's 5.4606 |
| **singular edges / arcs** | `base_complex.singular_edges` | v5 308 edges / 12 arcs vs v11 164 / 6 |
| **arc endpoints per surface** | census script | all four runs: exactly 2 on the outlet — killed the "cutting causes it" hypothesis |
| **boundary components** | `clean_blocks.boundary_components` | 3 instead of 1 on v9 = two internal cavities, 28 and 22 quads |
| **faces per block** | `clean_blocks.cuboid_status` | v11 {6: 115, 7: 2} vs v9 {6: 72, 7: 4, 8: 4, 9: 3, 10: 1, 12: 1} |
| **excess faces** Σ max(0, faces−6) | `clean_blocks._excess` | replaced a non-cuboid *count* that had built 17- and 27-faced blocks while the percentage rose |
| **cuboid share vs count** | `report_structure` | 82→65 blocks looked like 78 %→91 % while the absolute cuboid count *fell* 64→59 |
| **Hausdorff to the input surface** | `HexBlockValidator.boundary_hausdorff` | v11 0.0341 vs v9 0.1851 — the label fix, not the postprocessing |
| **block-edge kink angle** | `clean_blocks.edge_kink_stats` | 63 kinks > 30°, of which 45 on 14 `shell_hub\|shell_blade` curves |
| **interface gap** | `reattach.interface_gap` | 0.0418 median — the non-conformity assumption checked, not assumed |

Each metric gets one figure panel showing the *distribution*, not just the
extremum — this branch has been burned five times by a mean or an extremum
moving the right way while the structure degraded (PROGRESS "Stage 5c").

---

## 4. What I would NOT do

- **No new renders of v9-based results as the headline.** v9 ran on the
  mislabelled input. Its figures belong in the analysis as the *before*, not
  as the result.
- **No re-run of v1–v8** to make the series complete. Their logs and metrics
  are on disk and go into `RUNS.md` as a table; regenerating them costs
  ~8 h and adds nothing.
- **No figures of the TFI stage** — it does not exist yet.

---

## 5. Open questions for you

1. **`pip install pyvista`** into the venv, or stay with matplotlib? This is
   the single biggest lever on figure quality.
2. **Which run is the showcase?** My recommendation: **v11** throughout, with
   v9 appearing only in steps 09/10 as the before-comparison, and v5/v10 only
   in step 03. The alternative — one full series per run — is 4x the figures
   for little gain.
3. **Slides too?** `slides/` has a working `presentation.tex` for the 2D
   pipeline. A parallel deck for the 3D pipeline is roughly a day on top; say
   if it should be in scope.
