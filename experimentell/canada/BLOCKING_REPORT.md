# Automatic hex blocking for the canada machine: diagnosis and method study

Date: 2026-11-12/13. Scratch outputs: `/work/canada/blocking_study/`.
Tags used below: **[V]** = verified here by measurement or experiment, **[R]** = reconstructed from logs and
earlier transcripts (numbers taken from the files named), **[J]** = literature or engineering judgement, not
tested.

---

## 1. What was tried and how it failed [R]

### 1.1 Runner (ru) and guide vane (gv) passages
The source for these is the dtOO passage templates (same entity layout as tistos). The blade O-grid is cut out
(label 7) and AlgoHex is run on the rest. Then beam/sheet collapse, then TFI.

| Variant | Dir | Result |
|---|---|---|
| ru, dtOO rectangular O-grid cut (85° TE corners) + corner rule | `s0_corner` | frame field OK, parametrisation 777 invalid tets, Hausdorff 0.35, rejected |
| same + `--remesh-ogrid 0.03` | `s0_corner_rm` | AlgoHex exit 245 after 57 s |
| ru, round mmg envelope δ=0.083 (reduced mode) | `envelope/e2` | AlgoHex OK, beam 12 blocks / 11 cuboids, **TFI 8 inverted** (hub LE) |
| ru fluid mode with TE ridge | `envelope/ru_geo` | 20 cuboids, TFI min sJ 0.11 |
| gv reduced, δ=0.078 | `envelope/gv1` | LE–inlet gap 0.021 → **collapse to 5 blocks (artefact)** |
| gv fluid δ=0.05, size-max 0.08 | `gv10_d050` | 12 blocks / 11 cuboids |
| gv fluid + geometry TE ridge | `gv12_ridges` | 12 cuboids, VALID, but **4 blocks of 4 cells** (the sliver between envelope and inlet), TFI min sJ 0.045 |
| gv size-max 0.12 | – | 86–116 singular edges or broken parametrisation |

### 1.2 Draft tube / blk domain (Saugrohr with lip + 2 diffuser rings, quarter model)
| Variant | Dir | Result |
|---|---|---|
| Blade-free domain, rings as real walls, n=2000 / 20000 | `dt/b2`, `dt/b2_n20k` | 129 / 159 invalid tets |
| Rings as walls, clean geometry curves | `dt/geo4` | n=20000, 10.9 h, **3016 flipped tets**, beam untangle stuck (67 inverted) |
| Mesh-level envelopes + ridge cuts | `dt/env1..6`, `envg1` | preprocessing failed: fragmented ridges, plane drift |
| CAD envelopes (OCC, δ=0.035, knife edges blunted) | `dt/cad1` | n=2000/8000/20000: 18137/8985/7156 flipped after rounding, **~510 flipped already in the seamless map (independent of n)**, holes in the hex mesh up to 1.2 deep, beam rejects all (Hausdorff 1.2–1.4 > 0.15) |
| cad1 resume, n=100000 | `dt/cad1/state`, `sweep/n100000` | slots filled (≈1 cell across), 2551 flipped, 148 inverted hexes; clean → 330 blocks, 10 inverted cells that cannot be untangled |
| beam collapse on n=100000 (`--guard struct`) | `dt/cad1/state/struct` | 14 h → **109 blocks, INVALID** (6 inverted, 14 non-cuboid); macro search 159; macro+cell 159 |
| n=300000 | – | worse (5606 flipped) |
| AlgoHex size-field patch (`/work/tools/algohex-size`) | – | no effect: the number of layers across a slot can only change at singularities, and the singularities are fixed before the sizing is applied |

A hand estimate of the 2D meridional layout gave **~15–20 blocks** (transcript 2026-10-27). The automatic
results are 5–10× above that.

### 1.3 "2D blocking + sweep over span"
The user reports that this failed with an older model. No artefacts of that attempt exist in `/work/canada`
or in the transcripts searched (`e6473e12…`, `d6303bce…`, `3ba50fb8…`). Section 3.4 below discusses likely
causes [J].

---

## 2. Geometry measurements [V]

Script: `blocking_study/gaps.py` (minimum distance between surface labels of the AlgoHex input, volume, and
the uniform AlgoHex hex edge h = (V/n)^(1/3)).

| Domain | V | h @ n=2000 | Narrowest gap after cut-out | gap / h |
|---|---|---|---|---|
| gv, fluid δ=0.05 (`gv12_ridges`) | 1.61 | 0.093 | envelope–inlet **0.049**; envelope–outlet 0.118 | 0.53 / 1.27 |
| gv, reduced δ=0.078 (`gv1`) | 1.58 | 0.093 | envelope–inlet **0.021** | 0.23 |
| ru, δ=0.083 (`e2`, the one that works) | 1.85 | 0.098 | outlet 0.135, inlet 0.149 | 1.4–1.5 |
| draft tube CAD envelopes (`dt/cad1`) | 254 | 0.50 (0.23 @ 20k, 0.137 @ 100k) | lip–ring2 **0.061**, ring2–ring3 **0.099** | 0.12 (0.26 @ 20k, 0.45 @ 100k) |

- Raw gv LE–inlet distance ≈ 0.099. The envelope offset δ eats 50 % (δ=0.05) or 79 % (δ=0.078) of it.
- Diffuser rings: mean thickness 2V/A = **0.069** (ring2) and **0.065** (ring3), from `cad/ring*_solid.stl`.
  The envelopes are 0.12 thick. Raw slot widths are 0.061+2·0.035 ≈ **0.13** and **0.17**.
  The distances between ring mid-surfaces are ≈ **0.20** and **0.24**.
- Scale ratio in the draft tube: domain ~14 long, slot 0.06 → ratio ≈ 230. To get one uniform cell into the
  slot you need V/g³ ≈ 1.1 M hexes.
- **The draft-tube quarter is a rotational sweep.** `blocking_study/sweep_sections.py` cuts `cad/domain_blk.brep`
  with the half planes φ = 0.5, 10, 20, 30, 45, 60, 70, 80, 89.5°. Every section is **one face with exactly 2
  holes** (ring2 and ring3 envelopes; the lip hangs on the outer boundary) and the same sequence of boundary
  labels. Only the CAD sub-edge count on the ring3 envelope varies (4 vs. 6 pieces), which is a
  face-splitting artefact of the offset, not a topology change. The section area shrinks from 88 to 35
  (elliptical far field), but the topology stays the same.

---

## 3. Root causes

### 3.1 Thin volume in front of and behind the blade (gv, ru) [V for numbers, J for interpretation]
The gap is short because three things add up:
1. **Domain truncation.** The dtOO passage interface sits only 0.099 (gv) / ~0.23 (ru) in front of the LE. In
   the assembled machine this plane is an artificial interface, not a geometric feature.
2. **The envelope offset subtracts δ from that gap.** With δ=0.05–0.078 only 0.02–0.05 remains, which is
   0.23–0.53 of the uniform hex size.
3. **AlgoHex has one global size.** It scales the frame field by cbrt(V/n) (`Parametrization3DT_impl.hh`,
   `scale_frames_to_target_complexity`). A gap below ~h/2 rounds to 0 or 1 layers. With 0 layers the
   envelope "touches" the inlet: this is the 5-block artefact of gv1. With 1 layer you get the 4-cell sliver
   blocks of gv12, which TFI then fills with min sJ 0.045.

The ru passage only works because its gap is 1.4 h.

So the thin volume is a consequence of **cutting a volume out** (offset δ > 0) on a domain that is
**artificially truncated**. It is not a property of the blade. If the blade stays a zero-thickness surface,
the full 0.099 is kept. If the interface plane moves upstream, or the gv+ru+dt are blocked jointly, the gap
disappears altogether.

### 3.2 Draft tube: inverted cells and 109–330 blocks [V/R]
1. **The thin rings are not the frame-field problem; the scale is.** In the fluid domain a ring knife edge
   (solid angle 5–8°) is a re-entrant edge of 352–355°. AlgoHex computes the sector valence as
   `round(angle / 90°)` (`SingularGraphExtractionT_impl.hh:745`), which gives 4 sectors. That is regular, like
   a crack tip. The "0 sectors" problem only arises where the *domain* is a sharp wedge, for example the
   runner solid TE. What actually breaks is ring thickness 0.065 and slots 0.06–0.13 against a uniform
   h = 0.14–0.5. The integer grid map must place several integer iso-surfaces within less than one cell, so
   rounding flips tets. A fixed rest of ~510 seamless flips comes from the frame field and its singularities,
   which are fixed before n acts. This is why n=300k gets worse again.
2. **The envelopes make the slot narrower** (0.13 → 0.061). They trade the knife edge, which was harmless
   for the fluid, for a smaller gap, which is fatal.
3. **The block count follows from the raw singularity graph.** The n=100k hex mesh has 634 singular edges,
   plus holes and inverted cells. Sheet collapse with depth 10 can only remove a part of them; the
   109-block run still improved in its last depth. The domain's true structure is that of a 2D section ×
   sweep (section 2), so the minimal block count is roughly the 2D block count of the meridional section
   plus a few blocks for the axis treatment: ~15–25.
4. **Re-meshing the rings afterwards (the user's doubt) is justified.** Every envelope is a closed hole
   whose O-layer has to be rebuilt and re-inserted at block level. For three rings plus the lip, and with
   slots narrower than the envelope offset, there is no room left for the O-layers.

### 3.3 Why AlgoHex is structurally a poor fit here [J]
AlgoHex optimises a smooth frame field globally and then quantises with one global size. Thin slots,
multi-scale geometry and near-sweep domains are exactly the cases where this is weakest: rounding errors,
"holes" from HexEx, and many small singularities that collapse cannot remove. It is excellent for compact
CAD parts (HexMe). It is not designed for long, thin, sweep-like flow domains.

### 3.4 Why "2D blocking + sweep" may have failed before [J]
Likely causes; none can be checked, because no artefacts were found:
- **Wrong sweep direction for the blade rows.** The span direction hub → shroud is right. The φ direction
  is not, because the blade is in the way.
- **A layout per section instead of one layout swept.** Separate 2D layouts at several spans are not
  topologically identical, so the blocks do not connect.
- **A raw 2D quad mesh used as the layout.** It has dozens of irregular vertices. In this study a gmsh
  quasi-structured quad mesh (cross field, Algorithm 11) of the draft-tube section has 66 irregular
  vertices, and naive separatrix tracing gives 751–2087 patches (`blocking_study/sweep2d`) **[V]**. The 2D
  layout needs the same simplification step as 3D: irregular-vertex pairing / chord collapse, or a
  motorcycle graph with minimal separatrices. The difference is that in 2D this step is cheap and well
  understood.
- **The axis.** Blocks touching the rotation axis degenerate (a face collapses to an edge) unless a
  quarter-butterfly core is used.
- **Interfaces between the components** (gv↔ru↔dt) were not made conforming.

---

## 4. Experiments run here

All runs used the patched AlgoHex build `/work/datagen_keepbl/algohex_fix/HexMeshing`, pinned with
`taskset -c 96-103`, `OMP_NUM_THREADS=2`, each run limited by `timeout 3600`.

### 4.1 Synthetic "diffuser slab": walls vs. envelope vs. zero-thickness baffle [V]
The test geometry is a box 4 × 1 × 2 (x axial, z radial, y = sweep direction, standing in for φ between the
symmetry planes). It contains two staggered plates that span the full width in y, like the rings:
- thickness t = 0.065 (same as the canada rings), blunt LE, 8° knife TE
- raw slot g = 0.13, inclination 8° / 12°
- uniform h = 0.16 at n=2000, so t/h ≈ 0.4 and g/h ≈ 0.8. The ratios are close to canada at n≈20k.

Generator: `blocking_study/slab/gen_slab.py <walls|env|baffle>`. It writes the AlgoHex VTK. Feature curves are
label boundaries, creases > 25° and free edges of internal surfaces; tangent seams are dropped.

- **walls:** plates cut out as holes.
- **env:** plates replaced by their δ=0.035 Minkowski envelope, cut out.
- **baffle:** the plates' zero-thickness mid surfaces are embedded in the tet mesh (gmsh fragment) and tagged
  as **interior feature faces** (`face_colors` on interior triangles). Their free LE/TE lines are feature
  edges.
- **crack:** the same baffle, but split into two coincident boundary sheets before AlgoHex (`crack.py`).

| Variant | Tets | Seamless map | Quantisation | Result |
|---|---|---|---|---|
| walls, n=2000 | 87k | 0 invalid tets | QGP3D 2296 hexes | **timeout at 1 h** in the final Ipopt stiffening (exit 124): inconclusive, slow |
| env, n=2000 | 90k | 0 invalid | QGP3D 2303 hexes | **timeout at 1 h** (same stage): inconclusive, slow |
| **baffle, n=1000** | 85k | valid (0/0) | robust pipeline | **1000 hexes, 0 inverted, min sJ 0.40, 0 interior singular edges**, every boundary quad on the input surface, baffles 100 % covered by hex faces, 47 min |
| **baffle, n=2000** | 85k | valid (0/0) | robust pipeline | **1800 hexes, 0 inverted, min sJ 0.82, 0 singular edges**, baffle coverage 100 % / 109 %, 56 min |
| crack, n=2000 | 85k | **175 invalid tets, 62 bad valencies** | "Parameterizing without quantization" | unusable. ~50k "non-meshable footprint at boundary singular edge" warnings: AlgoHex cannot handle a 360° *boundary* edge |

**Block count of the baffle result** (`slab/split_blocks.py`): the hex mesh is split along the hex faces on the
baffles; the tips stay shared. The 360° tip edges and interior valence ≠ 4 edges are taken as block-defining,
then the base complex is computed. The result is **13 blocks** for both n=1000 and n=2000, with identical
topology and sizes such as [360, 300, …, 12]. This is the expected near-minimal H-type layout for a
rectangle with two staggered slits, swept in y. No sheet collapse was needed.

What this shows:
1. **AlgoHex accepts interior zero-thickness feature surfaces with free edges inside the volume.** The tips are
   regular (4 hexes around each), and the result is clean. The 134 "non-manifold feature surface at edge"
   messages, at the baffle–wall junctions where three feature sheets meet, were harmless.
2. The same baffle given as an explicit crack (a boundary with 360° edges) breaks the parametrisation. **The
   baffle must stay interior** in the AlgoHex input and be split only afterwards, on the hex/block mesh.
3. Removing the thickness length scale removes the multi-scale problem. The baffle runs have *zero*
   singular edges and give the minimal block structure directly. The two thick-wall variants did not finish
   in 1 h at the same n.

Caveat: this is a synthetic, planar test. The real rings are curved, elliptical, and meet the symmetry planes
obliquely. One real canada run is needed (section 6, step 2).

### 4.2 Sweepability of the draft tube [V]
`blocking_study/sweep_sections.py` (section 2): every meridional section has the same topology (1 face, 2
holes, same label sequence). `sweep2d/section.py` writes the sections at φ = 10° and 80° as planar BREPs:
both have 25 boundary edges and the same corner count. The naive 2D quad-mesh layout (`sweep2d/quad2d.py`)
gives 751–2087 patches. That only shows a raw quad mesh is not a layout; it is not a block-count estimate.

### 4.3 Not run (would exceed the 1 h limit)
- canada draft tube with ring **mid-surfaces as interior baffles**. The frame-field phase alone took ~3.2 h
  on the 172k-tet cad1 input.
- gv passage with the **blade as a camber-surface baffle**.

---

## 5. Method comparison

| Method | Thin gaps (1-cell regions) | Thin rings / knife edges | Few blocks for TFI | Generality | Local availability | Judgement |
|---|---|---|---|---|---|---|
| **Current**: envelope cut-out → AlgoHex → beam collapse | ✗ the offset narrows the gap; uniform h | ✗ hole + O-layer re-insertion | ✗ 109–330 on dt | high in principle | yes | fails on canada **[V/R]** |
| **AlgoHex with zero-thickness baffles** (walls kept as interior feature surfaces, split + inflate + pillow O-layer at block level) | ✓ full gap kept (0.099 gv, 0.20–0.24 dt mid-surfaces) | ✓ crack tip = 4 sectors, regular | ✓ 13 = minimal on the slab, 0 singularities **[V]** | high: any thin wall with a definable mid-surface (blades, rings, splitters, struts) | yes: AlgoHex + gmsh fragment + small Python | **most promising for 3D** |
| Multi-domain decomposition (gv / ru / dt / inlet as subdomains with conforming interfaces) | ✓ when interfaces are placed where there is room | depends on the per-domain method | ✓ | high | partly: gmsh/OCC split, interface matching to write | needed anyway as the frame |
| **2D blocking + sweep (2.5D)**: span sweep for blade rows, φ sweep for the draft tube | ✓ 2D sizing is local and cheap | ✓ 2D corners force valence ≥ 1; thin rings are C/O-grids in 2D | ✓ dt ≈ 2D count + axis butterfly | medium: needs a sweep direction with constant section topology (verified for dt) | gmsh 2D + own layout simplifier | **robust for the canada components**; needs a 2D layout simplifier |
| dtOO parameter-space template (H-O-H in (m, θ, span)) | ✓ | ✓ | ✓ fixed | low: only dtOO passages | yes | good fallback / reference, teaches no topology variety |
| Polycube (PolyCut, Evocube) | ✗ the same global-scale problem; thin slots distort | ✗ knife edges and twisted blades distort | ~ often not minimal | high | ✗ not installed, no downloads | not recommended now **[J]** |
| Medial-axis decomposition | ✓ the medial axis captures slots | ~ unstable at knife edges (branches) | ✓ in 2D, hard in 3D | medium | ✗ nothing installed | useful only for 2D sections **[J]** |
| Dual / sheet simplification (Gao 2015/17, beam collapse) | – (post-processing) | – | ✓ only if the input is clean | high | yes (`beam_collapse.py`) | keep it as the last step, on clean input only |

---

## 6. Recommendation and next steps

**General approach: "walls as baffles" + subdomain decomposition, with 2.5D as the robust fallback for
sweepable parts.**

1. **Represent every thin wall as a zero-thickness interior feature surface, never as a cut-out hole.**
   - Rings: use the mid-surface between the inner and outer skin. dtOO already builds the rings from a
     meridional mean line plus thickness, so the mean line revolved is the exact mid-surface.
   - Blades: use the camber surface (dtOO has it).
   - Feed AlgoHex the tet mesh with these surfaces embedded (gmsh `fragment`, interior `face_colors`, free
     edges as feature edges, feature vertices where a free edge meets a wall).
   - After extraction: split the hex/block mesh along the baffle faces, keeping the tip edges shared. Then
     morph each split sheet to the real pressure/suction or inner/outer surface.
   - Insert the O-/C-layer as a **pillowed sheet at block level** before the TFI (fits the user rule "blade
     before TFI": the O-layer is part of the block structure).
   - Generic: it uses no labels, only "surface is thin wall → mid-surface".
   - Effort: 3–5 days. Mid-surfaces from dtOO, the gmsh embed writer (prototype `slab/gen_slab.py`), the split
     (prototype `slab/split_blocks.py`), and morph + pillow (new, the main work).
2. **Decisive test on the real draft tube (needs approval, ~4–5 h on 1–2 cores):**
   - ring and lip mid-surfaces as baffles in the cad1 domain (no envelopes)
   - n=20000; the mid-surface slots of 0.20–0.24 are ≈1 h
   - check the seamless validity first: abort on an invalid seamless map
   - then split + base complex; compare with the 15–25 hand estimate.

   Same for the gv with a camber baffle (~15–60 min, 1 core).
3. **Do not truncate at artificial interfaces close to the LE.** Either place the gv/ru interfaces in the gap
   centre between the rows and LE planes, or block gv-LE-to-ru-TE jointly. Rotor–stator interfaces may stay
   non-conforming (AMI/GGI) in CFD, which decouples the rows.
4. **2.5D path for the draft tube, as fallback and as reference.**
   - The φ-sweep is verified.
   - Block the meridional section in 2D. Rings become slits (baffle idea in 2D) or thin C-grids.
   - Use a quarter-butterfly core at the axis ahead of the nose and behind the hub end.
   - Sweep through φ with exact projection onto the elliptic surfaces.
   - Needed: a 2D layout simplifier that turns the quad mesh into a minimal layout (motorcycle graph +
     chord collapse, or a 2D IGM with 2D sheet collapse via the existing beam code).
   - Effort: 3–4 days. This also gives a reference minimal count for the dataset.
5. **Keep beam collapse only for clean inputs** (0 invalid seamless, few singular edges). On canada-type
   input, invest in the representation (baffles), not in the search depth.

Risks:
- Curved, oblique baffles meeting the symmetry planes. Three feature faces meet there; this worked on the
  slab.
- The baffle tip on a blunt LE (round LE): fine after inflation, the C-grid wraps around it.
- AlgoHex runtime on 172k tets (~3–4 h) remains. A coarser input tet mesh is possible because the slots
  are no longer the size driver.

---

## 7. Files
- `blocking_study/gaps.py`: label gap / volume / h measurement
- `blocking_study/sweep_sections.py`, `sweep_sections2.py`: meridional section topology of the dt domain
- `blocking_study/slab/gen_slab.py`, `crack.py`, `eval_hex.py`, `split_blocks.py`: synthetic test
- `blocking_study/slab/{walls,env,baffle,baffle1k,crack}/`: inputs, logs, metrics
  - `hex_n*.ovm`, `*_eval.json`, `*_split.json`, `*_split.vtk` (block ids in `*_split_blockid.npy`)
- `blocking_study/sweep2d/`: φ=10°/80° sections (BREP, edge labels, quad mesh) and the naive layout count
- `blocking_study/user_*.txt`, `asst_*.txt`: transcript extracts used for section 1

---

## 8. Draft tube baffle test (2026-11-13)

Directory: `/work/canada/blocking_study/dt_baffle/`. Tags as above: **[V]** measured here, **[J]** judgement.

### 8.1 Input [V]
- `build_parts.py` (pythonOCC, `/opt/venv-dtoo`, `LD_LIBRARY_PATH=/dtOO-install/lib`):
  - Thin walls are detected **generically** among the curved faces of `case/cad_blk.stp`: pairs whose mutual
    distance is small and nearly constant. Found (0,3) lip/ring 1, (1,4) ring 2, (2,5) ring 3; median gaps
    0.110 / 0.075 / 0.075. No label or face numbers are used for this.
  - Mid-surface: skin-A grid projected onto skin B, midpoints, each streamwise column resampled to uniform
    arc length. This step is needed: without it the knife-edge clamping bunches rows to 0.005 and gmsh hangs
    in MeshAdapt. The grid is mirrored across x=0/y=0 and fitted with a B-spline.
  - The lip mid-surface continues upstream as the mid-line between the shroud (line between the interface
    rims) and the mantle. It closes the hub-chain obstacle (bulb + gv/ru hole + dt hub; faces as in
    `cad/domain_blk.py`). The ring-1 part continues it as an interior baffle.
  - Outer body: `domain_blk.outer_body()`. Interface faces: extended radially and trimmed at the lip.
- `assemble_mesh.py in_a --h-baffle 0.08` (gmsh 4.15 OCC, `python3`). Steps: cut to the quarter, fragment,
  remove the hub-chain volume. Baffles stay **embedded interior faces**. Mesh size: 0.08 at the baffles,
  0.4 max, 0.025 min. Feature edges: label change, crease > 25°, free baffle edges, and > 2 faces.
- `check_baffle_input.py in_a` → **PASS** (`in_a/input_check.json`):
  - 28030 nodes, 152083 tets, all positive; volume 259.87
  - 11769 baffle triangles, **all interior (2 tets), 0 on the boundary, 0 duplicate nodes** (min node
    distance 0.023)
  - closed 2-manifold boundary
  - 346 free baffle edges + 223 baffle/wall junction edges, all feature edges
  - 942 feature edges, 29 feature vertices, 0 dangling ends, 0 unmarked kinks
  - `check_mid.py`: baffle nodes sit at the skin midpoint within 1.5e-3 (median 1e-4)
- Labels: `in_a/labels.txt` (1-3 = BAFFLE_LIP / BAFFLE_MID_1_4 (ring 2) / BAFFLE_MID_2_5 (ring 3)).

### 8.2 AlgoHex runs [V]
`run_ah.sh <dir> <n> <cpus> <threads> <timeout>` runs the patched `/work/datagen_keepbl/algohex_fix/HexMeshing`
with `-i tet.vtk -n <n>`:
- n=2000: cores 0-3,96-99, OMP 4
- n=20000: cores 4-7,100-103, OMP 4, `--save-param-state in_a/state_n20000.ovm`

Both were started 10:02/10:13 and both took ~2.6 h (n=20000: integrability 0.6 h, parametrisation 1.75 h).

| | n=2000 | n=20000 | cad1 envelopes n=20000 (ref.) |
|---|---|---|---|
| seamless invalid tets / valencies | **0 / 0** | **0 / 0** | 511 / 0 |
| final invalid tets, valid volume | **0, 1.0** | **0, 1.0** | 7156, 0.999 |
| quantisation | QGP3D, 2326 hexes | QGP3D, 20845 hexes | failed → holes |
| inverted hexes (sJ ≤ 0), min sJ, sJ p1 | 17, -0.99, 0.05 | **61**, -0.97, 0.49 | – |
| boundary: max vertex dist / quads off surface | 0.017 / 5 | 0.024 / **0** | holes up to 1.2 |
| baffle coverage lip / ring 2 / ring 3 | 1.00 / 1.00 / 0.90 | 0.98 / 1.00 / 0.94 | – |
| interior singular edges (curves) | 143 (9) | 314 (9) | 634 @ n=100k |
| base-complex blocks, unsplit | 158 | **218** | 109-330 after beam (n=100k) |
| blocks after split along baffles (`eval_dt.py`) | 469 (319 ≤ 4 cells) | **932** (292 ≤ 4 cells) | – |

Coverage counts an interior hex face as on a baffle when all 4 vertices lie within 0.02·h of it. The planar
quad centroids sag on the curved sheets. `eval_dt.py` still gives 13 blocks on the slab baffle mesh
(`slabtest/`).

### 8.3 Diagnosis [V for locations, J for causes]
`in_a/singular_curves_n*.txt`, `diag.py`, `diag2.py`.

The singular structure is the **same 9 curves at n=2000 and n=20000**. The frame field fixes it, not n. The
pairs have identical length (3/5 pairs):

1. **A 3/5 pair through the whole domain** (z -4…10) in the bypass flow, at φ ≈ 15-21°, r 2.8-7.8.
2. **A val-3 curve in the diffuser** at φ ≈ 17-20° (z 0.2…10, r 1.8-3.8), and a val-5 curve near the hub/axis
   (r 0.5-0.76, φ ≈ 35°).
3. **A val-3 curve on the axis behind the hub end** (r 0.2-0.27, z 2.28…10), starting at the rim of the
   non-axisymmetric hub-end disc.
4. **A 3/5 pair along the whole ring-3 trailing edge** (z 5.98-6.41, φ 0…90°). Here the ring-3 mid-surface
   turns from radial (φ=0, r→6.1) to axial (φ=90°).
5. **A val-3 curve from the inlet to the lip nose** (r ≈ 2, φ 40°), and a short curve around the bulb.

The 61 inverted hexes sit:
- ~40 in a chain along curve family 1-2 at φ ≈ 17°, through the whole length (z -3.9…9.9);
- ~17 at the ring-3 trailing edge (z 6.24-6.37);
- 4 near φ ≈ 2-8°, r 5.9, z 6.0 (ring-3 TE at the y=0 plane).

None are at the ring-1/ring-2 leading or trailing edges or at the lip/interface junction. The baffle tips
themselves are fine apart from the twisted ring-3 TE.

[J] The long 3/5 pairs follow from the domain shape, not the walls:
- The quarter section goes from round to strongly elliptical: far field r 3 at φ=90° vs. up to 11 at φ=0,
  a 30° cone at φ=0 and a cylinder at φ=90°.
- The rings flare mainly at φ=0.
- The frame field, held by both symmetry planes, twists about the circumferential direction and relaxes
  through full-length singularity pairs.

Because the pairs run the whole length and cross every sheet, the base complex (218) and the split (932) are
far from the 15-25 estimate. Beam collapse was **not run**: the output is not clean (61 inverted).

### 8.4 Result and judgement
- **[V]** The representation works for validity. Embedding the rings and lip as interior baffles removes the
  failure mode of all earlier draft-tube runs: 0 invalid tets in the seamless and final map, at both n,
  versus 511 / 7156. Quantisation succeeds and the boundary fit is 0.024. It also runs 1.5× faster than cad1
  (2.6 h vs 4.1 h).
- **[V]** It does **not** give a minimal block structure on this domain: 9 frame-field singular curves (3/5
  pairs) through the whole length, 61 inverted hexes, 218 / 932 blocks.
- **[J]** The slab result (0 singular edges, 13 blocks) does not carry over. The real domain twists in φ, and
  no feature change targeting the thin walls can remove that. So no corrective AlgoHex variant was run.
  Candidates for a next test, none of them verified:
  - a cylindrical far field, to test the twist hypothesis directly (~2.6 h);
  - a full 360° / half model without the y=0 plane;
  - the 2.5D φ-sweep (§6.4). The verified sweepability and this twist make the 2.5D sweep the stronger path
    for the draft tube. The baffle idea still holds inside the 2D layout (rings as slits).

### 8.5 Files (ParaView)
- input: `in_a/tet.vtk` (cell data `color` = label, see `labels.txt`), `in_a/domain.brep`, `parts/*.brep`
  (mid-surfaces, body, interfaces)
- hex: `in_a/hex_n{2000,20000}.ovm`, `metrics_n*.json`, `hexmeshing_n*.log`
- evaluation:
  - `in_a/hex_n*_eval.json`
  - `in_a/hex_n*_split.vtk` (cell data `block`, `sJ`; threshold sJ ≤ 0 shows the inverted hexes)
  - `in_a/hex_n*_baffle_quads.vtk`
  - `in_a/hex_n*_singular.vtk` (singular edges as lines)
- saved state for resuming at another n: `in_a/state_n20000.ovm`
- scripts: `build_parts.py`, `assemble_mesh.py`, `check_baffle_input.py`, `check_mid.py`, `run_ah.sh`,
  `eval_dt.py`, `diag.py`, `diag2.py`

---

## 9. 2.5D phi-sweep test (2026-11-13)

Directory: `/work/canada/blocking_study/sweep_phi/`. Tags as above: **[V]** measured here, **[J]** judgement.
Input: the baffle domain of §8 (`dt_baffle/in_a/domain.brep` = 1 solid + 3 embedded zero-thickness faces for the lip,
ring 2 and ring 3). So the rings and the lip are **zero-thickness slits** (their mid-surfaces). Slits worked, so no
C-grid fallback was needed. The real wall thickness (0.075 / 0.11) is not in this mesh.

### 9.1 Sections [V]
`section_rz.py` (pythonOCC) cuts the solid and the free faces with the half plane φ, for φ = 0, 2.5, …, 90°
(37 sections; at 0/90° the plane is offset by 1e-4° so it does not coincide with the symmetry face). Each section
edge is tagged with the 3D face it lies on, and edges with r≡0 are tagged as the axis. These tags are used only to
check correspondence, never for decisions. All 37 sections have the same sequence: 12 loop edges + 3 slits, one of
them (the lip) attached to the loop. The quantised corner angles are identical at every φ.

### 9.2 2D layout of the section [V]
`layout2d.py`, `collapse2d.py`. A generic rectilinear (2D polycube) abstraction plus full-extension separatrices
and a guarded chord collapse:
1. The attached slit is inserted into the outer loop (out and back). Free slits are their own loops. Every segment
   is averaged over φ (normalised arc length) to give a **mean section**. The layout is built once on that mean
   section.
2. Each segment gets the rounded tangent direction {±r, ±z} (10° hysteresis). Where the rounded direction changes
   inside one CAD edge, a **virtual corner** is inserted. Exactly one is found: on the bulb nose, which turns 90°
   smoothly. Without it the corner indices sum to 5 instead of 4.
3. Corner valences come from the turn of the rounded direction (v = 2 − Δ/90):
   - 9× v=1;
   - 2× v=2: hub cone kink, lip junction;
   - 3× v=3: interface/lip wall, hub-end/cone, bulb virtual corner;
   - 5× v=4: the lip tip and the 4 ring tips.

   Turning sum: outer +360°, each slit −360°. So **no interior singular point is needed**.
4. Parameter coordinates: vertices joined by an r-segment share z, and vertices joined by a z-segment share r
   (union-find; the value is the mean physical coordinate).
5. From every vertex with v ≥ 2, v−1 rays run until they hit the boundary or a slit, giving a T-junction-free
   arrangement whose faces are the blocks. This gives **36 blocks, 23 rays**.
6. Chord collapse: two neighbouring parallel coordinate groups are merged (the chord between them is removed) if all
   of the following hold:
   - the parameter polygon stays valid;
   - the block count drops;
   - the Winslow re-fit at φ = 0, 10, 20, 30, 45, 60, 90° (h = 0.15) has 0 inverted quads and min sJ ≥ 0.1.

   Accepted: hub-end z ≈ lip-tip z (36 → 34). Also accepted, without changing the count: the two axis groups, and
   two r-lines 0.012 apart. Rejected:
   - ring-3 LE ≈ ring-2 TE (33 blocks): inverted at φ=60° (min sJ −0.04), and 0.076 at φ=0;
   - ring-2 LE ≈ lip tip (33): min sJ 0.025–0.08;
   - ring-2 TE ≈ ring-3 TE (30): inverted, −0.36 at φ=0.

   r-merges are impossible because the slits overlap axially. A guard of 0.03 gives the same 34
   (`collapse_g003.log`).
7. **Result: 34 2D blocks, 0 interior singular points.** All irregular vertices are boundary corners or slit tips.
   Every CAD vertex and slit tip is a block corner. Re-fit per φ (`Mesh2D.solve`): boundary nodes by normalised arc
   length per segment, interior by Winslow in the index space of the parameter grid, warm-started from the
   neighbouring φ. **0 inverted quads at all 37 φ.** min sJ is 0.53–0.67 for φ ≥ 25° and drops to **0.059 at φ=0**.
   The worst cell is a single cell at the ring-3 trailing edge, where the ring flares radially (r 6.08, z 6.3); the
   next worst is 0.29.

### 9.3 Sweep and axis [V for the build, J for the necessity argument]
`sweep3d.py`:
- Every 2D block becomes a hex block spanning φ, with its nodes taken from the same 2D grid in all 37 sections.
- Block edges and faces along φ lie on the section curves, i.e. on the CAD faces, at every φ level.
- **Axis:** each of the 4 axis-adjacent 2D blocks becomes a **quarter butterfly**:
  - a core square block (side 0.45 of the block's radial extent);
  - two half blocks (φ 0–45° and 45–90°) between the square and the outer arc.

  The core has a 90° corner on the axis, so there are no degenerate cells. The positions are interpolated in (s, φ)
  from the section grid lines.
- [J] A conforming butterfly needs a block vertex on the outer arc at 45°. That split cannot end inside the domain
  without further singular edges, and every swept block connects the two symmetry planes. So **all swept blocks are
  split at φ = 45°**. dtOO does the same: its 94 hex blocks are 47 φ pairs, split at ~35°.
- **3D: 2·30 swept + 4 cores + 8 halves = 72 hex blocks** (dtOO: 94 hex + 4 unstructured tet/pyramid cores = 98).
  Untested alternatives:
  - 1:2 block faces at the 4 butterfly arcs: 42 blocks, but non-conforming at block level;
  - a full-φ sweep without a butterfly: 34 blocks, but degenerate axis wedges (sJ = 0).
- Fill: every block is re-filled by **3D TFI from its 6 faces**. Resolution: h2d = 0.1, 36 φ cells; 217 908 hexes,
  233 439 nodes; run time 264 s on 8 cores.

### 9.4 Evaluation [V]
| | sweep (stacked sections) | **TFI-filled** |
|---|---|---|
| inverted hexes (sJ ≤ 0) | 0 | **0** |
| min sJ / p1 / median | 0.059 / 0.438 / 0.946 | **0.059 / 0.435 / 0.946** |
| worst blocks | b32lo 0.059 (ring-3 TE, φ=0 face), b3core 0.20, b23lo 0.29 | same |

- Conformity:
  - 140 shared block-face pairs, 152 block faces on the boundary, **0 non-conforming block faces**;
  - quad faces: 638 383 interior (2 hexes), 30 682 boundary (1 hex);
  - **0 boundary edges with an odd face count** (watertight; the baffles are two coincident boundary sheets);
  - shared-node mismatch 3.6e-15.
- Hausdorff to the CAD (`hausdorff_occ.py` + `hausdorff_mesh.py`; 6000 nodes + 6000 quad centres → CAD, and
  61 487 CAD samples → mesh, exact point-triangle distance):
  - nodes max **2.4e-3**. The maximum is a butterfly node on the bulb nose, interpolated and not projected; all
    section nodes lie on the OCC section curves.
  - quad centres max 0.022, p99 6.5e-3.
  - CAD → mesh max **0.024**, p99 3.8e-3.
  - **Hausdorff 0.024** (0.17 % of the 14 domain length).
  - The maximum is on the far field at φ ≈ 21°, r ≈ 8.2, z ≈ 8.8, where the far field changes from cone to
    cylinder between two 2.5° sections. Finer φ steps would reduce it.

### 9.5 Judgement [J]
- **The 2.5D φ sweep is viable for the draft tube.** It is all-hex, conforming, 0 inverted cells, exact CAD faces at
  every section, and has 72 blocks. That is fewer than dtOO, without unstructured axis cores. It runs in minutes,
  not hours. The φ twist that broke AlgoHex (§8) does not matter here, because the layout is built in the section
  and only the geometry varies with φ.
- **Weak spot: the φ=0 plane.** Ring 3 flares to nearly radial there and the far field is a 60° cone. This is what
  limits the chord collapse (34 instead of 30–33) and gives the min sJ of 0.059. Fixes would be sliding boundary
  nodes on the slits (optimised, not arc-length), or a quality-optimising smoother instead of plain Winslow.
- **Not done:**
  - thickness of the rings and lip: the slits must be inflated to thin C/O-grids at block level, i.e. pillowing
    around the mid-surface slit (rule "blade before TFI");
  - interfaces to gv/ru;
  - projection of the butterfly end-face nodes.
- **Generality:**
  - needs a sweep direction with constant section topology;
  - the rectilinear abstraction works without interior singularities only when the boundary turning quantises
    consistently (true here after one virtual corner). Other sections may need interior 3/5 pairs, which are not
    implemented;
  - the chord collapse is greedy;
  - the mid-surfaces come from `dt_baffle/build_parts.py`, which uses one fixed CAD face index (the mantle,
    `faces[11]`) for the lip extension. That index is inherited here, not added.
- Fixed values used:
  - φ step 2.5°, split at 45° (the middle of the quarter);
  - h2d 0.1 (guard h 0.15), guard min sJ 0.1, guard φ set as above;
  - hysteresis 10°, butterfly core a = 0.45;
  - slit attachment tolerance 1e-4, symmetry-plane offset 1e-4°.

### 9.6 Files (ParaView)
- `sweep_tfi_blocks.vtk`: TFI hex mesh, cell data `block`, `block2d`, `kind` (0 swept, 1 core, 2 half), `sJ`
- `sweep_noTFI_blocks.vtk`: the same before TFI (stacked sections)
- `sweep_block_edges.vtk`: the 12 edges of every 3D block (cell data `block`)
- `layout2d_phi{00,45,90}.vtk`: 2D section meshes (cell data `block`, `sJ`)
- pictures:
  - `layout_final.png`: parameter layout and φ = 0/45/90;
  - `layout_initial.png`: the 36-block layout before collapse;
  - `sections.png`;
  - `blocks3d.png`.
- metrics: `sweep_metrics.json`, `hausdorff.json`, `sweep3d.log`, `collapse_g010.log` / `layout_final.json`
  (= `layout_g010.json`), `collapse_g003.log`
- scripts: `section_rz.py`, `layout2d.py`, `collapse2d.py`, `sweep3d.py`, `hausdorff_occ.py`, `hausdorff_mesh.py`,
  `plot2.py`, `plot3d.py`, `corners.py`, `dtoo_ref.py` (φ spans of the dtOO blocks)
- sections: `sections/sec_*.json`

---

## 11. Block-level sheet collapse (2026-11-13)

Directory: `/work/canada/blocking_study/sweep_phi/block_collapse/`. Test object: the 72-block sweep structure of §9
(`sweep_tfi_blocks.vtk`). Tags: **[V]** measured, **[J]** judgement.

### 11.1 Method (generic code) [V for what it does]
- `cad_export.py` (pythonOCC): reads any .brep/.step and writes its unique faces/edges/vertices (TopExp maps),
  dense face and edge samples, and the topology (edge→faces, vertex→edges). `domain.brep` has 15 faces,
  42 non-degenerate edges and 31 vertices.
- `bstruct.py`: reads a hex VTK with a per-cell `block` id. The (i,j,k) index of each block is recovered by walking
  face neighbours, with no assumption about the writer. Also contains the coarse mesh (one hex per block), the dual
  sheets (classes of parallel block edges), TFI and parametric resampling, and closest-point queries.
- `block_collapse.py --vtk --cad --out [--h --sj-floor --htol --width --depth]`:
  1. **Classification.** Each block face is assigned a CAD face by majority vote of its fine quads. Nodes and edges
     are assigned to a CAD vertex, CAD edge, CAD face or the interior by distance (tol 0.01) and by consistency
     with the CAD closure relation. The fine boundary quads of the input, grouped per CAD face, serve as the CAD
     surrogate for projection.
  2. **Topological collapse.** For each sheet, the two ends of every sheet edge are welded. The collapse is
     rejected if any of these apply:
     - the sheet intersects itself, or touches itself at a node;
     - a sheet block has boundary faces on both sides (the region would vanish);
     - welded nodes, edges or faces lie on CAD entities that are not in each other's closure (e.g. two CAD
       vertices welded, or a node on a CAD edge welded onto a vertex not on that edge);
     - two boundary faces are glued, which would make a CAD face interior;
     - a face is non-manifold, a doublet appears, the boundary develops an odd or new non-manifold edge, or an
       edge fan splits (all checked relative to the input, so the lip/wall junction is allowed);
     - a CAD vertex is lost, or the block edges on a CAD edge no longer form one chain.
  3. **Geometry.** Each welded node sits at a fraction τ of its sheet edge. The lower-dimensional CAD class wins;
     if a block face is on the boundary, that side wins; between equal classes τ = 0.5, and the variants lo/hi put
     the node at one end on a consistently oriented sheet. Geometry is then built as follows:
     - merged block edges = old edge + the sheet-edge part (the union of the old curves), or the τ-curve of the
       sheet face;
     - merged block faces = the τ-surface of the sheet block;
     - a block that absorbs a half sheet block is concatenated with it. If pieces attach on two different axes,
       the block falls back to plain TFI.

     Everything is then resampled in index units of the old fine cells. Interval count per new sheet = max over
     its edges, so the input resolution is kept; `--h` gives a uniform size instead. Faces and blocks are
     TFI-corrected to their new boundary, and boundary nodes are projected onto the CAD surrogate.
  4. **Guard.** 0 inverted hexes, min sJ ≥ **0.05**, and Hausdorff no more than 0.01 above the re-fit baseline in
     both directions:
     - CAD samples → boundary;
     - boundary quad centres → surrogate.
  5. **Beam search.** Width 4, depth up to 8, best state anywhere; every sheet is tried with each placement
     (mid/lo/hi).
- `replay.py` replays a path (any resolution) and checks conformity. `hausdorff_fwd_occ.py` gives the exact
  mesh→CAD distance. `plot_section.py` draws the φ=0 picture.
- **Re-fit sanity [V]:** the identity re-fit through the same pipeline reproduces the input: 217 908 hexes,
  0 inverted, min sJ 0.0588, max node change 2.4e-3 (the bulb-nose butterfly node is now projected onto the CAD).
  An earlier arc-length resampling variant, and the uniform-h mode at h=0.1 (2.6 M cells), gave 59 inverted hexes
  in a butterfly half *without any collapse*. So index-unit resampling is required, and uniform h is not reliable
  here.

### 11.2 Sheet inventory of the 72-block structure [V]
22 sheets (Σ members = 216 = 3·72):
- **2 φ layers** (lo/hi, 38 blocks each = 30 swept + 4 butterfly halves + 4 cores). These are the sheets of the
  φ=45° split; each core belongs to both.
- **14 swept 2D chords**, each = lo+hi copy of a 2D chord (16, 16, 14, 10, 8, 6, 6, 4, 4, 2, 2, 2, 2, 2 blocks).
- **4 swept 2D chords that continue through a butterfly** (core + 2 halves): 11, 11, 9 and 7 blocks.
- **2 butterfly-internal radial sheets** (halves only, 6 and 2 blocks).

### 11.3 Results [V]
At depth 1 (input resolution), out of 22 sheets:

**Topology: 15 rejected.**
- 5× a sheet block lies between two boundary faces. These are chords from wall to wall or wall to slit in one
  block.
- 3× two CAD vertices would be welded.
- 4× forced boundary side incompatible: the welded node would leave its CAD edge. Both φ-hi and the "r-chords"
  next to the hub-end disc fall here.
- 2× both sides of the ring-3 baffle would be glued.
- 1× φ-lo: the CAD edge 31 (hub-end disc rim) would lose its last block edge.

  So both φ-layer collapses (= dropping the 45° split) are rejected topologically. The rim vertex at φ=45° would
  be welded onto the symmetry plane; the degenerate axis (§9.3) would come on top of that.

**Geometry: 6 rejected** (best placement shown):

| sheet | inverted hexes | min sJ | Hausdorff |
|---|---|---|---|
| 2D 5/12/20/28 | 419 | -0.69 | – |
| 2D 1/8/16 (butterfly) | 85 | -0.43 | h_rev 0.043 |
| 2D 23/31 | 582 | -0.43 | h_rev 0.040 |
| 2D 14/22/30 | 5 | -0.03 | h_rev 0.038 |
| butterfly radial, 6 halves | 0 | **0.026** | h_fwd 0.0188 > 0.0183 |
| butterfly radial, 2 halves | 0 | **0.043** | h_fwd 0.0188 |

The two butterfly radial sheets remove the half blocks. The core then reaches the arc directly, with a flat 180°
corner at φ=45°, so they fail only on quality and Hausdorff.

**Admissible: 1.** The single-block chord 2D block 15 (b15lo/hi), i.e. the thin strip (z 2.08–2.38) between the
lip tip (face 12) and the ring-2 leading edge (face 13). Its collapse joins the lip-tip line and the ring-2 LE line
into one inclined block edge (`final/section_phi0.png`). This is the 2D merge "ring-2 LE ≈ lip tip" that §9.2 had
rejected at the Winslow guard 0.1 (0.025–0.08).

**Depth 2: nothing admissible. Final: 70 blocks (72 → 70).**

| | sweep re-fit (72) | **block collapse (70)** |
|---|---|---|
| hexes | 217 908 | 221 436 |
| inverted | 0 | **0** |
| min sJ / p1 / median | 0.0588 / 0.435 / 0.946 | **0.0588 / 0.435 / 0.944** |
| worst blocks | ring-3 TE 0.059, core 0.20 | ring-3 TE 0.059, **merged neighbours 0.088 / 0.110** (b16hi/lo) |
| conformity | – | shared block faces identical (0.0), quads {1: 29 594, 2: 649 511}, 0 odd boundary edges |
| CAD→mesh (all 268 262 CAD samples) | 0.0268 | 0.0268 |
| mesh→CAD exact OCC (6000 nodes / 6000 quad centres) | – | nodes 2.4e-3, centres **0.0225** (p99 7.2e-3) |

Hausdorff ≈ 0.027, the same as the sweep. The 0.0268 vs 0.024 of §9.4 comes only from the denser CAD sampling
here; the maximum is the same far-field cone/cylinder spot between two φ sections. The beam at a uniform h=0.15
(800 k cells) gives the same single admissible collapse and the same 70.

### 11.4 Judgement [J]
- On this structure, block-level collapse saves **only 2 of 72 blocks**. The sweep layout was already
  chord-collapsed in 2D (§9.2), so almost every remaining sheet either runs from wall to wall in one block or
  would weld distinct CAD vertices or edges. Those are hard topological limits, not quality issues. The 45° split
  cannot be removed by collapse; that would need a different axis treatment.
- The geometric rejections come from local TFI without smoothing. The near-misses (2D 14/22/30: 5 inverted,
  -0.03; butterfly radial sheets 0.026/0.043) might pass with a smoother or optimised boundary sliding, but the
  butterfly ones create a 180° corner by construction. Running a smoother would add at most 6 more blocks (66).
- **As a general step it is worth keeping, but as a cheap cleanup and verification pass, not as a major
  reducer.** It:
  - runs in about 45 s for the full beam;
  - is generic (any conforming block VTK + CAD export);
  - catches exactly the sliver blocks a layout generator leaves (here the 0.3-high lip/ring-2 strip);
  - explains *why* each remaining sheet is necessary.

  The large reductions must come from the layout itself (2D chord collapse, axis treatment). For AlgoHex base
  complexes, with many sliver sheets, the yield should be higher, but that is untested. Next candidate: the C-grid
  variant via `--vtk`.
- Limits:
  - classification uses distance tol 0.01 plus a majority vote. One edge conflict (the non-axisymmetric hub-end
    disc: CAD vertex 30 / edge 30 are not resolved by the input blocks) is inherited and reported;
  - fallback TFI for pieces on two axes;
  - no smoothing;
  - the Hausdorff guard is measured against the input boundary as surrogate, and the final check is the exact
    OCC check.

### 11.5 Files (ParaView)
- `final/best_h0.0_tfi.vtk`: 70-block TFI mesh, cell data `block`, `origin` (= input block id), `sJ`.
- `final/best_h0.0_block_edges.vtk`, `final/best_h0.0_coarse.vtk`: one hex per block, point data `node_dim` and
  `node_cad`.
- `final/section_phi0.png`; `final/best_h0.0.json` (path + metrics); `final/best_h0.0_hausdorff_fwd.json`;
  `final/replay_best.log`.
- Baseline re-fit: `run_hin/baseline_h0.0_{tfi,block_edges,coarse}.vtk`.
- Full log of every sheet/variant with reasons: `run_hin/block_collapse.log`, `run_hin/records.json`.
- `cad.npz` (CAD export), scripts as in 11.1.

---

## 10. C-grid ring thickness (2026-11-13)

Directory: `/work/canada/blocking_study/sweep_phi/cgrid/`. Tags as above: **[V]** measured, **[J]** judgement.
It builds on the §9 layout: 34 outside blocks, unchanged, with the same layout at all 37 φ. No block-level sheet
collapse was run.

### 10.1 Real walls, generic detection [V]
`wall_faces.py` (pythonOCC, `case/cad_blk.stp`):
- **Thin-wall pairs:** same criterion as `build_parts.py`. Found (0,3) lip/ring 1 with gap 0.110, (1,4) ring 2 with
  0.0745, (2,5) ring 3 with 0.075.
- **Wall component:**
  - It is the closure of the CAD faces touching a pair (distance < 1e-4).
  - Faces lying on the baffle-domain boundary are excluded (far field, hub chain, interfaces 18/19; ≥ 50 % of
    the quarter samples lie on a domain face).
  - Result: lip = {0, 3, 8, 11, 12}, ring 2 = {1, 4, 6, 7}, ring 3 = {2, 5, 9, 10}. That is skins, flat end caps
    and, for the lip, the mantle.
- **Mantle criterion, replacing `build_parts.py`'s fixed `faces[11]`:**
  - It is the component face that touches the pair, lies upstream of the pair's minimum z, and extends axially by
    more than 5× the gap (so bevels/caps are excluded).
  - It finds **11**, the same face as the fixed index. The mid-surfaces, and with them `in_a/domain.brep` and the
    §9 sections, are therefore unchanged, and I did not regenerate them.
- **Sections per φ:** each component is cut and chained into one polyline.
  - Rings are closed loops, the lip is an open chain (inner skin → TE cap → outer skin → mantle).
  - Face 12 (an interface annulus crossing the mantle at z ≈ −1.675) does not connect to the chain. It is dropped
    and reported.
  - At φ = 0/90° the wall is cut 0.25° off the symmetry plane, because face seams break the section exactly there.
    The wall is mirror-symmetric there, so the error is 2nd order [J].

### 10.2 Topology: one pillow layer per wall (O around free rings, C around the lip + strip) [V build, J choice]
`cgrid2d.py`. For every slit of the §9 layout, the chain of layout edges along it is collected:
- **Free rings (2):** the closed slit loop, both sides, gets a **closed O-layer**.
  - Every separatrix that ended on the slit continues as a radial layer edge to the wall.
  - At each free tip, the 3 separatrices meet the continuation edge to the wall nose.
  - So the layout stays **conforming and T-junction-free**, with no change outside.
  - Cost: each old slit tip becomes an interior vertex of **valence 6** (3 separatrices + 2 layer sides + 1 radial
    edge). In 3D that is a valence-6 singular edge along φ.
- **Lip (attached):** the chain is extended over collinear v = 2 boundary edges that lie inside the real wall (the
  LIP_WALL mid-curve → mantle). Its two ends:
  - **Inner end (DT_INLET):** convex corner (v = 1). The layer ends with its end edge on the interface.
  - **Outer end (old interface/lip-wall corner, v = 3):** a C-layer cannot end there conformingly.
    - Its end edge is interior, and both alternatives give a non-quad (pentagon) block or a T-junction.
    - So the layer continues as a **strip ("H") along the straight-ahead separatrix** down to the inlet: 2 thin
      strip blocks.
- **Why not C+H everywhere [J]:**
  - A wake/H strip from every tip (5 tips) would run to inlet/outlet and cut every band it crosses. Estimate:
    about +20 thin blocks of 1–2 cells across.
  - The pillow keeps the §9 outside layout as it is and adds only one block per wall-side block edge.
  - The strip is used only where it is topologically forced (the attached lip).
- **Block count:**
  - ring 2: 4 O-layer blocks; ring 3: 3; lip: 5 C-layer + 2 strip.
  - **2D: 34 + 14 = 48 blocks.**
  - **3D: 2·44 swept + 4 butterfly cores + 8 halves = 100 blocks** (§9: 72; dtOO: 98 incl. 4 unstructured cores).

### 10.3 Geometry per φ [V]
- **Layer:** thickness δ = 0.03 with N_L = 12 cells.
  - Geometric wall-normal grading from h1 = 1e-3 × the mean free-ring chord = **2.33e-3**, ratio 1.013. With
    δ = 0.03 and 12 cells this is nearly uniform.
  - δ is limited by the clear gap lip/ring 2 (~0.11–0.17).
- **Wall points:** wall point W_i = normalised arc length between the wall noses; outer node W_i + δ n_i.
  - n is the wall normal, smoothed over ±δ/2.
  - The nose is the farthest wall point along the mid-curve tip tangent within t/2 of it: the cap middle on flat
    caps, the apex on knife edges.
- **Outside grid:** the chain nodes, the duplicated junction and the strip separatrix are moved to these positions
  and become Dirichlet nodes. The outside grid is then re-fitted by Winslow.
- **Fixes that were needed** (each tried, measured, kept only if it helped):
  1. Free slits came from OCC in arbitrary orientation per φ. They are now oriented like the neighbouring section.
     Without this, 37 inverted quads at φ=0.
  2. **Nose finder.** "Nearest wall point to the mid tip" fails at φ≈0, where the ring-3 TE is a long knife that
     extends ~0.07 beyond the clamped mid-surface tip. The tangent-ray rule fixes it.
  3. Block orientation was taken from the summed signed area, not from the first cell. Before this, a whole layer
     block was wrongly flagged inverted.
  4. **Tip radial = outward bisector of the two adjacent wall chords.** At a thin or knife TE the nose layer cells
     must share ~355°, so they are close to 180° each. The bisector equalises them.
  5. **2D local untangler.** Nelder-Mead on free outside nodes of cells with sJ < 0.1, maximising the minimum
     corner sJ. It fixes the last folded outside cell at the pointed ring-3 TE (φ=20°).
- **Tried and dropped:**
  - Knife-edge blunting: cutting the wall where the local thickness < 0.25 t removed up to 0.28 of ring-3 chord
    at φ=0. With fixes 4+5 it is no longer needed. Thresholds 0.1, 0.05 and **0** (final) were all valid.
  - Snapping the tip-neighbour wall node onto the cap corner, to recover the corners. Wall-only snap: ~300
    inverted quads per φ. Full snap: 1 inverted at φ=0 (two coincident grid nodes). Not adopted; the code is kept
    as `cgrid2d_snap_experiment.py.bak`.
- **2D result:** **0 inverted quads at all 37 φ**, min sJ 0.035 (φ=0, ring-3 knife TE).

### 10.4 3D and evaluation [V]
`sweep3d_c.py`: the §9 sweep with the butterfly, φ split at 45°. TFI uses **arc-length-based blending
parameters**. Run time 121 s on 8 cores; 292 212 hexes.

| fill | inverted | min sJ | p1 | median | worst blocks |
|---|---|---|---|---|---|
| stacked sections (no TFI) | 0 | 0.035 | 0.416 | 0.950 | b38lo, b40lo 0.035 (ring-3 TE layer, φ=0 face), b25lo 0.048 |
| plain 3D TFI | **8** | −0.043 | 0.413 | 0.950 | b15hi (5×2-cell outside block at the lip TE / ring-2 LE, φ 70–90°) |
| **TFI + local hex untangling (final)** | **0** | **0.035** | 0.413 | 0.950 | b38lo/b40lo 0.035, b25lo 0.052, b15hi 0.073 |

- The local hex untangling moves block-interior nodes only (block faces are kept). It moved 133 nodes in
  4 blocks: b15hi 18, b25lo 3, b38lo 56, b40lo 56.
- **Plain TFI is not sufficient for one small block.** Index-based TFI also gave 8 inverted hexes there.
- **Conformity:**
  - 208 shared block-face pairs, 184 on the boundary, **0 non-conforming**;
  - quads: 858 799 interior, 35 674 boundary; **0 odd boundary edges** (watertight);
  - shared-node mismatch 1.8e-15.
- **First cell on the walls (after TFI):** 2.330e-3 everywhere (min 2.330e-3, max 2.330e-3), as targeted
  (1e-3 × ring chord).
- **Hausdorff against the real CAD** (`hausdorff_c_occ.py`, `hausdorff_c_mesh.py`). The CAD is all of
  `cad_blk.stp` (walls with thickness, far field, hubs, interfaces) plus the symmetry/inlet/outlet planes.
  - **Mesh → CAD:**
    - wall nodes **2.3e-6** (on the wall section curves);
    - all nodes 2.4e-3 (butterfly bulb node, as §9);
    - quad centres 0.024 overall, 0.021 on walls.
  - **CAD → mesh, walls: max 0.038, p99 0.034.** Cause: the **flat end caps' corners are not mesh nodes**. The
    chain nodes are 0.1 apart, the caps are only 0.075–0.11 wide, so the layer's wall chord from the cap middle
    to the next skin node cuts each cap corner (worst at the lip TE, φ≈11°).
  - Overall CAD → mesh max 0.058 (p99 0.033). The maximum is on the y=0 symmetry face at r 5.26, z 5.74, where
    ring 3 is thin. [J] It is probably a sample inside the real wall at exactly φ=0 that the 0.25°-offset wall
    section misses; this is not verified. Plane samples inside the wall sections were removed (108).

### 10.5 Judgement [J]
- **The 2.5D sweep takes the real wall thickness at block level:**
  - one pillow layer per wall (closed O around the free rings, C around the lip plus a 2-block strip where the lip
    is attached);
  - 48 2D / 100 3D blocks, conforming, watertight, 0 inverted after TFI + local untangling;
  - a BL-type layer with exact first-cell height;
  - wall nodes exactly on the CAD.
- **Weak points:**
  1. **Flat-cap corners are cut** (Hausdorff 0.038). Fix: refine the chain locally near the tips (graded
     parameter intervals at slit tips) so that cap corners become nodes. Snapping on the coarse grid did not work.
  2. **Ring-3 knife TE at φ≈0–20°:** min sJ 0.035 in the nose layer cells. This is intrinsic to a single radial
     line at a knife apex. A C-H (wake) cut there, or local tip refinement, would raise it.
  3. **One small outside block** (b15hi) needs untangling after TFI. A smoother block-edge geometry (block-edge
     smoothing before TFI) would avoid it.
- **Generic:**
  - Walls, mantle, chains and their ends (v=1 → end on the boundary, v=3 → strip along the straight separatrix)
    come from geometry and valences.
  - The valence-6 tip vertices are a fixed consequence of pillowing a zero-thickness slit.
  - Fixed values: δ = 0.03, N_L = 12, h1 = 1e-3·chord, normal smoothing ±δ/2, nose band t/2, untangle thresholds
    0.1 (2D) / 0.05 (3D), wall section offset 0.25° at the symmetry planes.

### 10.6 Files (ParaView / pictures)
- `sweep_tfi_blocks.vtk`: final hex mesh. Cell data:
  - `block`, `block2d`;
  - `kind` (0 outside, 1 core, 2 half, 3 C/O-layer, 4 strip);
  - `sJ`.
- `sweep_noTFI_blocks.vtk`: stacked sections, before TFI.
- `sweep_block_edges.vtk`: edges of all 100 blocks.
- `layout2d_c_phi{00,45,90}.vtk`: 2D section meshes with `block` and `sJ`.
- pictures:
  - `cgrid_sections.png`: φ = 0/45/90; full section, lip TE / ring-2 LE, ring-3 TE, mantle end + strip; walls in
    blue;
  - `sec2d_phi*.png`, `blocks3d_c.png`, `walls_zoom.png`.
- metrics: `sweep_metrics.json`, `hausdorff_c.json`, `sweep3d_c.log`, `wall_faces.log`, `walls.json`,
  `walls/wall_*.json`.
- scripts: `wall_faces.py`, `cgrid2d.py`, `sweep3d_c.py`, `hausdorff_c_occ.py`, `hausdorff_c_mesh.py`,
  `plot_summary.py`, `plot_c.py`, `test_c.py`.
- Change to the §9 code: `layout2d.py`'s `Mesh2D` got two optional hooks (`extra_cuts`, `extra_b`). These are
  backward-compatible; the §9 results are unaffected.

---

## 12. Block collapse on the C-grid variant (2026-11-14)

Directory: `/work/canada/blocking_study/sweep_phi/block_collapse/cgrid/`. Input: the 100-block C-grid structure of
§10 (`cgrid/sweep_tfi_blocks.vtk`, not modified). Code: §11, unchanged except for the additions below.
Tags: **[V]** measured, **[J]** judgement.

### 12.1 Setup [V]
- **CAD.**
  - `cad_samples_multi.py` samples all 20 faces of `case/cad_blk.stp` plus the 4 symmetry/inlet/outlet planes of
    the baffle domain. The selection rule is the same as in `cgrid/hausdorff_c_occ.py`; samples with x<0 or y<0
    are dropped.
  - `cad_blk.stp` consists of **unconnected faces** (no shared edges), so its edge/vertex topology cannot be used.
  - New option `--cad-topology labels`: CAD edges = connected chains of block edges whose boundary block faces
    carry the same set of ≥ 2 CAD-face labels; vertices = chain ends, junctions, and nodes with ≥ 3 labels.
    Result: 48 chains, 32 vertices.
  - Face labels come from the distance to the local tangent plane of each candidate face (PCA of the 8 nearest
    samples). Plain nearest-sample voting mislabelled thin layer faces on the symmetry planes as wall (majority
    0.50); the new rule gives min majority 0.94.
- **Wall layer (new, `--wall-faces 0..12` = the walls.json components).**
  - Topology: a sheet that is wall-normal at a wall face (transverse face on a wall) is rejected ("wall-layer").
  - Geometry: the first wall-normal cell must stay within ±1 % of the input [2.3297e-3, 2.3297e-3] at all 6118
    wall nodes.
- **Guard.**
  - **min sJ ≥ input min sJ = 0.0353** (`--sj-floor input`). The absolute 0.05 would reject everything, because
    the knife-TE layer of the input is already at 0.035.
  - 0 inverted; h_fwd ≤ baseline + 0.01.
  - CAD→mesh is checked **per sample relative to the baseline**, at most +0.01 (`--rev-relative`). Plane samples
    inside the real walls are never reached by the mesh (absolute max 0.074 in the input).
- **Fix in `bstruct.TriSurface`.** Candidate triangles are now also found via the nearest vertices. With centroids
  only, the long BL triangles were missed and projection slid nodes on the symmetry plane: 144 inverted in the
  baseline. After the fix, the **baseline re-fit equals the input**: 0 inverted, min sJ 0.0353, wall first cell
  2.3297e-3. The §11 baseline already reproduced its input, and its results stand [J].

### 12.2 Sheet inventory [V] (25 sheets, Σ = 300 = 3·100)
- **2 φ layers**, 52 blocks each (30 outside + 12 layer + 2 strip + 4 halves + 4 cores).
- **3 wall-normal layer sheets:**
  - 14 blocks: lip C-layer + strip;
  - 8 blocks: ring-2 O-layer;
  - 6 blocks: ring-3 O-layer.
- **18 swept 2D chords**, 4 of them through the butterfly. 11 of them now also run through layer blocks (cut
  tangentially). This includes the §11 strip chord: 2D 15 now continues through layer 2D 34/43.
- **2 butterfly radial sheets**, halves only (6 and 2 blocks).

### 12.3 Results [V]
**Admissible collapses: none (25 sheets × mid/lo/hi). Final: 100 blocks (unchanged).** Quality as input: 0 inverted,
min sJ 0.0353, wall first cell 2.330e-3.

Rejections (first reason; geometry = best placement):

**Wall layer: 3.** The three wall-normal layer sheets.

**Topology: 8.**
- 3× two CAD (label) vertices would be welded;
- 2× forced boundary side incompatible with the welded node class;
- 2× two blocks would get the same corners;
- 1× region would vanish.

**Geometry: 14.**

| sheet | inverted hexes | min sJ | other |
|---|---|---|---|
| 2 φ layers | 96 / 265 | ≈ 0 | degenerate axis, as expected |
| 10 chords between or around walls | 419 – 109 085 | –0.69 to –1.0 | first wall cell destroyed (down to 1e-4…6e-7 and up to 0.37); CAD→mesh +0.015…+0.07 |
| butterfly radial, 6 halves | 0 | **0.029** < 0.0353 | wall first cell 2.170e-3 (−7 %); CAD→mesh +0.023; h_fwd 0.0190 > 0.0188 |
| butterfly radial, 2 halves | 0 | **0.035** (passes) | rejected **only** on wall first cell (−7 %), CAD→mesh +0.023 and h_fwd 0.0190 |

In the 10 chord collapses, the chord sits between two pillowed walls. Collapsing it squeezes the wall layers onto
each other.

[J] The −7 % on the walls after the butterfly collapse is not a wall-layer edit. It is most likely the per-sheet
interval count (sheets merge and take the max count). Not verified further.

### 12.4 Judgement [J]
- The C-grid structure is **collapse-minimal under the BL constraints**. Every sheet either:
  - is a wall layer; or
  - would weld distinct CAD vertices or edges; or
  - would squeeze two walls' layers together; or
  - (butterfly) degrades the axis region and changes the wall grading via the shared interval counts.
- The one §11 gain (lip-tip/ring-2-LE strip) does not carry over. The ring-2 O-layer now runs through that chord,
  so the chord can no longer be collapsed in isolation.
- Compared with §11, the label-derived CAD topology is weaker: no CAD vertices inside a label chain. More sheets
  therefore reach the geometry stage, and the geometric guard still rejects them. Pillowing also removed the
  "wall-to-wall in one block" topological rejections. Those cases now show up as massive inversion instead.
  [J] A cheap extra topology rule would catch them earlier: reject a sheet whose two sides belong to the layers
  of different walls.
- The **block-level collapse is useful here as a certificate**: it shows that 100 is a local minimum for this
  layout under the wall-layer rules. It gives no reduction.

### 12.5 Files
- `cgrid/run/block_collapse.log`, `cgrid/run/records.json`: every sheet/variant with its reason and metrics.
- `cgrid/run/baseline_h0.0_{tfi,block_edges,coarse}.vtk`: re-fit of the input. Identical to the input; coarse
  mesh with point data `node_dim` / `node_cad`.
- `cgrid/run.out`, `cgrid/cad_samples.npz`.
- Scripts: `cad_samples_multi.py`, `block_collapse.py` (new options `--cad-topology labels --wall-faces --h1-tol
  --sj-floor input --rev-relative`), `bstruct.py`.
- Command:

      block_collapse.py --vtk ../cgrid/sweep_tfi_blocks.vtk --cad cgrid/cad_samples.npz --cad-topology labels
        --wall-faces 0,...,12 --sj-floor input --rev-relative --out cgrid/run
