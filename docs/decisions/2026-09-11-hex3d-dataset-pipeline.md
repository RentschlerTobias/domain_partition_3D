# hex3d dataset pipeline: what a training sample is, and how it is produced

Grilling session, 2026-09-11. Continues
[2026-09-08-hex3d-block-structure-objective.md](2026-09-08-hex3d-block-structure-objective.md),
which settled *that* the block structures are training data. This one settles
*what a sample is* and *how the pipeline produces it*.

Ten decisions, none of them implemented. The executable plan is
[`../../experimentell/hex3d_algohex/DATASET_PIPELINE.md`](../../experimentell/hex3d_algohex/DATASET_PIPELINE.md).

## Design tree

```
Goal: a data-generation pipeline whose samples train a 3D block-structure transformer
│
├── x A   Sample artifact: 3D hex complex, not a 2D slice
├── x I   One structure for both routes, subdivided
├── x D   Curved edges fitted from the fine-mesh polyline
├── x D2  Stored as two absolute cubic control points
├── x C   Thin seam: neutral geometry file, meshtron owns the ML format
├── x B   Sample volume settled empirically after the first cluster run
├── x H   Self-contained image, enroot on bwUniCluster 3.0
├── x H2  Measure the fit first, then use n as an augmentation axis
├── x G   Store quality fields, filter at load time
├── x J   Keep the sample file plus blocks.vtk and tet.vtk
├── -  E   Conditioning point cloud            DEFERRED to the transformer work
├── -  F   Canonicalisation, coordinate frame  DEFERRED to the transformer work
└── -  K   Frame-field data for a later GNN    PARKED, separate plan
```

---

## A. The sample is the 3D hex block complex

**Decision.** A training sample is the 3D block complex produced by
`clean_blocks.py` (`<name>_blocks.vtk`), carrying curved block edges. The
existing 2D Quadtron stack is not modified; porting the transformer to 3D is
separate, later work.

**Options considered.** (A) 3D-native hex complex, or equivalently its quad
shell, (B) 2D quad structures cut out of the 3D complex to augment the
existing 10k 2D set, (C) 2D only, with hex3d reduced to a CFD and thesis
artifact, (D) A plus a 2D slice exporter.

**Why.** B and C make the entire hex3d branch a by-product: the cheap 2D
generator (`domain_partition/data_generator.py`) already covers the 2D
distribution with 10 014 samples, so 21 sliced runners add nothing there.

A does not contradict the phrase "quad block structure". A hex complex is
bounded by quad faces, and predicting that quad shell is the smaller change to
plan B's pointer head than predicting 8-corner cells.

And `curved_refill.py` measured that curved block edges are worth 5-8x in the
generative case — the case a transformer creates. That is a 3D result and can
only be cashed in 3D.

**Consequences.** Fixes a public surface — the sample schema — that the
tokenizer, the pointer head, the reconstruction and the geometric metrics will
all hang off, across two repositories. Hard to reverse.

It also commits the project to the data-volume problem. The 2D model
overfitted at 100 samples (train ppl 2.8 against val 17.8,
`meshtron/docs/ho_quad_transformer/01_current_model_and_diagnosis.md`) and
needed 10 014, while `data/random_tistos/investigated/` holds 21 geometries,
six of which have run end to end.

## I. One block structure serves both routes, via subdivision

**Decision.** One artifact per geometry: the coarse block complex with curved
edges. Training data for the first-order standard route is derived from it by
subdividing each block 1x or 2x.

**Why.** Stated by the user, and the 2D implementation confirms it.
`meshtron/augment_subdivide.py` subdivides each higher-order block by
transfinite Coons interpolation into n x n sub-cells (face count 6 -> 6n^2) on
the SAME domain, and `export_augmented.py --ns 2 3 4` turns 7 988 six-face
meshes into roughly 31 900.

Conformity is what makes it work, and the mechanism is worth stating because
the 3D port depends on it: each block edge is divided evenly by ARC LENGTH of
the true curve, so a shared interior edge — the same curve traversed in
reverse — yields coincident division points, which a global tolerance dedup
then welds.

**Correction this replaces.** Earlier in the same session the assistant
claimed the two routes want different collapse depths, on the grounds that v11
spans 224 to 11 900 cells per block (a factor of 53) and the first-order route
prefers blocks of comparable size, so `--collapse-rounds` would have to be
tuned twice. That is wrong. Subdivision is a post-process on the finished
artifact, costs seconds, and produces the uniformity directly.

**Consequences.** `--collapse-rounds` is tuned for one objective, fewest
blocks. The 3D analogue of `augment_subdivide.py` becomes required rather than
optional, and needs the 3D Gordon-Hall map (`tfi.tfi`) where the 2D version
needs only a Coons patch. The 2D version also carries a C1 edge-smoothing pass
(`_c1_align`) because block-wise Coons leaves tangent kinks at block borders;
the 3D port will meet the same problem.

## D. Curved block edges are fitted from the fine-mesh polyline

**Decision.** The cubic block edges are fitted through the polyline of fine
hex-mesh vertices that already constitutes each block edge in
`<name>_blocks.vtk`. Not from AlgoHex's frame field.

**Options considered.** (A) cubic through the fine-mesh polyline, (B)
integrate tangents out of the frame field, (C) project the polyline onto the
input triangulation first and fit there, (D) C on boundary edges and A on
interior ones.

**Why B fails on three counts.** The frame field is never written out:
`run_algohex.py` saves only the seamless map (`--sm-out-path`, `.hexex`) and
the post-singularity tet mesh (`--final-tetmesh-out-path`, `.ovm`), and
`DATA_GENERATION.md` already measured that the latter is binary OVM which
AlgoHex's own reader rejects. Reaching the field means writing code inside
AlgoHex.

Beyond availability, a frame field is a per-tet DIRECTION, not a curve. Turning
it into a block edge means integrating it along that edge — a streamline
trace, which is what the 2D pipeline does with `edge_to_streamline`.

And it is the field that produced the hex mesh in the first place, so it is
the upstream source, not the better one.

**Why A.** The polyline is already there — roughly 25 points per edge at
n=60000 — and it tracks the geometry to 0.006. `curved_refill._cubic_axis`
fits it in one call, and the result is measured:

| divisions kept | chord p95 | coons p95 | gain |
|---|---|---|---|
| 100 % | 0.08172 | 0.01059 | 7.7x |
| 50 % | 0.08009 | 0.01103 | 7.3x |
| 25 % | 0.08056 | 0.01533 | 5.2x |

0.0106 is close to the 0.006 the full AlgoHex mesh itself achieves, which is
the ceiling for any method that resamples this mesh.

**Consequences.** Cheap to reverse: the fit reads `<name>_blocks.vtk` and runs
in seconds, so changing the source re-runs the exporter, not the 43-minute
pipeline. Option C stays open as a measurement — it needs
`clean_blocks.project_to_surface`, which already exists — and is defined only
for boundary edges, since interior edges have no input geometry to project
onto.

**Open risk, which gates H2.** Fit quality depends on the fine-mesh
resolution: about 25 points per block edge at n=60000, but roughly 4 at
n=2000, which is the bare minimum for a cubic.

## D2. Edges are stored as two absolute cubic control points

**Decision.** Each directed block edge stores the two interior control points
of a cubic Bezier as absolute coordinates, in the same system as the block
corners. Six scalars per edge. Chord-local encoding, quantization and any
other invariance is the tokenizer's business and is applied at load time.

**Options considered.** (1) absolute control points, (2) chord-local `s,h` per
control point as the 2D route chose, (3) Hermite end tangents, (4) the
polyline resampled to K fixed points.

**Why the 2D answer does not transfer.** `06_edge_geometry_study.md` picked
chord-local `cubic_bezier` on two grounds and neither survives the move to 3D.

Its length advantage was against HERMITE, not against absolute control points.
In 2D both cost four scalars; the six-token figure in that table is hermite's
two sincos-encoded angles plus two magnitudes. Counted in 3D:

| representation | scalars per edge |
|---|---|
| absolute control points | 6 |
| chord-local | 6 |
| Hermite | 10 |

Chord-local buys nothing in length.

And chord-local is ambiguous in 3D. The space perpendicular to a chord is a
PLANE, not a direction, so the encoding needs a canonical transverse frame
recoverable from the two endpoints alone. Every simple rule degenerates
somewhere — global z fails on axial edges, the local radial direction fails on
radial ones, an adjacent face normal is ambiguous between two faces. The
degeneracy is not hypothetical here: `DATA_GENERATION.md` measures exactly one
radial direction class threading every block in every sample, so a rule keyed
on the radial or the axial direction meets its bad case immediately.

**Why storing absolute forfeits nothing.** The map between absolute and
chord-local is a pure function of the two endpoints, which the sample already
carries. The tokenizer applies it while loading. That keeps the choice of
frame open until the 3D model exists, without regenerating a single sample.

**Consequences.** Cheap to reverse — a load-time transform, not a re-export.
Option 4 would have been the expensive one, changing the shape and size of the
schema.

**Two measurements owed before this counts as validated.** The fit residual of
the cubic against the true polyline; the 2D median was 0.74 % of chord length.
And an inflection count: the 2D study found zero inflections over 330 164
edges, which is what justifies a single cubic segment per edge. In 3D a space
curve can also twist out of plane; one cubic Bezier represents torsion, but
the number should be seen before committing to one segment.

## C. Thin seam: a neutral geometry file, meshtron owns the ML format

**Decision.** `hex3d_algohex` writes one neutral geometry file per sample.
Everything ML-shaped lives in meshtron and runs at load time. This repository
stays free of torch.

```
vertices        [M,3]     block corners, Cartesian
blocks          [K,8]     global corner indices per hex block
quad_faces      [F,4]     the quad shell
edges           [E,2]     directed block edges
edge_ctrl       [E,2,3]   the two cubic control points
edge_polyline   [E,*,3]   the raw polyline          <- insurance for re-fits
dir_class       [E]       direction-class id and solved cell count
surface_points  [N,3]     surface points from _tet.vtk
surface_label   [N]       which of the seven surfaces
params          dict      the ~50 dtOO design parameters from params.json
quality         dict      block count, cuboids, tiny, inverted, min sJ,
                          Hausdorff, validator verdict, refill outcome,
                          curve-fit residual
provenance      dict      geometry name, -n, collapse rounds, git sha
```

**Options considered.** (1) thin seam as above, (2) a ready-to-train `.pt`
dict written here, mirroring `meshtron/domain_extractor.py`, (3) a third
module owning both.

**Why.** It is the split the 2D route already uses and that worked there:
`domain_partition/` writes raw `checkpoint_mesh_*.pt`, and
`meshtron/domain_extractor.py` produces the ML dict. Two conventions for one
project would be a self-inflicted wound.

More importantly it puts the expensive side and the volatile side on opposite
banks. The geometry file costs 13 to 55 minutes; the ML format costs seconds.
Every quantity that will be tuned during training — point-cloud density,
normalisation, sort order, quantization — belongs on the cheap side, because
tuning it must never re-run AlgoHex.

**Consequences.** The geometry file must carry everything downstream could
need, since anything missing costs a full re-run. `edge_polyline` is the
insurance premium: while the raw polyline is stored, the curves can be
re-fitted by any method, which is what makes D and D2 cheap in retrospect.
Items E and F leave this session and become part of the transformer work.

## B. Sample volume is settled empirically, not now

**Decision.** No target sample count is fixed. Once the pipeline runs it is
tested on bwUniCluster's 30-minute `dev_cpu_il` slot, throughput in a 48-hour
window is measured, and the 3D transformer is then trained with subdivision
augmentation to see how badly it overfits. The number follows the measurement.

**Why.** Stated by the user, and the only honest option: throughput depends on
concurrency under AlgoHex's ~6 GB peak, on `-n`, and on the rejection rate,
none of which are known yet.

**Consequences.** Shifts the session's weight onto H. The pipeline must be
built for a batch scheduler and many concurrent samples from the start — the
target is 64 cores, not this 2-vCPU box.

## H. AlgoHex ships as a self-contained enroot image

**Decision.** AlgoHex is packaged as a self-contained Docker image and
imported on bwUniCluster 3.0 as an enroot squashfs, following the pattern in
`duty/eigenfrequencies/cluster/`. The same image runs under Docker on the
user's own 6-machine cluster, where Docker rights exist.

**Options considered.** (1) enroot, (2) Apptainer `.sif`, (3) both, (4) native
build on the cluster, (5) keep AlgoHex on the VPS and move only the
post-processing.

**Why enroot and not Apptainer.** The assistant initially recommended
Apptainer on the general claim that it is the HPC standard. That was a
generalisation and it is wrong for this account. enroot is proven on the exact
target machine: `cluster/enroot_dtoo_import.md` is titled "dtOO enroot import
guide — bwUniCluster 3.0", the export/import/submit scripts sit beside it, and
shell history from nodes `uc2n602`/`603`/`605` shows `enroot start --root -m
... dtOO.sqsh` against a 5.8 GB workspace image, with `enroot list` reporting
a pyxis container. Choosing Apptainer would mean re-deriving mounts, workspace
paths and SLURM integration that already work.

Option 5 was rejected outright: AlgoHex takes about 10 minutes per sample and
must not run concurrently (~6 GB peak), so leaving it on a throttled 2-vCPU
VPS would leave 64 cluster cores waiting on two.

**The prerequisite, which is the real work.** There is no self-contained
AlgoHex image today. `external/algohex-src/Dockerfile` ends after the `cmake`
configure step — its final two lines, `RUN cd /app/build && ninja -j1` and the
symlink into `/usr/local/bin`, are commented out. The compile was run by hand
into the Docker volume `algohex-build-cache`, presumably because a
non-resumable single-layer `ninja -j1` on a 2-vCPU box would have taken hours.
Consequently `docker save algohex-configured` produces an image WITHOUT the
`HexMeshing` binary, since volume contents are not part of an image. This
holds for any runtime, enroot or Apptainer alike.

Two things follow. The first deliverable of the whole plan is a self-contained
image: re-enable those two lines and build on a fast machine, where `ninja
-j64` is minutes rather than hours — the cluster makes the fragile coinbrew
chain (MUMPS, IPOPT, Bonmin from source) easier, not harder. And this is not
only a portability task: today the entire data-generation capability of the
project lives in one unversioned Docker volume on one VPS.

**What can and cannot be verified locally.** enroot 4.2.0 is installed on this
box (`/usr/local/bin/enroot`), with `squashfuse` and `mksquashfs`, and all its
prerequisites are met: KVM virtual machine rather than a container, uid 0,
`/dev/fuse` present, `unprivileged_userns_clone = 1`,
`max_user_namespaces = 31642`, and both `unshare -U -r` and `unshare -m`
succeed. So `enroot import` and `enroot start` can be exercised here. What
cannot: SLURM (no `sbatch`/`srun`, hence no pyxis, which is a SLURM plugin),
and the Lustre transfer rates that dominate cluster job time — measured at
5-10 MB/s in `eigenfrequencies/docs/cluster-dtoo-enroot-befund-v2.md`, so
`enroot create` from an 8 GB image costs 14-27 minutes per job there.

Note that the earlier claim "enroot cannot be tested here" was true when it
was made: `duty/quadmesh/scripts/dryrun_planb.slurm` records that its scripts
were written "in einer Umgebung ohne Cluster-Zugang (kein sbatch, kein enroot,
kein SSH zu uc3)". enroot has been installed since.

## H2. Measure the curve fit first, then treat n as an augmentation axis

**Decision.** Before fixing `-n`, measure the cubic fit residual and the
inflection count on block edges produced at n=2000, n=8000 and n=60000 — all
three already exist on disk for T1_9. If the fit holds at low n, the pipeline
runs two or three values of n per geometry and emits each as its OWN sample.

**Why `-n` is not an efficiency knob.** Measured on T1_9, same input:

| | cells | blocks | cuboids | classes | pinch | Hausdorff |
|---|---|---|---|---|---|---|
| n=2000 | 1872 | 16 | 14 | 11 | 2 | 0.0448 |
| n=8000 | 7028 | 19 | 18 | 9 | 0 | 0.0312 |
| n=60000 | 54460 | 16 | 14 | 11 | 2 | 0.0341 |

n=2000 reproduces v11's structure exactly from 30x fewer cells, while n=8000
produces a DIFFERENT valid structure. Several valid collapse endpoints exist
and the resolution selects one.

**Why low n is cheaper and yields more.** AlgoHex is nearly independent of n —
12.2 / 10.7 / 10.3 minutes — because the parametrization computes a map and
quantization scales it afterwards. The cost is on our side, where the sheet
collapse iterates over cells: ~1 / ~3 / ~45 minutes.

```
n=60000 alone      55 min  ->  1 structure
n=2000 and n=8000  27 min  ->  2 structures
```

They are not duplicates but distinct valid decompositions of the same domain,
which is the right lesson for a generative model: a geometry does not have
exactly one correct block structure.

**Why the measurement gates it.** A block edge carries about 25 points at
n=60000 and about 4 at n=2000 — the bare minimum for a cubic. The measurement
decides between two worlds: throughput doubles if the fit holds, and n=60000
becomes the floor at 55 minutes a sample if it does not. It also decides what
fits in the 30-minute `dev_cpu_il` slot: no complete sample at n=60000, two at
n=2000.

## G. Quality is stored, not gated

**Decision.** Every run that produces a block structure becomes a sample. Its
quality figures — block count, cuboid count, blocks under 10 cells, inverted
cells, min scaled Jacobian, Hausdorff to the input surface, validator verdict,
refill outcome, curve-fit residual — are stored as fields. Filtering happens
in meshtron at load time. The only hard rejection is a stage that produced no
artifact at all.

**Options considered.** (1) store and filter at load time, (2) hard rejection
at generation time as today, (3) two-tier, hard-rejecting only
`HexBlockValidator` INVALID.

**Why.** Same argument as the thin seam: keep the expensive side and the
volatile side apart. Producing a sample costs 13 to 55 minutes and a rejection
there is irreversible; skipping one at load time costs nothing, and which
criteria actually hurt the model is not yet known. With 21 geometries there is
no room to discard anything.

The current rule turns out to be advisory already. `generate_dataset.sh` runs
`tfi.py --require-lattices` AFTER `clean_blocks.py` has written
`<name>_blocks.vtk`, so the sample itself already exists when the rejection
fires; only the refill is discarded, and `status.json` records `"refill":
"rejected: block without lattice"`. The fields are collected there too, so
this decision is mostly a matter of carrying them into the sample file.

The genuine hard failure looks different: `cand_000` holds only `_tet.vtk` and
a 513-byte AlgoHex log. That is a crashed stage with no structure, not a
quality question.

## J. Keep the sample file, blocks.vtk and tet.vtk

**Decision.** Per sample, keep the neutral sample file (~1 MB), plus
`<name>_blocks.vtk` (4.8 MB) and `<name>_tet.vtk` (2.3 MB). Delete the refill,
the `.msh` exports, the visualisation VTKs and the full CFD mesh after
extraction, except for a showcase set of 5-10 samples.

**Why.** Measured on `cand_001`, a finished sample is ~90 MB, of which:

| | size | recoverable |
|---|---|---|
| `blocks_h0.05.vtk`, `blocks.msh`, `_nfaces`, `_quality`, `_edges` | 21.9 MB | in seconds |
| `full.vtk` / `.msh` / `_part` / `_edges` | 66.0 MB | in seconds |
| **`blocks.vtk`** | **4.8 MB** | **only with 13-55 min of compute** |
| **`tet.vtk`** | **2.3 MB** | **only with 13-55 min of compute** |

`blocks.vtk` and `tet.vtk` are the recovery boundary: while both exist, any
field of the sample file can be recomputed in seconds. Everything else is
derivable from them.

```
sample file only          ~10 GB at 10 000 samples
+ blocks.vtk + tet.vtk    ~80 GB
everything                ~900 GB
```

**Consequences.** Also fix a tidiness bug found while measuring:
`gen_hex_<name>.ovm` currently lands in `output/hex3d_algohex/` rather than
the sample directory, which does not scale to thousands of runs.

## K. Frame-field data for a later GNN (parked)

Storing frame-field quantities alongside each sample, so that a GNN can later
be trained to approximate the field, is wanted but explicitly out of scope
here. It belongs in a separate plan, after this one is executed. Note the
dependency it inherits from decision D: the field is not currently written out
in any readable form, so this item starts with the same AlgoHex-internals
problem that D avoided.
