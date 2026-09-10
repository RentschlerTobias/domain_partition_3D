# Dataset pipeline — the executable plan

Turns the hex3d block structures into training data for a 3D block-structure
transformer. Written so that a session starting with **no context** can pick up
at the current task and continue.

**Read first:** the ten decisions behind this plan, with their measurements and
the alternatives that were rejected, are in
[`../../docs/decisions/2026-09-11-hex3d-dataset-pipeline.md`](../../docs/decisions/2026-09-11-hex3d-dataset-pipeline.md).
Do not re-open a settled decision without reading it. `HANDOFF.md` is the entry
point for the branch as a whole.

---

## Resume here

**Current task:** T2
**Last verified:** 2026-09-11, T0 and T1 `done`. T1 passed on the first
attempt: `enroot import` from `dockerd://`, `create`, and `start` with a
**writable** bind mount, uid 0 inside, host file read and host file written.
The `enroot-unshare` question is settled as a non-issue. Next is T2 — and read
its status note before doing any work, because the premise it was written on
turned out to be wrong.

---

## How to work this document

A task is the unit of work. Nothing outside a task carries state, so an agent
needs to read only the task it is on, plus this section.

Each task has five fixed fields:

* **Status** — `todo` / `doing` / `done` / `blocked`
* **Why** — one sentence, so the task is not executed blindly
* **Do** — the concrete steps, with file paths
* **Done when** — the observable criterion: a command and its expected output
* **Writes** — which files the task creates or changes

Three rules:

1. **Update the status and the "Resume here" block in the SAME commit as the
   work.** A status that lags behind the code is worse than no status.
2. **A task that turns out to be wrong is marked `blocked` with the reason**,
   not silently skipped. The reason is what the next session needs.
3. **Measurements go in as numbers, not as "looks fine".** This branch has been
   fooled at least six times by a metric moving the right way while the
   structure got worse — see "Beware counts that look like progress" in
   `HANDOFF.md`.

**Compute note.** This box is a 2-vCPU VPS with a fair-use limit on SUSTAINED
load; twelve hours of two parallel `clean_blocks` runs got it throttled to
~10 % of its own cores. Check `vmstat 1 3` column 17 for steal before starting
anything long. Use `/root/bin/capped` for CPU-bound work and
`run_algohex.py --cpus N` for the container. Long batches belong on the
cluster.

---

## T0 — Secure what exists

**Status:** `done` 2026-09-11 — steps 1, 2 and 4 done, step 3 dropped
(superseded, see T2)

**Why.** Measured on 2026-09-11. Only one item here is irreplaceable, and it is
not the one that looks alarming:

| | asset | reproducible | cost if lost |
|---|---|---|---|
| 1 | **4 unpushed commits** | **no** | the work is gone |
| 2 | `output/deliverable/` (625 MB) | yes, by re-running | hours of compute; T6 needs the n-sweep structures |
| 3 | compiled `HexMeshing` | yes, from public source | rebuild the coinbrew chain — hours on 2 cores, minutes on 64 |
| 4 | Dockerfile diff | trivially | almost nothing; T2 reverts it anyway |

AlgoHex is open source and `external/algohex-src` is a clone of
`cgg-bern/AlgoHex.git` at commit `3519289`, so the binary costs BUILD TIME to
replace, not information. The Dockerfile diff is four lines that switch
upstream's build off, plus one that sets `Release` instead of
`RelWithDebInfo`; only that last line carries any value and it is recorded in
the decision log.

The full picture:

| asset | where | backed up |
|---|---|---|
| pipeline code | `git@github.com:RentschlerTobias/domain_partition_3D.git` | yes — ~~4 commits unpushed~~, pushed in step 1 |
| AlgoHex source | clone of `cgg-bern/AlgoHex.git` @ `3519289` | yes, upstream |
| Dockerfile modification | ~~working tree only~~ `external_patches/algohex-Dockerfile.patch` | yes, step 4 |
| compiled `HexMeshing` | Docker volume `algohex-build-cache`, 75 MB, **and the image `algohex:portable`** | no — rebuildable from public source; T2 has the recipe |
| v11, v11m, v16m, the n-sweep, the tet inputs | `fixtures/blocks_core.tar.zst`, 5.5 MiB | yes, step 2 |
| the rest of `output/` — 6 gen samples, refills, `*_hex*.ovm` | local only, gitignored | no, 1.8 GB after the step-2 prune |
| `data/random_tistos/` | local only, gitignored | no, 588 MB, regenerable via dtOO |

**Do.**

1. ~~`git push origin data_generation`.~~ **Done 2026-09-11**, `a59417a..1087f26`.
   The four commits that existed only on this disk — among them `eff3fa9
   feat(hex3d): curved block edges are worth 5-8x` and `40e75e1`, the session
   handoff — are on origin, together with this plan.
2. ~~Decide what happens to `output/deliverable/` (625 MB).~~ **Decided and
   executed 2026-09-11.** Two halves, because the 625 MB is not one asset.

   **Committed, compressed.** The subset that is irreplaceable AND small goes
   into git, which is the only off-box store this box has — no rclone, restic
   or borg, no git-lfs, no configured SSH host. Measured: **22.4 MiB of ASCII
   VTK → 5.5 MiB at `zstd -19`**, now
   `experimentell/hex3d_algohex/fixtures/blocks_core.tar.zst`. It holds v11,
   v11m, v16m, the n=2000/n=8000 structures, every `*.divisions.json`, and the
   two tet inputs `T1_9_tet_v5.vtk` (v11's input, and T3's) and
   `T1_9_tet_v11.vtk`. Paths inside the tar are repo-relative, so
   `tar -I zstd -xf …` from the repo root restores them where the code looks.
   This also fixes a live fragility: `tests/test_divisions.py:19` and
   `tests/test_refill.py:25` both hardcode
   `output/hex3d_algohex/deliverable/T1_9_blocks_v11.vtk`, which is gitignored
   — the suite was one `rm -rf output/` from being unrunnable.

   **Deleted, because it is provably dead.** Not "large", *dead* — nothing in
   the repo can read it:

   | group | apparent | why |
   |---|---|---|
   | `*.hexex`, 28 files | 1.6 GB | checkpoint reuse was measured to FAIL — `DATA_GENERATION.md`, "Checkpoint reuse does not work". AlgoHex writes the intermediate tet mesh as binary OVM and its own `-i` reader rejects it. No Python reads `.hexex`; `run_algohex.py:44` only writes it. |
   | `*final_tet*.ovm`, 28 of 29 | 449 MB | same broken reader. `T1_9_final_tet_v11.ovm` KEPT — `ANALYSIS_PLAN.md` step 04 uses it as the singular-graph input. |
   | 34 top-level `*.log` | 1.1 GB | per-face AlgoHex debug spew, 6.6 M lines in the worst one. Truncated in place to head 200 + tail 2000, 8.1 MB total. The run summaries were never in them: they are the 33 untouched `*_hex_metrics*.json`, 132 KB, which is where the `RUNS.md` timings come from. |

   **Kept on disk, not committed:** the 243 MB of `*_hex*.ovm`. That is
   AlgoHex's own output, 10-12 min of compute each, and it is a better recovery
   boundary than decision J's `blocks.vtk` for the ~35 structures here — from a
   `_hex.ovm`, a changed `clean_blocks.py` can be re-run without AlgoHex. J is
   still right at 10 000-sample scale, where 13 MB per sample is 130 GB.

   **Accepted loss:** everything else. Bounded at ~12 min of AlgoHex per
   structure plus the sheet collapse, 1 / 3 / 45 min at n=2000 / 8000 / 60000
   — cluster work in any case.

   One measurement to record because it will mislead the next session: `/` is
   **btrfs with `compress=zstd:1`**, so deleting 3.1 GB of apparent size
   returned only **9.4 → 10.3 GiB** of free space. Highly compressible text was
   never costing what `du` claimed. Do not size future cleanups off `du`.
3. ~~Optional hedge: snapshot the built binary from the volume.~~ **Dropped.**
   T2 found that the self-contained image `algohex:portable` already exists and
   verified it, which supersedes a tarball of the volume.
4. ~~Optional: record the Dockerfile diff.~~ **Done**, as
   `external_patches/algohex-Dockerfile.patch`, verified by
   `git -C external/algohex-src apply --check --reverse`. It carries two
   independent changes and the header says which one to keep: `Release` yes,
   the commented-out `ninja` no.

**Done when.** `git log origin/data_generation..HEAD` is empty, and a decision
about `output/deliverable/` is written here. **Both hold.**

**Writes.** `experimentell/hex3d_algohex/fixtures/{blocks_core.tar.zst,README.md}`,
`external_patches/algohex-Dockerfile.patch`, a pointer in
`experimentell/hex3d_algohex/README.md`, and the deletions above.

---

## T1 — Verify enroot on this machine

**Status:** `done` 2026-09-11 — it works, including a writable bind mount

**Why.** The whole cluster route depends on enroot, and the previous session's
scripts were written in an environment where it was absent —
`duty/quadmesh/scripts/dryrun_planb.slurm` records "kein sbatch, kein enroot,
kein SSH zu uc3". It is present now but has never been exercised here.

Checked on 2026-09-11, all prerequisites are met:

```
enroot 4.2.0                    /usr/local/bin/enroot
squashfuse, mksquashfs          /usr/bin/
machine                         KVM VM, not a container
uid                             0
/dev/fuse                       present
unprivileged_userns_clone       1
max_user_namespaces             31642
unshare -U -r / unshare -m      both succeed
```

~~One open detail: `enroot-unshare` was not found on `PATH`.~~ **Settled: it
is a non-issue.** The binary exists nowhere under `/usr`, `/bin`, `/sbin`,
`/opt` or `/lib`, while `enroot-nsenter`, `-mount`, `-switchroot`,
`-aufs2ovlfs`, `-mksquashovlfs` and `-makeself` are all in `/usr/local/bin` —
and `enroot start` works regardless, so 4.2.0 does not use it. It was a naming
change, not a packaging defect. Do not spend time on it again.

**Do.** `enroot import docker://alpine`, `enroot create`, then `enroot start`
with a bind mount and `cat` a file from the host side.

**Result, 2026-09-11.** Three things, all of them clean:

```
$ export ENROOT_DATA_PATH=/root/enroot/data \
         ENROOT_CACHE_PATH=/root/enroot/cache ENROOT_TEMP_PATH=/root/enroot/tmp
$ enroot import -o /root/enroot/images/alpine.sqsh dockerd://alpine:latest
   -> 8.0 MB squashfs, 520 inodes
$ enroot create /root/enroot/images/alpine.sqsh     # container name: alpine
$ enroot start --mount /tmp:/hostmnt alpine cat /hostmnt/host_probe.txt
hello-from-host
$ enroot start --root --rw --mount /tmp:/hostmnt alpine \
      sh -c 'echo written-inside-container > /hostmnt/enroot_probe_out.txt; id -u'
0                                   # and the file appeared on the host
```

The second run is the one T3 and T4 actually need: uid 0 inside, and a
**writable** bind mount that AlgoHex can drop its `.ovm` into.

Two details worth carrying to T3. Import from `dockerd://`, not `docker://`:
the image is already in the local daemon and `/root/enroot/import.log` records
"IPv4 forwarding is disabled. Networking will not work." And set the three
`ENROOT_*` paths explicitly — the default data path is
`/root/.local/share/enroot`, which is empty, so a bare `enroot list` does not
see the `/root/enroot/data` containers the earlier dtOO work created.

Cleaned up afterwards (`enroot remove alpine`, `.sqsh` deleted); `enroot list`
is back to just `dtOOtest`.

**Done when.** The container prints the host file. **It did.** T2 and T3 can be
verified locally; only SLURM and the Lustre rates cannot.

**Writes.** Nothing in the repo; the result is recorded in this task.

---

## T2 — Self-contained AlgoHex image

**Status:** `todo`

**Why.** `docker save algohex-configured` yields an image WITHOUT
`HexMeshing`, because the binary lives in the volume `algohex-build-cache` and
volume contents are not part of an image. Nothing can reach any cluster until
this is fixed.

**The change is a revert, not an addition.** Upstream's Dockerfile already
builds AlgoHex; the local working-tree diff switched it off:

```diff
-    -D CMAKE_BUILD_TYPE=RelWithDebInfo
+    -D CMAKE_BUILD_TYPE=Release
...
-RUN cd /app/build && ninja
-RUN ln -s /app/build/Build/bin/* /usr/local/bin
+#RUN cd /app/build && ninja -j1
+#RUN ln -s /app/build/Build/bin/* /usr/local/bin
```

The reason was almost certainly this box: a non-resumable single-layer `ninja`
on 2 vCPUs takes hours, so the build was run by hand into a persistent volume
instead.

**Do.** Keep `Release`. Uncomment the two lines and build with
`ninja -j$(nproc)`, not `-j1`. Prefer a 64-core machine or the cluster — there
this is minutes, and the fragile coinbrew chain (MUMPS, IPOPT, Bonmin from
source) gets easier with cores, not harder.

**Done when.** `docker run --rm <image> HexMeshing --help` succeeds with **no
volume mounted**. This command fails today; its passing is the proof.

**Writes.** `external/algohex-src/Dockerfile` and the patch from T0.

---

## T3 — enroot path for AlgoHex, locally

**Status:** `todo`

**Why.** Everything except the SLURM wrapper can be verified on this machine,
and verifying it here is far cheaper than debugging it in a queue.

What cannot be verified here: SLURM itself (no `sbatch`/`srun`, hence no pyxis,
which is a SLURM plugin), and the Lustre transfer rates that dominate cluster
job time — `eigenfrequencies/docs/cluster-dtoo-enroot-befund-v2.md` measures
5-10 MB/s there, so `enroot create` from an 8 GB image costs 14-27 minutes per
job.

**Do.** Export the T2 image using
`duty/eigenfrequencies/cluster/export_dtoo_enroot.sh` as the template, then
`enroot import`, then run one AlgoHex job at `-n 2000` through `enroot start`
on `data/T1_9/T1_9_tet_v5.vtk`. The import recipe and its pitfalls are in
`cluster/enroot_dtoo_import.md` — note in particular the naming trap recorded
there: the container name is derived from the `.sqsh` basename, and a mismatch
fails silently and late.

**Done when.** The enroot run and the Docker run of the same input produce hex
meshes that agree cell for cell.

**Writes.** `experimentell/hex3d_algohex/scripts/export_algohex_enroot.sh` and
an import guide beside it.

---

## T4 — Container backend seam in `run_algohex.py`

**Status:** `todo`

**Why.** Docker on the user's own 6-machine cluster, enroot on bwUniCluster,
one call site for both.

**Do.**

* Drop the `BUILD_VOLUME` mount — obsolete after T2.
* Put the container invocation behind a single function with `docker` and
  `enroot` as alternative backends, selected by a flag. Keep `--cpus`.
* Fix the output path: `gen_hex_<name>.ovm` currently lands in
  `output/hex3d_algohex/` rather than the sample directory, which does not
  scale to thousands of runs.
* `--prefix` already exists and removes the baked-in `T1_9_`; confirm it is
  used consistently.

**Done when.** The T3 job runs under both backends from one command, selected
by the flag, with identical output.

**Writes.** `experimentell/hex3d_algohex/run_algohex.py`.

---

## T5 — Block-edge extraction and cubic fit

**Status:** `todo`

**Why.** A block edge in `<name>_blocks.vtk` is already a polyline of fine-mesh
vertices — roughly 25 points at n=60000 — and it tracks the geometry to 0.006.
Fitting a cubic through it is what makes the sample higher-order. Decision D
rejected the frame field as the source: it is never written out, it is a
per-tet direction rather than a curve, and it is the stage that produced this
mesh in the first place.

**Do.** New `block_edges.py`.

* Extract the edge polylines. The topology is already available: `tfi.lattices`
  returns per-block lattices, `tfi.face_grid` returns face grids, and
  `base_complex.py` holds the purely topological partition. An edge polyline is
  a lattice axis line.
* Fit a cubic Bezier by least squares with fixed endpoints, lifting
  `meshtron/prototype_twostage.py::TwoStageTokenizer._fit_cubic_bezier` from 2D
  to 3D.
* Store two absolute control points per **directed** edge. Directed, not
  undirected: the 2D route found that half-edges keep parallel edges with the
  same endpoints apart, and the same will hold anywhere two block faces meet
  along different curves.

**Done when.** Fit residual as a percentage of chord length, reported as median
/ p95 / max, plus an inflection count. The 2D reference is 0.74 % median over
330 164 edges (`meshtron/docs/ho_quad_transformer/06_edge_geometry_study.md`).
Write the numbers into this task.

**Writes.** `experimentell/hex3d_algohex/block_edges.py`.

---

## T6 — The measurement that fixes `-n`

**Status:** `todo`

**Why.** `-n` is not an efficiency knob; it changes the target data. Measured
on T1_9:

| | cells | blocks | cuboids | classes | pinch |
|---|---|---|---|---|---|
| n=2000 | 1872 | 16 | 14 | 11 | 2 |
| n=8000 | 7028 | 19 | 18 | 9 | 0 |
| n=60000 | 54460 | 16 | 14 | 11 | 2 |

And the cost sits on our side, not AlgoHex's — 12.2 / 10.7 / 10.3 minutes for
AlgoHex against ~1 / ~3 / ~45 for the sheet collapse. So n=2000 plus n=8000
costs 27 minutes for **two** structures against 55 for one at n=60000, and the
two are distinct valid decompositions rather than duplicates.

The open risk: a block edge carries about 25 points at n=60000 and about 4 at
n=2000, the bare minimum for a cubic. If the fit collapses at low n, n=60000
becomes the floor and a sample costs 55 minutes.

**Do.** Run T5 over the three structures already on disk:
`output/hex3d_algohex/deliverable/T1_9_blocks_n2000.vtk`, `_n8000.vtk`, and the
v11 structure at n=60000.

**Done when.** Residual and inflection numbers exist for all three, and the
`-n` set for the dataset is recorded here as a decision with its numbers.

---

## T7 — Curved TFI

**Status:** `todo`

**Why.** `tfi.py` contains no `cubic`, `spline` or `coons`; `refill_block`
interpolates linearly through the fine grid. `curved_refill.py` has the pieces
(`_cubic_axis`, `_coons`) and the measurement — chord p95 0.0817 against coons
0.0106, a factor of 7.7 — but writes no mesh. It is a measurement script.

**Do.** Give `refill_block` a curved mode: boundary face grids resampled with
`_cubic_axis` and `_coons`; the interior keeps the existing 3D Gordon-Hall map
`tfi.tfi`. Note that resampling is parameterised by normalised INDEX, not arc
length — that is what makes the complex weld, and it must not change.

**Done when.** Two gates, both already used in this branch:

* Identity refill — refilling at the counts a block already has returns every
  boundary vertex to ~2.2e-15.
* The `curved_refill.py` table reproduced through the generator rather than the
  measurement script: at 25-100 % of divisions, coons p95 near 0.011 against a
  chord's 0.082.

**Writes.** `experimentell/hex3d_algohex/tfi.py`.

---

## T8 — 3D subdivision augmentation

**Status:** `todo`

**Why.** One artifact serves both transformer routes. The coarse higher-order
complex is the plan-B target; subdividing it 1x or 2x produces the first-order
standard route's data. `meshtron/augment_subdivide.py` does exactly this in 2D
and `export_augmented.py --ns 2 3 4` turns 7 988 meshes into ~31 900.

**Do.** New `augment_subdivide_3d.py`. The invariant that makes it conforming
carries over unchanged and is the thing to get right: divide every block edge
by **arc length of the true curve**, so a shared edge traversed in reverse
yields coincident division points, then dedup globally by tolerance. Faces use
a Coons patch, volumes use `tfi.tfi`.

Expect the tangent-kink problem that `_c1_align` solves in 2D: block-wise Coons
leaves C1 discontinuities at block borders, because the interior grid lines of
two neighbouring blocks meet at a shared sub-vertex with different tangents.

**Done when.** v11 subdivided at n=2 and n=3 gives 128 and 432 blocks,
watertight, zero inverted cells, and boundary error no worse than the
unsubdivided structure. Re-run `mesh_quality.mixed_metrics` on the result.

**Writes.** `experimentell/hex3d_algohex/augment_subdivide_3d.py`.

---

## T9 — The sample exporter

**Status:** `todo`

**Why.** The seam: this repository writes neutral geometry, meshtron owns the
ML format. Everything that will be tuned during training — coordinate frame,
ordering, point-cloud density, normalisation, quantization — belongs on the
cheap side, because tuning it must never re-run AlgoHex.

Fields:

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

`edge_polyline` is the insurance premium: while it is stored, the curves can be
re-fitted by any method without touching AlgoHex.

**Do.** New `export_sample.py`, reading `<name>_blocks.vtk`, `<name>_tet.vtk`,
`params.json` and `status.json`. Surface points and labels via
`clean_blocks.SurfaceLabeller` — never classify a surface by a coordinate
threshold, it cuts across the triangulation. Direction classes via
`tfi.direction_classes` and `tfi.solve_block_divisions`. Store quality as
fields; do not reject.

**Done when.** Round-trip: reconstruct the block complex from the exported file
alone and compare against `<name>_blocks.vtk`. Topology exact, geometry to
floating-point tolerance. The 2D equivalent reported 10 014/10 014 exact.

**Writes.** `experimentell/hex3d_algohex/export_sample.py`.

---

## T10 — SLURM batch driver

**Status:** `todo`

**Why.** 64 cores, not this 2-vCPU box. AlgoHex peaks around 6 GB and must not
run concurrently beyond what memory allows — two concurrent runs are what
killed v7 and v8.

**Do.** Rework `scripts/generate_dataset.sh` into a SLURM array job: one sample
per task, enroot backend, `.ovm` written into the sample directory.

Retention per decision J: keep the sample file plus `<name>_blocks.vtk` and
`<name>_tet.vtk` (~8 MB); delete the refill, the `.msh` exports, `_nfaces`,
`_quality` and `full.*` — 79 of 90 MB — except for a showcase set of 5-10
samples. Bound concurrency by memory, not by core count.

**Done when.** `sbatch --test-only` passes. Note that
`duty/quadmesh/scripts/dryrun_planb.slurm` carries an explicit warning that
queue and module names in this project were written without cluster access and
were never verified — check partitions with `sinfo` first. Then a 3-task array
on `dev_cpu_il` (30-minute walltime), with per-sample runtime, peak memory and
retained bytes matching the estimates (13-27 min, ~6 GB, ~8 MB).

**Writes.** `experimentell/hex3d_algohex/scripts/generate_dataset.slurm`.

---

## T11 — Full run over the 21 candidates

**Status:** `todo`

**Do.** Run the pipeline over `data/random_tistos/investigated/`, then load the
samples from meshtron and confirm the quality fields support the filter that
route wants. Measure achievable throughput in a 48-hour window.

**Done when.** The throughput number exists. That is what settles the dataset
size, which decision B deliberately left open: after this, train the 3D
transformer with subdivision augmentation and see how badly it overfits.

---

## T12 — The v11m question

**Status:** `todo`

**Why.** Whether 12 blocks give a CFD mesh as good as 16 is still open, and the
previous session found the cause had been misdiagnosed twice. v11m's boundary
error of 0.174 comes from a **failed weld**, not from direction-class spread:
counted on the refilled meshes, v11m has 12 482 boundary faces against v11's
9 864 while having FEWER cells, so about 2 618 are interior faces whose two
sides did not weld. Their vertices sit inside the domain, which is why they
measure a tenth of a radius from the surface.

Measured directly at the solved counts, v11m's boundary faces are BETTER than
v11's — p95 0.00493 against 0.00625.

**Do.** `tfi.check_watertight` misses this because it tests for faces used by
three or more cells, while a seam coming apart produces extra faces used by
ONE. Add a boundary-face-count balance, locate the 2 618 faces, then re-run the
12-block refill — now with curved edges from T7, which is exactly the regime
where few large blocks are supposed to pay off.

**Done when.** The 12-block refill is watertight under the strengthened check,
and its CFD metrics stand beside v11's in a table here.

**Writes.** `experimentell/hex3d_algohex/tfi.py`.

---

## Deferred and parked

* **E — conditioning point cloud** and **F — canonicalisation, coordinate
  frame, ordering.** Deferred to the transformer work. Both are load-time
  concerns and the exporter stores enough for either choice.
* **K — frame-field data for a later GNN.** Wanted, separate plan. It inherits
  the problem decision D avoided: the field is not written out in any readable
  form, so it starts inside AlgoHex.
* **Option C on curved edges** — projecting the polyline onto the input
  triangulation before fitting, using the existing
  `clean_blocks.project_to_surface`. Open as a measurement, and defined for
  boundary edges only, since interior edges have no input geometry.
