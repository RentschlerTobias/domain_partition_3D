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

**Current task:** T10 — throughput. A verified sample exists ON THE CLUSTER.
**Last verified:** 2026-09-11. T0, T1, T2, T3, T4, T5, T6 and T9 `done`.

* T0 — `fixtures/blocks_core.tar.zst` (5.5 MiB) holds what cannot be
  recomputed cheaply; 3.1 GB of provably dead `.hexex` / `final_tet` / log spew
  deleted; the Dockerfile diff is in `external_patches/`.
* T1 — enroot works here, including a **writable** bind mount at uid 0. The
  `enroot-unshare` worry was a non-issue.
* T2 — **done without a rebuild.** `algohex:portable` already existed and
  passes T2's own gate; the coinbrew rebuild on 64 cores was never needed. Its
  recipe is now `external_patches/Dockerfile.portable`. Read T2's status note
  before touching this area.
* T5 — `block_edges.py`. v11: **0.226 % median** residual over 107 edges, 3x
  better than the 2D reference, and 16x better than a straight chord.
* T6 — **`-n` is settled: 2000 and 8000 as the dataset pair.** The fit does
  not collapse at four points; it costs 1.4x. Read T6's warning that the
  residual column runs backwards before quoting any number from it.
* T9 — `export_sample.py`, round-trip PASS on v11 and cand_001. A sample is
  **0.42-0.45 MB** of `.npz` including the labelled surface triangulation.
* T4 — one call site, `--backend docker|enroot`, `--out-dir` fixes the shared
  `.ovm` path, and the useless checkpoint files are opt-in.
* T3 — `algohex.sqsh` (2.11 GB, zstd verified) and the two backends agree:
  **2210 of 2210 cells matched**, each congruent under a cube rotation. Read
  T3 before comparing meshes again — file comparison is NOT a valid gate here.

**The cluster route is done and measured.** Job 6860177 on bwUniCluster
produced a sample with a passing round-trip in 12:25, from an AlgoHex built
there from its pinned public commit — no Docker, no image copied from any
machine. It reproduces the VPS reference at `-n 2000` exactly: 16 blocks, 14
cuboids, 107 edges, fit modes 89/0/18, residual median 0.048 %, Hausdorff
0.04475. That closes the provenance question T2 left open.

Every cluster-side trap is written up in
[`../../docs/cluster-enroot-findings.md`](../../docs/cluster-enroot-findings.md)
— read it before touching the scripts. The two that cost the most: enroot's
`ENROOT_RUNTIME_PATH` defaults into `/run/user/$UID`, which does not exist on a
compute node, and `dev_cpu_il` allows only ~4 QUEUED jobs per user, which makes
dependency chains the wrong shape at any length.

T7, T8 and T12 improve the sample but block nothing; take them whenever the
cluster is queueing. Note that T12 now has a cheap new lead: T9 reports
`interior_faces` and the one-sided shell count on block topology, which is the
same signature T12's failed weld produces. **T13 is new and deliberately
deferred** until T11 runs — the embedding-optimisation idea, with the
measurement that justifies it.

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

**Status:** `done` 2026-09-11 — **the premise below was wrong.** A
self-contained image already existed and passes T2's own gate. What was
actually missing was its recipe, and that is now committed.

**Correction, and how the mistake was made.** This task and decision H both
assert that no transportable AlgoHex exists. Both were reasoned from
`external/algohex-src/Dockerfile` alone, without checking `docker images`.
`algohex:portable` has been on this box since ~2026-09-05, and
`README.md` — in the same directory — already called it "the one to export to a
cluster". Measured 2026-09-11:

```
$ docker run --rm algohex:portable HexMeshing --help      # NO volume mounted
AlgoHex
HexMeshing [OPTIONS] ...                                  exit 0
$ docker run --rm algohex:portable sh -c \
      'ldd "$(command -v HexMeshing)" | grep -c "not found"'
0
$ docker run --rm algohex:portable ls /opt/algohex/Build/bin
HexMeshing  LocalMeshabilityCheck
$ docker run --rm algohex:portable ls /app/build/Build/bin
(nothing -- the volume's path does not exist in the image)
```

`docker history` shows how it was made: `FROM algohex-configured`, a 24.2 MB
`COPY` of the volume's `Build/` tree into `/opt/algohex/Build`,
`LD_LIBRARY_PATH=/opt/algohex/Build/lib:/opt/coin-or/lib`, then `ln -sf` plus a
build-time `HexMeshing -h` self-test.

**So the real gap was never the binary — it was that this recipe existed
nowhere**, having been typed by hand into a gitignored tree. The image was as
unversioned as the volume it came from. It is now
`external_patches/Dockerfile.portable`, with the volume-extraction step, the
build, the verify and the export documented in its header.

**What is still open, and it is not blocking.** The image is a *from-volume*
shortcut: reproducible only while `algohex-build-cache` exists. The from-source
path — uncomment upstream's two final lines, `ninja -j$(nproc)` — remains
unverified and still belongs on a machine with cores. That is now an
improvement, not a prerequisite: T3 can proceed today.

**Original why, kept for the record.** `docker save algohex-configured` yields
an image WITHOUT `HexMeshing`, because the binary lives in the volume
`algohex-build-cache` and volume contents are not part of an image. True of
`algohex-configured`, and the reason `algohex:portable` was built in the first
place.

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

**Do** — reduced to the leftover, for whoever gets a fast machine. Keep
`Release`. Uncomment the two lines and build with `ninja -j$(nproc)`, not
`-j1`. Prefer a 64-core machine or the cluster — there this is minutes, and the
fragile coinbrew chain (MUMPS, IPOPT, Bonmin from source) gets easier with
cores, not harder. Then diff its `HexMeshing` against the one in
`algohex:portable`.

**Done when.** ~~`docker run --rm <image> HexMeshing --help` succeeds with no
volume mounted. This command fails today.~~ **It does not fail; it passes for
`algohex:portable`, shown above.** The lesson for the next task written in this
document: a "Done when" phrased as a command should be RUN before the task is
written, not only after. This one would have cost a full coinbrew rebuild on 64
cores to satisfy something already satisfied.

**Writes.** `external_patches/Dockerfile.portable`. Not
`external/algohex-src/Dockerfile` — that stays as it is, with its diff recorded
by T0 step 4.

---

## T3 — enroot path for AlgoHex, locally

**Status:** `done` 2026-09-11 — **2210 of 2210 cells match**, and the two
runtimes are NOT byte-identical, which turns out to be the more useful finding

**Why.** Everything except the SLURM wrapper can be verified on this machine,
and verifying it here is far cheaper than debugging it in a queue.

What cannot be verified here: SLURM itself (no `sbatch`/`srun`, hence no pyxis,
which is a SLURM plugin), and the Lustre transfer rates that dominate cluster
job time — `eigenfrequencies/docs/cluster-dtoo-enroot-befund-v2.md` measures
5-10 MB/s there, so `enroot create` from an 8 GB image costs 14-27 minutes per
job.

**Do.** Export **`algohex:portable`** (T2 verified it; 2.72 GB) using
`../../../eigenfrequencies/cluster/export_dtoo_enroot.sh` as the template, then
`enroot import`, then run one AlgoHex job at `-n 2000` through `enroot start`
on `data/T1_9/T1_9_tet_v5.vtk`. The import recipe and its pitfalls are in
`../../../eigenfrequencies/cluster/enroot_dtoo_import.md` — note in particular
the naming trap recorded there: the container name is derived from the `.sqsh`
basename, and a mismatch fails silently and late.

Carry over from T1, all measured here: import with `enroot import -o … 
dockerd://algohex:portable` rather than `docker://` (this box has IPv4
forwarding off), set `ENROOT_DATA_PATH`, `ENROOT_CACHE_PATH` and
`ENROOT_TEMP_PATH` to the `/root/enroot/*` tree that the dtOO work already
uses, and pass `--root --rw --mount` so AlgoHex can write its `.ovm` back to
the host. Inside the image the binary is on `PATH` as plain `HexMeshing` — no
`/app/build/Build/bin/` prefix and no volume, which is the difference from
`run_algohex.py`'s current invocation and the reason T4 exists.

Watch free space: 10.3 GiB, against ~2.7 GB for the `docker save` tar plus a
squashfs of similar size. Delete the tar before `enroot create`.

**Done when.** The enroot run and the Docker run of the same input produce hex
meshes that agree cell for cell. **They do.**

**Result, 2026-09-11.** `scripts/export_algohex_enroot.sh` produced
`algohex.sqsh`, 2.11 GB, compression verified as **zstd** — the lzo trap from
`cluster_env.sh` is the one that imports without complaint and is then
unreadable on every compute node, so the script sets
`ENROOT_SQUASH_OPTIONS` itself instead of trusting the environment. It also
self-tests the image with `HexMeshing --help` before spending the time.

Then the same input, `data/T1_9/T1_9_tet_v5.vtk` at `-n 2000`, through both
backends of T4's seam:

| | cells | verts | AlgoHex's own `time_total` |
|---|---|---|---|
| `--backend docker` | 2210 | 2856 | 928 s |
| `--backend enroot` | 2210 | 2856 | 960 s |

enroot is 3.4 % slower here, which is within the noise of a 2-vCPU box under a
fair-use limit and not a reason to prefer either.

**The comparison, and it needed three attempts to state correctly.**

* Byte comparison: **differs**, first at line 51. The difference is `0` against
  `1.2569e-17`.
* Vertex coordinates: **4 of 2856 differ**, by at most **2.776e-16** — 6.4e-17
  of the bounding-box diagonal.
* Cells: matched one-to-one by vertex set, **2210 of 2210**, none unmatched.
  Corner ORDER differs within 2184 of them, and every one of those is
  deckungsgleich under a proper **cube rotation** (checked against all 24).

So it is the same mesh. What differs is four coordinates at machine epsilon and
the corner-traversal order, which is a writer artifact.

**Do not use file comparison as a gate, and do not trust order-dependent
metrics across runs.** Two of this task's own intermediate numbers were wrong
for exactly that reason: a per-cell scaled-Jacobian diff of 1.097 and a total
volume difference of 3.5e-4 were both artifacts of comparing cell *i* of one
mesh against cell *i* of the other after the ordering had changed. Order-free
statements — min scaled Jacobian identical to 12 digits (−0.173380172954), 7
inverted cells in both, identical parametrisation energy to all 16 digits
(52.816477628172926) — agreed the whole time.

The likely cause of the ordering difference is thread count, not the runtime:
the docker run was capped with `--cpus 1.8` and the enroot run with
`/root/bin/capped --quota 90%`, and a different thread count changes reduction
and write order. Same binary, same libraries, same image in both cases. The
practical consequence for T11 is that **AlgoHex output is not bit-reproducible**,
so reproducibility has to be asserted on topology and on geometry-to-tolerance,
which is what `export_sample.py --check` already does.

**Writes.** `experimentell/hex3d_algohex/scripts/export_algohex_enroot.sh` and
an import guide beside it.

---

## T4 — Container backend seam in `run_algohex.py`

**Status:** `done` 2026-09-11 — one call site, two backends, and two bugs from
T0/T2 closed with it

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
by the flag, with identical output. **That is T3's gate and it is what
verifies this task** — the code is in place, the equality is measured there.

**Result, 2026-09-11.** `container_cmd(backend, …)` builds the whole command
line for either runtime; everything above it is shared. `--backend docker`
gives

```
docker run --rm --network=host [--cpus N] -v <repo>:/work algohex:portable HexMeshing …
```

and `--backend enroot` gives `enroot start --root --rw --mount <repo>:/work
<image>.sqsh HexMeshing …`, both verified to build correctly, the enroot form
from the invocation T1 measured.

Four things changed beyond the seam itself:

* **The `BUILD_VOLUME` mount is gone**, and with it the reason this script
  could not leave the VPS. The binary is `HexMeshing` on `PATH` now, not
  `/app/build/Build/bin/HexMeshing`, and the default image is
  `algohex:portable`, not `algohex-configured`.
* **`--out-dir`**, so `gen_hex_<name>.ovm` lands in the sample directory.
  `generate_dataset.sh` passes it and reads the `.ovm` from there.
* **`--cpus` with `--backend enroot` is refused, not ignored.** enroot is not
  a resource manager; on the cluster that is `--cpus-per-task`, on a bare box
  `/root/bin/capped`. A silently unlimited AlgoHex is how this box got
  throttled, so the failure is loud.
* **The checkpoint files are opt-in** (`--checkpoints`, default off). They were
  written unconditionally to make `-n` sweeps cheap by replaying them, which
  T0 established does not work and was already measured not to work — AlgoHex
  writes the tet mesh as binary OVM its own reader rejects. At ~56 MB per run
  that is 560 GB across 10 000 samples for nothing. Still available for the
  singular-graph figure.

`enroot start` runs straight from the `.sqsh`, so there is no `enroot create`
and no rootfs unpack — which is what makes it affordable on Lustre, where
`enroot create` from an 8 GB image costs 14-27 min per job. It also means the
naming trap from `enroot_dtoo_import.md` does not apply: no container name is
derived when the image path is given.

**Writes.** `experimentell/hex3d_algohex/run_algohex.py`,
`experimentell/hex3d_algohex/scripts/generate_dataset.sh`.

---

## T5 — Block-edge extraction and cubic fit

**Status:** `done` 2026-09-11 — v11: **0.226 % median** residual, 3x better
than the 2D reference, and the cubic beats a straight chord by 16x

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

**Result, 2026-09-11.** v11, 16 blocks, 107 distinct block edges (192 lattice
axis lines before dedup):

| | median | p95 | max |
|---|---|---|---|
| **residual / chord** | **0.226 %** | 2.33 % | 3.24 % |
| a straight chord instead | 3.74 % | 17.10 % | — |
| planarity / chord | 0.14 % | 2.24 % | 4.00 % |

Three things follow, and the third is a correction to this task's own premise:

1. **The cubic is worth it, by 16x at the median** (0.226 % against 3.74 %).
   That is the same order as the 5-8x already measured for the generative case
   and it settles the representation: two control points per edge, not a chord.
2. **3D is easier than 2D here, not harder** — 0.226 % against the 2D
   reference's 0.74 %. The edges are also nearly planar (0.14 % median), so a
   single cubic, which is planar iff its control points are coplanar, is not
   fighting the data.
3. **Inflections: 4, on 4 of 107 edges, one each.** A cubic Bezier has exactly
   one inflection of capacity, so one segment per edge suffices. But the 2D
   study's *method* does not transfer, and that took three attempts to
   establish — see `curvature_sign_changes` in the module, which is kept
   precisely so nobody re-derives it. Curvature differentiates twice, which
   amplifies the fine mesh's own 1 %-of-chord vertex jitter by ~1/h²; measured
   against ground truth, the noise (max 5.995 on the 85 edges that provably do
   not S-curve) EXCEEDS the signal (max 2.740 on the 4 that do), so no
   threshold separates them. The reported test integrates instead: how often
   the polyline crosses its own chord by more than 1 % of it.

**One thing the fit had to add, and it is not a 3D problem.** Two free control
vectors need two interior points. The 2D code, lifted verbatim, returns both
control points at the ORIGIN when given only the two endpoints — `A` is the
zero matrix — and the residual cannot see it, because the only points being
measured are the endpoints, which lie on any such curve. A perfect score for a
curve through the origin. It never bit in 2D because streamlines carry many
samples; a block edge at n=2000 carries about four points, so it is the normal
case here. `fit_cubic_bezier` therefore picks its model by available data and
reports which: `cubic` (>= 2 interior points), `quadratic` (exactly 1,
degree-elevated, exactly determined), `chord` (none, honestly labelled). On
v11: 89 / 10 / 8.

**Cross-checked against the existing implementation.** `clean_blocks.
block_edge_curves` already derives the same skeleton from the face partition
and writes it as `<name>_edges.vtk`. It finds 108 curves against the lattice
route's 107, and **101 match exactly** on point count and chord length. All 7
differences trace to a single **collapsed block edge** — mesh vertices 17 and
384, 7.6e-5 apart against a median chord of 0.7195, a factor of 10⁴. The
lattice route keeps it as an edge because the lattice says the corners are
distinct; the topological route absorbs it into its neighbour, because in the
edge graph alone that vertex is a pass-through. The lattice answer is the one
T9 wants, and the degeneracy itself belongs in the exported `quality` fields,
not silently in a fitted cubic.

**Writes.** `experimentell/hex3d_algohex/block_edges.py`.

---

## T6 — The measurement that fixes `-n`

**Status:** `done` 2026-09-11 — the fit does NOT collapse at low n. **Decision:
n = 2000 and 8000 as the dataset pair; n=60000 only as a reference sample.**

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

**Result, 2026-09-11.**

| | edges | pts/edge | cubic / quad / straight | residual median | p95 | max | inflect |
|---|---|---|---|---|---|---|---|
| n=2000 | 107 | 6 | 89 / 0 / **18** | 0.048 % | 1.48 % | 3.02 % | 2 |
| n=8000 | 128 | 9 | 98 / 0 / **30** | 0.117 % | 1.62 % | 7.45 % | 1 |
| n=60000 (v11) | 107 | 17 | 89 / 10 / 8 | 0.226 % | 2.33 % | 3.24 % | 4 |

**Read that table with care — the residual column is upside down.** It gets
BETTER as n gets smaller, which cannot be true of the geometry: fewer polyline
points mean fewer degrees of freedom left to disagree with, and a 2-point edge
falls back to a straight chord whose residual against its own 2 points is
exactly zero. Both effects flatter low n. That is why 18 and 30 straight
fallbacks matter more than the residual: at n=2000 **17 %** of edges, at n=8000
**23 %**, carry no curvature information at all, against 7 % at n=60000.

**The control that answers the question the task actually asks.** Fit v11's
edges from a SUBSAMPLE of their own polylines, then measure against all of the
points — same geometry, fewer samples, so the degrees-of-freedom artifact
cannot hide in it:

| fitted from | 4 pts | 6 pts | 9 pts | 13 pts | all 17 |
|---|---|---|---|---|---|
| residual median | 0.325 % | 0.261 % | 0.250 % | 0.270 % | 0.226 % |
| p95 | 3.45 % | 2.57 % | 2.31 % | 2.29 % | 2.33 % |

**Four points cost 1.4x in median residual and nothing beyond that.** 0.325 %
is still less than half the 2D reference of 0.74 %. The fit does not collapse,
so the risk this task was created to test does not materialise and n=60000 is
not the floor.

**Decision: the dataset runs n=2000 and n=8000, both emitted as their own
sample.** 27 minutes for two structures against 55 for one, and they are
genuinely different decompositions — 16 blocks / 107 edges against 19 / 128 —
so they are augmentation, not duplication (decision H2). n=60000 stays as a
reference and showcase structure, where the 2x cost buys the cleanest edges.

Two caveats to carry into T11. The short-edge problem is real even if the
residual hides it: a low-n structure has more block edges that are only one
cell long, and those are straight by necessity, not by measurement. And
`edge_polyline` is stored per decision D2 precisely so this decision can be
revisited without re-running AlgoHex — if the transformer turns out to need
richer edges, re-fitting is free and re-meshing is not.

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

**Status:** `done` 2026-09-11 — round-trip PASS on v11 and cand_001,
**0.42 / 0.45 MB** per sample including the labelled surface

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

**Result, 2026-09-11.** `export_sample.py … --check` PASSes on both:

| | v11 | cand_001 |
|---|---|---|
| corner vertices / blocks | 50 / 16 | 62 / 21 |
| quad shell / interior faces | 50 / 23 | 64 / 31 |
| directed edges | 214 = 2x107 | 274 = 2x137 |
| fit modes cubic/quad/chord | 89 / 10 / 8 | 117 / 14 / 6 |
| surface | 18 548 tris, 7 labels | 19 048 tris, 7 labels |
| **file size** | **0.42 MB** | **0.45 MB** |

The gate checks four things, and the last two are the ones that would have
caught a silent error: blocks as sets of corner COORDINATES (comparing ids
would only test the remap against itself), the undirected edge-endpoint set,
every stored control pair still reproducing its own stored polyline to the
residual the fit reported, and each directed pair being the reverse of its
twin rather than an independent fit.

Format is `.npz`, not the 2D route's `.pt`: decision C keeps this repository
free of torch, and npz needs numpy alone on both banks. The ragged
`edge_polyline` rides as one flat array plus an offset vector; `params`,
`quality` and `provenance` as JSON strings in the same archive, so a sample
stays one file. 0.45 MB is under decision J's ~1 MB estimate even with the
whole labelled surface triangulation in it.

**Two deviations from decision C's field list, both deliberate.**

1. **The surface goes in as a triangulation, not as `surface_points [N,3]` +
   `surface_label [N]`.** One label per point cannot be written honestly: a
   vertex on a ring between two surfaces belongs to both, and picking a winner
   would bake a tie-break into the EXPENSIVE side of the seam — the one thing
   decision C exists to prevent. Stored as `surface_points`, `surface_tris`,
   `surface_tri_label`, from which any point cloud at any density with any
   tie-break follows in seconds at load time, which is where item E was
   deferred to anyway.
2. **`dir_class_count` is optional and only written with `--target-h`.** The
   counts are a function of a cell size nobody has committed to; the class
   *ids* are structural and always written.

**A 3D-only finding, and it is now a stored quality field.** cand_001's worst
edge misses by **32.7 %** — and it is not the 2D story of a degenerate
mini-edge. It is a well-sampled 26-point edge at a third of the median chord,
whose **arc is 4.5x its chord**: an edge wound around the passage, which one
cubic cannot represent at any residual. That is a geometric limit, not a fit
failure. On cand_001 the four edges above 5 % residual are exactly the four
with the highest arc/chord (4.51, 3.06, 1.66, 1.55) and the lowest planarity,
against a median of 1.003. Chord length does not predict it; arc/chord does.
v11 has none: max 1.210, zero edges over 5 %.

So `quality` now carries `arc_over_chord_p95`, `arc_over_chord_max` and
`edges_over_5pct` — per decision G these are stored, not gated, and they are
exactly what the load-time filter T11 has to build needs. If those edges ever
matter, `edge_polyline` is stored and splitting them into two segments is a
load-time change, not a re-run.

**Writes.** `experimentell/hex3d_algohex/export_sample.py`.

---

## T10 — SLURM batch driver

**Status:** `doing` — written as ONE job with `xargs -P`, not an array. The
shape follows from two measurements, not from taste. Untested on the cluster.

**Why.** 64 cores, not this 2-vCPU box. AlgoHex peaks around 6 GB and must not
run concurrently beyond what memory allows — two concurrent runs are what
killed v7 and v8.

**Do.** Rework `scripts/generate_dataset.sh` into a SLURM array job: one sample
per task, enroot backend, `.ovm` written into the sample directory.

Retention per decision J: keep the sample file plus `<name>_blocks.vtk` and
`<name>_tet.vtk` (~8 MB); delete the refill, the `.msh` exports, `_nfaces`,
`_quality` and `full.*` — 79 of 90 MB — except for a showcase set of 5-10
samples. Bound concurrency by memory, not by core count.

**Done when.** `sbatch --test-only` passes, then a first real batch reports a
throughput number. `sinfo` is denied on this cluster — use
`scontrol show partition`, and see
[`../../docs/cluster-enroot-findings.md`](../../docs/cluster-enroot-findings.md)
for that and ten other traps.

**Not an array, and that is a measured decision.** Two numbers from
2026-09-11 rule the array out:

* `dev_cpu_il` allows about **4 QUEUED jobs per user**
  (`QOSMaxSubmitJobPerUserLimit`), shared with everything else running. An
  array of 30 tasks does not fit a budget of 4.
* AlgoHex ran at **18 % efficiency on 8 cores** — roughly 1.5 cores actually
  in use — and **1.87 GB** at `-n 2000` (job 6860177). One sample per node
  would leave 62 of 64 cores idle.

So: one job, one queue slot, `xargs -P` inside it.
`scripts/batch_samples.slurm` computes concurrency from the slot it was given
and reports which constraint binds:

```
concurrency 32   (memory allows 51, cores allow 32)
```

with `--mem=128G` and `--cpus-per-task=64`. Memory is the constraint that
matters at large `-n` — two concurrent AlgoHex runs at ~6 GB each killed v7
and v8 on a 7.7 GiB box — so `MEM_PER_SAMPLE_MB` must be raised for `-n 60000`
or the node OOMs.

**Resumability instead of walltime guessing.** 30 minutes is the cap and one
wave at `-n 2000` takes 12-13, so one wave fits and two do not.
`scripts/sample_one.sh` exits 0 in milliseconds if `sample.npz` already
exists, so a job killed mid-wave loses only what was still running: resubmit
and the survivors are skipped. The batch summary lists exactly which
directories are missing a sample and which log to read.

**The first test needs no new geometry.** With `CANDS` empty the manifest is
the `-n` sweep over the committed fixture `T1_9_tet_v5.vtk`, and per T6 those
are genuinely different structures rather than duplicates — real samples, not
a rehearsal:

```bash
sbatch --test-only experimentell/hex3d_algohex/scripts/batch_samples.slurm
sbatch            experimentell/hex3d_algohex/scripts/batch_samples.slurm
```

Then the real thing, once geometry is available:

```bash
NS="2000 8000" CANDS="data/random_tistos/investigated/cand_0*" \
    sbatch experimentell/hex3d_algohex/scripts/batch_samples.slurm
```

**Retention** per decision J is in `sample_one.sh`: the sample plus
`blocks.vtk` survive, the `.ovm` and the `_nfaces`/`_quality`/`.msh`
derivatives are deleted. `KEEP_OVM=1` for a showcase set.

**Writes.** `experimentell/hex3d_algohex/scripts/batch_samples.slurm` and
`scripts/sample_one.sh`.

---

## T14 — Geometry supply: dtOO via enroot, driven by eigenfrequencies' adapter

**Status:** `doing` — the sampler exists:
`scripts/generate_machines.py` (decisions Q1-Q4 in
`docs/decisions/2026-09-11-dataset-sampler-strategy.md`). Sobol, seeded,
incremental; state in the dataset root's `sampler.json`; layout
`data/dataset/<strategy>/machine_00NN/` with `params.json` (tracked) and
later `mesh.msh` (gitignored, 429 MB scale). `--preview` draws and writes
params only, verified locally against the live `design_bounds()` (64 draws,
per-dimension realized span 0.97-0.99 of bounds); `--export` calls
`adapter.export_mesh()` through the single eigenfrequencies import seam and
runs only inside the dtOO container. Batch size N is still open: it follows
from T10's throughput number and the measured dtOO yield rate, which the
first `--export` batch measures via `export_error.txt` per failed draw.
Still `todo`: the export batch on the cluster, then `batch_samples.slurm`
with `CANDS="data/dataset/sobol/machine_00*"`.

**Why.** 21 candidate geometries is not a dataset, and the 429 MB of
`mesh.msh` are gitignored, so a fresh clone has none. Decision B left the
sample count open deliberately: it follows from throughput, which T10 measures.

**Everything needed is already built, in the sibling repo.** Located
2026-09-11:

| piece | where | state |
|---|---|---|
| the 30 design parameters **with bounds** | `eigenfrequencies/adapters/machines/tistos.yaml` | **all 30 of our `params.json` keys**, checked key by key |
| parameter vector → `.msh` | `DtooAdapter.export_mesh({label: value})` | one call, tested |
| bounds for a sampler | `DtooAdapter.design_bounds()` → `{label: (min, max)}` | ditto |
| dtOO on the cluster | `dtOO.sqsh` + `cluster/enroot_dtoo_import.md` + `submit_dtoo_enroot_smoke.sh` | enroot, same pattern as fenicsx |
| the tistos case itself | `github.com/ihs-ustutt/dtOO`, `demo/tistos` | see `DATA_GENERATION.md` |

So the chain is:

```
design_bounds()  ->  sampler  ->  export_mesh()  ->  .msh  ->  tet_prep_v5  ->  batch_samples
     exists         MISSING        exists          exists      exists          exists (T10)
```

**dtOO runs IN the container, not beside it.** `run_dtoo_export` imports
`dtOOPythonSWIG` directly — its own docstring says "dtOO container only" — so
the whole Python process must run inside enroot, exactly as
`_run_fenicsx` does it. `env_notes.md` describes a native `~/pe` install; that
is the older state, and the enroot route is the one with an import guide, an
export script and a smoke test beside it.

Two traps are already recorded in `submit_dtoo_enroot_smoke.sh` and must be
carried over verbatim:

```bash
enroot start --root "$DTOO_SQSH" bash -c '
    source /usr/lib/openfoam/openfoam2606/etc/bashrc
    source /dtOO-install/bin/env.sh
    python3.13 ...'
```

* **Both** env files, before python3.13, or `libPstream.so` and
  `libTKFeat.so.7.9` are not found.
* **No `set -u` inside.** OpenFOAM's bashrc reads unset variables
  (`WM_PROJECT_SITE`), so `set -u` aborts before dtOO is even imported and the
  image gets blamed for a bug in the wrapper.

**Do.** A sampler script, perhaps 50 lines, that runs inside the dtOO
container: read `design_bounds()`, draw N parameter sets, call `export_mesh()`
for each, write `params.json` beside each `.msh` in our candidate layout.

**Two decisions it contains.**

1. **Sampling.** For an optimisation, DE draws a population along a gradient.
   For a DATASET the goal is coverage of a 30-dimensional box, where uniform
   random clumps — use Latin Hypercube or Sobol (`scipy.stats.qmc`, already in
   the venv).
2. **Yield.** Not every draw is buildable; `DTOO_FAIL_PENALTY` exists for
   exactly that. A hint from our own data: `data/random_tistos/` holds 8 loose
   meshes plus 21 under `investigated/`, which looks like filtering. The yield
   rate is unknown and multiplies straight onto the target count, so measure it
   on the first batch rather than assuming it.

**Done when.** N buildable geometries exist with their `params.json`, the yield
rate is recorded here as a number, and `batch_samples.slurm` runs over them
with `CANDS`.

**The two passes are split across machines:** the Sobol draw needs only
`scipy` (runs locally / on the login node with the repo venv), the mesh export
needs the dtOO container. `scripts/generate_machines.py --preview` draws state
durably; when every machine dir carries `params.json`, its `--export` pass
skips the draw entirely, so the container needs nothing beyond the adapter:

```bash
# local: draw + dimensions coverage proof (done 2026-09-11, count 64)
python experimentell/hex3d_algohex/scripts/generate_machines.py --count 64 --preview

# cluster: mesh.msh per machine, enroot container, idempotent resubmit
sbatch experimentell/hex3d_algohex/scripts/batch_dtoo_export.slurm
```

Inside the container the job runs `scripts/export_all_machines.sh`, one fresh
`python3.13` process per pending machine (dtOO SWIG state must not accumulate;
see the sibling repo's `dtoo_cfd_build.py` docstring) and skips machines whose
`mesh.msh` already exists.

`batch_dtoo_export.slurm` mounts the repo rw and the eigenfrequencies
checkout at `/ef` (the Q4 import seam), sources both env files before
`python3.13` and runs without `set -u` inside, exactly as
`eigenfrequencies/cluster/submit_dtoo_enroot_smoke.sh` prescribes. It ends
with a yield summary (exported / failed via `export_error.txt` / not
attempted); the yield number from the first batch is what unblocks the
dataset-size decision, and the next step after it is
`NS="2000 8000" CANDS="data/dataset/sobol/machine_00*" sbatch
scripts/batch_samples.slurm`.

**Which mesh — the fluid grid, not the runner solid.** `tistos.yaml`'s
`mech_volume: ruWithRounding_mechMesh` is the STRUCTURAL mesh of the runner
solid (quadratic tets; `turbine_runner` builds it into `runner.msh` for the
modal stage). The dataset pipeline needs the FLUID flow-channel grid instead —
`ru_gridGmsh`, the boundedVolume of type `map3dTo3dGmsh` that produced
`T1_9_ru_gridGmsh.msh` and the 21 candidates, whose geometric entity ids
`tet_prep_v5.classify` reads (hub 1, shroud 2, inlet 3, outlet 4, periodic
5/6). Exporting the wrong volume is not loud: the structural mesh parses
fine, but `reduced_boundary` finds no keep elements, so `tet_prep_v5` reports
`0 boundary triangles` and crashes. `export_all_machines.sh` therefore pins
`DTOO_MECH_VOLUME=ru_gridGmsh` (`export.py` documents the env override).
Despite its name that field is the adapter's export slot, not a mesh kind:
the sibling repo's own `naca.yaml` points it at `gridGmsh`, and `tistos.yaml`
ships the structural mesh only because the modal stage is its default
consumer. `scripts/msh_histogram.py` checks an export's element structure
before `tet_prep` reads it. The first export pass (63 machines, 2026-09-12)
held the solid mesh; those files are preserved beside their design as
`machine_00NN/mesh_mech.msh` — kept for later pipeline tests on the runner —
while `mesh.msh` is the grid this chain consumes.

**Why the chain stops at `sample.npz` — no reattach.** The dataset unit is
the sample file (block structure + tet + params + status); it fully serves
the block-structure transformer. `reattach.py` (core + blade O-grid + wall
layers = the CFD-ready `full.vtk`) stays OUT of the batch chain for three
measured reasons:

1. Economics (decision J, measured on cand_001): the `full.vtk` group is
   66 MB per sample and re-derives in seconds from the retained
   `blocks.vtk` + `tet.vtk` pair (7.1 MB, worth 13-55 min of compute) —
   ~80 GB at 10k samples versus ~900 GB with reattach kept.
2. Environment: reattach is the only stage needing torch (via `dp3d.tmesh`),
   and torch is measured blocked in the cluster environment.
3. Contracts: the CFD-ready consumer does not exist inside this repo.

**When the pipeline joins the eigenfrequencies optimization loop (planned,
see decision Q4 in `docs/decisions/2026-09-11-dataset-sampler-strategy.md`),
reattach runs THERE** — in the eigenfrequencies CFD consumer seam, invoked
per evaluated individual (torch available in that stack), not as a dataset
stage. Until then, showcase samples can get their `full.vtk` rebuilt on the
local box / VPS by rerunning `reattach.py` from the retained pair.

---

## T11 — Full run, and the throughput that settles the dataset size

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

## T13 — Optimise the embedding of the block faces

**Status:** `deferred` — **do not start before T11 runs on the cluster.** The
user's call, and the right one: this improves samples, it does not produce
them, and its value is only assessable once there are enough samples to see
how often the defects below actually occur.

**Why.** Every geometric defect measured in T5 and T9 is an **embedding**
defect, not a topology defect. Our block corners and edges sit wherever
AlgoHex's integer-grid map put them, and no stage ever improves that placement
— there is no such stage in the pipeline at all. Measured 2026-09-11:

| | cand_001 | v11 |
|---|---|---|
| block edges | 137 | 107 |
| **on the quad shell** | **128 (93 %)** | **100 (93 %)** |
| interior | 9 | 7 |
| edges over 5 % fit residual | 4 — **all on the shell** | 0 |
| worst: arc/chord 4.506, chord 0.2067 | 32.67 % residual, on the shell | — |
| degenerate edge, chord 7.58e-05 | — | **on the shell** |

**That table is the argument.** The obvious objection to a surface method —
"it cannot touch our interior block faces" — is answered by the data: 93 % of
block edges are on the shell, and **100 % of the measured defects are.** The
interior is 7 % of edges and has produced no defect worth fixing.

**The method.** Heuschling, Lim & Kobbelt 2026, *Embedding Optimization of
Layouts via Distortion Minimization* (EG 2026 / CGF 45(2), RWTH Aachen) — see
`../../docs/LITERATURE.md` §E. It takes a target surface, a layout
connectivity and an initial embedding, and improves the embedding
**geometrically while preserving connectivity strictly**: repositions layout
nodes, re-embeds arcs as piecewise geodesic curves, inserts extra nodes along
arcs where flexibility is needed. Multi-resolution, so it can be optimised on
a coarse surface and prolongated.

The mapping onto our data is direct. Their target surface is our input
triangulation (`<name>_tet.vtk`, which the exporter already stores as
`surface_points` / `surface_tris` / `surface_tri_label`). Their layout
connectivity is our quad shell (`quad_faces`) with its nodes and arcs. Their
initial embedding is what AlgoHex gave us — the polylines in `edge_polyline`.

**What it would and would not fix.**

* Would: the degenerate shell edge (two layout nodes 7.6e-5 apart, node
  repositioning is exactly that operation), and the winding edges, where a
  geodesic re-embedding shortens an arc that currently wraps 4.5x its chord.
  Both are the direct cause of our worst curve-fit residuals.
* Would not: the block COUNT. Connectivity is preserved by construction, so
  this is no help for T12 or for the fewness question — that stays the
  Gao/Xu/Duan family in `LITERATURE.md` §C.

**Open before starting.**

1. Their code was announced but is not published yet
   (`github.com/7-AlexH/layout-embedding-optimization`). Check first; a
   reimplementation is a different size of task.
2. Whether re-embedding the shell keeps the complex **weldable**. Moving a
   shell node moves the interior block faces that meet it, and the
   INDEX-parameterised resampling in `refill_block` is what currently makes
   blocks weld. This is the real risk and it is the same class of problem as
   T12's failed weld.
3. Whether it belongs before or after the curve fit. `edge_polyline` is stored,
   so re-fitting after re-embedding is free — which argues for after, and for
   treating this as a load-time-adjacent stage rather than a re-run.

**Writes.** Nothing yet.

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
