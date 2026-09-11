# enroot and SLURM on bwUniCluster 3.0 — what actually bit us

Measured 2026-09-10/11 while taking the hex3d pipeline from nothing to a
verified sample on the cluster. Every entry below cost at least one failed job,
so it is written as cause → evidence → fix rather than as advice.

The end state it produced: job **6860177**, COMPLETED, 12:48 wall, 1.87 GB,
one sample with a passing round-trip — and numbers identical to the VPS
reference at `-n 2000` (16 blocks, 14 cuboids, 107 edges, fit residual median
0.048 %, Hausdorff 0.04475). The scripts are in
`experimentell/hex3d_algohex/scripts/`.

Read `../experimentell/hex3d_algohex/DATASET_PIPELINE.md` for what the pipeline
does; this file is only about the cluster underneath it.

---

## 1. `ENROOT_RUNTIME_PATH` — the one that is genuinely non-obvious

**Symptom.** A job ran two `enroot start` calls successfully and the third died
with one line and no shell trace at all:

```
mkdir: cannot create directory '/run/user/985462': Permission denied
```

**Cause.** `ENROOT_RUNTIME_PATH` defaults to `${XDG_RUNTIME_DIR}/enroot`, and
`sbatch --export=ALL` carries `XDG_RUNTIME_DIR=/run/user/$UID` from the login
node — where logind created that directory — into a compute node where it does
not exist and `/run/user` is not user-writable. enroot `mkdir -p`s its runtime
path on **every** invocation, so the failure appears whenever the directory is
missing, which is why two calls worked while it still existed and the third
failed after the owning session ended.

**Fix.** Pin it to a path you own, next to the temp path:

```bash
export ENROOT_RUNTIME_PATH="/tmp/${USER}-enroot-run"
```

**Why it is worth its own section.** Nothing in the error names enroot, the
job dies in under a second, and `bash -eux` prints nothing — so the natural
reading is "the build broke", which sends you into the build instead of the
runtime.

## 2. The four enroot paths, and which filesystem each needs

| variable | put it on | why |
|---|---|---|
| `ENROOT_DATA_PATH` | **workspace** | the unpacked container. Shared across nodes, which is what makes "unpack once" work — see §3 |
| `ENROOT_RUNTIME_PATH` | **node-local** (`/tmp`) | §1 |
| `ENROOT_TEMP_PATH` | **node-local** (`/tmp`) | `enroot import` flattens layers through an overlayfs and a parallel filesystem cannot host one. Inherited from `eigenfrequencies/cluster/cluster_env.sh`, not re-measured here |
| `ENROOT_CACHE_PATH` | workspace | just downloads |

Plus one option that is not a path:

```bash
export ENROOT_SQUASH_OPTIONS="-comp zstd -noD"
```

Without it `mksquashfs` defaults to **lzo**, and an lzo image imports without
complaint and is then unreadable on every compute node. Also inherited from the
sibling repo; our export scripts set it themselves rather than trusting the
environment, and `smoke_algohex.slurm` verifies the result with
`unsquashfs -s | grep Compression`.

## 3. Unpack once, then pass the container NAME

`enroot start` accepts either a `.sqsh` path or the name of a container created
by `enroot create`. Which one you want depends on how many jobs will run:

* **one job** — pass the `.sqsh`. No `enroot create`, no rootfs unpack.
* **many jobs** — `enroot create` ONCE into `ENROOT_DATA_PATH` on the
  workspace and pass the name. Otherwise every task re-reads the squashfs
  through squashfuse across the parallel filesystem.

Measured here: the first pipeline job spent about 4 of its 12 minutes on
`enroot create` (35 955 files, 4 784 directories, 3 085 symlinks) on node
uc2n601. The next job landed on **uc2n605** and did not unpack at all, because
`ENROOT_DATA_PATH` is on the workspace and therefore visible from both. The
unpack is paid once, ever.

The sibling project measured the other direction on a bigger image: unpacking
9.6 GB of `.sqsh` from Lustre took ~20 of 30 budgeted minutes at 20 s of CPU,
pure IO wait (job 6822816). Ours is 2.1 GB, so scale that down — but pay it
once.

**And passing the `.sqsh` path is strictly worse than "slower": it is a
different, broken mount.** `--root <sqsh>` serves the squashfs through
squashfuse, a user-space FUSE process that decompresses on every read. Their
measurement, from `submit_hydroflow_opt.sh`: dtOO's `CreateStates` takes **8
seconds** against an unpacked image and **had not finished after 20 minutes**
under squashfuse, with squashfuse sitting at 37 % CPU the whole time. Copying
the `.sqsh` to node-local disk removed the Lustre latency but not the FUSE
layer. This is what our own job **6865391** was: 30 minutes, 24:24 CPU, zero
output, TIMEOUT. It was not hung, it was crawling — and from the outside the
two are indistinguishable.

So for the dtOO image here the rule is absolute: `enroot create` once into one
shared store (`$(ws_find eigenfreq)/enroot-data`, beside the image both
pipelines start), then `enroot start` the **container name** `dtOO`. The
race-safe install (private `$name.$SLURM_JOB_ID`, then `mv -T` into place;
loser deletes and uses the winner's) is in `batch_dtoo_export.slurm` and copied
from their `submit_hydroflow_opt.sh`, including the mismatch check: a rootfs
that exists but whose name `enroot list` does not show is a broken install,
not a slow one. The name must match the image basename — `dtOO.sqsh` → `dtOO`.

`enroot create` derives the container name from the `.sqsh` basename, so
`algohex.sqsh` → `algohex`. That naming trap only applies to this form.

## 4. enroot can BUILD from source — no Docker needed anywhere

bwUniCluster has no Docker, which is why the cluster route uses enroot at all.
That does not stop you from building a container there: a Dockerfile is a list
of shell commands, and `enroot start --root --rw` runs them.

Verified in miniature on the VPS before trusting it on the cluster:

* `apt-get install` inside `enroot start --root --rw` works. enroot even ships
  `/usr/local/etc/enroot/hooks.d/10-aptfix.sh`, which fires on a Debian rootfs
  — apt in a container is an anticipated case, not a trick.
* Changes **persist** across separate `enroot start` invocations.
* `enroot export -o out.sqsh NAME` captures them (87 → 174 MB in the test).

Then on the cluster, for real: `debian:trixie` imported from Docker Hub, 20 apt
packages, coinbrew fetch, AlgoHex cloned at its pinned commit, Ipopt and Bonmin
compiled, `ninja -j64`, `enroot export`. The resulting image passes the same
gate the Docker one does, and produces identical meshes.

**Consequence for provenance:** the image is a pinned public commit plus two
scripts in this repo. No Docker volume, no registry, no file copied from
anyone's workstation.

## 5. Compute nodes have internet

`eigenfrequencies/cluster` notes say setup must happen on the login node
because "nur dort gibt es Internet". **That is not true on this cluster**, at
least not for `dev_cpu_il`: AlgoHex's cmake downloads its own dependencies
(`ALGOHEX_DOWNLOAD_MISSING_DEPS` is on by default when building standalone,
per upstream's README) and it did so from inside job 6856932, which completed.

This one is recorded because it nearly cost a redesign: the build was being
split around a network boundary that does not exist, adding a manual login-node
step for nothing. Downloading on the login node is still reasonable — it does
not spend a 30-minute slot — but it is convenience, not necessity.

## 6. `sinfo` is denied; `scontrol` is not

```
$ sinfo -o "%P %l %c %m %D"
slurm_load_node: Access/permission denied
```

Use this instead, which works:

```bash
scontrol show partition | grep -E "PartitionName|MaxTime"
```

Result, 2026-09-11:

| partition | MaxTime |
|---|---|
| `dev_cpu`, `dev_cpu_il` | 00:30:00 |
| `cpu`, `cpu_il`, `highmem` | 3-00:00:00 |

And where even `scontrol` is denied, `sbatch --test-only` is the fallback: it
rejects a wrong partition immediately and costs nothing.

## 7. The recommended partition is the wrong one

`eigenfrequencies/docs/cluster-resource-sizing.md` measured `cpu_il` (64 cores,
256 GiB) against `cpu` (96/384) and recommends `cpu_il` on throughput.

In practice `cpu_il` has **about a week of queue time**, which makes it the
slowest route to a first result no matter how many cores it has. Everything
here therefore runs on `dev_cpu_il` and is cut to fit 30-minute slots.
`cpu_il` remains the escape hatch for a job that genuinely cannot be cut.

A throughput recommendation that ignores queue time is a recommendation about
the wrong quantity when you are trying to get the first thing working.

## 8. `QOSMaxSubmitJobPerUserLimit` — about 4 QUEUED jobs

```
sbatch: error: QOSMaxSubmitJobPerUserLimit
sbatch: error: Batch job submission failed: Job violates accounting/QOS policy
```

Measured at **4 queued jobs per user on `dev_cpu_il`**, and the budget is
shared with everything else you are running: with two foreign jobs already
queued there was room for exactly one of four submissions.

**Consequence, and it is a design consequence, not a workaround.** A
`--dependency=afterok` chain needs all of its links queued at the same time, so
it is the wrong shape at any length. `cluster_setup.sh` submits **one job per
run** and is idempotent, so re-running it is the loop:

```bash
bash cluster_setup.sh      # until it says "image exists"
```

Each run picks the next step from the actual state, and every build stage skips
in seconds if its artifact is already there.

## 9. Make long work resumable before you need it

30 minutes is a hard cap, and no estimate for this build existed. Every stage
therefore checks its own artifact (`libipopt.so`, `libbonmin.so`,
`build.ninja`) and exits immediately if present, and a stage killed at the wall
resumes on resubmission because coinbrew drives `make` and AlgoHex drives
`ninja`, both of which keep their object files in the container rootfs.

It paid off: the build took three slots, and the third one skipped Ipopt and
Bonmin in seconds before compiling AlgoHex.

This is the same reasoning that made the original VPS build resumable. The
difference is that it is now cut deliberately and the cuts are in the repo.

## 10. enroot does not carry the Dockerfile's `ENV`

A symlinked binary starts and cannot find its own `.so` files, because
`LD_LIBRARY_PATH` from the Dockerfile is not part of the container. The build
therefore installs a **wrapper**, not a symlink:

```sh
#!/bin/sh
export LD_LIBRARY_PATH=/app/build/Build/lib:/opt/coin-or/lib:${LD_LIBRARY_PATH:-}
exec /app/build/Build/bin/HexMeshing "$@"
```

That is what makes plain `HexMeshing` work from `run_algohex.py` with no
environment setup, on either backend.

Side effect worth knowing: `ldd "$(command -v HexMeshing)"` now inspects a
shell script and reports "not a dynamic executable". Check the real binary
path.

## 11. Traps in our own tooling, recorded so they are not re-derived

Not cluster problems, but they cost jobs:

* **`grep -c` exits 1 when it matches nothing**, and "no missing libraries" is
  the success case. Inside an `&&` chain that aborts the check on a healthy
  image — job 6858990 failed its own precheck for this reason. Use `| wc -l`,
  which always exits 0.
* **A command substitution is a subshell.** `J=$(sub ...)` cannot set a
  variable in its caller, so a failure flag assigned inside `sub` was invisible
  and the script reported success after four failed submits. Return an exit
  status instead; that survives.
* **`[ -d "$REPO/.git" ]` is false in a worktree**, where `.git` is a file
  holding a gitdir pointer. Use `git -C "$REPO" rev-parse --git-dir`.
* **The default branch is `master`** and none of this work is there. Clone with
  `--branch`, and check the branch rather than printing it — `git checkout
  origin/<branch>` leaves a detached HEAD where `rev-parse --abbrev-ref HEAD`
  says `HEAD` and a later `git pull` fails.

## 12. OpenFOAM's `eval` runs your command twice — pass `--rc`

The dtOO image ships a command script that sources OpenFOAM's
`etc/bashrc`, and OpenFOAM's `etc/config.sh/functions` runs `eval` on the
argv enroot handed over. Two measured consequences (sibling repo, 2026-09-04,
their commit `fae9a3d`):

* **Quoting is destroyed.** `bash -c 'source env.sh; checkMesh ...'` is re-read
  as two outer statements: `checkMesh` ran before its environment existed and
  failed with `libfiniteVolume.so: cannot open shared object file`, while the
  second, correct execution succeeded — and the **first** one decides the exit
  code.
* **Everything runs twice.** Sourcing OpenFOAM's bashrc from a shell that had
  already `cd`-ed somewhere turned the repeat into an endless loop that
  produced no output at all.

The fix is one flag: `--rc <path>/cluster/enroot_rc.sh`, whose whole body is
`exec "$@"`. The path must resolve **inside** the container, so it lives in
the mounted eigenfrequencies checkout (we mount it at `/ef` and pass
`--rc /ef/cluster/enroot_rc.sh`). The container environment is unaffected:
`LD_LIBRARY_PATH` is byte-identical with and without the script, OpenFOAM's
own libraries included.

## 13. Logs are committed to git, not left in scratch

The sibling repo keeps `cluster/logs/<jobid>/` with the SLURM `.out` and one
log per stage, commits new ones with `git add cluster/logs && git commit`, and
symlinks `latest`. This repo now does the same through
`experimentell/hex3d_algohex/logs/`: `batch_dtoo_export.slurm` tees its whole
output — guards included — to `batch_dtoo_export_<jobid>.log` and prints the
git command at the end. Rationale, in their words: a failed run must leave its
evidence somewhere greppable, not only in a `.out` file that the next
submission pushes out of `ls`. The checksum of that habit is that every guard
failure above was diagnosed from a log that survived the job.

## Measured costs, for planning

| step | time | resources |
|---|---|---|
| sources into the container (login node) | 10-20 min | download-bound |
| Ipopt + MUMPS | one 30-min slot | 64 cores, mostly serial |
| Bonmin | 345 s | 64 cores |
| AlgoHex cmake + ninja + export | 8:47 | **14.55 GB peak**, 45 % of 32 G |
| `enroot create` (2.1 GB image) | ~4 min | once, ever |
| one sample at `-n 2000` | **12:25** | 1.87 GB, 18 % of 8 cores |

AlgoHex is largely serial: 18 % efficiency on 8 cores is about 1.5 cores in
use. That number is what makes per-node concurrency the right lever for
throughput rather than more cores per run.
