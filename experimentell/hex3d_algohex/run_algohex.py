"""Stage 2: run AlgoHex's HexMeshing on a tet mesh from tet_prep*.py.

Not a host build: there is no cmake/ninja on this host and the IPOPT/Bonmin
chain is fragile, so AlgoHex runs in a container. TWO container runtimes are
supported behind one call site (decision H): Docker on the user's own
6-machine cluster, enroot on bwUniCluster 3.0.

**The image is self-contained** (T2). `algohex:portable` carries `HexMeshing`
at `/opt/algohex/Build/bin`, symlinked onto `PATH`, with
`LD_LIBRARY_PATH=/opt/algohex/Build/lib:/opt/coin-or/lib`; its recipe is
`external_patches/Dockerfile.portable`. The old `algohex-build-cache` volume
mount is gone, and with it the reason this script could not leave the VPS:
`docker save algohex-configured` produced an image WITHOUT the binary, because
volume contents are not part of an image.

Docker needs `--network=host`: the default bridge network cannot resolve DNS
in this sandbox (see PROGRESS.md).
"""

import shutil
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent.parent
IN_VTK = REPO / "data" / "T1_9" / "T1_9_tet.vtk"
OUT_DIR = REPO / "output" / "hex3d_algohex"

DOCKER_IMAGE = "algohex:portable"
# enroot takes EITHER a .sqsh path or the name of an already-created container,
# and which one you want depends on how many jobs will run -- a distinction
# worth stating, because getting it backwards is expensive on a cluster.
#
# One job: pass the .sqsh. `enroot start` reads it through squashfuse, no
# `enroot create`, no rootfs unpack.
#
# An ARRAY of jobs: create the container ONCE into a persistent store and pass
# its NAME. Otherwise every task pays the squashfuse read across the parallel
# filesystem again. Measured in the sibling project on bwUniCluster
# (2026-09-07, job 6822816, dev_cpu_il): unpacking 9.6 GB of .sqsh from Lustre
# took ~20 of 30 budgeted minutes at 20 s of CPU -- pure IO wait -- and an
# earlier "8 s cold, 3.8 s warm" figure did not reproduce. See the staging
# block in eigenfrequencies/cluster/submit_hydroflow_opt.sh, which is why T10
# uses a shared store and pays each image once ever. This image is 2.1 GB, so
# scale that cost down accordingly, but pay it once.
#
# The naming trap from cluster/enroot_dtoo_import.md applies to the second form
# only: `enroot create` derives the container name from the .sqsh basename, so
# algohex.sqsh -> "algohex".
ENROOT_IMAGE = REPO / "output" / "hex3d_algohex" / "algohex.sqsh"
HEXMESHING_BIN = "HexMeshing"          # on PATH inside the image, since T2
WORK = "/work"                          # the repo, inside the container


def container_cmd(backend, args, image=None, cpus=None, mounts=()):
    """One command line for either runtime. `mounts` is (host, container).

    The only real asymmetry is `--cpus`. Docker can cap a container directly;
    enroot cannot, because it is not a resource manager -- under SLURM that is
    `--cpus-per-task`, and on a bare box it is `/root/bin/capped` around this
    process. Passing `--cpus` with `--backend enroot` is therefore refused
    rather than ignored: a silently unlimited AlgoHex is how this box got
    throttled in the first place.
    """
    binds = list(mounts) or [(str(REPO), WORK)]
    if backend == "docker":
        cmd = ["docker", "run", "--rm", "--network=host"]
        if cpus:
            cmd += ["--cpus", str(cpus)]
        for h, c in binds:
            cmd += ["-v", f"{h}:{c}"]
        return cmd + [image or DOCKER_IMAGE, HEXMESHING_BIN, *args]
    if backend == "enroot":
        if cpus:
            raise SystemExit(
                "[run_algohex] --cpus is docker-only; enroot does not limit "
                "CPU. Use SLURM's --cpus-per-task, or wrap this process in "
                "/root/bin/capped on a bare box.")
        # A bare word is a container name (already `enroot create`d); anything
        # that looks like a path must exist as a .sqsh. Checking the file only
        # would reject the array-friendly form outright.
        if image and "/" not in image and not image.endswith(".sqsh"):
            img = image
        else:
            p = Path(image) if image else ENROOT_IMAGE
            if not p.exists():
                raise SystemExit(
                    f"[run_algohex] no enroot image at {p}. Build one with "
                    f"scripts/export_algohex_enroot.sh, or pass the name of "
                    f"an already-created container instead of a path.")
            img = str(p)
        # --root: uid 0 inside, --rw: the image's own root is writable, which
        # AlgoHex needs for temporaries. Both verified in T1 together with a
        # writable bind mount.
        cmd = ["enroot", "start", "--root", "--rw"]
        for h, c in binds:
            cmd += ["--mount", f"{h}:{c}"]
        return cmd + [img, HEXMESHING_BIN, *args]
    raise SystemExit(f"[run_algohex] unknown backend {backend!r}")


def _work_path(p):
    """Host path -> its path inside the container.

    Only the repository is mounted, at /work, so anything outside it is
    invisible to AlgoHex. Saying so here beats a `relative_to` ValueError from
    four minutes into a run, which is how this was found.
    """
    p = Path(p).resolve()
    try:
        rel = p.relative_to(REPO)
    except ValueError:
        raise SystemExit(
            f"[run_algohex] {p} is outside the repository, and only the "
            f"repository is mounted (at {WORK}). Put inputs and outputs under "
            f"{REPO}, or extend the mounts in container_cmd.")
    return f"{WORK}/{rel}"


def run_hexmeshing(extra_args=(), tag="", in_vtk=None, prefix="T1_9",
                   cpus=None, backend="docker", image=None, out_dir=None,
                   checkpoints=False):
    """``tag`` suffixes every output so concurrent runs cannot clobber each
    other -- though note that two AlgoHex runs must never overlap anyway
    (~6 GB peak each; concurrency is what killed v7 and v8).

    ``out_dir`` is the sample directory. It defaults to the shared
    `output/hex3d_algohex/` only for backwards compatibility with the older
    scripts; batch runs must pass it, because thousands of samples writing
    `gen_hex_<name>.ovm` into one directory is the tidiness bug decision J
    flagged.
    """
    out_dir = Path(out_dir).resolve() if out_dir else OUT_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    sfx = f"_{tag}" if tag else ""
    in_vtk = (Path(in_vtk).resolve() if in_vtk else IN_VTK)
    out_ovm = out_dir / f"{prefix}_hex{sfx}.ovm"
    json_out = out_dir / f"{prefix}_hex_metrics{sfx}.json"
    log_path = out_dir / f"hexmeshing{sfx}.log"

    args = ["-i", _work_path(in_vtk),
            "-o", _work_path(out_ovm),
            "-j", _work_path(json_out)]

    # The seamless map and the post-singularity tet mesh used to be written
    # unconditionally, to make -n sweeps cheap by replaying them. That does
    # not work and was measured not to work: AlgoHex writes the tet mesh as
    # binary OVM and its own -i reader rejects it ("Checkpoint reuse does not
    # work", DATA_GENERATION.md). They cost ~45 MB and ~11 MB per run for
    # nothing -- 560 GB at 10 000 samples -- so they are opt-in now. Still
    # wanted for the singular-graph figure (ANALYSIS_PLAN step 04).
    if checkpoints:
        args += ["--sm-out-path",
                 _work_path(out_dir / f"{prefix}_seamless{sfx}.hexex"),
                 "--final-tetmesh-out-path",
                 _work_path(out_dir / f"{prefix}_final_tet{sfx}.ovm")]

    cmd = container_cmd(backend, [*args, *extra_args], image=image, cpus=cpus)
    print(f"[run_algohex] backend={backend}")
    print(f"[run_algohex] {' '.join(cmd)}")
    t0 = time.time()
    with open(log_path, "w") as log:
        proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
    dt = time.time() - t0
    print(f"[run_algohex] exit={proc.returncode} in {dt:.1f}s, log: {log_path}")

    tail = log_path.read_text().splitlines()[-40:]
    print("[run_algohex] --- log tail ---")
    print("\n".join(tail))

    return proc.returncode, out_ovm, json_out, log_path


def check_backend(backend):
    """Fail early and with the reason, not mid-run."""
    need = {"docker": "docker", "enroot": "enroot"}[backend]
    if shutil.which(need) is None:
        raise SystemExit(f"[run_algohex] backend {backend!r} needs {need!r} "
                         f"on PATH and it is not there")
    return True


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tag", default="", help="suffix for all output files")
    ap.add_argument("--in-vtk", default=None, help="input tet mesh .vtk")
    ap.add_argument("--cpus", type=float, default=None,
                    help="limit the container to N CPUs (docker only)")
    ap.add_argument("--prefix", default="T1_9",
                    help="output filename prefix; T1_9 was hard-coded while "
                         "that was the only geometry")
    ap.add_argument("--backend", choices=("docker", "enroot"),
                    default="docker",
                    help="docker on the 6-machine cluster, enroot on "
                         "bwUniCluster 3.0 (decision H)")
    ap.add_argument("--image", default=None,
                    help=f"docker tag (default {DOCKER_IMAGE}) or path to a "
                         f".sqsh (default {ENROOT_IMAGE.name})")
    ap.add_argument("--out-dir", default=None,
                    help="write outputs here instead of output/hex3d_algohex/"
                         "; batch runs should pass the sample directory")
    ap.add_argument("--checkpoints", action="store_true",
                    help="also write the seamless map and the intermediate "
                         "tet mesh (~56 MB). Off by default: replaying them "
                         "does not work, see DATA_GENERATION.md")
    ap.add_argument("rest", nargs=argparse.REMAINDER,
                    help="extra flags passed straight to HexMeshing")
    a = ap.parse_args()
    check_backend(a.backend)
    extra = [x for x in a.rest if x != "--"]
    rc, out_ovm, json_out, log_path = run_hexmeshing(
        extra, tag=a.tag, in_vtk=a.in_vtk, prefix=a.prefix, cpus=a.cpus,
        backend=a.backend, image=a.image, out_dir=a.out_dir,
        checkpoints=a.checkpoints)
    print(f"[run_algohex] hex mesh: {out_ovm}")
    sys.exit(rc)
