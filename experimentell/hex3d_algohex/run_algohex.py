"""Stage 2 (feat/algohex-3d-frame-field): run AlgoHex's HexMeshing on the
T1_9 tet mesh from tet_prep.py (see PLAN.md).

Docker, not a host build (no cmake/ninja on this host, and the IPOPT/Bonmin
build chain is fragile -- the vendored Dockerfile in external/algohex-src is
the CI-validated build path). The container needs --network=host: the
default bridge network can't resolve DNS in this sandbox (see PROGRESS.md).
"""

import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent.parent
IN_VTK = REPO / "data" / "T1_9" / "T1_9_tet.vtk"
OUT_DIR = REPO / "output" / "hex3d_algohex"
# No baked "algohex" image with symlinked binaries: the ninja build was run
# resumably against a persistent volume (see PROGRESS.md) rather than as a
# single RUN layer, so the binary lives in that volume, not the image.
IMAGE = "algohex-configured"
BUILD_VOLUME = "algohex-build-cache"
HEXMESHING_BIN = "/app/build/Build/bin/HexMeshing"


def run_hexmeshing(extra_args=(), tag="", in_vtk=None, prefix="T1_9"):
    """``tag`` suffixes every output so concurrent runs (e.g. a fast
    --without-integrable-field-optimization probe alongside the full run)
    cannot clobber each other."""
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    sfx = f"_{tag}" if tag else ""
    in_vtk = (Path(in_vtk).resolve() if in_vtk else IN_VTK)
    out_ovm = OUT_DIR / f"{prefix}_hex{sfx}.ovm"
    json_out = OUT_DIR / f"{prefix}_hex_metrics{sfx}.json"
    log_path = OUT_DIR / f"hexmeshing{sfx}.log"

    # Checkpoint the seamless map + the post-singularity-optimization tetmesh.
    # AlgoHex's hexMeshing_from_seamless_map() path (triggered by passing BOTH
    # --hexex-in-path and -i) skips field generation, locally-meshable-field
    # generation and integrability optimization -- together 86% of the runtime
    # on T1_9 (6704s of 7777s) -- and goes straight to parametrization +
    # quantization + extraction. Saving these makes -n sweeps cheap.
    sm_out = OUT_DIR / f"{prefix}_seamless{sfx}.hexex"
    final_tet = OUT_DIR / f"{prefix}_final_tet{sfx}.ovm"

    cmd = [
        "docker", "run", "--rm", "--network=host",
        "-v", f"{REPO}:/work",
        "-v", f"{BUILD_VOLUME}:/app/build",
        IMAGE, HEXMESHING_BIN,
        "-i", f"/work/{in_vtk.relative_to(REPO)}",
        "-o", f"/work/{out_ovm.relative_to(REPO)}",
        "-j", f"/work/{json_out.relative_to(REPO)}",
        "--sm-out-path", f"/work/{sm_out.relative_to(REPO)}",
        "--final-tetmesh-out-path", f"/work/{final_tet.relative_to(REPO)}",
        *extra_args,
    ]
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


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tag", default="", help="suffix for all output files")
    ap.add_argument("--in-vtk", default=None, help="input tet mesh .vtk")
    ap.add_argument("--prefix", default="T1_9",
                    help="output filename prefix; T1_9 was hard-coded while "
                         "that was the only geometry")
    ap.add_argument("rest", nargs=argparse.REMAINDER,
                    help="extra flags passed straight to HexMeshing")
    a = ap.parse_args()
    extra = [x for x in a.rest if x != "--"]
    rc, out_ovm, json_out, log_path = run_hexmeshing(
        extra, tag=a.tag, in_vtk=a.in_vtk, prefix=a.prefix)
    sys.exit(rc)
