"""Stage: extract a hex mesh from an AlgoHex seamless map with C4HexMeshing.

AlgoHex's own extractor (HexEx) rejects some seamless maps: "Detected N Proper,
M Flipped ... Tets after the Sanitization" ends the run with "Invalid Input -
Flipped or Degenerate Tet", although the seamless map itself is fine. C4HexMeshing
(Brueckler et al. 2023) takes that same map, builds the motorcycle complex (MC3D),
optionally collapses its zero-extent cells and quantises the complex into an
integer-grid map of its own, then extracts hexahedra. So it can succeed where the
HexEx extraction cannot, and it is a second, independent extractor rather than a
replacement of AlgoHex's field computation.

Input: the ``.hexex`` seamless map that ``run_algohex.py`` writes with
``--checkpoints`` (HexMeshing ``--sm-out-path``). Output: an OpenVolumeMesh hex
mesh, ``<prefix>_<n>_hex.ovm`` -- the format the beam collapse reads directly.

The ``c4hex_cli`` binary is external (GPLv3, built at ``/opt/stack/ext/C4HexMeshing``,
no LD_LIBRARY_PATH needed); select another one with ``--bin`` or ``$HEXDP_C4HEX_BIN``.
"""

from __future__ import annotations

import argparse
import glob
import os
import subprocess
import sys
import time
from pathlib import Path

DEFAULT_BIN = "/opt/stack/ext/C4HexMeshing/build/Build/bin/c4hex_cli"


def c4hex_bin(explicit=None):
    """The c4hex_cli binary; a missing one is an error with the three ways to set it."""
    p = explicit or os.environ.get("HEXDP_C4HEX_BIN") or DEFAULT_BIN
    if not Path(p).is_file():
        sys.exit(f"[c4hex] error: {p} not found (set --bin, $HEXDP_C4HEX_BIN, or build "
                 f"C4HexMeshing - see docs/backends.md)")
    return p


def run(hexex, out_dir, prefix="hex", bin=None, flags=(), timeout=None, log_name="c4hex.log"):
    """Run c4hex_cli on a seamless map and return the .ovm it wrote.

    ``flags`` are passed through to the CLI (e.g. ``--collapse``, ``--num-hexes``);
    they are not interpreted here, so a new C4Hex option needs no change in hexdp.
    A run that produces no hex mesh is an error with the log tail, because c4hex_cli
    reports its own failures (IGM untangling, non-convergence) as exit codes 20-22.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    log = out_dir / log_name
    cmd = [str(c4hex_bin(bin)), "--input", str(hexex), "--output-hex-file", str(out_dir / prefix),
           *map(str, flags)]
    print(f"[c4hex] {' '.join(cmd)}", flush=True)
    t0 = time.time()
    with open(log, "w") as fh:
        try:
            r = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, timeout=timeout)
        except subprocess.TimeoutExpired:
            sys.exit(f"[c4hex] timeout after {timeout} s (log {log})")
    dt = time.time() - t0
    found = sorted(glob.glob(str(out_dir / f"{prefix}_*_hex.ovm")))
    if r.returncode != 0 or not found:
        tail = "".join(open(log).readlines()[-15:])
        print(f"[c4hex] FAILED after {dt:.0f} s, exit {r.returncode}, {len(found)} hex mesh(es) "
              f"(log {log}):\n{tail}", file=sys.stderr)
        sys.exit(3)
    if len(found) > 1:
        sys.exit(f"[c4hex] error: {len(found)} hex meshes for one run, cannot choose: {found}")
    print(f"[c4hex] hex mesh {found[0]} ({dt:.0f} s)", flush=True)
    return Path(found[0])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--input", required=True, help="seamless map (.hexex) from HexMeshing --sm-out-path")
    ap.add_argument("--out-dir", default=".", help="where the .ovm and the log are written")
    ap.add_argument("--prefix", default="hex", help="output name: <prefix>_<n>_hex.ovm")
    ap.add_argument("--bin", default=None, help=f"c4hex_cli (default $HEXDP_C4HEX_BIN or {DEFAULT_BIN})")
    ap.add_argument("--timeout", type=float, default=None, help="kill the run after this many seconds")
    ap.add_argument("rest", nargs=argparse.REMAINDER, help='flags for c4hex_cli, after "--"')
    a = ap.parse_args(argv)
    ovm = run(a.input, a.out_dir, a.prefix, a.bin, [x for x in a.rest if x != "--"], a.timeout)
    print(f"[c4hex] {ovm}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
