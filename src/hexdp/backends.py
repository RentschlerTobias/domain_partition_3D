"""Optional compute backends and the hardware hexdp can use - what is there, and what depends on it.

    status = environment()          # {"hardware": {...}, "backends": {name: {...}}}
    available("algohex")            # True / False (cached)

The core (gmsh, numpy, scipy, meshio) is required; the backends below are optional. A strategy whose backend is
missing says so in applicable() (Score.no(..., missing=("backend ...",))), so an agent sees it before running.
Install notes and facts per backend: docs/backends.md.
"""

from __future__ import annotations

import functools
import os
import shutil
import subprocess
from pathlib import Path

BACKENDS = {
    "gmsh": {"needed_for": "geometry loading, tet meshes, .msh export (core)", "required": True},
    "algohex": {"needed_for": "strategy 'algohex' (general topology: frame field + integer-grid map); the tistos "
                              "fluid channel", "required": False,
                "facts": "build from source with many dependencies (CoMISo, OpenVolumeMesh, HexEx, QGP3D, MC3D, "
                         "...), about 30-60 min; per geometry (tistos, n=2000, 8 threads) median 35 min, 90 % "
                         "under 76 min, single runs 3-6 h; about 9 % fail without a result (invalid parametrisation, the "
                         "4 h limit, crashes) and about 2 % give a mesh off the surface. Try sweep / polycube / "
                         "templates first."},
    "occ": {"needed_for": "strategies 'sweep', 'cgrid' (exact CAD sections), block 'dtoo_rotor' (with dtoo)",
            "required": False},
    "dtoo": {"needed_for": "parametric dtOO machines: block 'dtoo_rotor' (runner from a parameter file)",
             "required": False},
    "evocube": {"needed_for": "strategy 'polycube' (polycube labeling; GPL-3, used as an external program)",
                "required": False},
    "mmg": {"needed_for": "envelope cut before AlgoHex (blunt trailing edges)", "required": False},
}


def _occ_python():
    return os.environ.get("HEXDP_OCC_PYTHON", "/opt/venv-dtoo/bin/python")


@functools.lru_cache(maxsize=None)
def _probe(name: str) -> dict:
    if name == "gmsh":
        try:
            import gmsh
            return {"ok": True, "version": getattr(gmsh, "__version__", "?")}
        except Exception as exc:                                      # noqa: BLE001
            return {"ok": False, "why": str(exc)}
    if name == "algohex":
        given = os.environ.get("HEXMESHING_BIN", "")
        b = given or "HexMeshing"
        p = shutil.which(b) or (b if Path(b).exists() else None)
        if not p:
            why = (f"$HEXMESHING_BIN={given}: not found" if given
                   else "HexMeshing not found on PATH (set $HEXMESHING_BIN to the binary of a patched build)")
            return {"ok": False, "path": None, "why": why}
        out = {"ok": True, "path": os.path.abspath(p), "how": "$HEXMESHING_BIN" if given else "PATH"}
        if not given:
            out["notes"] = ("found on PATH; set $HEXMESHING_BIN to the binary of a patched build to use that "
                            "instead (e.g. a build with extra AlgoHex flags)")
        return out
    if name in ("occ", "dtoo"):
        py = _occ_python()
        if not Path(py).exists():
            return {"ok": False, "why": f"no OCC python at {py} ($HEXDP_OCC_PYTHON)"}
        mod = "OCC.Core.TopoDS" if name == "occ" else "dtOOPythonSWIG"
        env = dict(os.environ, LD_LIBRARY_PATH=os.environ.get("HEXDP_OCC_LD_LIBRARY_PATH", "/dtOO-install/lib")
                   + os.pathsep + os.environ.get("LD_LIBRARY_PATH", ""))
        try:
            r = subprocess.run([py, "-c", f"import {mod}"], capture_output=True, timeout=60, env=env)
            return {"ok": r.returncode == 0, "python": py,
                    "why": None if r.returncode == 0 else r.stderr.decode(errors="replace").strip()[-200:]}
        except Exception as exc:                                      # noqa: BLE001
            return {"ok": False, "why": str(exc)}
    if name == "evocube":
        d = Path(os.environ.get("HEXDP_EVOCUBE_BIN", "/opt/stack/ext/evocube/build"))
        ok = all((d / t).exists() for t in ("prep_cli", "evolabel_cli", "polycube_withHexEx"))
        return {"ok": ok, "path": str(d), "why": None if ok else f"evocube programs not in {d} ($HEXDP_EVOCUBE_BIN)"}
    if name == "mmg":
        p = os.environ.get("MMG3D", "/opt/stack/tools/mmg-install/bin/mmg3d_O3")
        p = shutil.which(p) or (p if Path(p).exists() else None)
        return {"ok": bool(p), "path": p, "why": None if p else "mmg3d not found ($MMG3D)"}
    raise KeyError(name)


def available(name: str) -> bool:
    return bool(_probe(name)["ok"])


def cpus_available() -> list[int]:
    """The CPUs this process may run on (its affinity: a launcher like taskset restricts it from outside)."""
    try:
        return sorted(os.sched_getaffinity(0))
    except AttributeError:                                           # not Linux
        return list(range(os.cpu_count() or 1))


def parse_cpus(spec) -> list[int]:
    """'4' -> the first 4 available CPUs; '0-3,8' -> exactly those (must be available)."""
    avail = cpus_available()
    s = str(spec).strip()
    if s.isdigit() and "-" not in s and "," not in s:
        n = int(s)
        if not 1 <= n <= len(avail):
            raise ValueError(f"cpus={n}: between 1 and {len(avail)} available")
        return avail[:n]
    out = set()
    for part in s.split(","):
        a, _, b = part.partition("-")
        out |= set(range(int(a), int(b or a) + 1))
    missing = sorted(out - set(avail))
    if missing:
        raise ValueError(f"cpus {missing} are not available to this process (available: {compact(avail)})")
    return sorted(out)


def compact(cpus) -> str:
    """[0, 1, 2, 5] -> '0-2,5' (a taskset list)."""
    cpus = sorted(cpus)
    out, i = [], 0
    while i < len(cpus):
        j = i
        while j + 1 < len(cpus) and cpus[j + 1] == cpus[j] + 1:
            j += 1
        out.append(str(cpus[i]) if i == j else f"{cpus[i]}-{cpus[j]}")
        i = j + 1
    return ",".join(out)


def environment() -> dict:
    avail = cpus_available()
    try:
        load = os.getloadavg()
    except OSError:
        load = None
    mem = None
    try:
        info = dict(line.split(":", 1) for line in open("/proc/meminfo"))
        mem = {"total_gb": round(int(info["MemTotal"].split()[0]) / 1e6, 1),
               "available_gb": round(int(info["MemAvailable"].split()[0]) / 1e6, 1)}
    except Exception:                                                  # noqa: BLE001
        pass
    return {"hardware": {"cpus_available": len(avail), "cpu_list": compact(avail), "machine_cpus": os.cpu_count(),
                         "load_1_5_15": load, "memory": mem},
            "backends": {n: {**meta, **_probe(n)} for n, meta in BACKENDS.items()}}
