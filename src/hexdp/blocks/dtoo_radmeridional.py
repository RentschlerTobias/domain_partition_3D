"""dtOO stage for the parametrised radial machine of the dtOO demo ``radMeridional`` (dtOO environment): build the
machine for a parameter file and write one gmsh mesh per bounded volume.

    python -m hexdp.blocks.dtoo_radmeridional PARAMS.json --out DIR [--demo DIR] [--config ConfigOptimization]
                                              [--parts gv_mesh,ru_mesh,meshLayers]

PARAMS.json: {"configMeas": {...}, "configMeridional": {...}, "configGuideVane": {...}, "configRunner": {...},
"configLayer": {...}} - each entry overrides the given keys of the demo configuration (``Config`` of --config);
``{}`` builds the defaults. Writes DIR/<part>.msh (gmsh 2.2, physical groups <part>_hub/_shroud/_blade/_inlet/_outlet
/_suction/_pressure) and DIR/<part>.vtk per part, DIR/config.json (the full configuration used + timings).

The guide vane is built only when gv_mesh is asked for, the runner (and the draft cone, which follows it) only for
ru_mesh / meshLayers. A dtOO built without its OpenFOAM interface still works: the demo's builder package imports
the OpenFOAM classes eagerly, so missing names are replaced by placeholders that fail when used.
"""

from __future__ import annotations

import argparse
import importlib
import json
import os
import re
import sys
import time
from pathlib import Path

DEMO = "/opt/stack/ext/dtOO/demo/radMeridional"
PARTS = ("gv_mesh", "ru_mesh", "meshLayers")


def import_builders():
    """Import dtOOPythonApp.builder, adding placeholders for dtOOPythonSWIG names this dtOO build lacks."""
    import dtOOPythonSWIG
    missing = []
    for _ in range(50):
        try:
            importlib.import_module("dtOOPythonApp.builder")
            return missing
        except ImportError as exc:
            m = re.search(r"cannot import name '(\w+)' from 'dtOOPythonSWIG'", str(exc))
            if not m:
                raise
            name = m.group(1)

            def _absent(*a, _n=name, **k):
                raise RuntimeError(f"dtOOPythonSWIG.{_n} is not part of this dtOO build (OpenFOAM interface?)")

            setattr(dtOOPythonSWIG, name, type(name, (), {"__init__": _absent}))
            missing.append(name)
            for mod in [k for k in sys.modules if k.startswith("dtOOPythonApp")]:
                del sys.modules[mod]
    raise ImportError("dtOOPythonApp.builder: too many missing dtOOPythonSWIG names")


def plain(x):
    if isinstance(x, (int, float, str, bool, type(None))):
        return x
    if isinstance(x, (list, tuple)):
        return [plain(y) for y in x]
    if isinstance(x, dict):
        return {k: plain(v) for k, v in x.items()}
    return str(x)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("params")
    ap.add_argument("--out", required=True)
    ap.add_argument("--demo", default=os.environ.get("HEXDP_DTOO_RADMERIDIONAL", DEMO))
    ap.add_argument("--config", default="ConfigOptimization")
    ap.add_argument("--parts", default=",".join(PARTS))
    a = ap.parse_args(argv)
    parts = [p for p in a.parts.split(",") if p]
    bad = set(parts) - set(PARTS)
    if bad:
        ap.error(f"unknown parts {sorted(bad)} (known: {', '.join(PARTS)})")
    params = json.loads(Path(a.params).read_text()) if a.params != "-" else {}
    out = Path(a.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(Path(a.demo).resolve()))
    t0 = time.time()
    times = {}
    missing = import_builders()
    CO = importlib.import_module(a.config)
    import radMeridional
    cfg = CO.Config()
    meas, mer, gv, ru, lay = cfg.getConfig()
    for name, d in (("configMeas", meas), ("configMeridional", mer), ("configGuideVane", gv), ("configRunner", ru),
                    ("configLayer", lay)):
        unknown = set(params.get(name, {})) - set(d)
        if unknown:
            raise KeyError(f"{name}: unknown keys {sorted(unknown)}")
        d.update(params.get(name, {}))
    os.chdir(out)                                   # dtOO writes build.log and debug meshes into the cwd
    g = radMeridional.radMeridional()
    g.createMeridional(mer, *CO.ConfigOptimization(meas).getCurves())
    times["meridional"] = time.time() - t0
    if "gv_mesh" in parts:
        g.createBlade(gv)
        times["guide_vane"] = time.time() - t0
    if "ru_mesh" in parts or "meshLayers" in parts:
        g.createBlade(ru)
        times["runner"] = time.time() - t0
        g.createLayerRegion(lay)
        times["layer"] = time.time() - t0
    bV = g.getContainer().cptr_bV()
    for p in parts:
        v = bV[p]
        v.makeGrid()
        gm = v.getModel()
        gm.writeMSH(str(out / f"{p}.msh"), 2.2)
        gm.writeVTK(str(out / f"{p}.vtk"))
        times[p] = time.time() - t0
        print(f"[dtoo_radmeridional] {p} written after {times[p]:.0f} s", flush=True)
    for f in out.glob("*_debug.msh"):
        f.unlink()
    for f in out.glob("dtMeshGRegion_*"):
        f.unlink()
    (out / "config.json").write_text(json.dumps(
        {"demo": str(a.demo), "config": a.config, "params": params, "missing_swig_names": missing,
         "configMeas": plain(meas), "configMeridional": plain(mer), "configGuideVane": plain(gv),
         "configRunner": plain(ru), "configLayer": plain(lay), "times_s": times}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
