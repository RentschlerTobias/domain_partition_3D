"""dtOO stage of the building block ``dtoo_rotor`` (dtOO environment): build a dtOO machine for a parameter file and sample the
blade and fillet surfaces into rotor.npz (see build.py for the format).

    python -m hexdp.blocks.dtoo_rotor PARAMS.json --case /dtOO/demo/tistos --out rotor.npz

The machine template is copied to <out dir>/dtoo_case (dtOO reads and writes relative to it). The defaults are
the names of the tistos machine: the blade + fillet surfaces are the inner faces (v = 0) of the volumes
PREFIX{Blade,Rounding}_meshBlock_closed_{patch}_rot, the patches in cyclic order around the blade; the number of
blades and the hub length come from constants of the machine.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
import time
from pathlib import Path

import numpy as np


def build(case_dir: Path, params: dict, adjust_plugin: str | None):
    from dtOOPythonSWIG import (baseContainer, dtXmlParser, labeledVectorHandlingAnalyticFunction,
                                labeledVectorHandlingAnalyticGeometry, labeledVectorHandlingBoundedVolume,
                                labeledVectorHandlingConstValue, labeledVectorHandlingDtCase,
                                labeledVectorHandlingDtPlugin, logMe)
    os.chdir(case_dir)
    logMe.initLog("build.log")
    parser = dtXmlParser.init("machine.xml", "machineSave.xml").reference()
    parser.parse()
    cV = labeledVectorHandlingConstValue()
    parser.createConstValue(cV)
    parser.loadStateToConst("templateState", cV)
    for k, v in params.items():
        cV.get(k).setValue(float(v))

    def once():
        bC, aF, aG = baseContainer(), labeledVectorHandlingAnalyticFunction(), labeledVectorHandlingAnalyticGeometry()
        bV, dC, dP = labeledVectorHandlingBoundedVolume(), labeledVectorHandlingDtCase(), labeledVectorHandlingDtPlugin()
        parser.createAnalyticFunction(bC, cV, aF)
        parser.createAnalyticGeometry(bC, cV, aF, aG)
        parser.createBoundedVolume(bC, cV, aF, aG, bV)
        parser.createPlugin(bC, cV, aF, aG, bV, dC, dP)
        return cV, aG, dP
    cV, aG, dP = once()
    if adjust_plugin:
        dP.get(adjust_plugin).apply()
        cV, aG, dP = once()
    return cV, aG


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("params")
    ap.add_argument("--case", default="/dtOO/demo/tistos", help="dtOO machine template directory")
    ap.add_argument("--out", required=True)
    ap.add_argument("--prefix", default="tV_ruWithRounding")
    ap.add_argument("--patches", default="0_0,le_0,le_1,0_2,te_0,te_1", help="patch names in cyclic order")
    ap.add_argument("--nblades-const", default="cV_ru_nBlades")
    ap.add_argument("--z-consts", default="cV_ruMech_z_hub_0,cV_ruMech_z_hub_1")
    ap.add_argument("--adjust-plugin", default="ru_adjustDomain", help="plugin applied once before the final "
                    "build ('' = none)")
    ap.add_argument("--nu", type=int, default=201, help="samples along a patch")
    ap.add_argument("--nw", type=int, default=61, help="samples from the hub to the tip")
    a = ap.parse_args(argv)
    t0 = time.time()
    out = Path(a.out).absolute()
    params = json.loads(Path(a.params).read_text())
    case_dir = out.parent / "dtoo_case"
    if not case_dir.exists():
        shutil.copytree(a.case, case_dir)
    cV, aG = build(case_dir, params, a.adjust_plugin or None)
    from dtOOPythonSWIG import map3dTo3d
    print(f"[rotor] dtOO build {time.time() - t0:.1f} s", flush=True)
    U, W = np.linspace(0, 1, a.nu), np.linspace(0, 1, a.nw)
    surf = {}
    for part, key in (("Rounding", "fillet"), ("Blade", "blade")):
        stack = []
        for p in a.patches.split(","):
            m = map3dTo3d.MustDownCast(aG.get(f"{a.prefix}{part}_meshBlock_closed_{p}_rot"))
            X = np.empty((a.nu, a.nw, 3))
            for i, u in enumerate(U):
                for j, w in enumerate(W):
                    q = m.getPointPercent(u, 0.0, w)
                    X[i, j] = (q.x(), q.y(), q.z())
            stack.append(X)
        surf[key] = np.array(stack)
    z0, z1 = (float(cV.get(c).getValue()) for c in a.z_consts.split(","))
    nb = int(round(cV.get(a.nblades_const).getValue()))
    r_hub = float(np.hypot(surf["fillet"][:, :, 0, 0], surf["fillet"][:, :, 0, 1]).mean())
    np.savez(out, fillet=surf["fillet"], blade=surf["blade"], r_hub=r_hub, z0=z0, z1=z1, nb=nb,
             params=json.dumps(params), patches=a.patches)
    print(f"[rotor] {out}: r_hub {r_hub:.4f}, z {z0}..{z1}, {nb} blades, patches {surf['blade'].shape} "
          f"({time.time() - t0:.1f} s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
