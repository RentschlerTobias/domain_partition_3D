"""Export the true walls from the dtOO machine.

Run in blade/machine_scratch (machine.xml, machineSave.xml, xml/, params.json):
  machine_blade.npz / .vtk   aS_ru_blade_rot, dense parametric samples + span cuts
  machine_walls.npz          aS_ru_hub and aS_ru_shroud, dense parametric samples
                             (keys hub_S/hub_U/hub_V, shroud_S/...), so hub and
                             shroud are as exact as the blade (the msh facets sit
                             up to ~0.35 mm inside the cylinders).
"""
import json, sys, time
import numpy as np
from dtOOPythonSWIG import *

params = json.load(open("params.json"))
logMe.initLog("extract.log")
parser = dtXmlParser.init("machine.xml", "machineSave.xml").reference()
parser.parse()
cV = labeledVectorHandlingConstValue()
parser.createConstValue(cV)
parser.loadStateToConst("templateState", cV)
for k, v in params.items():
    cV.get(k).setValue(float(v))


def build():
    bC = baseContainer(); aF = labeledVectorHandlingAnalyticFunction()
    aG = labeledVectorHandlingAnalyticGeometry(); bV = labeledVectorHandlingBoundedVolume()
    dC = labeledVectorHandlingDtCase(); dP = labeledVectorHandlingDtPlugin()
    parser.createAnalyticFunction(bC, cV, aF)
    parser.createAnalyticGeometry(bC, cV, aF, aG)
    parser.createBoundedVolume(bC, cV, aF, aG, bV)
    parser.createPlugin(bC, cV, aF, aG, bV, dC, dP)
    return bC, aF, aG, bV, dC, dP


t0 = time.time()
c1 = build()
c1[5].get("ru_adjustDomain").apply()
bC, aF, aG, bV, dC, dP = build()
print(f"[time] build {time.time()-t0:.1f}s", flush=True)


def pt(g, *uv):
    p = g.getPoint(vectorReal(list(uv)))
    return (p.x(), p.y(), p.z())


g = aG.get("aS_ru_blade_rot")
umin, umax = g.getMin(0), g.getMax(0)
vmin, vmax = g.getMin(1), g.getMax(1)
print("range", umin, umax, vmin, vmax, flush=True)

NU, NV = 591, 81
U = np.linspace(umin, umax, NU)
V = np.linspace(vmin, vmax, NV)
S = np.empty((NU, NV, 3))
for i, u in enumerate(U):
    for j, v in enumerate(V):
        S[i, j] = pt(g, u, v)
print("sampled surface", S.shape, flush=True)

curves = {}
for name in ["aS_ru_blade_rot_v0.0", "aS_ru_blade_rot_v0.5", "aS_ru_blade_rot_v1.0"]:
    cg = aG.get(name)
    lo, hi = cg.getMin(0), cg.getMax(0)
    cu = np.linspace(lo, hi, 400)
    cs = []
    for u in cu:
        try:
            cs.append(pt(cg, u))
        except Exception:
            cs.append(pt(cg, u, 0.0))
    curves[name] = dict(U=cu, S=np.array(cs))
    print("curve", name, np.array(cs).shape, flush=True)

np.savez("machine_blade.npz", U=U, V=V, S=S,
         umin=umin, umax=umax, vmin=vmin, vmax=vmax, **curves)

walls = {}
for name, key in (("aS_ru_hub", "hub"), ("aS_ru_shroud", "shroud")):
    wg = aG.get(name)
    lo0, hi0, lo1, hi1 = wg.getMin(0), wg.getMax(0), wg.getMin(1), wg.getMax(1)
    NWU, NWV = 1441, 161          # 0.25 deg: chord error <= 5e-6 at r = 1.9
    WU, WV = np.linspace(lo0, hi0, NWU), np.linspace(lo1, hi1, NWV)
    WS = np.empty((NWU, NWV, 3))
    for i, u in enumerate(WU):
        for j, v in enumerate(WV):
            WS[i, j] = pt(wg, u, v)
    walls.update({f"{key}_S": WS, f"{key}_U": WU, f"{key}_V": WV})
    print("sampled", name, WS.shape, "range", lo0, hi0, lo1, hi1, flush=True)
np.savez("machine_walls.npz", **walls)

# legacy VTK structured grid
out = ["# vtk DataFile Version 3.0",
       "dtOO machine blade wall aS_ru_blade_rot (rotated)", "ASCII",
       "DATASET STRUCTURED_GRID", f"DIMENSIONS {NU} {NV} 1",
       f"POINTS {NU*NV} double"]
for i in range(NU):
    for j in range(NV):
        x, y, z = S[i, j]
        out.append(f"{x:.12g} {y:.12g} {z:.12g}")
open("machine_blade.vtk", "w").write("\n".join(out) + "\n")
print("wrote machine_blade.vtk", time.time() - t0, flush=True)
