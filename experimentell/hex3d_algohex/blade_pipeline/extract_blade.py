"""Export the true blade wall from the dtOO machine (aS_ru_blade_rot).

Writes machine_blade.npz (dense parametric samples + span cut curves) and
machine_blade.vtk (structured surface) into blade/.
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
