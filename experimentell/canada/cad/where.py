import sys
sys.argv=['x','0.0354']; exec(open('/work/canada/cad_envelope.py').read().split("for name, idx in RINGS.items():")[0])
from OCC.Core.STEPControl import STEPControl_Reader
from OCC.Core.BOPAlgo import BOPAlgo_ArgumentAnalyzer
for name in ("ring2","ring3"):
    r=STEPControl_Reader(); r.ReadFile(f"/work/canada/cad/{name}_envelope.stp"); r.TransferRoots(); s=r.OneShape()
    an=BOPAlgo_ArgumentAnalyzer(); an.SetShape1(s); an.SetSelfInterMode(True); an.Perform()
    res=an.GetCheckResult()
    print(name,"faulty entries:",res.Size())
    for it in res:
        st=it.GetCheckStatus(); fs=it.GetFaultyShapes1()
        for sh in fs:
            x0,y0,z0,x1,y1,z1=bbox(sh); print(f"   status {st}: shape type {sh.ShapeType()} bbox z[{z0:.2f},{z1:.2f}] x[{x0:.2f},{x1:.2f}] y[{y0:.2f},{y1:.2f}]")
