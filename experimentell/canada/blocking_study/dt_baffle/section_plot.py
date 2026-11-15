import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from OCC.Core.STEPControl import STEPControl_Reader
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_FACE, TopAbs_EDGE
from OCC.Core.TopoDS import topods
from OCC.Core.BRep import BRep_Tool
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Section
from OCC.Core.gp import gp_Pln, gp_Pnt, gp_Dir
r = STEPControl_Reader(); r.ReadFile("/work/canada/case/cad_blk.stp"); r.TransferRoots(); s = r.OneShape()
fs=[]; ex = TopExp_Explorer(s, TopAbs_FACE)
while ex.More(): fs.append(topods.Face(ex.Current())); ex.Next()
fig, axs = plt.subplots(2,1, figsize=(16,12))
for ax, phi in zip(axs, (0, 90)):
    d = gp_Dir(-np.sin(np.radians(phi)), np.cos(np.radians(phi)), 0)
    for i,f in enumerate(fs):
        sec = BRepAlgoAPI_Section(f, gp_Pln(gp_Pnt(0,0,0), d)); sec.Build()
        e = TopExp_Explorer(sec.Shape(), TopAbs_EDGE)
        while e.More():
            c,a,b = BRep_Tool.Curve(topods.Edge(e.Current()))
            P = np.array([[*c.Value(t).Coord()] for t in np.linspace(a,b,60)])
            rr = np.hypot(P[:,0],P[:,1])
            if (P[:,0]*np.cos(np.radians(phi))+P[:,1]*np.sin(np.radians(phi))).mean() > -1e-9:
                ax.plot(P[:,2], rr, '-', lw=1); ax.text(P[30,2], rr[30], str(i), fontsize=8)
            e.Next()
    ax.set_aspect('equal'); ax.set_title(f"phi={phi}"); ax.set_xlim(-2.5,7); ax.set_ylim(0,7); ax.grid()
plt.savefig("cad_sections.png", dpi=80)
