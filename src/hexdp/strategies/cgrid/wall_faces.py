"""Real thin walls (rings, lip) of the draft-tube CAD (case/cad_blk.stp), detected generically, and their meridional
sections at every phi of ../sections/.

1. thin-wall skin pairs: same criterion as dt_baffle/build_parts.py (two curved faces whose mutual distance is
   small and nearly constant compared with their size).
2. wall component of a pair = transitive closure of faces that touch it (distance < 1e-4), excluding faces that
   are part of the baffle-domain boundary (faces lying on dt_baffle/in_a/domain.brep: far field, hub chain,
   interfaces). This collects bevels / blunt edges and, for the lip, the mantle and the lip front.
3. geometric replacement of the fixed index faces[11] (mantle) of build_parts.py: the mantle is the face of a
   wall component that touches the upstream (min z) rim of the skin pair and is not one of the pair's skins.
4. per phi: section of each component with the half plane -> polylines in (r,z), chained (closed for free rings).
Writes walls.json (components) and walls/wall_<phi>.json.
run: LD_LIBRARY_PATH=/dtOO-install/lib PYTHONPATH=<repo>/src /opt/venv-dtoo/bin/python wall_faces.py \
       [--step cad_blk.stp] [--brep domain.brep] [--out DIR] [--sections DIR]
     (defaults: out = the working directory, sections = out/../sections)
"""
import sys, os, json, math, glob
import numpy as np
from OCC.Core.STEPControl import STEPControl_Reader
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_FACE, TopAbs_EDGE
from OCC.Core.TopoDS import topods
from OCC.Core.BRep import BRep_Tool
from OCC.Core.BRepTools import breptools
from OCC.Core.BRepAdaptor import BRepAdaptor_Surface, BRepAdaptor_Curve
from OCC.Core.GeomAbs import GeomAbs_Plane
from OCC.Core.BRepExtrema import BRepExtrema_DistShapeShape
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeVertex, BRepBuilderAPI_MakeFace
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Section
from OCC.Core.gp import gp_Pnt, gp_Dir, gp_Ax3, gp_Pln
from OCC.Core.Bnd import Bnd_Box
from OCC.Core.BRepBndLib import brepbndlib
from hexdp.strategies.sweep import section_rz



def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description="thin walls of the CAD and their meridional sections")
    ap.add_argument("--step", default=os.environ.get("HEXDP_CAD_STEP", "cad_blk.stp"))
    ap.add_argument("--brep", default=section_rz.SRC, help="domain solid (for the baffle-domain faces)")
    ap.add_argument("--out", default=os.getcwd(), help="writes walls.json and walls/ here")
    ap.add_argument("--sections", default=None, help="default: OUT/../sections")
    a = ap.parse_args(argv)
    HERE = a.out
    SECTIONS = a.sections or os.path.join(os.path.dirname(os.path.abspath(HERE)), "sections")
    section_rz.SRC = a.brep
    if not os.path.isfile(a.step):                      # OCC crashes (segfault) on a missing file
        sys.exit(f"[walls] no such STEP file: {a.step}")
    r = STEPControl_Reader(); r.ReadFile(a.step); r.TransferRoots(); shape = r.OneShape()
    faces = []; ex = TopExp_Explorer(shape, TopAbs_FACE)
    while ex.More():
        faces.append(topods.Face(ex.Current())); ex.Next()
    _, dom_faces, _ = section_rz.load()


    def bbox(s):
        b = Bnd_Box(); brepbndlib.Add(s, b); return np.array(b.Get())


    def dist(a, b):
        d = BRepExtrema_DistShapeShape(a, b); d.Perform(); return d.Value()


    def dist_pt(p, s):
        return dist(BRepBuilderAPI_MakeVertex(gp_Pnt(*p)).Vertex(), s)


    def face_samples(f, n=7):
        u0, u1, v0, v1 = breptools.UVBounds(f); s = BRep_Tool.Surface(f)
        return np.array([[*s.Value(u, v).Coord()] for u in np.linspace(u0, u1, n + 2)[1:-1] for v in np.linspace(v0, v1, n + 2)[1:-1]])


    # ---- 1. pairs (as build_parts.py)
    curved = [i for i, f in enumerate(faces) if BRepAdaptor_Surface(f).GetType() != GeomAbs_Plane]
    S = {i: face_samples(faces[i]) for i in curved}; B = {i: bbox(faces[i]) for i in curved}
    diag = {i: np.linalg.norm(B[i][3:] - B[i][:3]) for i in curved}
    pairs = []
    for ia, a in enumerate(curved):
        for b in curved[ia + 1:]:
            lo = np.maximum(B[a][:3], B[b][:3]); hi = np.minimum(B[a][3:], B[b][3:])
            if np.any(lo > hi + 0.2):
                continue
            d = np.r_[[dist_pt(p, faces[b]) for p in S[a]], [dist_pt(p, faces[a]) for p in S[b]]]
            size = min(diag[a], diag[b])
            if 1e-6 < np.median(d) < 0.05 * size and np.percentile(d, 95) < 0.1 * size:
                pairs.append((a, b, float(np.median(d)), float(d.max())))
    print("[walls] thin-wall pairs:", [(a, b, round(m, 4)) for a, b, m, _ in pairs], flush=True)

    # ---- 2. domain-boundary faces (lie on the baffle domain) and touching closure
    on_domain = set()
    for i, f in enumerate(faces):
        P = face_samples(f, 5); P = P[(P[:, 0] >= -1e-9) & (P[:, 1] >= -1e-9)]    # samples inside the quarter
        if len(P) and np.mean([min(dist_pt(p, g) for g in dom_faces) < 1e-5 for p in P]) > 0.5:
            on_domain.add(i)
    print("[walls] CAD faces on the baffle-domain boundary:", sorted(on_domain), flush=True)
    # faces lying on the (untrimmed) surface of a domain-boundary face continue that boundary (e.g. an interface
    # annulus inside the lip) and are not wall faces
    from OCC.Core.GeomAPI import GeomAPI_ProjectPointOnSurf
    from OCC.Core.gp import gp_Pnt as _P
    def on_surface_ext(i):
        P = face_samples(faces[i], 3)
        for j in on_domain:
            srf = BRep_Tool.Surface(faces[j])
            if all(GeomAPI_ProjectPointOnSurf(_P(*p), srf).NbPoints() and
                   GeomAPI_ProjectPointOnSurf(_P(*p), srf).LowerDistance() < 1e-4 for p in P):
                return j
        return None
    ext = {i: on_surface_ext(i) for i in range(len(faces)) if i not in on_domain}
    ext = {i: j for i, j in ext.items() if j is not None}
    print("[walls] faces on the extension of a domain-boundary face (excluded):", ext, flush=True)
    cand = [i for i in range(len(faces)) if i not in on_domain and i not in ext]
    touch = {i: {j for j in cand if j != i and np.all(bbox(faces[i])[:3] - 1e-3 <= bbox(faces[j])[3:]) and
                 np.all(bbox(faces[j])[:3] - 1e-3 <= bbox(faces[i])[3:]) and dist(faces[i], faces[j]) < 1e-4} for i in cand}
    comps = []
    for a, b, gap, gmax in pairs:
        comp = {a, b}; stack = [a, b]
        while stack:
            x = stack.pop()
            for y in touch[x]:
                if y not in comp:
                    comp.add(y); stack.append(y)
        # mantle (replacement for build_parts' fixed faces[11]): face of the component touching the upstream rim
        # (minimum z) of the pair that is not one of the two skins
        zmin = min(bbox(faces[a])[2], bbox(faces[b])[2])
        up = [j for j in comp - {a, b} if bbox(faces[j])[5] <= zmin + 1e-3]       # lies upstream of the pair
        up_touch = [j for j in up if (dist(faces[j], faces[a]) < 1e-4 or dist(faces[j], faces[b]) < 1e-4)
                    and bbox(faces[j])[5] - bbox(faces[j])[2] > 5 * gap]                  # not a bevel / edge face
        mantle = max(up_touch, key=lambda j: bbox(faces[j])[5] - bbox(faces[j])[2]) if up_touch else None
        comps.append(dict(pair=[a, b], gap=gap, max_gap=gmax, faces=sorted(comp), mantle=mantle))
        print(f"[walls] pair ({a},{b}) gap {gap:.4f}: component {sorted(comp)}, mantle {mantle}", flush=True)
    json.dump(dict(pairs=pairs, components=comps, on_domain=sorted(on_domain)), open(f"{HERE}/walls.json", "w"), indent=1)


    # ---- 4. sections
    def chain(polys, tol=5e-4):
        polys = [np.asarray(p) for p in polys]
        out = [polys.pop(0)]
        while polys:
            head, tail = out[0][0], out[-1][-1]; hit = False
            for k, p in enumerate(polys):
                for P, at_tail in ((p, True), (p[::-1], True), (p, False), (p[::-1], False)):
                    if at_tail and np.linalg.norm(P[0] - tail) < tol:
                        out.append(P); polys.pop(k); hit = True; break
                    if not at_tail and np.linalg.norm(P[-1] - head) < tol:
                        out.insert(0, P); polys.pop(k); hit = True; break
                if hit:
                    break
            if not hit:
                break
        dropped = [(np.round(p[0], 4).tolist(), np.round(p[-1], 4).tolist()) for p in polys]
        pts = np.vstack([out[0]] + [p[1:] for p in out[1:]])
        closed = np.linalg.norm(pts[0] - pts[-1]) < tol
        return pts, closed, dropped


    os.makedirs(f"{HERE}/walls", exist_ok=True)
    for fn in sorted(glob.glob(f"{SECTIONS}/sec_*.json")):
        phi = json.load(open(fn))["phi"]; p = math.radians(min(max(phi, 0.25), 90 - 0.25))   # 0.25 deg off the symmetry planes (face seams; the wall is symmetric there, error 2nd order)
        er = np.array([math.cos(p), math.sin(p), 0.0])
        ax = gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(-math.sin(p), math.cos(p), 0), gp_Dir(*er))
        pl = BRepBuilderAPI_MakeFace(gp_Pln(ax), 0.0, 30.0, -20.0, 30.0).Face()
        out = []
        for c in comps:
            polys = []
            for i in c["faces"]:
                sc = BRepAlgoAPI_Section(faces[i], pl); sc.Build()
                e = TopExp_Explorer(sc.Shape(), TopAbs_EDGE)
                while e.More():
                    Q = section_rz.sample(e.Current(), ds=0.002); polys.append(np.c_[Q @ er, Q[:, 2]]); e.Next()
            try:
                # start the chain with a piece of the first skin of the pair
                pts, closed, dropped = chain(polys)
            except RuntimeError as err:
                print(f"   phi {phi} pair {c['pair']}: {err}; per-face pieces:", flush=True)
                raise
            out.append(dict(pair=c["pair"], closed=bool(closed), pts=pts.tolist(), dropped=dropped))
        json.dump(dict(phi=phi, walls=out), open(f"{HERE}/walls/wall_{phi:05.1f}.json", "w"))
        print(f"[walls] phi {phi}: " + ", ".join(f"{w['pair']} {'closed' if w['closed'] else 'open'} {len(w['pts'])} pts dropped {w['dropped']}" for w in out), flush=True)


if __name__ == "__main__":
    main()
