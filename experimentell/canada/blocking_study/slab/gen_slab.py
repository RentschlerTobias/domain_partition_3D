"""Synthetic 'diffuser slab' for the blocking study: box (axial x, radial z, sweep y) with two staggered
thin rings (plates spanning the full sweep width, like the canada diffuser rings between the two
symmetry planes).  Variants:
  walls  : plates as holes with real thickness t, blunt LE, 8 deg knife TE
  env    : plates replaced by their delta-envelope (Minkowski sum with a disk), as holes
  baffle : plates replaced by their zero-thickness mid surface, embedded in the volume (internal feature faces)
Writes the AlgoHex VTK (vertex/line/triangle/tet cells, color = feature tag / surface label)."""
import sys, math, numpy as np, gmsh
var = sys.argv[1]; out = sys.argv[2]
L, Wy, Hz = 4.0, 1.0, 2.0
t, g, delta, te_deg = 0.065, 0.13, 0.035, 8.0
hmin, hmax = (float(sys.argv[3]) if len(sys.argv) > 3 else 0.03), 0.15
def plate(x0, z0, chord, slope_deg):
    a = math.radians(slope_deg); c, s = math.cos(a), math.sin(a)
    d = np.array([c, s]); n = np.array([-s, c]); p0 = np.array([x0, z0])
    lt = (t / 2) / math.tan(math.radians(te_deg / 2))      # taper length for the knife TE
    le_u, le_l = p0 + n * t / 2, p0 - n * t / 2
    tp_u, tp_l = p0 + d * (chord - lt) + n * t / 2, p0 + d * (chord - lt) - n * t / 2
    te = p0 + d * chord
    poly = [le_l, tp_l, te, tp_u, le_u]                    # convex, counter-clockwise
    mid = (p0, te)
    return poly, mid
A = plate(1.0, 0.80, 1.6, 8.0)
zA = lambda x: 0.80 + (x - 1.0) * math.tan(math.radians(8.0))
B = plate(2.2, zA(2.2) + g + t, 1.5, 12.0)
plates = [A, B]
gmsh.initialize(); gmsh.option.setNumber("General.Verbosity", 2)
occ = gmsh.model.occ
box = occ.addBox(0, 0, 0, L, Wy, Hz)
def face2d(poly):
    pts = [occ.addPoint(p[0], 0, p[1]) for p in poly]
    ls = [occ.addLine(pts[i], pts[(i + 1) % len(pts)]) for i in range(len(pts))]
    return occ.addPlaneSurface([occ.addCurveLoop(ls)])
tools = []
if var in ("walls", "env"):
    for poly, _ in plates:
        fs = [face2d(poly)]
        if var == "env":
            P = np.array(poly)
            for i in range(len(P)):
                a, b = P[i], P[(i + 1) % len(P)]
                e = b - a; nrm = np.array([e[1], -e[0]]) / np.linalg.norm(e)   # outward (ccw polygon)
                fs.append(face2d([a, b, b + nrm * delta, a + nrm * delta]))
                fs.append(occ.addDisk(a[0], 0, a[1], delta, delta, zAxis=[0, 1, 0]))
            f, _ = occ.fuse([(2, fs[0])], [(2, x) for x in fs[1:]]); occ.synchronize()
            fs = [x[1] for x in f]
        for f in fs:
            ex = occ.extrude([(2, f)], 0, Wy, 0)
            tools += [e for e in ex if e[0] == 3]
    occ.synchronize()
    if len(tools) > 1:
        tools, _ = occ.fuse(tools[:1], tools[1:])
    occ.cut([(3, box)], tools); occ.synchronize()
else:
    surfs = []
    for _, (p, q) in plates:
        a, b = occ.addPoint(p[0], 0, p[1]), occ.addPoint(q[0], 0, q[1])
        ex = occ.extrude([(1, occ.addLine(a, b))], 0, Wy, 0)
        surfs += [e for e in ex if e[0] == 2]
    occ.fragment([(3, box)], surfs); occ.synchronize()
vols = gmsh.model.getEntities(3); assert len(vols) == 1, vols
# mesh size: fine near the plates
fld = gmsh.model.mesh.field
pc = [c for c in gmsh.model.getEntities(1)]
# surfaces belonging to the plates: those not on the box planes
def on_box(tag):
    bb = gmsh.model.getBoundingBox(2, tag)
    return any(abs(bb[i] - bb[i + 3]) < 1e-5 and (abs(bb[i]) < 1e-5 or abs(bb[i] - [L, Wy, Hz][i]) < 1e-5) for i in range(3))
psurf = [s for _, s in gmsh.model.getEntities(2) if not on_box(s)]
d = fld.add("Distance"); fld.setNumbers(d, "SurfacesList", psurf); fld.setNumber(d, "Sampling", 60)
th = fld.add("Threshold"); fld.setNumber(th, "InField", d); fld.setNumber(th, "SizeMin", hmin)
fld.setNumber(th, "SizeMax", hmax); fld.setNumber(th, "DistMin", 0.02); fld.setNumber(th, "DistMax", 0.5)
fld.setAsBackgroundMesh(th)
gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.Algorithm3D", 10)
gmsh.model.mesh.generate(3)
nt, xyz, _ = gmsh.model.mesh.getNodes()
P = xyz.reshape(-1, 3); idx = {int(k): i for i, k in enumerate(nt)}
def elems(dim, tag, typ):
    ty, _, nn = gmsh.model.mesh.getElements(dim, tag)
    out = []
    for a, b in zip(ty, nn):
        if a == typ: out.append(np.array([idx[int(v)] for v in b]).reshape(-1, {1: 2, 2: 3, 4: 4, 15: 1}[typ]))
    return np.vstack(out) if out else np.zeros((0, {1: 2, 2: 3, 4: 4, 15: 1}[typ]), int)
tets = np.vstack([elems(3, v, 4) for _, v in vols])
def surf_label(s):
    bb = gmsh.model.getBoundingBox(2, s)
    for i in range(3):
        if abs(bb[i] - bb[i + 3]) < 1e-5:
            if abs(bb[i]) < 1e-5: return 1 + 2 * i
            if abs(bb[i] - [L, Wy, Hz][i]) < 1e-5: return 2 + 2 * i
    c = np.array([(bb[0] + bb[3]) / 2, (bb[2] + bb[5]) / 2])
    dA = min(np.linalg.norm(c - np.array(p)) for p in A[0]); dB = min(np.linalg.norm(c - np.array(p)) for p in B[0])
    return 7 if dA < dB else 8
tris, lab = [], []
for _, s in gmsh.model.getEntities(2):
    T = elems(2, s, 2); tris.append(T); lab += [surf_label(s)] * len(T)
tris = np.vstack(tris); lab = np.array(lab)
nrm = np.cross(P[tris[:, 1]] - P[tris[:, 0]], P[tris[:, 2]] - P[tris[:, 0]]); nrm /= np.linalg.norm(nrm, axis=1)[:, None]
e2t = {}
for i, tr in enumerate(tris):
    for a_, b_ in ((tr[0], tr[1]), (tr[1], tr[2]), (tr[2], tr[0])):
        e2t.setdefault((min(a_, b_), max(a_, b_)), []).append(i)
edges, kept, dropped = [], 0, 0
fvs = set()
for _, c in gmsh.model.getEntities(1):
    E = elems(1, c, 1)
    if not len(E): continue
    labs, ang, nadj = set(), 0.0, []
    for a_, b_ in E:
        ts = e2t.get((min(a_, b_), max(a_, b_)), [])
        nadj.append(len(ts)); labs |= {int(lab[x]) for x in ts}
        if len(ts) == 2:
            ang = max(ang, np.degrees(np.arccos(np.clip(abs(nrm[ts[0]] @ nrm[ts[1]]) if False else nrm[ts[0]] @ nrm[ts[1]], -1, 1))))
    # crease: unsigned normal deviation; consistent orientation is not guaranteed -> use min(ang, 180-ang)
    crease = min(ang, 180 - ang) > 25
    free = min(nadj) == 1                       # free edge of an internal surface
    if len(labs) > 1 or crease or free or max(nadj) > 2:
        edges.append(E); kept += 1; fvs |= {int(E[0, 0]), int(E[-1, 1])}
    else:
        dropped += 1
edges = np.vstack(edges)
from collections import Counter
val = Counter(int(x) for e in edges for x in e)
fv_end = {v_ for v_, k in val.items() if k != 2}
plate_lab = {7, 8}
inc = {}
for a_, b_ in edges:
    inc.setdefault(int(a_), []).append(int(b_)); inc.setdefault(int(b_), []).append(int(a_))
def turn(v_):
    a_, b_ = inc[v_]; d1 = P[a_] - P[v_]; d2 = P[b_] - P[v_]
    return 180 - np.degrees(np.arccos(np.clip(d1 @ d2 / np.linalg.norm(d1) / np.linalg.norm(d2), -1, 1)))
fv = sorted(fv_end | {v_ for v_ in fvs if val[v_] == 2 and turn(v_) > 30})
print(f'[slab] curves kept {kept}, smooth seams dropped {dropped}')
print(f"[slab] {var}: {len(P)} pts, {len(tets)} tets, {len(tris)} tris ({sum(np.isin(lab, list(plate_lab)))} on plates), "
      f"{len(edges)} feature edges, {len(fv)} feature vertices, plate surfaces {sorted(plate_lab)}")
v = P[tets]; vol = np.abs(np.einsum('ij,ij->i', np.cross(v[:, 1] - v[:, 0], v[:, 2] - v[:, 0]), v[:, 3] - v[:, 0])).sum() / 6
print(f"[slab] volume {vol:.4f}; uniform hex edge n=1000 {(vol/1000)**(1/3):.4f} n=2000 {(vol/2000)**(1/3):.4f} n=8000 {(vol/8000)**(1/3):.4f}")
gmsh.write(out.replace(".vtk", ".msh"))
gmsh.finalize()
with open(out, "w") as f:
    cells = [(1, [v_], 1) for v_ in fv] + [(3, list(e), 1) for e in edges] + \
            [(5, list(tr), int(l)) for tr, l in zip(tris, lab)] + [(10, list(tt), 0) for tt in tets]
    f.write("# vtk DataFile Version 2.0\nslab\nASCII\nDATASET UNSTRUCTURED_GRID\n")
    f.write(f"POINTS {len(P)} double\n" + "".join(f"{p[0]:.12g} {p[1]:.12g} {p[2]:.12g}\n" for p in P))
    f.write(f"CELLS {len(cells)} {sum(1 + len(c[1]) for c in cells)}\n")
    f.write("".join(f"{len(c[1])} " + " ".join(map(str, c[1])) + "\n" for c in cells))
    f.write(f"CELL_TYPES {len(cells)}\n" + "".join(f"{c[0]}\n" for c in cells))
    f.write(f"CELL_DATA {len(cells)}\nSCALARS color int 1\nLOOKUP_TABLE default\n" + "".join(f"{c[2]}\n" for c in cells))
