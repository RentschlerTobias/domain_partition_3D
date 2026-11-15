"""Stage B (gmsh OCC): assemble the draft-tube quarter domain with the thin walls as interior zero-thickness
baffles, tet-mesh it and write the AlgoHex input VTK (vertex = feature vertices, line = feature edges,
triangle = boundary + baffle faces with a surface label, tetra).

Steps: mid-surfaces cut to the quarter -> lip mid-surface fragmented with the two interface faces (the stub
upstream of the gv interface and the interface parts outside the lip are dropped) -> fragment of the outer
body with the hub-chain faces + all mid-surfaces -> the hub-chain volume (bulb + gv/ru hole + dt hub) is
removed; the parts of the mid-surfaces inside the domain stay as embedded interior faces.
Feature edges: gmsh curves with a label change, a crease > 25 deg, a free side (baffle LE/TE) or > 2 faces.

usage: python3 assemble_mesh.py <outdir> [--h-baffle 0.06] [--h-max 0.4]
"""
import argparse, json, os, collections
import numpy as np
import gmsh

ap = argparse.ArgumentParser()
ap.add_argument("out")
ap.add_argument("--h-baffle", type=float, default=0.06)
ap.add_argument("--h-max", type=float, default=0.4)
ap.add_argument("--dist-max", type=float, default=1.5)
ap.add_argument("--threads", type=int, default=8)
ap.add_argument("--geo-only", action="store_true")
a = ap.parse_args()
PARTS = "/work/canada/blocking_study/dt_baffle/parts"
os.makedirs(a.out, exist_ok=True)
info = json.load(open(f"{PARTS}/parts.json"))

gmsh.initialize(); gmsh.option.setNumber("General.Verbosity", 2)
gmsh.option.setNumber("General.NumThreads", a.threads)
gmsh.option.setNumber("Geometry.ToleranceBoolean", 1e-6)
occ = gmsh.model.occ
imp = lambda n: occ.importShapes(f"{PARTS}/{n}.brep", highestDimOnly=True)
bb = lambda d, t: np.array(gmsh.model.getBoundingBox(d, t))

body = imp("body"); assert len(body) == 1 and body[0][0] == 3, body
quarter_box = lambda: occ.addBox(0, 0, -100, 100, 100, 200)
mids = {}
for n, m in info["mid"].items():
    s = imp(n)
    q, _ = occ.intersect(s, [(3, quarter_box())])          # mid-surface cut to x >= 0, y >= 0
    mids[n] = q
occ.synchronize()

# ---- lip: trim with the interface faces
d_gv, d_dt = imp("if_gv_ext"), imp("if_dt_ext")
z_gv, z_dt, r_lip = info["lip_z_gv"], info["lip_z_dt"], info["lip_mid_r_at_dt_interface"]
out, omap = occ.fragment(mids["mid_lip"], d_gv + d_dt)
occ.synchronize()
lip_pieces, keep_if = [], {}
drop = []
for t in omap[0] if len(omap[0]) else []:
    pass
for k, src in enumerate(["lip"] * len(mids["mid_lip"]) + ["if_gv"] * len(d_gv) + ["if_dt"] * len(d_dt)):
    for d, t in omap[k]:
        if src == "lip":
            zc = occ.getCenterOfMass(d, t)[2]
            if zc < z_gv + 0.02:
                drop.append((d, t)); print(f"   drop lip stub {t}: centroid z {zc:.3f}")
            elif zc < z_dt:
                keep_if[t] = "LIP_WALL"
            else:
                keep_if[t] = "BAFFLE_LIP"
        else:
            area = occ.getMass(d, t); rim_area = np.pi / 4 * (info["r_ext"] ** 2 - r_lip ** 2)
            if abs(area - rim_area) < 0.1 * rim_area:
                drop.append((d, t)); print(f"   drop {src} rim {t}: area {area:.3f} (rim {rim_area:.3f})")
            else:
                keep_if[t] = "IN_OUTLET" if src == "if_gv" else "DT_INLET"
occ.remove(list(set(drop)), recursive=True)
occ.synchronize()
print("[assemble] lip/interface pieces:", keep_if)

hub = {}
for k, lab in (("bulb", "IN_HUB"), ("hub0", "DT_HUB"), ("hub1", "DT_HUB"), ("hubEnd", "DT_HUB")):
    for d, t in imp(k):
        hub[t] = lab
tools = [(2, t) for t in keep_if] + [(2, t) for t in hub] + [x for n, v in mids.items() if n != "mid_lip" for x in v]
src_lab = dict(keep_if); src_lab.update(hub)
for n, v in mids.items():
    if n != "mid_lip":
        for _, t in v:
            src_lab[t] = "BAFFLE_" + n.upper()
occ.synchronize()
out, omap = occ.fragment(body, tools)
occ.synchronize()
lab_of = {}
for (d, t), res in zip(body + tools, omap):
    if d == 2:
        for _, r in res:
            lab_of[r] = src_lab[t]
vols = gmsh.model.getEntities(3)
print("[assemble] volumes after fragment:", [(v, bb(3, v).round(2).tolist()) for _, v in vols])
# domain = the volume touching the far field (largest bbox)
dom = max(vols, key=lambda v: np.prod(bb(3, v[1])[3:] - bb(3, v[1])[:3]))
for v in vols:
    if v != dom:
        occ.remove([v])
occ.synchronize()
# free surfaces (outside the domain) are removed
for _, s in gmsh.model.getEntities(2):
    print("   surface", s, "adj", gmsh.model.getAdjacencies(2, s)[0], "src", lab_of.get(s), "com", np.round(occ.getCenterOfMass(2, s), 3).tolist(), flush=True)
# free surfaces that are not baffles: symmetry-plane pieces that bounded the removed hub-chain volume
free = [(2, s) for _, s in gmsh.model.getEntities(2) if len(gmsh.model.getAdjacencies(2, s)[0]) == 0 and s not in lab_of]
if free:
    print("[assemble] removing free surfaces:", free)
    occ.remove(free, recursive=True); occ.synchronize()
labels = {}
for _, s in gmsh.model.getEntities(2):
    if s in lab_of:
        labels[s] = lab_of[s]; continue
    c = occ.getCenterOfMass(2, s)
    if abs(c[1]) < 1e-7: labels[s] = "SYMNY"
    elif abs(c[0]) < 1e-7: labels[s] = "SYMNX"
    elif abs(c[2] + 4.0) < 1e-7: labels[s] = "IN_INLET"
    elif abs(c[2] - 10.0) < 1e-7: labels[s] = "OUT_OUTLET"
    else: labels[s] = "FAR"
cnt = collections.Counter(labels.values())
print("[assemble] surfaces per label:", dict(cnt))
names = sorted(set(labels.values()))
lab_id = {n: i + 1 for i, n in enumerate(names)}
baffle_ids = sorted(lab_id[n] for n in names if n.startswith("BAFFLE"))
gmsh.write(f"{a.out}/domain.brep")
if a.geo_only:
    raise SystemExit(0)

# ---- mesh
baf_s = [s for s, l in labels.items() if l.startswith("BAFFLE") or l == "LIP_WALL"]
f_d = gmsh.model.mesh.field.add("Distance"); gmsh.model.mesh.field.setNumbers(f_d, "SurfacesList", baf_s)
gmsh.model.mesh.field.setNumber(f_d, "Sampling", 100)
f_t = gmsh.model.mesh.field.add("Threshold")
gmsh.model.mesh.field.setNumber(f_t, "InField", f_d)
gmsh.model.mesh.field.setNumber(f_t, "SizeMin", a.h_baffle); gmsh.model.mesh.field.setNumber(f_t, "SizeMax", a.h_max)
gmsh.model.mesh.field.setNumber(f_t, "DistMin", 0.05); gmsh.model.mesh.field.setNumber(f_t, "DistMax", a.dist_max)
gmsh.model.mesh.field.setAsBackgroundMesh(f_t)
gmsh.option.setNumber("Mesh.MeshSizeMax", a.h_max)
gmsh.option.setNumber("Mesh.MeshSizeMin", 0.025)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 16)
gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.Algorithm", 6); gmsh.option.setNumber("Mesh.Algorithm3D", 10)
gmsh.model.mesh.generate(3)
gmsh.model.mesh.optimize("")
nt, xyz, _ = gmsh.model.mesh.getNodes()
P = xyz.reshape(-1, 3); idx = {int(k): i for i, k in enumerate(nt)}
NN = {1: 2, 2: 3, 4: 4}
def elems(dim, tag, typ):
    ty, _, nn = gmsh.model.mesh.getElements(dim, tag)
    out = [np.array([idx[int(v)] for v in b]).reshape(-1, NN[typ]) for t_, b in zip(ty, nn) if t_ == typ]
    return np.vstack(out) if out else np.zeros((0, NN[typ]), int)
tets = np.vstack([elems(3, v, 4) for _, v in gmsh.model.getEntities(3)])
tris, tlab, tsurf = [], [], []
for _, s in gmsh.model.getEntities(2):
    T = elems(2, s, 2); tris.append(T); tlab += [lab_id[labels[s]]] * len(T); tsurf += [s] * len(T)
tris = np.vstack(tris); tlab = np.array(tlab); tsurf = np.array(tsurf)
used = np.unique(np.r_[tets.ravel(), tris.ravel()])
nrm = np.cross(P[tris[:, 1]] - P[tris[:, 0]], P[tris[:, 2]] - P[tris[:, 0]]); nrm /= np.linalg.norm(nrm, axis=1)[:, None]
e2t = collections.defaultdict(list)
for i, tr in enumerate(tris):
    for x, y in ((tr[0], tr[1]), (tr[1], tr[2]), (tr[2], tr[0])):
        e2t[(min(x, y), max(x, y))].append(i)
edges, why = [], collections.Counter()
for _, c in gmsh.model.getEntities(1):
    E = elems(1, c, 1)
    if not len(E):
        continue
    labs, surfs, ang, nadj = set(), set(), 0.0, []
    for x, y in E:
        ts = e2t.get((min(x, y), max(x, y)), [])
        nadj.append(len(ts)); labs |= {int(tlab[t]) for t in ts}; surfs |= {int(tsurf[t]) for t in ts}
        if len(ts) == 2:
            ca = np.degrees(np.arccos(np.clip(nrm[ts[0]] @ nrm[ts[1]], -1, 1)))
            ang = max(ang, min(ca, 180 - ca))
    reason = ("label" if len(labs) > 1 else "free" if min(nadj) == 1 else "nonmanifold" if max(nadj) > 2
              else "crease" if ang > 25 else None)
    if reason:
        edges.append(E); why[reason] += 1
    else:
        why["dropped"] += 1
edges = np.vstack(edges)
val = collections.Counter(int(x) for e in edges for x in e)
inc = collections.defaultdict(list)
for x, y in edges:
    inc[int(x)].append(int(y)); inc[int(y)].append(int(x))
def turn(v):
    p, q = inc[v]; d1 = P[p] - P[v]; d2 = P[q] - P[v]
    return 180 - np.degrees(np.arccos(np.clip(d1 @ d2 / np.linalg.norm(d1) / np.linalg.norm(d2), -1, 1)))
fv = sorted({v for v, k in val.items() if k != 2} | {v for v, k in val.items() if k == 2 and turn(v) > 30})
print(f"[assemble] curves: {dict(why)}; feature edges {len(edges)}, feature vertices {len(fv)}")
print(f"[assemble] {len(used)} nodes, {len(tets)} tets, {len(tris)} triangles "
      f"({int(np.isin(tlab, baffle_ids).sum())} baffle)")
gmsh.option.setNumber("Mesh.SaveAll", 1); gmsh.write(f"{a.out}/tet.msh")
gmsh.finalize()
np.savez(f"{a.out}/tet_arrays.npz", P=P, tets=tets, tris=tris, tlab=tlab, edges=edges, fv=np.array(fv))
json.dump(dict(labels=lab_id, baffle_labels=baffle_ids, h_baffle=a.h_baffle, h_max=a.h_max,
               n_nodes=int(len(P)), n_tets=int(len(tets)), n_tris=int(len(tris)), n_feature_edges=int(len(edges)),
               n_feature_vertices=len(fv), curves=dict(why), surfaces_per_label=dict(cnt)),
          open(f"{a.out}/input_info.json", "w"), indent=1)
with open(f"{a.out}/tet.vtk", "w") as f:
    cells = [(1, [v], 1) for v in fv] + [(3, list(e), 1) for e in edges] + \
            [(5, list(tr), int(l)) for tr, l in zip(tris, tlab)] + [(10, list(tt), 0) for tt in tets]
    f.write("# vtk DataFile Version 2.0\ndt_baffle\nASCII\nDATASET UNSTRUCTURED_GRID\n")
    f.write(f"POINTS {len(P)} double\n" + "".join(f"{p[0]:.15g} {p[1]:.15g} {p[2]:.15g}\n" for p in P))
    f.write(f"CELLS {len(cells)} {sum(1 + len(c[1]) for c in cells)}\n")
    f.write("".join(f"{len(c[1])} " + " ".join(map(str, c[1])) + "\n" for c in cells))
    f.write(f"CELL_TYPES {len(cells)}\n" + "".join(f"{c[0]}\n" for c in cells))
    f.write(f"CELL_DATA {len(cells)}\nSCALARS color int 1\nLOOKUP_TABLE default\n" + "".join(f"{c[2]}\n" for c in cells))
with open(f"{a.out}/labels.txt", "w") as f:
    f.write("".join(f"{i} {n}\n" for n, i in sorted(lab_id.items(), key=lambda x: x[1])))
print("[assemble] wrote", f"{a.out}/tet.vtk")
