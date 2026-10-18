"""Cut a thin/sharp body out of a flow domain by a SMOOTH offset envelope.

Generic replacement for "cut the dtOO O-grid out": the envelope is the level
set {d(x) = delta} of the unsigned distance d to the body wall. An offset of a
surface is smooth wherever delta is below the concave curvature radius, so a
sharp or blunt trailing edge becomes a cylinder of radius delta and the rings
where the envelope meets the neighbouring cut surfaces have no corners. The
dtOO O-grid interface, by contrast, inherits the corners of a blunt TE
(canada: 2 x ~85 deg per ring, AlgoHex overflowed on them).

Pipeline (prototype, dtOO ru_gridGmsh layout):
  1. reduced domain + labels exactly as tet_prep_v5 (O-grid and BL prisms out)
  2. body wall = exterior faces of the removed O-grid that neither touch the
     remaining volume nor the real walls (caps at hub/shroud)
  3. tet-mesh the reduced domain, level set ls = d - delta on its nodes
  4. mmg3d -ls discretises {ls = 0} conformingly; keep the ls > 0 part
  5. relabel: the new iso surface becomes the O-grid interface label, all other
     labels survive through mmg's surface references
  6. feature graph from label boundaries, write the AlgoHex VTK

usage: envelope_cut.py <msh> <out_dir> [--delta-factor 1.25] [--size-max 0.06]
"""
import argparse, collections, os, subprocess, sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

HEX = Path(__file__).resolve().parent
sys.path.insert(0, str(HEX.parent.parent)); sys.path.insert(0, str(HEX))
from dp3d.extraction import parse_msh                       # noqa: E402
import tet_prep as tp                                       # noqa: E402
import tet_prep_v2 as v2                                    # noqa: E402
import tet_prep_v5 as v5                                    # noqa: E402

# mmg is not in the container: git clone https://github.com/MmgTools/mmg, cmake
# -DUSE_SCOTCH=OFF -DUSE_ELAS=OFF -DUSE_VTK=OFF, make install; point MMG3D at mmg3d_O3
MMG = os.environ.get("MMG3D", "/work/tools/mmg-install/bin/mmg3d_O3")
ISO_REF = 10                                                # mmg's MG_ISO


# ---------------------------------------------------------------- distance
def point_triangle_dist(p, a, b, c):
    """Exact distances of points p (n,3) to triangles (a,b,c) (n,3 each)
    (Ericson, Real-Time Collision Detection 5.1.5, vectorised)."""
    ab, ac, ap = b - a, c - a, p - a
    d1, d2 = (ab * ap).sum(1), (ac * ap).sum(1)
    bp = p - b
    d3, d4 = (ab * bp).sum(1), (ac * bp).sum(1)
    cp = p - c
    d5, d6 = (ab * cp).sum(1), (ac * cp).sum(1)
    va, vb, vc = d3 * d6 - d5 * d4, d5 * d2 - d1 * d6, d1 * d4 - d3 * d2
    den = np.where(np.abs(va + vb + vc) < 1e-300, 1e-300, va + vb + vc)
    v, w = vb / den, vc / den
    q = a + ab * v[:, None] + ac * w[:, None]                 # interior
    sd = lambda x, y: np.where(np.abs(y) < 1e-300, 0.0, x / np.where(np.abs(y) < 1e-300, 1, y))
    # regions, later assignments win (vertex > edge > face, as in Ericson)
    m = (vc <= 0) & (d1 >= 0) & (d3 <= 0)
    q[m] = (a + ab * sd(d1, d1 - d3)[:, None])[m]
    m = (vb <= 0) & (d2 >= 0) & (d6 <= 0)
    q[m] = (a + ac * sd(d2, d2 - d6)[:, None])[m]
    m = (va <= 0) & ((d4 - d3) >= 0) & ((d5 - d6) >= 0)
    q[m] = (b + (c - b) * sd(d4 - d3, (d4 - d3) + (d5 - d6))[:, None])[m]
    q[(d1 <= 0) & (d2 <= 0)] = a[(d1 <= 0) & (d2 <= 0)]
    q[(d3 >= 0) & (d4 <= d3)] = b[(d3 >= 0) & (d4 <= d3)]
    q[(d6 >= 0) & (d5 <= d6)] = c[(d6 >= 0) & (d5 <= d6)]
    return np.linalg.norm(p - q, axis=1)


def surface_distance(V, F, pts, k=24):
    """Unsigned distance of pts to the triangle mesh (V,F): exact over the k
    triangles with the nearest centroids (enough for a fine, regular wall)."""
    tree = cKDTree(V[F].mean(1))
    _, cand = tree.query(pts, k=k)
    best = np.full(len(pts), np.inf)
    for j in range(k):
        f = F[cand[:, j]]
        best = np.minimum(best, point_triangle_dist(pts, V[f[:, 0]], V[f[:, 1]], V[f[:, 2]]))
    return best


# ------------------------------------------------------------- body wall
def body_wall(nodes, els):
    """Exterior faces of the removed O-grid hexes that face neither the
    remaining volume nor a real wall: the blade surface itself."""
    v1f = set()
    for e in els:
        if e[0] in (4, 6, 7):
            for loc in tp._VOL_FACES[e[0]]:
                v1f.add(frozenset(e[3][i] for i in loc))
    cnt, fo = collections.Counter(), {}
    for e in els:
        if e[0] == 5 and e[2] in v5.OGRID_GEOM:
            for loc in tp._VOL_FACES[5]:
                f = tuple(e[3][i] for i in loc)
                cnt[frozenset(f)] += 1
                fo[frozenset(f)] = (f, e[3])
    look = {frozenset(e[3][:4]): e[2] for e in els if e[0] == 3}
    walln = {v for e in els if e[0] in (2, 3) and e[2] in (v5.GEOM_HUB, v5.GEOM_SHROUD)
             for v in e[3]}
    ent = collections.defaultdict(set)
    for k, c in cnt.items():
        if c == 1 and k not in v1f:
            ent[look.get(k)].add(k)
    body = [g for g, fs in ent.items() if not any(v in walln for f in fs for v in f)]
    bt = np.array(tp.orient_and_triangulate([fo[k] for g in body for k in ent[g]], nodes))
    u = np.unique(bt)
    m = {v: i for i, v in enumerate(u)}
    V = np.array([nodes[v] for v in u], float)
    F = np.vectorize(m.get)(bt)
    print(f"[env] body wall: entities {sorted(body)}, {len(F)} triangles")
    return V, F


# ------------------------------------------------------------- medit io
def write_mesh(path, P, tets, tris, refs):
    with open(path, "w") as f:
        f.write("MeshVersionFormatted 2\nDimension 3\n")
        f.write(f"Vertices\n{len(P)}\n")
        for p in P:
            f.write(f"{p[0]:.17g} {p[1]:.17g} {p[2]:.17g} 0\n")
        f.write(f"Triangles\n{len(tris)}\n")
        for t, r in zip(tris, refs):
            f.write(f"{t[0]+1} {t[1]+1} {t[2]+1} {r}\n")
        f.write(f"Tetrahedra\n{len(tets)}\n")
        for t in tets:
            f.write(f"{t[0]+1} {t[1]+1} {t[2]+1} {t[3]+1} 1\n")
        f.write("End\n")


def write_sol(path, vals):
    with open(path, "w") as f:
        f.write("MeshVersionFormatted 2\nDimension 3\n")
        f.write(f"SolAtVertices\n{len(vals)}\n1 1\n")
        f.writelines(f"{v:.17g}\n" for v in vals)
        f.write("End\n")


def read_mesh(path):
    tok = open(path).read().split()
    out, i = {}, 0
    while i < len(tok):
        key = tok[i]
        if key in ("Vertices", "Triangles", "Tetrahedra", "Edges", "Ridges",
                   "RequiredVertices", "Corners", "RequiredEdges", "RequiredTriangles",
                   "Normals", "NormalAtVertices", "Tangents", "TangentAtVertices", "Quadrilaterals"):
            n = int(tok[i + 1])
            w = {"Vertices": 4, "Triangles": 4, "Tetrahedra": 5, "Edges": 3, "Ridges": 1,
                 "RequiredVertices": 1, "Corners": 1, "RequiredEdges": 1, "RequiredTriangles": 1,
                 "Normals": 3, "NormalAtVertices": 2, "Tangents": 3, "TangentAtVertices": 2,
                 "Quadrilaterals": 5}[key]
            out[key] = np.array(tok[i + 2:i + 2 + n * w], float).reshape(n, w)
            i += 2 + n * w
        else:
            i += 1
    return out


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("msh")
    ap.add_argument("out_dir")
    ap.add_argument("--delta-factor", type=float, default=1.25,
                    help="delta = factor * max distance of the removed O-grid's outer face")
    ap.add_argument("--delta", type=float, default=None, help="absolute delta, overrides factor")
    ap.add_argument("--size-max", type=float, default=0.06)
    ap.add_argument("--hausd", type=float, default=0.002)
    a = ap.parse_args()
    out = Path(a.out_dir); out.mkdir(parents=True, exist_ok=True)

    nodes, els = parse_msh(Path(a.msh))
    tris = v5.reduced_boundary(nodes, els)
    ids = v5.classify(nodes, tris, v5.tagged_2d_lookup(els), v5.removed_face_kind(els))
    BV, BF = body_wall(nodes, els)

    T = np.array(tris)
    dist = {}
    for sid in np.unique(ids):
        pts = np.array([nodes[v] for v in np.unique(T[ids == sid])])
        d = surface_distance(BV, BF, pts)
        dist[sid] = (d.min(), d.max())
        print(f"[env]   {v5.NAMES[sid]:20s} distance to body [{d.min():.4f}, {d.max():.4f}]")
    delta = a.delta or a.delta_factor * dist[v5.S_OGRID_IFACE][1]
    clear = min(dist[s][0] for s in (v5.S_INLET, v5.S_OUTLET, v5.S_PER_A, v5.S_PER_B) if s in dist)
    print(f"[env] delta = {delta:.4f}; closest inlet/outlet/periodic = {clear:.4f}")
    if delta >= clear:
        raise SystemExit("[env] envelope would hit inlet/outlet/periodic; lower delta")

    # tet mesh of the reduced domain (O-grid hole still inside, it lies in ls < 0)
    P, tets, idx, used = v2.mesh_interior(nodes, tris, a.size_max)
    remap = {ot: idx[nt] for nt, ot in zip(range(1, len(used) + 1), used) if nt in idx}
    btris = [[remap[t] for t in tri] for tri in tris]
    ls = surface_distance(BV, BF, P) - delta
    write_mesh(out / "in.mesh", P, tets, btris, ids)
    write_sol(out / "in.sol", ls)
    cmd = [MMG, "-in", str(out / "in.mesh"), "-sol", str(out / "in.sol"), "-ls", "0",
           "-hmax", str(a.size_max), "-hausd", str(a.hausd), "-out", str(out / "cut.mesh")]
    print("[env]", " ".join(cmd), flush=True)
    with open(out / "mmg.log", "w") as lf:
        rc = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT).returncode
    if rc:
        raise SystemExit(f"[env] mmg3d failed rc={rc}, see {out/'mmg.log'}")

    M = read_mesh(out / "cut.mesh")
    Pm = M["Vertices"][:, :3]
    tet = M["Tetrahedra"]
    # mmg tags the two sides of the iso surface with different tet refs; keep
    # the side that carries the original inlet (ref S_INLET), not a fixed ref
    inlet_nodes = {int(x) - 1 for t in M["Triangles"] if int(t[3]) == v5.S_INLET for x in t[:3]}
    refs = collections.Counter(tet[:, 4].astype(int).tolist())
    side = {r: sum(1 for t in tet[tet[:, 4] == r][:, :4].astype(int) - 1
                   if inlet_nodes.intersection(t.tolist())) for r in refs}
    keep_ref = max(side, key=side.get)
    print(f"[env] tet refs {dict(refs)}, tets touching the inlet {side} -> keep ref {keep_ref}")
    keep = tet[tet[:, 4] == keep_ref][:, :4].astype(int) - 1
    print(f"[env] mmg: {len(tet)} tets, kept {len(keep)} outside the envelope")
    # boundary of the kept part, labels from mmg's surface refs
    tri_ref = {frozenset(t[:3].astype(int) - 1): int(t[3]) for t in M["Triangles"]}
    cnt = collections.Counter()
    for t in keep:
        for f in ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)):
            cnt[frozenset(t[list(f)])] += 1
    faces = [k for k, c in cnt.items() if c == 1]
    # outward orientation
    cen = {}
    for t in keep:
        c = Pm[t].mean(0)
        for f in ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)):
            k = frozenset(t[list(f)])
            if cnt[k] == 1:
                cen[k] = (t[list(f)], c)
    used_n = np.unique(keep)
    newi = -np.ones(len(Pm), int); newi[used_n] = np.arange(len(used_n))
    bt, bid, miss = [], [], 0
    for k in faces:
        f, c = cen[k]
        n = np.cross(Pm[f[1]] - Pm[f[0]], Pm[f[2]] - Pm[f[0]])
        if np.dot(n, Pm[f].mean(0) - c) < 0:
            f = f[::-1]
        r = tri_ref.get(k)
        if r is None:
            miss += 1; r = ISO_REF
        bid.append(v5.S_OGRID_IFACE if r in (ISO_REF, v5.S_OGRID_IFACE) else r)
        bt.append([newi[x] for x in f])
    bid = np.array(bid)
    print(f"[env] boundary {len(bt)} triangles ({miss} without mmg ref -> iso), "
          f"labels {dict(sorted(collections.Counter(bid.tolist()).items()))}")
    hist = tp.check_manifold([tuple(t) for t in bt])
    print(f"[env] edge multiplicity {dict(hist)}")
    Pk = Pm[used_n]
    tk = newi[keep]
    nd = {i: tuple(p) for i, p in enumerate(Pk)}
    feat, fverts = v2.feature_graph(nd, [tuple(t) for t in bt], bid)
    tp.write_algohex_vtk_analytic(Pk, tk, bt, bid, feat, fverts, out / "tet.vtk")
    print(f"[env] done: delta {delta:.4f}")


if __name__ == "__main__":
    main()
