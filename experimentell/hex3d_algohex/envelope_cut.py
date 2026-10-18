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

--domain fluid (no dependence on the dtOO O-grid thickness, delta free):
  1. domain = tets + pyramids + O-grid hexes (only the BL prisms out)
  2. three level-set cuts: hub BL (d_hub < t_h, label 5), shroud BL (label 6),
     envelope (d_body < delta, label 7); t_h/t_s = top of the prism layer,
     delta = 1/2 of the body's distance to inlet/outlet/periodic
  3. each cut pre-refines around its iso surface (mmg -ls never re-evaluates
     the level set), a final pass coarsens the volume to size-max again

usage: envelope_cut.py <msh> <out_dir> [--domain reduced|fluid] [--size-max 0.06]
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
    caps = {g: [fo[k][0] for k in fs] for g, fs in ent.items() if g not in body}
    return V, F, body, caps


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


# --------------------------------------------------------------- walls
def wall_surfaces(nodes, els, caps):
    """Hub and shroud wall as triangle meshes: the tagged hub/shroud elements
    plus the caps the removed O-grid puts on them (`caps` from body_wall: the
    O-grid exterior faces that touch neither the remaining volume nor the
    body), each assigned to the wall whose nodes it shares."""
    out = {}
    for g in (v5.GEOM_HUB, v5.GEOM_SHROUD):
        faces = [e[3] for e in els if e[0] in (2, 3) and e[2] == g]
        wn = {v for f in faces for v in f}
        mine = [c for c, fs in caps.items() if any(v in wn for f in fs for v in f)]
        faces += [f for c in mine for f in caps[c]]
        tri = []
        for f in faces:
            tri.append(f[:3])
            if len(f) == 4:
                tri.append((f[0], f[2], f[3]))
        tri = np.array(tri)
        u = np.unique(tri)
        m = {v: i for i, v in enumerate(u)}
        out[g] = (np.array([nodes[v] for v in u], float), np.vectorize(m.get)(tri))
        print(f"[env] wall {g}: {len(tri)} triangles, caps {sorted(mine)}")
    return out[v5.GEOM_HUB], out[v5.GEOM_SHROUD]


def fluid_boundary(nodes, els, body_ents):
    """Boundary of tets + pyramids + O-grid hexes (only the BL prisms out).
    Labels: inlet/outlet/periodic as in tet_prep_v5; everything else is a
    temporary label that the level-set cuts must remove completely."""
    keep = [el for el in els if el[0] in (4, 7) or (el[0] == 5 and el[2] in v5.OGRID_GEOM)]
    tris = tp.orient_and_triangulate(tp.volume_exterior_faces(keep), nodes)
    hist = tp.check_manifold(tris)
    if set(hist) - {2}:
        raise RuntimeError(f"fluid boundary is not a closed 2-manifold: {dict(hist)}")
    look = v5.tagged_2d_lookup(els)
    to_sid = {v5.GEOM_INLET: v5.S_INLET, v5.GEOM_OUTLET: v5.S_OUTLET,
              v5.GEOM_PER_A: v5.S_PER_A, v5.GEOM_PER_B: v5.S_PER_B}
    ids = np.array([to_sid.get(look.get(frozenset(t)), TMP) for t in tris])
    print(f"[env] fluid boundary {len(tris)} triangles, labels "
          f"{dict(sorted(collections.Counter(ids.tolist()).items()))}")
    return tris, ids


TMP = 21                                                   # cut away completely
HGRAD = None                                               # None = mmg's own default


def hgrad_opt():
    # pass -hgrad only when set: an explicit "-hgrad 1.3" does not reproduce
    # mmg's default 1.3 (16388 vs 16424 vertices on the canada runner)
    return ["-hgrad", str(HGRAD)] if HGRAD else []


# ---------------------------------------------------------------- ls cut
def ls_cut(P, tets, btris, ids, ls_fun, new_label, out, tag, size_max, hausd,
           iso_hausd=None, h_iso=None):
    """mmg3d -ls on (P, tets) with surface refs `ids` and level set ls_fun(P);
    keep the ls > 0 side, label its new iso triangles `new_label`.
    `iso_hausd` refines only the new iso surface (mmg local parameters on its
    ref): a small envelope radius needs it at the rounded TE, a global hausd
    would refine every surface (canada guide vane: 318k instead of ~100k tets)."""
    if h_iso:
        # mmg -ls interpolates the level set linearly on the INPUT mesh and
        # never re-evaluates it, so the input must already resolve the iso
        # surface: a coarse mesh (0.12) around a radius-0.033 TE rounding left
        # 60 deg kinks in the rings. Pre-pass: size map h = h_iso near the iso,
        # growing with |ls| up to size_max.
        write_mesh(out / f"{tag}_pre.mesh", P, tets, btris, ids)
        h = np.clip(h_iso + 0.5 * np.abs(ls_fun(P)), h_iso, size_max)
        write_sol(out / f"{tag}_pre.sol", h)
        cmd = [MMG, "-in", str(out / f"{tag}_pre.mesh"), "-sol", str(out / f"{tag}_pre.sol"),
               "-hausd", str(hausd), *hgrad_opt(), "-out", str(out / f"{tag}_pre_out.mesh")]
        with open(out / f"{tag}_pre_mmg.log", "w") as lf:
            rc = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT).returncode
        if rc:
            raise SystemExit(f"[env] mmg3d pre-refinement failed rc={rc}")
        M0 = read_mesh(out / f"{tag}_pre_out.mesh")
        P = M0["Vertices"][:, :3]
        tets = M0["Tetrahedra"][:, :4].astype(int) - 1
        btris = (M0["Triangles"][:, :3].astype(int) - 1).tolist()
        ids = M0["Triangles"][:, 3].astype(int)
        print(f"[env] {tag}: pre-refined to h {h_iso:.4f} near the iso surface: "
              f"{len(P)} nodes, {len(tets)} tets", flush=True)
    write_mesh(out / f"{tag}_in.mesh", P, tets, btris, ids)
    write_sol(out / f"{tag}_in.sol", ls_fun(P))
    if iso_hausd:
        # the rings lie on the iso surface AND on the earlier cut surfaces
        refs = [ISO_REF] + [r for r in (v5.S_BL_IFACE_HUB, v5.S_BL_IFACE_SHROUD) if r in set(ids)]
        with open(out / f"{tag}_in.mmg3d", "w") as f:
            f.write(f"Parameters\n{len(refs)}\n")
            for r in refs:
                f.write(f"{r} Triangles {iso_hausd / 10:.6g} {size_max:.6g} {iso_hausd:.6g}\n")
    cmd = [MMG, "-in", str(out / f"{tag}_in.mesh"), "-sol", str(out / f"{tag}_in.sol"),
           "-ls", "0", "-hmax", str(size_max), "-hausd", str(hausd), *hgrad_opt(),
           "-out", str(out / f"{tag}_cut.mesh")]
    print("[env]", " ".join(cmd), flush=True)
    with open(out / f"{tag}_mmg.log", "w") as lf:
        rc = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT).returncode
    if rc:
        raise SystemExit(f"[env] mmg3d failed rc={rc}, see {out / (tag + '_mmg.log')}")
    M = read_mesh(out / f"{tag}_cut.mesh")
    Pm = M["Vertices"][:, :3]
    tet = M["Tetrahedra"]
    # mmg tags the two sides with different tet refs; keep the one whose
    # centroids lie on the positive side, not a fixed ref number. The .sol
    # mmg writes back is its size metric, not the level set: evaluate anew.
    refs = sorted(set(tet[:, 4].astype(int).tolist()))
    cen_ls = ls_fun(Pm[tet[:, :4].astype(int) - 1].mean(1))
    mean_ls = {r: float(np.median(cen_ls[tet[:, 4] == r])) for r in refs}
    keep_ref = max(mean_ls, key=mean_ls.get)
    keep = tet[tet[:, 4] == keep_ref][:, :4].astype(int) - 1
    print(f"[env] {tag}: mean ls per tet ref {mean_ls} -> keep ref {keep_ref}, "
          f"{len(keep)} of {len(tet)} tets")
    tri_ref = {frozenset(t[:3].astype(int) - 1): int(t[3]) for t in M["Triangles"]}
    FACES = ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))
    cnt, cen = collections.Counter(), {}
    for t in keep:
        for f in FACES:
            cnt[frozenset(t[list(f)])] += 1
    for t in keep:
        c = Pm[t].mean(0)
        for f in FACES:
            k = frozenset(t[list(f)])
            if cnt[k] == 1:
                cen[k] = (t[list(f)], c)
    used_n = np.unique(keep)
    newi = -np.ones(len(Pm), int)
    newi[used_n] = np.arange(len(used_n))
    bt, bid, miss = [], [], 0
    for k, (f, c) in cen.items():
        n = np.cross(Pm[f[1]] - Pm[f[0]], Pm[f[2]] - Pm[f[0]])
        if np.dot(n, Pm[f].mean(0) - c) < 0:
            f = f[::-1]
        r = tri_ref.get(k)
        if r is None:
            miss += 1
            r = ISO_REF
        bid.append(new_label if r == ISO_REF else r)
        bt.append([newi[x] for x in f])
    bid = np.array(bid)
    hist = tp.check_manifold([tuple(t) for t in bt])
    print(f"[env] {tag}: boundary {len(bt)} triangles ({miss} without ref), labels "
          f"{dict(sorted(collections.Counter(bid.tolist()).items()))}, edge multiplicity {dict(hist)}")
    return Pm[used_n], newi[keep], bt, bid


def coarsen(P, tets, btris, ids, out, size_max, hausd, cut_hausd):
    """Final mmg3d pass without level set: the cuts leave their fine near-iso
    sizing in the whole volume (canada guide vane ~170-230k tets); remesh to
    size_max inside while the three cut surfaces keep cut_hausd, so their
    rings stay resolved."""
    write_mesh(out / "final_in.mesh", P, tets, btris, ids)
    with open(out / "final_in.mmg3d", "w") as f:
        refs = [r for r in (v5.S_BL_IFACE_HUB, v5.S_BL_IFACE_SHROUD, v5.S_OGRID_IFACE) if r in set(ids)]
        f.write(f"Parameters\n{len(refs)}\n")
        for r in refs:
            f.write(f"{r} Triangles {cut_hausd / 10:.6g} {size_max:.6g} {cut_hausd:.6g}\n")
    cmd = [MMG, "-in", str(out / "final_in.mesh"), "-hmax", str(size_max), "-hausd", str(hausd),
           *hgrad_opt(), "-out", str(out / "final_out.mesh")]
    with open(out / "final_mmg.log", "w") as lf:
        rc = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT).returncode
    if rc:
        raise SystemExit(f"[env] mmg3d final coarsening failed rc={rc}")
    M = read_mesh(out / "final_out.mesh")
    tet = M["Tetrahedra"][:, :4].astype(int) - 1
    used = np.unique(tet)
    newi = -np.ones(len(M["Vertices"]), int)
    newi[used] = np.arange(len(used))
    # boundary = triangles on exactly one tet; mmg keeps their refs
    cnt = collections.Counter(frozenset(t[list(f)]) for t in tet
                              for f in ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)))
    bt, bid = [], []
    for t in M["Triangles"]:
        k = frozenset(t[:3].astype(int) - 1)
        if cnt.get(k) == 1:
            bt.append([newi[x] for x in t[:3].astype(int) - 1])
            bid.append(int(t[3]))
    # orient outward
    Pn = M["Vertices"][used, :3]
    tn = newi[tet]
    cen = {}
    for t in tn:
        c = Pn[t].mean(0)
        for f in ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)):
            cen[frozenset(t[list(f)])] = c
    for i, f in enumerate(bt):
        n = np.cross(Pn[f[1]] - Pn[f[0]], Pn[f[2]] - Pn[f[0]])
        if np.dot(n, Pn[f].mean(0) - cen[frozenset(f)]) < 0:
            bt[i] = f[::-1]
    hist = tp.check_manifold([tuple(t) for t in bt])
    print(f"[env] final: {len(Pn)} nodes, {len(tn)} tets, {len(bt)} boundary triangles, "
          f"edge multiplicity {dict(hist)}")
    return Pn, tn, bt, np.array(bid)


def distances(BV, BF, nodes, T, ids, names):
    dist = {}
    for sid in np.unique(ids):
        pts = np.array([nodes[v] for v in np.unique(T[ids == sid])])
        d = surface_distance(BV, BF, pts)
        dist[int(sid)] = (d.min(), d.max())
        print(f"[env]   {names.get(int(sid), sid)!s:20s} distance to body [{d.min():.4f}, {d.max():.4f}]")
    return dist


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("msh")
    ap.add_argument("out_dir")
    ap.add_argument("--domain", choices=("reduced", "fluid"), default="reduced",
                    help="reduced: envelope around the cut-out dtOO O-grid (delta must exceed "
                         "its thickness); fluid: O-grid kept in the domain, hub/shroud BL and "
                         "envelope cut by three level sets, delta free")
    ap.add_argument("--delta-factor", type=float, default=1.25,
                    help="reduced: delta = factor * max distance of the removed O-grid's outer face")
    ap.add_argument("--delta-clear", type=float, default=0.5,
                    help="fluid: delta = this fraction of the body's distance to inlet/outlet/periodic")
    ap.add_argument("--delta", type=float, default=None, help="absolute delta, overrides both")
    ap.add_argument("--bl-factor", type=float, default=1.05,
                    help="fluid: hub/shroud cut height = factor * highest prism-layer top")
    ap.add_argument("--size-max", type=float, default=None,
                    help="default 0.06 reduced, 0.08 fluid (canada guide vane at 0.12: "
                         "19-38k tets, 86-116 singular edges or a broken parametrisation; "
                         "at 0.08: 44-50k tets, the 52 singular edges / 12 raw blocks of tistos)")
    ap.add_argument("--hausd", type=float, default=0.002)
    ap.add_argument("--hgrad", type=float, default=None,
                    help="mmg size gradation (default: mmg's own for reduced, 2.0 for fluid)")
    a = ap.parse_args()
    global HGRAD
    HGRAD = a.hgrad or (2.0 if a.domain == "fluid" else None)
    if a.size_max is None:
        a.size_max = 0.08 if a.domain == "fluid" else 0.06
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    nodes, els = parse_msh(Path(a.msh))
    BV, BF, body_ents, caps = body_wall(nodes, els)
    names = dict(v5.NAMES)
    names[TMP] = "tmp (to cut)"
    if a.domain == "reduced":
        tris = v5.reduced_boundary(nodes, els)
        ids = v5.classify(nodes, tris, v5.tagged_2d_lookup(els), v5.removed_face_kind(els))
        dist = distances(BV, BF, nodes, np.array(tris), ids, names)
        delta = a.delta or a.delta_factor * dist[v5.S_OGRID_IFACE][1]
    else:
        tris, ids = fluid_boundary(nodes, els, body_ents)
        dist = distances(BV, BF, nodes, np.array(tris), ids, names)
        clear0 = min(dist[s][0] for s in (v5.S_INLET, v5.S_OUTLET, v5.S_PER_A, v5.S_PER_B))
        delta = a.delta or a.delta_clear * clear0
    clear = min(dist[s][0] for s in (v5.S_INLET, v5.S_OUTLET, v5.S_PER_A, v5.S_PER_B) if s in dist)
    print(f"[env] delta = {delta:.4f}; closest inlet/outlet/periodic = {clear:.4f}")
    if delta >= clear:
        raise SystemExit("[env] envelope would hit inlet/outlet/periodic; lower delta")

    P, tets, idx, used = v2.mesh_interior(nodes, tris, a.size_max)
    remap = {ot: idx[nt] for nt, ot in zip(range(1, len(used) + 1), used) if nt in idx}
    btris = [[remap[t] for t in tri] for tri in tris]

    if a.domain == "fluid":
        (HV, HF), (SV, SF) = wall_surfaces(nodes, els, caps)
        # cut height = top of the BL prism layer, measured on the prisms themselves
        pn = np.array([nodes[v] for v in {v for e in els if e[0] == 6 for v in e[3]}], float)
        dh, ds = surface_distance(HV, HF, pn), surface_distance(SV, SF, pn)
        t_h = a.bl_factor * dh[dh < ds].max()
        t_s = a.bl_factor * ds[ds < dh].max()
        print(f"[env] BL cut heights: hub {t_h:.4f}, shroud {t_s:.4f}")
        for tag, (V, F, t, lab) in (("hub", (HV, HF, t_h, v5.S_BL_IFACE_HUB)),
                                    ("shroud", (SV, SF, t_s, v5.S_BL_IFACE_SHROUD))):
            P, tets, btris, ids = ls_cut(P, tets, btris, ids,
                                         lambda X, V=V, F=F, t=t: surface_distance(V, F, X) - t,
                                         lab, out, tag, a.size_max, a.hausd)
    # turn angle along the rounded TE <= ~15 deg: chord sag of a radius-delta arc
    iso_hausd = min(a.hausd, delta * (1 - np.cos(np.radians(15))))
    P, tets, btris, ids = ls_cut(P, tets, btris, ids,
                                 lambda X: surface_distance(BV, BF, X) - delta,
                                 v5.S_OGRID_IFACE, out, "body", a.size_max, a.hausd,
                                 iso_hausd if a.domain == "fluid" else None,
                                 # ~15 deg per segment around the radius-delta rounding
                                 delta * np.radians(15) if a.domain == "fluid" else None)
    if a.domain == "fluid":
        P, tets, btris, ids = coarsen(P, tets, btris, ids, out, a.size_max, a.hausd, iso_hausd)
    left = set(ids.tolist()) - set(range(1, 8))
    if left:
        raise SystemExit(f"[env] labels {sorted(left)} survived the cuts")
    nd = {i: tuple(p) for i, p in enumerate(P)}
    feat, fverts = v2.feature_graph(nd, [tuple(t) for t in btris], ids)
    tp.write_algohex_vtk_analytic(P, tets, btris, ids, feat, fverts, out / "tet.vtk")
    print(f"[env] done: domain {a.domain}, delta {delta:.4f}")


if __name__ == "__main__":
    main()
