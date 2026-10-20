"""AlgoHex input for a domain whose boundary carries named physical groups.

Generic counterpart of tet_prep_v5 (which hard-codes the dtOO ru_gridGmsh
entity ids): labels come from the MSH's own 2D physical groups, the volume is
re-tetrahedralised and made isotropic with mmg (wall-resolving hex layers give
needle-like boundary quads), feature curves = label boundaries + sharp edges
inside one label (dihedral above --sharp).

usage: tet_prep_named.py <msh> <out_dir> [--size-max H] [--hausd D] [--sharp DEG]
"""
import argparse, collections, subprocess, sys
from scipy.spatial import cKDTree
from pathlib import Path

import numpy as np
import meshio

HEX = Path(__file__).resolve().parent
sys.path.insert(0, str(HEX.parent.parent)); sys.path.insert(0, str(HEX))
import tet_prep as tp                                       # noqa: E402
import tet_prep_v2 as v2                                    # noqa: E402
import envelope_cut as ec                                   # noqa: E402

VOLFACES = {"tetra": ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)),
            "hexahedron": ((0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)),
            "pyramid": ((0, 1, 2, 3), (0, 1, 4), (1, 2, 4), (2, 3, 4), (3, 0, 4)),
            "wedge": ((0, 1, 2), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5))}


def split_quad(f):
    """Diagonal through the smallest vertex index: the same rule on both sides
    of every shared quad keeps the tet split conforming (orientation kept)."""
    i = int(np.argmin(f))
    a, b, c, d = (f[(i + k) % 4] for k in range(4))
    return [(a, b, c), (a, c, d)]


def to_tets(m):
    """Conforming tet split of the original volume mesh: tets as they are,
    pyramids into 2 (base diagonal by split_quad), hexes into 12 around a
    centroid (each face split by split_quad). No re-meshing of the interior,
    which gmsh cannot do from needle-like wall-layer faces (hung > 10 min)."""
    P = [m.points]
    n = len(m.points)
    tets = []
    for c in m.cells:
        if c.type == "tetra":
            tets.append(c.data)
        elif c.type == "pyramid":
            for cell in c.data:
                for t in split_quad(tuple(int(x) for x in cell[:4])):
                    tets.append(np.array([[*t, cell[4]]]))
        elif c.type == "hexahedron":
            cen = m.points[c.data].mean(1)
            ids = n + np.arange(len(c.data))
            P.append(cen)
            n += len(c.data)
            for cell, ci in zip(c.data, ids):
                for loc in VOLFACES["hexahedron"]:
                    for t in split_quad(tuple(int(cell[i]) for i in loc)):
                        tets.append(np.array([[*t, ci]]))
    T = np.vstack(tets)
    P = np.vstack(P)
    v = np.einsum("ij,ij->i", np.cross(P[T[:, 1]] - P[T[:, 0]], P[T[:, 2]] - P[T[:, 0]]), P[T[:, 3]] - P[T[:, 0]])
    T[v < 0] = T[v < 0][:, [0, 2, 1, 3]]
    print(f"[named] split into {len(T)} tets ({int((v == 0).sum())} flat)", flush=True)
    return P, T


def boundary(m):
    """Exterior faces of all volume cells, outward, triangulated, with the
    physical group of the matching tagged 2D element."""
    look = {}
    for c, ph in zip(m.cells, m.cell_data["gmsh:physical"]):
        if c.type in ("triangle", "quad"):
            for f, p in zip(c.data, ph):
                look[frozenset(f.tolist())] = int(p)
    cnt, own = collections.Counter(), {}
    for c in m.cells:
        if c.type not in VOLFACES:
            continue
        for cell in c.data:
            for loc in VOLFACES[c.type]:
                f = tuple(int(cell[i]) for i in loc)
                k = frozenset(f)
                cnt[k] += 1
                own[k] = (f, cell)
    P = m.points
    tris, ids, miss = [], [], 0
    for k, n in cnt.items():
        if n != 1:
            continue
        f, cell = own[k]
        pts = P[list(f)]
        nrm = np.cross(pts[1] - pts[0], pts[2] - pts[0])
        if np.dot(nrm, pts.mean(0) - P[cell].mean(0)) < 0:
            f = f[::-1]
        lab = look.get(k)
        if lab is None:
            miss += 1
            lab = 0
        for t in ([f] if len(f) == 3 else split_quad(f)):
            tris.append(t)
            ids.append(lab)
    print(f"[named] boundary {len(tris)} triangles, {miss} faces without a physical group")
    return tris, np.array(ids)


def sharp_edges(P, tris, ids, deg):
    e2t = collections.defaultdict(list)
    for i, t in enumerate(tris):
        for a in range(3):
            e2t[tuple(sorted((t[a], t[(a + 1) % 3])))].append(i)
    T = np.asarray(tris)
    n = np.cross(P[T[:, 1]] - P[T[:, 0]], P[T[:, 2]] - P[T[:, 0]])
    n /= np.linalg.norm(n, axis=1)[:, None]
    out = []
    for e, inc in e2t.items():
        if len(inc) == 2 and ids[inc[0]] == ids[inc[1]]:
            if np.degrees(np.arccos(np.clip(n[inc[0]] @ n[inc[1]], -1, 1))) > deg:
                out.append(e)
    return out


def anchored_chains(sharp, label_feat):
    """Keep a sharp-edge chain only if it is closed or both its ends sit on a
    label-boundary curve. A tightly rounded but smooth edge (canada duct lip)
    is detected only piecewise at a coarse surface size; its fragments end
    nowhere and would hand AlgoHex dozens of dangling curve ends."""
    adj = collections.defaultdict(list)
    for a, b in sharp:
        adj[a].append(b); adj[b].append(a)
    anchor = {v for e in label_feat for v in e}
    keep, seen = [], set()
    for s0 in adj:
        if s0 in seen:
            continue
        comp, st = [], [s0]
        seen.add(s0)
        while st:
            u = st.pop(); comp.append(u)
            for w in adj[u]:
                if w not in seen:
                    seen.add(w); st.append(w)
        ends = [v for v in comp if len(adj[v]) == 1]
        if all(v in anchor for v in ends):
            cs = set(comp)
            keep += [e for e in sharp if e[0] in cs]
    print(f"[named] sharp chains: kept {len(keep)} of {len(sharp)} edges (closed or anchored at both ends)")
    return keep


_E2T = {}


def edge_tris(btris, e):
    key = id(btris)
    if key not in _E2T:
        m = collections.defaultdict(list)
        for i, t in enumerate(btris):
            for x, y in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
                m[tuple(sorted((int(x), int(y))))].append(i)
        _E2T.clear(); _E2T[key] = m
    return _E2T[key][tuple(sorted((int(e[0]), int(e[1]))))]


def planar_labels(P, T, ids, tol=1e-9):
    """labels whose surface is a plane (symmetry planes, flat inlets): {label: (normal, offset)}"""
    out = {}
    span = np.linalg.norm(np.ptp(P, 0))
    for l in np.unique(ids):
        X = P[np.unique(T[ids == l])]
        c = X.mean(0)
        _, sv, vt = np.linalg.svd(X - c, full_matrices=False)
        n = vt[-1]
        if np.abs((X - c) @ n).max() < tol * span:
            out[int(l)] = (n, float(n @ c))
    print(f"[named] planar labels: {sorted(out)}")
    return out


def restore_planes(P, tets, T, ids, planes):
    """mmg slides points along the cut curves off a plane by up to ~1e-2; put the vertices of a
    planar label back onto it (a vertex on two planes onto their line), then check the tets"""
    P = P.copy()
    on = collections.defaultdict(set)
    for t, l in zip(T, ids):
        if int(l) in planes:
            for v in t:
                on[int(v)].add(int(l))
    moved = 0.0
    for v, ls in on.items():
        A = np.array([planes[l][0] for l in ls]); b = np.array([planes[l][1] for l in ls])
        x = P[v]
        # least-norm correction onto the intersection of the planes
        dx = A.T @ np.linalg.lstsq(A @ A.T, b - A @ x, rcond=None)[0]
        moved = max(moved, float(np.linalg.norm(dx)))
        P[v] = x + dx
    vol = np.einsum("ij,ij->i", np.cross(P[tets[:, 1]] - P[tets[:, 0]], P[tets[:, 2]] - P[tets[:, 0]]),
                    P[tets[:, 3]] - P[tets[:, 0]])
    print(f"[named] restored {len(on)} plane vertices, max move {moved:.2e}, non-positive tets after: {(vol <= 0).sum()}")
    if (vol <= 0).any():
        raise SystemExit("[named] plane restoration inverted tets")
    return P


def curve_corners(P, feat, deg):
    """valence-2 feature vertices where the curve turns by more than deg: real corners
    (hub end cap, blunt trailing edge) that AlgoHex must see as feature vertices"""
    adj = collections.defaultdict(list)
    for x, y in feat:
        adj[x].append(y); adj[y].append(x)
    out = []
    for v, nb in adj.items():
        if len(nb) == 2:
            d1, d2 = P[v] - P[nb[0]], P[nb[1]] - P[v]
            c = np.dot(d1, d2) / (np.linalg.norm(d1) * np.linalg.norm(d2) + 1e-300)
            if np.degrees(np.arccos(np.clip(c, -1, 1))) > deg:
                out.append(v)
    print(f"[named] curve corners (turn > {deg} deg) added as feature vertices: {len(out)}")
    return out


def remesh_surface(P, T, ids, out, h, hausd, hgrad):
    MMGS = ec.MMG.replace("mmg3d_O3", "mmgs_O3")
    u = np.unique(T)
    mp = -np.ones(len(P), int); mp[u] = np.arange(len(u))
    with open(out / "surf_in.mesh", "w") as f:
        f.write(f"MeshVersionFormatted 2\nDimension 3\nVertices\n{len(u)}\n")
        for p in P[u]:
            f.write(f"{p[0]:.17g} {p[1]:.17g} {p[2]:.17g} 0\n")
        f.write(f"Triangles\n{len(T)}\n")
        for t, r in zip(mp[T], ids):
            f.write(f"{t[0]+1} {t[1]+1} {t[2]+1} {r}\n")
        f.write("End\n")
    cmd = [MMGS, "-in", str(out / "surf_in.mesh"), "-hmax", str(h), "-hausd", str(hausd),
           "-hgrad", str(hgrad), "-out", str(out / "surf_out.mesh")]
    with open(out / "mmgs.log", "w") as lf:
        rc = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT).returncode
    if rc:
        raise SystemExit(f"[named] mmgs failed rc={rc}")
    M = ec.read_mesh(out / "surf_out.mesh")
    S = M["Vertices"][:, :3]
    ST = M["Triangles"][:, :3].astype(int) - 1
    SR = M["Triangles"][:, 3].astype(int)
    hist = tp.check_manifold([tuple(t) for t in ST])
    print(f"[named] surface remeshed: {len(T)} -> {len(ST)} triangles, edge multiplicity {dict(hist)}", flush=True)
    return S, ST, SR


def envelope(P, tets, btris, ids, a, out):
    """Cut the named thin walls out by offset envelopes {d_i = delta}, ONE BODY AT A TIME
    (edge-connected wall components), each with its own label. Where two envelopes meet, mmg
    cuts the second iso surface exactly into the first, so the crease is a complete label
    boundary instead of a piecewise-detected 'sharp' chain."""
    names = {int(l.split()[0]): l.split()[1] for l in open(out / "labels.txt")}
    inv = {v: k for k, v in names.items()}
    cut = [inv[n] for n in a.envelope.split(",")]
    B = np.asarray(btris)
    sel = np.isin(ids, cut)
    e2t = collections.defaultdict(list)
    for i in np.flatnonzero(sel):
        t = B[i]
        for x, y in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
            e2t[tuple(sorted((int(x), int(y))))].append(i)
    adj = collections.defaultdict(set)
    for inc in e2t.values():
        for x in inc:
            adj[x].update(inc)
    comp, c = {}, 0
    for s0 in np.flatnonzero(sel):
        if s0 in comp:
            continue
        st = [s0]; comp[s0] = c
        while st:
            u = st.pop()
            for w in adj[u]:
                if w not in comp:
                    comp[w] = c; st.append(w)
        c += 1
    bodies = [B[[j for j in comp if comp[j] == i]] for i in range(c)]
    gaps = []
    for i in range(c):
        for j in range(i + 1, c):
            gaps.append(ec.surface_distance(P, bodies[i], P[np.unique(bodies[j])]).min())
    wall_v = set(np.unique(B[sel]).tolist())
    for l in sorted(set(ids.tolist()) - set(cut)):
        lv = np.unique(B[ids == l])
        if wall_v.intersection(lv.tolist()):
            continue
        gaps.append(ec.surface_distance(P, B[sel], P[lv]).min())
    delta = a.delta or a.delta_gap * min(gaps)
    h_iso = delta * np.radians(15)
    sag = delta * (1 - np.cos(np.radians(15)))
    print(f"[named] envelope around {a.envelope}: {c} bodies, gaps {np.round(sorted(gaps)[:4], 3)}, "
          f"delta {delta:.4f}", flush=True)
    WV = P.copy()
    curved_all = np.unique(np.array(sharp_edges(P, [tuple(t) for t in btris], np.where(sel, 1, 0), 20.0)
                                    or [(0, 0)]).ravel())
    new0 = max(names) + 1
    ec.HGRAD = a.hgrad
    for i, F in enumerate(bodies):
        ls_fun = lambda X, F=F: ec.surface_distance(WV, F, X) - delta
        bv = set(np.unique(F).tolist())
        curved = [v for v in curved_all if v in bv] or list(bv)[:1]
        tc = cKDTree(WV[curved])
        hmap = lambda X: np.clip(np.minimum(delta, h_iso + 0.5 * tc.query(X)[0]) + 0.5 * np.abs(ls_fun(X)),
                                 h_iso, a.size_max)
        # pre-refine around this body's iso surface (mmg -ls never re-evaluates the level set)
        ec.write_mesh(out / f"envpre{i}_in.mesh", P, tets, btris, ids)
        ec.write_sol(out / f"envpre{i}_in.sol", hmap(P))
        fine = [new0 + k for k in range(i)]
        if fine:
            with open(out / f"envpre{i}_in.mmg3d", "w") as f:
                f.write(f"Parameters\n{len(fine)}\n")
                for r in fine:
                    f.write(f"{r} Triangles {h_iso:.6g} {a.size_max:.6g} {sag:.6g}\n")
        cmd = [ec.MMG, "-in", str(out / f"envpre{i}_in.mesh"), "-sol", str(out / f"envpre{i}_in.sol"),
               "-hausd", str(a.hausd), "-hgrad", str(a.hgrad), "-out", str(out / f"envpre{i}_out.mesh")]
        with open(out / f"envpre{i}_mmg.log", "w") as lf:
            if subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT).returncode:
                raise SystemExit(f"[named] mmg pre-refinement {i} failed")
        M = ec.read_mesh(out / f"envpre{i}_out.mesh")
        P = M["Vertices"][:, :3]; tets = M["Tetrahedra"][:, :4].astype(int) - 1
        btris = (M["Triangles"][:, :3].astype(int) - 1).tolist(); ids = M["Triangles"][:, 3].astype(int)
        P, tets, btris, ids = ec.ls_cut(P, tets, btris, ids, ls_fun, new0 + i, out, f"env{i}", a.size_max,
                                        a.hausd, iso_hausd=sag, fine_refs=fine, iso_hmin=h_iso)
        print(f"[named] body {i}: {len(tets)} tets after its cut", flush=True)
    left = set(ids.tolist()) & set(cut)
    if left:
        raise SystemExit(f"[named] wall labels {left} survived the envelope cuts")
    env = [new0 + i for i in range(c)]
    P, tets, btris, ids = ec.coarsen(P, tets, btris, ids, out, a.size_max, a.hausd, sag, refs=env, cut_hmin=h_iso)
    with open(out / "labels.txt", "a") as f:
        for i in range(c):
            f.write(f"{new0 + i} ENVELOPE_{i} (body {i} of {a.envelope})\n")
    return np.asarray(P), tets, btris, ids, env


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("msh")
    ap.add_argument("out_dir")
    ap.add_argument("--size-max", type=float, default=0.3)
    ap.add_argument("--hausd", type=float, default=0.01)
    ap.add_argument("--hgrad", type=float, default=1.5)
    ap.add_argument("--envelope", default=None,
                    help="comma-separated physical names of thin walls to cut out by an offset envelope")
    ap.add_argument("--corner", type=float, default=45.0,
                    help="feature-curve turn (deg) above which a curve vertex becomes a feature vertex")
    ap.add_argument("--delta-gap", type=float, default=0.25,
                    help="envelope distance = this fraction of the narrowest gap around the thin walls")
    ap.add_argument("--delta", type=float, default=None)
    ap.add_argument("--sharp", type=float, default=40.0,
                    help="dihedral (deg, 0 = flat) above which an edge inside one label is a feature")
    a = ap.parse_args()
    out = Path(a.out_dir); out.mkdir(parents=True, exist_ok=True)

    m = meshio.read(a.msh)
    names = {v[0]: k for k, v in m.field_data.items() if v[1] == 2}
    tris, ids = boundary(m)
    hist = tp.check_manifold(tris)
    if set(hist) - {2}:
        raise SystemExit(f"[named] boundary not a closed 2-manifold: {dict(hist)}")
    # compact label ids 1..k (AlgoHex colours), keep the mapping
    labs = sorted(set(ids.tolist()))
    lid = {l: i + 1 for i, l in enumerate(labs)}
    with open(out / "labels.txt", "w") as f:
        for l in labs:
            line = f"{lid[l]} {names.get(l, 'unnamed')} ({int((ids == l).sum())} tris)"
            f.write(line + "\n")
            print("[named]  ", line)
    ids = np.array([lid[l] for l in ids])

    # isotropic SURFACE remesh with mmgs (refs and ridges kept), then a fresh
    # interior tet mesh from it. Re-meshing the volume instead (dtOO hexes
    # split into 1.9M tets, mmg3d) did not finish in 28 min; gmsh cannot mesh
    # the interior from the original needle-like wall-layer faces either.
    S, ST, SR = remesh_surface(m.points, np.asarray(tris), ids, out, a.size_max, a.hausd, a.hgrad)
    nodes = {i: tuple(p) for i, p in enumerate(S)}
    P, tets, idx, used = v2.mesh_interior(nodes, [tuple(t) for t in ST], a.size_max)
    remap = {ot: idx[nt] for nt, ot in zip(range(1, len(used) + 1), used) if nt in idx}
    btris = [[remap[t] for t in tri] for tri in ST]
    ids = SR
    P = np.asarray(P)
    env = []
    planes = planar_labels(P, np.asarray(btris), ids)
    if a.envelope:
        P, tets, btris, ids, env = envelope(P, tets, btris, ids, a, out)
        P = restore_planes(P, np.asarray(tets), np.asarray(btris), ids, planes)
    nd = {i: tuple(p) for i, p in enumerate(P)}
    feat, _ = v2.feature_graph(nd, [tuple(t) for t in btris], ids)
    # sharp edges only on surfaces that are not envelopes (an envelope is smooth by construction;
    # its creases with other envelopes are label boundaries already)
    sh = [e for e in sharp_edges(P, btris, ids, a.sharp) if not set(env) & {int(x) for x in ids[edge_tris(btris, e)]}]
    sh = anchored_chains(sh, feat)
    feat = sorted(set(feat) | set(sh))
    fverts = sorted(set(ec.v5.feature_vertices(feat)) | set(curve_corners(P, feat, a.corner)))
    print(f"[named] feature edges {len(feat)} (sharp inside a label: {len(sh)}), feature vertices {len(fverts)}")
    tp.write_algohex_vtk_analytic(P, tets, btris, ids, feat, fverts, out / "tet.vtk")


if __name__ == "__main__":
    main()
