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

HEX = Path("/opt/stack/domain_partition_3D/experimentell/hex3d_algohex")
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


BOUNDARY_ENT = [None]
EXTRA_FEAT = [[]]          # explicit feature edges (crease ridges) from envelope()     # dtOO 2D entity of every boundary triangle (set by boundary)


def boundary(m):
    """Exterior faces of all volume cells, outward, triangulated, with the
    physical group of the matching tagged 2D element."""
    look, looke = {}, {}
    for c, ph, ge in zip(m.cells, m.cell_data["gmsh:physical"], m.cell_data["gmsh:geometrical"]):
        if c.type in ("triangle", "quad"):
            for f, p, g in zip(c.data, ph, ge):
                look[frozenset(f.tolist())] = int(p)
                looke[frozenset(f.tolist())] = int(g)
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
    tris, ids, miss, ents = [], [], 0, []
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
            ents.append(looke.get(k, -1))
    print(f"[named] boundary {len(tris)} triangles, {miss} faces without a physical group")
    BOUNDARY_ENT[0] = np.array(ents)
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


def extend_chains(P, btris, ids, sharp, label_feat, low=5.0, cosmin=0.8, steps=40, gap=0.15):
    """hysteresis: continue every sharp chain end that is not on a label boundary along weaker
    edges (dihedral > low, same label, continuing the direction) until it reaches one. The ring-3
    trailing edge (dihedral 90) faded below the threshold just before the symmetry plane and was
    dropped as 'unanchored'."""
    T = np.asarray(btris)
    n = np.cross(P[T[:, 1]] - P[T[:, 0]], P[T[:, 2]] - P[T[:, 0]])
    n /= np.linalg.norm(n, axis=1)[:, None]
    e2t = collections.defaultdict(list)
    vadj = collections.defaultdict(set)
    for i, t in enumerate(T):
        for x, y in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
            e2t[tuple(sorted((int(x), int(y))))].append(i)
            vadj[int(x)].add(int(y)); vadj[int(y)].add(int(x))
    def dih(e):
        inc = e2t[e]
        if len(inc) != 2 or ids[inc[0]] != ids[inc[1]]:
            return -1
        return np.degrees(np.arccos(np.clip(n[inc[0]] @ n[inc[1]], -1, 1)))
    anchor = {v for e in label_feat for v in e}
    S = {tuple(sorted(map(int, e))) for e in sharp}
    deg = collections.Counter(v for e in S for v in e)
    added, joined = 0, 0
    for v0 in [v for v, c in deg.items() if c == 1 and v not in anchor]:
        prev = next(u for e in S if v0 in e for u in e if u != v0)
        path, v, ok = [], v0, False
        for _ in range(steps):
            d = P[v] - P[prev]; d /= np.linalg.norm(d)
            best = None
            for w in vadj[v]:
                if w == prev:
                    continue
                e = tuple(sorted((v, w)))
                g = P[w] - P[v]; c = d @ g / np.linalg.norm(g)
                a = dih(e)
                if c > cosmin and a > low and (best is None or a > best[0]):
                    best = (a, w, e)
            if best is None:
                break
            path.append(best[2]); prev, v = v, best[1]
            if v in anchor or deg.get(v, 0) >= 1:
                ok = True
                break
        if ok:
            S.update(path); added += len(path)
            continue
        # gap closing: a fading edge runs INTO another curve (canada ring-3 bevel edge merging into
        # the trailing knife edge): join the end to the nearest curve vertex within `gap` along the
        # shortest path on the surface (same label only)
        # exclude the chain's own component (else the search 'joins' it back onto itself)
        sadj = collections.defaultdict(list)
        for e in S:
            sadj[e[0]].append(e[1]); sadj[e[1]].append(e[0])
        own, st = {v0}, [v0]
        while st:
            u = st.pop()
            for w in sadj[u]:
                if w not in own:
                    own.add(w); st.append(w)
        targets = (anchor | {u for e in S for u in e}) - own
        import heapq
        lab0 = {int(ids[i]) for e in S if v0 in e for i in e2t[e]}
        dist, prevv, hq, hit = {v0: 0.0}, {}, [(0.0, v0)], None
        while hq:
            dcur, u = heapq.heappop(hq)
            if dcur > gap:
                break
            if u in targets:
                hit = u; break
            for w in vadj[u]:
                e = tuple(sorted((u, w)))
                if not ({int(ids[i]) for i in e2t[e]} & lab0):
                    continue
                nd = dcur + np.linalg.norm(P[w] - P[u])
                if nd < dist.get(w, np.inf):
                    dist[w] = nd; prevv[w] = u; heapq.heappush(hq, (nd, w))
        if hit is not None:
            u = hit
            while u != v0:
                S.add(tuple(sorted((u, prevv[u])))); added += 1; u = prevv[u]
            joined += 1
    print(f"[named] hysteresis: extended chain ends by {added} edges ({joined} joined to a nearby curve)")
    return sorted(S)


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


def restore_planes(P, tets, T, ids, planes, fatal=True):
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
        if fatal:
            raise SystemExit("[named] plane restoration inverted tets")
        return None
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


def remesh_surface(P, T, ids, out, h, hausd, hgrad, hmap=None, ridge_angle=None):
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
    if hmap is not None:
        ec.write_sol(out / "surf_in.sol", hmap(P[u]))
        cmd[3:3] = ["-sol", str(out / "surf_in.sol")]
    if ridge_angle is not None:
        cmd[3:3] = ["-ar", str(ridge_angle)]
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


def final_coarsen(P, tets, btris, ids, out, a, env, sag, h_iso, planes):
    """like envelope_cut.coarsen, but planar labels get a near-zero Hausdorff tolerance: with the
    global one mmg slid plane vertices up to 0.04 off the symmetry planes (197 tets inverted
    when put back)"""
    ec.write_mesh(out / "final_in.mesh", P, tets, btris, ids)
    present = set(np.asarray(ids).tolist())
    rows = [f"{r} Triangles {h_iso:.6g} {a.size_max:.6g} {sag:.6g}" for r in env if r in present]
    rows += [f"{r} Triangles {a.size_max / 10:.6g} {a.size_max:.6g} 1e-04" for r in (planes or {}) if r in present]
    with open(out / "final_in.mmg3d", "w") as f:
        f.write(f"Parameters\n{len(rows)}\n" + "\n".join(rows) + "\n")
    cmd = [ec.MMG, "-in", str(out / "final_in.mesh"), "-hmax", str(a.size_max), "-hausd", str(a.hausd),
           "-hgrad", str(a.hgrad), "-out", str(out / "final_out.mesh")]
    with open(out / "final_mmg.log", "w") as lf:
        if subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT).returncode:
            raise SystemExit("[named] final mmg pass failed")
    M = ec.read_mesh(out / "final_out.mesh")
    tet = M["Tetrahedra"][:, :4].astype(int) - 1
    used = np.unique(tet)
    newi = -np.ones(len(M["Vertices"]), int); newi[used] = np.arange(len(used))
    cnt = collections.Counter(frozenset(t[list(f)]) for t in tet for f in ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)))
    bt, bid = [], []
    for t in M["Triangles"]:
        k = frozenset(t[:3].astype(int) - 1)
        if cnt.get(k) == 1:
            bt.append([newi[x] for x in t[:3].astype(int) - 1]); bid.append(int(t[3]))
    Pn = M["Vertices"][used, :3]; tn = newi[tet]
    cen = {}
    for t in tn:
        c = Pn[t].mean(0)
        for f in ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)):
            cen[frozenset(t[list(f)])] = c
    for i, f in enumerate(bt):
        nn = np.cross(Pn[f[1]] - Pn[f[0]], Pn[f[2]] - Pn[f[0]])
        if np.dot(nn, Pn[f].mean(0) - cen[frozenset(f)]) < 0:
            bt[i] = f[::-1]
    print(f"[named] final: {len(Pn)} nodes, {len(tn)} tets, {len(bt)} boundary triangles, "
          f"edge multiplicity {dict(tp.check_manifold([tuple(t) for t in bt]))}", flush=True)
    return Pn, tn, bt, np.array(bid)


def envelope(P, tets, btris, ids, a, out, fine=None, planes=None):
    """Cut the named thin walls out by offset envelopes {d_i = delta}, ONE BODY AT A TIME
    (edge-connected wall components), each with its own label. Where two envelopes meet, mmg
    cuts the second iso surface exactly into the first, so the crease is a complete label
    boundary instead of a piecewise-detected 'sharp' chain."""
    names = {int(l.split()[0]): l.split()[1] for l in open(out / "labels.txt")}
    inv = {v: k for k, v in names.items()}
    cut = [inv[n] for n in a.envelope.split(",")]
    B = np.asarray(btris)
    sel = np.isin(ids, cut)
    P_work, ids_work = P, np.asarray(ids)
    if fine is not None:
        # walls from the ORIGINAL dtOO triangulation: the coarse mmgs surface resolves the tightly
        # rounded duct lip only piecewise, its 'sharp' edges split the mantle into 131 faces
        P0, T0, I0 = fine[:3]
        E0 = fine[3] if len(fine) > 3 else None
        P, B, ids, sel = P0, np.asarray(T0), np.asarray(I0), np.isin(I0, cut)
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
    body_idx = [[j for j in comp if comp[j] == i] for i in range(c)]
    bodies = [B[ix] for ix in body_idx]
    body_ent = [E0[ix] if fine is not None and E0 is not None else None for ix in body_idx]
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
    curved_all = np.unique(np.array(sharp_edges(P, [tuple(t) for t in B.tolist()], np.where(sel, 1, 0), 20.0)
                                    or [(0, 0)]).ravel())
    if fine is not None:                       # back to the working (remeshed) mesh for the cuts
        P = P_work
        ids = np.asarray(ids_work)
    new0 = max(names) + 1
    sub_labels, ridge_bodies = [], []
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
        if planes:
            # put plane vertices back after every mmg pass: drift accumulated over all cuts
            # reached 0.039 and the final restoration inverted 283 tets
            Q = restore_planes(np.asarray(P), np.asarray(tets), np.asarray(btris), np.asarray(ids), planes, fatal=False)
            P = Q if Q is not None else P
        if a.ridges:
            # feature curves on this body's envelope at the offsets of its sharp edges; the
            # extra face labels 100+k of ec.ridge_split are renumbered per body
            # faces of this body from its dtOO 2D entities (geometry), not from dihedral detection
            ec.BODY_ENT[0] = body_ent[i]
            P, tets, btris, ids, fmap = ec.ridge_split(np.asarray(P), tets, btris, ids, WV, F, new0 + i, out,
                                                       delta=delta, size_max=a.size_max)
            body_sub = {new0 + i}
            ids = np.asarray(ids).copy()
            for k, lab in fmap.items():
                if lab != new0 + i:
                    ids[ids == lab] = 200 + 10 * i + k
                    sub_labels.append(200 + 10 * i + k)
                    body_sub.add(200 + 10 * i + k)
            ridge_bodies.append((body_sub, new0 + i, F, body_ent[i]))
            tets = np.asarray(tets)
            if planes:
                Q = restore_planes(np.asarray(P), tets, np.asarray(btris), ids, planes, fatal=False)
                P = Q if Q is not None else P
    left = set(ids.tolist()) & set(cut)
    if left:
        raise SystemExit(f"[named] wall labels {left} survived the envelope cuts")
    env = [new0 + i for i in range(c)] + sub_labels
    P, tets, btris, ids = final_coarsen(P, tets, btris, ids, out, a, env, sag, h_iso, planes)
    # keep only the ridge curves over real creases (whole curves), envelopes back to one label each
    extra = []
    for body_sub, main, F, ent in ridge_bodies:
        ids, ex = ec.select_crease_ridges(np.asarray(P), btris, np.asarray(ids), body_sub, main, WV, F, ent, delta)
        extra += ex
    EXTRA_FEAT[0] = extra
    print(f"[named] crease ridge curves kept: {len(extra)} edges")
    with open(out / "labels.txt", "a") as f:
        for i in range(c):
            f.write(f"{new0 + i} ENVELOPE_{i} (body {i} of {a.envelope})\n")
        for l in sub_labels:
            f.write(f"{l} ENVELOPE_{(l - 200) // 10}_face{(l - 200) % 10} (ridge split)\n")
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
    ap.add_argument("--fine-walls", default=None, help="comma-separated names of walls meshed finer")
    ap.add_argument("--wall-h", type=float, default=0.04)
    ap.add_argument("--edge-h", type=float, default=0.02)
    ap.add_argument("--ridge-angle", type=float, default=None, help="mmgs -ar: edges above this stay ridges")
    ap.add_argument("--corner", type=float, default=45.0,
                    help="feature-curve turn (deg) above which a curve vertex becomes a feature vertex")
    ap.add_argument("--delta-gap", type=float, default=0.25,
                    help="envelope distance = this fraction of the narrowest gap around the thin walls")
    ap.add_argument("--delta", type=float, default=None)
    ap.add_argument("--ridges", action="store_true",
                    help="feature curves on every envelope at the offsets of its wall's sharp edges")
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
    ids_fine = ids.copy()
    inv_names = {names.get(l, 'unnamed'): lid[l] for l in labs}

    # isotropic SURFACE remesh with mmgs (refs and ridges kept), then a fresh
    # interior tet mesh from it. Re-meshing the volume instead (dtOO hexes
    # split into 1.9M tets, mmg3d) did not finish in 28 min; gmsh cannot mesh
    # the interior from the original needle-like wall-layer faces either.
    hmap = None
    if a.fine_walls:
        # local size at thin walls: a.wall_h on them, a.edge_h at their strongly curved parts
        # (ring ends, duct lip), growing with distance up to size_max
        Tn = np.asarray(tris)
        fw = [inv_names[n] for n in a.fine_walls.split(",")]
        selw = np.isin(ids, fw)
        curved = np.unique(np.array(sharp_edges(m.points, [tuple(t) for t in Tn.tolist()],
                                                np.where(selw, 1, 0), 15.0) or [(0, 0)]).ravel())
        wv = np.unique(Tn[selw])
        curved = np.intersect1d(curved, wv)
        tw, tc = cKDTree(m.points[wv]), cKDTree(m.points[curved])
        hmap = lambda X: np.clip(np.minimum(a.wall_h + 0.5 * tw.query(X)[0], a.edge_h + 0.5 * tc.query(X)[0]),
                                 a.edge_h, a.size_max)
        print(f"[named] fine walls {a.fine_walls}: {len(wv)} wall vertices, {len(curved)} on curved parts", flush=True)
    S, ST, SR = remesh_surface(m.points, np.asarray(tris), ids, out, a.size_max, a.hausd, a.hgrad,
                               hmap=hmap, ridge_angle=a.ridge_angle)
    nodes = {i: tuple(p) for i, p in enumerate(S)}
    P, tets, idx, used = v2.mesh_interior(nodes, [tuple(t) for t in ST], a.size_max)
    remap = {ot: idx[nt] for nt, ot in zip(range(1, len(used) + 1), used) if nt in idx}
    btris = [[remap[t] for t in tri] for tri in ST]
    ids = SR
    P = np.asarray(P)
    env = []
    planes = planar_labels(P, np.asarray(btris), ids)
    if a.envelope:
        P, tets, btris, ids, env = envelope(P, tets, btris, ids, a, out,
                                            fine=(m.points, np.asarray(tris), np.asarray(ids_fine),
                                                  BOUNDARY_ENT[0]),
                                            planes=planes)
        P = restore_planes(P, np.asarray(tets), np.asarray(btris), ids, planes)
    nd = {i: tuple(p) for i, p in enumerate(P)}
    feat, _ = v2.feature_graph(nd, [tuple(t) for t in btris], ids)
    # sharp edges only on surfaces that are not envelopes (an envelope is smooth by construction;
    # its creases with other envelopes are label boundaries already)
    sh = [e for e in sharp_edges(P, btris, ids, a.sharp) if not set(env) & {int(x) for x in ids[edge_tris(btris, e)]}]
    sh = anchored_chains(extend_chains(P, btris, ids, sh, feat), feat)
    feat = sorted(set(feat) | set(sh) | {tuple(sorted(map(int, e))) for e in EXTRA_FEAT[0]})
    fverts = sorted(set(ec.v5.feature_vertices(feat)) | set(curve_corners(P, feat, a.corner)))
    print(f"[named] feature edges {len(feat)} (sharp inside a label: {len(sh)}), feature vertices {len(fverts)}")
    tp.write_algohex_vtk_analytic(P, tets, btris, ids, feat, fverts, out / "tet.vtk")


if __name__ == "__main__":
    main()
