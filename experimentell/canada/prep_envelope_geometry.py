"""AlgoHex input: named thin walls cut out by offset envelopes, ALL feature curves from the geometry.

 1. boundary, labels and geometric LINES exactly as prep_from_geometry (1D entities joined where they
    continue tangentially between the same labels; a line is a feature if it separates labels or
    creases somewhere above --crease deg)
 2. feature lines on the walls that get cut out are set aside as CREASE LINES of their body
 3. surface remesh (mmgs) and interior tet mesh; every later mmg3d pass gets the feature curves as
    required ridges (Edges + Ridges, -nr), so creases between 20 and 45 deg are not smoothed away
 4. one envelope per body (edge-connected wall component): level set d_body - delta, pre-refined,
    each with its own label; symmetry/planar labels restored after every pass
 5. ridge curves: for every crease line L of a body, phi_L(x) = (x - c_L(x)) . m(c), c_L the closest
    point on L and m the unit (n_A - n_B) of the two wall faces at L, evaluated only in the band of
    points whose closest wall point lies on L. All lines in ONE mmg -ls split (no insert/move/swap),
    after a pre-refinement around the ridges. The ridge = envelope edges with |phi| ~ 0 in the band:
    one chain per crease line, no relabelling, no detection.

usage: prep_envelope_geometry.py <full.msh> <physical.msh> <out_dir> --envelope A,B [--delta-gap 0.25]
"""
import argparse, collections, subprocess, sys
from pathlib import Path

import numpy as np
import meshio
import gmsh
from scipy.spatial import cKDTree

HEX = Path("/opt/stack/domain_partition_3D/experimentell/hex3d_algohex")
sys.path.insert(0, str(HEX.parent.parent)); sys.path.insert(0, str(HEX))
import tet_prep as tp                                       # noqa: E402
import tet_prep_v2 as v2                                    # noqa: E402
import envelope_cut as ec                                   # noqa: E402

MMG, MMGS = ec.MMG, ec.MMG.replace("mmg3d_O3", "mmgs_O3")
FACES4 = ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))


# ------------------------------------------------------------------ mesh io with ridges
def write_mesh(path, P, tets, tris, refs, edges=()):
    edges = list(edges)
    with open(path, "w") as f:
        f.write(f"MeshVersionFormatted 2\nDimension 3\nVertices\n{len(P)}\n")
        for p in P:
            f.write(f"{p[0]:.17g} {p[1]:.17g} {p[2]:.17g} 0\n")
        f.write(f"Triangles\n{len(tris)}\n")
        for t, r in zip(tris, refs):
            f.write(f"{t[0]+1} {t[1]+1} {t[2]+1} {int(r)}\n")
        if edges:
            f.write(f"Edges\n{len(edges)}\n" + "".join(f"{a+1} {b+1} 1\n" for a, b in edges))
            f.write(f"Ridges\n{len(edges)}\n" + "".join(f"{i+1}\n" for i in range(len(edges))))
            val = collections.Counter(v for e in edges for v in e)
            cor = [v for v, c in val.items() if c != 2]
            f.write(f"Corners\n{len(cor)}\n" + "".join(f"{v+1}\n" for v in cor))
        if len(tets):
            f.write(f"Tetrahedra\n{len(tets)}\n")
            for t in tets:
                f.write(f"{t[0]+1} {t[1]+1} {t[2]+1} {t[3]+1} 1\n")
        f.write("End\n")


def read_back(path, keep_ref=None):
    """mmg3d output -> (P, tets, boundary tris oriented outward, refs, edges); only tets of keep_ref"""
    M = ec.read_mesh(path)
    tet = M["Tetrahedra"]
    if keep_ref is not None:
        tet = tet[tet[:, 4] == keep_ref]
    tet = tet[:, :4].astype(int) - 1
    used = np.unique(tet)
    newi = -np.ones(len(M["Vertices"]), int); newi[used] = np.arange(len(used))
    Pn = M["Vertices"][used, :3]; tn = newi[tet]
    cnt = collections.Counter(frozenset(t[list(f)]) for t in tn for f in FACES4)
    tri_ref = {frozenset(newi[t[:3].astype(int) - 1]): int(t[3]) for t in M["Triangles"]}
    bt, bid = [], []
    for t in tn:
        c = Pn[t].mean(0)
        for f in FACES4:
            k = frozenset(t[list(f)])
            if cnt[k] == 1:
                ff = t[list(f)]
                n = np.cross(Pn[ff[1]] - Pn[ff[0]], Pn[ff[2]] - Pn[ff[0]])
                if np.dot(n, Pn[ff].mean(0) - c) < 0:
                    ff = ff[::-1]
                bt.append(list(ff)); bid.append(tri_ref.get(k, ec.ISO_REF))
    edges = []
    if "Edges" in M:
        for e in M["Edges"][:, :2].astype(int) - 1:
            a, b = newi[e[0]], newi[e[1]]
            if a >= 0 and b >= 0:
                edges.append((int(a), int(b)))
    return Pn, tn, np.array(bt), np.array(bid), edges


def run(cmd, log):
    with open(log, "w") as lf:
        if subprocess.run([str(c) for c in cmd], stdout=lf, stderr=subprocess.STDOUT).returncode:
            raise SystemExit(f"[envg] failed: {' '.join(map(str, cmd))} (see {log})")


# ------------------------------------------------------------------ geometry lines
def geometry(full, phys, crease):
    m = meshio.read(phys)
    pname = {v[0]: k for k, v in m.field_data.items()}
    g2p = {}
    for c, ph, ge in zip(m.cells, m.cell_data["gmsh:physical"], m.cell_data["gmsh:geometrical"]):
        if c.type in ("triangle", "quad"):
            for p, g in zip(np.asarray(ph), np.asarray(ge)):
                g2p[int(g)] = pname[int(p)]
    labels = sorted(set(g2p.values()))
    lid = {n: i + 1 for i, n in enumerate(labels)}
    gmsh.initialize(); gmsh.option.setNumber("General.Verbosity", 1); gmsh.open(full)
    tags, xyz, _ = gmsh.model.mesh.getNodes()
    P = xyz.reshape(-1, 3); row = {int(t): i for i, t in enumerate(tags)}
    face_ent = {}
    for _, s in gmsh.model.getEntities(2):
        if s not in g2p:
            continue
        et, _, en = gmsh.model.mesh.getElements(2, s)
        for t, nd in zip(et, en):
            k = gmsh.model.mesh.getElementProperties(t)[3]
            for e in np.asarray(nd).reshape(-1, k):
                e = [row[int(x)] for x in e]
                face_ent[frozenset(e)] = s
                if k == 4:
                    for q in ((e[0], e[1], e[2]), (e[0], e[2], e[3]), (e[0], e[1], e[3]), (e[1], e[2], e[3])):
                        face_ent[frozenset(q)] = s
    els = []
    for _, v in gmsh.model.getEntities(3):
        et, _, en = gmsh.model.mesh.getElements(3, v)
        for t, nd in zip(et, en):
            k = gmsh.model.mesh.getElementProperties(t)[3]
            for e in np.asarray(nd).reshape(-1, k):
                els.append((int(t), 0, v, [row[int(x)] for x in e]))
    tris = tp.orient_and_triangulate(tp.volume_exterior_faces(els), {i: tuple(p) for i, p in enumerate(P)})
    T = np.array(tris)
    TE = np.array([face_ent[frozenset(t)] for t in tris])
    TL = np.array([lid[g2p[s]] for s in TE])
    e2t = collections.defaultdict(list)
    for i, tr in enumerate(T):
        for x, y in ((tr[0], tr[1]), (tr[1], tr[2]), (tr[2], tr[0])):
            e2t[tuple(sorted((int(x), int(y))))].append(i)
    nrm = np.cross(P[T[:, 1]] - P[T[:, 0]], P[T[:, 2]] - P[T[:, 0]])
    nrm /= np.linalg.norm(nrm, axis=1)[:, None]
    cur = {}
    for _, c in gmsh.model.getEntities(1):
        et, _, en = gmsh.model.mesh.getElements(1, c)
        segs = [tuple(sorted((row[int(x)], row[int(y)]))) for t, nd in zip(et, en) for x, y in np.asarray(nd).reshape(-1, 2)]
        segs = [s_ for s_ in segs if len(e2t.get(s_, [])) == 2]
        if segs:
            cur[c] = segs
    gmsh.finalize()
    clab = {c: frozenset(int(TL[i]) for s_ in segs for i in e2t[s_]) for c, segs in cur.items()}
    ends = {}
    for c, segs in cur.items():
        deg_ = collections.Counter(v for s_ in segs for v in s_)
        out = {}
        for v in (v for v, k in deg_.items() if k == 1):
            w = next(u for s_ in segs if v in s_ for u in s_ if u != v)
            d = P[w] - P[v]; out[v] = d / np.linalg.norm(d)
        ends[c] = out
    at = collections.defaultdict(list)
    for c, ed in ends.items():
        for v, d in ed.items():
            at[v].append((c, d))
    par = {c: c for c in cur}
    def find(c):
        while par[c] != c:
            c = par[c]
        return c
    for v, lst in at.items():
        best = {}
        for i, (ci, di) in enumerate(lst):
            for j, (cj, dj) in enumerate(lst):
                if i < j and clab[ci] == clab[cj]:
                    turn = np.degrees(np.arccos(np.clip(-di @ dj, -1, 1)))
                    if turn < 30:
                        for x, y in ((i, j), (j, i)):
                            if x not in best or turn < best[x][0]:
                                best[x] = (turn, y)
        for x, (turn, y) in best.items():
            if best.get(y, (None, None))[1] == x:
                par[find(lst[x][0])] = find(lst[y][0])
    lines = collections.defaultdict(list)
    for c in cur:
        lines[find(c)].append(c)
    feature_lines = []        # (segments, kind, labels)
    for root, cs in lines.items():
        segs = [s_ for c in cs for s_ in cur[c]]
        labs = {int(TL[i]) for s_ in segs for i in e2t[s_]}
        ang = max(np.degrees(np.arccos(np.clip(nrm[e2t[s_][0]] @ nrm[e2t[s_][1]], -1, 1))) for s_ in segs)
        if len(labs) > 1:
            feature_lines.append((segs, "label", labs))
        elif ang > crease:
            feature_lines.append((segs, "crease", labs))
    return P, T, TL, TE, labels, lid, feature_lines, e2t, nrm


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("full"); ap.add_argument("phys"); ap.add_argument("out_dir")
    ap.add_argument("--envelope", required=True)
    ap.add_argument("--delta-gap", type=float, default=0.25)
    ap.add_argument("--size-max", type=float, default=0.3)
    ap.add_argument("--hausd", type=float, default=0.01)
    ap.add_argument("--hgrad", type=float, default=2.0)
    ap.add_argument("--crease", type=float, default=20.0)
    a = ap.parse_args()
    out = Path(a.out_dir); out.mkdir(parents=True, exist_ok=True)

    P0, T0, L0, E0, labels, lid, flines, e2t0, nrm0 = geometry(a.full, a.phys, a.crease)
    cut = {lid[n] for n in a.envelope.split(",")}
    keep_lines, crease_lines = [], []
    for segs, kind, labs in flines:
        if labs <= cut:
            if kind == "crease":
                crease_lines.append(segs)        # creases of walls that get cut out
        elif not (labs & cut):
            keep_lines.append(segs)              # curves of surfaces that stay
        # curves between a wall and a staying surface vanish under the envelope; the envelope's own
        # intersection with that surface becomes the new label boundary
    print(f"[envg] {len(flines)} feature lines: {len(keep_lines)} kept on staying surfaces, "
          f"{len(crease_lines)} crease lines on the cut walls", flush=True)

    # bodies, gaps, delta (on the fine dtOO wall)
    wsel = np.isin(L0, list(cut))
    wt = np.flatnonzero(wsel)
    adj = collections.defaultdict(set)
    ed = collections.defaultdict(list)
    for i in wt:
        t = T0[i]
        for x, y in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
            ed[tuple(sorted((int(x), int(y))))].append(i)
    for inc in ed.values():
        for x in inc:
            adj[x].update(inc)
    comp, nb = {}, 0
    for s0 in wt:
        if s0 in comp:
            continue
        st = [s0]; comp[s0] = nb
        while st:
            u = st.pop()
            for w in adj[u]:
                if w not in comp:
                    comp[w] = nb; st.append(w)
        nb += 1
    bodies = [T0[[j for j in comp if comp[j] == i]] for i in range(nb)]
    gaps = [ec.surface_distance(P0, bodies[i], P0[np.unique(bodies[j])]).min()
            for i in range(nb) for j in range(i + 1, nb)]
    wall_v = set(np.unique(T0[wsel]).tolist())
    for l in set(L0.tolist()) - cut:
        lv = np.unique(T0[L0 == l])
        if not wall_v.intersection(lv.tolist()):
            gaps.append(ec.surface_distance(P0, T0[wsel], P0[lv]).min())
    delta = a.delta_gap * min(gaps)
    h_iso = delta * np.radians(15)
    sag = delta * (1 - np.cos(np.radians(15)))
    body_of_line = []
    for segs in crease_lines:
        v = segs[0][0]
        body_of_line.append(next(i for i in range(nb) if v in set(np.unique(bodies[i]).tolist())))
    print(f"[envg] {nb} bodies, gaps {np.round(sorted(gaps)[:4], 3)}, delta {delta:.4f}; crease lines per body "
          f"{dict(collections.Counter(body_of_line))}", flush=True)

    # surface remesh with mmgs, the kept feature curves as required ridges; fine near the walls
    keep_segs = sorted({s_ for segs in keep_lines for s_ in segs})
    # the surface remesh gets ALL feature lines (incl. the knife edges of the walls): without them
    # mmgs smooths an 8-deg wedge and its two sides cross (gmsh: "segment and facet intersect")
    all_segs = sorted({s_ for segs, _, _ in flines for s_ in segs})
    u = np.unique(T0); mp = -np.ones(len(P0), int); mp[u] = np.arange(len(u))
    write_mesh(out / "surf_in.mesh", P0[u], [], mp[T0], L0, [tuple(mp[list(e)]) for e in all_segs])
    # same surface sizing as the input that meshed cleanly (geo4): walls 0.04, their curves 0.01,
    # hausd 0.005; coarser settings made mmgs produce a self-intersecting surface in the slots
    tw = cKDTree(P0[np.unique(T0[wsel])])
    cv = np.unique([v for segs in crease_lines for s_ in segs for v in s_] or [int(T0[wsel][0, 0])])
    tc = cKDTree(P0[cv])
    hs = np.clip(np.minimum(0.04 + 0.5 * tw.query(P0[u])[0], 0.01 + 0.5 * tc.query(P0[u])[0]), 0.01, a.size_max)
    ec.write_sol(out / "surf_in.sol", hs)
    run([MMGS, "-in", out / "surf_in.mesh", "-sol", out / "surf_in.sol", "-hausd", 0.005, "-hgrad", 1.5,
         "-nr", "-out", out / "surf_out.mesh"], out / "mmgs.log")
    M = ec.read_mesh(out / "surf_out.mesh")
    S = M["Vertices"][:, :3]; ST = M["Triangles"][:, :3].astype(int) - 1; SR = M["Triangles"][:, 3].astype(int)
    SE = [tuple(e) for e in (M["Edges"][:, :2].astype(int) - 1)] if "Edges" in M else []
    print(f"[envg] surface: {len(ST)} triangles, {len(SE)} required curve edges", flush=True)
    nodes = {i: tuple(p) for i, p in enumerate(S)}
    P, tets, idx, used = v2.mesh_interior(nodes, [tuple(t) for t in ST], a.size_max)
    P = np.asarray(P); tets = np.asarray(tets)
    remap = {ot: idx[nt] for nt, ot in zip(range(1, len(used) + 1), used) if nt in idx}
    btris = np.array([[remap[t] for t in tri] for tri in ST]); ids = SR.copy()
    edges = [(remap[x], remap[y]) for x, y in SE]
    # mmg3d passes: only curves of surfaces that stay (both incident labels outside the cut walls)
    e2l_ = collections.defaultdict(set)
    for t, l in zip(btris, ids):
        for x, y in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
            e2l_[tuple(sorted((int(x), int(y))))].add(int(l))
    edges = [e for e in edges if not (e2l_[tuple(sorted(e))] & cut)]

    # planar labels (restored after every pass)
    planes = {}
    span = np.linalg.norm(np.ptp(P, 0))
    for l in np.unique(ids):
        X = P[np.unique(btris[ids == l])]
        c = X.mean(0); n = np.linalg.svd(X - c, full_matrices=False)[2][-1]
        if np.abs((X - c) @ n).max() < 1e-9 * span:
            planes[int(l)] = (n, float(n @ c))
    def restore(P, tets, btris, ids):
        on = collections.defaultdict(set)
        for t, l in zip(btris, ids):
            if int(l) in planes:
                for v in t:
                    on[int(v)].add(int(l))
        Q = P.copy()
        for v, ls in on.items():
            A = np.array([planes[l][0] for l in ls]); b = np.array([planes[l][1] for l in ls])
            Q[v] = P[v] + A.T @ np.linalg.lstsq(A @ A.T, b - A @ P[v], rcond=None)[0]
        vol = np.einsum("ij,ij->i", np.cross(Q[tets[:, 1]] - Q[tets[:, 0]], Q[tets[:, 2]] - Q[tets[:, 0]]),
                        Q[tets[:, 3]] - Q[tets[:, 0]])
        if (vol <= 0).any():
            print(f"[envg] plane restoration would invert {(vol <= 0).sum()} tets, skipped", flush=True)
            return P
        return Q

    def local_params(path, fine_refs, plane_refs):
        rows = [f"{r} Triangles {h_iso:.6g} {a.size_max:.6g} {sag:.6g}" for r in fine_refs]
        rows += [f"{r} Triangles {a.size_max / 10:.6g} {a.size_max:.6g} 1e-04" for r in plane_refs]
        with open(path, "w") as f:
            f.write(f"Parameters\n{len(rows)}\n" + "\n".join(rows) + "\n")

    # envelopes, one body at a time
    new0 = max(lid.values()) + 1
    for i, F in enumerate(bodies):
        lab = new0 + i
        ls = lambda X, F=F: ec.surface_distance(P0, F, X) - delta
        tag = f"env{i}"
        shift = 1000 if ec.ISO_REF in set(ids.tolist()) else 0
        # pre-refine around the iso surface
        h = np.clip(h_iso + 0.5 * np.abs(ls(P)), h_iso, a.size_max)
        write_mesh(out / f"{tag}_pre.mesh", P, tets, btris, ids + shift, edges)
        ec.write_sol(out / f"{tag}_pre.sol", h)
        local_params(out / f"{tag}_pre.mmg3d", [r + shift for r in range(new0, lab)], [r + shift for r in planes])
        run([MMG, "-in", out / f"{tag}_pre.mesh", "-sol", out / f"{tag}_pre.sol", "-hausd", a.hausd,
             "-hgrad", a.hgrad, "-nr", "-out", out / f"{tag}_pre_out.mesh"], out / f"{tag}_pre.log")
        P, tets, btris, ids, edges = read_back(out / f"{tag}_pre_out.mesh")
        ids = ids - shift
        # cut
        write_mesh(out / f"{tag}_in.mesh", P, tets, btris, ids + shift, edges)
        ec.write_sol(out / f"{tag}_in.sol", ls(P))
        local_params(out / f"{tag}_in.mmg3d", [ec.ISO_REF] + [r + shift for r in range(new0, lab)],
                     [r + shift for r in planes])
        run([MMG, "-in", out / f"{tag}_in.mesh", "-sol", out / f"{tag}_in.sol", "-ls", 0, "-hmax", a.size_max,
             "-hausd", a.hausd, "-hgrad", a.hgrad, "-nr", "-out", out / f"{tag}_cut.mesh"], out / f"{tag}_cut.log")
        Mc = ec.read_mesh(out / f"{tag}_cut.mesh")
        tet = Mc["Tetrahedra"]
        cen = Mc["Vertices"][tet[:, :4].astype(int) - 1, :3].mean(1)
        cls = ls(cen)
        keep_ref = max(set(tet[:, 4].astype(int)), key=lambda r: np.median(cls[tet[:, 4] == r]))
        P, tets, btris, ids, edges = read_back(out / f"{tag}_cut.mesh", keep_ref)
        ids = np.where(ids == ec.ISO_REF, lab + shift, ids) - shift
        left = set(ids.tolist()) & cut
        P = restore(P, tets, btris, ids)
        print(f"[envg] body {i}: {len(tets)} tets, wall labels left {sorted(left)}", flush=True)
    env_labels = list(range(new0, new0 + nb))
    if set(ids.tolist()) & cut:
        raise SystemExit("[envg] wall labels survived the envelope cuts")

    # final coarsening, kept curves required, envelopes fine, planes frozen
    write_mesh(out / "final_in.mesh", P, tets, btris, ids, edges)
    local_params(out / "final_in.mmg3d", env_labels, list(planes))
    run([MMG, "-in", out / "final_in.mesh", "-hmax", a.size_max, "-hausd", a.hausd, "-hgrad", a.hgrad, "-nr",
         "-out", out / "final_out.mesh"], out / "final.log")
    P, tets, btris, ids, edges = read_back(out / "final_out.mesh")
    P = restore(P, tets, btris, ids)
    print(f"[envg] after final pass: {len(P)} nodes, {len(tets)} tets, {len(edges)} required curve edges", flush=True)

    # ridge curves over the crease lines
    ridge_edges = []
    if crease_lines:
        Bmid, Bm, Bline, Bbody = [], [], [], []
        for li, segs in enumerate(crease_lines):
            for s_ in segs:
                inc = e2t0[s_]
                n_a, n_b = nrm0[inc[0]], nrm0[inc[1]]
                m = n_a - n_b
                m = m / max(np.linalg.norm(m), 1e-300)
                Bmid.append(s_); Bm.append(m); Bline.append(li); Bbody.append(body_of_line[li])
        segA = np.array([P0[s_[0]] for s_ in Bmid]); segB = np.array([P0[s_[1]] for s_ in Bmid])
        Bm = np.array(Bm); Bline = np.array(Bline)
        segtree = cKDTree(0.5 * (segA + segB))

        def line_phi(X):
            """phi for every point: (x - c) . m on the crease line whose band contains x, else +1"""
            phi = np.ones(len(X)); which = -np.ones(len(X), int)
            _, cand = segtree.query(X, k=min(8, len(segA)))
            best = np.full(len(X), np.inf); bc = np.zeros_like(X); bj = np.zeros(len(X), int)
            for j in range(cand.shape[1]):
                k = cand[:, j]
                A, B = segA[k], segB[k]
                t = np.clip(((X - A) * (B - A)).sum(1) / np.maximum(((B - A) ** 2).sum(1), 1e-300), 0, 1)
                c = A + t[:, None] * (B - A)
                d = np.linalg.norm(X - c, axis=1)
                mk = d < best
                best[mk], bc[mk], bj[mk] = d[mk], c[mk], k[mk]
            # band: the closest point on the body's wall IS on this line (within a small tolerance)
            dw = np.full(len(X), np.inf)
            for b in set(body_of_line):
                sel = np.array([Bbody[j] == b for j in bj])
                if sel.any():
                    dw[sel] = ec.surface_distance(P0, bodies[b], X[sel])
            band = (best - dw < 0.05 * delta) & (best < 2 * delta)
            phi[band] = ((X[band] - bc[band]) * Bm[bj[band]]).sum(1)
            which[band] = Bline[bj[band]]
            return phi, which

        phi_v, _ = line_phi(P)
        dE = np.full(len(P), np.inf)
        for b in range(nb):
            dE = np.minimum(dE, np.abs(ec.surface_distance(P0, bodies[b], P) - delta))
        h = np.clip(h_iso + 0.5 * (np.abs(phi_v) + dE), h_iso, a.size_max)
        h[phi_v == 1] = np.clip(h_iso + 0.5 * dE[phi_v == 1] + 0.5 * delta, h_iso, a.size_max)
        write_mesh(out / "ridge_pre.mesh", P, tets, btris, ids, edges)
        ec.write_sol(out / "ridge_pre.sol", h)
        local_params(out / "ridge_pre.mmg3d", env_labels, list(planes))
        run([MMG, "-in", out / "ridge_pre.mesh", "-sol", out / "ridge_pre.sol", "-hausd", a.hausd, "-hgrad", a.hgrad,
             "-nr", "-out", out / "ridge_pre_out.mesh"], out / "ridge_pre.log")
        P, tets, btris, ids, edges = read_back(out / "ridge_pre_out.mesh")
        P = restore(P, tets, btris, ids)
        phi_v, which = line_phi(P)
        shift = 1000 if ec.ISO_REF in set(ids.tolist()) else 0
        write_mesh(out / "ridge_in.mesh", P, tets, btris, ids + shift, edges)
        ec.write_sol(out / "ridge_in.sol", phi_v)
        run([MMG, "-in", out / "ridge_in.mesh", "-sol", out / "ridge_in.sol", "-ls", 0, "-noinsert", "-noswap",
             "-nomove", "-nr", "-out", out / "ridge_out.mesh"], out / "ridge.log")
        Mr = ec.read_mesh(out / "ridge_out.mesh")
        Pn = Mr["Vertices"][:, :3]; tn = Mr["Tetrahedra"][:, :4].astype(int) - 1
        cnt = collections.Counter(frozenset(t[list(f)]) for t in tn for f in FACES4)
        ref = {frozenset(t[:3].astype(int) - 1): int(t[3]) for t in Mr["Triangles"]}
        bt, bid = [], []
        for t in tn:
            c = Pn[t].mean(0)
            for f in FACES4:
                k = frozenset(t[list(f)])
                if cnt[k] == 1:
                    ff = t[list(f)]
                    nn = np.cross(Pn[ff[1]] - Pn[ff[0]], Pn[ff[2]] - Pn[ff[0]])
                    if np.dot(nn, Pn[ff].mean(0) - c) < 0:
                        ff = ff[::-1]
                    bt.append(list(ff)); bid.append(ref.get(k, -1) - shift)
        P, tets, btris, ids = Pn, tn, np.array(bt), np.array(bid)
        edges = [tuple(e) for e in (Mr["Edges"][:, :2].astype(int) - 1)] if "Edges" in Mr else []
        phi_v, which = line_phi(P)
        hloc = np.median(np.linalg.norm(P[btris[:, 1]] - P[btris[:, 0]], axis=1))
        on = np.abs(phi_v) < 0.02 * h_iso
        env_set = set(env_labels)
        e2l = collections.defaultdict(set)
        for t, l in zip(btris, ids):
            for x, y in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
                e2l[tuple(sorted((int(x), int(y))))].add(int(l))
        for e, ls_ in e2l.items():
            if ls_ <= env_set and on[e[0]] and on[e[1]] and which[e[0]] == which[e[1]] >= 0:
                ridge_edges.append(e)
        # one chain per crease line?
        g = collections.defaultdict(list)
        for x, y in ridge_edges:
            g[x].append(y); g[y].append(x)
        seen, pieces = set(), collections.Counter()
        for s0 in g:
            if s0 in seen:
                continue
            st = [s0]; seen.add(s0)
            while st:
                uu = st.pop()
                for w in g[uu]:
                    if w not in seen:
                        seen.add(w); st.append(w)
            pieces[int(which[s0])] += 1
        print(f"[envg] ridge curves: {len(ridge_edges)} edges; pieces per crease line {dict(pieces)} "
              f"(of {len(crease_lines)} lines)", flush=True)

    # features: required curves + label boundaries + ridges
    nd = {i: tuple(p) for i, p in enumerate(P)}
    lab_feat, _ = v2.feature_graph(nd, [tuple(t) for t in btris], ids)
    feat = sorted({tuple(sorted(map(int, e))) for e in edges} | {tuple(sorted(e)) for e in lab_feat}
                  | {tuple(sorted(e)) for e in ridge_edges})
    fv = ec.v5.feature_vertices(feat)
    tp.write_algohex_vtk_analytic(P, tets, btris, ids, feat, fv, out / "tet.vtk")
    with open(out / "labels.txt", "w") as f:
        for n in labels:
            f.write(f"{lid[n]} {n}\n")
        for i in range(nb):
            f.write(f"{new0 + i} ENVELOPE_{i}\n")


if __name__ == "__main__":
    main()
