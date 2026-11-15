"""AlgoHex input whose feature curves come from the GEOMETRY (gmsh 1D entities of the dtOO model),
not from dihedral detection on a mesh.

 1. boundary = 2D elements of the 2D entities that carry a physical name (taken from the physical
    export of the same dtOO model; entity tags are identical)
 2. feature curve = 1D entity whose adjacent boundary entities have different physical names, or
    meet at a crease (mean dihedral along the curve > --crease deg); block seams on a smooth
    surface are dropped
 3. isotropic surface remesh with mmgs; the feature curves are handed over as ridges (Edges +
    Ridges + Corners), so they come back as complete mesh edges
 4. interior tet mesh from the remeshed surface (gmsh, boundary fixed)

usage: prep_from_geometry.py <full.msh (saveAll)> <physical.msh> <out_dir> [--size-max H] [--hausd D]
       [--fine-walls A,B --wall-h h --edge-h h] [--gap-cells N]

 gap sizing: at a boundary vertex the gap is its distance to the nearest boundary of a label that is
 neither its own nor adjacent to it (sharing a mesh edge); h = gap / N there, spread with gradient 0.5.
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

MMGS = ec.MMG.replace("mmg3d_O3", "mmgs_O3")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("full"); ap.add_argument("phys"); ap.add_argument("out_dir")
    ap.add_argument("--size-max", type=float, default=0.3)
    ap.add_argument("--hausd", type=float, default=0.005)
    ap.add_argument("--hgrad", type=float, default=1.5)
    ap.add_argument("--crease", type=float, default=20.0)
    ap.add_argument("--fine-walls", default=None)
    ap.add_argument("--wall-h", type=float, default=0.04)
    ap.add_argument("--edge-h", type=float, default=0.01)
    ap.add_argument("--gap-cells", type=float, default=0,
                    help="cells across narrow gaps between non-adjacent boundary labels (0 = off)")
    ap.add_argument("--gap-h-min", type=float, default=0.01)
    a = ap.parse_args()
    out = Path(a.out_dir); out.mkdir(parents=True, exist_ok=True)

    # physical name per 2D entity tag
    m = meshio.read(a.phys)
    pname = {v[0]: k for k, v in m.field_data.items()}
    g2p = {}
    for c, ph, ge in zip(m.cells, m.cell_data["gmsh:physical"], m.cell_data["gmsh:geometrical"]):
        if c.type in ("triangle", "quad"):
            for p, g in zip(np.asarray(ph), np.asarray(ge)):
                g2p[int(g)] = pname[int(p)]
    labels = sorted(set(g2p.values()))
    lid = {n: i + 1 for i, n in enumerate(labels)}

    gmsh.initialize(); gmsh.option.setNumber("General.Verbosity", 1); gmsh.open(a.full)
    tags, xyz, _ = gmsh.model.mesh.getNodes()
    P = xyz.reshape(-1, 3)
    row = {int(t): i for i, t in enumerate(tags)}
    # boundary = exterior faces of the volume elements, oriented outward (tet_prep helpers); label and
    # 2D entity from the 2D element on that face
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
    nodes0 = {i: tuple(p) for i, p in enumerate(P)}
    tris = tp.orient_and_triangulate(tp.volume_exterior_faces(els), nodes0)
    T = np.array(tris)
    TE = np.array([face_ent.get(frozenset(t), -1) for t in tris])
    if (TE < 0).any():
        raise SystemExit(f"[geo] {(TE < 0).sum()} boundary triangles without a 2D entity")
    TL = np.array([lid[g2p[s]] for s in TE])
    hist = tp.check_manifold([tuple(t) for t in T])
    print(f"[geo] boundary {len(T)} triangles from {len(set(TE.tolist()))} entities, labels {len(labels)}, "
          f"edge multiplicity {dict(hist)}", flush=True)
    # feature curves from 1D entities
    e2t = collections.defaultdict(list)
    for i, tr in enumerate(T):
        for x, y in ((tr[0], tr[1]), (tr[1], tr[2]), (tr[2], tr[0])):
            e2t[tuple(sorted((int(x), int(y))))].append(i)
    nrm = np.cross(P[T[:, 1]] - P[T[:, 0]], P[T[:, 2]] - P[T[:, 0]])
    nrm /= np.linalg.norm(nrm, axis=1)[:, None]
    # 1D entities -> geometric LINES: dtOO splits one edge into several curve entities (canada ring 3:
    # curve 230 (8 -> 30 deg) continues as curve 162 (8 deg); judging entities one by one dropped 162
    # and left a loose end). Join entities at shared end points where they continue tangentially
    # (mutual best match, turn < 30 deg), then decide per line.
    cur = {}
    for _, c in gmsh.model.getEntities(1):
        et, _, en = gmsh.model.mesh.getElements(1, c)
        segs = [tuple(sorted((row[int(x)], row[int(y)]))) for t, nd in zip(et, en) for x, y in np.asarray(nd).reshape(-1, 2)]
        segs = [s_ for s_ in segs if len(e2t.get(s_, [])) == 2]
        if segs:
            cur[c] = segs
    def ends_dir(segs):
        deg_ = collections.Counter(v for s_ in segs for v in s_)
        out = {}
        for v in (v for v, k in deg_.items() if k == 1):
            w = next(u for s_ in segs if v in s_ for u in s_ if u != v)
            d = P[w] - P[v]; out[v] = d / np.linalg.norm(d)
        return out
    ends = {c: ends_dir(segs) for c, segs in cur.items()}
    # adjacent surface labels of a curve: only curves between the same labels can be one line
    # (a seam on the symmetry plane continued a ring edge tangentially and was swallowed into it)
    clab = {c: frozenset(int(TL[i]) for s_ in segs for i in e2t[s_]) for c, segs in cur.items()}
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
                if i < j:
                    turn = np.degrees(np.arccos(np.clip(-di @ dj, -1, 1)))   # outgoing directions
                    if turn < 30 and clab[ci] == clab[cj]:
                        for x, y in ((i, j), (j, i)):
                            if x not in best or turn < best[x][0]:
                                best[x] = (turn, y)
        for x, (turn, y) in best.items():
            if best.get(y, (None, None))[1] == x:
                par[find(lst[x][0])] = find(lst[y][0])
    lines = collections.defaultdict(list)
    for c in cur:
        lines[find(c)].append(c)
    feat, kinds = [], collections.Counter()
    for root, cs in lines.items():
        segs = [s_ for c in cs for s_ in cur[c]]
        labs = {int(TL[i]) for s_ in segs for i in e2t[s_]}
        ang = max(np.degrees(np.arccos(np.clip(nrm[e2t[s_][0]] @ nrm[e2t[s_][1]], -1, 1))) for s_ in segs)
        if len(labs) > 1 or ang > a.crease:
            feat += segs
            kinds["label boundary" if len(labs) > 1 else "crease"] += 1
        else:
            kinds["smooth seam (dropped)"] += 1
        if len(cs) > 1 and ang > a.crease and len(labs) == 1:
            print(f"[geo] line from entities {sorted(cs)}: max dihedral {ang:.1f} -> crease")
    gmsh.finalize()
    print(f"[geo] {len(cur)} curve entities -> {len(lines)} lines: {dict(kinds)}; feature segments {len(feat)}", flush=True)

    # mmgs surface remesh with the feature curves as ridges
    u = np.unique(T); mp = -np.ones(len(P), int); mp[u] = np.arange(len(u))
    F = np.array(feat)
    val = collections.Counter(int(x) for e in F for x in e)
    corners = sorted(v for v, c in val.items() if c != 2)
    with open(out / "surf_in.mesh", "w") as f:
        f.write(f"MeshVersionFormatted 2\nDimension 3\nVertices\n{len(u)}\n")
        for p in P[u]:
            f.write(f"{p[0]:.17g} {p[1]:.17g} {p[2]:.17g} 0\n")
        f.write(f"Triangles\n{len(T)}\n")
        for t, l in zip(mp[T], TL):
            f.write(f"{t[0]+1} {t[1]+1} {t[2]+1} {l}\n")
        f.write(f"Edges\n{len(F)}\n")
        for e in mp[F]:
            f.write(f"{e[0]+1} {e[1]+1} 1\n")
        f.write(f"Ridges\n{len(F)}\n" + "".join(f"{i+1}\n" for i in range(len(F))))
        f.write(f"Corners\n{len(corners)}\n" + "".join(f"{mp[v]+1}\n" for v in corners))
        f.write("End\n")
    cmd = [MMGS, "-in", str(out / "surf_in.mesh"), "-hmax", str(a.size_max), "-hausd", str(a.hausd),
           "-hgrad", str(a.hgrad), "-nr", "-out", str(out / "surf_out.mesh")]
    h = np.full(len(u), a.size_max)
    if a.fine_walls:
        fw = {lid[n] for n in a.fine_walls.split(",")}
        wv = np.unique(T[np.isin(TL, list(fw))])
        cv = np.intersect1d(np.unique(F), wv)
        tw, tc = cKDTree(P[wv]), cKDTree(P[cv])
        h = np.clip(np.minimum(a.wall_h + 0.5 * tw.query(P[u])[0], a.edge_h + 0.5 * tc.query(P[u])[0]),
                    a.edge_h, a.size_max)
    if a.gap_cells:
        adj = collections.defaultdict(set)
        for ts_ in e2t.values():
            ls_ = {int(TL[i]) for i in ts_}
            for x in ls_:
                adj[x] |= ls_
        vl = collections.defaultdict(set)
        for tr, l in zip(T, TL):
            for v in tr:
                vl[int(v)].add(int(l))
        gap = np.full(len(P), np.inf)
        for L in sorted(set(TL.tolist())):
            mine = np.array([v for v, ls_ in vl.items() if L in ls_])
            other = np.array([v for v, ls_ in vl.items() if not (ls_ & adj[L])])
            if len(other):
                gap[mine] = np.minimum(gap[mine], cKDTree(P[other]).query(P[mine])[0])
        hg = np.clip(gap / a.gap_cells, a.gap_h_min, a.size_max)
        src = np.where(hg < a.size_max)[0]
        if len(src):
            k = min(16, len(src))
            d, j = cKDTree(P[src]).query(P[u], k=k)
            d, j = d.reshape(len(u), -1), j.reshape(len(u), -1)
            h = np.minimum(h, (hg[src][j] + 0.5 * d).min(1))
            print(f"[geo] gap sizing: {len(src)} vertices in gaps < {a.gap_cells * a.size_max:.2f}, "
                  f"smallest gap {gap[src].min():.4f} -> h {hg[src].min():.4f}", flush=True)
    if a.fine_walls or a.gap_cells:
        ec.write_sol(out / "surf_in.sol", h)
        cmd[3:3] = ["-sol", str(out / "surf_in.sol")]
    with open(out / "mmgs.log", "w") as lf:
        if subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT).returncode:
            raise SystemExit("[geo] mmgs failed")
    M = ec.read_mesh(out / "surf_out.mesh")
    S = M["Vertices"][:, :3]; ST = M["Triangles"][:, :3].astype(int) - 1; SR = M["Triangles"][:, 3].astype(int)
    SE = M["Edges"][:, :2].astype(int) - 1 if "Edges" in M else np.zeros((0, 2), int)
    print(f"[geo] surface remeshed: {len(T)} -> {len(ST)} triangles, {len(SE)} feature edges kept, "
          f"edge multiplicity {dict(tp.check_manifold([tuple(t) for t in ST]))}", flush=True)
    nodes = {i: tuple(p) for i, p in enumerate(S)}
    P2, tets, idx, used = v2.mesh_interior(nodes, [tuple(t) for t in ST], a.size_max)
    remap = {ot: idx[nt] for nt, ot in zip(range(1, len(used) + 1), used) if nt in idx}
    bt = [[remap[t] for t in tri] for tri in ST]
    fe = [(remap[int(x)], remap[int(y)]) for x, y in SE]
    # label boundaries are feature curves too (they are in SE already if they were 1D entities)
    nd = {i: tuple(p) for i, p in enumerate(P2)}
    lab_feat, _ = v2.feature_graph(nd, [tuple(t) for t in bt], SR)
    fe = sorted({tuple(sorted(e)) for e in fe} | {tuple(sorted(e)) for e in lab_feat})
    fv = ec.v5.feature_vertices(fe)
    tp.write_algohex_vtk_analytic(P2, tets, bt, SR, fe, fv, out / "tet.vtk")
    with open(out / "labels.txt", "w") as f:
        for n in labels:
            f.write(f"{lid[n]} {n}\n")


if __name__ == "__main__":
    main()
