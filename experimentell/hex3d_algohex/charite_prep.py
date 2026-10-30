"""Stage 1 for the Charite left-ventricle (LV) geometries: STL -> AlgoHex input.

The tistos route starts from a dtOO gmsh grid (`tet_prep_v5.py`); the LV cases
come as three surface STLs per case directory:

    LV_inner.stl       endocardium   (label 1)
    LV_outer.stl       epicardium    (label 2)
    LV_connection.stl  basal ring that closes the wall at the valve plane (label 3)

The surface is `LV_all.stl` (closed 2-manifold, exactly the union of the
three parts); the parts only supply the labels, matched per triangle on its
vertex positions. The domain is the myocardium between inner and outer wall.

Feature curves = label boundaries, as in tet_prep_v2.feature_graph: the two
basal rims (inner/ring, outer/ring). Both are closed loops, so there are no
feature vertices -- the same situation as the closed O-ring curves on tistos.

Self-intersections (one test case has one at the basal rim) are repaired by a local
umbrella step on the vertices involved before meshing; the move is logged and
stored in transform.json.

Scale: the STLs are in metres (bounding-box diagonal ~0.15). AlgoHex and
the downstream tools have been tuned on tistos (diagonal ~4.3), so the
geometry is centred and scaled (default x25) and the transform is written to
<out>/transform.json; divide by `scale` and add `center` to go back.

    python charite_prep.py <case_dir> --out-dir <out_dir>
"""

import json
from pathlib import Path

import meshio
import numpy as np

import tet_prep as tp
import tet_prep_v2 as v2

PARTS = (("LV_inner", 1), ("LV_outer", 2), ("LV_connection", 3))
NAMES = {1: "endocardium", 2: "epicardium", 3: "basal_ring"}


def _merge(P, T, tol):
    key = np.round(P / tol).astype(np.int64)
    _, first, inv = np.unique(key, axis=0, return_index=True, return_inverse=True)
    return P[first], inv.reshape(-1)[T]


def load_case(case_dir, tol=1e-9):
    """LV_all.stl is the surface; every triangle gets the label of the part STL
    that contains the same triangle (matched on its three vertex positions).
    Unmatched triangles fall back to the part with the nearest triangle
    centroid; the number of fallbacks is reported."""
    from scipy.spatial import cKDTree
    m = meshio.read(Path(case_dir) / "LV_all.stl")
    nodes, T = _merge(m.points.astype(np.float64), m.cells_dict["triangle"], tol)
    keep = (T[:, 0] != T[:, 1]) & (T[:, 1] != T[:, 2]) & (T[:, 0] != T[:, 2])
    T = T[keep]

    def tri_keys(P, tris):
        k = np.round(P[tris] / tol).astype(np.int64)          # [n,3,3]
        k = np.sort(k.view([("", k.dtype)] * 3).reshape(len(tris), 3), axis=1)
        return [tuple(map(tuple, row.tolist())) for row in k]

    lab_of, cents, cent_lab = {}, [], []
    for name, sid in PARTS:
        p = meshio.read(Path(case_dir) / f"{name}.stl")
        Pp, Tp = p.points.astype(np.float64), p.cells_dict["triangle"]
        for k in tri_keys(Pp, Tp):
            lab_of[k] = sid
        cents.append(Pp[Tp].mean(1)); cent_lab.append(np.full(len(Tp), sid))
    ids = np.array([lab_of.get(k, 0) for k in tri_keys(nodes, T)])
    miss = ids == 0
    if miss.any():
        tree = cKDTree(np.vstack(cents))
        _, j = tree.query(nodes[T[miss]].mean(1))
        ids[miss] = np.concatenate(cent_lab)[j]
    print(f"[charite] LV_all: {len(T)} triangles, labels from parts: "
          f"{int((~miss).sum())} exact, {int(miss.sum())} by nearest centroid")
    return nodes, T, ids


def _seg_hits_tri(p, q, a, b, c, eps=1e-12):
    """Moller-Trumbore, open segment p-q against triangle abc."""
    d, e1, e2 = q - p, b - a, c - a
    h = np.cross(d, e2); det = e1 @ h
    if abs(det) < eps:
        return False
    f = 1.0 / det; s = p - a; u = f * (s @ h)
    if u < 0 or u > 1:
        return False
    qq = np.cross(s, e1); v = f * (d @ qq)
    if v < 0 or u + v > 1:
        return False
    t = f * (e2 @ qq)
    return 1e-9 < t < 1 - 1e-9


def self_intersections(nodes, tris):
    """Pairs of non-adjacent triangles that intersect (edge through triangle)."""
    from scipy.spatial import cKDTree
    C = nodes[tris].mean(1)
    r = np.linalg.norm(nodes[tris] - C[:, None], axis=2).max(1)
    hits = []
    for i, j in cKDTree(C).query_pairs(2 * r.max()):
        if set(tris[i]) & set(tris[j]) or np.linalg.norm(C[i] - C[j]) > r[i] + r[j]:
            continue
        A, B = nodes[tris[i]], nodes[tris[j]]
        if any(_seg_hits_tri(A[k], A[(k + 1) % 3], *B) for k in range(3)) or \
                any(_seg_hits_tri(B[k], B[(k + 1) % 3], *A) for k in range(3)):
            hits.append((i, j))
    return hits


def repair_self_intersections(nodes, tris, ids, max_iter=10):
    """Move the vertices of intersecting triangles to the mean of their edge
    neighbours (umbrella step) until the surface is free of self-intersections.
    One test case has one epicardium/basal-ring pair at the rim, which TetGen
    rejects ("A segment and a facet intersect"): a sliver ring triangle whose
    three corners all sit on the outer rim, pierced at a rim zigzag. Vertices
    on a label boundary (the feature curves) are therefore smoothed ALONG the
    curve only (mean of their two curve neighbours); the free umbrella step put
    86 deg kinks into the rim. Returns the max displacement."""
    from collections import defaultdict
    nb = defaultdict(set)
    vlab = defaultdict(set)
    for t, sid in zip(tris, ids):
        for k in range(3):
            nb[t[k]].update((t[(k + 1) % 3], t[(k + 2) % 3]))
            vlab[t[k]].add(int(sid))
    fixed = {v for v, ls in vlab.items() if len(ls) > 1}
    e2l = defaultdict(set)
    for t, sid in zip(tris, ids):
        for k in range(3):
            e2l[frozenset((t[k], t[(k + 1) % 3]))].add(int(sid))
    curve_nb = defaultdict(set)
    for e, ls in e2l.items():
        if len(ls) > 1:
            a, b = tuple(e)
            curve_nb[a].add(b); curve_nb[b].add(a)
    orig = nodes.copy()
    for it in range(max_iter):
        hits = self_intersections(nodes, tris)
        if not hits:
            moved = float(np.linalg.norm(nodes - orig, axis=1).max())
            if it:
                print(f"[charite] self-intersections repaired in {it} umbrella steps, "
                      f"max vertex move {moved:.2e}")
            return moved
        vs = {v for i, j in hits for v in (*tris[i], *tris[j])}
        for v in vs:
            ring = curve_nb[v] if v in fixed else nb[v]
            nodes[v] = np.mean([nodes[w] for w in ring], axis=0)
    raise RuntimeError(f"{len(hits)} self-intersections left after {max_iter} steps")


def load_surface_vtu(path):
    """A closed, already labelled triangle surface (cell data `label` 1/2/3),
    e.g. the myocardium reconstructed from surface patches."""
    m = meshio.read(path)
    T = m.cells_dict["triangle"]
    ids = np.asarray(m.cell_data_dict["label"]["triangle"]).astype(int)
    return np.asarray(m.points, dtype=np.float64), T, ids


def main(case_dir, out_dir, scale=25.0, size_max=0.12, surface_vtu=None):
    nodes, tris, ids = load_surface_vtu(surface_vtu) if surface_vtu else load_case(case_dir)
    hits = self_intersections(nodes, tris)
    print(f"[charite] self-intersecting triangle pairs: {len(hits)}")
    repair_move = repair_self_intersections(nodes, tris, ids) if hits else 0.0
    hist = tp.check_manifold(tris.tolist())
    print(f"[charite] {Path(surface_vtu or case_dir).name}: {len(nodes)} nodes, {len(tris)} triangles, "
          f"edge multiplicity {dict(hist)}, labels "
          f"{ {NAMES[k]: int((ids == k).sum()) for k in NAMES} }")
    if set(hist) - {2}:
        raise RuntimeError("merged surface is not a closed 2-manifold")
    center = 0.5 * (nodes.min(0) + nodes.max(0))
    nodes = (nodes - center) * scale

    feat, fverts = v2.feature_graph(nodes, tris.tolist(), ids)
    P, tets, idx, used = v2.mesh_interior(nodes, tris.tolist(), size_max)
    remap = {ot: idx[nt] for nt, ot in zip(range(1, len(used) + 1), used) if nt in idx}
    btris = [[remap[t] for t in tri] for tri in tris.tolist()]
    fe = [(remap[a], remap[b]) for a, b in feat]
    fv = [remap[v] for v in fverts]
    got = tp.tet_boundary_faces(tets, P)
    ok = {frozenset(t) for t in btris} == {frozenset(t) for t in got}
    print(f"[charite] boundary preserved: {'OK' if ok else 'MISMATCH'}")

    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    tp.write_algohex_vtk_analytic(P, tets, btris, ids, fe, fv, out / "tet.vtk")
    json.dump({"case": str(surface_vtu or case_dir), "scale": scale, "center": center.tolist(),
               "size_max": size_max, "labels": NAMES,
               "self_intersections_repaired": len(hits), "repair_max_move": repair_move,
               "n_tets": len(tets), "n_feature_edges": len(fe)},
              open(out / "transform.json", "w"), indent=1)
    print(f"[charite] wrote {out / 'tet.vtk'}: {len(P)} points, {len(tets)} tets")
    return out / "tet.vtk"


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("case_dir", nargs="?", default=None)
    ap.add_argument("--surface-vtu", default=None,
                    help="labelled closed surface instead of the STL case dir")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--scale", type=float, default=25.0)
    ap.add_argument("--size-max", type=float, default=0.12)
    a = ap.parse_args()
    if not (a.case_dir or a.surface_vtu):
        ap.error("give a case dir or --surface-vtu")
    main(a.case_dir, a.out_dir, a.scale, a.size_max, a.surface_vtu)
