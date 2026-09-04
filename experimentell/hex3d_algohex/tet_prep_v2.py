"""Stage 1 v2: keep the ORIGINAL boundary mesh, mesh only the interior.

Motivated by inspection feedback on the v1 pipeline:
  * the gmsh `classifySurfaces` + `createGeometry` remesh throws away the
    original structured CFD surface mesh and replaces it with an isotropic
    Delaunay triangulation. Geometric fidelity is fine (max deviation 0.007,
    0.15% of the domain diagonal) but the structure and the alignment with
    hub/shroud/blade are lost.
  * feature curves derived from that remesh zigzag: the classification
    boundary cuts arbitrarily across an unstructured triangulation.
    Measured on the v1 curves: 18 of 393 curve nodes kink by more than 30
    deg, 7 by more than 60 deg (max 152 deg), and segment lengths spread
    over a factor of 146.

This version instead:
  1. takes the MSH's own tagged 2D elements and keeps only the geom ids that
     lie on the true outer skin (the interior O-grid block interfaces are
     tagged as 2D elements too, which is what makes the raw set
     non-manifold: edge multiplicity {2: 53664, 3: 724});
  2. groups those ids into the 7 physical surfaces;
  3. derives feature curves as boundaries between physical surfaces IN THE
     ORIGINAL MESH, so they follow the structured mesh lines and are smooth
     and complete by construction;
  4. meshes only the interior with gmsh, keeping the given surface
     triangulation fixed (discrete entities, no classifySurfaces).
"""

import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "experimentell" / "hex3d_algohex"))

from dp3d.extraction import parse_msh                                 # noqa: E402
import tet_prep as tp                                                 # noqa: E402

MSH = REPO / "data" / "T1_9" / "T1_9_ru_gridGmsh.msh"
OUT = REPO / "output" / "hex3d_algohex"

# geom ids confirmed to lie on the true outer skin (see export_vtk.py report)
OUTER_GEOM = [1, 2, 3, 4, 5, 6, 7, 8, 9, 13, 14, 15, 18, 19, 20, 23, 24, 25,
              28, 29, 30]


def outer_boundary(nodes, elements):
    """Outer-skin 2D elements, triangulated, with their geom id."""
    outer = {frozenset(f) for f, _c in tp.volume_exterior_faces(elements)}
    tris, gid = [], []
    for et, _p, geom, nds in elements:
        if et == 2:
            nd = nds[:3]
            if frozenset(nd) in outer:
                tris.append(tuple(nd))
                gid.append(geom)
        elif et == 3:
            # volume_exterior_faces returns quad faces as 4-node sets, so the
            # test must use the whole quad; splitting first and testing the
            # two triangles never matches and silently drops every quad patch
            # -- which is where the blade lives.
            nd = nds[:4]
            if frozenset(nd) in outer:
                tris.append((nd[0], nd[1], nd[2]))
                gid.append(geom)
                tris.append((nd[0], nd[2], nd[3]))
                gid.append(geom)
    return tris, np.asarray(gid)


def to_physical(nodes, tris, gid):
    """Map geom ids to the 7 physical surfaces, by geometry of each patch."""
    P = {t: np.asarray(nodes[t], float) for tri in tris for t in tri}
    cen = np.array([np.mean([P[t] for t in tri], axis=0) for tri in tris])
    nrm = np.zeros((len(tris), 3))
    for i, tri in enumerate(tris):
        n = np.cross(P[tri[1]] - P[tri[0]], P[tri[2]] - P[tri[0]])
        ln = np.linalg.norm(n)
        nrm[i] = n / ln if ln > 0 else n
    r = np.hypot(cen[:, 0], cen[:, 1])
    rad = np.zeros_like(nrm)
    rad[:, 0] = cen[:, 0] / np.maximum(r, 1e-12)
    rad[:, 1] = cen[:, 1] / np.maximum(r, 1e-12)
    radial = np.abs(np.einsum("ij,ij->i", nrm, rad))
    axial = np.abs(nrm[:, 2])

    ids = np.zeros(len(tris), int)
    # decide per GEOM PATCH by majority, so a patch is never split
    for g in np.unique(gid):
        m = gid == g
        rr, zz = r[m], cen[m, 2]
        if np.median(np.abs(rr - tp.R_HUB)) < 0.02 and np.median(radial[m]) > 0.7:
            ids[m] = tp.SURF_HUB
        elif np.median(np.abs(rr - tp.R_SHROUD)) < 0.02 and np.median(radial[m]) > 0.7:
            ids[m] = tp.SURF_SHROUD
        elif np.median(np.abs(zz - tp.Z_INLET)) < 1e-6 and np.median(axial[m]) > 0.7:
            ids[m] = tp.SURF_INLET
        elif np.median(np.abs(zz - tp.Z_OUTLET)) < 1e-6 and np.median(axial[m]) > 0.7:
            ids[m] = tp.SURF_OUTLET
        else:
            ids[m] = -1                      # blade or periodic, split below

    rest = ids == -1
    comp, ncomp = tp._connected_components(tris, rest)
    th_mid = {}
    for c in range(ncomp):
        sel = [i for i, cc in comp.items() if cc == c]
        th_mid[c] = np.median(np.arctan2(cen[sel, 1], cen[sel, 0]))
    order = sorted(th_mid, key=lambda c: th_mid[c])
    per_a, per_b = order[0], order[-1]
    for i, c in comp.items():
        ids[i] = (tp.SURF_PER_A if c == per_a else
                  tp.SURF_PER_B if c == per_b else tp.SURF_BLADE)
    print(f"[v2] remainder split into {ncomp} components "
          f"{sorted(Counter(comp.values()).values(), reverse=True)}")
    for sid in sorted(tp.SURF_NAMES):
        print(f"[v2]   surf {sid} {tp.SURF_NAMES[sid]:12s}: "
              f"{int((ids == sid).sum()):6d} tris")
    return ids


def feature_graph(nodes, tris, ids):
    """Feature curves = edges between different PHYSICAL surfaces. In the
    original structured mesh these follow the mesh lines, so they are smooth
    and closed by construction -- no threshold, no pruning needed."""
    e2t = defaultdict(list)
    for i, t in enumerate(tris):
        for a in range(3):
            p, q = t[a], t[(a + 1) % 3]
            e2t[(min(p, q), max(p, q))].append(i)
    feat = []
    for e, inc in e2t.items():
        if len(inc) == 2 and ids[inc[0]] != ids[inc[1]]:
            feat.append(e)
    val = Counter()
    for a, b in feat:
        val[a] += 1
        val[b] += 1
    fverts = sorted(v for v, c in val.items() if c != 2)
    hist = dict(sorted(Counter(val.values()).items()))
    print(f"[v2] feature edges: {len(feat)}, valence histogram {hist}, "
          f"feature vertices {len(fverts)}")

    P = {t: np.asarray(nodes[t], float) for tri in tris for t in tri}
    adj = defaultdict(list)
    for a, b in feat:
        adj[a].append(b)
        adj[b].append(a)
    ang = []
    for v, nb in adj.items():
        if len(nb) != 2:
            continue
        u, w = P[nb[0]] - P[v], P[nb[1]] - P[v]
        nu, nw = np.linalg.norm(u), np.linalg.norm(w)
        if nu > 1e-12 and nw > 1e-12:
            ang.append(180 - np.degrees(np.arccos(
                np.clip(np.dot(u, w) / (nu * nw), -1, 1))))
    ang = np.array(ang)
    if len(ang):
        print(f"[v2] curve smoothness: p95 {np.percentile(ang, 95):.1f} deg, "
              f"max {ang.max():.1f} deg, >30deg: {(ang > 30).sum()}/{len(ang)}")
    return feat, fverts


def mesh_interior(nodes, tris, size_max=0.12):
    """Tet-mesh the interior, keeping the given surface triangulation."""
    import gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("t19")

    used = sorted({t for tri in tris for t in tri})
    nmap = {t: i + 1 for i, t in enumerate(used)}
    coords = np.array([nodes[t] for t in used], float).ravel()

    s = gmsh.model.addDiscreteEntity(2)
    gmsh.model.mesh.addNodes(2, s, list(range(1, len(used) + 1)), coords)
    conn = [nmap[t] for tri in tris for t in tri]
    gmsh.model.mesh.addElementsByType(s, 2, [], conn)

    # A geo surface loop + volume over the discrete surface is what actually
    # meshes the interior while leaving the given triangulation untouched.
    # gmsh.model.addDiscreteEntity(3, -1, [s]) silently produces 0 tets, and
    # classifySurfaces/createGeometry (the v1 route) would remesh the surface.
    loop = gmsh.model.geo.addSurfaceLoop([s])
    gmsh.model.geo.addVolume([loop])
    gmsh.model.geo.synchronize()
    gmsh.option.setNumber("Mesh.MeshSizeMax", size_max)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.model.mesh.generate(3)

    ntag, ncoord, _ = gmsh.model.mesh.getNodes()
    P = np.asarray(ncoord, float).reshape(-1, 3)
    idx = {int(t): i for i, t in enumerate(ntag)}
    tets = []
    et, _tg, en = gmsh.model.mesh.getElements(3)
    for t, nd in zip(et, en):
        if t != 4:
            continue
        for row in np.asarray(nd, np.int64).reshape(-1, 4):
            tets.append([idx[int(x)] for x in row])
    gmsh.finalize()
    print(f"[v2] interior meshed: {len(P)} nodes, {len(tets)} tets "
          f"(boundary kept fixed)")
    return P, tets, idx, used


def main(size_max=0.12, out_name="T1_9_tet_v2.vtk"):
    nodes, elements = parse_msh(MSH)
    tris, gid = outer_boundary(nodes, elements)
    print(f"[v2] outer skin: {len(tris)} triangles over geom ids "
          f"{sorted(set(gid.tolist()))}")
    hist = tp.check_manifold(tris)
    print(f"[v2] edge multiplicity: {dict(hist)}")
    ids = to_physical(nodes, tris, gid)
    feat, fverts = feature_graph(nodes, tris, ids)

    P, tets, idx, used = mesh_interior(nodes, tris, size_max)
    # original node tag -> new index (gmsh keeps our node numbering 1..N for
    # the boundary nodes we fed in, and appends interior nodes after them)
    remap = {ot: idx[nt] for nt, ot in
             zip(range(1, len(used) + 1), used) if nt in idx}
    missing = [t for tri in tris for t in tri if t not in remap]
    if missing:
        raise RuntimeError(f"{len(set(missing))} boundary nodes lost by gmsh")
    btris = [[remap[t] for t in tri] for tri in tris]
    fe = [(remap[a], remap[b]) for a, b in feat]
    fv = [remap[v] for v in fverts]

    # gate: the tet mesh's own boundary must be exactly the surface we fed in
    got = tp.tet_boundary_faces(tets, P)
    want = {frozenset(t) for t in btris}
    have = {frozenset(t) for t in got}
    print(f"[v2] boundary preserved: {len(want & have)}/{len(want)} triangles "
          f"({'OK' if want == have else 'MISMATCH'})")

    out = MSH.parent / out_name
    tp.write_algohex_vtk_analytic(P, tets, btris, ids, fe, fv, out)
    return out


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--size-max", type=float, default=0.12)
    ap.add_argument("--out", default="T1_9_tet_v2.vtk")
    a = ap.parse_args()
    main(a.size_max, a.out)
