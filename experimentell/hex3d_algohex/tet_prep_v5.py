"""Stage 1 v5: reduced domain with PROPER surface labels.

The v3/v4 reduced domains (no blade O-grid, no prism boundary layer) ran the
frame-field stages beautifully -- 8 min to local meshability instead of 2 h
-- but their quantization blew up (OOM at the same point for -n 60000 and
-n 15000 alike, so it is not driven by the target cell count).

The likely cause fits a pattern seen three times now: fewer feature
constraints give a WORSE integer-grid map, not a better one.

    run   feature edges   surfaces   quantization
    v1    853 (427 flat)  21         completed
    v5    630             7          0 invalid tets
    v3    437 (clean)     7          IGM diverged
    v7/8  166             3          OOM

In the v3/v4 domains almost the whole boundary (16 502 of 18 548 triangles)
was one unconstrained "shell": once the boundary layers are removed nothing
touches hub or shroud any more, and the connectivity test that identified
the periodic walls no longer fires, so they were swallowed by the shell.

This version labels the reduced boundary EXACTLY, by matching each face
against the source MSH's own tagged 2D elements (geom 3 = inlet, 4 = outlet,
5/6 = the periodic pair) instead of guessing geometrically. Faces matching
nothing are the newly exposed interfaces to the removed cells, and are split
further into hub-side / shroud-side / blade-side shells so the field is
constrained there too.
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
import tet_prep_v2 as v2                                              # noqa: E402

MSH = REPO / "data" / "T1_9" / "T1_9_ru_gridGmsh.msh"
OGRID_GEOM = {2, 3, 4, 5, 6}          # volume geom regions of the blade O-grid
GEOM_INLET, GEOM_OUTLET = 3, 4        # tagged 2D element geom ids
GEOM_PER_A, GEOM_PER_B = 6, 5
GEOM_HUB, GEOM_SHROUD = 1, 2          # the real walls (r = 0.5 / 1.9)

S_INLET, S_OUTLET, S_PER_A, S_PER_B = 1, 2, 3, 4
S_BL_IFACE_HUB, S_BL_IFACE_SHROUD, S_OGRID_IFACE = 5, 6, 7
S_HUB, S_SHROUD = 8, 9                # only present with --keep-prisms
NAMES = {S_INLET: "inlet", S_OUTLET: "outlet", S_PER_A: "periodic_A",
         S_PER_B: "periodic_B",
         # These three are INTERFACES to cells cut out of the domain, not
         # walls. The name "shell_blade" cost three wrong diagnoses: the
         # blade wall is not in this mesh at all, the surface is the cut
         # face towards the removed O-grid.
         S_BL_IFACE_HUB: "bl_interface_hub",
         S_BL_IFACE_SHROUD: "bl_interface_shroud",
         S_OGRID_IFACE: "ogrid_interface",
         S_HUB: "hub", S_SHROUD: "shroud"}


def reduced_boundary(nodes, elements, drop_prisms=True):
    keep_types = {4, 7} if drop_prisms else {4, 6, 7}
    keep = [el for el in elements if el[0] in keep_types
            and el[2] not in OGRID_GEOM]
    ext = tp.volume_exterior_faces(keep)
    tris = tp.orient_and_triangulate(ext, nodes)
    hist = tp.check_manifold(tris)
    print(f"[v5] reduced domain: {len(tris)} boundary triangles, "
          f"edge multiplicity {dict(hist)}")
    if set(hist) - {2}:
        raise RuntimeError("boundary is not a closed 2-manifold")
    return tris


def tagged_2d_lookup(elements):
    """node-set -> geom id, for the MSH's own tagged 2D elements."""
    look = {}
    for et, _p, geom, nds in elements:
        if et == 2:
            look[frozenset(nds[:3])] = geom
        elif et == 3:
            nd = nds[:4]
            look[frozenset(nd)] = geom
            # BOTH diagonals. orient_and_triangulate reverses a quad whose
            # normal points inward and only then splits on 0-2, which for a
            # reversed quad is the OTHER diagonal of the original -- so half
            # the tagged quads silently failed to match. Visible as an
            # asymmetry: inlet matched r [0.614,1.787] while outlet matched
            # the full [0.499,1.899].
            look[frozenset((nd[0], nd[1], nd[2]))] = geom
            look[frozenset((nd[0], nd[2], nd[3]))] = geom
            look[frozenset((nd[0], nd[1], nd[3]))] = geom
            look[frozenset((nd[1], nd[2], nd[3]))] = geom
    return look


def removed_face_kind(elements):
    """face-key -> which kind of removed cell it belonged to.

    'ogrid' = the hexahedral O-grid around the blade (geom regions 2-6),
    'bl'    = the prismatic hub/shroud boundary layer."""
    kind = {}
    for et, _p, geom, nds in elements:
        if et not in tp._VOL_FACES:
            continue
        if geom in OGRID_GEOM:
            k = "ogrid"
        elif et == 6:
            k = "bl"
        else:
            continue
        for loc in tp._VOL_FACES[et]:
            face = [nds[i] for i in loc]
            kind[frozenset(face)] = k
            if len(face) == 4:
                # both diagonals, same reason as in tagged_2d_lookup
                kind[frozenset(face[:3])] = k
                kind[frozenset((face[0], face[2], face[3]))] = k
                kind[frozenset((face[0], face[1], face[3]))] = k
                kind[frozenset((face[1], face[2], face[3]))] = k
    return kind


def classify(nodes, tris, look, removed_kind):
    ids = np.zeros(len(tris), int)
    P = {t: np.asarray(nodes[t], float) for tri in tris for t in tri}
    cen = np.array([np.mean([P[t] for t in tri], axis=0) for tri in tris])
    r = np.hypot(cen[:, 0], cen[:, 1])

    n_match = 0
    for i, tri in enumerate(tris):
        g = look.get(frozenset(tri))
        if g == GEOM_INLET:
            ids[i] = S_INLET
        elif g == GEOM_OUTLET:
            ids[i] = S_OUTLET
        elif g == GEOM_PER_A:
            ids[i] = S_PER_A
        elif g == GEOM_PER_B:
            ids[i] = S_PER_B
        elif g == GEOM_HUB:
            ids[i] = S_HUB
        elif g == GEOM_SHROUD:
            ids[i] = S_SHROUD
        else:
            continue
        n_match += 1
    print(f"[v5] matched against tagged 2D elements: {n_match}/{len(tris)}")

    # Remainder = interfaces to the removed cells. Split them by WHICH kind of
    # cell was removed on the other side (prism boundary layer vs O-grid hex),
    # not by a radius threshold: a threshold cuts arbitrarily across the
    # triangulation and produces exactly the zigzag it is meant to avoid
    # (measured on a first attempt: p95 kink 174.8 deg, 318 of 679 nodes
    # above 30 deg). The cell-type boundary follows real mesh lines.
    rest = np.where(ids == 0)[0]
    for i in rest:
        ids[i] = S_OGRID_IFACE if removed_kind.get(frozenset(tris[i])) == "ogrid" \
            else S_BL_IFACE_HUB
    # the boundary-layer shell falls apart into hub side and shroud side
    bl = ids == S_BL_IFACE_HUB
    comp, ncomp = tp._connected_components(tris, bl)
    if ncomp >= 2:
        sizes = Counter(comp.values())
        rmed = {c: np.median(r[[i for i, cc in comp.items() if cc == c]])
                for c in range(ncomp)}
        shroud_c = max(rmed, key=lambda c: rmed[c])
        for i, c in comp.items():
            if c == shroud_c:
                ids[i] = S_BL_IFACE_SHROUD
        print(f"[v5] BL shell split into {ncomp} components "
              f"{sorted(sizes.values(), reverse=True)}")
    for sid in sorted(NAMES):
        n = int((ids == sid).sum())
        sel = ids == sid
        extra = ""
        if n:
            extra = (f"  r [{r[sel].min():.3f},{r[sel].max():.3f}]"
                     f" z [{cen[sel, 2].min():.2f},{cen[sel, 2].max():.2f}]")
        print(f"[v5]   surf {sid} {NAMES[sid]:14s}: {n:6d} tris{extra}")
    if (ids == 0).any():
        raise RuntimeError(f"{(ids == 0).sum()} triangles unclassified")
    return ids


def remesh_patch(nodes, tris, ids, sid, h):
    """Re-mesh one labelled surface isotropically, boundary node-for-node fixed.

    The O-grid interface is the one surface of the input whose triangulation
    is bad -- 24 % of its triangles below quality 0.2, because it is the
    O-grid's anisotropic quads cut in half (aspect median 4.6, up to 17).
    Every other surface is clean. It is a cut face, not geometry, so nothing
    stops us from meshing it better.

    Its boundary must survive exactly: those 224 edges are shared with the
    two boundary-layer interfaces, and they are the feature curves. Two gmsh
    routes fail here -- `classifySurfaces` invents its own boundary curves and
    re-meshes them (224 -> 93 edges), and a discrete surface without declared
    boundary meshes to nothing. The route that works is to add the boundary
    loops as discrete CURVES carrying their own mesh and hand them to
    `addDiscreteEntity(2, tag, boundary=...)`. `MeshSizeExtendFromBoundary`
    and `MeshSizeFromPoints` must be off, or the surface inherits the fine
    boundary spacing and comes out at 22k triangles instead of 6k."""
    import gmsh
    from collections import defaultdict
    keep = [t for t, i in zip(tris, ids) if i != sid]
    keep_ids = [i for i in ids if i != sid]
    patch = [t for t, i in zip(tris, ids) if i == sid]
    used = sorted({v for t in patch for v in t})
    tag = {v: i + 1 for i, v in enumerate(used)}
    ec = defaultdict(int)
    for t in patch:
        for k in range(3):
            u, v = t[k], t[(k + 1) % 3]
            ec[(min(u, v), max(u, v))] += 1
    bnd = [g for g, c in ec.items() if c == 1]
    adj = defaultdict(list)
    for a, b in bnd:
        adj[a].append(b)
        adj[b].append(a)
    loops, seen = [], set()
    for s0 in adj:
        if s0 in seen:
            continue
        lp, cur, prev = [s0], s0, None
        seen.add(s0)
        while True:
            nx = [w for w in adj[cur] if w != prev]
            if not nx or nx[0] == s0:
                break
            prev, cur = cur, nx[0]
            lp.append(cur)
            seen.add(cur)
        loops.append(lp)

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("patch")
    ct = []
    for li, lp in enumerate(loops, 1):
        c = gmsh.model.addDiscreteEntity(1, li)
        ct.append(c)
        gmsh.model.mesh.addNodes(1, c, [tag[v] for v in lp],
                                 np.array([nodes[v] for v in lp], float).ravel().tolist())
        seg = []
        for k in range(len(lp)):
            seg += [tag[lp[k]], tag[lp[(k + 1) % len(lp)]]]
        gmsh.model.mesh.addElementsByType(c, 1, [], seg)
    sf = gmsh.model.addDiscreteEntity(2, 1, ct)
    onb = {x for lp in loops for x in lp}
    inner = [v for v in used if v not in onb]
    gmsh.model.mesh.addNodes(2, sf, [tag[v] for v in inner],
                             np.array([nodes[v] for v in inner], float).ravel().tolist())
    gmsh.model.mesh.addElementsByType(sf, 2, [],
                                      [tag[int(v)] for t in patch for v in t])
    gmsh.model.mesh.createGeometry([(2, sf)])
    gmsh.model.geo.synchronize()
    for k, v in (("Mesh.MeshSizeMin", h), ("Mesh.MeshSizeMax", h),
                 ("Mesh.MeshSizeFromCurvature", 0),
                 ("Mesh.MeshSizeExtendFromBoundary", 0),
                 ("Mesh.MeshSizeFromPoints", 0), ("Mesh.Algorithm", 6)):
        gmsh.option.setNumber(k, v)
    gmsh.model.mesh.generate(2)
    nt, nc, _ = gmsh.model.mesh.getNodes()
    Q = np.array(nc).reshape(-1, 3)
    et, _i, ev = gmsh.model.mesh.getElements(2)
    T = next(np.array(vv).reshape(-1, 3) for t, vv in zip(et, ev) if t == 2)
    row = {int(t): i for i, t in enumerate(nt)}
    T = np.array([[row[int(x)] for x in t] for t in T])
    gmsh.finalize()

    # Map back by POSITION, not by tag. gmsh assigns its own node tags when
    # it regenerates the surface, and they collide with the ones handed in --
    # trusting them silently mismapped the boundary and left 447 dangling
    # edges plus one edge with three faces, which the manifold check caught.
    from scipy.spatial import cKDTree
    ref = np.array([nodes[v] for v in used], float)
    tree = cKDTree(ref)
    d, j = tree.query(Q)
    back = {}
    nxt = max(nodes) + 1
    for r in range(len(Q)):
        if d[r] < 1e-9:
            back[r] = used[j[r]]
        else:
            nodes[nxt] = tuple(Q[r])
            back[r] = nxt
            nxt += 1
    reused = sum(1 for r in back if back[r] in set(used))
    print(f"[v5]   {reused} of {len(Q)} output nodes coincide with input nodes")
    new = [[back[int(x)] for x in t] for t in T]
    print(f"[v5] remeshed surface {NAMES[sid]}: {len(patch)} -> {len(new)} "
          f"triangles at h={h}, {len(bnd)} boundary edges kept")
    return keep + new, np.array(keep_ids + [sid] * len(new), int)


def main(size_max=0.12, out_name="T1_9_tet_v5.vtk", drop_prisms=True,
         remesh_ogrid=None, merge_interfaces=False):
    nodes, elements = parse_msh(MSH)
    tris = reduced_boundary(nodes, elements, drop_prisms)
    look = tagged_2d_lookup(elements)
    removed_kind = removed_face_kind(elements)
    ids = classify(nodes, tris, look, removed_kind)
    if merge_interfaces:
        # Give the three artificial cut surfaces ONE label. Their mutual
        # boundaries -- two closed rings of 112 edges each around the blade
        # root and tip -- then stop being label boundaries and disappear from
        # the feature graph. They are artefacts of where we cut: in the
        # uncut domain the hub runs on and the blade runs on, there is no
        # edge there. The pinch sits on the hub ring, 0.89 away from the
        # nearest branch point, i.e. on its smooth part.
        # Whole closed curves vanish, so no dangling ends are created -- the
        # failure that gave run v2 a SIGSEGV.
        n0 = int(((ids == S_BL_IFACE_SHROUD) | (ids == S_OGRID_IFACE)).sum())
        ids[ids == S_BL_IFACE_SHROUD] = S_BL_IFACE_HUB
        ids[ids == S_OGRID_IFACE] = S_BL_IFACE_HUB
        print(f"[v5] merged the 3 cut surfaces into one: {n0} triangles "
              f"relabelled to {NAMES[S_BL_IFACE_HUB]}")
    if remesh_ogrid:
        tris, ids = remesh_patch(nodes, tris, ids, S_OGRID_IFACE, remesh_ogrid)
        hist = tp.check_manifold(tris)
        print(f"[v5] after remeshing, edge multiplicity {dict(hist)}")
        if set(hist) - {2}:
            raise RuntimeError("remeshing broke the closed 2-manifold")
    feat, fverts = v2.feature_graph(nodes, tris, ids)

    P, tets, idx, used = v2.mesh_interior(nodes, tris, size_max)
    remap = {ot: idx[nt] for nt, ot in
             zip(range(1, len(used) + 1), used) if nt in idx}
    btris = [[remap[t] for t in tri] for tri in tris]
    fe = [(remap[a], remap[b]) for a, b in feat]
    fv = [remap[v] for v in fverts]
    got = tp.tet_boundary_faces(tets, P)
    ok = {frozenset(t) for t in btris} == {frozenset(t) for t in got}
    print(f"[v5] boundary preserved: {'OK' if ok else 'MISMATCH'}")

    out = MSH.parent / out_name
    tp.write_algohex_vtk_analytic(P, tets, btris, ids, fe, fv, out)
    return out


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--size-max", type=float, default=0.12)
    ap.add_argument("--out", default="T1_9_tet_v5.vtk")
    ap.add_argument("--merge-interfaces", action="store_true",
                    help="one label for all three artificial cut surfaces, "
                         "which removes their mutual feature rings")
    ap.add_argument("--remesh-ogrid", type=float, default=None,
                    metavar="H",
                    help="re-mesh the O-grid interface isotropically at edge "
                         "length H, keeping its boundary node-for-node")
    ap.add_argument("--keep-prisms", action="store_true",
                    help="keep the hub/shroud boundary layer in the domain")
    a = ap.parse_args()
    main(a.size_max, a.out, drop_prisms=not a.keep_prisms,
         remesh_ogrid=a.remesh_ogrid, merge_interfaces=a.merge_interfaces)
