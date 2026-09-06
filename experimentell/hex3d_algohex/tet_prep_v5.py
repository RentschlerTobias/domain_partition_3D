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


def main(size_max=0.12, out_name="T1_9_tet_v5.vtk", drop_prisms=True):
    nodes, elements = parse_msh(MSH)
    tris = reduced_boundary(nodes, elements, drop_prisms)
    look = tagged_2d_lookup(elements)
    removed_kind = removed_face_kind(elements)
    ids = classify(nodes, tris, look, removed_kind)
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
    ap.add_argument("--keep-prisms", action="store_true",
                    help="keep the hub/shroud boundary layer in the domain")
    a = ap.parse_args()
    main(a.size_max, a.out, drop_prisms=not a.keep_prisms)
