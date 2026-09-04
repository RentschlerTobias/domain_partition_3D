"""Stage 1 v3: reduced domain -- passage WITHOUT the blade O-grid shell.

Rationale (Tobias): do not make the frame field resolve the thin, twisted
blade boundary layer at all. Cut the existing structured O-grid around the
blade out of the domain, run AlgoHex only on the remaining volume, extract
its blocks, and later assemble them with the O-grid's own block structure.

Conformity is not an issue: TFI regenerates the interior of every block from
its boundary curves, so the two parts only have to agree at BLOCK level, and
the conforming-division MILP (dp3d/tmesh.py:743) handles that. The cell-level
mismatch between AlgoHex's quads and the O-grid's quads is irrelevant.

The removed region is the hex O-grid, geom regions 2-6 of the source MSH
(44 800 hexes). Its outer shell becomes a new interior-facing boundary of
the reduced domain, tagged as its own surface (SURF_SHELL); the blade itself
is no longer part of this domain.
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
OGRID_GEOM = {2, 3, 4, 5, 6}          # hex O-grid regions around the blade
SURF_SHELL = 8                        # new: O-grid outer shell
NAMES = dict(tp.SURF_NAMES)
NAMES[SURF_SHELL] = "ogrid_shell"
NAMES.pop(tp.SURF_BLADE, None)        # blade is not on this domain any more


def reduced_boundary(nodes, elements):
    """Boundary of the passage minus the O-grid region."""
    keep = [el for el in elements
            if not (el[0] in tp._VOL_FACES and el[2] in OGRID_GEOM)]
    ext = tp.volume_exterior_faces(keep)
    tris = tp.orient_and_triangulate(ext, nodes)
    hist = tp.check_manifold(tris)
    print(f"[v3] reduced domain: {len(tris)} boundary triangles, "
          f"edge multiplicity {dict(hist)}")
    if set(hist) - {2}:
        raise RuntimeError("reduced domain boundary is not a closed 2-manifold")
    return tris


def classify(nodes, tris):
    """hub / shroud / inlet / outlet / periodic analytically; everything left
    over is the O-grid shell (it replaces the blade on this domain)."""
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
    ids[(np.abs(r - tp.R_HUB) < 0.02) & (radial > 0.7)] = tp.SURF_HUB
    ids[(np.abs(r - tp.R_SHROUD) < 0.02) & (radial > 0.7)] = tp.SURF_SHROUD
    ids[(np.abs(cen[:, 2] - tp.Z_INLET) < 1e-6) & (axial > 0.7)] = tp.SURF_INLET
    ids[(np.abs(cen[:, 2] - tp.Z_OUTLET) < 1e-6) & (axial > 0.7)] = tp.SURF_OUTLET

    rest = ids == 0
    comp, ncomp = tp._connected_components(tris, rest)
    sizes = Counter(comp.values())
    print(f"[v3] remainder: {int(rest.sum())} tris in {ncomp} components "
          f"{sorted(sizes.values(), reverse=True)}")
    # the shell is the component that does NOT touch the outer walls; the two
    # periodic walls are the two whose theta is extreme
    th_mid, touches_outer = {}, {}
    for c in range(ncomp):
        sel = [i for i, cc in comp.items() if cc == c]
        th_mid[c] = np.median(np.arctan2(cen[sel, 1], cen[sel, 0]))
        rr = r[sel]
        touches_outer[c] = (np.abs(rr - tp.R_HUB).min() < 0.03 or
                            np.abs(rr - tp.R_SHROUD).min() < 0.03)
    outer = [c for c in range(ncomp) if touches_outer[c]]
    inner = [c for c in range(ncomp) if not touches_outer[c]]
    order = sorted(outer, key=lambda c: th_mid[c])
    per_a = order[0] if order else None
    per_b = order[-1] if len(order) > 1 else None
    for i, c in comp.items():
        if c == per_a:
            ids[i] = tp.SURF_PER_A
        elif c == per_b:
            ids[i] = tp.SURF_PER_B
        else:
            ids[i] = SURF_SHELL
    for sid in sorted(set(ids.tolist())):
        print(f"[v3]   surf {sid} {NAMES.get(sid, '?'):12s}: "
              f"{int((ids == sid).sum()):6d} tris")
    if inner:
        sel = [i for i, cc in comp.items() if cc in inner]
        print(f"[v3]   shell extent: r [{r[sel].min():.3f},{r[sel].max():.3f}] "
              f"z [{cen[sel, 2].min():.3f},{cen[sel, 2].max():.3f}]")
    return ids


def main(size_max=0.12, out_name="T1_9_tet_v3.vtk"):
    nodes, elements = parse_msh(MSH)
    tris = reduced_boundary(nodes, elements)
    ids = classify(nodes, tris)
    feat, fverts = v2.feature_graph(nodes, tris, ids)

    P, tets, idx, used = v2.mesh_interior(nodes, tris, size_max)
    remap = {ot: idx[nt] for nt, ot in
             zip(range(1, len(used) + 1), used) if nt in idx}
    btris = [[remap[t] for t in tri] for tri in tris]
    fe = [(remap[a], remap[b]) for a, b in feat]
    fv = [remap[v] for v in fverts]
    got = tp.tet_boundary_faces(tets, P)
    want = {frozenset(t) for t in btris}
    have = {frozenset(t) for t in got}
    print(f"[v3] boundary preserved: {len(want & have)}/{len(want)} "
          f"({'OK' if want == have else 'MISMATCH'})")

    out = MSH.parent / out_name
    tp.write_algohex_vtk_analytic(P, tets, btris, ids, fe, fv, out)
    return out


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--size-max", type=float, default=0.12)
    ap.add_argument("--out", default="T1_9_tet_v3.vtk")
    a = ap.parse_args()
    main(a.size_max, a.out)
