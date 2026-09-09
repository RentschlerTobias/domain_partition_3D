"""Stage 1 for a SOLID: the runner as a structural body, not a flow passage.

The complement of `tet_prep_v5`. That one cuts a flow passage out of a hybrid
mesh -- drop the blade O-grid, drop the prism layer, label seven surfaces,
hand the remainder to AlgoHex, glue the removed parts back afterwards. Here
nothing is removed and nothing is glued back: the body is the domain.

What carries over is the part that actually matters, and it is worth being
explicit about why. `feature_graph` puts a feature edge wherever two boundary
triangles carry DIFFERENT surface labels -- no angle criterion anywhere -- so
the labelling IS the constraint set for the frame field. The runner solid
arrives already labelled: RUHUB, RU_HUB_FIX and RUBLADE as physical groups, so
their mutual boundaries become the feature curves without any classification
step. Compare `tet_prep_v5.classify`, which needs 60 lines and an exact
lookup against the source MSH's tagged 2D elements to recover the same thing.

The mesh is second order -- 10-node tets, 6-node triangles, which is what a
structural solver wants -- and AlgoHex needs a linear one, so the mid-side
nodes are dropped and the point list compacted. That is a projection, not an
approximation: the corner nodes are unchanged.

    PY tet_prep_solid.py --msh runner_best.msh --out runner_tet.vtk
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

# gmsh element type -> number of CORNER nodes (the rest are mid-side)
CORNERS = {2: 3, 9: 3, 4: 4, 11: 4}
SURFACE_TYPES = (2, 9)
VOLUME_TYPES = (4, 11)


def read_solid(msh, verbose=True):
    """(points, tets, boundary triangles, surface id per triangle, names).

    Surfaces are taken from the physical tags as they are: no geometric
    classification, no thresholds. Whatever the mesher called a group is a
    group here."""
    nodes, elements = parse_msh(msh)

    tris, tri_phys, tets = [], [], []
    for et, phys, _geom, nds in elements:
        if et in SURFACE_TYPES:
            tris.append([int(x) for x in nds[:CORNERS[et]]])
            tri_phys.append(int(phys))
        elif et in VOLUME_TYPES:
            tets.append([int(x) for x in nds[:CORNERS[et]]])
    if not tets:
        raise RuntimeError(f"no volume elements in {msh}")

    used = sorted({v for t in tets for v in t} | {v for t in tris for v in t})
    idx = {v: i for i, v in enumerate(used)}
    P = np.array([nodes[v] for v in used], float)
    tets = np.array([[idx[v] for v in t] for t in tets], np.int64)
    tris = np.array([[idx[v] for v in t] for t in tris], np.int64)
    sids = np.array(tri_phys, int)
    if verbose:
        print(f"[solid] {len(P)} points, {len(tets)} tets, {len(tris)} "
              f"boundary triangles")
        for s, n in sorted(Counter(sids.tolist()).items()):
            print(f"[solid]   physical {s}: {n} triangles")
    return P, tets, tris, sids


def physical_names(msh):
    """{tag: name} from the $PhysicalNames block.

    Read for reporting only. The pipeline keys on the tags, so a mesh without
    names still works -- it just prints numbers."""
    out, on = {}, False
    for line in Path(msh).read_text().splitlines():
        if line.startswith("$PhysicalNames"):
            on = True
            continue
        if line.startswith("$EndPhysicalNames"):
            break
        if on and '"' in line:
            parts = line.split()
            out[int(parts[1])] = parts[2].strip('"')
    return out


def sharp_edges(P, tris, sids, deg=40.0, verbose=True):
    """Boundary edges that are a geometric crease, whatever the labels say.

    `tet_prep_v5` gets away without any angle criterion because on the flow
    passage the label boundaries ARE the geometric edges: inlet meets hub in a
    real crease, and the rings between the cut surfaces kink by 81 degrees.
    That coincidence does not hold for a solid. RUHUB, RU_HUB_FIX and RUBLADE
    are material and boundary-condition groups, not a decomposition by
    geometry -- the blade's trailing edge is knife-sharp and sits in the
    MIDDLE of RUBLADE, so a label-boundary criterion cannot see it.

    Measured on runner_best: 526 boundary edges kink by more than 40 degrees,
    and only 300 of them lie on a label boundary. The other 226 -- 43 % -- had
    no constraint on them, and AlgoHex's integrability optimisation diverged
    after 90 minutes with the Newton solver unable to resolve its KKT system.

    Returns the extra edges, i.e. creases that are NOT already label
    boundaries."""
    e2t = defaultdict(list)
    for i, t in enumerate(tris):
        for a in range(3):
            u, v = int(t[a]), int(t[(a + 1) % 3])
            e2t[(min(u, v), max(u, v))].append(i)
    nrm = {}

    def normal(i):
        if i not in nrm:
            a, b, c = P[np.asarray(tris[i])]
            n = np.cross(b - a, c - a)
            ln = np.linalg.norm(n)
            nrm[i] = n / ln if ln > 1e-15 else n
        return nrm[i]

    extra, angs = [], []
    for e, inc in e2t.items():
        if len(inc) != 2:
            continue
        d = float(np.clip(np.dot(normal(inc[0]), normal(inc[1])), -1, 1))
        ang = float(np.degrees(np.arccos(d)))
        angs.append(ang)
        if ang > deg and sids[inc[0]] == sids[inc[1]]:
            extra.append(e)
    if verbose:
        angs = np.array(angs)
        print(f"[solid] creases over {deg} deg: {int((angs > deg).sum())} of "
              f"{len(angs)} boundary edges, {len(extra)} of them inside a "
              f"single surface and therefore missed by the labels alone")
    return extra


def feature_vertices(feat):
    """Curve endpoints and branch points: valence != 2.

    Angle-based creases fray -- a crease can stop mid-surface where the
    geometry rounds off -- and an undeclared dangling end is what earned run
    v2 a SIGSEGV. Recomputing the set from the final edges declares them."""
    val = Counter()
    for a, b in feat:
        val[a] += 1
        val[b] += 1
    return sorted(v for v, c in val.items() if c != 2)


def check_closed(tris, verbose=True):
    """The boundary must be a closed 2-manifold, or AlgoHex has nothing to
    align to. Same check `tet_prep_v5` runs, and it has caught real damage
    before -- a re-meshing that left 447 dangling edges passed every other
    test."""
    hist = tp.check_manifold([list(t) for t in tris])
    if verbose:
        print(f"[solid] boundary edge multiplicity: {dict(hist)}")
    if set(hist) - {2}:
        raise RuntimeError("boundary is not a closed 2-manifold")


def main(msh, out_name, sharp=40.0, verbose=True):
    P, tets, tris, sids = read_solid(msh, verbose)
    names = physical_names(msh)
    if verbose:
        print(f"[solid] surfaces: "
              f"{ {names.get(s, s): int((sids == s).sum()) for s in sorted(set(sids.tolist()))} }")

    # orient consistently, then the feature graph from the label boundaries
    tri_list = [list(map(int, t)) for t in tris]
    check_closed(tri_list, verbose)
    nodes = {i: tuple(p) for i, p in enumerate(P)}
    feat, fverts = v2.feature_graph(nodes, tri_list, sids)
    if sharp:
        extra = sharp_edges(P, tri_list, sids, sharp, verbose)
        feat = list(feat) + extra
        fverts = feature_vertices(feat)
        if verbose:
            print(f"[solid] feature graph: {len(feat)} edges "
                  f"({len(feat) - len(extra)} from labels + {len(extra)} "
                  f"creases), {len(fverts)} feature vertices")

    out = Path(out_name)
    out.parent.mkdir(parents=True, exist_ok=True)
    tp.write_algohex_vtk_analytic(P, tets, tri_list, sids, feat, fverts, out)
    return out


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--msh", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--sharp", type=float, default=40.0, metavar="DEG",
                    help="also make a feature edge of every boundary crease "
                         "above DEG that the labels do not already catch; "
                         "0 disables (which is what tet_prep_v5 does, and "
                         "what diverged on this solid)")
    a = ap.parse_args()
    main(a.msh, a.out, sharp=a.sharp)
