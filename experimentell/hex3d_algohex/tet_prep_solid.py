"""Stage 1 for a SOLID: the runner as a structural body, not a flow passage.

The complement of `tet_prep_v5`. That one cuts a flow passage out of a hybrid
mesh -- drop the blade O-grid, drop the prism layer, label seven surfaces,
hand the remainder to AlgoHex, glue the removed parts back afterwards. Here
nothing is removed and nothing is glued back: the body is the domain.

`feature_graph` puts a feature edge wherever two boundary triangles carry
DIFFERENT surface labels -- no angle criterion anywhere -- so the labelling
IS the constraint set for the frame field. The runner solid's PHYSICAL groups
(RUHUB, RU_HUB_FIX, RUBLADE) are too coarse for this: they are material/BC
groups, not a geometric decomposition, and the blade's trailing edge sits in
the MIDDLE of RUBLADE. A first version patched this with a global dihedral
threshold (`--sharp`), which produced a badly conditioned frame field --
AlgoHex diverged, then hit an apparent integer overflow after a 12-hour run
(see SESSION_2026-09-10.md).

The fix needs no new geometry and no re-export. Gmsh's own element tags carry
TWO ids per triangle, physical and geometrical, and the geometrical one is
never collapsed: RUBLADE alone spans 13 distinct geometrical faces (measured:
ids 7-19), because the CAD model already builds the blade from separate
patches (leading edge, trailing edge, fillets) and only the physical-group
step lumps them into one label. Checked directly on `runner_best.msh`: all
226 of the sharp creases (>40 deg) that the physical labels missed are
already a boundary between two geometrical ids. So the feature graph is built
from the geometrical id -- the same tag `tet_prep_v5.classify` reads for the
CFD passage -- while the physical id still becomes the exported surface
label, since RUHUB/RU_HUB_FIX/RUBLADE is what a structural BC wants to see.
`--sharp` remains available as an opt-in top-up, off by default.

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
    """(points, tets, boundary triangles, physical id, geometrical id).

    The physical id is the surface LABEL as exported (RUHUB/RU_HUB_FIX/
    RUBLADE) -- whatever the mesher called a group is a group. The
    geometrical id is the finer, un-collapsed CAD face id gmsh always
    carries alongside it; nothing here classifies or thresholds anything,
    both are read straight off the element tags."""
    nodes, elements = parse_msh(msh)

    tris, tri_phys, tri_geom, tets = [], [], [], []
    for et, phys, geom, nds in elements:
        if et in SURFACE_TYPES:
            tris.append([int(x) for x in nds[:CORNERS[et]]])
            tri_phys.append(int(phys))
            tri_geom.append(int(geom))
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
    gids = np.array(tri_geom, int)
    if verbose:
        print(f"[solid] {len(P)} points, {len(tets)} tets, {len(tris)} "
              f"boundary triangles")
        for s, n in sorted(Counter(sids.tolist()).items()):
            print(f"[solid]   physical {s}: {n} triangles")
        print(f"[solid]   {len(set(gids.tolist()))} geometrical faces "
              f"under those {len(set(sids.tolist()))} physical groups")
    return P, tets, tris, sids, gids


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


# The hub is a body of revolution over the WHOLE axial range, while the
# blade only occupies a middle band and grows outward from the hub wall --
# most of the tet volume is a plain revolve that needs no frame field.
# `--crop-hub` drops it so AlgoHex only has to solve the blade, not the
# whole runner; the dropped part is meant to be swept/revolved and glued
# back on separately (not done here -- this only produces AlgoHex's input).
#
# The cut is grown TOPOLOGICALLY, from the RUBLADE-labelled faces outward by
# `rings` layers of tet face-adjacency -- not by a coordinate threshold. A
# first attempt thresholded on radius and z and produced TWO independently
# shaped cuts (an inner shaft tube and two outer hub-wall caps) that met each
# other awkwardly near the fillet and crashed AlgoHex's singular-graph fixer
# (complex edges, invalid singular vertices). One seed, one BFS gives a
# single connected collar and a single connected cut surface instead -- the
# same lesson HANDOFF.md already has on record: a threshold cuts across the
# triangulation, growing from a real label does not.
CUT_SID = -1


def _tet_faces(t):
    a, b, c, d = (int(x) for x in t)
    return ((a, b, c), (a, b, d), (a, c, d), (b, c, d))


def blade_collar(tets, tris, sids, seed_sid, rings=3):
    """Mask selecting the tets within `rings` face-adjacency layers of a
    seed_sid-labelled boundary face -- the blade and a collar of the hub
    around it, grown outward one tet layer at a time."""
    face_to_tets = defaultdict(list)
    for ti, t in enumerate(tets):
        for f in _tet_faces(t):
            face_to_tets[frozenset(f)].append(ti)

    seed_faces = {frozenset(t) for t, s in zip(tris.tolist(), sids.tolist())
                  if s == seed_sid}
    frontier = {ti for f in seed_faces for ti in face_to_tets.get(f, [])}
    keep = set(frontier)
    for _ in range(rings):
        nxt = set()
        for ti in frontier:
            for f in _tet_faces(tets[ti]):
                nxt.update(face_to_tets.get(frozenset(f), []))
        nxt -= keep
        keep |= nxt
        frontier = nxt
    mask = np.zeros(len(tets), bool)
    mask[list(keep)] = True
    return mask


def heal_pinch(tets, keep, P, max_rounds=20, verbose=True):
    """Absorb tets around any non-manifold edge until the boundary closes.

    A growing collar can still touch ITSELF: two lobes meeting along a
    single edge without sharing a tet layer between them (the pinch measured
    on `--crop-rings 3`: 312 edges at multiplicity 4, at r~0.98 z~1.47 --
    the hub fillet). The fix is local, not a bigger uniform collar (that
    oscillated: 3, 3, 14 bad edges at rings 15/20/30, never zero): pull in
    every tet from the FULL mesh that touches a bad edge, which thickens
    exactly the join that is too thin, and recheck."""
    keep = keep.copy()
    for it in range(max_rounds):
        tris_k = tp.tet_boundary_faces(tets[keep], P)
        hist = tp.check_manifold(tris_k)
        if set(hist) <= {2}:
            if verbose:
                print(f"[solid] pinch heal: closed after {it} repair round(s)")
            return keep
        edge_count = Counter()
        for t in tris_k:
            for i in range(3):
                a, b = t[i], t[(i + 1) % 3]
                edge_count[(a, b) if a < b else (b, a)] += 1
        bad_edges = {e for e, c in edge_count.items() if c != 2}
        absorb = {ti for ti, t in enumerate(tets)
                  if any(a in t and b in t for a, b in bad_edges)}
        new = absorb - set(np.where(keep)[0])
        if verbose:
            print(f"[solid] pinch heal round {it + 1}: {len(bad_edges)} bad "
                  f"edges, absorbing {len(new)} more tets")
        if not new:
            raise RuntimeError("pinch heal stalled: bad edges remain but no "
                                "further tets touch them")
        keep[list(new)] = True
    raise RuntimeError(f"pinch heal did not converge in {max_rounds} rounds")


def apply_crop(tets, tris, sids, gids, seed_sid, rings=3, P=None, verbose=True):
    """Drop the hub/shaft tets and rebuild the boundary of what's left.

    The new boundary is a mix of ORIGINAL exterior faces (matched back by
    node set, so they keep their real sid/gid) and newly exposed CUT faces
    where a neighbour tet was removed -- those get `CUT_SID`, the same
    unconstrained-shell treatment `tet_prep_v5` gives its own cut
    interfaces. A cut face bordering a real surface still becomes a feature
    edge (the ids differ), which is exactly the intended constraint: the cut
    is where the swept-back hub piece will be welded on later."""
    orig = {frozenset(t): (int(s), int(g))
            for t, s, g in zip(tris.tolist(), sids.tolist(), gids.tolist())}
    keep = blade_collar(tets, tris, sids, seed_sid, rings)
    keep = heal_pinch(tets, keep, P, verbose=verbose)
    n0 = len(tets)
    tets = tets[keep]
    new_tris = tp.tet_boundary_faces(tets, P)
    new_sids, new_gids, n_cut = [], [], 0
    for t in new_tris:
        s, g = orig.get(frozenset(t), (CUT_SID, CUT_SID))
        if s == CUT_SID:
            n_cut += 1
        new_sids.append(s)
        new_gids.append(g)
    # recompact, since cropping orphans most of the original point set
    used = sorted({v for t in tets for v in t} | {v for t in new_tris for v in t})
    ridx = {v: i for i, v in enumerate(used)}
    P = P[used]
    tets = np.array([[ridx[int(v)] for v in t] for t in tets], np.int64)
    tris = np.array([[ridx[int(v)] for v in t] for t in new_tris], np.int64)
    sids = np.array(new_sids, int)
    gids = np.array(new_gids, int)
    if verbose:
        print(f"[solid] crop-hub: kept {len(tets)} of {n0} tets "
              f"({100 * len(tets) / n0:.1f}%), {len(tris)} boundary "
              f"triangles, {n_cut} of them newly exposed cut faces")
    return P, tets, tris, sids, gids


def main(msh, out_name, sharp=0.0, crop=False, crop_rings=3,
         crop_seed="RUBLADE", verbose=True):
    P, tets, tris, sids, gids = read_solid(msh, verbose)
    names = physical_names(msh)
    if crop:
        seed_sid = {v: k for k, v in names.items()}[crop_seed]
        P, tets, tris, sids, gids = apply_crop(
            tets, tris, sids, gids, seed_sid, crop_rings, P, verbose)
    if verbose:
        print(f"[solid] surfaces: "
              f"{ {names.get(s, s): int((sids == s).sum()) for s in sorted(set(sids.tolist()))} }")

    # orient consistently, then the feature graph from the GEOMETRICAL label
    # boundaries -- finer than the physical groups and, measured on
    # runner_best, a strict superset of the angle-threshold creases `--sharp`
    # was invented to catch (see module docstring).
    tri_list = [list(map(int, t)) for t in tris]
    check_closed(tri_list, verbose)
    nodes = {i: tuple(p) for i, p in enumerate(P)}
    feat, fverts = v2.feature_graph(nodes, tri_list, gids)
    if sharp:
        extra = sharp_edges(P, tri_list, sids, sharp, verbose)
        feat = list(set(feat) | {tuple(sorted(e)) for e in extra})
        fverts = feature_vertices(feat)
        if verbose:
            print(f"[solid] feature graph after --sharp top-up: {len(feat)} "
                  f"edges, {len(fverts)} feature vertices")

    # the exported surface label is still the coarse physical group -- that
    # is what a structural BC or a human wants to see, not the CAD face id
    out = Path(out_name)
    out.parent.mkdir(parents=True, exist_ok=True)
    tp.write_algohex_vtk_analytic(P, tets, tri_list, sids, feat, fverts, out)
    return out


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--msh", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--sharp", type=float, default=0.0, metavar="DEG",
                    help="also make a feature edge of every boundary crease "
                         "above DEG that the geometrical-id feature graph "
                         "does not already catch; off by default since it "
                         "measures empty on runner_best (see module "
                         "docstring) and a blanket angle threshold is what "
                         "produced the badly conditioned frame field before")
    ap.add_argument("--crop-hub", action="store_true",
                    help="drop the plain-revolution shaft/hub disc and hand "
                         "AlgoHex only a topological collar around the "
                         "blade; the dropped part still needs to be swept "
                         "and glued back on separately, which this flag "
                         "does not do")
    ap.add_argument("--crop-rings", type=int, default=3, metavar="N",
                    help="collar width in tet face-adjacency layers grown "
                         "from the blade surface (only with --crop-hub)")
    ap.add_argument("--crop-seed", default="RUBLADE", metavar="NAME",
                    help="physical group to grow the collar from")
    a = ap.parse_args()
    main(a.msh, a.out, sharp=a.sharp, crop=a.crop_hub,
         crop_rings=a.crop_rings, crop_seed=a.crop_seed)
