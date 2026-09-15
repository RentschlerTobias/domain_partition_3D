"""Stage 1 for the blade alone: cap the open RUBLADE patch and mesh it fresh.

The complement of cropping tets out of the existing volume mesh
(`tet_prep_solid.apply_crop`): that approach carves through an existing
tessellation, and the cut inevitably zigzags across tet edges that were
never aligned to it -- confirmed visually (a "chaotic" cut face) and in
AlgoHex's own crashes on the result. RUBLADE itself needs none of that: it
is already a clean, given 2D surface, and its open boundary is a SINGLE loop
sitting exactly on the hub cylinder (measured on runner_best: 182 edges, all
at r=0.5000, no branch-cut wraparound in theta). Capping that loop and
tetrahedralizing fresh gives a mesh with no leftover fragments at all.

The cap is built by unrolling the loop onto the (theta*r, z) plane -- exact,
since a cylinder is developable -- meshing that flat polygon in 2D, and
mapping the result back onto r=0.5. RUBLADE's own triangles are kept
untouched; only the cap is new. The two share their boundary nodes exactly
(matched by position, not by gmsh's own renumbering -- the same lesson
`tet_prep_v5.remesh_patch` already has on record).

    PY blade_volume.py --msh runner_best.msh --out runner_blade_tet.vtk
"""

import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "experimentell" / "hex3d_algohex"))

from dp3d.extraction import parse_msh                                 # noqa: E402
import tet_prep as tp                                                 # noqa: E402
import tet_prep_v2 as v2                                              # noqa: E402

CAP_SID = -1


def read_blade_patch(msh, phys_name="RUBLADE"):
    """(points dict, triangle list, geom-id per triangle, physical id used).

    Points stay keyed by their ORIGINAL msh node id -- everything downstream
    matches by that id or by position, never by a freshly assigned index."""
    nodes, elements = parse_msh(msh)
    names = {}
    for line in Path(msh).read_text().splitlines():
        if line.startswith("$EndPhysicalNames"):
            break
    for line in Path(msh).read_text().splitlines():
        if '"' in line and line.strip() and line.strip()[0].isdigit():
            parts = line.split()
            if len(parts) == 3:
                names[int(parts[1])] = parts[2].strip('"')
    phys_id = {v: k for k, v in names.items()}[phys_name]

    tris, gids = [], []
    for et, phys, geom, nds in elements:
        if et in (2, 9) and phys == phys_id:
            tris.append([int(x) for x in nds[:3]])
            gids.append(int(geom))
    used = sorted({v for t in tris for v in t})
    print(f"[blade] {phys_name}: {len(tris)} triangles, {len(used)} points, "
          f"{len(set(gids))} geometrical sub-faces")
    return nodes, tris, gids, phys_id


def boundary_loop(tris):
    """The single ordered vertex loop bounding an open triangle patch.

    Raises if the patch has zero or more than one open boundary -- either
    means the seed selection is not what `blade_volume` assumes."""
    from collections import defaultdict
    edge_count = defaultdict(int)
    for t in tris:
        for i in range(3):
            a, b = t[i], t[(i + 1) % 3]
            e = (a, b) if a < b else (b, a)
            edge_count[e] += 1
    open_edges = [e for e, c in edge_count.items() if c == 1]
    adj = defaultdict(list)
    for a, b in open_edges:
        adj[a].append(b)
        adj[b].append(a)
    if not open_edges:
        raise RuntimeError("patch has no open boundary -- already closed?")
    s0 = open_edges[0][0]
    loop, seen, cur, prev = [s0], {s0}, s0, None
    while True:
        nxt = [w for w in adj[cur] if w != prev]
        if not nxt:
            raise RuntimeError("boundary walk dead-ended -- not a simple loop")
        nx = nxt[0]
        if nx == s0:
            break
        loop.append(nx)
        seen.add(nx)
        prev, cur = cur, nx
    if len(seen) != len(open_edges):
        raise RuntimeError(
            f"boundary has {len(open_edges)} edges but the walk only found "
            f"{len(seen)} vertices -- more than one loop, not handled")
    return loop


def cap_loop_on_cylinder(nodes, loop, r, h=None, verbose=True):
    """Triangulate the disk bounded by `loop`, constrained to sit on the
    cylinder of radius `r` the loop already lies on.

    Unrolling is exact for a cylinder (it is developable), so the flat 2D
    mesh gmsh produces maps back with no distortion. Returns (new_points
    dict keyed by fresh ids starting after `max(nodes)`, cap triangle list
    using original loop ids on the boundary and the fresh ids inside)."""
    import gmsh
    pts3 = np.array([nodes[v] for v in loop], float)
    theta = np.arctan2(pts3[:, 1], pts3[:, 0])
    if theta.max() - theta.min() > np.pi:
        raise RuntimeError("loop spans more than pi in theta -- the flat "
                            "unroll would wrap; not handled")
    x2 = theta * r
    y2 = pts3[:, 2]
    if h is None:
        d = np.hypot(np.diff(x2, append=x2[0]), np.diff(y2, append=y2[0]))
        h = float(np.median(d))

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("cap")
    pt_tags = [gmsh.model.geo.addPoint(x, y, 0, h) for x, y in zip(x2, y2)]
    n = len(pt_tags)
    line_tags = [gmsh.model.geo.addLine(pt_tags[i], pt_tags[(i + 1) % n])
                 for i in range(n)]
    cl = gmsh.model.geo.addCurveLoop(line_tags)
    sf = gmsh.model.geo.addPlaneSurface([cl])
    gmsh.model.geo.synchronize()
    # ExtendFromBoundary must be ON here -- unlike `remesh_patch`, which
    # turns it off to stop a FINE boundary flooding a surface that should
    # stay coarse, this cap has no resolution of its own yet and needs the
    # loop's spacing carried inward, or the interior falls back to a coarse
    # default (measured: 9 interior points for a 182-point boundary).
    for k, v in (("Mesh.MeshSizeExtendFromBoundary", 1),
                 ("Mesh.MeshSizeFromPoints", 1),
                 ("Mesh.MeshSizeFromCurvature", 0)):
        gmsh.option.setNumber(k, v)
    gmsh.model.mesh.generate(2)
    nt, nc, _ = gmsh.model.mesh.getNodes()
    Q = np.array(nc).reshape(-1, 3)[:, :2]
    et, _i, ev = gmsh.model.mesh.getElements(2)
    T = next(np.array(vv).reshape(-1, 3) for t, vv in zip(et, ev) if t == 2)
    row = {int(t): i for i, t in enumerate(nt)}
    T = np.array([[row[int(x)] for x in t] for t in T])
    gmsh.finalize()

    # map gmsh's output nodes back onto the loop by POSITION (gmsh assigns
    # its own tags), then fold theta*r back to 3D on the r-cylinder
    from scipy.spatial import cKDTree
    tree = cKDTree(np.column_stack([x2, y2]))
    d, j = tree.query(Q)
    back = {}
    nxt_id = max(nodes) + 1
    new_pts = {}
    for k in range(len(Q)):
        if d[k] < 1e-9:
            back[k] = loop[j[k]]
        else:
            th, z = Q[k, 0] / r, Q[k, 1]
            new_pts[nxt_id] = (r * np.cos(th), r * np.sin(th), z)
            back[k] = nxt_id
            nxt_id += 1
    cap_tris = [[back[int(a)] for a in t] for t in T]
    if verbose:
        print(f"[blade] cap: {len(loop)} boundary pts, {len(new_pts)} new "
              f"interior pts, {len(cap_tris)} triangles, h={h:.4f}")
    return new_pts, cap_tris


def tetrahedralize_shell(nodes, tris, verbose=True):
    """Fresh constrained-Delaunay tetrahedralization of the volume enclosed
    by a closed, watertight triangle shell, via TetGen.

    gmsh's documented STL-to-volume route (`classifySurfaces` +
    `createGeometry`, then `generate(3)`) turns the discrete mesh into
    smooth parametric patches first -- and on this shell that silently
    threw almost the whole 15558-triangle RUBLADE surface away (measured:
    455 tets back, an order of magnitude too few; some patches most likely
    failed automatic reparametrization and were dropped without error).
    TetGen's PLC mode is built for exactly this instead: keep the given
    boundary facets EXACTLY, fill the interior only -- no reparametrization
    step to lose triangles in.

    Returns (points, tets, idx) where idx maps an ORIGINAL node id (as used
    in `tris`) to its row in `points`; TetGen appends any new interior
    (Steiner) points after the input ones, so the mapping stays valid."""
    import tetgen
    used = sorted({v for t in tris for v in t})
    idx = {v: i for i, v in enumerate(used)}
    P_in = np.array([nodes[v] for v in used], float)
    F_in = np.array([[idx[v] for v in t] for t in tris], np.int64)
    tg = tetgen.TetGen(P_in, F_in)
    tg.tetrahedralize(order=1, nobisect=True, quality=False)
    P = np.asarray(tg.node, float)
    Tt = np.asarray(tg.elem, np.int64)
    if verbose:
        print(f"[blade] tetrahedralized: {len(P)} points, {len(Tt)} tets "
              f"({len(P) - len(used)} Steiner points added)")
    return P, Tt, idx


def main(msh, out_name, phys_name="RUBLADE", verbose=True):
    nodes, tris, gids, phys_id = read_blade_patch(msh, phys_name)
    loop = boundary_loop(tris)
    r = float(np.hypot(*np.array(nodes[loop[0]])[:2]))
    new_pts, cap_tris = cap_loop_on_cylinder(nodes, loop, r, verbose=verbose)
    nodes = dict(nodes)
    nodes.update(new_pts)

    all_tris = tris + cap_tris
    all_gids = gids + [CAP_SID] * len(cap_tris)
    P, tets, idx = tetrahedralize_shell(nodes, all_tris, verbose)

    # re-derive the boundary from the FRESH tet mesh and carry gids across
    # by node-set lookup against the shell we fed TetGen (`idx` maps an
    # original shell node id to its row in `P`, and TetGen with
    # nobisect=True does not move or split boundary facets)
    boundary = tp.tet_boundary_faces(tets, P)
    orig = {frozenset(idx[v] for v in t): g for t, g in zip(all_tris, all_gids)}
    sid_list = [orig[frozenset(t)] for t in boundary]

    tris_out = np.array(boundary, np.int64)
    sids_out = np.array(sid_list, int)
    hist = tp.check_manifold([list(t) for t in boundary])
    print(f"[blade] boundary edge multiplicity: {dict(hist)}")
    if set(hist) - {2}:
        raise RuntimeError("boundary is not a closed 2-manifold")

    nodes_c = {i: tuple(p) for i, p in enumerate(P)}
    feat, fverts = v2.feature_graph(nodes_c, [list(t) for t in boundary], sids_out)
    print(f"[blade] feature graph: {len(feat)} edges, {len(fverts)} vertices")

    out = Path(out_name)
    out.parent.mkdir(parents=True, exist_ok=True)
    tp.write_algohex_vtk_analytic(P, tets, [list(t) for t in boundary],
                                   sids_out, feat, fverts, out)
    print(f"[blade] wrote {out} ({len(P)} points, {len(tets)} tets, "
          f"{len(boundary)} boundary triangles)")
    return out


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--msh", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--phys", default="RUBLADE")
    a = ap.parse_args()
    main(a.msh, a.out, phys_name=a.phys)
