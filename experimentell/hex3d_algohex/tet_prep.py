"""Stage 1 (feat/algohex-3d-frame-field): T1_9 pure-tet volume + AlgoHex
feature-tag VTK, the input to AlgoHex's ``HexMeshing`` (see PLAN.md).

The source MSH (``data/T1_9/T1_9_ru_gridGmsh.msh``) is a HYBRID volume mesh
(tets/hexes/prisms/pyramids), not something AlgoHex can consume directly.
Its explicit tagged 2D "surface" elements turned out to be an unreliable
guide to the true boundary: the same physical surface (e.g. the hub) is
represented redundantly by both a coarse triangle patch AND separate small
quad O-grid patches at a different geom id -- merging them naively produces
a non-manifold soup (edges shared by 3 triangles). Instead this module
derives the boundary purely topologically: a face referenced by exactly one
volume cell is on the true exterior (generalizes
``dp3d.extraction.hex_exterior_quads`` from hex-only to the hybrid mix).
"""

import math
import sys
from collections import Counter
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))
from dp3d.extraction import parse_msh, write_stl                          # noqa: E402

MSH = REPO / "data" / "T1_9" / "T1_9_ru_gridGmsh.msh"
OUT = REPO / "output" / "hex3d_algohex"

# VTK legacy cell type codes used by AlgoHex's VtkColorReader.
_VTK_VERTEX, _VTK_LINE, _VTK_TRIANGLE, _VTK_TETRA = 1, 3, 5, 10

# T1_9 analytic geometry (verified against the source MSH, see PROGRESS.md).
R_HUB, R_SHROUD = 0.5, 1.9
Z_INLET, Z_OUTLET = 0.0, 2.5

# Physical surface ids used as AlgoHex face colors. Every boundary face must
# get a nonzero color (the field has to stay tangent to ALL walls, not just
# sharp ones); the ids additionally separate the physical patches so that
# only inter-patch boundaries become feature curves.
SURF_HUB, SURF_SHROUD, SURF_INLET, SURF_OUTLET = 1, 2, 3, 4
SURF_BLADE, SURF_PER_A, SURF_PER_B = 5, 6, 7
SURF_NAMES = {SURF_HUB: "hub", SURF_SHROUD: "shroud", SURF_INLET: "inlet",
              SURF_OUTLET: "outlet", SURF_BLADE: "blade",
              SURF_PER_A: "periodic_A", SURF_PER_B: "periodic_B"}

# A blade LE/TE is a genuine sharp crease INSIDE the blade patch, so sharp
# dihedral has to be honoured in addition to patch boundaries.
SHARP_DIHEDRAL_DEG = 40.0

# Gmsh low-order volume element faces (node-index SETS; winding is fixed up
# geometrically in orient_and_triangulate, so the order listed here doesn't
# matter). etype: 4=tet, 5=hex, 6=prism/wedge, 7=pyramid.
_VOL_FACES = {
    4: [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)],
    5: [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4),
        (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)],
    6: [(0, 1, 2), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)],
    7: [(0, 1, 2, 3), (0, 1, 4), (1, 2, 4), (2, 3, 4), (3, 0, 4)],
}


def volume_exterior_faces(elements):
    """True exterior boundary faces of the hybrid volume: faces referenced
    by exactly one volume cell. Returns [(face_node_tags, cell_node_tags)]
    -- the cell's nodes are kept alongside for the outward-orientation fix."""
    count = Counter()
    faces = {}
    cell_of = {}
    for etype, _phys, _geom, nds in elements:
        if etype not in _VOL_FACES:
            continue
        for loc in _VOL_FACES[etype]:
            face = tuple(nds[i] for i in loc)
            key = frozenset(face)
            count[key] += 1
            faces[key] = face
            cell_of[key] = nds
    return [(faces[k], cell_of[k]) for k, c in count.items() if c == 1]


def orient_and_triangulate(ext_faces, nodes):
    """Outward-orient each exterior face (normal points away from its cell's
    centroid) and split quads into 2 triangles with a fixed diagonal."""
    tris = []
    for face, cell_nds in ext_faces:
        pts = np.array([nodes[t] for t in face])
        cc = np.mean([nodes[t] for t in cell_nds], axis=0)
        fc = pts.mean(axis=0)
        n = np.cross(pts[1] - pts[0], pts[2] - pts[0])
        if np.dot(n, fc - cc) < 0:
            face = face[::-1]
        if len(face) == 3:
            tris.append(face)
        else:
            tris.append((face[0], face[1], face[2]))
            tris.append((face[0], face[2], face[3]))
    return tris


def check_manifold(tris):
    edge_count = Counter()
    for t in tris:
        for i in range(3):
            a, b = t[i], t[(i + 1) % 3]
            edge_count[(a, b) if a < b else (b, a)] += 1
    hist = Counter(edge_count.values())
    return hist


def extract_boundary_stl():
    """Stage 1a: topological exterior boundary of T1_9 -> STL, manifold-
    checked. Cheap (pure numpy), re-run standalone for inspection."""
    OUT.mkdir(parents=True, exist_ok=True)
    print(f"[tet_prep] parsing {MSH} ...")
    nodes, elements = parse_msh(MSH)
    print(f"[tet_prep] {len(nodes)} nodes, {len(elements)} elements")

    ext = volume_exterior_faces(elements)
    print(f"[tet_prep] {len(ext)} topological exterior faces")
    tris = orient_and_triangulate(ext, nodes)
    print(f"[tet_prep] {len(tris)} boundary triangles after quad split")

    hist = check_manifold(tris)
    print(f"[tet_prep] edge-multiplicity histogram: {dict(hist)}")
    if set(hist) - {2}:
        bad = sum(v for k, v in hist.items() if k != 2)
        print(f"[tet_prep] WARNING: {bad} non-manifold edge instances "
              f"(expected: only multiplicity 2)")

    stl_path = OUT / "T1_9_boundary_ext.stl"
    write_stl(tris, nodes, stl_path)
    return stl_path


# --------------------------------------------------------------------------
# Stage 1b: gmsh volume tetrahedralization + AlgoHex feature-colored VTK
# --------------------------------------------------------------------------

def build_tet_mesh(stl_path, elem_size_max=0.08, elem_size_min=0.02,
                    surface_angle=math.pi / 6, curve_angle=math.pi / 3):
    """gmsh: merge the boundary STL, recover surfaces/curves/points from
    dihedral angle (classifySurfaces/createGeometry -- same recipe as
    experimentell/gmsh_pipeline/remesh_step1.py, proven on this same T1_9
    geometry), mesh the enclosed volume with tets at a coarsened uniform
    size (the STL-native boundary resolution is far too fine for a first
    AlgoHex pass). Returns the gmsh module left initialized so the caller
    can pull element data before finalizing."""
    import gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.merge(str(stl_path))
    gmsh.model.mesh.classifySurfaces(surface_angle, True, True, curve_angle)
    gmsh.model.mesh.createGeometry()

    surfs = [t for d, t in gmsh.model.getEntities(2)]
    sl = gmsh.model.geo.addSurfaceLoop(surfs)
    gmsh.model.geo.addVolume([sl])
    gmsh.model.geo.synchronize()

    n_pts = len(gmsh.model.getEntities(0))
    n_curves = len(gmsh.model.getEntities(1))
    print(f"[tet_prep] classified: {n_pts} points, {n_curves} curves, "
          f"{len(surfs)} surfaces")

    gmsh.option.setNumber("Mesh.MeshSizeMax", elem_size_max)
    gmsh.option.setNumber("Mesh.MeshSizeMin", elem_size_min)
    gmsh.model.mesh.generate(3)

    n_nodes = len(gmsh.model.mesh.getNodes()[0])
    n_tets = sum(len(t) for t in gmsh.model.mesh.getElements(3)[1])
    print(f"[tet_prep] volume mesh: {n_nodes} nodes, {n_tets} tets")
    return gmsh


def tet_boundary_faces(tets, points):
    """Boundary triangles of a pure tet mesh: faces used by exactly one tet,
    wound so their normals point OUTWARD (away from the owning tet's apex).

    Consistent winding is essential: without it adjacent boundary triangles
    can come out anti-parallel, and the dihedral test used for sharp-crease
    detection then reports ~180 deg on perfectly flat surfaces (observed:
    154 bogus 'creases' on the planar inlet)."""
    count = Counter()
    rep = {}
    for t in tets:
        for loc in ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)):
            face = tuple(int(t[i]) for i in loc)
            key = frozenset(face)
            count[key] += 1
            # remember the face together with the tet's 4th (apex) vertex
            apex = [int(v) for v in t if int(v) not in key]
            rep[key] = (face, apex[0] if apex else None)
    out = []
    for k, c in count.items():
        if c != 1:
            continue
        face, apex = rep[k]
        n = np.cross(points[face[1]] - points[face[0]],
                     points[face[2]] - points[face[0]])
        if apex is not None and np.dot(n, points[face[0]] - points[apex]) < 0:
            face = (face[0], face[2], face[1])
        out.append(face)
    return out


def _connected_components(tris, seed_mask):
    """Edge-connected components among the triangles selected by seed_mask."""
    from collections import defaultdict
    idxs = [i for i, m in enumerate(seed_mask) if m]
    e2t = defaultdict(list)
    for i in idxs:
        t = tris[i]
        for a in range(3):
            p, q = t[a], t[(a + 1) % 3]
            e2t[(min(p, q), max(p, q))].append(i)
    comp = {}
    cid = 0
    for start in idxs:
        if start in comp:
            continue
        stack = [start]
        comp[start] = cid
        while stack:
            i = stack.pop()
            t = tris[i]
            for a in range(3):
                p, q = t[a], t[(a + 1) % 3]
                for j in e2t[(min(p, q), max(p, q))]:
                    if j not in comp:
                        comp[j] = cid
                        stack.append(j)
        cid += 1
    return comp, cid


def classify_boundary(points, btris):
    """Assign each boundary triangle its PHYSICAL surface id.

    Replaces gmsh `classifySurfaces` ids, which split the smooth hub/shroud
    cylinders into many parametrization patches; every such artificial patch
    boundary became a feature curve, over-constraining the octahedral field
    (measured: 427 of 853 tagged feature edges were geometrically flat) and
    destroying integrability. Hub/shroud/inlet/outlet are recognised
    analytically (position + normal direction); the remainder splits by
    connectivity into the blade and the two periodic side walls.
    """
    P = points
    ids = np.zeros(len(btris), np.int64)
    cen = np.array([P[list(t)].mean(axis=0) for t in btris])
    nrm = np.zeros((len(btris), 3))
    for i, t in enumerate(btris):
        n = np.cross(P[t[1]] - P[t[0]], P[t[2]] - P[t[0]])
        ln = np.linalg.norm(n)
        nrm[i] = n / ln if ln > 0 else n

    r = np.hypot(cen[:, 0], cen[:, 1])
    rad = np.zeros((len(btris), 3))
    rad[:, 0] = cen[:, 0] / np.maximum(r, 1e-12)
    rad[:, 1] = cen[:, 1] / np.maximum(r, 1e-12)
    radial = np.abs(np.einsum("ij,ij->i", nrm, rad))
    axial = np.abs(nrm[:, 2])

    ids[(np.abs(r - R_HUB) < 0.02) & (radial > 0.7)] = SURF_HUB
    ids[(np.abs(r - R_SHROUD) < 0.02) & (radial > 0.7)] = SURF_SHROUD
    ids[(np.abs(cen[:, 2] - Z_INLET) < 1e-6) & (axial > 0.7)] = SURF_INLET
    ids[(np.abs(cen[:, 2] - Z_OUTLET) < 1e-6) & (axial > 0.7)] = SURF_OUTLET

    # remainder = blade + the two periodic side walls, separated by connectivity
    rest = ids == 0
    comp, ncomp = _connected_components(btris, rest)
    sizes = Counter(comp.values())
    print(f"[tet_prep] unclassified remainder: {rest.sum()} tris in "
          f"{ncomp} components, sizes={sorted(sizes.values(), reverse=True)}")

    # the periodic pair maps onto each other under a pitch rotation about z;
    # the blade does not -- identify it as the component whose theta-extent
    # is widest (it wraps with the blade twist) is unreliable, so use the
    # rotation test directly against the known 90 deg pitch (4 blades).
    theta_mid = {}
    for c in range(ncomp):
        sel = [i for i, cc in comp.items() if cc == c]
        th = np.arctan2(cen[sel, 1], cen[sel, 0])
        theta_mid[c] = np.median(th)
    order = sorted(theta_mid, key=lambda c: theta_mid[c])
    # blade sits between the two periodic walls in theta
    per_a, blade_c, per_b = order[0], order[len(order) // 2], order[-1]
    if ncomp != 3:
        print(f"[tet_prep] WARNING: expected 3 remainder components "
              f"(blade + 2 periodic), got {ncomp}")
    for i, c in comp.items():
        ids[i] = (SURF_PER_A if c == per_a else
                  SURF_PER_B if c == per_b else SURF_BLADE)

    for sid in sorted(SURF_NAMES):
        n = int((ids == sid).sum())
        print(f"[tet_prep]   surf {sid} {SURF_NAMES[sid]:12s}: {n:6d} tris")
    if (ids == 0).any():
        print(f"[tet_prep] WARNING: {(ids == 0).sum()} triangles unclassified")
    return ids


def feature_edges_and_vertices(points, btris, tri_ids):
    """Feature curves = boundary between DIFFERENT physical surfaces, plus
    genuinely sharp creases inside a patch (blade LE/TE). Feature vertices =
    junctions where 3+ feature curves meet."""
    from collections import defaultdict
    P = points
    e2t = defaultdict(list)
    for i, t in enumerate(btris):
        for a in range(3):
            p, q = t[a], t[(a + 1) % 3]
            e2t[(min(p, q), max(p, q))].append(i)

    def tri_normal(t):
        n = np.cross(P[t[1]] - P[t[0]], P[t[2]] - P[t[0]])
        ln = np.linalg.norm(n)
        return n / ln if ln > 0 else n

    cos_thr = math.cos(math.radians(SHARP_DIHEDRAL_DEG))
    feat, n_patch, n_sharp = [], 0, 0
    for e, inc in e2t.items():
        if len(inc) != 2:
            continue
        i, j = inc
        if tri_ids[i] != tri_ids[j]:
            feat.append(e)
            n_patch += 1
        else:
            d = float(np.dot(tri_normal(btris[i]), tri_normal(btris[j])))
            if d < cos_thr:
                feat.append(e)
                n_sharp += 1
    print(f"[tet_prep] feature edges: {len(feat)} "
          f"({n_patch} patch boundaries + {n_sharp} sharp creases)")

    inc_cnt = Counter()
    for a, b in feat:
        inc_cnt[a] += 1
        inc_cnt[b] += 1
    fverts = sorted(v for v, c in inc_cnt.items() if c >= 3)
    print(f"[tet_prep] feature vertices (3+ incident feature edges): "
          f"{len(fverts)}")
    return feat, fverts


def write_algohex_vtk_analytic(points, tets, btris, tri_ids, feat_edges,
                               feat_verts, out_path):
    """AlgoHex VTK v2.0 ASCII: vertex/line/triangle/tet cells over one shared
    point list, with an int CELL_DATA color per cell (read back by AlgoHex as
    vertex_colors / edge_colors / face_colors; nonzero = feature)."""
    cells = []
    for v in feat_verts:
        cells.append((_VTK_VERTEX, [int(v)], 1))
    for a, b in feat_edges:
        cells.append((_VTK_LINE, [int(a), int(b)], 1))
    for t, sid in zip(btris, tri_ids):
        cells.append((_VTK_TRIANGLE, [int(x) for x in t], int(sid)))
    for t in tets:
        cells.append((_VTK_TETRA, [int(x) for x in t], 0))

    with open(out_path, "w") as f:
        f.write("# vtk DataFile Version 2.0\n")
        f.write("dp3d T1_9 AlgoHex tetmesh (analytic feature tags)\n")
        f.write("ASCII\nDATASET UNSTRUCTURED_GRID\n")
        f.write(f"POINTS {len(points)} float\n")
        for p in points:
            f.write(f"{p[0]} {p[1]} {p[2]}\n")
        total = sum(1 + len(c[1]) for c in cells)
        f.write(f"CELLS {len(cells)} {total}\n")
        for _, pts, _ in cells:
            f.write(f"{len(pts)} " + " ".join(str(p) for p in pts) + "\n")
        f.write(f"CELL_TYPES {len(cells)}\n")
        for vt, _, _ in cells:
            f.write(f"{vt}\n")
        f.write(f"CELL_DATA {len(cells)}\n")
        f.write("SCALARS color int 1\nLOOKUP_TABLE default\n")
        for _, _, col in cells:
            f.write(f"{col}\n")
    print(f"[tet_prep] wrote {out_path} ({len(points)} points, "
          f"{len(cells)} cells: {len(feat_verts)} vertex, "
          f"{len(feat_edges)} line, {len(btris)} triangle, {len(tets)} tetra)")


def write_algohex_vtk(gmsh, out_path):
    """Assemble AlgoHex's expected VTK v2.0 ASCII input: vertex/line/
    triangle/tet cells sharing one point list, colored by a single
    CELL_DATA int scalar (feature point/edge/face id; tets uncolored).
    Mirrors the hand-rolled legacy-VTK writer style of
    dp3d.extraction._write_quad_vtk."""
    tag, coords, _ = gmsh.model.mesh.getNodes()
    tag = np.asarray(tag, dtype=np.int64)
    coords = np.asarray(coords, dtype=float).reshape(-1, 3)
    idx = {int(t): i for i, t in enumerate(tag)}

    cells = []   # (vtk_type, [point-indices], color)

    for d, p in gmsh.model.getEntities(0):
        etypes, etags, enodes = gmsh.model.mesh.getElements(0, p)
        for nds in enodes:
            for i in range(0, len(nds), 1):
                cells.append((_VTK_VERTEX, [idx[int(nds[i])]], p))

    for d, c in gmsh.model.getEntities(1):
        etypes, etags, enodes = gmsh.model.mesh.getElements(1, c)
        for et, nds in zip(etypes, enodes):
            if et != 1:            # 2-node line
                continue
            nds = np.asarray(nds, dtype=np.int64).reshape(-1, 2)
            for a, b in nds:
                cells.append((_VTK_LINE, [idx[int(a)], idx[int(b)]], c))

    for d, s in gmsh.model.getEntities(2):
        etypes, etags, enodes = gmsh.model.mesh.getElements(2, s)
        for et, nds in zip(etypes, enodes):
            if et != 2:            # 3-node triangle
                continue
            nds = np.asarray(nds, dtype=np.int64).reshape(-1, 3)
            for a, b, c2 in nds:
                cells.append((_VTK_TRIANGLE,
                              [idx[int(a)], idx[int(b)], idx[int(c2)]], s))

    etypes, etags, enodes = gmsh.model.mesh.getElements(3)
    for et, nds in zip(etypes, enodes):
        if et != 4:                 # 4-node tet
            continue
        nds = np.asarray(nds, dtype=np.int64).reshape(-1, 4)
        for a, b, c2, d2 in nds:
            cells.append((_VTK_TETRA,
                          [idx[int(a)], idx[int(b)], idx[int(c2)], idx[int(d2)]],
                          0))

    n_verts = sum(1 for c in cells if c[0] == _VTK_VERTEX)
    n_lines = sum(1 for c in cells if c[0] == _VTK_LINE)
    n_tris = sum(1 for c in cells if c[0] == _VTK_TRIANGLE)
    n_tets = sum(1 for c in cells if c[0] == _VTK_TETRA)
    print(f"[tet_prep] VTK cells: {n_verts} vertex, {n_lines} line, "
          f"{n_tris} triangle, {n_tets} tetra")

    with open(out_path, "w") as f:
        f.write("# vtk DataFile Version 2.0\n")
        f.write("dp3d T1_9 AlgoHex tetmesh (feature-colored)\n")
        f.write("ASCII\nDATASET UNSTRUCTURED_GRID\n")
        f.write(f"POINTS {len(coords)} float\n")
        for p in coords:
            f.write(f"{p[0]} {p[1]} {p[2]}\n")
        total_ints = sum(1 + len(pts) for _, pts, _ in cells)
        f.write(f"CELLS {len(cells)} {total_ints}\n")
        for _, pts, _ in cells:
            f.write(f"{len(pts)} " + " ".join(str(p) for p in pts) + "\n")
        f.write(f"CELL_TYPES {len(cells)}\n")
        for vt, _, _ in cells:
            f.write(f"{vt}\n")
        f.write(f"CELL_DATA {len(cells)}\n")
        f.write("SCALARS color int 1\nLOOKUP_TABLE default\n")
        for _, _, col in cells:
            f.write(f"{col}\n")
    print(f"[tet_prep] wrote {out_path} ({len(coords)} points, "
          f"{len(cells)} cells)")


def main():
    stl_path = extract_boundary_stl()
    gmsh = build_tet_mesh(stl_path)

    tag, coords, _ = gmsh.model.mesh.getNodes()
    tag = np.asarray(tag, np.int64)
    points = np.asarray(coords, float).reshape(-1, 3)
    idx = {int(t): i for i, t in enumerate(tag)}

    tets = []
    etypes, _etags, enodes = gmsh.model.mesh.getElements(3)
    for et, nds in zip(etypes, enodes):
        if et != 4:
            continue
        nds = np.asarray(nds, np.int64).reshape(-1, 4)
        for row in nds:
            tets.append([idx[int(v)] for v in row])
    gmsh.finalize()
    print(f"[tet_prep] {len(points)} points, {len(tets)} tets")

    btris = tet_boundary_faces(tets, points)
    print(f"[tet_prep] {len(btris)} tet-mesh boundary triangles")
    tri_ids = classify_boundary(points, btris)
    feat_edges, feat_verts = feature_edges_and_vertices(points, btris, tri_ids)

    out_vtk = MSH.parent / "T1_9_tet.vtk"
    write_algohex_vtk_analytic(points, tets, btris, tri_ids, feat_edges,
                               feat_verts, out_vtk)
    return out_vtk


if __name__ == "__main__":
    main()
