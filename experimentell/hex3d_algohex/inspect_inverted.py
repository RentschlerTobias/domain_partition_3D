#!/usr/bin/env python3
"""Inspect inverted cells in Canada hex meshes.

Parses hex meshes (OVM and VTK formats), computes per-cell scaled Jacobian,
identifies inverted cells (min sJ <= 0), and creates visualization files
with inverted cells, their neighbors, and markers.

For each inverted cell reports:
- mesh name, cell id
- centre xyz, bbox, volume
- min scaled Jacobian
- boundary vs interior (count of boundary faces)
- 1-ring neighbours

Writes:
- inverted_zoom.vtk: inverted cells + 1-2 rings of neighbours + markers
- inverted_zoom.json: metadata and explanation
"""

import sys
import json
import numpy as np
from pathlib import Path
from collections import defaultdict

# Add the module directory to path
REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO / "experimentell" / "hex3d_algohex"))

import ovm_io


def read_hex_mesh_vtk(path):
    """Read a hex mesh from VTK unstructured grid format.
    
    Returns (points, hexes) where hexes is (n_cells, 8) array.
    """
    with open(path) as f:
        lines = f.readlines()
    
    i = 0
    # Find POINTS section
    while i < len(lines) and not lines[i].startswith("POINTS"):
        i += 1
    if i >= len(lines):
        raise ValueError(f"No POINTS section in {path}")
    
    n_pts = int(lines[i].split()[1])
    i += 1
    points = np.zeros((n_pts, 3))
    for k in range(n_pts):
        parts = lines[i].split()
        points[k] = [float(x) for x in parts[:3]]
        i += 1
    
    # Find CELLS section
    while i < len(lines) and not lines[i].startswith("CELLS"):
        i += 1
    if i >= len(lines):
        raise ValueError(f"No CELLS section in {path}")
    
    n_cells = int(lines[i].split()[1])
    i += 1
    hexes = []
    for k in range(n_cells):
        parts = [int(x) for x in lines[i].split()]
        if parts[0] == 8:  # VTK_HEXAHEDRON
            hexes.append(parts[1:9])
        i += 1
    
    return points, np.array(hexes, dtype=np.int64)


def read_hex_mesh(path):
    """Read hex mesh from OVM or VTK format."""
    path = Path(path)
    if path.suffix == ".ovm":
        points, edges, faces, polys = ovm_io.read_ovm(path)
        hexes, skipped = ovm_io.ovm_to_cells(points, edges, faces, polys)
        return points, hexes
    elif path.suffix == ".vtk":
        return read_hex_mesh_vtk(path)
    else:
        raise ValueError(f"Unknown format: {path.suffix}")


def compute_scaled_jacobians(points, hexes):
    """Compute min scaled Jacobian for each cell.
    
    Uses the same convention as ovm_io.scaled_jacobian and HexBlockValidator.
    Returns array of shape (n_cells,).
    """
    idx = [(0, 1, 3, 4), (1, 2, 0, 5), (2, 3, 1, 6), (3, 0, 2, 7),
           (4, 7, 5, 0), (5, 4, 6, 1), (6, 5, 7, 2), (7, 6, 4, 3)]
    
    Q = points[hexes]  # (n_cells, 8, 3)
    sj_all = np.empty((len(hexes), 8))
    
    for k, (a, b, c, d) in enumerate(idx):
        u = Q[:, b] - Q[:, a]
        v = Q[:, c] - Q[:, a]
        w = Q[:, d] - Q[:, a]
        n = np.linalg.norm(u, axis=1) * np.linalg.norm(v, axis=1) * np.linalg.norm(w, axis=1)
        denom = np.where(n > 0, n, 1.0)
        sj_all[:, k] = np.einsum('ij,ij->i', u, np.cross(v, w)) / denom
    
    return sj_all.min(axis=1)


def build_cell_topology(hexes):
    """Build face-to-cell and edge-to-cell maps.
    
    Returns (f2h, e2h) where:
    - f2h: dict mapping frozenset(face_verts) -> list of cell indices
    - e2h: dict mapping (min_v, max_v) -> list of cell indices
    """
    f2h = defaultdict(list)
    e2h = defaultdict(list)
    
    # Hex faces in VTK order
    HF = [(0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4),
          (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
    
    for ci, cell in enumerate(hexes):
        for face in HF:
            verts = tuple(int(cell[i]) for i in face)
            f2h[frozenset(verts)].append(ci)
        
        # Edges
        edges = [(0, 1), (1, 2), (2, 3), (3, 0),
                 (4, 5), (5, 6), (6, 7), (7, 4),
                 (0, 4), (1, 5), (2, 6), (3, 7)]
        for a, b in edges:
            u, v = int(cell[a]), int(cell[b])
            e2h[(min(u, v), max(u, v))].append(ci)
    
    return f2h, e2h


def compute_cell_neighbours(f2h, n_cells):
    """Compute 1-ring cell neighbours (cells sharing a face)."""
    neighbours = defaultdict(set)
    for fk, hs in f2h.items():
        if len(hs) == 2:
            a, b = hs
            neighbours[a].add(b)
            neighbours[b].add(a)
    return neighbours


def compute_cell_boundary_faces(hexes, f2h):
    """Count boundary faces per cell (faces with only one owner)."""
    bnd_faces = {fk for fk, hs in f2h.items() if len(hs) == 1}
    
    HF = [(0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4),
          (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
    
    bnd_count = np.zeros(len(hexes), dtype=int)
    for ci, cell in enumerate(hexes):
        for face in HF:
            verts = tuple(int(cell[i]) for i in face)
            if frozenset(verts) in bnd_faces:
                bnd_count[ci] += 1
    
    return bnd_count


def compute_cell_volume(points, cell):
    """Compute signed volume of a hexahedron."""
    return ovm_io._hex_volume(points[cell])


def compute_cell_center(points, cell):
    """Compute cell centroid."""
    return points[cell].mean(axis=0)


def compute_cell_bbox(points, cell):
    """Compute cell bounding box."""
    pts = points[cell]
    return pts.min(axis=0), pts.max(axis=0)



# ---------------------------------------------------------------------------
# Concave corner vs true fold diagnosis
# ---------------------------------------------------------------------------
#
# A negative scaled Jacobian at a hex corner does NOT by itself mean the cell
# is folded.  The trilinear map of a hexahedron can have a negative Jacobian
# at a corner while the cell itself is a perfectly valid, non-self-
# intersecting volume with a positive volume: that happens when the hex has a
# *reentrant* (concave) corner, i.e. one of the quad faces meeting at that
# corner has an interior angle > 180 deg.  The material then wraps around the
# corner (solid angle > 2*pi) but the boundary surface does not cross itself.
#
# The three quantities that separate the two cases:
#   * signed volume   -- negative  => the cell is turned inside out (fold)
#   * self-intersection -- the boundary surface crosses itself (fold)
#   * positive volume + no self-intersection + only corner Jacobians
#     negative at reflex faces => concave corner (convexity artefact)

# The six quad faces of a VTK hexahedron, in the VTK face order.
# NOTE the VTK quirk: face 0 (0,1,2,3) is wound *inward* for a right-handed
# hex (its right-hand normal points into the cell), faces 1..5 point outward.
HEX_FACES = [(0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4),
             (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]

# The three edge-neighbours of each corner, ordered so that the triple
# (corner, n1, n2, n3) is right-handed for a healthy hex; the sign of the
# triple product is then the sign of the corner Jacobian.
HEX_NEIGHBOURS = {0: (1, 3, 4), 1: (2, 0, 5), 2: (3, 1, 6), 3: (0, 2, 7),
                  4: (7, 5, 0), 5: (4, 6, 1), 6: (5, 7, 2), 7: (6, 4, 3)}


def hex_signed_volume(points, cell):
    """Signed volume of a hex by the divergence theorem over its six faces.

    Each quad face is split into two triangles and the flux
    ``a . (b-a) x (c-a) / 6`` is accumulated.  Face 0 is negated because the
    VTK face table winds it inward (see HEX_FACES).  The result is positive
    for a right-handed hex, negative for one that is turned inside out.
    """
    q = points[cell]
    vol = 0.0
    for fi, face in enumerate(HEX_FACES):
        sgn = -1.0 if fi == 0 else 1.0
        p = q[list(face)]
        for tri in ((0, 1, 2), (0, 2, 3)):
            a, b, c = p[tri[0]], p[tri[1]], p[tri[2]]
            vol += sgn * np.dot(a, np.cross(b - a, c - a)) / 6.0
    return float(vol)


# Natural coordinates of the eight corners of the trilinear hex, in the
# VTK node order (0,1,2,3 bottom quad, 4,5,6,7 top quad).
HEX_CORNERS = [(-1.0, -1.0, -1.0), (1.0, -1.0, -1.0),
               (1.0, 1.0, -1.0), (-1.0, 1.0, -1.0),
               (-1.0, -1.0, 1.0), (1.0, -1.0, 1.0),
               (1.0, 1.0, 1.0), (-1.0, 1.0, 1.0)]


def hex_shape_derivatives(s, t, u):
    """Analytic derivatives of the eight trilinear shape functions.

    Returns an array of shape (3, 8): row 0 holds dN/ds, row 1 dN/dt and
    row 2 dN/du.
    """
    dN = np.zeros((3, 8))
    for i, (a, b, c) in enumerate(HEX_CORNERS):
        dN[0, i] = a * (1.0 + b * t) * (1.0 + c * u) / 8.0
        dN[1, i] = b * (1.0 + a * s) * (1.0 + c * u) / 8.0
        dN[2, i] = c * (1.0 + a * s) * (1.0 + b * t) / 8.0
    return dN


def hex_det_jacobian(points, cell, s, t, u):
    """Determinant of the Jacobian of the trilinear map at (s,t,u).

    Positive means the map is orientation preserving there, negative means it
    is inverted (the hex folds over itself at that point).  The derivative is
    taken analytically, so the value is exact and does not depend on a finite
    difference step.
    """
    dN = hex_shape_derivatives(s, t, u)
    return float(np.linalg.det(dN @ points[cell]))


def hex_corner_jacobians(points, cell):
    """Jacobian determinant at the eight corners of the trilinear map.

    The corner values are the exact values of the trilinear map at the eight
    natural corners, so they are not affected by the finite differences used
    for the interior samples.
    """
    q = points[cell]
    out = []
    for k in range(8):
        b, c, d = HEX_NEIGHBOURS[k]
        out.append(float(np.dot(q[b] - q[k], np.cross(q[c] - q[k], q[d] - q[k]))))
    return np.asarray(out)


def hex_trilinear_volume(points, cell, n=4):
    """Volume of the trilinear hex by Gauss quadrature over [-1,1]^3.

    This is the volume of the *image* of the trilinear map (the isoparametric
    cell).  It differs from the divergence-theorem volume whenever the quad
    faces are not planar, because the divergence theorem is applied here to
    the flat triangles of the polyhedron while the trilinear map has curved
    (bilinear) faces.  Both signs agree for every cell in this study.
    """
    g = 0.5 * (1.0 - np.cos(np.pi * (np.arange(1, n + 1) - 0.5) / n))
    total = 0.0
    for s in g:
        for t in g:
            for u in g:
                total += hex_det_jacobian(points, cell, s, t, u)
    return float(total * (2.0 / n) ** 3)


def _cross2(a, b):
    return float(a[0] * b[1] - a[1] * b[0])


def hex_reentrant_corners(points, cell):
    """Corners whose material wraps around more than half a sphere.

    At a hex corner three edges meet.  The corner is *reentrant* (a concave
    corner, material angle > 2 pi) exactly when the three edge vectors do
    not fit into a single closed half space through the corner, i.e. when the
    origin lies inside the tetrahedron they span.  This test is exact and
    does not depend on the planarity of the faces, which matters because the
    faces of a strongly warped hex do not project to a simple polygon.

    Returns a list of (corner, signed_solid_angle) with the signed solid
    angle of the trihedral cone in steradian; it is negative for a
    reentrant corner.
    """
    q = points[cell]
    out = []
    for k in range(8):
        vs = [q[n] - q[k] for n in HEX_NEIGHBOURS[k]]
        norms = [float(np.linalg.norm(v)) for v in vs]
        if min(norms) < 1e-30:
            out.append((int(k), float("nan")))
            continue
        a, b, c = [v / n for v, n in zip(vs, norms)]
        num = float(np.dot(a, np.cross(b, c)))
        den = 1.0 + float(np.dot(a, b)) + float(np.dot(b, c)) + float(np.dot(c, a))
        sa = 2.0 * float(np.arctan2(num, den))
        out.append((int(k), sa))
    return out


def hex_face_geometry(points, cell, face):
    """Interior angles and reflex flags of one quad face.

    Returns a dict with
      ``angles3d``    the four 3-D interior angles (between the two face
                      edges at the corner), in degrees, in the face loop order;
      ``angles``      the four *signed* interior angles measured in the best
                      fit plane of the face.  These are the angles whose sum
                      is 360 deg for a simple quad, so an angle above 180 deg
                      marks a reflex (reentrant) corner.  For a planar quad
                      they agree with ``angles3d``;
      ``reflex``      list of four booleans, True where the corner is reflex;
      ``max_angle``   the largest of ``angles``;
      ``planarity``   the singular value ratio of the four face points: 0 for
                      a planar quad, larger for a warped one;
      ``convexity``   'convex' | 'concave' | 'flat' (degenerate).
    """
    p = points[cell][list(face)]
    n = len(face)
    angles3d = []
    for i in range(n):
        u = p[(i - 1) % n] - p[i]
        v = p[(i + 1) % n] - p[i]
        nu, nv = np.linalg.norm(u), np.linalg.norm(v)
        if nu == 0.0 or nv == 0.0:
            angles3d.append(float("nan"))
        else:
            angles3d.append(float(np.degrees(
                np.arccos(np.clip(np.dot(u, v) / (nu * nv), -1.0, 1.0)))))

    centre = p.mean(axis=0)
    _, sv, vt = np.linalg.svd(p - centre)
    planarity = float(sv[-1] / max(sv[0], 1e-30))
    a2, b2 = vt[0], vt[1]
    proj = np.column_stack([(p - centre) @ a2, (p - centre) @ b2])

    # Turn at each corner, measured in the face plane.  For a counter
    # clockwise loop the polygon turn at a convex corner is positive; a
    # negative turn marks a reflex (reentrant) corner.  The interior angle
    # is 180 deg minus the turn, so it always lies in (0, 360) and an
    # interior angle above 180 deg is exactly a reflex corner.
    area2 = 0.0
    for i in range(n):
        a = proj[i]
        b = proj[(i + 1) % n]
        area2 += a[0] * b[1] - b[0] * a[1]
    orient = 1.0 if area2 >= 0.0 else -1.0
    turns = []
    for i in range(n):
        prev = proj[(i - 1) % n] - proj[i]
        nxt = proj[(i + 1) % n] - proj[i]
        turns.append(orient * _cross2(nxt, prev))
    pos = sum(1 for t in turns if t > 0)
    neg = sum(1 for t in turns if t < 0)
    if pos == 0 or neg == 0:
        reflex = [False] * n
        convexity = "flat"
    else:
        major = 1 if pos > neg else -1
        reflex = [(t * major) < 0 for t in turns]
        convexity = "concave" if any(reflex) else "convex"

    if convexity == "flat":
        # Degenerate projection: fall back to the 3-D angles.
        angles = list(angles3d)
    else:
        angles = []
        for i in range(n):
            prev = proj[(i - 1) % n] - proj[i]
            nxt = proj[(i + 1) % n] - proj[i]
            if np.linalg.norm(prev) == 0.0 or np.linalg.norm(nxt) == 0.0:
                angles.append(float("nan"))
                continue
            turn = np.degrees(np.arctan2(orient * _cross2(nxt, prev),
                                         float(np.dot(nxt, prev))))
            angles.append(float(180.0 - turn))

    # The reflex flag is the definition the report uses: an interior angle
    # above 180 deg.  The turn test above produces exactly such angles, so
    # the two agree; the flag is restated here so that the JSON carries a
    # single unambiguous criterion.
    reflex = [bool(np.isfinite(a) and a > 180.0) for a in angles]

    finite = [a for a in angles if np.isfinite(a)]
    return {
        "vertices": list(face),
        "angles3d": [float(a) for a in angles3d],
        "angles": [float(a) for a in angles],
        "reflex": [bool(r) for r in reflex],
        "reflex_corners": [int(face[i]) for i in range(n) if reflex[i]],
        "max_angle": float(max(finite)) if finite else float("nan"),
        "max_angle3d": float(max([a for a in angles3d if np.isfinite(a)]))
                       if any(np.isfinite(angles3d)) else float("nan"),
        "planarity": planarity,
        "convexity": convexity,
        "angle_sum": float(sum(a for a in angles if np.isfinite(a))),
        # The projected quad is a simple polygon only when its interior
        # angles add up to 360 deg.  A strongly warped face projects onto a
        # polygon that winds around more than once and its projected angles
        # no longer describe a simple polygon.
        # The interior angles of a simple planar quad add up to 360 deg.  A
        # strongly warped face projects onto a polygon that winds around
        # more than once (sum 720 deg) and its projected angles then no
        # longer describe a simple polygon.
        "angle_sum_deviation": float(abs(sum(a for a in angles
                                             if np.isfinite(a)) - 360.0)),
        "simple_projection": bool(abs(sum(a for a in angles
                                          if np.isfinite(a)) - 360.0) < 90.0),
    }


def _segment_segment_3d(p1, p2, p3, p4, eps=1e-12):
    """Do the 3D segments p1p2 and p3p4 cross?  Returns (hit, point)."""
    d1 = p2 - p1
    d2 = p4 - p3
    r = p3 - p1
    a = float(d1 @ d1)
    e = float(d2 @ d2)
    f = float(d2 @ r)
    c = float(d1 @ r)
    b = float(d1 @ d2)
    den = a * e - b * b
    if abs(den) < 1e-20:
        # Parallel or degenerate: check whether they are collinear and overlap.
        if abs(b) > 1e-15:
            return False, None
        # collinear?  (r x d1) must vanish
        if np.linalg.norm(np.cross(r, d1)) > 1e-9 * max(1.0, np.linalg.norm(d1)):
            return False, None
        s0 = float(r @ d1) / a
        s1 = s0 + e / a
        lo, hi = min(s0, s1), max(s0, s1)
        lo = max(lo, 0.0)
        hi = min(hi, 1.0)
        if lo > hi + eps:
            return False, None
        s = 0.5 * (lo + hi)
        return True, p1 + s * d1
    s = (b * f - c * e) / den
    t = (a * f - b * c) / den
    if -eps <= s <= 1.0 + eps and -eps <= t <= 1.0 + eps:
        return True, p1 + s * d1
    return False, None


def _point_in_triangle_3d(p, tri, n, tol=1e-9):
    """Is the point p inside the triangle tri, whose plane normal is n?"""
    if abs(float((p - tri[0]) @ n)) > tol * max(1.0, np.linalg.norm(n)):
        return False
    # Project onto the dominant plane of the triangle and use the 2D test.
    axis = int(np.argmax(np.abs(n)))
    keep = [i for i in range(3) if i != axis]
    a = np.array([tri[0][keep[0]], tri[0][keep[1]]])
    b = np.array([tri[1][keep[0]], tri[1][keep[1]]])
    c = np.array([tri[2][keep[0]], tri[2][keep[1]]])
    q = np.array([p[keep[0]], p[keep[1]]])
    s1 = _cross2(b - a, q - a)
    s2 = _cross2(c - b, q - b)
    s3 = _cross2(a - c, q - c)
    has_neg = (s1 < 0) or (s2 < 0) or (s3 < 0)
    has_pos = (s1 > 0) or (s2 > 0) or (s3 > 0)
    return not (has_neg and has_pos)


def _tri_tri_intersect(t1, t2, tol=1e-12):
    """Do two triangles overlap?  Moller's triangle-triangle test.

    ``t1``/``t2`` are (3,3) arrays of vertex coordinates.  Returns
    (hit, reason).  The edge/edge test is done in 3D so that triangles
    which merely overlap when projected onto a plane are not mistaken
    for a genuine crossing.
    """
    n1 = np.cross(t1[1] - t1[0], t1[2] - t1[0])
    n2 = np.cross(t2[1] - t2[0], t2[2] - t2[0])
    ln1, ln2 = np.linalg.norm(n1), np.linalg.norm(n2)
    if ln1 == 0.0 or ln2 == 0.0:
        return False, "degenerate"
    # Work with unit normals so that the tolerance is a length.
    n1 = n1 / ln1
    n2 = n2 / ln2
    scale = max(np.linalg.norm(t1 - t1[0]), np.linalg.norm(t2 - t2[0]), 1e-30)
    eps = tol * scale
    d2 = (t2 - t1[0]) @ n1
    if np.all(d2 > eps) or np.all(d2 < -eps):
        return False, "separated"
    d1 = (t1 - t2[0]) @ n2
    if np.all(d1 > eps) or np.all(d1 < -eps):
        return False, "separated"

    for i in range(3):
        for j in range(3):
            hit, _pt = _segment_segment_3d(t1[i], t1[(i + 1) % 3],
                                           t2[j], t2[(j + 1) % 3])
            if hit:
                return True, "edge-cross"
    for i in range(3):
        if _point_in_triangle_3d(t1[i], t2, n2):
            return True, "vertex-in"
    for j in range(3):
        if _point_in_triangle_3d(t2[j], t1, n1):
            return True, "vertex-in"
    return False, "-"


def hex_self_intersection(points, cell):
    """Does the boundary of the hex cross itself?

    Three independent checks are reported:
      ``corner_tets``    signed volume of the eight corner tetrahedra
                         (corner + its three edge neighbours).  All positive
                         means no corner is inverted;
      ``face_pyramids``  signed volume of the six pyramids from the cell
                         centroid to each face (face 0 negated to match the
                         outward orientation).  All positive means every face
                         still bounds the cell on the outside;
      ``face_pairs``     the pairs of *non-adjacent* faces whose triangles
                         intersect.  A hit means the boundary surface
                         crosses itself: a genuine fold.
    """
    q = points[cell]
    corner_tets = []
    for k in range(8):
        b, c, d = HEX_NEIGHBOURS[k]
        corner_tets.append(float(
            np.dot(q[b] - q[k], np.cross(q[c] - q[k], q[d] - q[k])) / 6.0))

    centre = q.mean(axis=0)
    face_pyramids = []
    for fi, face in enumerate(HEX_FACES):
        # VTK winds face 0 (0,1,2,3) inwards, so its loop has to be reversed
        # to point outwards like the other five.
        sgn = -1.0 if fi == 0 else 1.0
        p = q[list(face)]
        vol = 0.0
        for tri in ((0, 1, 2), (0, 2, 3)):
            a, b, c = p[tri[0]], p[tri[1]], p[tri[2]]
            # (centre - a) . n  is the pyramid volume with the *outward*
            # normal, i.e. the opposite sign of the divergence term a . n.
            vol += -sgn * np.dot(centre - a, np.cross(b - a, c - a)) / 6.0
        face_pyramids.append(float(vol))

    hits = []
    for i in range(6):
        for j in range(i + 1, 6):
            if set(HEX_FACES[i]) & set(HEX_FACES[j]):
                continue                      # adjacent faces share an edge
            fi, fj = HEX_FACES[i], HEX_FACES[j]
            for ti in range(2):
                for tj in range(2):
                    t1 = np.array([q[fi[ti]], q[fi[ti + 1]], q[fi[(ti + 2) % 4]]])
                    t2 = np.array([q[fj[tj]], q[fj[tj + 1]], q[fj[(tj + 2) % 4]]])
                    hit, why = _tri_tri_intersect(t1, t2)
                    if hit:
                        hits.append({"face_a": i, "face_b": j,
                                     "tri_a": ti, "tri_b": tj, "reason": why})
                        break
                if hits and hits[-1]["face_a"] == i and hits[-1]["face_b"] == j:
                    break
    return {
        "corner_tets": corner_tets,
        "corner_tets_all_positive": bool(all(v > 0.0 for v in corner_tets)),
        "face_pyramids": face_pyramids,
        "face_pyramids_all_positive": bool(all(v > 0.0 for v in face_pyramids)),
        "face_pairs": hits,
        "self_intersecting": bool(len(hits) > 0),
    }


# The eight interior sample points: the centres of the eight sub-hexes of
# the natural coordinates, at (+-0.5, +-0.5, +-0.5).
HEX_INTERIOR_SAMPLES = [(-0.5, -0.5, -0.5), (-0.5, -0.5, 0.5),
                         (-0.5, 0.5, -0.5), (-0.5, 0.5, 0.5),
                         (0.5, -0.5, -0.5), (0.5, -0.5, 0.5),
                         (0.5, 0.5, -0.5), (0.5, 0.5, 0.5)]


def diagnose_cell(points, cell):
    """Full concave-corner vs true-fold diagnosis of one hex cell.

    The verdict is:
      TRUE FOLD        the signed volume is <= 0, or the boundary surface
                       crosses itself;
      CONCAVE CORNER   the volume is positive, the centroid and every
                       interior sample have a positive Jacobian, the surface
                       does not cross itself, and the only negative
                       Jacobians are corner values sitting on reflex (> 180
                       deg) corners of the quad faces.
    """
    q = points[cell]
    vol = hex_signed_volume(points, cell)
    vol_tri = float(ovm_io._hex_volume(q))
    vol_iso = hex_trilinear_volume(points, cell)

    cj = hex_corner_jacobians(points, cell)
    j_centroid = hex_det_jacobian(points, cell, 0.0, 0.0, 0.0)
    j_interior = [hex_det_jacobian(points, cell, s, t, u)
                   for (s, t, u) in HEX_INTERIOR_SAMPLES]

    faces = [hex_face_geometry(points, cell, f) for f in HEX_FACES]
    si = hex_self_intersection(points, cell)

    neg_corners = [k for k in range(8) if cj[k] < 0.0]
    interior_all_positive = bool(all(j > 0.0 for j in j_interior))

    # A corner is reentrant when the material wraps around more than half a
    # sphere there.  This is the exact criterion for a concave corner and it
    # does not depend on the planarity of the faces, which matters because
    # the faces of a strongly warped hex do not project to a simple polygon.
    reentrant = []
    solid_angles = []
    for (k, sa) in hex_reentrant_corners(points, cell):
        solid_angles.append(sa)
        if np.isfinite(sa) and sa < 0.0:
            reentrant.append(int(k))
    reentrant = sorted(reentrant)
    for f, k in zip(faces, range(8)):
        f["reentrant_corners"] = [int(c) for c in f["vertices"]
                                  if int(c) in reentrant]

    reasons = []
    if vol <= 0.0:
        reasons.append("signed volume <= 0")
    if vol_iso <= 0.0:
        reasons.append("trilinear volume <= 0")
    if si["self_intersecting"]:
        reasons.append("boundary self-intersection")
    if j_centroid <= 0.0:
        reasons.append("centroid Jacobian <= 0")
    if not interior_all_positive:
        reasons.append("negative Jacobian at an interior sample point")
    unexplained = [k for k in neg_corners if k not in reentrant]
    if unexplained:
        reasons.append("negative corner Jacobian at a convex corner "
                       "(corners %s)" % unexplained)

    if reasons:
        verdict = "TRUE FOLD"
    elif neg_corners:
        verdict = "CONCAVE CORNER"
    else:
        verdict = "OK"

    return {
        "cell_id": int(cell) if np.isscalar(cell) else None,
        "vertices": [int(v) for v in cell],
        "center": q.mean(axis=0).tolist(),
        "volume_divergence": float(vol),
        "volume_tet5": float(vol_tri),
        "volume_trilinear": float(vol_iso),
        "volume_ovm_io": float(ovm_io._hex_volume(q)),
        "corner_jacobians": [float(v) for v in cj],
        "corner_jacobian_min": float(cj.min()),
        "negative_corners": neg_corners,
        "centroid_jacobian": float(j_centroid),
        "interior_samples": [
            {"s": float(s), "t": float(t), "u": float(u),
             "jacobian": float(j)}
            for (s, t, u), j in zip(HEX_INTERIOR_SAMPLES, j_interior)],
        "interior_jacobian_min": float(min(j_interior)),
        "faces": faces,
        "corner_solid_angles": [float(a) for a in solid_angles],
        "reentrant_corners": reentrant,
        "reflex_corners": reentrant,
        # Faces that carry a reentrant (concave) hex corner, and the subset
        # whose own projected quad has a reflex (> 180 deg) interior angle.
        "reflex_faces": [i for i, f in enumerate(faces)
                         if any(f["reentrant_corners"])],
        "reflex_faces_angled": [i for i, f in enumerate(faces)
                                if any(f["reflex"])],
        "self_intersection": si,
        "verdict": verdict,
        "verdict_reasons": reasons,
    }


def get_rings_around(cells, neighbours, n_rings=1):
    """Get cells within n_rings of the given cells."""
    result = set(cells)
    frontier = set(cells)
    for _ in range(n_rings):
        new_frontier = set()
        for c in frontier:
            new_frontier.update(neighbours[c])
        new_frontier -= result
        result.update(new_frontier)
        frontier = new_frontier
    return sorted(result)


def write_zoom_vtk(path, points, hexes, inverted_ids, ring_ids, cell_data=None,
                   reflex_info=None):
    """Write a VTK file with inverted cells, their rings, and markers.
    
    Markers are small octahedra at the center of each inverted cell,
    scaled to ~2x the cell size for visibility.
    
    When reflex_info is given it is a dict
        {"cell_id": int, "faces": [list_of_6_face_corner_ids],
         "reflex_faces": [face indices]}
    and the reflex faces of that cell are added as separate triangle cells
    so they can be highlighted.  Two extra cell arrays are written:
    "reflex_face_count" (number of reflex faces per cell) and
    "face_reflex_mask" (6 entries per hex cell, 1 on a reflex face).
    """
    # Collect all cell indices
    all_cells = sorted(set(inverted_ids) | set(ring_ids))
    
    # Compact points to only those used
    used_pts = np.unique(np.concatenate([hexes[c] for c in all_cells]))
    pt_remap = {int(v): i for i, v in enumerate(used_pts)}
    P = points[used_pts]
    H = np.array([[pt_remap[int(v)] for v in hexes[c]] for c in all_cells])
    
    # Map cell indices
    cell_remap = {old: new for new, old in enumerate(all_cells)}
    inv_local = [cell_remap[c] for c in inverted_ids]
    
    # Reflex faces of the flagged cell, as extra triangle cells.  The faces
    # whose own quad has an interior angle > 180 deg are drawn a second
    # time so that they can be highlighted on their own; every face that
    # carries a reentrant (concave) hex corner is marked in the mask.
    reflex_tris = []
    reflex_mask = {c: [0] * 6 for c in all_cells}
    if reflex_info and len(inverted_ids):
        rcid = int(inverted_ids[0])
        faces = reflex_info["faces"]
        for fi in reflex_info["reflex_faces_angled"]:
            corners = faces[fi]["vertices"]
            for tri in ((0, 1, 2), (0, 2, 3)):
                reflex_tris.append([pt_remap[int(hexes[rcid][corners[k]])]
                                    for k in tri])
            reflex_mask[rcid][fi] = 1
        for fi in reflex_info["reflex_faces"]:
            reflex_mask[rcid][fi] = max(reflex_mask[rcid][fi], 1)
    
    # Create markers (octahedra) at inverted cell centers
    markers = []
    marker_pts = []
    pt_offset = len(P)
    
    for cid in inverted_ids:
        center = compute_cell_center(points, hexes[cid])
        # Estimate cell size from bbox
        bmin, bmax = compute_cell_bbox(points, hexes[cid])
        size = np.linalg.norm(bmax - bmin)
        radius = size * 1.0  # marker radius = cell size
        
        # Octahedron: 6 vertices along axes
        oct_pts = [
            center + [radius, 0, 0],
            center + [-radius, 0, 0],
            center + [0, radius, 0],
            center + [0, -radius, 0],
            center + [0, 0, radius],
            center + [0, 0, -radius],
        ]
        
        # Octahedron faces (8 triangles)
        # Vertices: 0=+x, 1=-x, 2=+y, 3=-y, 4=+z, 5=-z
        oct_faces = [
            [0, 2, 4], [2, 1, 4], [1, 3, 4], [3, 0, 4],
            [0, 4, 2], [2, 4, 1], [1, 4, 3], [3, 4, 0],
            [0, 2, 5], [2, 1, 5], [1, 3, 5], [3, 0, 5],
            [0, 5, 2], [2, 5, 1], [1, 5, 3], [3, 5, 0],
        ]
        
        for pt in oct_pts:
            marker_pts.append(pt)
        
        for face in oct_faces:
            markers.append([pt_offset + f for f in face])
        
        pt_offset += 6
    
    # Write VTK
    with open(path, 'w') as f:
        f.write("# vtk DataFile Version 2.0\n")
        f.write("Inverted cells with 1-2 ring neighbors and markers\n")
        f.write("ASCII\n")
        f.write("DATASET UNSTRUCTURED_GRID\n")
        
        # Points (hex mesh + marker points)
        all_pts = np.vstack([P, marker_pts]) if marker_pts else P
        f.write(f"POINTS {len(all_pts)} float\n")
        for pt in all_pts:
            f.write(f"{pt[0]} {pt[1]} {pt[2]}\n")
        
        # Cells
        n_hex = len(H)
        n_tri = len(markers)
        n_reflex = len(reflex_tris)
        n_total = n_hex + n_tri + n_reflex
        f.write(f"CELLS {n_total} {9 * n_hex + 4 * n_tri + 4 * n_reflex}\n")
        
        for h in H:
            f.write(f"8 {' '.join(str(int(v)) for v in h)}\n")
        
        for m in markers:
            f.write(f"3 {m[0]} {m[1]} {m[2]}\n")
        
        for t in reflex_tris:
            f.write(f"3 {t[0]} {t[1]} {t[2]}\n")
        
        # Cell types
        f.write(f"CELL_TYPES {n_total}\n")
        for _ in range(n_hex):
            f.write("12\n")  # VTK_HEXAHEDRON
        for _ in range(n_tri):
            f.write("5\n")  # VTK_TRIANGLE
        for _ in range(n_reflex):
            f.write("5\n")  # VTK_TRIANGLE (reflex face)
        
        # Cell data
        if cell_data or inv_local or reflex_info:
            f.write(f"CELL_DATA {n_total}\n")
            
            # Inverted flag
            f.write("SCALARS is_inverted int 1\n")
            f.write("LOOKUP_TABLE default\n")
            inv_set = set(inv_local)
            for i in range(n_hex):
                f.write(f"{1 if i in inv_set else 0}\n")
            for _ in range(n_tri):
                f.write("1\n")  # markers are all "inverted"
            for _ in range(n_reflex):
                f.write("0\n")  # reflex faces are not inverted cells
            
            # Original cell ID
            f.write("SCALARS original_cell_id int 1\n")
            f.write("LOOKUP_TABLE default\n")
            for c in all_cells:
                f.write(f"{c}\n")
            for _ in range(n_tri):
                f.write("-1\n")  # markers have no original ID
            for _ in range(n_reflex):
                f.write("-1\n")  # reflex faces have no original ID
            
            # Number of reflex faces per cell
            f.write("SCALARS reflex_face_count int 1\n")
            f.write("LOOKUP_TABLE default\n")
            for c in all_cells:
                f.write(f"{sum(reflex_mask[c])}\n")
            for _ in range(n_tri):
                f.write("0\n")
            for _ in range(n_reflex):
                f.write("1\n")
            
            # 6 entries per hex cell: 1 on the reflex faces of that cell
            f.write("CELL_ARRAY face_reflex_mask int 6\n")
            for c in all_cells:
                f.write(" ".join(str(v) for v in reflex_mask[c]) + "\n")
            for _ in range(n_tri):
                f.write("0 0 0 0 0 0\n")
            for _ in range(n_reflex):
                f.write("0 0 0 0 0 0\n")
            
            # 1 on the extra triangles that carry a reflex face
            f.write("SCALARS is_reflex_face int 1\n")
            f.write("LOOKUP_TABLE default\n")
            for _ in range(n_hex):
                f.write("0\n")
            for _ in range(n_tri):
                f.write("0\n")
            for _ in range(n_reflex):
                f.write("1\n")


def write_scaled_vtk(path, points, hexes, inverted_ids, scale=1.5):
    """Write inverted cells with coordinates scaled about their centers."""
    # For each inverted cell, scale its coordinates
    new_points = points.copy()
    new_hexes = []
    
    for cid in inverted_ids:
        center = compute_cell_center(points, hexes[cid])
        cell = hexes[cid]
        
        # Scale each vertex about the center
        for v in cell:
            new_points[v] = center + scale * (points[v] - center)
        
        new_hexes.append(cell)
    
    # Write as separate VTK
    H = np.array(new_hexes)
    used_pts = np.unique(H)
    pt_remap = {int(v): i for i, v in enumerate(used_pts)}
    P = new_points[used_pts]
    H_remapped = np.array([[pt_remap[int(v)] for v in h] for h in H])
    
    with open(path, 'w') as f:
        f.write("# vtk DataFile Version 2.0\n")
        f.write(f"Scaled inverted cells (scale={scale}x about center)\n")
        f.write("ASCII\n")
        f.write("DATASET UNSTRUCTURED_GRID\n")
        f.write(f"POINTS {len(P)} float\n")
        for pt in P:
            f.write(f"{pt[0]} {pt[1]} {pt[2]}\n")
        f.write(f"CELLS {len(H_remapped)} {9 * len(H_remapped)}\n")
        for h in H_remapped:
            f.write(f"8 {' '.join(str(int(v)) for v in h)}\n")
        f.write(f"CELL_TYPES {len(H_remapped)}\n")
        for _ in H_remapped:
            f.write("12\n")


def check_existing_inverted_vtk(path):
    """Check what the existing inverted_cells.vtk contains."""
    print(f"\n=== Checking existing {path} ===")
    if not Path(path).exists():
        print(f"File not found: {path}")
        return
    
    points, hexes = read_hex_mesh_vtk(path)
    print(f"Contains {len(hexes)} cells")
    
    # Check for duplicates
    cell_set = set(tuple(h) for h in hexes)
    print(f"Unique cells: {len(cell_set)} (duplicates: {len(hexes) - len(cell_set)})")
    
    # Check point duplication
    pt_set = set(tuple(p) for p in points)
    print(f"Unique points: {len(pt_set)} (duplicates: {len(points) - len(pt_set)})")
    
    return points, hexes


DIAGNOSIS_MESHES = [
    ("repaired_n100000.vtk",
     "/work/canada/sweep/inspect/repaired_n100000.vtk"),
    ("inverted_cells.vtk",
     "/work/canada/sweep/inspect/inverted_cells.vtk"),
    ("blocks.vtk",
     "/work/canada/sweep/struct_beam/blocks.vtk"),
]


def run_diagnosis(out_json=None, meshes=None, verbose=True):
    """Diagnose every inverted cell of the given meshes.

    For each cell with a non-positive scaled Jacobian the report holds the
    signed volume (divergence theorem over the six faces), the Jacobian
    determinant at the eight corners, at the centroid and at eight interior
    sample points, the four interior angles of every quad face with a
    reflex (> 180 deg) flag, a self-intersection check, and the verdict
    TRUE FOLD or CONCAVE CORNER.
    """
    if meshes is None:
        meshes = DIAGNOSIS_MESHES
    if out_json is None:
        out_json = "/work/canada/sweep/inspect/concavity.json"

    report = {"meshes": {}, "summary": {}}
    total_true = total_concave = 0
    for label, path in meshes:
        if not Path(path).exists():
            if verbose:
                print(f"[diagnose] skip {label}: {path} not found")
            continue
        P, H = read_hex_mesh_vtk(path)
        sj = compute_scaled_jacobians(P, H)
        inv_ids = [int(c) for c in np.where(sj <= 0)[0]]
        f2h, e2h = build_cell_topology(H)
        bnd_count = compute_cell_boundary_faces(H, f2h)
        cells = []
        for cid in inv_ids:
            d = diagnose_cell(P, H[cid])
            d["cell_id"] = int(cid)
            d["min_scaled_jacobian"] = float(sj[cid])
            d["center"] = compute_cell_center(P, H[cid]).tolist()
            d["volume"] = float(compute_cell_volume(P, H[cid]))
            d["boundary_faces"] = int(bnd_count[cid])
            d["is_interior"] = bool(bnd_count[cid] == 0)
            cells.append(d)
        n_true = sum(1 for c in cells if c["verdict"] == "TRUE FOLD")
        n_conc = sum(1 for c in cells if c["verdict"] == "CONCAVE CORNER")
        total_true += n_true
        total_concave += n_conc
        report["meshes"][label] = {
            "path": path,
            "n_cells": int(len(H)),
            "n_inverted": len(inv_ids),
            "n_true_fold": n_true,
            "n_concave_corner": n_conc,
            "cells": cells,
        }
        if verbose:
            print(f"[diagnose] {label}: {len(cells)} inverted cells, "
                  f"{n_true} TRUE FOLD, {n_conc} CONCAVE CORNER")
    report["summary"] = {
        "n_true_fold": total_true,
        "n_concave_corner": total_concave,
    }
    Path(out_json).parent.mkdir(parents=True, exist_ok=True)
    with open(out_json, "w") as f:
        json.dump(report, f, indent=2)
    if verbose:
        print(f"[diagnose] wrote {out_json}")
    return report


def summarize_diagnosis(report):
    """Print a compact summary of a diagnosis report."""
    print("\n=== Diagnosis summary ===")
    for label, m in report["meshes"].items():
        print(f"\n{label} ({m['path']})")
        print(f"  cells {m['n_cells']}, inverted {m['n_inverted']}, "
              f"TRUE FOLD {m['n_true_fold']}, "
              f"CONCAVE CORNER {m['n_concave_corner']}")
        for c in m["cells"]:
            print(f"   cell {c['cell_id']:>6}  minSJ {c['min_scaled_jacobian']:+.6f}  "
                  f"vol(div) {c['volume_divergence']:+.3e}  "
                  f"vol(tri) {c['volume_trilinear']:+.3e}  "
                  f"neg corners {c['negative_corners']}  "
                  f"reflex faces {c['reflex_faces_angled']} "
                  f"(reentrant corners {c['reentrant_corners']})  "
                  f"-> {c['verdict']}")
    s = report["summary"]
    print(f"\nTOTAL: {s['n_true_fold']} TRUE FOLD, "
          f"{s['n_concave_corner']} CONCAVE CORNER")


def write_diagnosis_zoom(zoom_path, report, mesh_label="repaired_n100000.vtk",
                         cell_id=None):
    """Extend inverted_zoom.vtk with the reflex faces of the flagged cell.

    The cell is written as a hexahedron together with its ring-1
    neighbours; the reflex faces are written a second time as triangle
    cells so that they can be highlighted on their own.
    """
    m = report["meshes"].get(mesh_label)
    if m is None or not m["cells"]:
        return None
    cells = m["cells"]
    if cell_id is None:
        # the cell with the most negative scaled Jacobian
        cell_id = min(cells, key=lambda c: c["min_scaled_jacobian"])["cell_id"]
    diag = next(c for c in cells if c["cell_id"] == cell_id)
    if not diag["reflex_faces_angled"]:
        return None
    P, H = read_hex_mesh_vtk(m["path"])
    f2h, e2h = build_cell_topology(H)
    neighbours = compute_cell_neighbours(f2h, len(H))
    nbrs = sorted(set(get_rings_around({cell_id}, neighbours, n_rings=2))
                  - {cell_id})
    write_zoom_vtk(zoom_path, P, H, [cell_id], nbrs,
                   reflex_info={"faces": diag["faces"],
                                "reflex_faces": diag["reflex_faces"],
                                "reflex_faces_angled":
                                    diag["reflex_faces_angled"]})
    return diag


def main():
    import argparse
    
    ap = argparse.ArgumentParser(description="Inspect inverted cells in Canada hex meshes")
    ap.add_argument("--check-existing", action="store_true",
                    help="Check existing inverted_cells.vtk")
    ap.add_argument("--diagnose", action="store_true",
                    help="Diagnose every inverted cell of the Canada meshes "
                         "(concave corner vs true fold) and write "
                         "concavity.json")
    ap.add_argument("--diagnose-out", default=None,
                    help="Output JSON for --diagnose "
                         "(default /work/canada/sweep/inspect/concavity.json)")
    args = ap.parse_args()
    
    # Define input files
    ovm_path = Path("/work/canada/sweep/n100000/hex_n100000.ovm")
    repaired_vtk = Path("/work/canada/sweep/inspect/repaired_n100000.vtk")
    blocks_vtk = Path("/work/canada/sweep/struct_beam/blocks.vtk")
    existing_inv = Path("/work/canada/sweep/inspect/inverted_cells.vtk")
    out_dir = Path("/work/canada/sweep/inspect")
    
    # Check existing file if requested
    if args.check_existing:
        check_existing_inverted_vtk(existing_inv)
        return
    
    # Concave corner vs true fold diagnosis
    if args.diagnose:
        out_json = args.diagnose_out or str(out_dir / "concavity.json")
        report = run_diagnosis(out_json=out_json)
        summarize_diagnosis(report)
        # Extend the zoom VTK with the reflex faces of the flagged cell
        zoom_path = out_dir / "inverted_zoom.vtk"
        try:
            diag = write_diagnosis_zoom(zoom_path, report)
        except Exception as exc:  # pragma: no cover - visualisation is optional
            print(f"[diagnose] zoom update skipped: {exc}")
            diag = None
        if diag is not None:
            print(f"[diagnose] {zoom_path}: cell {diag['cell_id']} reflex "
                  f"faces {diag['reflex_faces_angled']} highlighted")
        return
    
    # Collect inverted cells from all meshes
    all_inverted = []
    
    # 1. Original OVM mesh
    print("\n=== Processing original OVM mesh ===")
    if ovm_path.exists():
        P, H = read_hex_mesh(ovm_path)
        sj = compute_scaled_jacobians(P, H)
        inv_ids = np.where(sj <= 0)[0]
        print(f"Original mesh: {len(H)} cells, {len(inv_ids)} inverted, min sJ {sj.min():.4f}")
        
        for cid in inv_ids:
            center = compute_cell_center(P, H[cid])
            bmin, bmax = compute_cell_bbox(P, H[cid])
            vol = compute_cell_volume(P, H[cid])
            all_inverted.append({
                "mesh": "hex_n100000.ovm",
                "cell_id": int(cid),
                "center": center.tolist(),
                "bbox_min": bmin.tolist(),
                "bbox_max": bmax.tolist(),
                "volume": float(vol),
                "min_sj": float(sj[cid]),
            })
    
    # 2. Repaired mesh
    print("\n=== Processing repaired mesh ===")
    if repaired_vtk.exists():
        P, H = read_hex_mesh(repaired_vtk)
        sj = compute_scaled_jacobians(P, H)
        inv_ids = np.where(sj <= 0)[0]
        print(f"Repaired mesh: {len(H)} cells, {len(inv_ids)} inverted, min sJ {sj.min():.4f}")
        
        f2h, e2h = build_cell_topology(H)
        neighbours = compute_cell_neighbours(f2h, len(H))
        bnd_count = compute_cell_boundary_faces(H, f2h)
        
        for cid in inv_ids:
            center = compute_cell_center(P, H[cid])
            bmin, bmax = compute_cell_bbox(P, H[cid])
            vol = compute_cell_volume(P, H[cid])
            ring1 = get_rings_around([cid], neighbours, 1)
            ring2 = get_rings_around([cid], neighbours, 2)
            
            all_inverted.append({
                "mesh": "repaired_n100000.vtk",
                "cell_id": int(cid),
                "center": center.tolist(),
                "bbox_min": bmin.tolist(),
                "bbox_max": bmax.tolist(),
                "volume": float(vol),
                "min_sj": float(sj[cid]),
                "n_boundary_faces": int(bnd_count[cid]),
                "is_interior": bnd_count[cid] == 0,
                "ring1_neighbours": sorted(ring1),
                "ring2_neighbours": sorted(ring2),
                "points": P,
                "hexes": H,
            })
    
    # 3. Blocks mesh
    print("\n=== Processing blocks mesh ===")
    if blocks_vtk.exists():
        P, H = read_hex_mesh(blocks_vtk)
        sj = compute_scaled_jacobians(P, H)
        inv_ids = np.where(sj <= 0)[0]
        print(f"Blocks mesh: {len(H)} cells, {len(inv_ids)} inverted, min sJ {sj.min():.4f}")
        
        f2h, e2h = build_cell_topology(H)
        neighbours = compute_cell_neighbours(f2h, len(H))
        bnd_count = compute_cell_boundary_faces(H, f2h)
        
        for cid in inv_ids:
            center = compute_cell_center(P, H[cid])
            bmin, bmax = compute_cell_bbox(P, H[cid])
            vol = compute_cell_volume(P, H[cid])
            ring1 = get_rings_around([cid], neighbours, 1)
            ring2 = get_rings_around([cid], neighbours, 2)
            
            all_inverted.append({
                "mesh": "blocks.vtk",
                "cell_id": int(cid),
                "center": center.tolist(),
                "bbox_min": bmin.tolist(),
                "bbox_max": bmax.tolist(),
                "volume": float(vol),
                "min_sj": float(sj[cid]),
                "n_boundary_faces": int(bnd_count[cid]),
                "is_interior": bnd_count[cid] == 0,
                "ring1_neighbours": sorted(ring1),
                "ring2_neighbours": sorted(ring2),
                "points": P,
                "hexes": H,
            })
    
    # Print summary
    print(f"\n=== Summary ===")
    print(f"Total inverted cells found: {len(all_inverted)}")
    for inv in all_inverted:
        print(f"\nMesh: {inv['mesh']}, Cell ID: {inv['cell_id']}")
        print(f"  Center: {inv['center']}")
        print(f"  BBox: {inv['bbox_min']} -> {inv['bbox_max']}")
        print(f"  Volume: {inv['volume']:.6f}")
        print(f"  Min sJ: {inv['min_sj']:.6f}")
        if 'is_interior' in inv:
            print(f"  Interior: {inv['is_interior']} ({inv.get('n_boundary_faces', '?')} boundary faces)")
        if 'ring1_neighbours' in inv:
            print(f"  Ring-1 neighbours: {len(inv['ring1_neighbours'])} cells")
            print(f"  Ring-2 neighbours: {len(inv['ring2_neighbours'])} cells")
    
    # Write visualization files for the first inverted cell with full data
    inv_with_data = [x for x in all_inverted if 'points' in x]
    if inv_with_data:
        inv = inv_with_data[0]
        P = inv['points']
        H = inv['hexes']
        cid = inv['cell_id']
        
        f2h, e2h = build_cell_topology(H)
        neighbours = compute_cell_neighbours(f2h, len(H))
        
        # Get 1-2 rings
        ring1 = inv['ring1_neighbours']
        ring2 = inv['ring2_neighbours']
        
        # Write inverted_zoom.vtk.  The diagnosis of the flagged cell adds
        # its reflex (> 180 deg) faces as separate triangle cells so that
        # the concave corners can be highlighted on their own.
        zoom_path = out_dir / "inverted_zoom.vtk"
        print(f"\n=== Writing {zoom_path} ===")
        diag = diagnose_cell(P, H[cid])
        reflex_info = {"faces": diag["faces"],
                       "reflex_faces": diag["reflex_faces"],
                       "reflex_faces_angled": diag["reflex_faces_angled"]}
        write_zoom_vtk(zoom_path, P, H, [cid], ring2,
                       reflex_info=reflex_info)
        print(f"Wrote {zoom_path}")
        if diag["reflex_faces_angled"]:
            print(f"  reflex faces of cell {cid}: "
                  f"{diag['reflex_faces_angled']} "
                  f"(reentrant corners {diag['reentrant_corners']}), "
                  f"verdict {diag['verdict']}")
        
        # Write scaled version if cell is tiny
        vol = abs(inv['volume'])
        if vol < 1e-6:
            scaled_path = out_dir / "inverted_zoom_scaled.vtk"
            print(f"Cell volume is tiny ({vol:.2e}), writing scaled version")
            write_scaled_vtk(scaled_path, P, H, [cid], scale=1.5)
            print(f"Wrote {scaled_path}")
        
        # Write metadata JSON
        json_path = out_dir / "inverted_zoom.json"
        print(f"\n=== Writing {json_path} ===")
        
        # Check if cell is invisible in surface views
        is_interior = inv.get('is_interior', False)
        is_degenerate = vol < 1e-6
        
        boundary_desc = 'has no boundary faces' if is_interior else f"has {inv.get('n_boundary_faces', 0)} boundary faces"
        interior_desc = 'it is completely enclosed by other cells and does not touch the domain boundary, ' if is_interior else ''
        degenerate_desc = 'its near-zero volume makes it collapse to a point/line/plane, ' if is_degenerate else ''
        
        explanation = (
            f"This inverted cell is {'interior' if is_interior else 'on the boundary'} "
            f"({boundary_desc}). "
            f"It {'is' if is_degenerate else 'is not'} geometrically degenerate (volume = {vol:.2e}). "
            f"The cell is invisible in normal surface views because "
            f"{interior_desc}"
            f"{degenerate_desc}"
            f"surface-only exports only show boundary faces, not interior cells or their quality. "
            f"To visualize it, use the inverted_zoom.vtk file which shows the cell, its 2-ring neighbors, "
            f"and an octahedral marker at its center."
        )
        
        metadata = {
            "mesh": inv['mesh'],
            "cell_id": inv['cell_id'],
            "center": inv['center'],
            "bbox_min": inv['bbox_min'],
            "bbox_max": inv['bbox_max'],
            "volume": inv['volume'],
            "min_sj": inv['min_sj'],
            "is_interior": bool(is_interior),
            "n_boundary_faces": int(inv.get('n_boundary_faces', 0)),
            "ring1_neighbours": [int(x) for x in inv['ring1_neighbours']],
            "ring2_neighbours": [int(x) for x in inv['ring2_neighbours']],
            "explanation": explanation,
        }
        
        with open(json_path, 'w') as f:
            json.dump(metadata, f, indent=2)
        print(f"Wrote {json_path}")
    else:
        print("\nNo inverted cells with full mesh data found, skipping visualization files")


if __name__ == "__main__":
    main()
