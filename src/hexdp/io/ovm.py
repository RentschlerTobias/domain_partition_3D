"""OpenVolumeMesh (.ovm ASCII) reader + VTK writer for AlgoHex output.

AlgoHex writes its hex mesh in OVM's half-edge/half-face format, which is
not directly readable by meshio. Layout (see demo/HexMeshing/cylinder.ovm):

    OVM ASCII
    Vertices / <n> / "x y z" per line
    Edges    / <n> / "from_vertex to_vertex" per line
    Faces    / <n> / "<k> he0 he1 ... he(k-1)" per line
    Polyhedra/ <n> / "<k> hf0 hf1 ... hf(k-1)" per line

Index conventions (OpenVolumeMesh):
    halfedge 2*e   = edges[e] forward (from -> to)
    halfedge 2*e+1 = edges[e] reversed (to -> from)
    halfface 2*f   = faces[f] as listed
    halfface 2*f+1 = faces[f] with reversed orientation

Hex connectivity is recovered topologically rather than by trusting a
particular half-face ordering convention: take one face as "bottom", find
the disjoint opposite face as "top", then pair each bottom vertex with its
unique top neighbour along a cell edge. Final winding is fixed by the sign
of the cell volume, so the result is valid VTK regardless of how AlgoHex
ordered the half-faces.
"""

import numpy as np


def read_ovm(path):
    """Parse an OVM ASCII file. Returns (points, cells) where cells is a
    list of vertex-index lists (length 8 for hexes, 4 for tets)."""
    with open(path) as f:
        tokens = f.read().split("\n")

    i = 0

    def next_nonempty():
        nonlocal i
        while i < len(tokens) and not tokens[i].strip():
            i += 1
        val = tokens[i].strip()
        i += 1
        return val

    header = next_nonempty()
    if not header.startswith("OVM"):
        raise ValueError(f"not an OVM file: {header!r}")

    def expect_section(name):
        sec = next_nonempty()
        if sec.lower() != name.lower():
            raise ValueError(f"expected section {name!r}, got {sec!r}")
        return int(next_nonempty())

    n_verts = expect_section("Vertices")
    points = np.empty((n_verts, 3), float)
    for k in range(n_verts):
        points[k] = [float(x) for x in next_nonempty().split()]

    n_edges = expect_section("Edges")
    edges = np.empty((n_edges, 2), np.int64)
    for k in range(n_edges):
        edges[k] = [int(x) for x in next_nonempty().split()]

    n_faces = expect_section("Faces")
    faces = []
    for k in range(n_faces):
        parts = [int(x) for x in next_nonempty().split()]
        faces.append(parts[1:1 + parts[0]])

    n_polys = expect_section("Polyhedra")
    polys = []
    for k in range(n_polys):
        parts = [int(x) for x in next_nonempty().split()]
        polys.append(parts[1:1 + parts[0]])

    return points, edges, faces, polys


def halfedge_verts(he, edges):
    e, d = he // 2, he % 2
    a, b = int(edges[e][0]), int(edges[e][1])
    return (a, b) if d == 0 else (b, a)


def halfface_verts(hf, faces, edges):
    """Ordered vertex loop of a half-face."""
    f, d = hf // 2, hf % 2
    hes = faces[f]
    if d == 1:
        hes = [he ^ 1 for he in reversed(hes)]
    return [halfedge_verts(he, edges)[0] for he in hes]


def _cell_edges(hfs, faces, edges):
    """Undirected vertex-adjacency within one cell."""
    adj = {}
    for hf in hfs:
        loop = halfface_verts(hf, faces, edges)
        for k in range(len(loop)):
            a, b = loop[k], loop[(k + 1) % len(loop)]
            adj.setdefault(a, set()).add(b)
            adj.setdefault(b, set()).add(a)
    return adj


def _hex_volume(pts):
    """Signed volume of a trilinear hex via the standard 5-tet decomposition.

    NOTE: the variant in experimentell/3d_extrapolation/hexa_interpolation.py
    (hexa_cell_volumes) has two tets with their last two vertices swapped --
    (0,5,4,7) instead of (0,5,7,4) and (2,7,6,5) instead of (2,7,5,6) -- so
    those two contribute with the wrong sign and the total comes out at
    exactly 1/3 of the true volume. Verified: unit cube -> 0.3333, 2x3x4 box
    -> 8.0 instead of 24.0. That bug is pre-existing in this repo; it is
    harmless there because the value is only used for volume statistics
    (cell validity is judged by corner Jacobians), but it is load-bearing
    here, so this copy is corrected.
    """
    c = [np.asarray(pts[k], float) for k in range(8)]

    def tet(a, b, cc, d):
        return np.dot(b - a, np.cross(cc - a, d - a)) / 6.0

    return (tet(c[0], c[1], c[2], c[5]) + tet(c[0], c[2], c[3], c[7])
            + tet(c[0], c[5], c[7], c[4]) + tet(c[2], c[7], c[5], c[6])
            + tet(c[0], c[2], c[7], c[5]))


def hex_cell_vertices(hfs, faces, edges, points):
    """8 vertices of a hexahedron in VTK_HEXAHEDRON order."""
    if len(hfs) != 6:
        return None
    bottom = halfface_verts(hfs[0], faces, edges)
    if len(bottom) != 4:
        return None
    bset = set(bottom)
    top_hf = None
    for hf in hfs[1:]:
        loop = halfface_verts(hf, faces, edges)
        if len(loop) == 4 and not (set(loop) & bset):
            top_hf = hf
            break
    if top_hf is None:
        return None
    top = set(halfface_verts(top_hf, faces, edges))

    adj = _cell_edges(hfs, faces, edges)
    partner = []
    for v in bottom:
        up = [w for w in adj.get(v, ()) if w in top]
        if len(up) != 1:
            return None
        partner.append(up[0])

    cell = list(bottom) + partner
    if _hex_volume(points[cell]) < 0:
        # flip: swap bottom/top so the cell has positive orientation
        cell = partner + list(bottom)
    return cell


def ovm_to_cells(points, edges, faces, polys):
    """Convert OVM polyhedra to VTK cells. Returns (hexes, skipped)."""
    hexes, skipped = [], 0
    for hfs in polys:
        cell = hex_cell_vertices(hfs, faces, edges, points)
        if cell is None:
            skipped += 1
        else:
            hexes.append(cell)
    return np.asarray(hexes, np.int64) if hexes else np.zeros((0, 8), np.int64), skipped


def write_hex_vtk(path, points, hexes):
    """Legacy VTK unstructured grid of hexahedra (same writer style as
    dp3d.extraction._write_quad_vtk)."""
    with open(path, "w") as f:
        f.write("# vtk DataFile Version 2.0\n")
        f.write("AlgoHex hex mesh\n")
        f.write("ASCII\nDATASET UNSTRUCTURED_GRID\n")
        f.write(f"POINTS {len(points)} float\n")
        for p in points:
            f.write(f"{p[0]} {p[1]} {p[2]}\n")
        f.write(f"CELLS {len(hexes)} {9 * len(hexes)}\n")
        for c in hexes:
            f.write("8 " + " ".join(str(int(v)) for v in c) + "\n")
        f.write(f"CELL_TYPES {len(hexes)}\n")
        for _ in hexes:
            f.write("12\n")          # VTK_HEXAHEDRON


def write_hex_msh(path, points, hexes, cell_data=None, lines=None,
                  line_tags=None):
    """Gmsh 2.2 ASCII hexahedral mesh (element type 5). Optional per-cell
    integer tag is written as the physical tag so it survives the round trip
    into gmsh/ParaView.

    ``lines`` adds 2-node line elements (type 1) alongside the hexes, with
    ``line_tags`` as their physical tag. That is how the block edges are
    carried: a solid hex mesh shows nothing of the block structure, and the
    edge wireframe is what makes it visible (and is what the TFI division
    MILP will need to tag anyway)."""
    n_line = 0 if lines is None else len(lines)
    with open(path, "w") as f:
        f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
        f.write(f"$Nodes\n{len(points)}\n")
        for i, p in enumerate(points, 1):
            f.write(f"{i} {p[0]} {p[1]} {p[2]}\n")
        f.write("$EndNodes\n")
        f.write(f"$Elements\n{n_line + len(hexes)}\n")
        eid = 0
        for k in range(n_line):
            eid += 1
            tag = int(line_tags[k]) if line_tags is not None else 1
            a, b = lines[k]
            f.write(f"{eid} 1 2 {tag} 1 {int(a) + 1} {int(b) + 1}\n")
        for c in hexes:
            eid += 1
            phys = int(cell_data[eid - n_line - 1]) if cell_data is not None else 1
            nodes = " ".join(str(int(v) + 1) for v in c)
            f.write(f"{eid} 5 2 {phys} 1 {nodes}\n")
        f.write("$EndElements\n")


def scaled_jacobian(points, cell):
    """Minimum corner scaled Jacobian of one hex -- the correct validity
    test (a sheared hex can have positive signed volume yet inverted
    corners, cf. hexa_corner_jacobians in hexa_interpolation.py)."""
    P8 = points[cell]
    idx = [(0, 1, 3, 4), (1, 2, 0, 5), (2, 3, 1, 6), (3, 0, 2, 7),
           (4, 7, 5, 0), (5, 4, 6, 1), (6, 5, 7, 2), (7, 6, 4, 3)]
    out = []
    for a, b, c, d in idx:
        u, v, w = P8[b] - P8[a], P8[c] - P8[a], P8[d] - P8[a]
        n = np.linalg.norm(u) * np.linalg.norm(v) * np.linalg.norm(w)
        out.append(float(np.dot(u, np.cross(v, w)) / n) if n > 0 else 0.0)
    return min(out)


def convert(ovm_path, vtk_path):
    points, edges, faces, polys = read_ovm(ovm_path)
    print(f"[ovm_io] {ovm_path}: {len(points)} verts, {len(edges)} edges, "
          f"{len(faces)} faces, {len(polys)} polyhedra")
    hexes, skipped = ovm_to_cells(points, edges, faces, polys)
    print(f"[ovm_io] recovered {len(hexes)} hexahedra "
          f"({skipped} non-hex/failed cells skipped)")
    write_hex_vtk(vtk_path, points, hexes)
    print(f"[ovm_io] wrote {vtk_path}")
    return points, hexes


if __name__ == "__main__":
    import sys
    convert(sys.argv[1], sys.argv[2])
