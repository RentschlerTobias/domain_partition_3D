"""Load any boundary description into a ``Geometry`` through gmsh.

    .msh                      entities + physical groups as stored (dtOO, gmsh, HexMe)
    .step .stp .iges .igs .brep  CAD through gmsh's OpenCASCADE kernel; triangulated for queries,
                              exact closest points through ``closest_on_cad_faces``
    .stl                      no entities: surfaces are segmented by gmsh classifySurfaces (angle) and
                              the result is marked in ``Geometry.notes``
    .vtk                      a tet mesh with an integer cell field per entity (AlgoHex input format,
                              e.g. HexMe 'CellEntityIds' or tet_prep 'color'): converted to .msh with the
                              field as entity and physical tags

gmsh is a process-global library; every function here opens and closes its own session.
"""

from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path

import numpy as np

from hexdp.core.geometry import (
    Curve,
    Geometry,
    PatchRole,
    RoleSource,
    Source,
    Surface,
    classify_curves,
    role_from_name,
)

CAD_SUFFIXES = {".step", ".stp", ".iges", ".igs", ".brep"}


@contextmanager
def gmsh_session():
    import gmsh

    own = not gmsh.isInitialized()
    if own:
        # no SIGINT handler: gmsh may run in a worker thread (MCP server), where signals cannot be set
        gmsh.initialize(interruptible=False)
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.option.setNumber("General.Verbosity", 1)
    try:
        gmsh.clear()
        yield gmsh
    finally:
        gmsh.clear()
        if own:
            gmsh.finalize()


def load(path, mesh_size: float | None = None, stl_angle_deg: float = 40.0, sharp_deg: float = 30.0) -> Geometry:
    """Read ``path`` into a Geometry. ``mesh_size`` sets the triangulation size for CAD input
    (default: 1/50 of the bounding-box diagonal)."""
    path = Path(path)
    suf = path.suffix.lower()
    if suf == ".json":                       # parameters of a parametric design: no surfaces until a strategy builds it
        return Geometry(Source.PARAMS, str(path), {}, {}, {},
                        notes=["parameter file: the geometry is built by a strategy (e.g. 'runner' from dtOO)"])
    with gmsh_session() as gmsh:
        if suf == ".msh":
            gmsh.open(str(path))
            source, notes = Source.MSH, []
            _ensure_boundary_mesh(gmsh, mesh_size)
        elif suf in CAD_SUFFIXES:
            gmsh.model.occ.importShapes(str(path))
            gmsh.model.occ.synchronize()
            source, notes = Source.CAD, []
            _mesh_surfaces(gmsh, mesh_size)
        elif suf == ".vtk":
            tmp = _vtk_to_msh(path)
            gmsh.open(str(tmp))
            source, notes = Source.MSH, [f"VTK tet mesh: entities from cell field '{tmp.stem.split('__')[-1]}'"]
        elif suf == ".stl":
            gmsh.merge(str(path))
            gmsh.model.mesh.classifySurfaces(np.radians(stl_angle_deg), True, False, np.pi)
            gmsh.model.mesh.createGeometry()
            source = Source.STL
            notes = [f"STL has no entities: surfaces segmented at {stl_angle_deg:g} deg, features inferred"]
        else:
            raise ValueError(f"unsupported geometry file '{path.name}'")
        geo = _collect(gmsh, source, str(path), notes)
    classify_curves(geo, sharp_deg=sharp_deg)
    return geo


def _vtk_to_msh(path):
    """Write a gmsh 2.2 copy of a VTK tet mesh whose integer cell field tags the entities."""
    import tempfile

    import meshio

    m = meshio.read(path)
    field = next((k for k, v in m.cell_data.items() if np.asarray(v[0]).dtype.kind in "iu"), None)
    if field is None:
        raise ValueError(f"{path}: no integer cell field to take the entities from")
    keep = [(c, np.asarray(d).ravel().astype(int)) for c, d in zip(m.cells, m.cell_data[field])
            if c.type in ("vertex", "line", "triangle", "tetra")]
    out = meshio.Mesh(m.points, [(c.type, c.data) for c, _ in keep],
                      cell_data={"gmsh:physical": [d for _, d in keep], "gmsh:geometrical": [d for _, d in keep]})
    tmp = Path(tempfile.mkdtemp(prefix="hexdp_vtk_")) / f"{Path(path).stem}__{field}.msh"
    meshio.write(tmp, out, file_format="gmsh22", binary=False)
    return tmp


def _mesh_surfaces(gmsh, mesh_size):
    lo, hi = _bbox(gmsh)
    h = mesh_size or float(np.linalg.norm(hi - lo)) / 50
    gmsh.option.setNumber("Mesh.MeshSizeMax", h)
    gmsh.option.setNumber("Mesh.MeshSizeMin", h / 20)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 12)
    gmsh.model.mesh.generate(2)


def _ensure_boundary_mesh(gmsh, mesh_size):
    """A .msh normally carries the surface triangles; mesh the surfaces only if they are empty
    (a geometry-only .msh)."""
    for dim, tag in gmsh.model.getEntities(2):
        _, elem_tags, _ = gmsh.model.mesh.getElements(dim, tag)
        if any(len(t) for t in elem_tags):
            return
    _mesh_surfaces(gmsh, mesh_size)


def _bbox(gmsh):
    b = np.array(gmsh.model.getBoundingBox(-1, -1))
    return b[:3], b[3:]


def _names(gmsh) -> dict[tuple[int, int], str]:
    """Entity -> name: physical group names first, then the entity's own name (CAD)."""
    out = {}
    for dim, ptag in gmsh.model.getPhysicalGroups():
        name = gmsh.model.getPhysicalName(dim, ptag)
        if not name:
            continue
        for t in gmsh.model.getEntitiesForPhysicalGroup(dim, ptag):
            out.setdefault((dim, int(t)), name)
    for dim in (1, 2):
        for _, t in gmsh.model.getEntities(dim):
            name = gmsh.model.getEntityName(dim, t)
            if name:
                out.setdefault((dim, int(t)), name.split("/")[-1])
    return out


def _surface_mesh(gmsh, tag):
    """Triangles of one surface (quads split), as local node array + index array."""
    types, _, enodes = gmsh.model.mesh.getElements(2, tag)
    tris = []
    for et, en in zip(types, enodes):
        _, _, _, nn, _, _ = gmsh.model.mesh.getElementProperties(et)
        corners = 3 if nn in (3, 6) else 4 if nn in (4, 8, 9) else None
        if corners is None:
            continue
        E = np.asarray(en, np.int64).reshape(-1, nn)[:, :corners]
        if corners == 4:
            E = np.vstack([E[:, [0, 1, 2]], E[:, [0, 2, 3]]])
        tris.append(E)
    if not tris:
        return np.zeros((0, 3)), np.zeros((0, 3), np.int64), np.zeros(0, np.int64)
    T = np.vstack(tris)
    used, inv = np.unique(T, return_inverse=True)
    coords = _node_coords(gmsh, used)
    return coords, inv.reshape(T.shape), used


def _node_coords(gmsh, node_tags):
    tags, coord, _ = gmsh.model.mesh.getNodes(-1, -1, includeBoundary=True)
    lookup = dict(zip(np.asarray(tags, np.int64).tolist(), range(len(tags))))
    C = np.asarray(coord).reshape(-1, 3)
    return C[[lookup[int(t)] for t in node_tags]]


def _curve_points(gmsh, tag):
    """Curve polyline in curve order: the chain of the curve's line elements (works for any .msh, also
    without stored parametric coordinates), else a parametric sampling when the curve has no mesh."""
    types, _, enodes = gmsh.model.mesh.getElements(1, tag)
    segs = [np.asarray(en, np.int64).reshape(-1, 2 if t == 1 else 3)[:, :2] for t, en in zip(types, enodes) if t in (1, 8)]
    if segs:
        order = _chain(np.vstack(segs))
        if order is not None:
            return _node_coords(gmsh, order)
    try:
        lo, hi = gmsh.model.getParametrizationBounds(1, tag)
        u = np.linspace(lo[0], hi[0], 33)
        return np.asarray(gmsh.model.getValue(1, tag, u)).reshape(-1, 3)
    except Exception:
        _, coord, _ = gmsh.model.mesh.getNodes(1, tag, includeBoundary=True)
        return np.asarray(coord).reshape(-1, 3)


def _chain(S):
    """Order the segments (k, 2) of one open or closed polyline into a node sequence; None if branched."""
    nbr: dict[int, list[int]] = {}
    for a, b in S.tolist():
        nbr.setdefault(a, []).append(b)
        nbr.setdefault(b, []).append(a)
    if any(len(v) > 2 for v in nbr.values()):
        return None
    ends = [n for n, v in nbr.items() if len(v) == 1]
    start = ends[0] if ends else next(iter(nbr))
    seq, prev = [start], None
    while True:
        nxt = [n for n in nbr[seq[-1]] if n != prev]
        if not nxt or (nxt[0] == start):
            if nxt and nxt[0] == start:
                seq.append(start)      # closed curve: repeat the start node
            break
        prev = seq[-1]
        seq.append(nxt[0])
        if len(seq) > len(nbr) + 1:
            return None
    return np.asarray(seq, np.int64)


def _collect(gmsh, source, path, notes) -> Geometry:
    names = _names(gmsh)
    surfaces, surface_nodes = {}, {}
    for _, t in gmsh.model.getEntities(2):
        nodes, tris, tags = _surface_mesh(gmsh, t)
        surface_nodes[int(t)] = set(tags.tolist())
        name = names.get((2, t))
        role = role_from_name(name)
        surfaces[int(t)] = Surface(int(t), nodes, tris, name, role,
                                   RoleSource.NAME if role is not PatchRole.UNKNOWN else RoleSource.NONE)
    curves, derived = {}, 0
    for _, t in gmsh.model.getEntities(1):
        up, _ = gmsh.model.getAdjacencies(1, t)
        adj = tuple(int(x) for x in up if int(x) in surfaces)
        if not adj:
            adj = _surfaces_by_nodes(gmsh, t, surface_nodes)
            derived += bool(adj)
        curves[int(t)] = Curve(int(t), _curve_points(gmsh, t), adj, name=names.get((1, t)))
    if derived:
        notes.append(f"{derived} curve-surface adjacencies derived from shared mesh nodes (no entity topology stored)")
    points, point_curves = {}, {}
    for _, t in gmsh.model.getEntities(0):
        points[int(t)] = np.asarray(gmsh.model.getValue(0, t, []), float)
        up, _ = gmsh.model.getAdjacencies(0, t)
        point_curves[int(t)] = tuple(int(x) for x in up)
    surface_volumes = {}
    for _, t in gmsh.model.getEntities(2):
        up, _ = gmsh.model.getAdjacencies(2, t)
        surface_volumes[int(t)] = tuple(int(x) for x in up)
    empty = [t for t, s in surfaces.items() if not len(s.tris)]
    if empty:
        notes.append(f"{len(empty)} surfaces without triangles: {empty[:10]}")
        for t in empty:
            surfaces.pop(t)
    return Geometry(source, path, surfaces, curves, points, point_curves, surface_volumes,
                    exact=source is Source.CAD, notes=notes)


def _surfaces_by_nodes(gmsh, tag, surface_nodes):
    """Surfaces whose triangles use (nearly) all nodes of the curve's line elements."""
    types, _, enodes = gmsh.model.mesh.getElements(1, tag)
    cn = set(np.concatenate([np.asarray(e, np.int64) for e in enodes]).tolist()) if len(enodes) else set()
    if not cn:
        return ()
    return tuple(s for s, nodes in surface_nodes.items() if len(cn & nodes) >= 0.9 * len(cn))


def closest_on_cad_faces(path, Q, faces, fallback):
    """Exact closest points on the trimmed CAD faces ``faces`` (one face tag per point).

    gmsh's getClosestPoint projects onto the untrimmed carrier surface; a result outside the face's
    trimming is replaced by ``fallback`` (the triangulation projection) for that point. Returns
    (points, distances, exact_mask)."""
    Q = np.asarray(Q, float)
    R = np.array(fallback, float).reshape(-1, 3)
    exact = np.zeros(len(Q), bool)
    with gmsh_session() as gmsh:
        gmsh.model.occ.importShapes(str(path))
        gmsh.model.occ.synchronize()
        for f in np.unique(faces):
            sel = np.where(faces == f)[0]
            xyz, uv = gmsh.model.getClosestPoint(2, int(f), Q[sel].ravel().tolist())
            xyz = np.asarray(xyz).reshape(-1, 3)
            uv = np.asarray(uv).reshape(-1, 2)
            # isInside returns the number of the given points inside the entity: ask point by point
            inside = np.array([gmsh.model.isInside(2, int(f), p.tolist(), parametric=True) > 0 for p in uv], bool)
            R[sel[inside]] = xyz[inside]
            exact[sel[inside]] = True
    return R, np.linalg.norm(R - Q, axis=1), exact
