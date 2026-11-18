"""AlgoHex-format tet meshes (legacy VTK: vertex/line/triangle/tetra cells over one point list, an int id per cell;
nonzero ids on vertices/lines/triangles = features). HexMe ships them; tet_prep writes them.

    m = read(path)                          # TetVTK(points, tets, tet_ids, tris, tri_ids, lines, line_ids, verts, vert_ids)
    m.volumes()                             # {volume id: tet count}
    volume_roles(m)                         # {"outer": [...], "embedded": [...]} (box vs. bodies inside it)
    sub = extract_volume(m, 1)              # one volume with its own boundary and the features on it
    write(sub, path)

extract_volume turns a body embedded in a meshed box (HexMe 'b' models: volume 1 = box minus body, 2 = body) into
the domain around the body (an external-flow case) or the body alone. Faces the volume shares with the removed
volumes become its boundary; boundary faces without a triangle id get `new_face_id`.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

_TYPES = {"vertex": 1, "line": 3, "triangle": 5, "tetra": 10}


@dataclass
class TetVTK:
    points: np.ndarray
    tets: np.ndarray
    tet_ids: np.ndarray
    tris: np.ndarray
    tri_ids: np.ndarray
    lines: np.ndarray
    line_ids: np.ndarray
    verts: np.ndarray
    vert_ids: np.ndarray

    def volumes(self) -> dict[int, int]:
        u, n = np.unique(self.tet_ids, return_counts=True)
        return {int(a): int(b) for a, b in zip(u, n)}


def read(path) -> TetVTK:
    import meshio
    m = meshio.read(str(path))
    key = next((k for k in m.cell_data if k.lower() in ("cellentityids", "color", "gmsh:physical")), None)
    got = {t: (np.zeros((0, k), int), np.zeros(0, int)) for t, k in (("vertex", 1), ("line", 2), ("triangle", 3),
                                                                        ("tetra", 4))}
    for i, c in enumerate(m.cells):
        if c.type in got:
            ids = np.asarray(m.cell_data[key][i], int).ravel() if key else np.zeros(len(c.data), int)
            a, b = got[c.type]
            got[c.type] = (np.r_[a, np.asarray(c.data, int).reshape(len(c.data), -1)], np.r_[b, ids])
    return TetVTK(np.asarray(m.points, float), got["tetra"][0], got["tetra"][1], got["triangle"][0], got["triangle"][1],
                  got["line"][0], got["line"][1], got["vertex"][0][:, 0], got["vertex"][1])


def _boundary_faces(tets):
    F = np.r_[tets[:, [1, 2, 3]], tets[:, [0, 3, 2]], tets[:, [0, 1, 3]], tets[:, [0, 2, 1]]]   # outward if tet positive
    key = np.sort(F, 1)
    _, inv, cnt = np.unique(key, axis=0, return_inverse=True, return_counts=True)
    return F[cnt[inv.ravel()] == 1]


def volume_roles(m: TetVTK) -> dict[str, list[int]]:
    """{'outer': volume ids touching the outer boundary of the whole mesh (e.g. the box around a body),
    'embedded': the others (bodies inside it)} - independent of how the ids are numbered."""
    outer = set(np.unique(_boundary_faces(m.tets).ravel()).tolist())
    roles = {"outer": [], "embedded": []}
    for v in m.volumes():
        pts = np.unique(m.tets[m.tet_ids == v])
        roles["outer" if outer.intersection(pts.tolist()) else "embedded"].append(v)
    return roles


def extract_volume(m: TetVTK, volume_ids, new_face_id: int | None = None) -> TetVTK:
    vols = {int(volume_ids)} if np.isscalar(volume_ids) else {int(v) for v in volume_ids}
    keep = np.isin(m.tet_ids, list(vols))
    if not keep.any():
        raise ValueError(f"no tets with volume id {sorted(vols)} (volumes: {m.volumes()})")
    T = m.tets[keep]
    B = _boundary_faces(T)
    # boundary faces carry the id of the matching input triangle, else new_face_id
    tri_key = {tuple(sorted(t)): int(i) for t, i in zip(m.tris, m.tri_ids)}
    nid = (max(tri_key.values(), default=0) + 1) if new_face_id is None else int(new_face_id)
    B_ids = np.array([tri_key.get(tuple(sorted(f)), nid) for f in B], int)
    on_bnd = np.zeros(len(m.points), bool)
    on_bnd[B.ravel()] = True
    bkeys = {tuple(sorted(e)) for f in B for e in ((f[0], f[1]), (f[1], f[2]), (f[0], f[2]))}
    lk = np.array([tuple(sorted(e)) in bkeys for e in m.lines], bool) if len(m.lines) else np.zeros(0, bool)
    vk = on_bnd[m.verts] if len(m.verts) else np.zeros(0, bool)
    used = np.unique(T)
    remap = -np.ones(len(m.points), int)
    remap[used] = np.arange(len(used))
    return TetVTK(m.points[used], remap[T], m.tet_ids[keep], remap[B], B_ids, remap[m.lines[lk]], m.line_ids[lk],
                  remap[m.verts[vk]], m.vert_ids[vk])


def write(m: TetVTK, path, title: str = "hexdp tet mesh (AlgoHex format)"):
    cells = [("vertex", m.verts[:, None], m.vert_ids), ("line", m.lines, m.line_ids),
             ("triangle", m.tris, m.tri_ids), ("tetra", m.tets, m.tet_ids)]
    n = sum(len(c) for _, c, _ in cells)
    with open(path, "w") as f:
        f.write(f"# vtk DataFile Version 2.0\n{title}\nASCII\nDATASET UNSTRUCTURED_GRID\n")
        f.write(f"POINTS {len(m.points)} double\n")
        np.savetxt(f, m.points, fmt="%.17g")
        f.write(f"CELLS {n} {sum(len(c) * (c.shape[1] + 1) for _, c, _ in cells if len(c))}\n")
        for _, c, _ in cells:
            if len(c):
                np.savetxt(f, np.c_[np.full(len(c), c.shape[1]), c], fmt="%d")
        f.write(f"CELL_TYPES {n}\n")
        for t, c, _ in cells:
            if len(c):
                np.savetxt(f, np.full(len(c), _TYPES[t]), fmt="%d")
        f.write(f"CELL_DATA {n}\nSCALARS CellEntityIds int 1\nLOOKUP_TABLE default\n")
        for _, c, ids in cells:
            if len(c):
                np.savetxt(f, ids, fmt="%d")
    return Path(path)
