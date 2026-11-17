"""Gmsh .msh (4.1) export of a block structure: the fine hex mesh with physical groups, e.g. for FEniCSx
(``dolfinx.io.gmshio``).

Physical groups: one volume group ("domain", tag 1) or one per block, and one surface group per geometry
surface of the boundary faces (from ``face_info``: name "<role>_<surface tag>", tag = surface tag), faces
without that information in "boundary" (tag 0 is not used by gmsh, so it gets max tag + 1). Elements are
linear (hex8, quad4); the element order of a FEniCS function space is independent of the mesh.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from hexdp.core.blockstructure import OFF, grid_cells


def _grids(bs):
    """Fine node coordinates and per-block node-id grids (the coarse hexes when there are none)."""
    if bs.grids:
        return np.asarray(bs.fine_points, float), [np.asarray(G) for G in bs.grids]
    grids = [np.array(c8)[[np.flatnonzero((OFF == o).all(1))[0] for o in np.ndindex(2, 2, 2)]].reshape(2, 2, 2)
             for c8 in bs.blocks]
    return np.asarray(bs.points, float), grids


def boundary_groups(bs) -> dict[tuple[str, int], list[tuple[int, int, int]]]:
    """(name, tag) -> boundary faces (block, axis, side), grouped by the geometry surface of face_info."""
    groups: dict[tuple[str, int], list] = {}
    rest = []
    for key, users in bs.faces().items():
        if len(users) != 1:
            continue
        info = bs.face_info.get(key)
        if info is None or info.surface is None:
            rest.append(users[0])
        else:
            groups.setdefault((f"{info.role.value}_{info.surface}", int(info.surface)), []).append(users[0])
    if rest:
        groups[("boundary", max([t for _, t in groups] + [0]) + 1)] = rest
    return groups


def write_msh(bs, path, block_groups: bool = False) -> Path:
    import gmsh

    P, grids = _grids(bs)
    path = Path(path).with_suffix(".msh")
    own = not gmsh.isInitialized()
    if own:
        gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("hexdp")
        used = np.unique(np.concatenate([G.ravel() for G in grids]))
        vol_tags = list(range(1, len(grids) + 1)) if block_groups else [1]
        for t in vol_tags:
            gmsh.model.addDiscreteEntity(3, t)
        gmsh.model.mesh.addNodes(3, vol_tags[0], (used + 1).tolist(), P[used].ravel().tolist())
        etag = 1
        cells_of = [grid_cells(G) for G in grids]
        if block_groups:
            for b, cells in enumerate(cells_of):
                gmsh.model.mesh.addElementsByType(b + 1, 5, list(range(etag, etag + len(cells))), (cells + 1).ravel().tolist())
                etag += len(cells)
        else:
            cells = np.vstack(cells_of)
            gmsh.model.mesh.addElementsByType(1, 5, list(range(etag, etag + len(cells))), (cells + 1).ravel().tolist())
            etag += len(cells)
        for t in vol_tags:
            gmsh.model.addPhysicalGroup(3, [t], t, name="domain" if not block_groups else f"block_{t - 1}")
        for (name, tag), faces in sorted(boundary_groups(bs).items(), key=lambda kv: kv[0][1]):
            quads = []
            for b, axis, side in faces:
                G = grids[b]
                idx = [slice(None)] * 3
                idx[axis] = -1 if side else 0
                F = G[tuple(idx)]
                quads.append(np.stack([F[:-1, :-1], F[1:, :-1], F[1:, 1:], F[:-1, 1:]], -1).reshape(-1, 4))
            quads = np.vstack(quads)
            gmsh.model.addDiscreteEntity(2, tag)
            gmsh.model.mesh.addElementsByType(tag, 3, list(range(etag, etag + len(quads))), (quads + 1).ravel().tolist())
            etag += len(quads)
            gmsh.model.addPhysicalGroup(2, [tag], tag, name=name)
        gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
        gmsh.option.setNumber("Mesh.SaveAll", 0)
        gmsh.write(str(path))
        gmsh.model.remove()
    finally:
        if own:
            gmsh.finalize()
    return path

