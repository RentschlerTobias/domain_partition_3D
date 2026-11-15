"""Block structures from and to hex VTK files with a per-cell block id.

``read_blocks_vtk`` is the reader of the existing pipelines (beam_collapse, the canada block collapse):
the fine hex mesh plus a cell field (default ``block``) naming the block of every cell. The (i, j, k)
indexing of each block is recovered by walking face neighbours, so nothing is assumed about the writer.
"""

from __future__ import annotations

import collections

import numpy as np

from hexdp.core.blockstructure import HEX_EDGES, OFF, PARALLEL, BlockStructure, grid_cells


def _hex_signed(X):
    """sign of the volume of the fine hexes of a structured grid (n1+1, n2+1, n3+1, 3), median."""
    a = X[1:, :-1, :-1] - X[:-1, :-1, :-1]
    b = X[:-1, 1:, :-1] - X[:-1, :-1, :-1]
    c = X[:-1, :-1, 1:] - X[:-1, :-1, :-1]
    return float(np.median(np.einsum("...i,...i->...", np.cross(a, b), c)))


def structured_grids(P, H, bid):
    """Per block id the structured node-id grid G (n1+1, n2+1, n3+1), right-handed."""
    nbr_local: dict[int, list[int]] = {}
    for a, b in HEX_EDGES:
        nbr_local.setdefault(a, []).append(b)
        nbr_local.setdefault(b, []).append(a)
    face_loc = {(ax, s): [n for n in range(8) if OFF[n][ax] == s] for ax in range(3) for s in (0, 1)}
    grids, ids = [], []
    for b in np.unique(bid):
        cells = np.where(bid == b)[0]
        f2c = collections.defaultdict(list)
        for c in cells:
            for loc in face_loc.values():
                f2c[frozenset(H[c, loc].tolist())].append(c)
        c0 = cells[0]
        ijk = {int(H[c0, n]): OFF[n].copy() for n in range(8)}
        seen, dq = {c0}, collections.deque([c0])
        while dq:
            c = dq.popleft()
            hc = H[c]
            lc = {int(v): i for i, v in enumerate(hc)}
            for loc in face_loc.values():
                key = frozenset(hc[loc].tolist())
                for g in f2c[key]:
                    if g in seen:
                        continue
                    hg = H[g]
                    lg = {int(v): i for i, v in enumerate(hg)}
                    for n in key:
                        m_ = [int(hc[x]) for x in nbr_local[lc[n]] if int(hc[x]) not in key][0]
                        p_ = [int(hg[x]) for x in nbr_local[lg[n]] if int(hg[x]) not in key][0]
                        ijk[p_] = 2 * ijk[n] - ijk[m_]
                    seen.add(g)
                    dq.append(g)
        if len(seen) != len(cells):
            raise ValueError(f"block {b}: cells are not face-connected")
        K = np.array(list(ijk.values()))
        lo = K.min(0)
        dims = K.max(0) - lo + 1
        G = -np.ones(dims, np.int64)
        for n, v in ijk.items():
            G[tuple(v - lo)] = n
        if not ((G >= 0).all() and G.size == len(ijk) and np.prod(dims - 1) == len(cells)):
            raise ValueError(f"block {b}: cells do not form a structured grid")
        if _hex_signed(P[G]) < 0:
            G = G[:, :, ::-1]
        grids.append(G)
        ids.append(int(b))
    return grids, ids


BLOCK_FIELDS = ("block", "block_id", "blocks", "BlockId")


def merge_coincident(P, H):
    """Some writers store the nodes on block interfaces once per block. Coincident points (exactly equal
    coordinates) become one node, so neighbouring blocks share their faces again. Returns (P, H, merged)."""
    _, first, inv = np.unique(P, axis=0, return_index=True, return_inverse=True)
    if len(first) == len(P):
        return P, H, 0
    keep = np.sort(first)                       # keep the first occurrence, in file order
    new_id = np.empty(len(first), np.int64)
    new_id[np.argsort(first)] = np.arange(len(first))
    old_to_new = new_id[inv.reshape(-1)]
    return P[keep], old_to_new[H], len(P) - len(keep)


def read_blocks_vtk(path, field: str | None = None) -> BlockStructure:
    """Hex VTK with a per-cell block id -> BlockStructure with coarse topology, cell counts and fine grids.
    ``field`` defaults to the first of BLOCK_FIELDS present."""
    import meshio

    m = meshio.read(path)
    if field is None:
        field = next((f for f in BLOCK_FIELDS if f in m.cell_data), BLOCK_FIELDS[0])
    if field not in m.cell_data:
        raise ValueError(f"{path}: no cell field '{field}' (have {sorted(m.cell_data)})")
    H = np.vstack([c.data for c in m.cells if c.type == "hexahedron"]).astype(np.int64)
    bid = np.concatenate([d for d, c in zip(m.cell_data[field], m.cells) if c.type == "hexahedron"]).astype(int)
    P = np.asarray(m.points, float)
    P, H, merged = merge_coincident(P, H)
    grids, _ = structured_grids(P, H, bid)
    bs = from_grids(P, grids)
    bs.meta["source"] = str(path)
    if merged:
        bs.meta["merged_coincident_points"] = merged
    bs.meta["block_field"] = field
    return bs


def from_grids(P, grids) -> BlockStructure:
    """Coarse structure from structured fine grids: block corners become coarse vertices (shared by
    fine node id), edge cell counts come from the grid dimensions."""
    corner_ids, blocks = {}, []
    for G in grids:
        n = np.array(G.shape) - 1
        c8 = []
        for o in OFF:
            fid = int(G[o[0] * n[0], o[1] * n[1], o[2] * n[2]])
            c8.append(corner_ids.setdefault(fid, len(corner_ids)))
        blocks.append(c8)
    fine_of = np.empty(len(corner_ids), np.int64)
    for fid, cid in corner_ids.items():
        fine_of[cid] = fid
    bs = BlockStructure(points=P[fine_of], blocks=np.array(blocks), fine_points=P, grids=list(grids))
    for b, G in enumerate(grids):
        n = np.array(G.shape) - 1
        for ax, pairs in enumerate(PARALLEL):
            for i, j in pairs:
                bs.edge_cells[frozenset((int(bs.blocks[b][i]), int(bs.blocks[b][j])))] = int(n[ax])
    return bs


def write_blocks_vtk(bs: BlockStructure, path, field: str = "block") -> None:
    """Fine hex mesh with a per-cell block id (needs ``grids``), or the coarse hexes when there are none."""
    import meshio

    if not bs.grids:
        meshio.write(path, meshio.Mesh(bs.points, [("hexahedron", bs.blocks)],
                                       cell_data={field: [np.arange(bs.n_blocks)]}))
        return
    cells, ids = [], []
    for b, G in enumerate(bs.grids):
        c = grid_cells(G)
        cells.append(c)
        ids.append(np.full(len(c), b))
    meshio.write(path, meshio.Mesh(bs.fine_points, [("hexahedron", np.vstack(cells))],
                                   cell_data={field: [np.concatenate(ids)]}))
