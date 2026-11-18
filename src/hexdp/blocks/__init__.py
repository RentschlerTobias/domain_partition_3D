"""Building blocks of a composite blocking: steps that produce intermediate results (inputs, surface layouts) or
parts of the block structure, combined by the ``composite`` strategy from a workflow's ``steps``.

A step is a plugin of kind ``blocks`` (entry-point group hexdp.blocks) with ``params`` like every plugin and

    run(inputs: dict[str, object], params: dict, ctx) -> object

``inputs`` are the results of earlier steps named in the step's ``from``. A result that is a ``Part`` becomes part of
the final block structure; ``compose`` merges the parts at coincident fine nodes (conforming interfaces) and labels
the boundary faces with the names the parts give them.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class Part:
    """Blocks of one sub-domain: fine node positions, per block a structured node-id grid (local ids), and a
    classifier naming a boundary face from its fine node positions (None = no name; a face that becomes interior
    after composing is not named)."""
    name: str
    points: np.ndarray
    grids: list
    label: object = None                    # callable(F (n, 3)) -> str | None
    info: dict = field(default_factory=dict)


def compose(parts, tol: float = 1e-9, groups: dict | None = None):
    """One BlockStructure from the parts: fine nodes closer than tol x size are merged (the interfaces of the parts
    must have coincident nodes - conforming by construction), every block right-handed, boundary faces labelled
    by the part they belong to. Returns (bs, info)."""
    from scipy.spatial import cKDTree

    from hexdp.core.blockstructure import FaceInfo
    from hexdp.core.geometry import PatchRole
    from hexdp.io.blocks_vtk import from_grids

    P = np.vstack([p.points for p in parts])
    off = np.cumsum([0] + [len(p.points) for p in parts])
    scale = max(float(np.ptp(P, axis=0).max()), 1e-300)
    pairs = cKDTree(P).query_pairs(tol * scale, output_type="ndarray")
    root = np.arange(len(P))

    def find(x):
        while root[x] != x:
            root[x] = root[root[x]]
            x = root[x]
        return x
    for a, b in pairs:
        ra, rb = find(a), find(b)
        if ra != rb:
            root[max(ra, rb)] = min(ra, rb)
    root = np.array([find(i) for i in range(len(P))])
    uniq, new = np.unique(root, return_inverse=True)
    new = new.ravel()
    Pm = P[uniq]
    grids, owner = [], []
    for k, part in enumerate(parts):
        for G in part.grids:
            G = new[np.asarray(G) + off[k]]
            X = Pm[G]
            d = np.dot(np.cross(X[1, 0, 0] - X[0, 0, 0], X[0, 1, 0] - X[0, 0, 0]), X[0, 0, 1] - X[0, 0, 0])
            grids.append(G if d > 0 else G[:, :, ::-1])
            owner.append(k)
    bs = from_grids(Pm, grids)
    groups = dict(groups or {})
    for key, users in bs.faces().items():
        if len(users) != 1:
            continue
        b, ax, sd = users[0]
        lab = parts[owner[b]].label
        name = lab(bs._face_points(b, ax, sd)) if lab else None
        if name is None:
            continue
        tag = groups.setdefault(name, len(groups) + 1)
        bs.face_info[key] = FaceInfo(PatchRole.WALL, tag)
    bs.meta["surface_names"] = {t: n for n, t in groups.items()}
    info = {"parts": {p.name: len(p.grids) for p in parts}, "merged_nodes": int(len(P) - len(Pm))}
    return bs, info
