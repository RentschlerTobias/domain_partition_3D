"""BlockStructure: a conforming coarse hex mesh, one coarse hex per block, plus per-edge and per-face data.

Corner order is the VTK hexahedron order, (i, j, k) offsets:

    c0=(0,0,0) c1=(1,0,0) c2=(1,1,0) c3=(0,1,0) c4=(0,0,1) c5=(1,0,1) c6=(1,1,1) c7=(0,1,1)

Edges and faces are keyed by their vertex sets (``frozenset``), so two blocks that share an edge or a
face share the key. On top of the topology a structure carries, all optional:

    edge_cells    cell count per edge (all edges of one sheet must agree)
    edge_grading  grading per edge (ratio last/first cell, along the edge's vertex order)
    edge_curves   polyline geometry of an edge (curved block edges)
    face_roles    patch role and geometry surface of a boundary face
    grids         per block, the structured fine-node grid (n1+1, n2+1, n3+1) into ``fine_points``

``check()`` reports what is wrong instead of raising, so a strategy can show its intermediate state.
"""

from __future__ import annotations

import collections
import json
from dataclasses import dataclass, field

import numpy as np

from hexdp.core.geometry import PatchRole, closest_on_tris

OFF = np.array([(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)])
# the four parallel corner pairs per block axis (0 = i, 1 = j, 2 = k)
PARALLEL = (((0, 1), (3, 2), (7, 6), (4, 5)), ((0, 3), (1, 2), (5, 6), (4, 7)), ((0, 4), (1, 5), (2, 6), (3, 7)))
HEX_EDGES = ((0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7))
# block faces as (axis, side) -> the 4 corners in cyclic order
FACE_CORNERS = {
    (0, 0): (0, 3, 7, 4), (0, 1): (1, 2, 6, 5),
    (1, 0): (0, 1, 5, 4), (1, 1): (3, 2, 6, 7),
    (2, 0): (0, 1, 2, 3), (2, 1): (4, 5, 6, 7),
}


@dataclass(frozen=True)
class Issue:
    # orientation | inverted_cells | nonmanifold_face | hanging_vertex | t_junction | sheet_counts | degenerate
    kind: str
    message: str
    where: tuple = ()


# Issues that make a block structure unusable. A t_junction (1:n block layout over a conforming fine mesh)
# is only a warning: the mesh itself is valid.
INVALID_KINDS = frozenset({"orientation", "inverted_cells", "nonmanifold_face", "hanging_vertex",
                           "sheet_counts", "degenerate"})


def invalid(issues) -> list:
    """The issues among ``issues`` that make the structure invalid."""
    return [i for i in issues if i.kind in INVALID_KINDS]


@dataclass
class FaceInfo:
    role: PatchRole = PatchRole.UNKNOWN
    surface: int | None = None     # geometry surface tag


@dataclass
class Interface:
    kind: str                       # conformal | periodic | ami
    faces_a: list[frozenset]
    faces_b: list[frozenset]
    transform: dict | None = None   # e.g. {"rotation_axis": [...], "angle_deg": 90}


@dataclass
class BlockStructure:
    points: np.ndarray                                   # (n, 3) coarse vertices
    blocks: np.ndarray                                   # (m, 8) vertex indices, VTK order
    edge_cells: dict[frozenset, int] = field(default_factory=dict)
    edge_grading: dict[frozenset, float] = field(default_factory=dict)
    edge_curves: dict[frozenset, np.ndarray] = field(default_factory=dict)
    face_info: dict[frozenset, FaceInfo] = field(default_factory=dict)
    interfaces: list[Interface] = field(default_factory=list)
    fine_points: np.ndarray | None = None
    grids: list[np.ndarray] | None = None
    meta: dict = field(default_factory=dict)   # provenance: source file, inputs a later stage needs

    def __post_init__(self):
        self.points = np.asarray(self.points, float).reshape(-1, 3)
        self.blocks = np.asarray(self.blocks, np.int64).reshape(-1, 8)

    # ---------------------------------------------------------------- topology
    @property
    def n_blocks(self) -> int:
        return len(self.blocks)

    def block_faces(self, b: int) -> dict[tuple[int, int], frozenset]:
        c8 = self.blocks[b]
        return {k: frozenset(int(c8[i]) for i in v) for k, v in FACE_CORNERS.items()}

    def faces(self) -> dict[frozenset, list[tuple[int, int, int]]]:
        """face key -> [(block, axis, side)] of every block that has it."""
        out: dict[frozenset, list] = collections.defaultdict(list)
        for b in range(self.n_blocks):
            for (a, s), key in self.block_faces(b).items():
                out[key].append((b, a, s))
        return dict(out)

    def boundary_faces(self) -> list[frozenset]:
        return [k for k, v in self.faces().items() if len(v) == 1]

    def edges(self) -> dict[frozenset, list[tuple[int, int]]]:
        """edge key -> [(block, axis)]."""
        out: dict[frozenset, list] = collections.defaultdict(list)
        for b, c8 in enumerate(self.blocks):
            for ax, grp in enumerate(PARALLEL):
                for i, j in grp:
                    out[frozenset((int(c8[i]), int(c8[j])))].append((b, ax))
        return dict(out)

    def sheets(self) -> list[dict]:
        """Dual sheets: classes of block edges linked by being parallel in some block.
        Each: {"edges": set of edge keys, "members": [(block, axis)]}, largest first."""
        parent: dict[frozenset, frozenset] = {}

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for c8 in self.blocks:
            for grp in PARALLEL:
                es = [frozenset((int(c8[i]), int(c8[j]))) for i, j in grp]
                for e in es:
                    parent.setdefault(e, e)
                for e in es[1:]:
                    ra, rb = find(es[0]), find(e)
                    if ra != rb:
                        parent[rb] = ra
        cls: dict[frozenset, set] = collections.defaultdict(set)
        for e in parent:
            cls[find(e)].add(e)
        out = []
        for es in cls.values():
            mem = [(b, ax) for b, c8 in enumerate(self.blocks) for ax, grp in enumerate(PARALLEL)
                   if frozenset((int(c8[grp[0][0]]), int(c8[grp[0][1]]))) in es]
            out.append({"edges": es, "members": mem})
        out.sort(key=lambda s: (-len(s["members"]), min(min(e) for e in s["edges"])))
        return out

    # ---------------------------------------------------------------- checks
    def corner_jacobians(self) -> np.ndarray:
        """(m, 8) scaled Jacobian at each corner of the straight coarse hexes."""
        return hex_corner_jacobians(self.points[self.blocks])

    def fine_min_jacobian(self) -> np.ndarray | None:
        """Per block, the minimum corner scaled Jacobian over its fine cells (None without grids)."""
        if not self.grids:
            return None
        out = []
        for G in self.grids:
            out.append(hex_corner_jacobians(self.fine_points[grid_cells(G)]).min())
        return np.array(out)

    def check(self, tol: float = 1e-9) -> list[Issue]:
        issues: list[Issue] = []
        for b, c8 in enumerate(self.blocks):
            if len(set(c8.tolist())) < 8:
                issues.append(Issue("degenerate", f"block {b} repeats a vertex", (b,)))
        fine = self.fine_min_jacobian()
        if fine is not None:
            # the fine grid is the block: straight coarse hexes of curved blocks may fold without harm
            for b in np.where(fine <= 0)[0]:
                issues.append(Issue("inverted_cells", f"block {b}: fine cells with scaled Jacobian {fine[b]:.3g} <= 0",
                                    (int(b),)))
        else:
            J = self.corner_jacobians()
            for b in np.where(J.min(1) <= 0)[0]:
                issues.append(Issue("orientation", f"block {b}: min corner Jacobian {J[b].min():.3g} <= 0", (int(b),)))
        for key, users in self.faces().items():
            if len(users) > 2:
                issues.append(Issue("nonmanifold_face", f"face {sorted(key)} shared by {len(users)} blocks",
                                    tuple(u[0] for u in users)))
        issues += self._hanging_vertices(tol)
        if self.edge_cells:
            for s in self.sheets():
                counts = {self.edge_cells[e] for e in s["edges"] if e in self.edge_cells}
                if len(counts) > 1:
                    issues.append(Issue("sheet_counts", f"edges of one sheet have cell counts {sorted(counts)}",
                                        tuple(sorted({m[0] for m in s["members"]}))))
        return issues

    def _hanging_vertices(self, tol) -> list[Issue]:
        """A vertex that lies on a boundary face of the structure without being one of its corners is a
        T-junction: the faces next to it cannot be conforming. With fine grids the test uses the real
        (curved) face, otherwise the straight coarse quad."""
        scale = max(float(np.ptp(self.points, axis=0).max()), 1e-300)
        faces = self.faces()
        bfaces = [k for k, v in faces.items() if len(v) == 1]
        users = {v for f in bfaces for v in f}
        hits: dict[int, list[list[int]]] = {}
        for f in bfaces:
            b, a, side = faces[f][0]
            surf = self._face_points(b, a, side)
            lo, hi = surf.min(0) - tol * scale, surf.max(0) + tol * scale
            cand = [v for v in users if v not in f and np.all(self.points[v] >= lo) and np.all(self.points[v] <= hi)]
            if not cand:
                continue
            A, B, C = self._face_triangles(b, a, side)
            for v in cand:
                p = self.points[v][None, None]
                r = closest_on_tris(p, A[None], B[None], C[None])[0]
                if np.linalg.norm(r - self.points[v], axis=1).min() <= tol * scale:
                    on_node = bool(self.grids) and np.linalg.norm(surf - self.points[v], axis=1).min() <= tol * scale
                    hits.setdefault(int(v), []).append((sorted(int(x) for x in f), on_node))
        out = []
        for v, fs in sorted(hits.items()):
            faces_ = [f for f, _ in fs]
            if all(n for _, n in fs):
                # the fine meshes conform (a fine node sits there); only the block layout is 1:n
                out.append(Issue("t_junction", f"vertex {v} splits a block edge of boundary face(s) {faces_}; "
                                               "fine nodes coincide", (v,)))
            else:
                out.append(Issue("hanging_vertex", f"vertex {v} lies on boundary face(s) {faces_}", (v,)))
        return out

    def _face_grid(self, b, axis, side):
        """Node coordinates (p, q, 3) of block b's face (axis, side): fine grid if present, else the 2x2 corners."""
        if self.grids:
            G = self.grids[b]
            idx = [slice(None)] * 3
            idx[axis] = -1 if side else 0
            return self.fine_points[G[tuple(idx)]]
        C = self.points[self.blocks[b]][[np.flatnonzero((OFF[:, axis] == side) & (OFF[:, (axis + 1) % 3] == u) &
                                                        (OFF[:, (axis + 2) % 3] == w))[0]
                                         for u in (0, 1) for w in (0, 1)]]
        return C.reshape(2, 2, 3)

    def _face_points(self, b, axis, side):
        return self._face_grid(b, axis, side).reshape(-1, 3)

    def _face_triangles(self, b, axis, side):
        X = self._face_grid(b, axis, side)
        p00, p10, p11, p01 = X[:-1, :-1], X[1:, :-1], X[1:, 1:], X[:-1, 1:]
        A = np.concatenate([p00.reshape(-1, 3), p00.reshape(-1, 3)])
        B = np.concatenate([p10.reshape(-1, 3), p11.reshape(-1, 3)])
        C = np.concatenate([p11.reshape(-1, 3), p01.reshape(-1, 3)])
        return A, B, C

    def _face_cycle(self, key):
        for b in range(self.n_blocks):
            for (a, s), k in self.block_faces(b).items():
                if k == key:
                    return [self.blocks[b][i] for i in FACE_CORNERS[(a, s)]]
        raise KeyError(key)

    # ---------------------------------------------------------------- serialisation
    def to_json(self) -> dict:
        ek = lambda e: "-".join(map(str, sorted(e)))  # noqa: E731
        return {
            "points": self.points.tolist(),
            "blocks": self.blocks.tolist(),
            "edge_cells": {ek(e): int(n) for e, n in self.edge_cells.items()},
            "edge_grading": {ek(e): float(g) for e, g in self.edge_grading.items()},
            "edge_curves": {ek(e): np.asarray(c).tolist() for e, c in self.edge_curves.items()},
            "face_info": {ek(f): {"role": i.role.value, "surface": i.surface} for f, i in self.face_info.items()},
            "interfaces": [{"kind": i.kind, "faces_a": [ek(f) for f in i.faces_a],
                            "faces_b": [ek(f) for f in i.faces_b], "transform": i.transform}
                           for i in self.interfaces],
            "meta": {k: v for k, v in self.meta.items() if isinstance(v, (str, int, float, bool))},
        }

    @classmethod
    def from_json(cls, d: dict) -> "BlockStructure":
        key = lambda s: frozenset(int(x) for x in s.split("-"))  # noqa: E731
        return cls(
            points=np.asarray(d["points"], float),
            blocks=np.asarray(d["blocks"], np.int64),
            edge_cells={key(k): int(v) for k, v in d.get("edge_cells", {}).items()},
            edge_grading={key(k): float(v) for k, v in d.get("edge_grading", {}).items()},
            edge_curves={key(k): np.asarray(v, float) for k, v in d.get("edge_curves", {}).items()},
            face_info={key(k): FaceInfo(PatchRole(v["role"]), v.get("surface"))
                       for k, v in d.get("face_info", {}).items()},
            interfaces=[Interface(i["kind"], [key(f) for f in i["faces_a"]], [key(f) for f in i["faces_b"]],
                                  i.get("transform")) for i in d.get("interfaces", [])],
            meta=dict(d.get("meta", {})),
        )

    def save(self, path) -> None:
        with open(path, "w") as fh:
            json.dump(self.to_json(), fh)

    @classmethod
    def load(cls, path) -> "BlockStructure":
        with open(path) as fh:
            return cls.from_json(json.load(fh))

    def summary(self) -> dict:
        issues = self.check()
        kinds: dict[str, int] = {}
        for i in issues:
            kinds[i.kind] = kinds.get(i.kind, 0) + 1
        J = self.corner_jacobians()
        fine = self.fine_min_jacobian()
        return {"blocks": self.n_blocks, "vertices": len(self.points), "edges": len(self.edges()),
                "boundary_faces": len(self.boundary_faces()), "sheets": len(self.sheets()),
                "fine_cells": int(sum(np.prod(np.array(G.shape) - 1) for G in self.grids)) if self.grids else 0,
                "min_sj_fine": float(fine.min()) if fine is not None else None,
                "min_sj_coarse": float(J.min()) if len(J) else None, "issues": kinds,
                "valid": not invalid(issues)}


def grid_cells(G: np.ndarray) -> np.ndarray:
    """(n, 8) VTK-ordered node ids of the fine hexes of a structured grid G (n1+1, n2+1, n3+1)."""
    return np.stack([G[:-1, :-1, :-1], G[1:, :-1, :-1], G[1:, 1:, :-1], G[:-1, 1:, :-1],
                     G[:-1, :-1, 1:], G[1:, :-1, 1:], G[1:, 1:, 1:], G[:-1, 1:, 1:]], -1).reshape(-1, 8)


def hex_corner_jacobians(X: np.ndarray) -> np.ndarray:
    """Scaled Jacobian at the 8 corners of hexes X (n, 8, 3) in VTK order."""
    nb = {0: (1, 3, 4), 1: (2, 0, 5), 2: (3, 1, 6), 3: (0, 2, 7),
          4: (7, 5, 0), 5: (4, 6, 1), 6: (5, 7, 2), 7: (6, 4, 3)}
    J = np.empty(X.shape[:2])
    for c, (a, b, d) in nb.items():
        u, v, w = X[:, a] - X[:, c], X[:, b] - X[:, c], X[:, d] - X[:, c]
        det = np.einsum("ij,ij->i", np.cross(u, v), w)
        norm = np.linalg.norm(u, axis=1) * np.linalg.norm(v, axis=1) * np.linalg.norm(w, axis=1)
        J[:, c] = det / np.maximum(norm, 1e-300)
    return J
