"""Geometry: boundary entities with patch roles, feature curves and closest-point queries.

One model for every source. A gmsh .msh, a CAD file (STEP/IGES/BREP) and an STL are all loaded
by ``hexdp.io.gmsh_model.load`` into the same ``Geometry``:

    surfaces   dim-2 entities, each with a triangulation, a name and a patch role
    curves     dim-1 entities with a polyline, the surfaces they bound and a kind
               (sharp / smooth / seam / open / free)
    points     dim-0 entities
    topology   curve -> surfaces, point -> curves, surface -> volumes

Roles come from names (physical groups, CAD names) where they exist; otherwise they stay
UNKNOWN until a detector or the user sets them. Every role records where it came from.
The model never decides a role from a label number.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from enum import Enum

import numpy as np
from scipy.spatial import cKDTree


class PatchRole(str, Enum):
    WALL = "wall"
    INLET = "inlet"
    OUTLET = "outlet"
    PERIODIC = "periodic"
    INTERFACE = "interface"
    SYMMETRY = "symmetry"
    UNKNOWN = "unknown"


class RoleSource(str, Enum):
    NAME = "name"            # physical group / CAD name
    HEURISTIC = "heuristic"  # a detector inferred it
    MANUAL = "manual"        # set by the user or an agent
    NONE = "none"


class CurveKind(str, Enum):
    SHARP = "sharp"    # between two surfaces, normals jump: a real feature
    SMOOTH = "smooth"  # between two surfaces, tangent-continuous: a modelling patch boundary
    SEAM = "seam"      # the same surface on both sides (periodic seam)
    OPEN = "open"      # one surface only: boundary of an open shell
    FREE = "free"      # no surface


class Source(str, Enum):
    MSH = "msh"
    CAD = "cad"
    STL = "stl"


# Ordered: the first rule that matches a name token wins. Tokens are the name split at non-alphanumerics,
# lower case; a rule matches a token that starts with one of its prefixes.
ROLE_RULES: tuple[tuple[PatchRole, tuple[str, ...]], ...] = (
    (PatchRole.PERIODIC, ("perio", "cyclic", "cycl")),
    (PatchRole.INTERFACE, ("interface", "ami", "mixing")),
    (PatchRole.INLET, ("inlet", "inflow", "einlass", "eintritt")),
    (PatchRole.OUTLET, ("outlet", "outflow", "auslass", "austritt")),
    (PatchRole.SYMMETRY, ("symm",)),
    (PatchRole.WALL, ("wall", "hub", "shroud", "blade", "nabe", "kranz", "schaufel", "casing")),
)


def role_from_name(name: str | None) -> PatchRole:
    """Patch role suggested by an entity or group name, UNKNOWN if nothing matches."""
    if not name:
        return PatchRole.UNKNOWN
    tokens = [t for t in re.split(r"[^0-9a-zA-Z]+", name.lower()) if t]
    for role, prefixes in ROLE_RULES:
        if any(t.startswith(p) for t in tokens for p in prefixes):
            return role
    return PatchRole.UNKNOWN


@dataclass
class Surface:
    tag: int
    nodes: np.ndarray            # (n, 3)
    tris: np.ndarray             # (m, 3) indices into nodes
    name: str | None = None
    role: PatchRole = PatchRole.UNKNOWN
    role_source: RoleSource = RoleSource.NONE

    @property
    def area(self) -> float:
        a, b, c = (self.nodes[self.tris[:, i]] for i in range(3))
        return float(0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1).sum())


@dataclass
class Curve:
    tag: int
    points: np.ndarray           # (n, 3) polyline in curve order
    surfaces: tuple[int, ...] = ()
    kind: CurveKind = CurveKind.FREE
    angle_deg: float | None = None   # max normal jump along the curve (two-surface curves)
    name: str | None = None

    @property
    def length(self) -> float:
        return float(np.linalg.norm(np.diff(self.points, axis=0), axis=1).sum()) if len(self.points) > 1 else 0.0


@dataclass
class Geometry:
    source: Source
    path: str | None
    surfaces: dict[int, Surface]
    curves: dict[int, Curve]
    points: dict[int, np.ndarray]
    point_curves: dict[int, tuple[int, ...]] = field(default_factory=dict)
    surface_volumes: dict[int, tuple[int, ...]] = field(default_factory=dict)
    exact: bool = False          # closest-point queries can go to the exact CAD surfaces
    notes: list[str] = field(default_factory=list)

    # -------------------------------------------------------------- summary
    def bbox(self) -> tuple[np.ndarray, np.ndarray]:
        P = np.vstack([s.nodes for s in self.surfaces.values()])
        return P.min(0), P.max(0)

    def roles(self) -> dict[PatchRole, list[int]]:
        out: dict[PatchRole, list[int]] = {}
        for s in self.surfaces.values():
            out.setdefault(s.role, []).append(s.tag)
        return out

    def set_role(self, tags, role: PatchRole, source: RoleSource = RoleSource.MANUAL) -> None:
        for t in np.atleast_1d(tags):
            s = self.surfaces[int(t)]
            s.role, s.role_source = PatchRole(role), source

    def feature_curves(self, include_smooth: bool = False) -> list[Curve]:
        keep = {CurveKind.SHARP, CurveKind.OPEN} | ({CurveKind.SMOOTH} if include_smooth else set())
        return [c for c in self.curves.values() if c.kind in keep]

    def summary(self) -> dict:
        lo, hi = self.bbox()
        kinds: dict[str, int] = {}
        for c in self.curves.values():
            kinds[c.kind.value] = kinds.get(c.kind.value, 0) + 1
        return {
            "source": self.source.value,
            "exact": self.exact,
            "bbox": [lo.tolist(), hi.tolist()],
            "surfaces": len(self.surfaces),
            "curves": len(self.curves),
            "points": len(self.points),
            "curve_kinds": kinds,
            "roles": {r.value: sorted(t) for r, t in self.roles().items()},
            "notes": list(self.notes),
        }

    # -------------------------------------------------------------- closest-point queries
    def _surface_index(self, tags):
        key = tuple(sorted(tags))
        cache = self.__dict__.setdefault("_tri_cache", {})
        if key not in cache:
            A, B, C, owner = [], [], [], []
            for t in key:
                s = self.surfaces[t]
                A.append(s.nodes[s.tris[:, 0]]); B.append(s.nodes[s.tris[:, 1]]); C.append(s.nodes[s.tris[:, 2]])
                owner.append(np.full(len(s.tris), t))
            cache[key] = TriSurface(np.vstack(A), np.vstack(B), np.vstack(C)), np.concatenate(owner)
        return cache[key]

    def project(self, Q, surfaces=None):
        """Closest point on the triangulated surfaces (all, or the given tags).

        Returns (points, distances, surface tag per point). The triangulation is the reference; for a CAD
        source ``project_exact`` refines the result onto the exact surfaces."""
        tags = list(self.surfaces) if surfaces is None else [int(t) for t in np.atleast_1d(surfaces)]
        tri, owner = self._surface_index(tags)
        R, d, j = tri.project(Q, return_index=True)
        return R, d, owner[j]

    def project_curve(self, Q, curves=None):
        """Closest point on curve polylines. Returns (points, distances, curve tag per point)."""
        tags = list(self.curves) if curves is None else [int(t) for t in np.atleast_1d(curves)]
        P = [self.curves[t].points for t in tags if len(self.curves[t].points) > 1]
        ids = [np.full(len(self.curves[t].points), t) for t in tags if len(self.curves[t].points) > 1]
        return Polylines(np.vstack(P), np.concatenate(ids)).project(Q)

    def project_exact(self, Q, surfaces=None):
        """Closest point on the exact CAD surfaces (trimmed faces); needs ``exact``."""
        if not self.exact:
            raise ValueError("this geometry has no exact surfaces (source is not CAD)")
        from hexdp.io.gmsh_model import closest_on_cad_faces
        R0, _, owner = self.project(Q, surfaces)
        return closest_on_cad_faces(self.path, np.asarray(Q, float).reshape(-1, 3), owner, R0)


# ------------------------------------------------------------------ point-to-triangle / polyline helpers
def closest_on_tris(p, a, b, c):
    """Closest point on triangles (a, b, c) to points p, vectorised (Ericson, Real-Time Collision Detection)."""
    ab, ac, ap = b - a, c - a, p - a
    d1, d2 = (ab * ap).sum(-1), (ac * ap).sum(-1)
    bp = p - b; d3, d4 = (ab * bp).sum(-1), (ac * bp).sum(-1)
    cp = p - c; d5, d6 = (ab * cp).sum(-1), (ac * cp).sum(-1)
    va = d3 * d6 - d5 * d4; vb = d5 * d2 - d1 * d6; vc = d1 * d4 - d3 * d2
    denom = np.where(va + vb + vc == 0, 1e-300, va + vb + vc)
    v = vb / denom; w = vc / denom
    out = a + ab * v[..., None] + ac * w[..., None]

    def put(mask, val):
        nonlocal out
        out = np.where(mask[..., None], val, out)

    with np.errstate(divide="ignore", invalid="ignore"):
        e_bc = (va <= 0) & (d4 - d3 >= 0) & (d5 - d6 >= 0)
        put(e_bc, b + (c - b) * np.nan_to_num((d4 - d3) / ((d4 - d3) + (d5 - d6)))[..., None])
        e_ac = (vb <= 0) & (d2 >= 0) & (d6 <= 0)
        put(e_ac, a + ac * np.nan_to_num(d2 / (d2 - d6))[..., None])
        e_ab = (vc <= 0) & (d1 >= 0) & (d3 <= 0)
        put(e_ab, a + ab * np.nan_to_num(d1 / (d1 - d3))[..., None])
    put((d6 >= 0) & (d5 <= d6), c)
    put((d3 >= 0) & (d4 <= d3), b)
    put((d1 <= 0) & (d2 <= 0), a)
    return out


class TriSurface:
    """Triangle set with closest-point queries: KD trees on centroids and vertices give candidate
    triangles, the exact point-triangle distance picks the closest (vertices too, so long thin
    boundary-layer triangles are not missed)."""

    def __init__(self, A, B, C, k=12):
        self.A, self.B, self.C = A, B, C
        self.ntri = len(A)
        self.k = min(k, self.ntri)
        self.tree = cKDTree((A + B + C) / 3)
        self.vtree = cKDTree(np.vstack([A, B, C]))

    def project(self, Q, return_index=False):
        Q = np.asarray(Q, float)
        sh = Q.shape
        Q = Q.reshape(-1, 3)
        _, nn = self.tree.query(Q, k=self.k)
        _, nv = self.vtree.query(Q, k=min(self.k, 3 * self.ntri))
        nn = np.concatenate([nn.reshape(len(Q), -1), nv.reshape(len(Q), -1) % self.ntri], 1)
        R = closest_on_tris(Q[:, None], self.A[nn], self.B[nn], self.C[nn])
        d = np.linalg.norm(R - Q[:, None], axis=-1)
        j = d.argmin(1)
        r = np.arange(len(Q))
        out, dist = R[r, j].reshape(sh), d[r, j].reshape(sh[:-1])
        return (out, dist, nn[r, j]) if return_index else (out, dist)


class Polylines:
    """Closest point on polylines (segments between consecutive samples with the same id)."""

    def __init__(self, P, ids):
        ok = ids[1:] == ids[:-1]
        self.a, self.b, self.id = P[:-1][ok], P[1:][ok], ids[:-1][ok]
        self.tree = cKDTree((self.a + self.b) / 2)

    def project(self, Q, k=8):
        Q = np.asarray(Q, float).reshape(-1, 3)
        _, nn = self.tree.query(Q, k=min(k, len(self.a)))
        nn = nn.reshape(len(Q), -1)
        a, b = self.a[nn], self.b[nn]
        ab = b - a
        t = np.clip(((Q[:, None] - a) * ab).sum(-1) / np.maximum((ab * ab).sum(-1), 1e-300), 0, 1)
        R = a + ab * t[..., None]
        d = np.linalg.norm(R - Q[:, None], axis=-1)
        j = d.argmin(1)
        r = np.arange(len(Q))
        return R[r, j], d[r, j], self.id[nn[r, j]]


# ------------------------------------------------------------------ curve classification
def classify_curves(geo: Geometry, sharp_deg: float = 30.0) -> None:
    """Set ``kind`` and ``angle_deg`` of every curve from the surfaces next to it.

    This sorts the curves the source already has; it does not detect features. A curve between two
    surfaces is SHARP when the surface normals across it differ by more than ``sharp_deg`` somewhere
    along it, otherwise SMOOTH. A curve with one adjacent surface is a SEAM when that surface has
    triangles on both sides of it, else OPEN."""
    for c in geo.curves.values():
        surf = tuple(dict.fromkeys(c.surfaces))
        if len(surf) == 0:
            c.kind, c.angle_deg = CurveKind.FREE, None
        elif len(surf) == 1:
            c.kind = CurveKind.SEAM if _two_sided(geo.surfaces[surf[0]], c.points) else CurveKind.OPEN
            c.angle_deg = None
        else:
            ang = _normal_jump(geo, c, surf)
            c.angle_deg = ang
            c.kind = CurveKind.SHARP if ang > sharp_deg else CurveKind.SMOOTH


def _tri_normals(s: Surface):
    a, b, c = (s.nodes[s.tris[:, i]] for i in range(3))
    n = np.cross(b - a, c - a)
    return n / np.maximum(np.linalg.norm(n, axis=1, keepdims=True), 1e-300), (a + b + c) / 3


def _edge_midpoints(c: Curve):
    return 0.5 * (c.points[1:] + c.points[:-1]) if len(c.points) > 1 else c.points


def _normal_jump(geo: Geometry, c: Curve, surf) -> float:
    """Largest angle between the normals of the two nearest triangles of each pair of adjacent surfaces,
    sampled at the curve's segment midpoints. Orientation-free (|cos|), so it measures the crease."""
    M = _edge_midpoints(c)
    if not len(M):
        return 0.0
    normals = []
    for t in surf[:2]:
        n, cen = _tri_normals(geo.surfaces[t])
        _, j = cKDTree(cen).query(M)
        normals.append(n[j])
    cosang = np.abs((normals[0] * normals[1]).sum(1)).clip(0, 1)
    return float(np.degrees(np.arccos(cosang.min())))


def _two_sided(s: Surface, pts: np.ndarray, tol_rel: float = 1e-6) -> bool:
    """Does the surface triangulation have two triangles on most segments of the polyline?"""
    if len(pts) < 2:
        return False
    tree = cKDTree(s.nodes)
    scale = max(np.ptp(s.nodes, axis=0).max(), 1e-300)
    d, idx = tree.query(pts)
    if (d > tol_rel * scale + 1e-12).mean() > 0.5:
        return False     # curve nodes are not surface nodes: cannot tell, call it open
    E = np.sort(np.vstack([s.tris[:, [0, 1]], s.tris[:, [1, 2]], s.tris[:, [2, 0]]]), axis=1)
    _, counts = np.unique(E, axis=0, return_counts=True)
    count = {tuple(e): k for e, k in zip(np.unique(E, axis=0), counts)}
    segs = np.sort(np.stack([idx[:-1], idx[1:]], 1), axis=1)
    seen = [count.get(tuple(e), 0) for e in segs if e[0] != e[1]]
    return bool(seen) and np.mean([k >= 2 for k in seen]) > 0.5
