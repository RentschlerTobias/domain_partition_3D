"""Exporter ``dataset``: one training sample (sample.npz) of a block structure for learning blockings from a
boundary point cloud.

Input of a model (the closed, oriented boundary of the target volume):
    surf_points (n, 3), surf_normals (n, 3) outward unit normals, surf_labels (n,) patch label ids
    -- sampled area-uniformly (seeded) on the boundary triangulation, which is stored as well:
    bnd_vertices (v, 3), bnd_triangles (t, 3) oriented outward, bnd_tri_labels (t,)
Target (the coarse block structure):
    corners (c, 3), blocks (b, 8) VTK hex order, edge_curves (e, k, 3) block edges sampled k times (curved):
    boundary edges are resampled along their chord and projected onto the input surface, because the fine grid keeps
    only 2-5 nodes per edge after the collapse and interpolation alone cannot bend a 2-node cylinder rim,
    edge_ends (e, 2) corner ids, face_labels: boundary faces of blocks (f, 6) = block, axis,
    side, label id, then 2 x unused (kept rectangular)
Meta: label_names (json), blocks/valid/min_sj, case parameters and the source (json strings).
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from hexdp.contracts import ParamSpec


def _boundary_triangles(bs):
    """Outward-oriented triangles of the fine boundary quads, with the face label id per triangle."""
    from hexdp.core.blockstructure import grid_cells  # noqa: F401  (fine grids required)
    if not bs.grids:
        raise ValueError("dataset export needs the fine grids (a filled block structure)")
    P = np.asarray(bs.fine_points, float)
    tris, labs = [], []
    for key, users in bs.faces().items():
        if len(users) != 1:
            continue
        b, ax, sd = users[0]
        G = bs.grids[b]
        idx = [slice(None)] * 3
        idx[ax] = -1 if sd else 0
        F = G[tuple(idx)]
        q = np.stack([F[:-1, :-1], F[1:, :-1], F[1:, 1:], F[:-1, 1:]], -1).reshape(-1, 4)
        info = bs.face_info.get(key)
        lab = int(info.surface) if info is not None and info.surface is not None else 0
        # orientation per quad: outward = from the next inner grid layer to the face
        jdx = [slice(None)] * 3
        jdx[ax] = -2 if sd else 1
        Fi = G[tuple(jdx)]
        qi = np.stack([Fi[:-1, :-1], Fi[1:, :-1], Fi[1:, 1:], Fi[:-1, 1:]], -1).reshape(-1, 4)
        out_dir = P[q].mean(1) - P[qi].mean(1)
        t = np.r_[q[:, [0, 1, 2]], q[:, [0, 2, 3]]]
        od = np.r_[out_dir, out_dir]
        n = np.cross(P[t[:, 1]] - P[t[:, 0]], P[t[:, 2]] - P[t[:, 0]])
        flip = np.einsum("ij,ij->i", n, od) < 0
        t[flip] = t[flip][:, [0, 2, 1]]
        tris.append(t)
        labs.append(np.full(len(t), lab))
    T = np.vstack(tris)
    used, inv = np.unique(T, return_inverse=True)
    return P[used], inv.reshape(-1, 3), np.concatenate(labs)


def sample_boundary(V, T, labels, n, seed=0):
    rng = np.random.default_rng(seed)
    A, B, C = V[T[:, 0]], V[T[:, 1]], V[T[:, 2]]
    N = np.cross(B - A, C - A)
    area = 0.5 * np.linalg.norm(N, axis=1)
    k = rng.choice(len(T), size=n, p=area / area.sum())
    u, v = rng.random(n), rng.random(n)
    swap = u + v > 1
    u[swap], v[swap] = 1 - u[swap], 1 - v[swap]
    pts = A[k] + u[:, None] * (B[k] - A[k]) + v[:, None] * (C[k] - A[k])
    nrm = N[k] / (np.linalg.norm(N[k], axis=1, keepdims=True) + 1e-300)
    return pts, nrm, labels[k]


def _project_edge(project, labeller, Q, a, b):
    """Project the interior samples of one edge onto the surface.

    Refused in two cases: a point that would move by more than half the chord (the projection probably jumped
    across a gap), and a projected line whose points fall on different patches. The second case is a crease or a
    rim - the nearest point then zigzags between the two faces, and the edge needs the seam curve instead (not
    built yet), so the straight chord is the better answer there."""
    P = project(labeller, Q)
    moved = np.linalg.norm(P - Q, axis=1)
    if moved.max() > 0.5 * float(np.linalg.norm(b - a)):
        return Q
    ids = labeller.label(P)          # the two points next to the fixed corners may already touch another patch
    if len(np.unique(ids[1:-1] if len(ids) > 2 else ids)) > 1:
        return Q
    return P


def _seam_polylines(labeller):
    """The label boundaries of the input surface as polylines (points, closed, patch labels)."""
    from hexdp.fill.curve_postprocess import seam_polylines
    return seam_polylines(np.asarray(labeller.P, float), np.asarray(labeller.tris, int),
                          np.asarray(labeller.ids, int))


def _closest_on_polyline(poly, q):
    """(arclength of the closest point, its distance) for one point on a polyline."""
    a, d = poly[:-1], np.diff(poly, axis=0)
    L2 = np.einsum("ij,ij->i", d, d)
    L2[L2 <= 0] = 1e-30
    t = np.clip(np.einsum("ij,ij->i", q - a, d) / L2, 0.0, 1.0)
    p = a + t[:, None] * d
    dist = np.linalg.norm(p - q, axis=1)
    i = int(np.argmin(dist))
    s = np.r_[0.0, np.cumsum(np.sqrt(L2))]
    return float(s[i] + t[i] * np.sqrt(L2[i])), float(dist[i])


def _seam_curve(seams, labels, q, k):
    """The seam curve between an edge's two patches, resampled between the projections of its endpoints.

    A block edge between two patches is a crease or a rim: its nearest point on either face alternates and the
    curve zigzags, so the edge follows the seam of the input triangulation instead. That is what the old tistos
    curve_postprocess did by patch pair. None when no seam with that pair runs near both endpoints."""
    pair = tuple(sorted(labels))
    best = None
    for poly, _closed, labs in seams:
        if tuple(sorted(labs)) != pair or len(poly) < 2:
            continue
        s0, d0 = _closest_on_polyline(poly, q[0])
        s1, d1 = _closest_on_polyline(poly, q[-1])
        score = max(d0, d1)
        if best is None or score < best[0]:
            best = (score, poly, s0, s1)
    if best is None:
        return None
    score, poly, s0, s1 = best
    if score > 0.5 * float(np.linalg.norm(q[-1] - q[0])) or abs(s1 - s0) < 1e-9:
        return None
    s = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(poly, axis=0), axis=1))]
    t = np.linspace(s0, s1, k)
    out = np.c_[tuple(np.interp(t, s, poly[:, c]) for c in range(3))]
    out[0], out[-1] = q[0], q[-1]          # the block corners stay where they are
    return out


def _boundary_edge_labels(bs, labeller=None):
    """Corner pair -> the surface patches of the boundary quads at that edge (two labels = a crease or a rim).

    The surface tag comes from the block structure when the strategy recorded it (algohex does); a structure built
    from a hex mesh without face labels (polycube, or a read blocks.vtk) asks the labeller instead."""
    P = np.asarray(bs.fine_points, float)
    out = {}
    for key, users in bs.faces().items():
        if len(users) != 1:
            continue
        b, ax, sd = users[0]
        info = bs.face_info.get(key)
        lab = -1 if info is None or info.surface is None else int(info.surface)
        G = bs.grids[b]
        idx = [slice(None)] * 3
        idx[ax] = -1 if sd else 0
        F = G[tuple(idx)]
        quad = [F[0, 0], F[-1, 0], F[-1, -1], F[0, -1]]
        cs = [int(np.argmin(np.linalg.norm(bs.points - P[v], axis=1))) for v in quad]
        if lab < 0 and labeller is not None:
            lab = int(labeller.label(np.asarray([P[v] for v in quad]).mean(0, keepdims=True))[0])
        for i in range(4):
            a, c = cs[i], cs[(i + 1) % 4]
            if a != c:
                out.setdefault((min(a, c), max(a, c)), set()).add(lab)
    return out


def _boundary_corner_pairs(bs):
    """Corner pairs of the edges that belong to a boundary block face (a face used by one block).

    Only these edges lie on the geometry. A boundary edge with two fine nodes is a straight chord between two
    surface points - on a cylinder rim that chord cuts off the arc, which is why the edges are projected."""
    return set(_boundary_edge_labels(bs))


def block_edges(bs, k=9, labeller=None):
    """Every block edge once: its two corners and k points along the fine grid line.

    With a ``labeller`` (the input surface), the interior samples of a boundary edge are projected onto the
    geometry. The fine grid keeps only 2-5 nodes per edge after the sheet collapse, so an interpolation can only
    follow the polygon through them: a 4-node cylinder rim stays a square. Projecting the resampled chord puts the
    curvature back without touching the block structure. Endpoints never move (they are block corners), and a
    projection that would move a point more than half the chord is refused (it probably crossed a gap)."""
    seen, ends, curves = set(), [], []
    labels_at = _boundary_edge_labels(bs, labeller) if labeller is not None else {}
    bnd = set(labels_at)
    seams = _seam_polylines(labeller) if labeller is not None else []
    P = np.asarray(bs.fine_points, float)
    for b, G in enumerate(bs.grids):
        n = np.array(G.shape) - 1
        for ax in range(3):
            o = [a for a in range(3) if a != ax]
            for i in (0, n[o[0]]):
                for j in (0, n[o[1]]):
                    idx = [None] * 3
                    idx[o[0]], idx[o[1]], idx[ax] = i, j, slice(None)
                    line = P[G[tuple(idx)]]
                    c0 = int(np.argmin(np.linalg.norm(bs.points - line[0], axis=1)))
                    c1 = int(np.argmin(np.linalg.norm(bs.points - line[-1], axis=1)))
                    key = (min(c0, c1), max(c0, c1), round(float(line[len(line) // 2].sum()), 9))
                    if key in seen:
                        continue
                    seen.add(key)
                    s = np.r_[0, np.cumsum(np.linalg.norm(np.diff(line, axis=0), axis=1))]
                    t = np.linspace(0, s[-1], k)
                    q = np.c_[[np.interp(t, s, line[:, c]) for c in range(3)]].T
                    eid = (min(c0, c1), max(c0, c1))
                    if eid in bnd and k > 2:
                        seam = _seam_curve(seams, labels_at[eid], q, k) if len(labels_at[eid]) == 2 else None
                        if seam is None:
                            from hexdp.core.clean_blocks import project_to_surface
                            q[1:-1] = _project_edge(project_to_surface, labeller, q[1:-1], q[0], q[-1])
                        else:
                            q = seam
                    curves.append(q)
                    ends.append((c0, c1))
    return np.array(ends), np.array(curves)


def _surface_labeller(bs, ctx=None):
    """The input surface the beam worked on, else None.

    The strategies record it in the meta; a workspace that was written by hand (or by a stage that does not) still
    has the tet mesh as tet.vtk, which is the same file."""
    tet = (bs.meta or {}).get("tet_vtk")
    if not tet and ctx is not None:
        cand = Path(ctx.workspace) / "tet.vtk"
        tet = str(cand) if cand.exists() else None
    if not tet or not Path(tet).exists():
        return None
    from hexdp.core.clean_blocks import SurfaceLabeller, read_input_surface
    P0, tri, tid = read_input_surface(str(tet))
    return SurfaceLabeller(P0, tri, tid)


def _edge_off_surface(lab, curves, ends, bnd, part="inner"):
    """Worst distance from the surface of the sampled points of the *boundary* edges.

    "inner" measures the edge curve itself (the corners are fixed and belong to the block structure, not to the
    edge shape); "ends" measures the corners, i.e. how far the block structure sits off the geometry."""
    from hexdp.core.clean_blocks import project_to_surface
    keep = [c for c, (a, b) in zip(curves, ends) if (min(a, b), max(a, b)) in bnd and len(c) > 2]
    if not keep:
        return 0.0
    Q = np.concatenate([c[1:-1] if part == "inner" else c[[0, -1]] for c in keep], axis=0)
    P = project_to_surface(lab, Q)
    return float(np.linalg.norm(Q - P, axis=1).max())


class DatasetExporter:
    """One training sample (sample.npz): boundary point cloud with outward normals and patch labels (input), the
    coarse block structure with curved block edges and boundary face labels (target), metrics and parameters."""

    name = "dataset"
    params = {
        "n_points": ParamSpec(20000, "boundary points of the input point cloud", lower=100),
        "seed": ParamSpec(0, "sampling seed", lower=0),
        "edge_samples": ParamSpec(9, "points per block edge curve", lower=2),
        "project_edges": ParamSpec(True, "project the boundary edge curves onto the input surface (the fine grid "
                                          "keeps only 2-5 nodes per edge, so interpolation alone cannot bend them)"),
    }
    formats = ("npz",)

    def write(self, bs, path, params, ctx):
        out = Path(path).with_suffix(".npz")
        V, T, tl = _boundary_triangles(bs)
        pts, nrm, pl = sample_boundary(V, T, tl, params["n_points"], params["seed"])
        lab = _surface_labeller(bs, ctx) if params.get("project_edges", True) else None
        ends, curves = block_edges(bs, params["edge_samples"], lab)
        faces = []
        for key, users in bs.faces().items():
            if len(users) == 1:
                b, ax, sd = users[0]
                info = bs.face_info.get(key)
                faces.append((b, ax, sd, int(info.surface) if info and info.surface is not None else 0, 0, 0))
        s = bs.summary()
        case = {}
        rep = Path(ctx.workspace) / "case.json"
        if rep.exists():
            case = json.loads(rep.read_text())
        np.savez_compressed(
            out, surf_points=pts.astype(np.float32), surf_normals=nrm.astype(np.float32), surf_labels=pl.astype(np.int16),
            bnd_vertices=V.astype(np.float32), bnd_triangles=T.astype(np.int32), bnd_tri_labels=tl.astype(np.int16),
            corners=np.asarray(bs.points, np.float32), blocks=np.asarray(bs.blocks, np.int32),
            edge_ends=ends.astype(np.int32), edge_curves=curves.astype(np.float32), face_labels=np.array(faces, np.int32),
            label_names=json.dumps({int(t): n for t, n in bs.meta.get("surface_names", {}).items()}),
            metrics=json.dumps({"blocks": s["blocks"], "valid": s["valid"], "min_sj": s["min_sj_fine"],
                                "fine_cells": s["fine_cells"],
                                **({"edge_off_surface_max": round(_edge_off_surface(
                                        lab, curves, ends, _boundary_corner_pairs(bs)), 6),
                                    "corner_off_surface_max": round(_edge_off_surface(
                                        lab, curves, ends, _boundary_corner_pairs(bs), "ends"), 6)}
                                   if lab is not None else {})}),
            case=json.dumps(case))
        return [out]
