"""Exporter ``dataset``: one training sample (sample.npz) of a block structure for learning blockings from a
boundary point cloud.

Input of a model (the closed, oriented boundary of the target volume):
    surf_points (n, 3), surf_normals (n, 3) outward unit normals, surf_labels (n,) patch label ids
    -- sampled area-uniformly (seeded) on the boundary triangulation, which is stored as well:
    bnd_vertices (v, 3), bnd_triangles (t, 3) oriented outward, bnd_tri_labels (t,)
Target (the coarse block structure):
    corners (c, 3), blocks (b, 8) VTK hex order, edge_curves (e, k, 3) block edges sampled k times along the fine
    grid line (curved), edge_ends (e, 2) corner ids, face_labels: boundary faces of blocks (f, 6) = block, axis,
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


def block_edges(bs, k=9):
    """Every block edge once: its two corners and k points along the fine grid line."""
    seen, ends, curves = set(), [], []
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
                    curves.append(np.c_[[np.interp(t, s, line[:, c]) for c in range(3)]].T)
                    ends.append((c0, c1))
    return np.array(ends), np.array(curves)


class DatasetExporter:
    """One training sample (sample.npz): boundary point cloud with outward normals and patch labels (input), the
    coarse block structure with curved block edges and boundary face labels (target), metrics and parameters."""

    name = "dataset"
    params = {
        "n_points": ParamSpec(20000, "boundary points of the input point cloud", lower=100),
        "seed": ParamSpec(0, "sampling seed", lower=0),
        "edge_samples": ParamSpec(9, "points per block edge curve", lower=2),
    }
    formats = ("npz",)

    def write(self, bs, path, params, ctx):
        out = Path(path).with_suffix(".npz")
        V, T, tl = _boundary_triangles(bs)
        pts, nrm, pl = sample_boundary(V, T, tl, params["n_points"], params["seed"])
        ends, curves = block_edges(bs, params["edge_samples"])
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
                                "fine_cells": s["fine_cells"]}),
            case=json.dumps(case))
        return [out]
