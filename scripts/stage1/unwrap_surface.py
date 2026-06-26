#!/usr/bin/env python3
"""
Stage 1 / Step A: Unwrap a cylindrical turbine surface (hub or shroud) into a
flat 2D domain.

Hub and shroud of T1_9 are exact cylinders (constant radius). The mapping

    theta = atan2(y, x)
    s     = r * theta        (circumferential arc length)
    t     = z                (axial)

is isometric (distortion-free), so the 2D triangulation is identical to the 3D
one - only the vertex coordinates change. The result is a rectangular-ish domain
with a blade-shaped hole, i.e. exactly the domain_partition 2D case.

Outputs (returned by ``unwrap``):
    points3d   (N,3)  welded 3D coordinates
    st         (N,2)  unwrapped (s, t) coordinates
    tris       (M,3)  triangle connectivity (0-based, into welded points)
    r          float  cylinder radius
    loops      list[list[int]]  ordered boundary node loops (outer first)
    node_dim   (N,)   2 = interior, 1 = boundary, 0 = corner

Run directly to dump a diagnostic PNG/JSON for the hub.
"""

import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import meshio


def _weld(points, tris, decimals=6):
    """Merge coincident STL vertices, remap triangles."""
    key = np.round(points, decimals)
    uniq, inv = np.unique(key, axis=0, return_inverse=True)
    return uniq, inv[tris]


def _boundary_edges(tris):
    """Edges incident to exactly one triangle."""
    cnt = defaultdict(int)
    for a, b, c in tris:
        for u, v in ((a, b), (b, c), (c, a)):
            cnt[(min(u, v), max(u, v))] += 1
    return [e for e, n in cnt.items() if n == 1]


def _trace_loops(bnd_edges):
    """Order boundary edges into closed loops (assumes manifold deg-2 boundary)."""
    adj = defaultdict(list)
    for a, b in bnd_edges:
        adj[a].append(b)
        adj[b].append(a)
    seen = set()
    loops = []
    for start in adj:
        if start in seen:
            continue
        loop = [start]
        seen.add(start)
        prev, cur = None, start
        while True:
            nxts = [n for n in adj[cur] if n != prev]
            if not nxts:
                break
            nxt = nxts[0]
            if nxt == start:
                break
            loop.append(nxt)
            seen.add(nxt)
            prev, cur = cur, nxt
        loops.append(loop)
    return loops


def _polygon_area(st):
    """Signed shoelace area of an ordered (n,2) polygon."""
    x, y = st[:, 0], st[:, 1]
    return 0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)


def _chord_extremes(loop, st):
    """Return the loop node-pair that is farthest apart (blade LE/TE)."""
    pts = st[loop]
    # O(n^2) on a ~100-node loop is fine.
    d2 = np.sum((pts[:, None, :] - pts[None, :, :]) ** 2, axis=-1)
    i, j = np.unravel_index(np.argmax(d2), d2.shape)
    return [loop[i], loop[j]]


def _detect_corners(loop, st, angle_thresh_deg=40.0):
    """Flag loop vertices whose turning angle exceeds threshold as corners."""
    pts = st[loop]
    n = len(loop)
    corners = []
    for i in range(n):
        a = pts[(i - 1) % n]
        b = pts[i]
        c = pts[(i + 1) % n]
        v1 = b - a
        v2 = c - b
        n1, n2 = np.linalg.norm(v1), np.linalg.norm(v2)
        if n1 < 1e-12 or n2 < 1e-12:
            continue
        cosang = np.clip(np.dot(v1, v2) / (n1 * n2), -1.0, 1.0)
        turn = np.degrees(np.arccos(cosang))
        if turn > angle_thresh_deg:
            corners.append(loop[i])
    return corners


def unwrap(stl_path, corner_angle_deg=40.0):
    mesh = meshio.read(str(stl_path))
    raw_pts = mesh.points
    tris_raw = np.vstack([c.data for c in mesh.cells if c.type == "triangle"])
    pts, tris = _weld(raw_pts, tris_raw)

    r_per = np.sqrt(pts[:, 0] ** 2 + pts[:, 1] ** 2)
    r = float(np.mean(r_per))
    if r_per.std() > 1e-3 * r:
        print(f"WARNING: radius not constant (std={r_per.std():.4g}); "
              f"surface may not be a true cylinder.")

    theta = np.arctan2(pts[:, 1], pts[:, 0])
    # Guard against the +/-pi branch cut (safe for T1_9: theta well inside (-pi,pi)).
    if theta.max() - theta.min() > 1.9 * np.pi:
        raise ValueError("theta spans the branch cut; unwrap needs a cut shift.")
    s = r * theta
    t = pts[:, 2]
    st = np.column_stack([s, t])

    loops = _trace_loops(_boundary_edges(tris))
    # Outer loop = largest absolute enclosed area in (s,t).
    loops.sort(key=lambda lp: abs(_polygon_area(st[lp])), reverse=True)

    node_dim = np.full(len(pts), 2, dtype=np.int64)
    loop_corners = []
    for lp in loops:
        for v in lp:
            node_dim[v] = 1
        corners = list(dict.fromkeys(_detect_corners(lp, st, corner_angle_deg)))
        # Every loop must carry >=2 corners so it can be split into segments
        # whose endpoints are termination nodes (StreamlineMerging needs this).
        # A smooth blade hole has no angular corners -> anchor LE/TE instead.
        if len(corners) < 2:
            for v in _chord_extremes(lp, st):
                if v not in corners:
                    corners.append(v)
        loop_corners.append(corners)
        for v in corners:
            node_dim[v] = 0

    return {
        "points3d": pts,
        "st": st,
        "tris": tris,
        "r": r,
        "loops": loops,
        "loop_corners": loop_corners,
        "node_dim": node_dim,
    }


def _diagnostic(stl_path, out_dir):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    data = unwrap(stl_path)
    st, loops, node_dim = data["st"], data["loops"], data["node_dim"]
    print(f"r={data['r']:.4f}  N={len(st)}  tris={len(data['tris'])}")
    print(f"loops: {[len(l) for l in loops]}  (outer first)")
    print(f"corners: {(node_dim == 0).sum()}  boundary: {(node_dim == 1).sum()}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.triplot(st[:, 0], st[:, 1], data["tris"], lw=0.2, color="0.7")
        colors = ["C0", "C1", "C2", "C3"]
        for i, lp in enumerate(loops):
            ring = st[lp + [lp[0]]]
            ax.plot(ring[:, 0], ring[:, 1], colors[i % 4], lw=1.5,
                    label=f"loop {i} (n={len(lp)})")
        cm = node_dim == 0
        ax.scatter(st[cm, 0], st[cm, 1], c="red", s=40, zorder=5, label="corners")
        ax.set_aspect("equal")
        ax.set_xlabel("s = r*theta")
        ax.set_ylabel("t = z")
        ax.legend()
        ax.set_title(f"Unwrapped {Path(stl_path).name}")
        png = out_dir / "unwrap_diagnostic.png"
        fig.savefig(png, dpi=130, bbox_inches="tight")
        print(f"wrote {png}")
    except Exception as e:  # noqa: BLE001
        print(f"plot skipped: {e}")

    meta = {
        "r": data["r"],
        "n_nodes": int(len(st)),
        "n_tris": int(len(data["tris"])),
        "loop_sizes": [len(l) for l in loops],
        "n_corners": int((node_dim == 0).sum()),
    }
    (out_dir / "unwrap_meta.json").write_text(json.dumps(meta, indent=2))


if __name__ == "__main__":
    import sys
    stl = sys.argv[1] if len(sys.argv) > 1 else \
        "/root/repos/block_structured_meshing/T1_9_hub_raw.stl"
    out = sys.argv[2] if len(sys.argv) > 2 else \
        "/root/repos/block_structured_meshing/output/T1_9/hub_stage1"
    _diagnostic(stl, out)
