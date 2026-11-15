"""B5: conformal fine mesh of the 12 core blocks, BL at hub and shroud.

    core_tfi.py <run_dir> [--target-h 0.05] [--first 0.0009] [--growth 1.2]

Conformal by construction instead of by repair:
  * each block is evaluated as a tensor-product UNIFORM Catmull-Rom spline over
    its own coarse lattice (blade/blocks_b3.vtk), in cylinder coordinates
    (r, theta, z) about the machine axis. The spline interpolates, so on
    a block face it depends on that face's nodes only, and uniform CR is
    symmetric under index reversal / transposition: two blocks evaluate the face
    they share identically (welded at 1e-9 afterwards, checked watertight);
  * boundary faces are then projected onto their exact surface (geometry.py:
    dtOO hub / shroud / blade, msh planes), points on a face edge shared by two
    boundary faces onto both surfaces (alternating); the projection displacement
    is carried into the block by TFI (tfi.tfi on the displacement field);
  * divisions: tfi.direction_classes + tfi.solve_block_divisions at --target-h;
    a class whose axes end on hub / shroud gets BL layers there (first cell
    --first, geometric --growth until the core spacing). The distribution is
    oriented by RADIUS, not by the block's index direction, so both blocks of a
    shared face see the same one (the reason tfi.clustered_fractions was never
    used: a one-sided distribution in index space tears the seam).
The ring blocks (B3, ids >= core count) are not meshed here (B6); their faces
towards the core are the cutout faces, evaluated by the same rule.
"""
import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))

import tfi                                                            # noqa: E402
import quality as qual                                                # noqa: E402
from geometry import Geometry                                         # noqa: E402

WALLS = ("hub", "shroud")


def lin_weights(n_coarse, s):
    """(len(s), n_coarse + 1) piecewise-linear weights (no overshoot)."""
    W = np.zeros((len(s), n_coarse + 1))
    for r, x in enumerate(s):
        i = min(int(np.floor(x)), n_coarse - 1)
        t = x - i
        W[r, i], W[r, i + 1] = 1 - t, t
    return W


def cr_weights(n_coarse, s):
    """(len(s), n_coarse + 1) uniform Catmull-Rom weights at index params s in
    [0, n_coarse]; reflected ghost nodes at both ends (symmetric under reversal)."""
    m = n_coarse
    W = np.zeros((len(s), m + 1))
    if m == 1:
        W[:, 0], W[:, 1] = 1 - s, s
        return W
    for r, x in enumerate(s):
        i = min(int(np.floor(x)), m - 1)
        t = x - i
        b = np.array([-t ** 3 + 2 * t ** 2 - t, 3 * t ** 3 - 5 * t ** 2 + 2,
                      -3 * t ** 3 + 4 * t ** 2 + t, t ** 3 - t ** 2]) / 2.0
        for k, idx in enumerate((i - 1, i, i + 1, i + 2)):
            if idx < 0:            # ghost = 2 X0 - X1
                W[r, 0] += 2 * b[k]
                W[r, 1] -= b[k]
            elif idx > m:          # ghost = 2 Xm - Xm-1
                W[r, m] += 2 * b[k]
                W[r, m - 1] -= b[k]
            else:
                W[r, idx] += b[k]
    return W


def graded(n, f_lo, f_hi, growth):
    """n cells on [0, 1]: geometric from f_lo (if > 0) at 0 and f_hi at 1 towards
    the uniform size, uniform in between; exactly n cells, normalised."""
    h = np.full(n, 1.0 / n)
    for _ in range(50):
        u = h.copy()
        core = 1.0 / n
        if f_lo > 0:
            k, x = 0, f_lo
            while k < n and x < core:
                u[k] = x
                x *= growth
                k += 1
        if f_hi > 0:
            k, x = n - 1, f_hi
            while k >= 0 and x < core:
                u[k] = min(u[k], x) if f_lo > 0 else x
                x *= growth
                k -= 1
        u /= u.sum()
        if np.allclose(u, h, rtol=1e-12, atol=1e-15):
            break
        h = u
    return np.r_[0, np.cumsum(h)] / np.sum(h)


def tfi_param(X, u, v, w):
    """tfi.tfi with the blending parameters of the actual node distribution
    (u, v, w in [0, 1] per axis) instead of the index: a boundary displacement
    then reaches a BL layer in proportion to its physical distance from the
    wall, so the 0.9 mm first cell moves with the wall instead of being sheared."""
    u = np.asarray(u)[:, None, None, None]
    v = np.asarray(v)[None, :, None, None]
    w = np.asarray(w)[None, None, :, None]
    F = ((1 - u) * X[0][None] + u * X[-1][None]
         + (1 - v) * X[:, 0][:, None] + v * X[:, -1][:, None]
         + (1 - w) * X[:, :, 0][:, :, None] + w * X[:, :, -1][:, :, None])
    bu, bv, bw, iu = [(1 - u), u], [(1 - v), v], [(1 - w), w], [0, -1]
    E = np.zeros_like(F)
    for a in (0, 1):
        for b in (0, 1):
            E += bu[a] * bv[b] * X[iu[a], iu[b], :][None, None]
            E += bu[a] * bw[b] * X[iu[a], :, iu[b]][None, :, None]
            E += bv[a] * bw[b] * X[:, iu[a], iu[b]][:, None, None]
    V = np.zeros_like(F)
    for a in (0, 1):
        for b in (0, 1):
            for c in (0, 1):
                V += bu[a] * bv[b] * bw[c] * X[iu[a], iu[b], iu[c]]
    return F - E + V


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run")
    ap.add_argument("--blocks", default=None)
    ap.add_argument("--target-h", type=float, default=0.05)
    ap.add_argument("--first", type=float, default=0.0009)
    ap.add_argument("--growth", type=float, default=1.2)
    ap.add_argument("--interp", choices=("spline", "linear"), default="spline",
                    help="coarse-lattice interpolation in (r, theta, z)")
    ap.add_argument("--no-project", action="store_true", help="diagnostic: skip the surface projection")
    ap.add_argument("--periodic", choices=("circumferential", "normal"), default="circumferential",
                    help="how points are put on the periodic surfaces")
    ap.add_argument("--deviation-vtk", default=None,
                    help="diagnostic: write the spline-vs-surface distance of every boundary face")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    run = Path(a.run)
    blocks = Path(a.blocks) if a.blocks else run / "blade" / "blocks_b3.vtk"
    out = Path(a.out) if a.out else run / "blade" / "tfi_b5_core.vtk"
    b3 = json.loads((run / "blade" / "b3_ring.json").read_text())
    vs = {int(k): set(v) for k, v in b3["vertex_surfaces"].items()}
    ncore = int(b3["core_blocks"])

    P, H, B, f2h = tfi.load_blocks(str(blocks))
    lat, missing = tfi.lattices(P, H, B, f2h, verbose=False)
    classes = list(tfi.direction_classes(lat, f2h, H, B, verbose=False))
    classes = [list(c) for c in (classes.values() if isinstance(classes, dict) else classes)]
    counts = tfi.solve_block_divisions(lat, classes, P, a.target_h, verbose=False)
    cof = tfi.class_of_axis(classes)
    G = Geometry(run)

    def face_walls(vert, ax, side):
        g = tfi.face_grid(vert, ax, side)
        common = set.intersection(*[vs.get(int(v), set()) for v in g.ravel()])
        return common

    # per class: is there a wall at the low-r / high-r end, and the physical length
    cls_info = {}
    for ci, c in enumerate(classes):
        lo = hi = False
        L = []
        for r, ax in c:
            vert = lat[r][1]
            for side in (0, 1):
                w = face_walls(vert, ax, side) & set(WALLS)
                if w:
                    g = tfi.face_grid(vert, ax, side)
                    rr = np.hypot(*P[g.ravel()][:, :2].T).mean()
                    lo |= "hub" in w
                    hi |= "shroud" in w
            L.append(tfi.axis_length(P, vert, ax))
        cls_info[ci] = {"hub": lo, "shroud": hi, "len_min": float(min(L)), "len_max": float(max(L))}
    # BL layers: cells added so the geometric run reaches the core spacing
    fracs = {}
    for ci, info in cls_info.items():
        n = int(counts[ci])
        if info["hub"] or info["shroud"]:
            h_core = info["len_min"] / n
            nbl = int(np.ceil(np.log(max(h_core / a.first, 1.0)) / np.log(a.growth)))
            n += nbl * (int(info["hub"]) + int(info["shroud"]))
            f = a.first / info["len_min"]
            fracs[ci] = graded(n, f if info["hub"] else 0, f if info["shroud"] else 0, a.growth)
        else:
            fracs[ci] = np.linspace(0, 1, n + 1)
        info["cells"] = len(fracs[ci]) - 1

    pts, hexes, bids = [], [], []
    dev_pts, dev_val, dev_sid, dev_blk, dev_q = [], [], [], [], []
    base = 0
    proj_rep = defaultdict(float)
    proj_blk = {}
    for r in sorted(lat):
        if r >= ncore:
            continue
        dims, vert = lat[r]
        X = P[vert]                                            # (n0+1, n1+1, n2+1, 3)
        Ws, ss = [], []
        for ax in range(3):
            f = fracs[cof[(r, ax)]].copy()
            # orient by radius: the class fractions run low r -> high r
            sl = [0, 0, 0]
            sl_end = [0, 0, 0]
            sl_end[ax] = -1
            r0 = np.hypot(*X[tuple(sl)][:2])
            r1 = np.hypot(*X[tuple(sl_end)][:2])
            if r0 > r1 + 1e-12:
                f = 1 - f[::-1]
            nc = vert.shape[ax] - 1
            ss.append(f)
            Ws.append((cr_weights if a.interp == "spline" else lin_weights)(nc, f * nc))
        # spline in cylinder coordinates about the machine axis (z): a spline through
        # hub nodes then stays on r = r_hub instead of sagging 13 mm into the hub
        # between coarse nodes 0.3 m apart (that sag sheared the BL cells)
        th0 = np.arctan2(*X[..., :2].reshape(-1, 2).mean(0)[::-1])
        Cy = np.stack([np.hypot(X[..., 0], X[..., 1]),
                       np.angle(np.exp(1j * (np.arctan2(X[..., 1], X[..., 0]) - th0))),
                       X[..., 2]], -1)
        Fc = np.einsum("ai,bj,ck,ijkx->abcx", Ws[0], Ws[1], Ws[2], Cy)
        F = np.stack([Fc[..., 0] * np.cos(Fc[..., 1] + th0), Fc[..., 0] * np.sin(Fc[..., 1] + th0),
                      Fc[..., 2]], -1)
        # boundary faces of this block: constraint per fine boundary point
        cons = {}
        for ax in range(3):
            for side in (0, 1):
                w = face_walls(vert, ax, side) - {"cutout"}
                if not w:
                    continue
                name = sorted(w)[0] if len(w) == 1 else None
                if name is None:
                    continue
                sl = [slice(None)] * 3
                sl[ax] = 0 if side == 0 else -1
                idx = np.argwhere(np.ones(F.shape[:3], bool)[tuple(sl)])
                for q in idx:
                    full = list(q)
                    full.insert(ax, 0 if side == 0 else F.shape[ax] - 1)
                    cons.setdefault(tuple(full), set()).add(name)
        D = np.zeros_like(F)
        for key, names in cons.items():
            x = F[key]
            if len(names) >= 3:
                continue                                       # block corner: coarse node, exact
            y = x.copy()
            # periodic last: the other surface first (hub / shroud: radial, inlet /
            # outlet: axial), then the circumferential move onto the periodic one
            order = sorted(names, key=lambda n: n.startswith("periodic"))
            for _ in range(6 if len(names) == 2 else 1):
                for n in order:
                    y = (G.project_any if a.periodic == "circumferential" else G.project)(n, y)[0]
            D[key] = y - x
            kk_ = ",".join(sorted(names))
            proj_rep[kk_] = max(proj_rep[kk_], float(np.linalg.norm(y - x)))
            proj_blk[(r, kk_)] = max(proj_blk.get((r, kk_), 0.0), float(np.linalg.norm(y - x)))
        if a.deviation_vtk:
            # spline vs exact surface on every constrained boundary face of the block
            for ax in range(3):
                for side in (0, 1):
                    w = face_walls(vert, ax, side) - {"cutout"}
                    if len(w) != 1:
                        continue
                    sl = [slice(None)] * 3
                    sl[ax] = 0 if side == 0 else -1
                    Xf, Df = F[tuple(sl)], D[tuple(sl)]
                    n0, n1 = Xf.shape[:2]
                    off = sum(len(x) for x in dev_pts)
                    gi = np.arange(n0 * n1).reshape(n0, n1) + off
                    dev_pts.append(Xf.reshape(-1, 3))
                    dev_val.append(np.linalg.norm(Df, axis=-1).ravel())
                    dev_sid.append(np.full(n0 * n1, sorted(G.surfaces).index(next(iter(w))), float))
                    dev_blk.append(np.full((n0 - 1) * (n1 - 1), r, float))
                    dev_q.append(np.c_[gi[:-1, :-1].ravel(), gi[1:, :-1].ravel(), gi[1:, 1:].ravel(), gi[:-1, 1:].ravel()])
        D = tfi_param(D, *ss)                                  # carry the displacement inside
        if not a.no_project:
            F = F + D
        ni, nj, nk = F.shape[:3]
        idx = np.arange(ni * nj * nk).reshape(ni, nj, nk) + base
        pts.append(F.reshape(-1, 3))
        hexes.append(tfi.block_cells(idx))
        bids += [r] * ((ni - 1) * (nj - 1) * (nk - 1))
        base += ni * nj * nk
    Pf = np.vstack(pts)
    Hf = np.vstack(hexes)
    Pw, remap = tfi.weld(Pf, 1e-9)
    Hw = remap[Hf]
    Bw = np.array(bids)
    watertight, nbnd = tfi.check_watertight(Pw, Hw, verbose=False)
    qs = qual.summary(Pw, Hw)
    cl = qual.classify(Pw, Hw)
    # first cell heights at hub / shroud
    fc = {}
    rr = np.hypot(Pw[:, 0], Pw[:, 1])
    for name, r0 in (("hub", G.r_hub), ("shroud", G.r_shroud)):
        on = np.abs(rr - r0) < 1e-5
        hs = []
        for c in Hw:
            o = on[c]
            if o.sum() == 4:
                top = Pw[c[~o]].mean(0)
                hs.append(abs(np.hypot(*top[:2]) - r0))        # next layer's radial distance
        fc[name] = {"cells": len(hs), "min": float(np.min(hs)), "median": float(np.median(hs)),
                    "max": float(np.max(hs))} if hs else None
    rep = {"blocks": str(blocks), "out": str(out), "core_blocks": ncore,
           "cells": int(len(Hw)), "points": int(len(Pw)), "watertight": bool(watertight),
           "classes": [{"class": ci, "axes": len(classes[ci]), **cls_info[ci]} for ci in cls_info],
           "projection_max": dict(proj_rep),
           "projection_max_per_block": {f"{b}:{k}": v for (b, k), v in sorted(proj_blk.items())},
           "quality": qs, "first_cell_height": fc,
           "gate": {"watertight": bool(watertight), "inverted_0": qs["inverted"] == 0}}
    import meshio
    if a.deviation_vtk and dev_pts:
        import meshio
        meshio.write(a.deviation_vtk, meshio.Mesh(np.vstack(dev_pts), [("quad", np.vstack(dev_q))],
                                                  point_data={"deviation": np.concatenate(dev_val),
                                                              "surface": np.concatenate(dev_sid)},
                                                  cell_data={"block_id": [np.concatenate(dev_blk)]}))
        rep["surface_ids"] = {i: n for i, n in enumerate(sorted(G.surfaces))}
        rep["deviation_per_surface"] = {}
        sid, val = np.concatenate(dev_sid), np.concatenate(dev_val)
        for i, n in enumerate(sorted(G.surfaces)):
            if (sid == i).any():
                v = val[sid == i]
                rep["deviation_per_surface"][n] = {"max": float(v.max()), "p95": float(np.percentile(v, 95)),
                                                   "median": float(np.median(v))}
        out.with_suffix(".json").write_text(json.dumps(rep, indent=1))
    import meshio
    meshio.write(str(out), meshio.Mesh(Pw, [("hexahedron", Hw)],
                                       cell_data={"block_id": [Bw], "sj": [cl["sj"]],
                                                  "inverted": [cl["inverted"].astype(float)],
                                                  "concave": [cl["concave"].astype(float)]}))
    out.with_suffix(".json").write_text(json.dumps(rep, indent=1))
    print(json.dumps({k: rep[k] for k in ("cells", "watertight", "projection_max", "quality",
                                          "first_cell_height", "gate")}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
