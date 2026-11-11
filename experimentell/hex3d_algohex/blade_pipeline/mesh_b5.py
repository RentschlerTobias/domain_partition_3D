"""B5: fine mesh of the core blocks, edge -> face -> volume, every entity once.

    mesh_b5.py <run_dir> [--target-h 0.05] [--first 0.0009] [--growth 1.2]

The coarse lattice (blade/blocks_b3.vtk) only provides TOPOLOGY. Geometry comes
from the exact surfaces and curves:
  edges   the B4 curves (blade/edges_b4.npz: on the dtOO feature curves /
          surfaces), resampled by arc length; span edges that start or end on
          hub / shroud get the BL distribution with the absolute first cell
          --first at that end (per edge, so y1 is the same everywhere);
  faces   Coons patch of the 4 discretised edges in cylinder coordinates
          (r, theta, z), then put onto the face's surface: periodic faces by
          theta = theta_per(r, z) (circumferential, r and z untouched), the other
          surfaces by closest point (a small correction: hub / shroud / inlet /
          outlet are revolution surfaces / planes, exact in (r, theta, z) up to
          the Coons blend); interior faces stay the Coons patch;
  volume  Gordon-Hall TFI of the block's 6 finished faces.
Each edge and face is computed once and handed to every block that owns it
(in that block's orientation), so the mesh is conformal by construction.
Divisions: tfi.direction_classes + tfi.solve_block_divisions, BL classes get
the extra layers (core_tfi.class_cells). Ring blocks (ids >= core count) are B6;
the faces they share with the core (cutout) are stored in blade/b5_entities.npz.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))

import tfi                                                            # noqa: E402
import quality as qual                                                # noqa: E402
from core_tfi import graded, tfi_param                                # noqa: E402
from geometry import Geometry                                         # noqa: E402
from winslow import winslow                                           # noqa: E402

WALLS = ("hub", "shroud")


def cyl(X):
    return np.stack([np.hypot(X[..., 0], X[..., 1]), np.arctan2(X[..., 1], X[..., 0]), X[..., 2]], -1)


def cart(C, th0=0.0):
    return np.stack([C[..., 0] * np.cos(C[..., 1] + th0), C[..., 0] * np.sin(C[..., 1] + th0), C[..., 2]], -1)


def resample(poly, f):
    """Points at arc-length fractions f along a polyline."""
    L = np.r_[0, np.cumsum(np.linalg.norm(np.diff(poly, axis=0), axis=1))]
    t = f * L[-1]
    return np.stack([np.interp(t, L, poly[:, k]) for k in range(3)], -1), float(L[-1])


def coons(B, T, L, R):
    """Coons patch from edges B (s, j=0), T (s, j=m), L (i=0, t), R (i=n, t), all as
    arrays of points; blending parameters = mean arc-length fractions of the two
    opposite edges (the boundary is reproduced exactly for any choice)."""
    def frac(E):
        d = np.r_[0, np.cumsum(np.linalg.norm(np.diff(E, axis=0), axis=1))]
        return d / d[-1]
    u = 0.5 * (frac(B) + frac(T))[:, None, None]
    v = 0.5 * (frac(L) + frac(R))[None, :, None]
    P = ((1 - v) * B[:, None] + v * T[:, None] + (1 - u) * L[None] + u * R[None]
         - (1 - u) * (1 - v) * B[0] - u * (1 - v) * B[-1] - (1 - u) * v * T[0] - u * v * T[-1])
    return P


def class_cells(lat, classes, counts, vs, P, first, growth):
    """Cells per class incl. BL layers at hub / shroud ends (as in core_tfi)."""
    out = {}
    for ci, c in enumerate(classes):
        hub = shroud = False
        L = []
        for r, ax in c:
            vert = lat[r][1]
            for side in (0, 1):
                g = tfi.face_grid(vert, ax, side)
                common = set.intersection(*[vs.get(int(v), set()) for v in g.ravel()])
                hub |= "hub" in common
                shroud |= "shroud" in common
            L.append(tfi.axis_length(P, vert, ax))
        n = int(counts[ci])
        if hub or shroud:
            nbl = int(np.ceil(np.log(max(min(L) / n / first, 1.0)) / np.log(growth)))
            n += nbl * (int(hub) + int(shroud))
        out[ci] = {"cells": n, "hub": hub, "shroud": shroud}
    return out


def match(stored, wanted):
    """The dihedral transform (as a function on 2D arrays) taking the stored
    coarse id grid onto the wanted one."""
    for tr in (False, True):
        A = stored.T if tr else stored
        for fa in (False, True):
            for fb in (False, True):
                B = A[::-1] if fa else A
                B = B[:, ::-1] if fb else B
                if B.shape == wanted.shape and np.array_equal(B, wanted):
                    def f(X, tr=tr, fa=fa, fb=fb):
                        Y = np.swapaxes(X, 0, 1) if tr else X
                        Y = Y[::-1] if fa else Y
                        return Y[:, ::-1] if fb else Y
                    return f
    raise RuntimeError("face orientation not found")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run")
    ap.add_argument("--target-h", type=float, default=0.05)
    ap.add_argument("--first", type=float, default=0.0009)
    ap.add_argument("--growth", type=float, default=1.2)
    ap.add_argument("--winslow", type=int, default=300, help="Winslow iterations per block (0 = off)")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    run = Path(a.run)
    bd = run / "blade"
    out = Path(a.out) if a.out else bd / "mesh_b5_core.vtk"
    b3 = json.loads((bd / "b3_ring.json").read_text())
    vs = {int(k): set(v) for k, v in b3["vertex_surfaces"].items()}
    ncore = int(b3["core_blocks"])
    P, H, B, f2h = tfi.load_blocks(str(bd / "blocks_b3.vtk"))
    lat, missing = tfi.lattices(P, H, B, f2h, verbose=False)
    classes = tfi.direction_classes(lat, f2h, H, B, verbose=False)
    classes = [list(c) for c in (classes.values() if isinstance(classes, dict) else classes)]
    counts = tfi.solve_block_divisions(lat, classes, P, a.target_h, verbose=False)
    cof = tfi.class_of_axis(classes)
    cc = class_cells(lat, classes, counts, vs, P, a.first, a.growth)
    G = Geometry(run)

    # B4 curves by their coarse node sequence
    e4 = np.load(bd / "edges_b4.npz")
    b4 = {}
    for i in range(len(e4["nodes"])):
        ids = e4["edge_ids"][i]
        ids = tuple(int(x) for x in ids[ids >= 0])
        b4[tuple(sorted(ids))] = (ids, e4["polyline"][e4["offset"][i]:e4["offset"][i + 1]], str(e4["kind"][i]))

    # ---------------------------------------------------------------- edges
    edges = {}
    erep = {"max_proj": 0.0}

    def edge(ids, n):
        """Fine points of the edge with coarse ids `ids` (that orientation), n cells."""
        key = tuple(sorted(int(x) for x in ids))
        if key not in edges:
            ids0, poly, kind = b4[key]
            ends = [vs.get(ids0[0], set()), vs.get(ids0[-1], set())]
            common = set.intersection(*[vs.get(int(v), set()) for v in ids0])
            L_ = float(np.sum(np.linalg.norm(np.diff(poly, axis=0), axis=1)))
            lo = any(w in ends[0] and w not in common for w in WALLS)
            hi = any(w in ends[1] and w not in common for w in WALLS)
            f = graded(n, a.first / L_ if lo else 0, a.first / L_ if hi else 0, a.growth) \
                if (lo or hi) else np.linspace(0, 1, n + 1)
            X, _ = resample(poly, f)
            X[0], X[-1] = P[ids0[0]], P[ids0[-1]]
            # back onto the constraint (the polyline chords between B4 samples)
            names = [x for x in kind.split(":")[1].split("+")] if ":" in kind else []
            Y = X.copy()
            if names and "blade" not in names:
                for _ in range(6 if len(names) == 2 else 1):
                    for nm in sorted(names, key=lambda s: s.startswith("periodic")):
                        Y[1:-1] = G.project_any(nm, Y[1:-1])
                erep["max_proj"] = max(erep["max_proj"], float(np.abs(Y - X).max()))
            edges[key] = (tuple(ids0), Y, f)
        ids0, Y, f = edges[key]
        return (Y, f) if tuple(int(x) for x in ids) == ids0 else (Y[::-1], 1 - f[::-1])

    # ---------------------------------------------------------------- faces
    faces = {}
    frep = {}

    def face(g, n0, n1):
        """Fine (n0+1, n1+1, 3) points of the face with coarse id grid g."""
        key = frozenset(int(x) for x in g.ravel())
        if key not in faces:
            Bt, _ = edge(g[:, 0], n0)
            Tp, _ = edge(g[:, -1], n0)
            Lf, _ = edge(g[0, :], n1)
            Rt, _ = edge(g[-1, :], n1)
            th0 = float(np.angle(np.exp(1j * np.arctan2(P[g.ravel()][:, 1], P[g.ravel()][:, 0])).mean()))
            rel = lambda X: (lambda C: np.stack([C[..., 0], np.angle(np.exp(1j * (C[..., 1] - th0))), C[..., 2]], -1))(cyl(X))
            F = cart(coons(rel(Bt), rel(Tp), rel(Lf), rel(Rt)), th0)
            common = set.intersection(*[vs.get(int(v), set()) for v in g.ravel()])
            # the cutout is a surface for the core (its outer face towards the ring)
            if not common:
                cut = all(np.linalg.norm(G.project("cutout", P[g.ravel()]) - P[g.ravel()], axis=1) < 1e-6)
                common = {"cutout"} if cut else set()
            name = sorted(common)[0] if len(common) == 1 else None
            inner = F[1:-1, 1:-1]
            if name is not None and inner.size:
                Y = G.project_any(name, inner.reshape(-1, 3)).reshape(inner.shape)
                frep[name] = max(frep.get(name, 0.0), float(np.linalg.norm(Y - inner, axis=-1).max()))
                F[1:-1, 1:-1] = Y
            faces[key] = (g.copy(), F, name)
        g0, F, name = faces[key]
        return match(g0, g)(F)

    # ---------------------------------------------------------------- volumes
    pts, hexes, bids = [], [], []
    wrep = {}
    base = 0
    for r in sorted(lat):
        if r >= ncore:
            continue
        dims, vert = lat[r]
        n = [cc[cof[(r, ax)]]["cells"] for ax in range(3)]
        X = np.full((n[0] + 1, n[1] + 1, n[2] + 1, 3), np.nan)
        for ax in range(3):
            o = [x for x in range(3) if x != ax]
            for side in (0, 1):
                g = tfi.face_grid(vert, ax, side)
                F = face(g, n[o[0]], n[o[1]])
                sl = [slice(None)] * 3
                sl[ax] = 0 if side == 0 else -1
                X[tuple(sl)] = F
        # blending parameters per axis: mean fractions of the 4 parallel edges
        par = []
        for ax in range(3):
            fs = []
            for i in (0, -1):
                for j in (0, -1):
                    sl = [i, j]
                    sl.insert(ax, slice(None))
                    _, f = edge(vert[tuple(sl)], n[ax])
                    fs.append(f)
            par.append(np.mean(fs, axis=0))
        Xi = tfi_param(np.nan_to_num(X), *par)
        bnd = ~np.isnan(X[..., 0])
        Xi[bnd] = X[bnd]
        if a.winslow > 0:
            Xi, wrep[int(r)] = winslow(Xi, iters=a.winslow)
        ni, nj, nk = Xi.shape[:3]
        idx = np.arange(ni * nj * nk).reshape(ni, nj, nk) + base
        pts.append(Xi.reshape(-1, 3))
        hexes.append(tfi.block_cells(idx))
        bids += [r] * ((ni - 1) * (nj - 1) * (nk - 1))
        base += ni * nj * nk
    Pf, Hf = np.vstack(pts), np.vstack(hexes)
    Pw, remap = tfi.weld(Pf, 1e-9)
    Hw = remap[Hf]
    watertight, _ = tfi.check_watertight(Pw, Hw, verbose=False)
    qs = qual.summary(Pw, Hw)
    cl = qual.classify(Pw, Hw)
    fc = {}
    rr = np.hypot(Pw[:, 0], Pw[:, 1])
    for nm, r0 in (("hub", G.r_hub), ("shroud", G.r_shroud)):
        on = np.abs(rr - r0) < 1e-5
        cnt = on[Hw].sum(1)
        hs = [abs(rr[c][~on[c]].mean() - r0) for c in Hw[cnt == 4]]
        fc[nm] = {"cells": len(hs), "min": float(np.min(hs)), "median": float(np.median(hs)),
                  "max": float(np.max(hs))} if hs else None
    rep = {"out": str(out), "cells": int(len(Hw)), "points": int(len(Pw)), "watertight": bool(watertight),
           "edges": len(edges), "faces": len(faces), "edge_reprojection_max": erep["max_proj"],
           "face_correction_max": frep, "winslow": wrep, "quality": qs, "first_cell_height": fc,
           "classes": [{"class": ci, "axes": len(classes[ci]), **cc[ci]} for ci in cc],
           "gate": {"watertight": bool(watertight), "inverted_0": qs["inverted"] == 0}}
    import meshio
    meshio.write(str(out), meshio.Mesh(Pw, [("hexahedron", Hw)],
                                       cell_data={"block_id": [np.array(bids)], "sj": [cl["sj"]],
                                                  "inverted": [cl["inverted"].astype(float)],
                                                  "concave": [cl["concave"].astype(float)]}))
    out.with_suffix(".json").write_text(json.dumps(rep, indent=1))
    # entities for B6 (the ring meshes against the same cutout faces)
    np.savez(bd / "b5_entities.npz",
             face_keys=np.array([np.array(sorted(k)) for k in faces], dtype=object),
             face_grids=np.array([v[0] for v in faces.values()], dtype=object),
             face_points=np.array([v[1] for v in faces.values()], dtype=object),
             face_surface=np.array([str(v[2]) for v in faces.values()]))
    print(json.dumps({k: rep[k] for k in ("cells", "watertight", "edge_reprojection_max", "face_correction_max",
                                          "quality", "first_cell_height", "gate")}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
