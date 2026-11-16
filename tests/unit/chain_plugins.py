"""Minimal plugins of every kind, enough to run the whole chain on an axis-aligned box. They show the
contracts in use; the real plugins arrive with the migration (step 1.4)."""

import math
import os

import numpy as np

from hexdp.contracts import Feature, ParamSpec, Score
from hexdp.core.blockstructure import OFF, PARALLEL, BlockStructure, FaceInfo
from hexdp.io.blocks_vtk import from_grids, write_blocks_vtk


class BBoxDetector:
    """Reports the bounding box and whether the boundary is closed."""
    name = "bbox"
    params = {}

    def detect(self, geo, params, ctx):
        lo, hi = geo.bbox()
        closed = all(len(c.surfaces) >= 2 or c.kind.value == "seam" for c in geo.curves.values())
        return [Feature("bbox", {"lo": lo.tolist(), "hi": hi.tolist(), "closed": closed},
                        reasons=[f"extent {np.round(hi - lo, 3).tolist()}", "closed" if closed else "open"])]


class BoxStrategy:
    """One block spanning the bounding box; only for boxes."""
    name = "box"
    params = {"flat_tol": ParamSpec(1e-6, "max distance of a box corner to the geometry", unit="m", lower=0),
              "ref_file": ParamSpec("", "any file; must arrive as an existing absolute path", kind="path")}

    def applicable(self, geo, features):
        bb = next((f for f in features if f.kind == "bbox"), None)
        if bb is None:
            return Score.no("no bounding box", missing=("bbox detector",))
        if len(geo.surfaces) != 6:
            return Score.no(f"{len(geo.surfaces)} surfaces, a box has 6")
        return Score(0.5, ["six surfaces, closed boundary: a single block fits"])

    def run(self, geo, features, params, ctx):
        f = params.get("ref_file")
        assert not f or (os.path.isabs(f) and os.path.exists(f)), f
        lo, hi = (np.array(v) for v in (features[0].data["lo"], features[0].data["hi"]))
        P = np.array([np.where(o, hi, lo) for o in OFF])
        bs = BlockStructure(P, np.arange(8)[None])
        for key in bs.faces():
            c = P[sorted(key)].mean(0)
            _, _, tag = geo.project(c[None])
            s = geo.surfaces[int(tag[0])]
            bs.face_info[key] = FaceInfo(s.role, s.tag)
        return bs


class NeverStrategy:
    """Needs a blade; never applies to a box."""
    name = "never"
    params = {}

    def applicable(self, geo, features):
        return Score.no("no blade found", missing=("blade feature",))

    def run(self, geo, features, params, ctx):
        raise AssertionError("must not run")


class SplitOp:
    """Split every block in two along its i axis (straight blocks only)."""
    name = "split_i"
    params = {}

    def applicable(self, bs, geo):
        return Score(1.0, ["test op"]) if bs.n_blocks == 1 else Score.no("only splits a single block")

    def apply(self, bs, geo, params, ctx):
        X = bs.points[bs.blocks[0]]
        mid = 0.5 * (X[[0, 3, 7, 4]] + X[[1, 2, 6, 5]])
        P = np.vstack([X, mid])          # 0..7 old corners, 8..11 mids of edges 0-1, 3-2, 7-6, 4-5
        a = [0, 8, 9, 3, 4, 11, 10, 7]
        b = [8, 1, 2, 9, 11, 5, 6, 10]
        return BlockStructure(P, np.array([a, b]))


class UniformDiscretiser:
    """Cells per sheet from a target size: n = ceil(longest edge of the sheet / target_h)."""
    name = "uniform"
    params = {"target_h": ParamSpec(0.25, "target cell size", unit="m", lower=1e-9)}

    def apply(self, bs, geo, params, ctx):
        for sh in bs.sheets():
            L = max(np.linalg.norm(np.subtract(*bs.points[sorted(e)])) for e in sh["edges"])
            n = max(1, math.ceil(L / params["target_h"] - 1e-9))
            for e in sh["edges"]:
                bs.edge_cells[e] = n
        return bs


class LinearFiller:
    """Trilinear fine grids from the 8 corners of each block; shared nodes merged."""
    name = "linear"
    params = {}

    def fill(self, bs, geo, params, ctx):
        P, index, grids = [], {}, []
        for c8 in bs.blocks:
            X = bs.points[c8]
            n = [bs.edge_cells[frozenset((int(c8[i]), int(c8[j])))] for i, j in (p[0] for p in PARALLEL)]
            u, v, w = (np.linspace(0, 1, k + 1) for k in n)
            U, V, W = np.meshgrid(u, v, w, indexing="ij")
            G = np.empty(U.shape, np.int64)
            for idx in np.ndindex(U.shape):
                t = np.array([U[idx], V[idx], W[idx]])
                wts = np.prod(np.where(OFF, t, 1 - t), axis=1)
                p = wts @ X
                key = tuple(np.round(p, 10))
                G[idx] = index.setdefault(key, len(P))
                if G[idx] == len(P):
                    P.append(p)
            grids.append(G)
        out = from_grids(np.array(P), grids)
        out.face_info = bs.face_info
        return out


class VtkExporter:
    """Fine hex mesh with block ids."""
    name = "vtk"
    params = {}
    formats = ("vtk",)

    def write(self, bs, path, params, ctx):
        p = path.with_suffix(".vtk")
        write_blocks_vtk(bs, p)
        return [p]
