"""Op plugin ``block_collapse``: collapse whole block sheets of a block structure, guarded by the CAD.

Wraps ops/block_collapse.py (canada blocking study): CAD classification of nodes, edges and faces
(cad_export, pythonOCC), topological sheet collapse with closure-relation checks, TFI re-fit, quality
and Hausdorff guards, beam search over collapse sequences.
"""

from __future__ import annotations

import os

from hexdp.contracts import ParamSpec, Score
from hexdp.io.blocks_vtk import read_blocks_vtk, write_blocks_vtk
from hexdp.pipeline.legacy import run_module


class BlockCollapseOp:
    """Remove block sheets while every CAD entity keeps its blocks and quality / Hausdorff stay in bounds."""

    name = "block_collapse"
    params = {
        "h": ParamSpec(0.0, "TFI fill cell size (<= 0: keep the input resolution)", unit="m"),
        "tol": ParamSpec(0.01, "CAD classification distance tolerance", unit="m", lower=0),
        "sj_floor": ParamSpec("0.05", "min scaled Jacobian: a number, or 'input' = not below the input"),
        "htol": ParamSpec(0.01, "allowed Hausdorff increase over the baseline re-fit", unit="m", lower=0),
        "width": ParamSpec(4, "beam width", lower=1),
        "depth": ParamSpec(8, "beam depth", lower=1),
        "workers": ParamSpec(8, "worker processes", lower=1),
        "cad": ParamSpec("", "CAD solid for the classification (default: the geometry's file)"),
    }

    def applicable(self, bs, geo):
        if not bs.grids:
            return Score.no("needs the fine block grids (a filled block structure)", missing=("fine grids",))
        if not (getattr(geo, "exact", False) or bs.meta.get("cad")):
            return Score.no("classifies against exact CAD faces", missing=("CAD geometry",))
        return Score(0.6, [f"{bs.n_blocks} blocks with fine grids; CAD available for the guards"])

    def apply(self, bs, geo, params, ctx):
        ws = ctx.workspace
        cad_src = params["cad"] or bs.meta.get("cad") or geo.path
        cad = ws / "cad.npz"
        if not cad.exists():
            run_module(ctx, "hexdp.io.cad_export", [cad_src, cad], "cad_export.log", occ=True)
        src = bs.meta.get("source")
        if not (src and os.path.exists(src)):           # write the structure only if it has no file yet
            src = ws / "blocks_in.vtk"
            write_blocks_vtk(bs, src)
        out = ws / "block_collapse"
        best = out / f"best_h{float(params['h'])}_tfi.vtk"   # block_collapse names it after str(h)
        if not best.exists():
            run_module(ctx, "hexdp.ops.block_collapse",
                       ["--vtk", src, "--cad", cad, "--out", out, "--field", bs.meta.get("block_field", "block"),
                        "--h", params["h"], "--tol", params["tol"], "--sj-floor", params["sj_floor"],
                        "--htol", params["htol"], "--width", params["width"], "--depth", params["depth"],
                        "--workers", params["workers"]], "block_collapse.log")
        res = read_blocks_vtk(best)
        res.meta.update({k: v for k, v in bs.meta.items() if k not in ("source", "block_field")})
        res.meta["source"] = str(best)
        return res
