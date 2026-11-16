"""Filler plugins ``tfi`` and ``b5``, discretiser ``keep``.

``tfi`` wraps fill/tfi.py: conforming cell counts per direction class (MILP), curved edges (Coons),
TFI fill of every block, labelled by the input surface. ``b5`` wraps the blade pipeline's fine mesh
(strategies/blade_ring/mesh_b5.py): every edge and face meshed once on the exact dtOO surfaces, then
Gordon-Hall TFI and Winslow smoothing; it needs a prepared blade run directory (B1-B4).
"""

from __future__ import annotations

import os
from pathlib import Path

from hexdp.contracts import ParamSpec
from hexdp.io.blocks_vtk import read_blocks_vtk, write_blocks_vtk
from hexdp.pipeline.legacy import run_module


class KeepDiscretiser:
    """Keep the cell counts the structure already has (the filler may solve its own)."""

    name = "keep"
    params = {}

    def apply(self, bs, geo, params, ctx):
        return bs


def _input_file(bs, ctx, name):
    src = bs.meta.get("source")
    if src and os.path.exists(src):
        return Path(src)
    p = ctx.path(name)
    write_blocks_vtk(bs, p)
    return p


class TfiFiller:
    """Conforming cell counts from a target size and TFI fill with curved edges, snapped to the input surface."""

    name = "tfi"
    params = {
        "target_h": ParamSpec(0.05, "target cell size", unit="m", lower=0),
        "require_lattices": ParamSpec(True, "fail if a block is not a lattice"),
        "solve_divisions": ParamSpec(True, "solve conforming per-class cell counts (MILP)"),
        "curved_edges": ParamSpec("coons", "block edge curves", choices=("coons", "off")),
        "input_vtk": ParamSpec("", "tet input mesh that labels the surface (default: from the strategy)", kind="path"),
    }

    def fill(self, bs, geo, params, ctx):
        src = _input_file(bs, ctx, "tfi_in.vtk")
        tet = params["input_vtk"] or bs.meta.get("tet_vtk")
        out = ctx.workspace / "tfi.vtk"
        args = [src, "--target-h", params["target_h"], "--out", out, "--curved-edges", params["curved_edges"]]
        if tet:
            args += ["--input-vtk", tet]
        if params["require_lattices"]:
            args.append("--require-lattices")
        if params["solve_divisions"]:
            args += ["--solve-divisions", "--apply-divisions"]
        run_module(ctx, "hexdp.fill.tfi", args, "tfi.log")
        res = read_blocks_vtk(out)
        res.meta.update({k: v for k, v in bs.meta.items() if k not in ("source", "block_field")})
        return res


class B5Filler:
    """Fine mesh on the exact dtOO surfaces (blade pipeline B5, with --ring B6), every entity once."""

    name = "b5"
    params = {
        "run_dir": ParamSpec("", "prepared blade run directory (blade/blocks_b3.vtk, edges_b4.npz, ...); required", kind="path"),
        "target_h": ParamSpec(0.05, "target cell size", unit="m", lower=0),
        "first": ParamSpec(0.0009, "first cell height at hub / shroud", unit="m", lower=0),
        "growth": ParamSpec(1.2, "boundary-layer growth ratio", lower=1),
        "winslow": ParamSpec(300, "Winslow iterations per block (0 = off)", lower=0),
        "ring": ParamSpec(False, "mesh the blade ring blocks too (B6)"),
    }

    NEEDS = ("blade/blocks_b3.vtk", "blade/b3_ring.json", "blade/edges_b4.npz", "blade/machine_scratch/machine_blade.npz")

    def fill(self, bs, geo, params, ctx):
        run = Path(params["run_dir"] or "")
        missing = [n for n in self.NEEDS if not (run / n).exists()]
        if not params["run_dir"] or missing:
            raise ValueError(f"b5 needs a prepared blade run directory; missing in '{run}': {missing}")
        out = ctx.workspace / ("mesh_b6.vtk" if params["ring"] else "mesh_b5_core.vtk")
        args = [run, "--target-h", params["target_h"], "--first", params["first"], "--growth", params["growth"],
                "--winslow", params["winslow"], "--out", out]
        if params["ring"]:
            args.append("--ring")
        run_module(ctx, "hexdp.strategies.blade_ring.mesh_b5", args, out.stem + ".log")
        res = read_blocks_vtk(out)
        res.meta.update({k: v for k, v in bs.meta.items() if k not in ("source", "block_field")})
        return res
