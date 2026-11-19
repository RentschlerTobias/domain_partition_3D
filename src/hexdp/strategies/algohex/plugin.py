"""Strategy plugin ``algohex``: tet preparation -> AlgoHex field-aligned hex mesh -> beam sheet collapse.

The stages are the keep-BL datagen chain (datagen/run_sample_keepbl.sh), run as subprocesses in the
workspace. A stage whose output already exists in the workspace is reused, so a run can be resumed
(AlgoHex alone takes minutes to hours).
"""

from __future__ import annotations

import shutil

from hexdp.contracts import ParamSpec, Score
from hexdp.core.geometry import Source
from hexdp.io.blocks_vtk import read_blocks_vtk
from hexdp.pipeline.legacy import run_module


class AlgoHexStrategy:
    """Unstructured field-aligned blocking: AlgoHex frame field and integer-grid map, then the beam search
    collapses separatrix sheets to the fewest blocks that stay valid and close to the input surface."""

    name = "algohex"
    params = {
        "size_max": ParamSpec(0.12, "maximum tet size of the AlgoHex input", unit="m", lower=0),
        "keep_prisms": ParamSpec(True, "keep the hub/shroud boundary layer in the AlgoHex domain"),
        "remesh_all": ParamSpec(0.0, "remesh the whole input at this size (0 = no remesh)", unit="m", lower=0),
        "n": ParamSpec(2000, "AlgoHex sizing parameter -n", lower=1),
        "threads": ParamSpec(0, "AlgoHex/MUMPS threads (0 = the job's core count)", lower=0),
        "hexmeshing_bin": ParamSpec("", "HexMeshing binary (default: $HEXMESHING_BIN)"),
        "beam_width": ParamSpec(32, "beam width of the sheet collapse", lower=1),
        "beam_depth": ParamSpec(20, "beam depth (max collapses)", lower=1),
        "hausdorff_tol": ParamSpec(0.0, "absolute Hausdorff tolerance of the collapse (0 = relative, hausdorff_rel)",
                                   lower=0),
        "hausdorff_rel": ParamSpec(0.02, "Hausdorff tolerance as a fraction of the bounding-box diagonal", lower=0),
        "max_raw_rel": ParamSpec(0.05, "reject the raw hex mesh above this Hausdorff distance (fraction of the "
                                 "diagonal; boundary vertices are snapped back onto the surface below it)", lower=0),
        "beam_workers": ParamSpec(4, "worker processes of the beam search", lower=1),
    }

    def applicable(self, geo, features):
        from hexdp.backends import BACKENDS, available
        if not available("algohex"):
            return Score.no("AlgoHex (HexMeshing) is not installed - " + BACKENDS["algohex"]["facts"],
                            missing=("backend algohex",))
        if geo.path and str(geo.path).lower().endswith(".vtk"):
            return Score(0.4, ["general topology: no sweep direction or symmetry needed",
                               "input is already an AlgoHex tet mesh with feature tags (no tet_prep)"])
        if geo.source is not Source.MSH or not geo.path:
            return Score.no("tet_prep reads a gmsh .msh (the dtOO hybrid mesh)", missing=("dtOO .msh input",))
        return Score(0.4, ["general topology: no sweep direction or symmetry needed",
                           "input is a .msh with entities (tet_prep v5 cuts the blade O-grid)"])

    def _warn_ignored(self, params, ctx):
        """size_max and remesh_all act only in tet_prep; with a reused tet.vtk they change nothing."""
        ignored = [f"{k}={params[k]}" for k in ("size_max", "remesh_all")
                   if params[k] != self.params[k].default]
        if not ignored:
            return
        msg = (f"{', '.join(ignored)}: size_max/remesh_all act only in tet_prep; tet.vtk is reused here "
               "(an already tagged VTK tet input, or a resumed run), so these parameters have no effect")
        print(f"[algohex] WARNING: {msg}", flush=True)
        ctx.decide("algohex", "WARNING: size_max/remesh_all ignored", [msg])

    def run(self, geo, features, params, ctx):
        ws = ctx.workspace
        tet, ovm, blocks = ws / "tet.vtk", ws / "hex_hex_x.ovm", ws / "blocks.vtk"
        if not tet.exists() and str(geo.path).lower().endswith(".vtk"):
            shutil.copyfile(geo.path, tet)               # an AlgoHex input already (HexMe): used as it is
            ctx.decide("algohex", "input is an AlgoHex tet mesh", [str(geo.path)])
        if not tet.exists():
            args = ["--msh", geo.path, "--out-dir", ws, "--out", tet.name, "--size-max", params["size_max"]]
            if params["keep_prisms"]:
                args.append("--keep-prisms")
            args += ["--remesh-all", params["remesh_all"]]
            run_module(ctx, "hexdp.strategies.algohex.tet_prep", args, "tet_prep.log")
        else:
            self._warn_ignored(params, ctx)
            ctx.decide("algohex", "reuse tet.vtk", [str(tet)])
        if not ovm.exists():
            th = params["threads"] or ctx.cores
            env = {"HEXMESHING_BIN": params["hexmeshing_bin"]} if params["hexmeshing_bin"] else {}
            run_module(ctx, "hexdp.strategies.algohex.run_algohex",
                       ["--backend", "native", "--prefix", "hex", "--tag", "x", "--in-vtk", tet, "--out-dir", ws,
                        "--", "-n", params["n"]], "algohex.log", threads=th, env=env)
        else:
            ctx.decide("algohex", "reuse hex_hex_x.ovm", [str(ovm)])
        if not blocks.exists():
            run_module(ctx, "hexdp.ops.beam_collapse",
                       [ovm, "--input-vtk", tet, "--out", blocks, "--width", params["beam_width"],
                        "--depth", params["beam_depth"], "--guard", "struct",
                        *(["--hausdorff-tol", params["hausdorff_tol"]] if params["hausdorff_tol"] > 0 else
                          ["--relative", "--hausdorff-tol", params["hausdorff_rel"],
                           "--max-raw-hausdorff", params["max_raw_rel"]]),
                        "--workers", params["beam_workers"]],
                       "beam.log")
        bs = read_blocks_vtk(blocks)
        bs.meta["tet_vtk"] = str(tet)
        bs.meta["msh"] = geo.path
        return bs

