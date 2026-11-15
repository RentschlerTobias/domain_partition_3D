"""Strategy plugin ``algohex``: tet preparation -> AlgoHex field-aligned hex mesh -> beam sheet collapse.

The stages are the keep-BL datagen chain (datagen/run_sample_keepbl.sh), run as subprocesses in the
workspace. A stage whose output already exists in the workspace is reused, so a run can be resumed
(AlgoHex alone takes minutes to hours).
"""

from __future__ import annotations

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
        "hausdorff_tol": ParamSpec(0.1, "max Hausdorff distance raw hex mesh -> input surface", lower=0),
        "beam_workers": ParamSpec(4, "worker processes of the beam search", lower=1),
    }

    def applicable(self, geo, features):
        if geo.source is not Source.MSH or not geo.path:
            return Score.no("tet_prep reads a gmsh .msh (the dtOO hybrid mesh)", missing=("dtOO .msh input",))
        return Score(0.4, ["general topology: no sweep direction or symmetry needed",
                           "input is a .msh with entities (tet_prep v5 cuts the blade O-grid)"])

    def run(self, geo, features, params, ctx):
        ws = ctx.workspace
        tet, ovm, blocks = ws / "tet.vtk", ws / "hex_hex_x.ovm", ws / "blocks.vtk"
        if not tet.exists():
            args = ["--msh", geo.path, "--out-dir", ws, "--out", tet.name, "--size-max", params["size_max"]]
            if params["keep_prisms"]:
                args.append("--keep-prisms")
            args += ["--remesh-all", params["remesh_all"]]
            run_module(ctx, "hexdp.strategies.algohex.tet_prep", args, "tet_prep.log")
        else:
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
                        "--hausdorff-tol", params["hausdorff_tol"], "--workers", params["beam_workers"]],
                       "beam.log")
        bs = read_blocks_vtk(blocks)
        bs.meta["tet_vtk"] = str(tet)
        bs.meta["msh"] = geo.path
        return bs

