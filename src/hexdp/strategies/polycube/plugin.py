"""Strategy plugin ``polycube``: polycube labeling (evocube, external GPL-3 tool, called as a program) -> HexEx hex
mesh of the polycube map -> beam sheet collapse to blocks (as for AlgoHex).

Input: a CAD solid (STEP/BREP/IGES: tet-meshed with gmsh here) or a tet mesh in the AlgoHex VTK format (e.g. HexMe).
The evocube programs (prep_cli, evolabel_cli, polycube_withHexEx; built locally in /opt/stack/ext/evocube) are found in
``evocube_bin`` or $HEXDP_EVOCUBE_BIN. Stage outputs in the workspace are reused.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import numpy as np

from hexdp.contracts import ParamSpec, Score
from hexdp.core.geometry import Source
from hexdp.io.blocks_vtk import read_blocks_vtk
from hexdp.pipeline.legacy import StageError, run_module

TOOLS = ("prep_cli", "evolabel_cli", "polycube_withHexEx")


def _bin_dir(param):
    d = Path(param or os.environ.get("HEXDP_EVOCUBE_BIN", "/opt/stack/ext/evocube/build"))
    missing = [t for t in TOOLS if not (d / t).exists()]
    if missing:
        raise StageError(f"evocube programs {missing} not found in {d} (set evocube_bin or $HEXDP_EVOCUBE_BIN)")
    return d


def _label_note(ws):
    """Which of the six cube directions the evocube labeling uses (for the error message, not a check).

    A label is ``2*axis + dir`` over the three axes (evocube ``src/flagging_utils.cpp:axesMatrix``), so the six
    valid values are 0..5 - 0 is a direction like any other. Only ``labeling_on_tets.txt`` uses -1 ("no label").
    """
    vals = set()
    try:
        for line in open(Path(ws) / "labeling.txt"):
            tok = line.split()
            if len(tok) == 1:            # the file is one label per line
                try:
                    vals.add(int(tok[0]))
                except ValueError:
                    pass
    except OSError:
        return ""
    vals.discard(-1)
    if not vals:
        return ""
    used = sorted(vals & set(range(6)))
    bad = sorted(vals - set(range(6)))
    note = f" (labeling uses {len(used)} of 6 directions: {used}"
    return note + (f", invalid labels {bad})" if bad else ")")


def _lib_dirs(bindir):
    """The library directories shipped with the evocube build: ``lib`` (or ``lib64``) and every
    subdirectory that holds shared libraries (e.g. ``lib/gco``, ``lib/OpenVolumeMesh/src``).
    The CLIs need them on LD_LIBRARY_PATH."""
    roots = [Path(bindir), *([Path(bindir).parent] if Path(bindir).name == "bin" else [])]
    out = []
    for root in roots:
        for name in ("lib", "lib64"):
            lib = root / name
            if not lib.is_dir():
                continue
            for d in sorted({p.parent for p in lib.rglob("*.so*") if p.is_file()}):
                if d not in out:
                    out.append(d)
    return out


def _run(ctx, cmd, log, threads, libs=()):
    full = [str(c) for c in cmd]
    if ctx.cpu_list and shutil.which("taskset"):
        full = ["taskset", "-c", ctx.cpu_list, *full]
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    if libs:
        old = env.get("LD_LIBRARY_PATH")
        env["LD_LIBRARY_PATH"] = os.pathsep.join([*(str(p) for p in libs), *([old] if old else [])])
    with open(ctx.path(log), "w") as fh:
        try:
            r = subprocess.run(full, cwd=str(ctx.workspace), stdout=fh, stderr=subprocess.STDOUT, env=env,
                               timeout=ctx.timeout_s)
        except subprocess.TimeoutExpired as exc:
            raise StageError(f"{cmd[0]}: timeout after {ctx.timeout_s} s") from exc
    if r.returncode != 0:
        tail = "".join(open(ctx.path(log)).readlines()[-12:])
        raise StageError(f"{Path(cmd[0]).name} failed with exit code {r.returncode}:\n{tail}")


def tet_inputs(geo, ws, size):
    """tetra.mesh (evocube) and tet.vtk (AlgoHex format: surface ids for the beam's labeller) of the geometry."""
    import gmsh

    from hexdp.strategies.algohex.tet_prep_v1 import write_algohex_vtk
    if str(geo.path).lower().endswith(".vtk"):
        import meshio
        shutil.copyfile(geo.path, ws / "tet.vtk")
        m = meshio.read(geo.path)
        T = np.vstack([c.data for c in m.cells if c.type == "tetra"])
        meshio.write(ws / "tetra.mesh", meshio.Mesh(m.points, [("tetra", T)]), file_format="medit")
        return
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.occ.importShapes(str(geo.path))
        gmsh.model.occ.synchronize()
        lo, hi = geo.bbox()
        gmsh.option.setNumber("Mesh.MeshSizeMax", size or float(np.linalg.norm(hi - lo)) / 30)
        gmsh.model.mesh.generate(3)
        gmsh.write(str(ws / "tetra.mesh"))
        write_algohex_vtk(gmsh, str(ws / "tet.vtk"))
    finally:
        gmsh.finalize()


class PolycubeStrategy:
    """Polycube blocking: evocube labels the boundary with +-x/y/z charts (genetic search, valid polycube), HexEx
    extracts the hex mesh of the polycube map, the beam collapses its sheets to few blocks. Good for compact genus-0
    CAD parts (L 3, T 4, cylinder 1 block); a part with a through hole fails for now (the polycube mesh has no interior
    singularity, a block wraps around the hole); weak for thin twisted shapes (runner blades)."""

    name = "polycube"
    params = {
        "evocube_bin": ParamSpec("", "directory of the evocube programs (default $HEXDP_EVOCUBE_BIN)"),
        "size": ParamSpec(0.0, "tet size for a CAD input (0 = bbox diagonal / 30)", unit="m", lower=0),
        "hex_scale": ParamSpec(1.0, "HexEx scale of the polycube map (larger = more hexes)", lower=0.01),
        "threads": ParamSpec(0, "threads of the evocube stages (0 = the job's core count)", lower=0),
        "label_threads": ParamSpec(1, "threads of the labeling search; 1 is deterministic (evocube shares one "
                                      "rand() across the OpenMP mutations, so every other count gives a different "
                                      "labeling per run) and costs ~2 s here", lower=1),
        "min_sj": ParamSpec(0.05, "reject a block structure whose worst fine cell is below this scaled Jacobian "
                                  "(the structural check alone accepts nearly degenerate cells)", lower=0),
        "beam_width": ParamSpec(16, "beam width of the sheet collapse", lower=1),
        "beam_depth": ParamSpec(16, "beam depth", lower=1),
        "hausdorff_tol": ParamSpec(0.05, "max boundary deviation of the collapse as a fraction of the bounding-box "
                                         "diagonal (relative, see beam)", lower=0),
        "beam_workers": ParamSpec(4, "worker processes of the beam search", lower=1),
    }

    def applicable(self, geo, features):
        from hexdp.backends import available
        if not available("evocube"):
            return Score.no("evocube is not installed (polycube labeling; see docs/backends.md)",
                            missing=("backend evocube",))
        p = str(geo.path or "").lower()
        if geo.surfaces:
            topo = geo.boundary_topology()
            if topo["closed"] and topo["components"] > 1:
                return Score.no(f"the boundary has {topo['components']} components (a hollow domain, e.g. a "
                                "box-embedded part): a polycube labeling needs one closed boundary")
        if p.endswith(".vtk") or geo.source is Source.CAD:
            return Score(0.35, ["general topology without a frame field: a valid polycube map by construction",
                                "suits compact parts; few blocks where the shape is close to a polycube"])
        return Score.no("needs a CAD solid or an AlgoHex-format tet mesh (.vtk)", missing=("CAD solid or tet .vtk",))

    def run(self, geo, features, params, ctx):
        ws = ctx.workspace
        bindir = _bin_dir(params["evocube_bin"])
        libs = _lib_dirs(bindir)
        th = params["threads"] or ctx.cores
        if not (ws / "tet.vtk").exists():
            tet_inputs(geo, ws, params["size"])
        if not (ws / "labeling_on_tets.txt").exists():
            _run(ctx, [bindir / "prep_cli", "."], "evocube_prep.log", th, libs)
            _run(ctx, [bindir / "evolabel_cli", "boundary.obj", "."], "evocube_label.log", params["label_threads"],
                 libs)
        if not (ws / "hexes.vtk").exists():
            _run(ctx, [bindir / "polycube_withHexEx", ".", params["hex_scale"]], "evocube_hexex.log", th, libs)
            import meshio

            from hexdp.io.ovm import write_hex_vtk
            m = meshio.read(ws / "hexes.mesh")
            _cells = [c.data for c in m.cells if c.type == "hexahedron" and len(c.data)]
            if not _cells:
                raise RuntimeError(
                    f"polycube: evocube extracted no hexahedra{_label_note(ws)}. The polycube map is degenerate, "
                    f"so HexEx has nothing to extract; use the algohex strategy for this geometry "
                    f"(logs: {ws / 'evocube_label.log'}).")
            H = np.vstack(_cells)
            write_hex_vtk(str(ws / "hexes.vtk"), m.points, H)
            ctx.decide("polycube", f"{len(H)} hexes from the polycube map", [str(ws / "hexes.vtk")])
        if not (ws / "blocks.vtk").exists():
            run_module(ctx, "hexdp.ops.beam_collapse",
                       [ws / "hexes.vtk", "--input-vtk", ws / "tet.vtk", "--out", ws / "blocks.vtk",
                        "--width", params["beam_width"], "--depth", params["beam_depth"], "--guard", "struct",
                        "--relative",
                        "--hausdorff-tol", params["hausdorff_tol"], "--workers", params["beam_workers"],
                        "--min-sj", params["min_sj"]],
                       "beam.log")
        bs = read_blocks_vtk(ws / "blocks.vtk")
        bs.meta["source"] = str(ws / "blocks.vtk")
        return bs
