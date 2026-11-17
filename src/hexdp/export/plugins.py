"""Exporter plugins ``vtk`` (fine hex mesh with block ids), ``foam`` (OpenFOAM polyMesh) and ``msh`` (gmsh,
physical groups, e.g. for FEniCSx)."""

from __future__ import annotations

import os
import shutil
from pathlib import Path

from hexdp.contracts import ParamSpec
from hexdp.io.blocks_vtk import write_blocks_vtk
from hexdp.pipeline.legacy import run_module


class VtkExporter:
    """Fine hex mesh with a per-cell block id (legacy VTK)."""

    name = "vtk"
    params = {}
    formats = ("vtk",)

    def write(self, bs, path, params, ctx):
        p = Path(path).with_suffix(".vtk")
        src = bs.meta.get("source")
        if src and os.path.exists(src):
            shutil.copyfile(src, p)          # byte-identical to the stage output
        else:
            write_blocks_vtk(bs, p)
        return [p]


class FoamExporter:
    """OpenFOAM polyMesh with patches from the mesh's patch field (blade pipeline export_foam)."""

    name = "foam"
    params = {"check": ParamSpec(False, "run checkMesh afterwards (needs OpenFOAM in the environment)")}
    formats = ("openfoam",)

    def write(self, bs, path, params, ctx):
        src = bs.meta.get("source")
        if not (src and os.path.exists(src)):
            raise ValueError("foam export needs the stage's mesh file (with a patch field), e.g. from the b5 filler")
        case = Path(path)
        args = [src, "--case", case, "--json", case.with_suffix(".json")]
        if not params["check"]:
            args.append("--no-check")
        run_module(ctx, "hexdp.export.foam", args, "foam_export.log")
        return [case]


class MshExporter:
    """Gmsh .msh 4.1: fine hex mesh, a volume group (or one per block) and a surface group per geometry
    surface of the boundary (e.g. for FEniCSx dolfinx.io.gmshio)."""

    name = "msh"
    params = {"block_groups": ParamSpec(False, "one physical volume per block instead of one 'domain'")}
    formats = ("gmsh",)

    def write(self, bs, path, params, ctx):
        from hexdp.export.msh import write_msh

        return [write_msh(bs, path, block_groups=params["block_groups"])]
