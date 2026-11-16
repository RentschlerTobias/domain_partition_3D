"""Strategy plugins ``sweep`` and ``cgrid``: 2.5D phi sweep of a meridional (r-z) block layout.

``sweep``: r-z sections of the CAD domain at phi = 0, step, ..., phi_max (section_rz, pythonOCC) ->
2D layout with greedy chord collapse (collapse2d) -> sweep of every 2D block over [0, 45] and [45, 90]
deg plus a quarter butterfly on the axis, TFI fill (sweep3d).

``cgrid``: as ``sweep`` plus C-layers around thin walls (baffles): wall detection and sections in the
wall CAD (wall_faces, pythonOCC) -> C-grid layout and sweep (sweep3d_c).

Both sweep about the z axis over a 0..90 deg sector (the 45 deg split of the butterfly); applicable()
checks that with the result of the ``sweep`` detector. Stage outputs in the workspace are reused.
"""

from __future__ import annotations

from hexdp.contracts import ParamSpec, Score
from hexdp.core.geometry import Source
from hexdp.io.blocks_vtk import read_blocks_vtk
from hexdp.pipeline.legacy import run_module

SWEEP_PARAMS = {
    "phi_step": ParamSpec(2.5, "angle between r-z sections", unit="deg", lower=0.1),
    "phi_max": ParamSpec(90.0, "sector angle (the layout assumes a 45 deg mid section)", unit="deg", choices=(90.0,)),
    "h2d": ParamSpec(0.1, "2D cell size of the section mesh", unit="m", lower=0),
    "a_core": ParamSpec(0.45, "relative size of the butterfly core on the axis", lower=0.05, upper=0.95),
    "collapse_guard": ParamSpec(0.1, "min scaled Jacobian a chord collapse must keep", lower=0, upper=1),
    "collapse_h": ParamSpec(0.15, "cell size of the collapse guard re-fit", unit="m", lower=0),
    "workers": ParamSpec(8, "worker processes", lower=1),
}


def _sweep_score(geo, features, name):
    import numpy as np

    rej = next((f for f in features if f.kind == "sweep_rejected"), None)
    if rej is not None:                     # the geometric reason first: it holds for any source
        return Score.no("sweep test failed: " + "; ".join(rej.reasons[1:]))
    if geo.source is not Source.CAD:
        return Score.no(f"{name} slices exact CAD faces (STEP/BREP)", missing=("CAD domain solid",))
    sw = next((f for f in features if f.kind == "sweep"), None)
    if sw is None:
        return Score(0.1, ["no sweep test result: assumes a rotational 0-90 deg sector"],
                     missing=["sweep feature (add the 'sweep' detector)"])
    a, p = np.array(sw.data["axis_direction"]), np.array(sw.data["axis_point"])
    scale = float(np.linalg.norm(np.subtract(*geo.bbox())))
    if abs(abs(a[2]) - 1) > 1e-3 or np.linalg.norm(p[:2]) > 1e-3 * scale:
        return Score.no(f"{name} sweeps about the z axis through the origin; the axis is {np.round(a, 3).tolist()} "
                        f"through {np.round(p, 3).tolist()}", missing=("transform to the z axis",))
    if abs(sw.data["angular_span_deg"] - 90.0) > 0.5:
        return Score.no(f"{name} needs a 90 deg sector (45 deg mid section), the domain spans "
                        f"{sw.data['angular_span_deg']:.1f} deg")
    return Score(round(0.8 * sw.confidence + 0.1, 3), ["sweep test passed"] + sw.reasons)


def prepare_sections_and_layout(geo, params, ctx):
    """sections/ and layout_final.json in the workspace (reused when present)."""
    ws = ctx.workspace
    if not (ws / "sections").is_dir() or not any((ws / "sections").glob("sec_*.json")):
        n = int(round(params["phi_max"] / params["phi_step"]))
        phis = [f"{params['phi_step'] * i:g}" for i in range(n + 1)]
        run_module(ctx, "hexdp.strategies.sweep.section_rz", ["--brep", geo.path, *phis], "section_rz.log", occ=True)
    else:
        ctx.decide("sweep", "reuse sections/", [str(ws / "sections")])
    if not (ws / "layout_final.json").exists():
        run_module(ctx, "hexdp.strategies.sweep.collapse2d",
                   [params["collapse_guard"], params["collapse_h"], "--dir", ws, "--workers", params["workers"]],
                   "collapse2d.log")
    else:
        ctx.decide("sweep", "reuse layout_final.json", [str(ws / "layout_final.json")])


class SweepStrategy:
    """2.5D phi sweep of a collapsed r-z section layout, quarter butterfly on the axis, TFI fill."""

    name = "sweep"
    params = dict(SWEEP_PARAMS)

    def applicable(self, geo, features):
        return _sweep_score(geo, features, self.name)

    def run(self, geo, features, params, ctx):
        prepare_sections_and_layout(geo, params, ctx)
        out = ctx.workspace / "sweep_tfi_blocks.vtk"
        if not out.exists():
            run_module(ctx, "hexdp.strategies.sweep.sweep3d", [params["h2d"], params["a_core"], "--dir", ctx.workspace],
                       "sweep3d.log")
        bs = read_blocks_vtk(out)
        bs.meta["cad"] = geo.path
        return bs


class CGridStrategy:
    """Sweep with C-layers (graded wall layer) around thin walls, which stay inside the domain as baffles."""

    name = "cgrid"
    params = dict(SWEEP_PARAMS, wall_cad=ParamSpec("", "CAD with the thin walls (STEP); required", kind="path"))

    def applicable(self, geo, features):
        s = _sweep_score(geo, features, self.name)
        if not s.ok:
            return s
        walls = next((f for f in features if f.kind == "thin_wall"), None)
        if walls is None:
            return Score(s.value * 0.9, s.reasons + ["needs the wall CAD as parameter wall_cad"],
                         s.missing + ["thin_wall feature (detector, plan step 3)"])
        return Score(min(1.0, s.value + 0.1), s.reasons + walls.reasons, s.missing)

    def run(self, geo, features, params, ctx):
        if not params["wall_cad"]:
            raise ValueError("cgrid needs parameter wall_cad (the STEP with the thin walls)")
        prepare_sections_and_layout(geo, params, ctx)
        cg = ctx.workspace / "cgrid"
        cg.mkdir(exist_ok=True)
        if not (cg / "walls.json").exists():
            run_module(ctx, "hexdp.strategies.cgrid.wall_faces",
                       ["--step", params["wall_cad"], "--brep", geo.path, "--out", cg], "wall_faces.log", cwd=cg, occ=True)
        out = cg / "sweep_tfi_blocks.vtk"
        if not out.exists():
            run_module(ctx, "hexdp.strategies.cgrid.sweep3d_c",
                       [params["h2d"], params["a_core"], "--cgrid-dir", cg, "--sweep-dir", ctx.workspace,
                        "--workers", params["workers"]], "sweep3d_c.log", cwd=cg)
        bs = read_blocks_vtk(out)
        bs.meta["cad"] = geo.path
        return bs
