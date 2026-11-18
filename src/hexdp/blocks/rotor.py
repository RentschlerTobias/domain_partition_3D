"""Building block ``dtoo_rotor``: the input of a bladed rotor from a parametric dtOO design (the case geometry is
the parameter file). Runs the dtOO stage (dtoo_rotor.py, dtOO environment) into rotor.npz in the workspace (reused)
and returns the rotor description for band_layout and section_sweep:

    r, z0, z1, symmetry (blades), loops = [(footprint on the unrolled hub, pins)] x blades,
    stages = [fillet patches, blade patches] (or [blade] with the root extended to the hub when fillet is off)
"""

from __future__ import annotations

import numpy as np

from hexdp.blocks.geom import area2d
from hexdp.contracts import ParamSpec
from hexdp.pipeline.legacy import run_module


def footprint(patches, r_hub):
    """Loop on the unrolled hub (x = r phi, z) from the first stage's patches at w = 0, clockwise (a hole for the
    band orientation), the patch joints as pins."""
    curves = [np.asarray(f)[:, 0] for f in patches]
    pts3 = np.vstack([c[:-1] for c in curves])
    pts = np.c_[r_hub * np.unwrap(np.arctan2(pts3[:, 1], pts3[:, 0])), pts3[:, 2]]
    pins = np.cumsum([0] + [len(c) - 1 for c in curves[:-1]])
    if area2d(pts) > 0:
        pts = pts[::-1]
        pins = sorted((len(pts) - p) % len(pts) for p in pins)
    return pts, [int(p) for p in pins]


def extend_root(blade, r_hub, nw):
    """Patches from the hub to the blade root: the root profile moved radially (same phi and z) - replaces a
    fillet that runs out tangentially into the hub."""
    out = []
    for B in blade:
        root = np.asarray(B, float)[:, 0]
        r_root = np.hypot(root[:, 0], root[:, 1])
        F = np.empty((len(root), nw, 3))
        for j, t in enumerate(np.linspace(0, 1, nw)):
            r = r_hub + t * (r_root - r_hub)
            F[:, j, 0] = root[:, 0] * r / r_root
            F[:, j, 1] = root[:, 1] * r / r_root
            F[:, j, 2] = root[:, 2]
        out.append(F)
    return out


def rotor_description(z, fillet=True):
    r = float(z["r_hub"])
    nb = int(z["nb"])
    blade = list(np.asarray(z["blade"]))
    first = list(np.asarray(z["fillet"])) if fillet else extend_root(blade, r, np.asarray(z["fillet"]).shape[2])
    pts, pins = footprint(first, r)
    P = 2 * np.pi * r
    loops = [(pts + [k * P / nb, 0.0], pins) for k in range(nb)]
    return dict(r=r, z0=float(z["z0"]), z1=float(z["z1"]), symmetry=nb, loops=loops, stages=[first, blade],
                fillet=bool(fillet))


class DtooRotor:
    """Bladed rotor input from a dtOO parameter file (the case geometry): hub, blade + fillet patches."""

    name = "dtoo_rotor"
    params = {
        "dtoo_case": ParamSpec("/dtOO/demo/tistos", "dtOO machine template directory", kind="path"),
        "fillet": ParamSpec(True, "keep the blade fillet (off: blade root extended radially to the hub)"),
        "prefix": ParamSpec("tV_ruWithRounding", "dtOO label prefix of the blade O-grid volumes"),
        "patches": ParamSpec("0_0,le_0,le_1,0_2,te_0,te_1", "patch names in cyclic order around the blade"),
        "nblades_const": ParamSpec("cV_ru_nBlades", "dtOO constant: number of blades"),
        "z_consts": ParamSpec("cV_ruMech_z_hub_0,cV_ruMech_z_hub_1", "dtOO constants: hub z range"),
        "adjust_plugin": ParamSpec("ru_adjustDomain", "dtOO plugin applied before the final build ('' = none)"),
        "nu": ParamSpec(201, "surface samples along a patch", lower=3),
        "nw": ParamSpec(61, "surface samples from the hub to the tip", lower=3),
    }

    def run(self, inputs, params, ctx):
        geo = inputs["geometry"]
        f = ctx.workspace / "rotor.npz"
        if not f.exists():
            args = [geo.path, "--case", params["dtoo_case"], "--out", f, "--prefix", params["prefix"],
                    "--patches", params["patches"], "--nblades-const", params["nblades_const"],
                    "--z-consts", params["z_consts"], "--adjust-plugin", params["adjust_plugin"],
                    "--nu", params["nu"], "--nw", params["nw"]]
            run_module(ctx, "hexdp.blocks.dtoo_rotor", args, "dtoo_rotor.log", occ=True,
                       env={"OSLO_LOCK_PATH": "/tmp", "FOAM_SIGFPE": "0"})
        else:
            ctx.decide("dtoo_rotor", "reuse rotor.npz", [str(f)])
        return rotor_description(np.load(f), fillet=params["fillet"])
