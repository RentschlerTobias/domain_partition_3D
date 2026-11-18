"""Building block ``band_layout``: structured grid of a cylinder surface band r = R, z0 <= z <= z1 with embedded
footprints (e.g. where blades attach), from layout2d on the unrolled band (x = R phi, periodic).

Inputs (from earlier steps): a dict with r, z0, z1 and ``loops`` = [(points (n, 2) in (R phi, z), pin indices)],
optional ``symmetry`` (the loops are copies under a rotation by 2 pi / symmetry: merges are applied to all copies).

Method: the layout is made on the band sheared by the direction the loops are stretched in (best_shear), the chord
collapse merges symmetric line sets while the fine mesh (in the physical plane, the band ends met at right angles)
keeps min sJ >= guard. With ``core_pairing`` (a later solid fill with a core) every column and the opposite one a core
sheet continues share one cell count. The result is a C x R tensor grid; every embedded loop covers a rectangle of it.
"""

from __future__ import annotations

import numpy as np

from hexdp.contracts import ParamSpec
from hexdp.core.layout2d import Domain2D, LayoutError, Loop, best_shear, collapse, mesh, sheared


class BandError(RuntimeError):
    pass


def core_counts(chosen):
    """Counts rule for a core: a column and the opposite column (same position from the other side of a q0 x q1
    core rectangle) get one count, the geometric mean; the core start and q0 that change the counts least."""
    def paired(ncol, off, q0):
        n, C = list(ncol), len(ncol)
        L = [q0, C // 2 - q0, q0, C // 2 - q0]
        start = np.cumsum([0] + L)
        for sd in (0, 1):
            for j in range(L[sd]):
                a = (off + start[sd] + j) % C
                b = (off + start[sd + 2] + L[sd] - 1 - j) % C
                n[a] = n[b] = max(1, int(round(np.sqrt(ncol[a] * ncol[b]))))
        return n

    def rule(ncol, nrow):
        C = len(ncol)
        if C % 2 or C < 4:
            raise LayoutError(f"{C} grid columns: a core needs an even number >= 4")
        cands = [(float(np.sum(np.abs(np.log(np.array(paired(ncol, o, q0)) / np.array(ncol))))),
                  abs(2 * q0 - C // 2), o, q0)
                 for q0 in sorted({C // 4, (C + 3) // 4}) if 1 <= q0 < C // 2 for o in range(C)]
        _, _, off, q0 = min(cands)
        chosen.update(offset=off, q0=q0)
        return paired(ncol, off, q0), nrow
    return rule


class BandLayout:
    """Structured grid of a cylinder band with embedded footprints (unrolled, periodic layout2d)."""

    name = "band_layout"
    params = {
        "h": ParamSpec(0.04, "fine cell size", unit="m", lower=1e-6),
        "nmin": ParamSpec(2, "minimum cells per block edge", lower=1),
        "guard": ParamSpec(0.1, "min scaled Jacobian a chord collapse must keep", lower=0, upper=1),
        "core_pairing": ParamSpec(True, "equal counts on the columns a core joins (for revolve_fill without bore)"),
        "embedded": ParamSpec(True, "loops are embedded (both sides are surface) rather than holes"),
    }

    def run(self, inputs, params, ctx=None):
        src = next(v for v in inputs.values() if isinstance(v, dict) and "loops" in v)
        return band_layout(src, h=params["h"], nmin=params["nmin"], guard=params["guard"],
                           core_pairing=params["core_pairing"], embedded=params["embedded"])


def band_layout(src, h=0.04, nmin=2, guard=0.1, core_pairing=True, embedded=True, log=None):
    R, Z0, Z1 = float(src["r"]), float(src["z0"]), float(src["z1"])
    NB = int(src.get("symmetry", 1))
    P = 2 * np.pi * R
    loops = [Loop(np.asarray(pts, float), pins=pins, embedded=embedded) for pts, pins in src["loops"]]
    nbd = 200
    bot = Loop(np.c_[np.linspace(0, P, nbd, endpoint=False), np.full(nbd, Z0)], shift=P)
    top = Loop(np.c_[np.linspace(P, 0, nbd, endpoint=False), np.full(nbd, Z1)], shift=-P)
    dom0 = Domain2D([bot, top, *loops], period=P)
    SH = best_shear(dom0)
    chosen = {}
    counts = core_counts(chosen) if core_pairing else None
    first = 2
    n_loops = len(loops)

    def images(lay, merge):
        """A merge with its images at the symmetric copies of the loops (loop i <-> i + k n/NB)."""
        ax, a, b = merge
        (la, ea), (lb, eb) = lay.verts[a], lay.verts[b]
        if NB <= 1 or la < first or lb < first:
            return [merge]
        per = n_loops // NB
        vid = {v: i for i, v in enumerate(lay.verts)}
        out = []
        for d in range(NB):
            ia = vid.get((first + (la - first + d * per) % n_loops, ea))
            ib = vid.get((first + (lb - first + d * per) % n_loops, eb))
            if ia is not None and ib is not None:
                out.append((ax, ia, ib))
        return out
    kw = dict(nmin=nmin, shear=SH, slide_free=True, counts=counts)
    try:
        lay, _ = collapse(sheared(dom0, SH), h=h, guard=guard, log=log, expand=images, mesh_dom=dom0, mesh_kw=kw)
        m = mesh(lay, dom0, h, **kw)
    except LayoutError as exc:
        raise BandError(f"band layout: {exc}") from exc
    if m.min_sj <= 0:
        raise BandError(f"band surface mesh has inverted cells (min sJ {m.min_sj:.3f})")
    C = len(lay.X)
    NI, NJ = int(m.oi[-1]), int(m.oj[-1])
    G2 = -np.ones((NI, NJ + 1), int)
    for k, (I, J) in enumerate(m.IJ):
        G2[I, J] = k
    if (G2 < 0).any():
        raise BandError("the band layout is not a full tensor grid")
    xs = m.X[:, 0]
    S3 = np.c_[R * np.cos(xs / R), R * np.sin(xs / R), m.X[:, 1]]
    return dict(kind="band", r=R, z0=Z0, z1=Z1, period=P, symmetry=NB, dom0=dom0, shear=SH, lay=lay, mesh=m,
                columns=C, rows=len(lay.Y) - 1, G2=G2, points=S3, core=dict(chosen) if core_pairing else None,
                loops=src["loops"], min_sj=float(m.min_sj))
