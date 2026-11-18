"""Building block ``revolve_fill``: the solid of revolution under a band surface grid (band_layout): r <= R between
z0 and z1, filled from the surface grid radially inward, either to a bore cylinder r = r_bore (a ring of C x R
blocks, no core) or to a core rectangle (q0 x q1 core blocks per row, needs core_pairing in band_layout).

Per fine level of the surface grid: radial lines from every surface node to the bore (same angle) or to the core
boundary (a rectangle planar per level and turned with the grid columns; the waviness of a level is taken by the
radial lines), the core interior by TFI. Boundary names: <prefix>_z0, <prefix>_z1, <prefix>_mantle (r = R),
<prefix>_bore (r = r_bore).
"""

from __future__ import annotations

import numpy as np

from hexdp.blocks import Part
from hexdp.contracts import ParamSpec


class RevolveError(RuntimeError):
    pass


class RevolveFill:
    """Fill the solid of revolution under a band surface grid: radial ring to a bore, or ring + core."""

    name = "revolve_fill"
    params = {
        "r_bore": ParamSpec(0.0, "bore radius (0 = full solid with a core)", unit="m", lower=0),
        "nr": ParamSpec(6, "cells from the surface to the bore / core", lower=1),
        "a_core": ParamSpec(0.45, "core rectangle half diagonal / R (full solid only)", lower=0.1, upper=0.9),
        "prefix": ParamSpec("hub", "boundary name prefix"),
    }

    def run(self, inputs, params, ctx=None):
        band = next(v for v in inputs.values() if isinstance(v, dict) and v.get("kind") == "band")
        return revolve_fill(band, r_bore=params["r_bore"], nr=params["nr"], a_core=params["a_core"],
                            prefix=params["prefix"])


def revolve_fill(band, r_bore=0.0, nr=6, a_core=0.45, prefix="hub"):
    m, R, C = band["mesh"], band["r"], band["columns"]
    oi, oj, G2 = m.oi, m.oj, band["G2"]
    NI, NJ = G2.shape[0], G2.shape[1] - 1
    pts = list(band["points"])
    ncol = [int((oi[c + 1] if c + 1 < C else NI) - oi[c]) for c in range(C)]
    if r_bore > 0:
        if r_bore >= R:
            raise RevolveError(f"bore radius {r_bore} >= surface radius {R}")
        RING = np.zeros((NI, nr + 1, NJ + 1), int)
        for J in range(NJ + 1):
            for I in range(NI):
                S = np.asarray(pts[G2[I, J]])
                B = np.array([S[0] * r_bore / R, S[1] * r_bore / R, S[2]])
                RING[I, 0, J] = G2[I, J]
                for k in range(1, nr + 1):
                    t = k / nr
                    pts.append((1 - t) * S + t * B)
                    RING[I, k, J] = len(pts) - 1
        grids = []
        for r in range(band["rows"]):
            J0, J1 = int(oj[r]), int(oj[r + 1])
            for col in range(C):
                I0, I1 = int(oi[col]), int(oi[col + 1]) if col + 1 < C else NI
                grids.append(RING[[i % NI for i in range(I0, I1 + 1)]][:, :, J0:J1 + 1])
        info = {"blocks": len(grids), "core": None, "r_bore": r_bore}
    else:
        core = band.get("core")
        if not core:
            raise RevolveError("a full solid needs band_layout with core_pairing (or give r_bore > 0)")
        off, q0 = core["offset"], core["q0"]
        q1 = C // 2 - q0
        L = [q0, q1, q0, q1]
        start = np.cumsum([0] + L)
        side = [sum(ncol[(off + start[k] + j) % C] for j in range(L[k])) for k in range(4)]
        if side[0] != side[2] or side[1] != side[3]:
            raise RevolveError(f"core sides {side} do not match (band_layout counts)")
        sA, sB = side[0], side[1]
        per_pq = ([(p, 0) for p in range(sA)] + [(sA, k) for k in range(sB)] + [(sA - p, sB) for p in range(sA)]
                  + [(0, sB - k) for k in range(sB)])
        base = int(oi[off])
        RING = np.zeros((NI, nr + 1, NJ + 1), int)
        CORE = np.zeros((sA + 1, sB + 1, NJ + 1), int)
        half = a_core * R / np.sqrt(2)
        for J in range(NJ + 1):
            Sj = np.array([pts[i] for i in G2[:, J]])
            ang = np.arctan2(Sj[:, 1], Sj[:, 0])
            d = [ang[int(oi[(off + start[k]) % C])] - np.radians(-135 + 90 * k) for k in range(4)]
            theta = float(np.angle(np.mean(np.exp(1j * np.array(d)))))
            c, s = np.cos(theta), np.sin(theta)
            zc = Sj[:, 2].mean()
            for I in range(NI):
                p, qq = per_pq[(I - base) % NI]
                u, v = -half + 2 * half * p / sA, -half + 2 * half * qq / sB
                B = np.array([c * u - s * v, s * u + c * v, zc])
                RING[I, 0, J] = G2[I, J]
                for k in range(1, nr + 1):
                    t = k / nr
                    pts.append((1 - t) * Sj[I] + t * B)
                    RING[I, k, J] = len(pts) - 1
            Zb = {per_pq[(I - base) % NI]: RING[I, nr, J] for I in range(NI)}
            for p in range(sA + 1):
                for qq in range(sB + 1):
                    if (p, qq) in Zb:
                        CORE[p, qq, J] = Zb[(p, qq)]
                        continue
                    u, v = p / sA, qq / sB
                    X = lambda a, b: pts[Zb[(a, b)]]          # noqa: E731
                    val = ((1 - v) * X(p, 0) + v * X(p, sB) + (1 - u) * X(0, qq) + u * X(sA, qq)
                           - ((1 - u) * (1 - v) * X(0, 0) + u * (1 - v) * X(sA, 0) + u * v * X(sA, sB)
                              + (1 - u) * v * X(0, sB)))
                    pts.append(val)
                    CORE[p, qq, J] = len(pts) - 1
        grids = []
        pc = list(np.cumsum([0] + [ncol[(off + j) % C] for j in range(q0)]))
        qc = list(np.cumsum([0] + [ncol[(off + q0 + j) % C] for j in range(q1)]))
        for r in range(band["rows"]):
            J0, J1 = int(oj[r]), int(oj[r + 1])
            for col in range(C):
                I0, I1 = int(oi[col]), int(oi[col + 1]) if col + 1 < C else NI
                grids.append(RING[[i % NI for i in range(I0, I1 + 1)]][:, :, J0:J1 + 1])
            for a in range(q0):
                for b in range(q1):
                    grids.append(CORE[pc[a]:pc[a + 1] + 1, qc[b]:qc[b + 1] + 1, J0:J1 + 1])
        info = {"blocks": len(grids), "core": [q0, q1], "r_bore": 0.0}
    z0, z1 = band["z0"], band["z1"]

    def label(F):
        r = np.hypot(F[:, 0], F[:, 1])
        if np.ptp(F[:, 2]) < 1e-9 and abs(F[0, 2] - z0) < 1e-6:
            return f"{prefix}_z0"
        if np.ptp(F[:, 2]) < 1e-9 and abs(F[0, 2] - z1) < 1e-6:
            return f"{prefix}_z1"
        if np.abs(r - R).max() < 1e-6:
            return f"{prefix}_mantle"
        if r_bore > 0 and np.abs(r - r_bore).max() < 1e-6:
            return f"{prefix}_bore"
        return None
    return Part(prefix, np.array(pts), grids, label, info)
