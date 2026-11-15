"""Hex cell classification shared by every blade_pipeline report.

A negative min corner scaled Jacobian does NOT mean "inverted": a valid cell
with one corner angle above 180 deg (concave) has a negative corner Jacobian
too (seen on canada). The two cases are told apart the way a finite-volume
code sees a cell: split it into 24 tetrahedra (cell centroid, face centroid,
one face edge). OpenFOAM's volume / checkMesh works on this decomposition.

  inverted : some centroid tet has volume <= 0 (folded, or negative volume)
  concave  : all 24 tets positive, but a corner scaled Jacobian <= 0
  ok       : all corner scaled Jacobians > 0

Coarse blocks are gated on `ok` (a concave block corner turns into folded fine
cells in the TFI); fine meshes on `inverted == 0` plus checkMesh, with the
concave count reported separately.
"""
import sys
from pathlib import Path

import numpy as np

from hexdp.core import clean_blocks as cb

# VTK hex faces, outward winding
FACES = np.array([(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4),
                  (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)])


def centroid_tets(P, H):
    """(ncell, 24) signed volumes of the centroid decomposition."""
    X = P[H]                                   # (n, 8, 3)
    c = X.mean(1)                              # cell centroid
    F = X[:, FACES]                            # (n, 6, 4, 3)
    fc = F.mean(2)                             # face centroids
    a = F
    b = np.roll(F, -1, axis=2)
    # tet (c, fc, a, b) with outward face winding -> positive for a valid cell
    v = np.einsum("nfkx,nfkx->nfk",
                  a - c[:, None, None, :],
                  np.cross(b - c[:, None, None, :], fc[:, :, None, :] - c[:, None, None, :])) / 6.0
    return v.reshape(len(H), 24)


def classify(P, H):
    P = np.asarray(P, float)
    H = np.asarray(H, int)
    tets = centroid_tets(P, H)
    sj = cb.scaled_jacobians(P, H)
    vol = tets.sum(1)
    # orientation: the mesh convention is the sign of the majority of cells
    s = 1.0 if np.median(vol) >= 0 else -1.0
    inv = (s * tets <= 0).any(1)
    concave = ~inv & (sj <= 0)
    return {"volume": s * vol, "min_tet": (s * tets).min(1), "sj": sj,
            "inverted": inv, "concave": concave, "ok": ~inv & ~concave}


def summary(P, H):
    c = classify(P, H)
    return {"cells": int(len(H)), "inverted": int(c["inverted"].sum()),
            "concave": int(c["concave"].sum()), "ok": int(c["ok"].sum()),
            "min_sj": float(c["sj"].min()), "p05_sj": float(np.percentile(c["sj"], 5)),
            "mean_sj": float(c["sj"].mean()), "min_volume": float(c["volume"].min())}
