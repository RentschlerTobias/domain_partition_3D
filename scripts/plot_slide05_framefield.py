#!/usr/bin/env python3
"""
Stage 1 / Slide 05: Smoothed vector field + frame field (back-mapping).

Two-panel figure (saved as a single PNG):
  - Panel 1 (left): Smoothed vector field (after harmonic solve) — vector
    arrows showing the smoothed field direction. The field is the
    (cos4θ, sin4θ) representation field stored in ``mesh.u`` after
    ``FrameField(mesh)`` runs the harmonic solve + per-node normalization.
  - Panel 2 (right): Frame field (back-mapping: vector → cross) — many
    small crosses (4 directions) covering the domain. The cross directions
    are recovered from the representative vector via
    ``base = arctan2(u[:,1], u[:,0]) / 4`` and ``ang = base + k*π/2`` for
    k=0,1,2,3 (the full 4-RoSy cross).

Reuses the cross-field logic from scripts/stage1/field_diagnostic.py
without modifying it. The cross-field is computed by ``FrameField`` (which
sets ``mesh.u`` as the (cos4θ, sin4θ) representation field).

Visual recipe:
  - Panel 1: mesh triangulation background (light gray), outer + inner
    contours (black), smoothed vector field arrows at every node
    (light blue, scale=40, no arrow heads).
  - Panel 2: mesh triangulation background (light gray), outer + inner
    contours (black), 4 small orthogonal crosses per node (light blue,
    scale=40, no arrow heads) — the full 4-RoSy cross at every node.
  - No axes, no labels, no titles, no legends.
  - Transparent background (alpha channel).
  - Portrait orientation (figsize=(8, 6) for 2-panel).
  - DPI >= 200.

Output: plots/slide05_framefield.png
"""

import sys
from pathlib import Path

# Make scripts/stage1/ importable so we can reuse field_diagnostic.py logic
# (FrameField, build_dp_data, unwrap_surface) as-is.
SCRIPT_DIR = Path(__file__).resolve().parent
STAGE1_DIR = SCRIPT_DIR / "stage1"
if str(STAGE1_DIR) not in sys.path:
    sys.path.insert(0, str(STAGE1_DIR))

REPO_ROOT = SCRIPT_DIR.parent
STL_PATH = REPO_ROOT / "T1_9_hub_raw.stl"
OUT_A = REPO_ROOT / "plots" / "step05_vectorfield.png"
OUT_B = REPO_ROOT / "plots" / "step05_crossfield.png"

import matplotlib

matplotlib.use("Agg")  # headless rendering
import matplotlib.pyplot as plt
import numpy as np
import torch

from dp_adapter import build_dp_data  # noqa: E402
from tools import FrameField  # noqa: E402
from unwrap_surface import unwrap  # noqa: E402


def main() -> None:
    OUT_A.parent.mkdir(parents=True, exist_ok=True)

    # --- 1. Build the mesh + cross-field (reuses field_diagnostic.py logic) ---
    mesh, _transform = build_dp_data(STL_PATH)
    ff = FrameField(mesh)

    # --- 2. Smoothed vector field (after harmonic solve + normalization) ---
    # mesh.u is the (cos4θ, sin4θ) representation field, normalized per-node.
    u = mesh.u.numpy()  # shape (N, 2)

    # --- 3. Cross-field directions (same formula as field_diagnostic.py L37) ---
    base = np.arctan2(u[:, 1], u[:, 0]) / 4.0  # principal cross angle θ

    # --- 4. Unwrapped (s, t) coordinates for plotting ---
    # build_dp_data does NOT renumber nodes (see dp_adapter.py comment),
    # so node indices map directly into unwrap()["st"].
    st_data = unwrap(STL_PATH)
    st = st_data["st"]
    tris = st_data["tris"]
    loops = st_data["loops"]

    def _draw_frame(ax):
        """Mesh + outer/blade contours — shared background for both plots."""
        ax.triplot(st[:, 0], st[:, 1], tris, lw=0.1, color="0.93")
        outer_loop = loops[0]
        ring = st[outer_loop + [outer_loop[0]]]
        ax.plot(ring[:, 0], ring[:, 1], lw=2.0, color="black")
        for inner in loops[1:]:
            inner_ring = st[inner + [inner[0]]]
            ax.plot(inner_ring[:, 0], inner_ring[:, 1], lw=2.0, color="black")

    # ============================================================
    # Plot A: Smoothed vector field (after harmonic solve).
    # mesh.u is the (cos4θ, sin4θ) representation field — the smoothed
    # unit vectors after the harmonic solve + per-node normalization.
    # ============================================================
    figA, axA = plt.subplots(figsize=(4, 6))  # portrait
    _draw_frame(axA)
    axA.quiver(
        st[:, 0], st[:, 1],
        u[:, 0], u[:, 1],
        color="#87CEEB",  # light blue
        scale=26,
        width=0.0022,
        headwidth=0,
        headlength=0,
        pivot="mid",
        alpha=0.9,
    )
    axA.set_aspect("equal")
    axA.set_axis_off()
    figA.savefig(OUT_A, dpi=200, transparent=True, bbox_inches="tight")
    plt.close(figA)
    print(f"wrote {OUT_A}")

    # ============================================================
    # Plot B: Frame field (back-mapping vector → cross). 4 orthogonal
    # cross arms per node recovered from the representative vector:
    #   base = arctan2(u[:,1], u[:,0]) / 4 ; ang_k = base + k·π/2.
    # ============================================================
    figB, axB = plt.subplots(figsize=(4, 6))  # portrait
    _draw_frame(axB)
    for k in range(4):
        ang = base + k * (np.pi / 2)
        axB.quiver(
            st[:, 0], st[:, 1],
            np.cos(ang), np.sin(ang),
            color="#87CEEB",  # light blue
            scale=26,
            width=0.0022,
            headwidth=0,
            headlength=0,
            pivot="mid",
            alpha=0.9,
        )
    axB.set_aspect("equal")
    axB.set_axis_off()
    figB.savefig(OUT_B, dpi=200, transparent=True, bbox_inches="tight")
    plt.close(figB)
    print(f"wrote {OUT_B}")

    print(f"  total nodes:    {len(st)}")
    print(f"  triangles:      {len(tris)}")


if __name__ == "__main__":
    main()
