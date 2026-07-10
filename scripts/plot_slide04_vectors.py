#!/usr/bin/env python3
"""
Stage 1 / Slide 04: Representative vectors in a 4-RoSy field.

Two-panel figure (saved as a single PNG):
  - Panel 1 (conceptual): Large cross with 4 directions (light blue) + red
    representative vector at 22.5° showing 4-RoSy symmetry.
  - Panel 2 (realistic): Unwrapped cylinder with red 4 small orthogonal
    vectors + blue representative vector, mesh background.

Reuses the cross-field logic from scripts/stage1/field_diagnostic.py
without modifying it. The cross-field is computed by ``FrameField`` (which
sets ``mesh.u`` as the (cos4θ, sin4θ) representation field). The four
cross directions are then drawn as small orthogonal quiver arrows at every
boundary node.

Visual recipe:
  - Panel 1: 4 light-blue arrows at 22.5°, 112.5°, 202.5°, 292.5° (the
    4-RoSy cross), with the one at 22.5° overlaid in red as the
    representative vector.
  - Panel 2: mesh triangulation background (light gray), outer + inner
    contours (black), 4 small red orthogonal vectors per boundary node
    (finer than slide03: scale=60 vs scale=40), and one blue
    representative vector per boundary node (the principal direction θ).
  - No axes, no labels, no titles, no legends.
  - Transparent background (alpha channel).
  - Portrait orientation (figsize=(8, 6) for 2-panel).
  - DPI >= 200.

Output: plots/slide04_vectors.png
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
OUT_A = REPO_ROOT / "plots" / "step04_representatives.png"
OUT_B = REPO_ROOT / "plots" / "step04_mapping.png"

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

    # --- 2. Identify boundary nodes (same criterion as FrameField itself) ---
    mask_boundary_edges = mesh.edge_attr == 1
    boundary_node_ids = torch.unique(
        mesh.edge_index[0, mask_boundary_edges]
    ).numpy()

    # --- 3. Cross-field directions (same formula as field_diagnostic.py L37) ---
    u = mesh.u.numpy()  # representation field (cos4θ, sin4θ)
    base = np.arctan2(u[:, 1], u[:, 0]) / 4.0  # principal cross angle θ

    # --- 4. Unwrapped (s, t) coordinates for plotting ---
    # build_dp_data does NOT renumber nodes (see dp_adapter.py comment),
    # so boundary_node_ids index directly into unwrap()["st"].
    st_data = unwrap(STL_PATH)
    st = st_data["st"]
    tris = st_data["tris"]
    loops = st_data["loops"]

    # Boundary node positions in (s, t) space.
    bxy = st[boundary_node_ids]

    # ============================================================
    # Plot A (realistic): cross-field + representative vectors on
    # the unwrapped surface (boundary nodes). Same as step03 plus the
    # blue representative vector per node.
    # ============================================================
    figA, axA = plt.subplots(figsize=(4, 6))  # portrait

    axA.triplot(st[:, 0], st[:, 1], tris, lw=0.1, color="0.93")

    outer_loop = loops[0]
    ring = st[outer_loop + [outer_loop[0]]]
    axA.plot(ring[:, 0], ring[:, 1], lw=2.0, color="black")
    for inner in loops[1:]:
        inner_ring = st[inner + [inner[0]]]
        axA.plot(inner_ring[:, 0], inner_ring[:, 1], lw=2.0, color="black")

    # 4-RoSy cross per boundary node (light blue, enlarged like step03).
    for k in range(4):
        ang = base[boundary_node_ids] + k * (np.pi / 2)
        axA.quiver(
            bxy[:, 0], bxy[:, 1],
            np.cos(ang), np.sin(ang),
            color="#3aa6d0",   # blue cross arms (same as step03)
            scale=9,
            width=0.006,
            headwidth=0,
            headlength=0,
            pivot="mid",
            alpha=0.85,
        )

    # Representative vector per boundary node (principal direction θ), red,
    # with an arrow head so the "chosen" direction stands out from the cross.
    axA.quiver(
        bxy[:, 0], bxy[:, 1],
        np.cos(base[boundary_node_ids]), np.sin(base[boundary_node_ids]),
        color="red",
        scale=9,
        width=0.007,
        headwidth=4,
        headlength=5,
        pivot="tail",
        alpha=1.0,
    )

    axA.set_aspect("equal")
    axA.set_axis_off()
    figA.savefig(OUT_A, dpi=200, transparent=True, bbox_inches="tight")
    plt.close(figA)
    print(f"wrote {OUT_A}")

    # ============================================================
    # Plot B (conceptual): cross -> representative vector mapping.
    # A 4-RoSy cross at example angle θ. The representative vector is
    # the 4-RoSy image: angle = 4 * min(cross angles) = 4·θ, matching
    # frame_field.map_cross_vectors_to_reference_vector (ref = 4*min).
    # All arrows start at the origin (pivot="tail").
    # ============================================================
    theta_deg = 20.0                       # example smallest-to-x-axis angle
    theta = np.deg2rad(theta_deg)
    rep = 4.0 * theta                      # representative vector angle = 4·θ

    figB, axB = plt.subplots(figsize=(4, 6))
    axB.set_xlim(-1.4, 1.4)
    axB.set_ylim(-1.4, 1.4)
    axB.set_aspect("equal")

    # x-axis reference line (thin gray) — the axis θ is measured against.
    axB.axhline(0.0, color="0.75", lw=1.0, zorder=1)

    # The 4 orthogonal cross arms (blue) from the origin.
    for k in range(4):
        ang = theta + k * (np.pi / 2)
        axB.quiver(
            0.0, 0.0, np.cos(ang), np.sin(ang),
            color="#1f77ff", scale=2.4, width=0.012,
            headwidth=4, headlength=5, pivot="tail", zorder=3, alpha=0.95,
        )

    # Representative vector (red) at 4·θ from the origin.
    axB.quiver(
        0.0, 0.0, np.cos(rep), np.sin(rep),
        color="red", scale=2.4, width=0.014,
        headwidth=4, headlength=5, pivot="tail", zorder=4,
    )

    # Minimal schematic annotations for θ and 4·θ.
    axB.text(0.62 * np.cos(theta / 2), 0.62 * np.sin(theta / 2),
             r"$\theta$", color="0.35", fontsize=12, ha="center", va="center")
    axB.text(np.cos(rep) - 0.22, np.sin(rep) + 0.02,
             r"$4\theta$", color="red", fontsize=12, ha="right", va="center")

    axB.set_axis_off()
    figB.savefig(OUT_B, dpi=200, transparent=True, bbox_inches="tight")
    plt.close(figB)
    print(f"wrote {OUT_B}")

    print(f"  boundary nodes: {len(boundary_node_ids)}")
    print(f"  total nodes:    {len(st)}")
    print(f"  triangles:      {len(tris)}")


if __name__ == "__main__":
    main()
