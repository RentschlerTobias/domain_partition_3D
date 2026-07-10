#!/usr/bin/env python3
"""
Stage 1 / Slide 03: 4-RoSy cross-field on the unwrapped hub surface,
restricted to boundary nodes only.

Reuses the cross-field logic from scripts/stage1/field_diagnostic.py
(lines 44-62) without modifying it. The cross-field is computed by
``FrameField`` (which sets ``mesh.u`` as the (cos4θ, sin4θ) representation
field). The two line-field directions are then drawn as small orthogonal
quiver arrows at every boundary node.

Visual recipe (matches slide02_unwrapped.png styling):
  - Mesh triangulation background: very subtle, light gray (lw=0.1, color="0.93")
  - Outer contour: thick black (lw=2.0, color="black")
  - Cross-field arrows on boundary nodes only: light blue (#87CEEB),
    scale=40, no arrow heads
  - No axes, no labels, no titles, no legends
  - Transparent background (alpha channel)
  - Portrait orientation (figsize=(4, 6))
  - DPI >= 200

Output: plots/slide03_crossfield.png
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
OUT_PATH = REPO_ROOT / "plots" / "step03_crossfield.png"

import matplotlib

matplotlib.use("Agg")  # headless rendering
import matplotlib.pyplot as plt
import numpy as np
import torch

from dp_adapter import build_dp_data  # noqa: E402
from tools import FrameField  # noqa: E402
from unwrap_surface import unwrap  # noqa: E402


def main() -> None:
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)

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
    base = np.arctan2(u[:, 1], u[:, 0]) / 4.0  # principal cross angle

    # --- 4. Unwrapped (s, t) coordinates for plotting ---
    # build_dp_data does NOT renumber nodes (see dp_adapter.py comment),
    # so boundary_node_ids index directly into unwrap()["st"].
    st_data = unwrap(STL_PATH)
    st = st_data["st"]
    tris = st_data["tris"]
    loops = st_data["loops"]

    # Boundary node positions in (s, t) space.
    bxy = st[boundary_node_ids]

    # --- 5. Plot ---
    fig, ax = plt.subplots(figsize=(4, 6))  # portrait

    # Mesh triangulation background: very subtle, light gray.
    ax.triplot(st[:, 0], st[:, 1], tris, lw=0.1, color="0.93")

    # Outer contour: thick black (same as slide02).
    outer_loop = loops[0]
    ring = st[outer_loop + [outer_loop[0]]]
    ax.plot(ring[:, 0], ring[:, 1], lw=2.0, color="black")

    # Inner (blade) loops: also thick black so the blade-shaped hole is visible.
    for inner in loops[1:]:
        inner_ring = st[inner + [inner[0]]]
        ax.plot(inner_ring[:, 0], inner_ring[:, 1], lw=2.0, color="black")

    # 4-RoSy cross-field on boundary nodes only: two orthogonal line-field
    # directions per node (k=0 and k=1 give the full cross as a "+").
    for k in range(2):
        ang = base[boundary_node_ids] + k * (np.pi / 2)
        ax.quiver(
            bxy[:, 0], bxy[:, 1],
            np.cos(ang), np.sin(ang),
            color="#3aa6d0",   # blue, higher contrast than the pale field blue
            scale=18,          # smaller scale => longer arms; balanced size
            width=0.004,
            headwidth=0,
            headlength=0,
            pivot="mid",
            alpha=0.95,
        )

    ax.set_aspect("equal")
    ax.set_axis_off()

    fig.savefig(
        OUT_PATH,
        dpi=200,
        transparent=True,
        bbox_inches="tight",
    )
    plt.close(fig)

    print(f"wrote {OUT_PATH}")
    print(f"  boundary nodes: {len(boundary_node_ids)}")
    print(f"  total nodes:    {len(st)}")
    print(f"  triangles:      {len(tris)}")


if __name__ == "__main__":
    main()
