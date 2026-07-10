#!/usr/bin/env python3
"""
Stage 1 / Slide 02: 2D unwrapped hub surface visualization.

Reuses the unwrap logic from scripts/stage1/unwrap_surface.py and renders
the (s, t) domain as a clean line drawing suitable for a presentation slide:

  - Mesh triangulation in background: very subtle, light gray (lw=0.1, color="0.93")
  - Outer boundary: thick black lines (lw=2.0, color="black")
  - No filling, no axes, no labels, no titles
  - Transparent background (alpha channel)
  - Portrait orientation (figsize=(4, 6))
  - DPI >= 200

Output: plots/slide02_unwrapped.png
"""

import sys
from pathlib import Path

# Make scripts/stage1/ importable so we can reuse unwrap_surface.py as-is.
SCRIPT_DIR = Path(__file__).resolve().parent
STAGE1_DIR = SCRIPT_DIR / "stage1"
if str(STAGE1_DIR) not in sys.path:
    sys.path.insert(0, str(STAGE1_DIR))

REPO_ROOT = SCRIPT_DIR.parent
STL_PATH = REPO_ROOT / "T1_9_hub_raw.stl"
OUT_PATH = REPO_ROOT / "plots" / "slide02_unwrapped.png"

import matplotlib

matplotlib.use("Agg")  # headless rendering
import matplotlib.pyplot as plt

from unwrap_surface import unwrap  # noqa: E402


def main() -> None:
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)

    data = unwrap(STL_PATH)
    st = data["st"]
    tris = data["tris"]
    loops = data["loops"]

    # Portrait figure; equal aspect so the (s, t) domain is not stretched.
    fig, ax = plt.subplots(figsize=(4, 6))

    # Mesh triangulation background: very subtle, light gray.
    ax.triplot(st[:, 0], st[:, 1], tris, lw=0.1, color="0.93")

    # Outer boundary: thick black lines. loops[0] is the outer loop
    # (largest absolute enclosed area in (s, t) - see unwrap_surface.py).
    outer_loop = loops[0]
    ring = st[outer_loop + [outer_loop[0]]]
    ax.plot(ring[:, 0], ring[:, 1], lw=2.0, color="black")

    # Inner (blade) loops: also draw as thick black lines so the blade-shaped
    # hole is visible. Same styling as the outer boundary.
    for inner in loops[1:]:
        inner_ring = st[inner + [inner[0]]]
        ax.plot(inner_ring[:, 0], inner_ring[:, 1], lw=2.0, color="black")

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


if __name__ == "__main__":
    main()
