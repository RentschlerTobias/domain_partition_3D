"""Generate plots/slide01_3d_hub.png — 3D visualization of the T1_9 hub surface.

Loads T1_9_hub_raw.stl from the repo root and renders it with matplotlib's
plot_trisurf. No axes, no labels, no titles. Transparent background, portrait
orientation, DPI >= 200.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless rendering
import matplotlib.pyplot as plt
import meshio
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers 3d projection)


REPO_ROOT = Path(__file__).resolve().parent.parent
STL_PATH = REPO_ROOT / "T1_9_hub_raw.stl"
OUT_PATH = REPO_ROOT / "plots" / "slide01_3d_hub.png"


def main() -> None:
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)

    mesh = meshio.read(str(STL_PATH))
    points = mesh.points
    triangles = mesh.cells[0].data  # first cell block = triangle connectivity

    # Portrait orientation: height > width
    fig = plt.figure(figsize=(4, 6))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_trisurf(
        points[:, 0],
        points[:, 1],
        points[:, 2],
        triangles=triangles,
        color="lightgray",
        alpha=0.8,
        edgecolor="none",
    )
    ax.set_axis_off()

    spans = points.max(axis=0) - points.min(axis=0)
    ax.set_box_aspect((spans[0], spans[1], spans[2]))

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
