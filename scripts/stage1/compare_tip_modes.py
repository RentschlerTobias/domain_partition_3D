#!/usr/bin/env python3
"""Side-by-side: blade-tip corners ON vs OFF (field tip singularities either way).

Left column  = MARK_BLADE_TIP_CORNERS True  (LE/TE are boundary split points;
               field tip singularities emit; artificial tip emission skipped)
Right column = MARK_BLADE_TIP_CORNERS False (blade is one smooth closed loop;
               only field tip singularities drive the partition)

Top row    = separatrices, Bottom row = quad blocks. Title shows block count.
"""
import importlib
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, "/root/repos/domain_partition")
sys.path.insert(0, str(Path(__file__).resolve().parent))

STL = "/root/repos/block_structured_meshing/T1_9_hub_raw.stl"
OUT = Path("/root/repos/block_structured_meshing/output/T1_9/hub_stage1")


def run(flag):
    import unwrap_surface as us
    us.set_blade_tip_corners(flag)
    import partition_surface as ps
    block_mesh, mesh, tf = ps.partition(STL, verbose=False)
    bx = block_mesh.x.numpy()
    bf = block_mesh.faces.numpy().T
    irr, inv, hi = ps._block_annotations(bx, bf)
    return mesh, bx, bf, irr, inv


def draw_seps(ax, mesh, title):
    xy = mesh.x[:, 0:2].numpy()
    tris = mesh.faces.T.numpy()
    ax.triplot(xy[:, 0], xy[:, 1], tris, lw=0.1, color="0.92")
    n_b = len(mesh.streamlines) - len(mesh.separatrices)
    for i, s in enumerate(mesh.streamlines):
        s = np.asarray(s, float)
        if s.ndim != 2 or len(s) < 2:
            continue
        ax.plot(s[:, 0], s[:, 1], "0.55" if i < n_b else "C3",
                lw=0.8 if i < n_b else 1.4)
    sc = xy[tris[mesh.singularities.numpy() != 0]].mean(axis=1)
    ax.scatter(sc[:, 0], sc[:, 1], c="blue", s=70, zorder=6)
    cm = mesh.x[:, 2].numpy() == 0
    ax.scatter(xy[cm, 0], xy[cm, 1], c="green", marker="^", s=55, zorder=6)
    ax.set_aspect("equal")
    ax.set_title(title)


def draw_blocks(ax, mesh, bx, bf, irr, inv, title):
    xy = mesh.x[:, 0:2].numpy()
    tris = mesh.faces.T.numpy()
    ax.triplot(xy[:, 0], xy[:, 1], tris, lw=0.1, color="0.93")
    for qi, quad in enumerate(bf):
        ring = bx[list(quad) + [quad[0]]]
        ax.fill(ring[:, 0], ring[:, 1], alpha=0.3,
                color="red" if qi in inv else "C0")
        ax.plot(ring[:, 0], ring[:, 1], "C0", lw=1.3)
    if len(irr):
        ax.scatter(bx[irr, 0], bx[irr, 1], facecolors="none", edgecolors="red",
                   s=170, linewidths=2, zorder=8)
    ax.set_aspect("equal")
    ax.set_title(title)


def main():
    fig, axes = plt.subplots(2, 2, figsize=(16, 16))
    for col, flag in enumerate((True, False)):
        # fresh import so module-level flags re-read cleanly
        for mod in ("unwrap_surface", "dp_adapter", "partition_surface",
                    "clean_separatrix"):
            if mod in sys.modules:
                importlib.reload(sys.modules[mod])
        mesh, bx, bf, irr, inv = run(flag)
        tag = "ON (LE/TE split pts)" if flag else "OFF (smooth blade)"
        draw_seps(axes[0, col], mesh,
                  f"separatrices — blade-tip corners {tag}")
        draw_blocks(axes[1, col], mesh, bx, bf, irr, inv,
                    f"{bf.shape[0]} blocks, {len(irr)} irregular — corners {tag}")
    out = OUT / "compare_tip_modes.png"
    fig.savefig(out, dpi=120, bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
