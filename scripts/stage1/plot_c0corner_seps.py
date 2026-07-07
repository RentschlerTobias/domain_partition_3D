#!/usr/bin/env python3
"""Step-0 experiment: does emitting separatrices from the OUTER C0 corners alone
tile the outer channels into quads? Same config as plot_xiao_clip_noletE.py
(Xiao + own clipping, no LE/TE) but with EMANATE_OUTER_CORNERS on.

Two panels: separatrices (left), quad blocks (right). Saved to
output/T1_9/hub_stage1/xiao_clip_c0corners.png.
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, "/root/repos/domain_partition")
sys.path.insert(0, str(Path(__file__).resolve().parent))

import unwrap_surface as us
import clean_separatrix as cs
import partition_surface as ps

STL = "/root/repos/block_structured_meshing/T1_9_hub_raw.stl"
OUT = Path("/root/repos/block_structured_meshing/output/T1_9/hub_stage1")


def main():
    us.set_blade_tip_corners(False)
    cs.set_emanate_outer_corners(True)     # the experiment toggle
    block_mesh, mesh, tf = ps.partition(STL, verbose=True)

    xy = mesh.x[:, 0:2].numpy()
    tris = mesh.faces.T.numpy()
    fig, ax = plt.subplots(1, 2, figsize=(20, 10))

    ax[0].triplot(xy[:, 0], xy[:, 1], tris, lw=0.1, color="0.92")
    n_b = len(mesh.streamlines) - len(mesh.separatrices)
    for i, s in enumerate(mesh.streamlines):
        s = np.asarray(s, float)
        if s.ndim != 2 or len(s) < 2:
            continue
        ax[0].plot(s[:, 0], s[:, 1], "0.55" if i < n_b else "C3",
                   lw=0.8 if i < n_b else 1.4)
    sc = xy[tris[mesh.singularities.numpy() != 0]].mean(axis=1)
    ax[0].scatter(sc[:, 0], sc[:, 1], c="blue", s=70, zorder=6, label="singularity")
    cm = mesh.x[:, 2].numpy() == 0
    ax[0].scatter(xy[cm, 0], xy[cm, 1], c="green", marker="^", s=55, zorder=6,
                  label="c0 corner")
    ax[0].set_aspect("equal")
    ax[0].legend()
    ax[0].set_title("C0-corner separatrices: separatrices")

    ax[1].triplot(xy[:, 0], xy[:, 1], tris, lw=0.1, color="0.93")
    if block_mesh.faces is not None and block_mesh.faces.numel() > 0:
        bx = block_mesh.x.numpy()
        bf = block_mesh.faces.numpy().T
        irr, inv, hi = ps._block_annotations(bx, bf)
        for qi, quad in enumerate(bf):
            ring = bx[list(quad) + [quad[0]]]
            ax[1].fill(ring[:, 0], ring[:, 1], alpha=0.3,
                       color="red" if qi in inv else "C0")
            ax[1].plot(ring[:, 0], ring[:, 1], "C0", lw=1.3)
        if len(irr):
            ax[1].scatter(bx[irr, 0], bx[irr, 1], facecolors="none",
                          edgecolors="red", s=170, linewidths=2, zorder=8)
        ttl = f"{bf.shape[0]} blocks, {len(irr)} irregular, {len(inv)} inverted"
    else:
        ttl = "no blocks"
    ax[1].set_aspect("equal")
    ax[1].set_title(f"C0-corner separatrices: {ttl}")

    out = OUT / "xiao_clip_c0corners.png"
    fig.savefig(out, dpi=120, bbox_inches="tight")
    print(f"wrote {out}")
    print(f"RESULT: {ttl}")


if __name__ == "__main__":
    main()
