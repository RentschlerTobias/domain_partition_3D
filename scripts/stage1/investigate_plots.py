#!/usr/bin/env python3
"""Diagnostic plots for the periodic stage-1 partition.

Investigates the two open issues:
  (1) why the passage is not fully tiled (periodic seam treated as wall), and
  (2) why each blade tip carries 2 field singularities instead of 1.

Panels:
  A  unwrapped mesh + cross-field (4 arms/node) + singularities by index
     + periodic seam highlighted + blade-tip corners
  B  zoom on the left blade tip (LE)
  C  zoom on the right blade tip (TE)
  D  separatrices + quad blocks, seam highlighted, untiled regions visible

Run:  /root/venv/bin/python scripts/stage1/investigate_plots.py
"""

import sys
from pathlib import Path

import numpy as np
import torch
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, "/root/repos/domain_partition")
sys.path.insert(0, str(Path(__file__).resolve().parent))

import partition_surface as ps  # applies monkeypatches (periodic solver etc.)
from tools import FrameField
from dp_adapter import build_dp_data
from clean_separatrix import _poincare_indices, _singularity_coords

STL = "/root/repos/block_structured_meshing/T1_9_hub_raw.stl"
OUT = Path("/root/repos/block_structured_meshing/output/T1_9/hub_stage1")


def cross_arms(u):
    """4 unit cross directions per node from the representation field u=(cos4t,sin4t)."""
    base = np.arctan2(u[:, 1], u[:, 0]) / 4.0
    arms = []
    for k in range(4):
        a = base + k * np.pi / 2
        arms.append(np.column_stack([np.cos(a), np.sin(a)]))
    return arms  # list of (N,2)


def draw_field(ax, xy, tris, arms, scale, stride=3):
    ax.triplot(xy[:, 0], xy[:, 1], tris, lw=0.1, color="0.9")
    idx = np.arange(0, len(xy), stride)
    for arm in arms:
        for s in (1, -1):  # draw both ways -> looks like a cross
            ax.quiver(xy[idx, 0], xy[idx, 1],
                      s * arm[idx, 0], s * arm[idx, 1],
                      angles="xy", scale_units="xy", scale=1.0 / scale,
                      width=0.0015, headwidth=0, headlength=0, color="0.45",
                      alpha=0.7)


def main():
    mesh, tf = build_dp_data(STL)
    pp = mesh.periodic_pairs.numpy() if mesh.periodic_pairs.numel() else np.zeros((0, 2), int)
    ff = FrameField(mesh)
    m = ff.mesh
    xy = m.x[:, 0:2].numpy()
    tris = m.faces.T.numpy()
    u = m.u.numpy()
    arms = cross_arms(u)

    # singularities (representation index) on the periodic field
    sing = _poincare_indices(m)
    fids = torch.where(sing != 0)[0].tolist()
    sc = np.array([_singularity_coords(m, f).numpy() for f in fids]) if fids else np.zeros((0, 2))
    sidx = np.array([int(sing[f]) for f in fids])

    ct = m.corner_type.numpy()
    tip_xy = xy[ct == 1]
    outer_xy = xy[ct == 0]

    # mean edge length for arm scale
    e = np.linalg.norm(xy[tris[:, 0]] - xy[tris[:, 1]], axis=1).mean()
    scale = 0.9 * e

    def plot_field_panel(ax, title, xlim=None, ylim=None, stride=3):
        draw_field(ax, xy, tris, arms, scale, stride=stride)
        # seam edges
        for a, b in pp:
            ax.plot([xy[a, 0]], [xy[a, 1]], ".", color="orange", ms=3)
            ax.plot([xy[b, 0]], [xy[b, 1]], ".", color="orange", ms=3)
        if len(sc):
            pos = sidx > 0
            neg = sidx < 0
            ax.scatter(sc[neg, 0], sc[neg, 1], c="blue", s=90, zorder=6,
                       edgecolors="k", label="sing -1 (5-valent)")
            ax.scatter(sc[pos, 0], sc[pos, 1], c="red", s=90, zorder=6,
                       edgecolors="k", label="sing +1 (3-valent)")
        ax.scatter(tip_xy[:, 0], tip_xy[:, 1], c="magenta", marker="*", s=160,
                   zorder=7, label="blade tip (LE/TE)")
        ax.scatter(outer_xy[:, 0], outer_xy[:, 1], c="green", marker="^", s=70,
                   zorder=7, label="outer corner")
        ax.set_aspect("equal")
        ax.set_title(title)
        if xlim:
            ax.set_xlim(*xlim)
        if ylim:
            ax.set_ylim(*ylim)

    fig, axes = plt.subplots(2, 2, figsize=(18, 16))

    # Panel A: full field
    plot_field_panel(axes[0, 0], "A: cross field + singularities + seam (orange)")
    axes[0, 0].legend(loc="upper left", fontsize=8)

    # zoom windows around the two tip clusters
    if len(tip_xy) >= 2:
        order = np.argsort(tip_xy[:, 0])
        le = tip_xy[order[0]]
        te = tip_xy[order[-1]]
    else:
        le, te = np.array([0.2, 0.46]), np.array([0.77, 0.70])
    w = 0.18
    plot_field_panel(axes[0, 1], "B: zoom LE (left tip)",
                     xlim=(le[0] - w, le[0] + w), ylim=(le[1] - w, le[1] + w), stride=1)
    plot_field_panel(axes[1, 0], "C: zoom TE (right tip)",
                     xlim=(te[0] - w, te[0] + w), ylim=(te[1] - w, te[1] + w), stride=1)

    # Panel D: blocks + separatrices + seam
    block_mesh, mesh2, tf2 = ps.partition(STL, verbose=False)
    bx = block_mesh.x.numpy()
    bf = block_mesh.faces.numpy().T
    axD = axes[1, 1]
    axD.triplot(xy[:, 0], xy[:, 1], tris, lw=0.1, color="0.92")
    n_b = len(mesh2.streamlines) - len(mesh2.separatrices)
    for i, s in enumerate(mesh2.streamlines):
        s = np.asarray(s, float)
        if s.ndim != 2 or len(s) < 2:
            continue
        axD.plot(s[:, 0], s[:, 1], "0.6" if i < n_b else "C3",
                 lw=0.8 if i < n_b else 1.3)
    irr, inv, hi = ps._block_annotations(bx, bf)
    for qi, quad in enumerate(bf):
        ring = bx[list(quad) + [quad[0]]]
        axD.fill(ring[:, 0], ring[:, 1], alpha=0.25,
                 color="red" if qi in inv else "C0")
        axD.plot(ring[:, 0], ring[:, 1], "C0", lw=1.2)
    for a, b in pp:  # seam
        axD.plot([xy[a, 0], xy[b, 0]], [xy[a, 1], xy[b, 1]], color="orange",
                 lw=0.3, alpha=0.4)
    axD.scatter(xy[ct == 1][:, 0], xy[ct == 1][:, 1], c="magenta", marker="*", s=120, zorder=7)
    if len(irr):
        axD.scatter(bx[irr, 0], bx[irr, 1], facecolors="none", edgecolors="red",
                    s=180, linewidths=2, zorder=8, label=f"irregular interior ({len(irr)})")
        axD.legend(loc="upper left", fontsize=8)
    axD.set_aspect("equal")
    axD.set_title(f"D: {bf.shape[0]} blocks (seam=orange) — untiled outer channels visible")

    out = OUT / "investigate.png"
    fig.savefig(out, dpi=120, bbox_inches="tight")
    print(f"wrote {out}")
    print(f"singularities: {len(fids)}  indices={sidx.tolist()}")
    print(f"blade tips: LE~({le[0]:.2f},{le[1]:.2f}) TE~({te[0]:.2f},{te[1]:.2f})")


if __name__ == "__main__":
    main()
