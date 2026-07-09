#!/usr/bin/env python3
"""Plot the 2D domain partition (block structure) of hub and shroud side by
side. No hexa. The shroud partition is the hub T-a template morphed onto the
shroud surface (HubToShroudMorph): the two cross fields are topologically
different (hub 4x idx=-1, shroud 4x idx=+1), so a directly-run shroud field
gives an incompatible 6-block/4-reject partition -- the morphed hub template
is the usable shared block topology. The actual shroud blade loop is overlaid
(red dashed) to show the morph lands on it.
"""
import sys
import os
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, "/root/repos/domain_partition")
sys.path.insert(0, str(Path(__file__).resolve().parent))

import unwrap_surface as us
import clean_separatrix as cs
import tmesh_partition as tp
import hexa_interpolation as hx
from dp_adapter import build_dp_data

OUT = Path("/root/repos/block_structured_meshing/output/T1_9/"
           "domain_partition_2d")


def draw(ax, res, mesh, ttl, blade_actual=None):
    xy = mesh.x[:, 0:2].numpy()
    tris = mesh.faces.T.numpy()
    ax.triplot(xy[:, 0], xy[:, 1], tris, lw=0.08, color="0.94", zorder=0)
    for bi, blk in enumerate(res["blocks"]):
        ring = blk["ring"]
        ax.fill(ring[:, 0], ring[:, 1], alpha=0.22, color="C0", zorder=1)
        for s in blk["sides"]:
            ax.plot(s[:, 0], s[:, 1], "C0", lw=1.4, zorder=2)
        c = blk["centroid"]
        ax.text(c[0], c[1], str(bi), fontsize=9, color="darkblue",
                ha="center", va="center", zorder=4,
                bbox=dict(fc="white", ec="C0", alpha=0.8, pad=0.8))
    cn = np.unique(np.concatenate([b["corners"] for b in res["blocks"]]))
    P = res["nodes"][cn]
    ax.scatter(P[:, 0], P[:, 1], c="k", s=16, zorder=3)
    if blade_actual is not None:
        bl = np.asarray(blade_actual, float)
        ax.plot(bl[:, 0], bl[:, 1], "r--", lw=1.2, zorder=5,
                label="actual shroud blade loop")
        ax.legend(loc="upper left", fontsize=9)
    ax.set_aspect("equal")
    ax.set_title(ttl, fontsize=12)


def main():
    us.set_blade_tip_corners(False)
    cs.set_emanate_outer_corners(True)
    OUT.mkdir(parents=True, exist_ok=True)

    hub = tp.run_tmesh(stl=hx.HUB_STL, out_dir=OUT, make_plots=False,
                       tag="hub", verbose=True)
    shroud_mesh, tf_s = build_dp_data(hx.SHROUD_STL)
    morph = hx.HubToShroudMorph(hub["mesh"], shroud_mesh)
    sh = hx.morph_result(hub["result"], morph)

    fig, axs = plt.subplots(1, 2, figsize=(20, 9))
    draw(axs[0], hub["result"], hub["mesh"],
         f"HUB domain partition (r=0.5): {len(hub['result']['blocks'])} blocks")
    draw(axs[1], sh, shroud_mesh,
         f"SHROUD domain partition (r=1.9, hub template morphed): "
         f"{len(sh['blocks'])} blocks",
         blade_actual=shroud_mesh.blade_loops[0])
    out = OUT / "hub_shroud_blocking.png"
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
