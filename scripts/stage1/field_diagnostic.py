#!/usr/bin/env python3
"""Diagnostic: build the cross field on the unwrapped hub and plot it with the
detected singularities. Validates the foundation before building a separatrix
tracer on top of mesh.u."""

import sys
from pathlib import Path

import numpy as np
import torch

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT.parent / "domain_partition_2D"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

# pull in the boundary-BC fix + monkeypatches
import partition_surface as ps  # noqa: E402
from tools import FrameField  # noqa: E402
from tools.singularity_detector import detect_singularities  # noqa: E402
from dp_adapter import build_dp_data  # noqa: E402


def main(stl, out_dir):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    mesh, tf = build_dp_data(stl)
    ff = FrameField(mesh)
    m = detect_singularities(ff.mesh)

    xy = m.x[:, 0:2].numpy()
    u = m.u.numpy()                     # representation field (cos4θ, sin4θ)
    tris = m.faces.T.numpy()
    sing_faces = np.where(m.singularities.numpy() != 0)[0]
    sing_idx = m.singularities.numpy()[sing_faces]

    # cross directions: angle = atan2(uy,ux)/4 + k*pi/2
    base = np.arctan2(u[:, 1], u[:, 0]) / 4.0
    print(f"singularities: {len(sing_faces)}  "
          f"(+1: {int((sing_idx>0).sum())}, -1: {int((sing_idx<0).sum())})")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(9, 7))
        ax.triplot(xy[:, 0], xy[:, 1], tris, lw=0.15, color="0.85")
        for k in range(2):  # two of the four cross directions (line field)
            ang = base + k * (np.pi / 2)
            ax.quiver(xy[:, 0], xy[:, 1], np.cos(ang), np.sin(ang),
                      color="C0", scale=40, width=0.0015, headwidth=0,
                      headlength=0, pivot="mid", alpha=0.6)
        # singularity barycenters
        sc = xy[tris[sing_faces]].mean(axis=1)
        pos = sing_idx > 0
        ax.scatter(sc[pos, 0], sc[pos, 1], c="red", s=60, zorder=5,
                   label=f"+idx ({int(pos.sum())})")
        ax.scatter(sc[~pos, 0], sc[~pos, 1], c="blue", s=60, zorder=5,
                   marker="s", label=f"-idx ({int((~pos).sum())})")
        ax.set_aspect("equal")
        ax.legend()
        ax.set_title("Hub cross field + singularities (unwrapped)")
        png = out_dir / "field_diagnostic.png"
        fig.savefig(png, dpi=140, bbox_inches="tight")
        print(f"wrote {png}")
    except Exception as e:  # noqa: BLE001
        print(f"plot skipped: {e}")


if __name__ == "__main__":
    stl = sys.argv[1] if len(sys.argv) > 1 else \
        str(REPO_ROOT / "T1_9_hub_raw.stl")
    out = str(REPO_ROOT / "output" / "T1_9" / "hub_stage1")
    main(stl, out)
