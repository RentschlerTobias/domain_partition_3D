#!/usr/bin/env python3
"""Zoom plot for the problematic streamline region (indices 20-27 area).

Usage:
    cd /home/t1dde/Duty/projects/domain_partition/domain_partition_3D
    python scripts/stage1/run_xiao_only_zoom.py [hub|shroud]

Outputs to:
    output/T1_9/<part>_stage1/xiao_only/xiao_zoom_problematic.png
"""

import json
import sys
from pathlib import Path

import numpy as np

_DP2D = Path(__file__).resolve().parent.parent.parent.parent / "domain_partition_2D"
sys.path.insert(0, str(_DP2D))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import unwrap_surface as us                                    # noqa: E402
import clean_separatrix as cs                                  # noqa: E402
import partition_surface as ps                                 # noqa: E402
import tmesh_faces as tmf                                      # noqa: E402
from dp_adapter import build_dp_data                           # noqa: E402
from tools import FrameField, StreamlineGenerator_v2           # noqa: E402
from tools.singularity_detector import detect_singularities    # noqa: E402
from tools.streamline_merging import StreamlineMerging         # noqa: E402
from tools.streamline_intersection_splitter import (           # noqa: E402
    StreamlineIntersectionSplitter)

_ROOT = Path(__file__).resolve().parent.parent.parent

part = "shroud" if len(sys.argv) > 1 and sys.argv[1] in ("hub", "shroud") else "shroud"
STL = str(_ROOT / f"T1_9_{part}_raw.stl")
OUT = _ROOT / "output" / "T1_9" / f"{part}_stage1" / "xiao_only"


def run_xiao(stl=STL):
    ps.set_periodic(True)
    ps.set_tile_periodic(False)
    mesh, _tf = build_dp_data(stl)
    ff = FrameField(mesh)
    m = detect_singularities(ff.mesh)
    sl = StreamlineGenerator_v2(ff.mesh)
    ps._drop_degenerate_corner_seps(sl.mesh)
    ps._snap_separatrix_endpoints(sl.mesh, radius=0.045)
    merging = StreamlineMerging(sl.mesh, verbose=False)
    splitter = StreamlineIntersectionSplitter(offset_boundingBox=0.05,
                                              num_samples=5)
    updated = splitter.process_streamlines(merging.new_streamlines)
    gen = tmf.TMeshFaceGenerator(updated, blade_loops=list(mesh.blade_loops),
                                 flat_tol_deg=15.0, verbose=False)
    result = gen.get_blocks()
    return {"result": result, "mesh": mesh, "sl": sl}


def _get_singularity_labels(mesh):
    """Map singularity coords to user labels S1..S4.
    User's convention from conversation:
        S1 = (0.6766, 0.6071)   [my S5]
        S2 = (0.6589, 0.6647)   [my S4]
        S3 = (0.3688, 0.4631)   [my S3]
        S4 = (0.3258, 0.5707)   [my S2]
    """
    sings = {}
    for i, d in enumerate(mesh.separatrices):
        sc = d.get("singularity_coords")
        if sc is None:
            continue
        key = (round(float(sc[0]), 6), round(float(sc[1]), 6))
        if key not in sings:
            sings[key] = {"coords": np.asarray(sc, float), "seps": []}
        idx = len(mesh.streamlines) - len(mesh.separatrices) + i
        s = np.asarray(mesh.streamlines[idx], float)
        sings[key]["seps"].append((idx, s[-1].copy()))

    # Map to user labels
    mapping = {
        (0.6766, 0.6071): "S1",
        (0.6589, 0.6647): "S2",
        (0.3688, 0.4631): "S3",
        (0.3258, 0.5707): "S4",
    }
    labeled = {}
    for key, data in sings.items():
        rkey = (round(key[0], 4), round(key[1], 4))
        label = mapping.get(rkey, f"?{rkey}")
        labeled[label] = data
    return labeled


def main():
    us.set_blade_tip_corners(False)
    cs.set_emanate_outer_corners(True)
    OUT.mkdir(parents=True, exist_ok=True)

    data = run_xiao()
    mesh = data["mesh"]
    sl = data["sl"]
    sings = _get_singularity_labels(sl.mesh)
    n_b = len(sl.mesh.streamlines) - len(sl.mesh.separatrices)

    # Highlight indices mentioned by user
    highlight = {20, 22, 23, 25, 27}

    fig, ax = plt.subplots(figsize=(14, 12))
    xy = mesh.x[:, 0:2].numpy()
    tris = mesh.faces.T.numpy()
    ax.triplot(xy[:, 0], xy[:, 1], tris, lw=0.06, color="0.96")

    # Plot all streamlines very faint
    for i, s in enumerate(sl.mesh.streamlines):
        s = np.asarray(s, float)
        if s.ndim != 2 or len(s) < 2:
            continue
        color = "0.85"
        lw = 0.4
        alpha = 0.3
        zorder = 1
        if i in highlight:
            color = "red"
            lw = 3.0
            alpha = 1.0
            zorder = 10
        ax.plot(s[:, 0], s[:, 1], color=color, lw=lw, alpha=alpha, zorder=zorder)
        if i in highlight:
            mid = len(s) // 2
            ax.text(s[mid, 0], s[mid, 1], f"{i}", fontsize=12,
                    color="darkred", ha="center", va="center",
                    fontweight="bold",
                    bbox=dict(boxstyle="round,pad=0.2", facecolor="yellow",
                              edgecolor="darkred", alpha=0.9),
                    zorder=11)

    # Plot singularities
    colors_sing = {"S1": "blue", "S2": "green", "S3": "purple", "S4": "orange"}
    for label, data in sings.items():
        c = data["coords"]
        color = colors_sing.get(label, "gray")
        ax.scatter(c[0], c[1], c=color, s=300, marker="*",
                   zorder=12, edgecolors="black", linewidths=2)
        ax.text(c[0], c[1], f"{label}\n({c[0]:.3f},{c[1]:.3f})", fontsize=11,
                color="black", ha="center", va="bottom", fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                          edgecolor=color, linewidth=2, alpha=0.9),
                zorder=13)

    # Mark the "gap" angle at S2 where user expects another streamline
    if "S2" in sings and "S1" in sings:
        s2 = sings["S2"]["coords"]
        # Direction of sl 23 (S1→S2)
        s23 = np.asarray(sl.mesh.streamlines[23], float)
        v23 = s23[-1] - s23[0]
        v23 = v23 / (np.linalg.norm(v23) + 1e-30)
        # Direction of sl 25 (S2→boundary)
        s25 = np.asarray(sl.mesh.streamlines[25], float)
        v25 = s25[-1] - s25[0]
        v25 = v25 / (np.linalg.norm(v25) + 1e-30)
        # Direction of sl 20 (S1→boundary, ~parallel to expected)
        s20 = np.asarray(sl.mesh.streamlines[20], float)
        v20 = s20[-1] - s20[0]
        v20 = v20 / (np.linalg.norm(v20) + 1e-30)

        # Draw angle arc at S2
        ang23 = np.arctan2(v23[1], v23[0])
        ang25 = np.arctan2(v25[1], v25[0])
        ang20 = np.arctan2(v20[1], v20[0])
        theta1, theta2 = min(ang23, ang25), max(ang23, ang25)
        arc_r = 0.05
        arc_angles = np.linspace(theta1, theta2, 30)
        ax.plot(s2[0] + arc_r * np.cos(arc_angles),
                s2[1] + arc_r * np.sin(arc_angles),
                "magenta", lw=2.5, linestyle="--", zorder=14)
        ax.text(s2[0] + 0.08, s2[1] + 0.03, "GAP?", fontsize=12,
                color="magenta", fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                          edgecolor="magenta", alpha=0.9),
                zorder=15)

        # Draw expected direction (~parallel to sl 20) from S2
        expected_len = 0.15
        ax.annotate("", xy=(s2[0] + expected_len * v20[0], s2[1] + expected_len * v20[1]),
                    xytext=(s2[0], s2[1]),
                    arrowprops=dict(arrowstyle="->", color="magenta", lw=2.5,
                                    linestyle="--"),
                    zorder=14)

    # Legend
    for i in highlight:
        ax.plot([], [], color="red", lw=3.0, label=f"sl {i}")
    ax.legend(loc="upper left", fontsize=10, title="Highlighted streamlines")

    ax.set_aspect("equal")
    ax.set_xlim(0.25, 0.75)
    ax.set_ylim(0.35, 0.75)
    ax.set_title("Zoom: problematic region (sl 20, 22, 23, 25, 27)\n"
                 "S2 gap angle highlighted in magenta")

    out = OUT / "xiao_zoom_problematic.png"
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
