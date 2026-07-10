#!/usr/bin/env python3
"""Run ONLY the Xiao method and create detailed labeled plots.

Usage:
    cd /home/t1dde/Duty/projects/domain_partition/domain_partition_3D
    python scripts/stage1/run_xiao_only.py [hub|shroud]

Outputs to:
    output/T1_9/<part>_stage1/xiao_only/xiao_streamlines_labeled.png
    output/T1_9/<part>_stage1/xiao_only/xiao_singularity_graph.png
    output/T1_9/<part>_stage1/xiao_only/xiao_blocks.png
    output/T1_9/<part>_stage1/xiao_only/xiao_metrics.json
    output/T1_9/<part>_stage1/xiao_only/xiao_streamline_table.txt
"""

import json
import sys
import time
from collections import defaultdict
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

# CRITICAL: partition_surface monkey-patches StreamlineMerging.__init__ to
# call find_missed_streamline_endpoints (Xiao 2020 Alg.2 Case 2). Without
# this import, the merge conflict between streamlines like 22 and 27 is
# NOT caught and the O-grid ring never closes properly.
assert hasattr(StreamlineMerging.__init__, '__wrapped__') or \
       'find_missed_streamline_endpoints' in StreamlineMerging.__init__.__code__.co_names, \
       "StreamlineMerging monkeypatch missing! partition_surface must be imported."

_ROOT = Path(__file__).resolve().parent.parent.parent

part = "shroud" if len(sys.argv) > 1 and sys.argv[1] in ("hub", "shroud") else "shroud"
STL = str(_ROOT / f"T1_9_{part}_raw.stl")
OUT = _ROOT / "output" / "T1_9" / f"{part}_stage1" / "xiao_only"


def run_xiao(stl=STL):
    """Pure Xiao + custom clipping: pipeline up to the snap, then blocks."""
    ps.set_periodic(True)
    ps.set_tile_periodic(False)
    t0 = time.time()
    mesh, _tf = build_dp_data(stl)
    ff = FrameField(mesh)
    m = detect_singularities(ff.mesh)
    n_sing = int((m.singularities != 0).sum())
    sl = StreamlineGenerator_v2(ff.mesh)
    ps._drop_degenerate_corner_seps(sl.mesh)
    ps._snap_separatrix_endpoints(sl.mesh, radius=0.045)

    merging = StreamlineMerging(sl.mesh, verbose=False)
    splitter = StreamlineIntersectionSplitter(offset_boundingBox=0.05,
                                              num_samples=5)
    updated = splitter.process_streamlines(merging.new_streamlines)
    gen = tmf.TMeshFaceGenerator(updated, blade_loops=list(mesh.blade_loops),
                                 flat_tol_deg=15.0, verbose=True)
    result = gen.get_blocks()

    n_boundary = len(sl.mesh.streamlines) - len(sl.mesh.separatrices)
    boundary_ref = [np.asarray(s, float)
                    for s in sl.mesh.streamlines[:n_boundary]]
    _reg, tnodes, irregular = tmf.node_regularity(
        result, lambda p: ps._min_boundary_dist(p, boundary_ref))
    metrics = {
        "approach": "Xiao + custom clipping (no seam postprocessing)",
        "singularities": n_sing,
        "blocks": len(result["blocks"]),
        "nonquad_regions": len(result["rejects"]),
        "irregular_interior_nodes": len(irregular),
        "runtime_s": round(time.time() - t0, 1),
    }
    return {"metrics": metrics, "result": result, "mesh": mesh,
            "sl": sl, "irregular": irregular, "tnodes": tnodes,
            "merging": merging}


def _get_singularities(mesh):
    """Extract unique singularity positions from separatrix metadata."""
    sings = {}
    for i, d in enumerate(mesh.separatrices):
        sc = d.get("singularity_coords")
        if sc is None:
            continue
        key = (round(float(sc[0]), 6), round(float(sc[1]), 6))
        if key not in sings:
            sings[key] = {
                "coords": np.asarray(sc, float),
                "sep_indices": [],
                "sep_targets": [],
            }
        sings[key]["sep_indices"].append(len(mesh.streamlines) - len(mesh.separatrices) + i)
        # Target of this separatrix
        idx = len(mesh.streamlines) - len(mesh.separatrices) + i
        s = np.asarray(mesh.streamlines[idx], float)
        sings[key]["sep_targets"].append((idx, s[-1].copy()))
    # Sort by x-coordinate for consistent S1, S2, S3, S4 labeling
    items = sorted(sings.items(), key=lambda kv: kv[1]["coords"][0])
    labeled = {}
    for li, (key, data) in enumerate(items):
        name = f"S{li+1}"
        labeled[name] = data
        labeled[name]["name"] = name
    return labeled


def _find_target_singularity(target_point, singularities, tol=0.02):
    """Find which singularity (if any) the target point corresponds to."""
    for name, data in singularities.items():
        if np.linalg.norm(target_point - data["coords"]) < tol:
            return name
    return None


def _plot_streamlines_unlabeled(ax, mesh, singularities, highlight_indices=None):
    """Plot streamlines without index labels, only colors."""
    xy = mesh.x[:, 0:2].numpy()
    tris = mesh.faces.T.numpy()
    n_b = len(mesh.streamlines) - len(mesh.separatrices)

    ax.triplot(xy[:, 0], xy[:, 1], tris, lw=0.06, color="0.96")
    cmap = plt.cm.tab20

    for i, s in enumerate(mesh.streamlines):
        s = np.asarray(s, float)
        if s.ndim != 2 or len(s) < 2:
            continue
        color = cmap(i % 20 / 20.0)
        lw = 2.0 if i >= n_b else 0.7
        alpha = 1.0 if i >= n_b else 0.5
        zorder = 5 if i >= n_b else 2
        if highlight_indices and i in highlight_indices:
            color = "red"
            lw = 3.0
            zorder = 10
        ax.plot(s[:, 0], s[:, 1], color=color, lw=lw, alpha=alpha,
                zorder=zorder)

    sing_faces = mesh.singularities.numpy() != 0
    if sing_faces.any():
        sc = xy[tris[sing_faces]].mean(axis=1)
        ax.scatter(sc[:, 0], sc[:, 1], c="blue", s=120, marker="*",
                   zorder=8, edgecolors="navy", linewidths=1)

    ax.plot([], [], "-", color="0.5", lw=0.7, alpha=0.5, label="boundary")
    ax.plot([], [], "-", color="0.3", lw=2.0, label="separatrix")
    ax.set_aspect("equal")
    ax.legend(loc="upper left", fontsize=8)


def _plot_streamlines_labeled(ax, mesh, singularities, highlight_indices=None):
    """Plot every streamline in a unique color with numeric labels."""
    xy = mesh.x[:, 0:2].numpy()
    tris = mesh.faces.T.numpy()
    n_b = len(mesh.streamlines) - len(mesh.separatrices)

    # Triangular background mesh (very faint)
    ax.triplot(xy[:, 0], xy[:, 1], tris, lw=0.06, color="0.96")

    # Color map for streamlines
    cmap = plt.cm.tab20

    for i, s in enumerate(mesh.streamlines):
        s = np.asarray(s, float)
        if s.ndim != 2 or len(s) < 2:
            continue
        color = cmap(i % 20 / 20.0)
        lw = 2.0 if i >= n_b else 0.7
        alpha = 1.0 if i >= n_b else 0.5
        zorder = 5 if i >= n_b else 2
        ls = "-"
        if highlight_indices and i in highlight_indices:
            color = "red"
            lw = 3.0
            alpha = 1.0
            zorder = 10
        ax.plot(s[:, 0], s[:, 1], color=color, lw=lw, alpha=alpha,
                ls=ls, zorder=zorder)

        # Label at midpoint
        mid = len(s) // 2
        ax.text(s[mid, 0], s[mid, 1], str(i), fontsize=7,
                color="black" if i >= n_b else "0.4",
                ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.15", facecolor="white",
                          edgecolor="none", alpha=0.7),
                zorder=zorder + 1)

    # Singularities
    sing_faces = mesh.singularities.numpy() != 0
    if sing_faces.any():
        sc = xy[tris[sing_faces]].mean(axis=1)
        ax.scatter(sc[:, 0], sc[:, 1], c="blue", s=120, marker="*",
                   zorder=8, edgecolors="navy", linewidths=1,
                   label=f"singularities: {len(sc)}")

    # Singularities with S1, S2, S3, S4 labels and coordinates
    for name, data in singularities.items():
        c = data["coords"]
        ax.text(c[0], c[1], f"{name}\n({c[0]:.3f},{c[1]:.3f})", fontsize=9,
                color="navy", ha="center", va="bottom", fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.2", facecolor="lightyellow",
                          edgecolor="navy", alpha=0.8),
                zorder=9)

    # Separatrix/boundary legend helpers
    ax.plot([], [], "-", color="0.5", lw=0.7, alpha=0.5, label="boundary")
    ax.plot([], [], "-", color="0.3", lw=2.0, label="separatrix")
    ax.set_aspect("equal")
    ax.legend(loc="upper left", fontsize=8)


def _plot_merged_streamlines(ax, merged_streamlines, singularities):
    n_total = len(merged_streamlines)
    n_boundary = min(36, n_total)
    n_sep = n_total - n_boundary

    cmap = plt.cm.tab20
    for i, coords in enumerate(merged_streamlines):
        coords = np.asarray(coords, float)
        if len(coords) < 2:
            continue
        color = cmap(i % 20 / 20.0)
        lw = 2.0 if i >= n_boundary else 0.7
        alpha = 1.0 if i >= n_boundary else 0.5
        zorder = 5 if i >= n_boundary else 2
        ax.plot(coords[:, 0], coords[:, 1], color=color, lw=lw, alpha=alpha, zorder=zorder)

        mid = len(coords) // 2
        ax.text(coords[mid, 0], coords[mid, 1], str(i), fontsize=7,
                color="black" if i >= n_boundary else "0.4",
                ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.15", facecolor="white",
                          edgecolor="none", alpha=0.7),
                zorder=zorder + 1)

    for name, data in singularities.items():
        c = data["coords"]
        ax.scatter(c[0], c[1], c="blue", s=120, marker="*",
                   zorder=8, edgecolors="navy", linewidths=1)
        ax.text(c[0], c[1], f"{name}\n({c[0]:.3f},{c[1]:.3f})", fontsize=9,
                color="navy", ha="center", va="bottom", fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.2", facecolor="lightyellow",
                          edgecolor="navy", alpha=0.8),
                zorder=9)

    ax.plot([], [], "-", color="0.5", lw=0.7, alpha=0.5, label="boundary")
    ax.plot([], [], "-", color="0.3", lw=2.0, label="separatrix")
    ax.set_aspect("equal")
    ax.legend(loc="upper left", fontsize=8)


def _plot_singularity_graph(ax, mesh, singularities):
    """Plot only the singularity-to-singularity connections with distinct colors."""
    xy = mesh.x[:, 0:2].numpy()
    tris = mesh.faces.T.numpy()
    n_b = len(mesh.streamlines) - len(mesh.separatrices)

    # Faint background
    ax.triplot(xy[:, 0], xy[:, 1], tris, lw=0.04, color="0.97")

    # Colors for singularity-to-singularity connections
    connection_colors = {
        ("S1", "S2"): "red",
        ("S1", "S3"): "blue",
        ("S1", "S4"): "green",
        ("S2", "S3"): "purple",
        ("S2", "S4"): "orange",
        ("S3", "S4"): "brown",
    }

    # Track connections
    connections = defaultdict(list)
    connection_sl_indices = {}

    for i, d in enumerate(mesh.separatrices):
        idx = n_b + i
        s = np.asarray(mesh.streamlines[idx], float)
        if s.ndim != 2 or len(s) < 2:
            continue
        sc = d.get("singularity_coords")
        if sc is None:
            continue
        # Find source singularity
        src = None
        for name, data in singularities.items():
            if np.linalg.norm(data["coords"] - np.asarray(sc, float)) < 0.01:
                src = name
                break
        # Find target singularity (or boundary)
        tgt = _find_target_singularity(s[-1], singularities)
        if src and tgt and src != tgt:
            key = tuple(sorted([src, tgt]))
            connections[key].append(idx)
            connection_sl_indices[idx] = key

    # Plot connections
    for key, indices in connections.items():
        color = connection_colors.get(key, "gray")
        for idx in indices:
            s = np.asarray(mesh.streamlines[idx], float)
            ax.plot(s[:, 0], s[:, 1], color=color, lw=2.5, alpha=0.8, zorder=5)

    # Plot all streamlines very faint
    for i, s in enumerate(mesh.streamlines):
        s = np.asarray(s, float)
        if s.ndim != 2 or len(s) < 2:
            continue
        ax.plot(s[:, 0], s[:, 1], "0.7", lw=0.3, alpha=0.3, zorder=1)

    # Singularities
    for name, data in singularities.items():
        c = data["coords"]
        ax.scatter(c[0], c[1], c="blue", s=200, marker="*",
                   zorder=8, edgecolors="navy", linewidths=2)
        ax.text(c[0], c[1], f"{name}\n({c[0]:.3f},{c[1]:.3f})", fontsize=10,
                color="navy", ha="center", va="bottom", fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.2", facecolor="lightyellow",
                          edgecolor="navy", alpha=0.8),
                zorder=9)

    # Legend for connections
    for key, color in connection_colors.items():
        if key in connections:
            ax.plot([], [], color=color, lw=2.5,
                    label=f"{'↔'.join(key)} ({len(connections[key])} sl)")
    ax.set_aspect("equal")
    ax.legend(loc="upper left", fontsize=8, title="Singularity connections")

    return connections


def _plot_blocks(ax, result, mesh, irregular, title=""):
    """Plot block structure with irregular nodes highlighted."""
    xy = mesh.x[:, 0:2].numpy()
    tris = mesh.faces.T.numpy()
    nodes = result["nodes"]
    ax.triplot(xy[:, 0], xy[:, 1], tris, lw=0.08, color="0.94")
    for blk in result["blocks"]:
        ring = blk["ring"]
        ax.fill(ring[:, 0], ring[:, 1], alpha=0.2, color="C0")
        for s in blk["sides"]:
            ax.plot(s[:, 0], s[:, 1], "C0", lw=1.1)
        c = nodes[blk["corners"]]
        ax.scatter(c[:, 0], c[:, 1], c="k", s=12, zorder=6)
    for ri, rej in enumerate(result["rejects"]):
        if rej.get("ring") is not None:
            r = rej["ring"]
            ax.fill(r[:, 0], r[:, 1], color="red", alpha=0.3)
            ax.plot(r[:, 0], r[:, 1], "red", lw=1.8,
                    label="non-quad region" if ri == 0 else None)
            cen = r[:-1].mean(axis=0)
            ax.annotate(f"{rej['n_real']}-corner", cen, fontsize=8,
                        color="darkred", ha="center")
    if irregular:
        ax.scatter(nodes[irregular, 0], nodes[irregular, 1],
                   facecolors="none", edgecolors="red", s=170, linewidths=2,
                   zorder=8, label=f"irregular interior: {len(irregular)}")
    ax.set_aspect("equal")
    handles, labels = ax.get_legend_handles_labels()
    if labels:
        ax.legend(loc="upper left", fontsize=8)
    ax.set_title(title, fontsize=10)


def main():
    us.set_blade_tip_corners(False)
    cs.set_emanate_outer_corners(True)
    OUT.mkdir(parents=True, exist_ok=True)

    print(f"=== Running Xiao method on {part} ===")
    data = run_xiao()
    mesh = data["mesh"]
    sl = data["sl"]
    result = data["result"]
    irregular = data["irregular"]
    metrics = data["metrics"]

    singularities = _get_singularities(sl.mesh)
    n_b = len(sl.mesh.streamlines) - len(sl.mesh.separatrices)

    merging = data["merging"]

    # --- Unlabeled streamline plot (nur Farben, keine Zahlen) ---
    fig, ax = plt.subplots(figsize=(14, 10))
    _plot_streamlines_unlabeled(ax, sl.mesh, singularities)
    ax.set_title(
        f"Xiao method: {n_b} boundary streamlines, "
        f"{len(sl.mesh.separatrices)} separatrices, "
        f"{metrics['singularities']} singularities\n"
        f"(no labels; bold = separatrix)")
    out = OUT / "xiao_streamlines_unlabeled.png"
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")

    # --- Detailed streamline plot (farbcodiert + gelabelt) ---
    fig, ax = plt.subplots(figsize=(16, 12))
    _plot_streamlines_labeled(ax, sl.mesh, singularities)
    ax.set_title(
        f"Xiao method: {n_b} boundary streamlines, "
        f"{len(sl.mesh.separatrices)} separatrices, "
        f"{metrics['singularities']} singularities\n"
        f"(labels = streamline index; bold = separatrix)")
    out = OUT / "xiao_streamlines_labeled.png"
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")

    # --- Post-merge streamline plot (shows actual merged result) ---
    fig, ax = plt.subplots(figsize=(16, 12))
    _plot_merged_streamlines(ax, merging.new_streamlines, singularities)
    ax.set_title(
        f"Xiao method (after merge): {metrics['singularities']} singularities, "
        f"{len(merging.new_streamlines)} streamlines post-merge")
    out = OUT / "xiao_streamlines_merged.png"
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")

    # --- Singularity connection graph ---
    fig, ax = plt.subplots(figsize=(14, 11))
    connections = _plot_singularity_graph(ax, sl.mesh, singularities)
    ax.set_title(
        f"Singularity connection graph: {len(connections)} direct connections\n"
        f"(thick colored = sing↔sing; faint gray = all streamlines)")
    out = OUT / "xiao_singularity_graph.png"
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")

    # --- Block partition plot ---
    fig, ax = plt.subplots(figsize=(12, 10))
    _plot_blocks(ax, result, mesh, irregular,
                 title=f"Xiao method: {metrics['blocks']} blocks, "
                       f"{metrics['nonquad_regions']} non-quad regions, "
                       f"{len(irregular)} irregular nodes")
    out = OUT / "xiao_blocks.png"
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")

    # --- Streamline table for inspection ---
    lines = ["Streamline index table", "=" * 60, ""]
    lines.append(f"BOUNDARY streamlines (indices 0 .. {n_b-1}):")
    for i in range(n_b):
        s = np.asarray(sl.mesh.streamlines[i], float)
        lines.append(f"  {i:3d}:  len={len(s):4d}  "
                     f"start=({s[0,0]:.4f},{s[0,1]:.4f})  "
                     f"end=({s[-1,0]:.4f},{s[-1,1]:.4f})")
    lines.append("")
    lines.append(f"SEPARATRICES (indices {n_b} .. {len(sl.mesh.streamlines)-1}):")
    for i, d in enumerate(sl.mesh.separatrices):
        idx = n_b + i
        s = np.asarray(sl.mesh.streamlines[idx], float)
        sing = d.get("singularity_coords")
        sing_str = f"S=({sing[0]:.4f},{sing[1]:.4f})" if sing is not None else "S=?"
        lines.append(f"  {idx:3d}:  len={len(s):4d}  {sing_str}  "
                     f"start=({s[0,0]:.4f},{s[0,1]:.4f})  "
                     f"end=({s[-1,0]:.4f},{s[-1,1]:.4f})")
    txt = "\n".join(lines)
    (OUT / "xiao_streamline_table.txt").write_text(txt)
    print(txt)
    print(f"wrote {OUT}/xiao_streamline_table.txt")

    # --- Singularity analysis ---
    lines = ["Singularity analysis", "=" * 60, ""]
    for name, data in singularities.items():
        c = data["coords"]
        lines.append(f"{name}: ({c[0]:.4f}, {c[1]:.4f})")
        lines.append(f"  Separatrices ({len(data['sep_indices'])} total):")
        for idx, tgt in data['sep_targets']:
            tgt_name = _find_target_singularity(tgt, singularities)
            if tgt_name:
                lines.append(f"    sl {idx}: → {tgt_name} (direct sing↔sing)")
            else:
                lines.append(f"    sl {idx}: → boundary ({tgt[0]:.4f},{tgt[1]:.4f})")
        lines.append("")

    # Adjacency matrix
    lines.append("Singularity adjacency (direct connections):")
    names = [f"S{i}" for i in range(1, len(singularities) + 1)]
    for a in names:
        for b in names:
            if a >= b:
                continue
            key = tuple(sorted([a, b]))
            if key in connections:
                lines.append(f"  {a} ↔ {b}: YES ({len(connections[key])} streamlines)")
                for idx in connections[key]:
                    s = np.asarray(sl.mesh.streamlines[idx], float)
                    lines.append(f"      sl {idx}: len={len(s)}")
            else:
                lines.append(f"  {a} ↔ {b}: NO DIRECT CONNECTION")
        lines.append("")

    txt2 = "\n".join(lines)
    (OUT / "xiao_singularity_analysis.txt").write_text(txt2)
    print(txt2)
    print(f"wrote {OUT}/xiao_singularity_analysis.txt")

    # --- Metrics JSON ---
    (OUT / "xiao_metrics.json").write_text(json.dumps(metrics, indent=2))
    print(f"wrote {OUT}/xiao_metrics.json")
    print(f"\nMetrics: {metrics}")


if __name__ == "__main__":
    main()
