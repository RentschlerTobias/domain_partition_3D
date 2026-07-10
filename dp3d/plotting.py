"""Plots for the dp3d pipeline.

Two families:
- analysis plots called from tmesh.run_tmesh / xiao (block structure, TFI
  grid, tiled passages, streamline states, singularity graph)
- showcase step plots (step01..step09): the presentation-style walkthrough
  of every pipeline stage, written by run_showcase() when --plots is set
"""

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from collections import defaultdict
from pathlib import Path

from . import partition_surface as ps
from . import tmesh_faces as tmf


# --------------------------------------------------------------------------
# analysis plots (T-a/T-b pipeline)
# --------------------------------------------------------------------------

def plot_blocks(result, mesh, tnodes, irregular, seam_info, out_png, tag=""):
    xy = mesh.x[:, 0:2].numpy()
    tris = mesh.faces.T.numpy()
    nodes = result["nodes"]
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.triplot(xy[:, 0], xy[:, 1], tris, lw=0.1, color="0.92")
    for blk in result["blocks"]:
        ring = blk["ring"]
        ax.fill(ring[:, 0], ring[:, 1], alpha=0.22, color="0.9")
        for s in blk["sides"]:
            ax.plot(s[:, 0], s[:, 1], "k", lw=1.2)
        c = nodes[blk["corners"]]
        ax.scatter(c[:, 0], c[:, 1], c="k", s=14, zorder=6)
    for rej in result["rejects"]:
        if rej.get("ring") is not None:
            r = rej["ring"]
            ax.plot(r[:, 0], r[:, 1], "red", lw=1.8)
            ax.fill(r[:, 0], r[:, 1], color="red", alpha=0.25)
    if tnodes:
        ax.scatter(nodes[tnodes, 0], nodes[tnodes, 1], marker="s",
                   facecolors="none", edgecolors="darkorange", s=70,
                   linewidths=1.6, zorder=7,
                   label=f"T-node (hanging): {len(tnodes)}")
    if irregular:
        ax.scatter(nodes[irregular, 0], nodes[irregular, 1],
                   facecolors="none", edgecolors="red", s=150, linewidths=2.0,
                   zorder=7, label=f"irregular interior: {len(irregular)}")
    if seam_info:
        for W in (seam_info["WL"], seam_info["WR"]):
            ax.plot(W[:, 0], W[:, 1], "green", lw=0.8, alpha=0.7)
        hang = [c for c in seam_info["canon"]
                if c["srcL"] is None or c["srcR"] is None]
        if hang:
            P = []
            for c in hang:
                side = "L" if c["srcL"] is None else "R"
                P.append(c["coord"] if side == "L"
                         else c["coord"] + np.array([seam_info["pitch"], 0.0]))
            P = np.array(P)
            ax.scatter(P[:, 0], P[:, 1], marker="s", facecolors="none",
                       edgecolors="purple", s=90, linewidths=1.6, zorder=7,
                       label=f"hanging seam junction: {len(P)}")
    ax.set_aspect("equal")
    ax.axis("off")
    fig.savefig(out_png, dpi=200, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"wrote {out_png}")


def plot_tfi(tfi, result, pitch, out_png, tag=""):
    fig, ax = plt.subplots(figsize=(10, 8))
    for X in tfi["grids"]:
        for i in range(X.shape[0]):
            ax.plot(X[i, :, 0], X[i, :, 1], "0.4", lw=0.3)
        for j in range(X.shape[1]):
            ax.plot(X[:, j, 0], X[:, j, 1], "0.4", lw=0.3)
    for blk in result["blocks"]:
        for s in blk["sides"]:
            ax.plot(s[:, 0], s[:, 1], "black", lw=1.0)
    ax.set_aspect("equal")
    ax.set_axis_off()
    fig.savefig(out_png, dpi=200, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"wrote {out_png}")


def plot_tiled(result, tfi, seam_info, pitch, out_blocks, out_tfi, tag=""):
    """3 passages side by side (shifts -pitch, 0, +pitch): visual periodicity
    check -- at the interior seams of the trio the grid points of neighbouring
    copies must coincide."""
    shifts = [-pitch, 0.0, pitch]
    cols = ["0.55", "black", "0.55"]

    fig, ax = plt.subplots(figsize=(18, 7))
    for sh, col in zip(shifts, cols):
        for blk in result["blocks"]:
            for s in blk["sides"]:
                ax.plot(s[:, 0] + sh, s[:, 1], col, lw=1.1)
    ax.set_aspect("equal")
    ax.set_axis_off()
    fig.savefig(out_blocks, dpi=200, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"wrote {out_blocks}")

    fig, ax = plt.subplots(figsize=(18, 7))
    for sh, col in zip(shifts, cols):
        for X in tfi["grids"]:
            for i in range(X.shape[0]):
                ax.plot(X[i, :, 0] + sh, X[i, :, 1], col, lw=0.25)
            for j in range(X.shape[1]):
                ax.plot(X[:, j, 0] + sh, X[:, j, 1], col, lw=0.25)
    ax.set_aspect("equal")
    ax.set_axis_off()
    fig.savefig(out_tfi, dpi=200, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"wrote {out_tfi}")


def plot_streamlines_clean(mesh, out_png, tag=""):
    """Presentation plot of the current streamline/separatrix state: pure
    image, no title/axes/legend, transparent background. Mesh triangulation
    light gray, boundary curves black, separatrices light blue, singularities
    red circles (+index) / blue squares (-index)."""
    xy = mesh.x[:, 0:2].numpy()
    tris = mesh.faces.T.numpy()
    n_b = len(mesh.streamlines) - len(mesh.separatrices)

    fig, ax = plt.subplots(figsize=(4, 6))
    ax.triplot(xy[:, 0], xy[:, 1], tris, lw=0.1, color="0.93")
    for i, s in enumerate(mesh.streamlines):
        s = np.asarray(s, float)
        if s.ndim != 2 or len(s) < 2:
            continue
        if i < n_b:
            ax.plot(s[:, 0], s[:, 1], color="black", lw=1.5)
        else:
            ax.plot(s[:, 0], s[:, 1], color="#87CEEB", lw=1.3)

    sing = mesh.singularities.numpy()
    sfaces = np.where(sing != 0)[0]
    if len(sfaces):
        sc = xy[tris[sfaces]].mean(axis=1)
        pos = sing[sfaces] > 0
        ax.scatter(sc[pos, 0], sc[pos, 1], c="red", s=60, zorder=6)
        ax.scatter(sc[~pos, 0], sc[~pos, 1], c="blue", s=60, zorder=6,
                   marker="s")

    ax.set_aspect("equal")
    ax.set_axis_off()
    fig.savefig(out_png, dpi=200, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"wrote {out_png}")


# --------------------------------------------------------------------------
# xiao analysis plots
# --------------------------------------------------------------------------

def get_singularities(mesh):
    """Unique singularity positions from separatrix metadata, labeled S1..Sn
    by ascending x."""
    sings = {}
    n_b = len(mesh.streamlines) - len(mesh.separatrices)
    for i, d in enumerate(mesh.separatrices):
        sc = d.get("singularity_coords")
        if sc is None:
            continue
        key = (round(float(sc[0]), 6), round(float(sc[1]), 6))
        if key not in sings:
            sings[key] = {"coords": np.asarray(sc, float),
                          "sep_indices": [], "sep_targets": []}
        idx = n_b + i
        s = np.asarray(mesh.streamlines[idx], float)
        sings[key]["sep_indices"].append(idx)
        sings[key]["sep_targets"].append((idx, s[-1].copy()))
    items = sorted(sings.items(), key=lambda kv: kv[1]["coords"][0])
    labeled = {}
    for li, (key, data) in enumerate(items):
        name = f"S{li + 1}"
        labeled[name] = data
        labeled[name]["name"] = name
    return labeled


def find_target_singularity(target_point, singularities, tol=0.02):
    for name, data in singularities.items():
        if np.linalg.norm(target_point - data["coords"]) < tol:
            return name
    return None


def _label_singularities(ax, singularities):
    for name, data in singularities.items():
        c = data["coords"]
        ax.text(c[0], c[1], f"{name}\n({c[0]:.3f},{c[1]:.3f})", fontsize=9,
                color="navy", ha="center", va="bottom", fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.2", facecolor="lightyellow",
                          edgecolor="navy", alpha=0.8),
                zorder=9)


def plot_streamlines_indexed(ax, mesh, singularities, labels=True):
    """All streamlines in distinct colors; optional per-curve index labels."""
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
        ax.plot(s[:, 0], s[:, 1], color=color, lw=lw, alpha=alpha,
                zorder=zorder)
        if labels:
            mid = len(s) // 2
            ax.text(s[mid, 0], s[mid, 1], str(i), fontsize=7,
                    color="black" if i >= n_b else "0.4",
                    ha="center", va="center",
                    bbox=dict(boxstyle="round,pad=0.15", facecolor="white",
                              edgecolor="none", alpha=0.7),
                    zorder=zorder + 1)

    sing_faces = mesh.singularities.numpy() != 0
    if sing_faces.any():
        sc = xy[tris[sing_faces]].mean(axis=1)
        ax.scatter(sc[:, 0], sc[:, 1], c="blue", s=120, marker="*",
                   zorder=8, edgecolors="navy", linewidths=1,
                   label=f"singularities: {len(sc)}")
    if labels:
        _label_singularities(ax, singularities)

    ax.plot([], [], "-", color="0.5", lw=0.7, alpha=0.5, label="boundary")
    ax.plot([], [], "-", color="0.3", lw=2.0, label="separatrix")
    ax.set_aspect("equal")
    ax.legend(loc="upper left", fontsize=8)


def plot_merged_streamlines(ax, merged_streamlines, singularities,
                            n_boundary_max=36):
    n_total = len(merged_streamlines)
    n_boundary = min(n_boundary_max, n_total)

    cmap = plt.cm.tab20
    for i, coords in enumerate(merged_streamlines):
        coords = np.asarray(coords, float)
        if len(coords) < 2:
            continue
        color = cmap(i % 20 / 20.0)
        lw = 2.0 if i >= n_boundary else 0.7
        alpha = 1.0 if i >= n_boundary else 0.5
        zorder = 5 if i >= n_boundary else 2
        ax.plot(coords[:, 0], coords[:, 1], color=color, lw=lw, alpha=alpha,
                zorder=zorder)
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
    _label_singularities(ax, singularities)

    ax.plot([], [], "-", color="0.5", lw=0.7, alpha=0.5, label="boundary")
    ax.plot([], [], "-", color="0.3", lw=2.0, label="separatrix")
    ax.set_aspect("equal")
    ax.legend(loc="upper left", fontsize=8)


def plot_singularity_graph(ax, mesh, singularities):
    """Singularity-to-singularity separatrix connections in distinct colors.
    Returns the {(Sa, Sb): [streamline indices]} connection map."""
    xy = mesh.x[:, 0:2].numpy()
    tris = mesh.faces.T.numpy()
    n_b = len(mesh.streamlines) - len(mesh.separatrices)

    ax.triplot(xy[:, 0], xy[:, 1], tris, lw=0.04, color="0.97")

    connection_colors = {
        ("S1", "S2"): "red",
        ("S1", "S3"): "blue",
        ("S1", "S4"): "green",
        ("S2", "S3"): "purple",
        ("S2", "S4"): "orange",
        ("S3", "S4"): "brown",
    }

    connections = defaultdict(list)
    for i, d in enumerate(mesh.separatrices):
        idx = n_b + i
        s = np.asarray(mesh.streamlines[idx], float)
        if s.ndim != 2 or len(s) < 2:
            continue
        sc = d.get("singularity_coords")
        if sc is None:
            continue
        src = None
        for name, data in singularities.items():
            if np.linalg.norm(data["coords"] - np.asarray(sc, float)) < 0.01:
                src = name
                break
        tgt = find_target_singularity(s[-1], singularities)
        if src and tgt and src != tgt:
            connections[tuple(sorted([src, tgt]))].append(idx)

    for key, indices in connections.items():
        color = connection_colors.get(key, "gray")
        for idx in indices:
            s = np.asarray(mesh.streamlines[idx], float)
            ax.plot(s[:, 0], s[:, 1], color=color, lw=2.5, alpha=0.8, zorder=5)

    for i, s in enumerate(mesh.streamlines):
        s = np.asarray(s, float)
        if s.ndim != 2 or len(s) < 2:
            continue
        ax.plot(s[:, 0], s[:, 1], "0.7", lw=0.3, alpha=0.3, zorder=1)

    for name, data in singularities.items():
        c = data["coords"]
        ax.scatter(c[0], c[1], c="blue", s=200, marker="*",
                   zorder=8, edgecolors="navy", linewidths=2)
    _label_singularities(ax, singularities)

    for key, color in connection_colors.items():
        if key in connections:
            ax.plot([], [], color=color, lw=2.5,
                    label=f"{'-'.join(key)} ({len(connections[key])} sl)")
    ax.set_aspect("equal")
    ax.legend(loc="upper left", fontsize=8, title="Singularity connections")
    return connections


def plot_xiao_blocks(ax, result, mesh, irregular, title=""):
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


# --------------------------------------------------------------------------
# showcase steps (presentation-style walkthrough, --plots)
# --------------------------------------------------------------------------

# Left seam of the hub runs (0,0)->(0.407,1): s_seam(t) = SEAM_SLOPE * t.
SEAM_SLOPE = 0.407


def fold_curve(s, pitch, seam_slope=SEAM_SLOPE):
    """Fold a cover-coordinate polyline into the fundamental pitch, splitting
    it wherever it crosses a seam. The crossing point is interpolated exactly
    and appended to BOTH adjacent segments, so every winding is drawn
    seam-to-seam with no visual gap."""
    s = np.asarray(s, float)
    u = s[:, 0] - seam_slope * s[:, 1]
    k = np.floor(u / pitch).astype(int)
    segs, cur = [], [s[0]]
    for i in range(1, len(s)):
        if k[i] != k[i - 1]:
            ub = pitch * max(k[i], k[i - 1])
            f = (ub - u[i - 1]) / (u[i] - u[i - 1] + 1e-30)
            pc = s[i - 1] + np.clip(f, 0.0, 1.0) * (s[i] - s[i - 1])
            cur.append(pc)
            segs.append((k[i - 1], np.asarray(cur)))
            cur = [pc, s[i]]
        else:
            cur.append(s[i])
    segs.append((k[-1], np.asarray(cur)))
    out = []
    for kk, seg in segs:
        seg = seg.copy()
        seg[:, 0] -= kk * pitch
        if len(seg) >= 2:
            out.append(seg)
    return out


def _save(fig, out_png):
    fig.savefig(out_png, dpi=200, transparent=True, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_png}")


def _draw_domain(ax, st, tris, loops):
    """Shared showcase background: faint triangulation + black contours."""
    ax.triplot(st[:, 0], st[:, 1], tris, lw=0.1, color="0.93")
    for loop in loops:
        ring = st[loop + [loop[0]]]
        ax.plot(ring[:, 0], ring[:, 1], lw=2.0, color="black")


def showcase_3d_surface(stl_path, out_png):
    """step01: 3D input surface."""
    import meshio
    m = meshio.read(str(stl_path))
    points = m.points
    triangles = m.cells[0].data

    fig = plt.figure(figsize=(4, 6))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_trisurf(points[:, 0], points[:, 1], points[:, 2],
                    triangles=triangles, color="lightgray", alpha=0.8,
                    edgecolor="none")
    ax.set_axis_off()
    spans = points.max(axis=0) - points.min(axis=0)
    ax.set_box_aspect((spans[0], spans[1], spans[2]))
    _save(fig, out_png)


def showcase_unwrapped(st, tris, loops, out_png):
    """step02: unwrapped (s, t) domain."""
    fig, ax = plt.subplots(figsize=(4, 6))
    _draw_domain(ax, st, tris, loops)
    ax.set_aspect("equal")
    ax.set_axis_off()
    _save(fig, out_png)


def showcase_crossfield(mesh, st, tris, loops, out_png):
    """step03: 4-RoSy cross-field on boundary nodes."""
    import torch
    boundary_ids = torch.unique(
        mesh.edge_index[0, mesh.edge_attr == 1]).numpy()
    u = mesh.u.numpy()
    base = np.arctan2(u[:, 1], u[:, 0]) / 4.0
    bxy = st[boundary_ids]

    fig, ax = plt.subplots(figsize=(4, 6))
    _draw_domain(ax, st, tris, loops)
    for k in range(2):
        ang = base[boundary_ids] + k * (np.pi / 2)
        ax.quiver(bxy[:, 0], bxy[:, 1], np.cos(ang), np.sin(ang),
                  color="#3aa6d0", scale=18, width=0.004,
                  headwidth=0, headlength=0, pivot="mid", alpha=0.95)
    ax.set_aspect("equal")
    ax.set_axis_off()
    _save(fig, out_png)


def showcase_representatives(mesh, st, tris, loops, out_field, out_concept):
    """step04: cross-field + representative vectors (field view) and the
    conceptual cross -> 4*theta representative mapping."""
    import torch
    boundary_ids = torch.unique(
        mesh.edge_index[0, mesh.edge_attr == 1]).numpy()
    u = mesh.u.numpy()
    base = np.arctan2(u[:, 1], u[:, 0]) / 4.0
    bxy = st[boundary_ids]

    fig, ax = plt.subplots(figsize=(4, 6))
    _draw_domain(ax, st, tris, loops)
    for k in range(4):
        ang = base[boundary_ids] + k * (np.pi / 2)
        ax.quiver(bxy[:, 0], bxy[:, 1], np.cos(ang), np.sin(ang),
                  color="#3aa6d0", scale=18, width=0.004,
                  headwidth=0, headlength=0, pivot="mid", alpha=0.85)
    ax.quiver(bxy[:, 0], bxy[:, 1],
              np.cos(base[boundary_ids]), np.sin(base[boundary_ids]),
              color="red", scale=18, width=0.005,
              headwidth=4, headlength=5, pivot="tail", alpha=1.0)
    ax.set_aspect("equal")
    ax.set_axis_off()
    _save(fig, out_field)

    # conceptual: cross at theta, representative vector at 4*theta
    theta = np.deg2rad(20.0)
    rep = 4.0 * theta
    fig, ax = plt.subplots(figsize=(4, 6))
    ax.set_xlim(-1.4, 1.4)
    ax.set_ylim(-1.4, 1.4)
    ax.set_aspect("equal")
    ax.axhline(0.0, color="0.75", lw=1.0, zorder=1)
    for k in range(4):
        ang = theta + k * (np.pi / 2)
        ax.quiver(0.0, 0.0, np.cos(ang), np.sin(ang),
                  color="#1f77ff", scale=2.4, width=0.012,
                  headwidth=4, headlength=5, pivot="tail", zorder=3,
                  alpha=0.95)
    ax.quiver(0.0, 0.0, np.cos(rep), np.sin(rep),
              color="red", scale=2.4, width=0.014,
              headwidth=4, headlength=5, pivot="tail", zorder=4)
    ax.text(0.62 * np.cos(theta / 2), 0.62 * np.sin(theta / 2),
            r"$\theta$", color="0.35", fontsize=12, ha="center", va="center")
    ax.text(np.cos(rep) - 0.22, np.sin(rep) + 0.02,
            r"$4\theta$", color="red", fontsize=12, ha="right", va="center")
    ax.set_axis_off()
    _save(fig, out_concept)


def showcase_framefield(mesh, st, tris, loops, out_vectors, out_crosses):
    """step05: smoothed representation vector field and the back-mapped
    4-RoSy frame field on every node."""
    u = mesh.u.numpy()
    base = np.arctan2(u[:, 1], u[:, 0]) / 4.0

    fig, ax = plt.subplots(figsize=(4, 6))
    _draw_domain(ax, st, tris, loops)
    ax.quiver(st[:, 0], st[:, 1], u[:, 0], u[:, 1],
              color="#87CEEB", scale=26, width=0.0022,
              headwidth=0, headlength=0, pivot="mid", alpha=0.9)
    ax.set_aspect("equal")
    ax.set_axis_off()
    _save(fig, out_vectors)

    fig, ax = plt.subplots(figsize=(4, 6))
    _draw_domain(ax, st, tris, loops)
    for k in range(4):
        ang = base + k * (np.pi / 2)
        ax.quiver(st[:, 0], st[:, 1], np.cos(ang), np.sin(ang),
                  color="#87CEEB", scale=26, width=0.0022,
                  headwidth=0, headlength=0, pivot="mid", alpha=0.9)
    ax.set_aspect("equal")
    ax.set_axis_off()
    _save(fig, out_crosses)


def showcase_singularities(mesh, st, tris, loops, final_sing_st, out_png):
    """step06: frame field with the FINAL partition singularities marked
    (red rings). ``final_sing_st`` are positions already in raw (s, t)."""
    u = mesh.u.numpy()
    base = np.arctan2(u[:, 1], u[:, 0]) / 4.0

    fig, ax = plt.subplots(figsize=(4, 6))
    _draw_domain(ax, st, tris, loops)
    for k in range(4):
        ang = base + k * (np.pi / 2)
        ax.quiver(st[:, 0], st[:, 1], np.cos(ang), np.sin(ang),
                  color="#87CEEB", scale=40, width=0.0015,
                  headwidth=0, headlength=0, pivot="mid", alpha=0.9)
    if len(final_sing_st):
        ax.scatter(final_sing_st[:, 0], final_sing_st[:, 1],
                   facecolors="none", edgecolors="red", s=150,
                   linewidths=2.0, zorder=6)
    ax.set_aspect("equal")
    ax.set_axis_off()
    _save(fig, out_png)


def showcase_integration(mesh, out_png):
    """step07: raw streamline integration state (pre-snap, full spirals),
    folded into the fundamental pitch."""
    pitch = mesh.pitch_norm
    n_b = len(mesh.streamlines) - len(mesh.separatrices)
    boundary = [np.asarray(b, float) for b in mesh.streamlines[:n_b]]
    term = ps._termination_nodes(mesh)

    def _ends_ok(p):
        # judged mod pitch: ending on a wrap image of a target is fine
        for k in range(-ps.MAX_WRAPS, ps.MAX_WRAPS + 1):
            q = np.asarray(p, float) - [k * pitch, 0.0]
            if len(term) and np.min(np.linalg.norm(term - q, axis=1)) < 0.05:
                return True
            if ps._min_boundary_dist(q, boundary) < 0.05:
                return True
        return False

    fig, ax = plt.subplots(figsize=(4, 6))
    for i, s in enumerate(mesh.streamlines):
        s = np.asarray(s, float)
        if s.ndim != 2 or len(s) < 2:
            continue
        if i < n_b:
            if abs(s[-1, 0] - s[0, 0]) > 0.5 * pitch:  # ring spans >= 1 pitch
                for seg in fold_curve(s, pitch):
                    ax.plot(seg[:, 0], seg[:, 1], color="black", lw=1.5)
            else:
                ax.plot(s[:, 0], s[:, 1], color="black", lw=1.5)
        else:
            segs = fold_curve(s, pitch)
            # drop the trailing partial winding of a MAX_WRAPS-capped helix
            if segs and not _ends_ok(s[-1]):
                segs = segs[:-1]
            for seg in segs:
                ax.plot(seg[:, 0], seg[:, 1], color="#87CEEB", lw=1.3)
    ax.set_aspect("equal")
    ax.set_axis_off()
    _save(fig, out_png)


def showcase_simplification(mesh, merged_streamlines, wedge_arms, out_png):
    """step08: streamlines after the Xiao merge (Alg. 2 cases 1-3), with the
    wedge arms that the T-a seam postprocessing deletes shown dashed orange."""
    pitch = mesh.pitch_norm
    n_b = len(mesh.streamlines) - len(mesh.separatrices)
    boundary = [np.asarray(b, float) for b in mesh.streamlines[:n_b]]

    fig, ax = plt.subplots(figsize=(4, 6))
    for s in merged_streamlines:
        s = np.asarray(s, float)
        if s.ndim != 2 or len(s) < 2:
            continue
        for seg in fold_curve(s, pitch):
            ax.plot(seg[:, 0], seg[:, 1], color="#87CEEB", lw=1.3, zorder=2)
    for b in boundary:
        b = np.asarray(b, float)
        if b.ndim != 2 or len(b) < 2:
            continue
        if abs(b[-1, 0] - b[0, 0]) > 0.5 * pitch:
            for seg in fold_curve(b, pitch):
                ax.plot(seg[:, 0], seg[:, 1], color="black", lw=1.5, zorder=4)
        else:
            ax.plot(b[:, 0], b[:, 1], color="black", lw=1.5, zorder=4)
    for arm in wedge_arms:
        arm = np.asarray(arm, float)
        if arm.ndim != 2 or len(arm) < 2:
            continue
        for seg in fold_curve(arm, pitch):
            ax.plot(seg[:, 0], seg[:, 1], "--", color="darkorange",
                    lw=2.0, zorder=7)
    ax.set_aspect("equal")
    ax.set_axis_off()
    _save(fig, out_png)


def showcase_postprocessing(ta_run, out_png):
    """step09: final T-a block structure; wedge streamlines deleted at the
    seam drawn red, their source singularities ringed red."""
    result = ta_run["result"]
    mesh = ta_run["mesh"]
    boundary_ref = ta_run["boundary_ref"]
    nodes = result["nodes"]
    wedge_arms = getattr(mesh, "dropped_wedge_arms", []) or []

    def _bdist(p):
        return ps._min_boundary_dist(p, boundary_ref)

    _regular, _tnodes, irregular = tmf.node_regularity(result, _bdist)

    xy = mesh.x[:, 0:2].numpy()
    tris = mesh.faces.T.numpy()
    fig, ax = plt.subplots(figsize=(4, 6))
    ax.triplot(xy[:, 0], xy[:, 1], tris, lw=0.1, color="0.93")

    for blk in result["blocks"]:
        ring = blk["ring"]
        ax.fill(ring[:, 0], ring[:, 1], alpha=0.12, color="#87CEEB")
        for s in blk["sides"]:
            ax.plot(s[:, 0], s[:, 1], color="#3aa6d0", lw=1.1, zorder=3)
        c = nodes[blk["corners"]]
        ax.scatter(c[:, 0], c[:, 1], c="black", s=10, zorder=5)

    for b in boundary_ref:
        b = np.asarray(b, float)
        if b.ndim == 2 and len(b) >= 2:
            ax.plot(b[:, 0], b[:, 1], color="black", lw=1.6, zorder=4)

    sing_xy = nodes[irregular] if irregular else np.zeros((0, 2))
    ring_sings = set()
    for arm in wedge_arms:
        arm = np.asarray(arm, float)
        if arm.ndim != 2 or len(arm) < 2:
            continue
        ax.plot(arm[:, 0], arm[:, 1], color="red", lw=2.2, zorder=8)
        if len(sing_xy):
            d0 = np.min(np.linalg.norm(sing_xy - arm[0], axis=1))
            d1 = np.min(np.linalg.norm(sing_xy - arm[-1], axis=1))
            tip = arm[0] if d0 <= d1 else arm[-1]
            ring_sings.add(
                int(np.argmin(np.linalg.norm(sing_xy - tip, axis=1))))
    for si in ring_sings:
        ax.scatter(sing_xy[si, 0], sing_xy[si, 1], facecolors="none",
                   edgecolors="red", s=170, linewidths=2.2, zorder=9)

    ax.set_aspect("equal")
    ax.set_axis_off()
    _save(fig, out_png)


def run_showcase(stl, out_dir, ta_run):
    """Write the step01..step09 walkthrough plots for one surface.

    ``ta_run`` is a completed tmesh.run_tmesh(..., tag="ta") result; it
    supplies the final singularities (step06), the deleted wedge arms
    (step08) and the final block structure (step09). Assumes the module
    globals (periodic field, corner modes) are already configured exactly
    like the run that produced ``ta_run``.
    """
    from .dp_adapter import build_dp_data
    from .field import FrameField, StreamlineGenerator_v2
    from .field.singularity_detector import detect_singularities
    from .field.streamline_merging import StreamlineMerging
    from .unwrap_surface import unwrap

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    showcase_3d_surface(stl, out_dir / "step01_3d_surface.png")

    data = unwrap(stl)
    st, tris, loops = data["st"], data["tris"], data["loops"]
    showcase_unwrapped(st, tris, loops, out_dir / "step02_unwrapped.png")

    # fresh field state for the early steps (same config as the ta run)
    mesh, _tf = build_dp_data(stl)
    ff = FrameField(mesh)
    detect_singularities(ff.mesh)

    showcase_crossfield(mesh, st, tris, loops,
                        out_dir / "step03_crossfield.png")
    showcase_representatives(mesh, st, tris, loops,
                             out_dir / "step04_representatives.png",
                             out_dir / "step04_mapping.png")
    showcase_framefield(mesh, st, tris, loops,
                        out_dir / "step05_vectorfield.png",
                        out_dir / "step05_crossfield.png")

    # final partition singularities, mapped normalized [0,1]^2 -> raw (s,t)
    # (both frames share node ordering; per-axis linear fit recovers the map)
    sings = get_singularities(ta_run["mesh"])
    if sings:
        pos_norm = np.array([d["coords"] for d in sings.values()])
        mx = mesh.x[:, 0:2].numpy()
        cx = np.polyfit(mx[:, 0], st[:, 0], 1)
        cy = np.polyfit(mx[:, 1], st[:, 1], 1)
        final_sing_st = np.column_stack([np.polyval(cx, pos_norm[:, 0]),
                                         np.polyval(cy, pos_norm[:, 1])])
    else:
        final_sing_st = np.zeros((0, 2))
    showcase_singularities(mesh, st, tris, loops, final_sing_st,
                           out_dir / "step06_singularities.png")

    # raw integration state (pre-snap), then the Xiao-merged state
    sl = StreamlineGenerator_v2(ff.mesh)
    ps._drop_degenerate_corner_seps(sl.mesh)
    showcase_integration(sl.mesh, out_dir / "step07_integration.png")

    ps._snap_separatrix_endpoints(sl.mesh, radius=0.045)
    merging = StreamlineMerging(sl.mesh, verbose=False)
    wedge_arms = getattr(ta_run["mesh"], "dropped_wedge_arms", []) or []
    showcase_simplification(
        sl.mesh, [np.asarray(s, float) for s in merging.new_streamlines],
        wedge_arms,
        out_dir / "step08_streamline_simplification_merge.png")

    showcase_postprocessing(ta_run, out_dir / "step09_postprocessing.png")
