"""Showcase plots: one consistent presentation-style series (transparent
background, no axes/titles) walking through every pipeline stage.

Steps 1-7 are method-independent (all methods share field and integration up
to the endpoint snap) and are written once per part; steps 8-11 are written
per method from the completed run dicts. Entry point: write_plots().
"""

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from collections import defaultdict
from pathlib import Path

from . import partition_surface as ps
from . import tmesh_faces as tmf

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


def singularity_connections(mesh, singularities):
    """{(Sa, Sb): [streamline indices]} of direct singularity-singularity
    separatrix connections."""
    n_b = len(mesh.streamlines) - len(mesh.separatrices)
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
    return connections


def _save(fig, out_png):
    fig.savefig(out_png, dpi=200, transparent=True, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_png}")


def _draw_domain(ax, st, tris, loops):
    """Shared background: faint triangulation + black boundary contours."""
    ax.triplot(st[:, 0], st[:, 1], tris, lw=0.1, color="0.93")
    for loop in loops:
        ring = st[loop + [loop[0]]]
        ax.plot(ring[:, 0], ring[:, 1], lw=2.0, color="black")


# --------------------------------------------------------------------------
# steps 1-7: method-independent
# --------------------------------------------------------------------------

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


def showcase_representatives(mesh, st, tris, loops, out_png):
    """step04: boundary cross-field with the representative vector per node."""
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
    _save(fig, out_png)


def showcase_mapping(out_png):
    """step04 companion: conceptual cross -> 4*theta representative mapping."""
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
    _save(fig, out_png)


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


def _ends_ok_fn(mesh):
    """Predicate: does a curve endpoint land on a legitimate target (boundary
    or termination node), judged modulo the pitch?"""
    pitch = mesh.pitch_norm
    n_b = len(mesh.streamlines) - len(mesh.separatrices)
    boundary = [np.asarray(b, float) for b in mesh.streamlines[:n_b]]
    term = ps._termination_nodes(mesh)

    def _ends_ok(p):
        for k in range(-ps.MAX_WRAPS, ps.MAX_WRAPS + 1):
            q = np.asarray(p, float) - [k * pitch, 0.0]
            if len(term) and np.min(np.linalg.norm(term - q, axis=1)) < 0.05:
                return True
            if ps._min_boundary_dist(q, boundary) < 0.05:
                return True
        return False

    return _ends_ok


def _fold_boundary(ax, s, pitch, **kw):
    """Boundary curve: rings spanning >= 1 pitch are folded, walls drawn
    as-is."""
    if abs(s[-1, 0] - s[0, 0]) > 0.5 * pitch:
        for seg in fold_curve(s, pitch):
            ax.plot(seg[:, 0], seg[:, 1], **kw)
    else:
        ax.plot(s[:, 0], s[:, 1], **kw)


def showcase_integration(mesh, out_png):
    """step07: raw streamline integration state (pre-snap, full spirals),
    folded into the fundamental pitch."""
    pitch = mesh.pitch_norm
    n_b = len(mesh.streamlines) - len(mesh.separatrices)
    _ends_ok = _ends_ok_fn(mesh)

    fig, ax = plt.subplots(figsize=(4, 6))
    for i, s in enumerate(mesh.streamlines):
        s = np.asarray(s, float)
        if s.ndim != 2 or len(s) < 2:
            continue
        if i < n_b:
            _fold_boundary(ax, s, pitch, color="black", lw=1.5)
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


def showcase_integration_labeled(mesh, out_png):
    """step07 labeled: every curve in its own color with an index label on
    its longest folded segment; singularities as color-coded stars with
    'S1 (x, y)' legend entries instead of in-plot text boxes."""
    pitch = mesh.pitch_norm
    n_b = len(mesh.streamlines) - len(mesh.separatrices)
    _ends_ok = _ends_ok_fn(mesh)
    xy = mesh.x[:, 0:2].numpy()
    tris = mesh.faces.T.numpy()

    fig, ax = plt.subplots(figsize=(6, 9))
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
        segs = fold_curve(s, pitch)
        if i >= n_b and segs and not _ends_ok(s[-1]):
            segs = segs[:-1]
        for seg in segs:
            ax.plot(seg[:, 0], seg[:, 1], color=color, lw=lw, alpha=alpha,
                    zorder=zorder)
        if segs:
            longest = max(segs, key=len)
            mid = longest[len(longest) // 2]
            ax.text(mid[0], mid[1], str(i), fontsize=7,
                    color="black" if i >= n_b else "0.4",
                    ha="center", va="center",
                    bbox=dict(boxstyle="round,pad=0.15", facecolor="white",
                              edgecolor="none", alpha=0.7),
                    zorder=zorder + 1)

    sing_cmap = plt.cm.tab10
    for si, (name, data) in enumerate(get_singularities(mesh).items()):
        c = data["coords"]
        ax.scatter(c[0], c[1], c=[sing_cmap(si % 10)], s=180, marker="*",
                   zorder=8, edgecolors="black", linewidths=0.8,
                   label=f"{name} ({c[0]:.3f}, {c[1]:.3f})")

    ax.plot([], [], "-", color="0.5", lw=0.7, alpha=0.5, label="boundary")
    ax.plot([], [], "-", color="0.3", lw=2.0, label="separatrix")
    ax.set_aspect("equal")
    ax.set_axis_off()
    ax.legend(loc="upper left", fontsize=8)
    _save(fig, out_png)


# --------------------------------------------------------------------------
# steps 8-11: per method
# --------------------------------------------------------------------------

def showcase_simplification(mesh, merged_streamlines, wedge_arms, out_png):
    """step08: streamlines after the Xiao merge (Alg. 2 cases 1-3). For ta/tb
    the wedge arms that the seam postprocessing deletes are shown dashed
    orange; pass an empty ``wedge_arms`` for xiao."""
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
        _fold_boundary(ax, b, pitch, color="black", lw=1.5, zorder=4)
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


def showcase_postprocessing(run, out_png):
    """step09: final block structure of one method; wedge streamlines deleted
    at the seam (ta/tb) drawn red with their source singularities ringed."""
    result = run["result"]
    mesh = run["mesh"]
    boundary_ref = run["boundary_ref"]
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


def showcase_tfi(tfi, result, out_png):
    """step10: TFI grid (ta/tb only)."""
    fig, ax = plt.subplots(figsize=(4, 6))
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
    _save(fig, out_png)


def showcase_tiled(result, tfi, pitch, out_png):
    """step11: TFI grid tiled over 3 passages (shifts -pitch, 0, +pitch), the
    visual periodicity check -- grid points of neighbouring copies must
    coincide at the interior seams."""
    shifts = [-pitch, 0.0, pitch]
    cols = ["0.55", "black", "0.55"]

    fig, ax = plt.subplots(figsize=(18, 7))
    for sh, col in zip(shifts, cols):
        for X in tfi["grids"]:
            for i in range(X.shape[0]):
                ax.plot(X[i, :, 0] + sh, X[i, :, 1], col, lw=0.25)
            for j in range(X.shape[1]):
                ax.plot(X[:, j, 0] + sh, X[:, j, 1], col, lw=0.25)
        for blk in result["blocks"]:
            for s in blk["sides"]:
                ax.plot(s[:, 0] + sh, s[:, 1], col, lw=0.8)
    ax.set_aspect("equal")
    ax.set_axis_off()
    _save(fig, out_png)


# --------------------------------------------------------------------------
# orchestrator
# --------------------------------------------------------------------------

def write_plots(stl, part, runs, out_dir):
    """Write the full showcase series for one surface into ``out_dir``.

    ``runs`` maps method name (ta/tb/xiao) to its completed run dict from
    tmesh.run_tmesh / xiao.run_xiao. Steps 1-7 are written once per part,
    steps 8-11 per method. Assumes the module globals (periodic field,
    corner modes) are still configured like the runs.
    """
    from .dp_adapter import build_dp_data
    from .field import FrameField, StreamlineGenerator_v2
    from .field.singularity_detector import detect_singularities
    from .field.streamline_merging import StreamlineMerging
    from .unwrap_surface import unwrap

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    showcase_3d_surface(stl, out_dir / f"step01_3d_surface_{part}.png")

    data = unwrap(stl)
    st, tris, loops = data["st"], data["tris"], data["loops"]
    showcase_unwrapped(st, tris, loops,
                       out_dir / f"step02_unwrapped_{part}.png")

    # fresh field state, same config as the runs
    mesh, _tf = build_dp_data(stl)
    ff = FrameField(mesh)
    detect_singularities(ff.mesh)

    showcase_crossfield(mesh, st, tris, loops,
                        out_dir / f"step03_crossfield_{part}.png")
    showcase_representatives(mesh, st, tris, loops,
                             out_dir / f"step04_representatives_{part}.png")
    showcase_mapping(out_dir / "step04_mapping.png")
    showcase_framefield(mesh, st, tris, loops,
                        out_dir / f"step05_vectorfield_{part}.png",
                        out_dir / f"step05_framefield_{part}.png")

    # final partition singularities (identical across methods), mapped
    # normalized [0,1]^2 -> raw (s,t): both frames share node ordering, so a
    # per-axis linear fit recovers the map
    ref = next((runs[m] for m in ("ta", "tb", "xiao") if m in runs), None)
    final_sing_st = np.zeros((0, 2))
    if ref is not None:
        sings = get_singularities(ref["mesh"])
        if sings:
            pos_norm = np.array([d["coords"] for d in sings.values()])
            mx = mesh.x[:, 0:2].numpy()
            cx = np.polyfit(mx[:, 0], st[:, 0], 1)
            cy = np.polyfit(mx[:, 1], st[:, 1], 1)
            final_sing_st = np.column_stack(
                [np.polyval(cx, pos_norm[:, 0]),
                 np.polyval(cy, pos_norm[:, 1])])
    showcase_singularities(mesh, st, tris, loops, final_sing_st,
                           out_dir / f"step06_singularities_{part}.png")

    # raw integration state (pre-snap), then the Xiao-merged state
    sl = StreamlineGenerator_v2(ff.mesh)
    ps._drop_degenerate_corner_seps(sl.mesh)
    showcase_integration(
        sl.mesh, out_dir / f"step07_streamline_integration_{part}.png")
    showcase_integration_labeled(
        sl.mesh,
        out_dir / f"step07_streamline_integration_{part}_labeled.png")

    ps._snap_separatrix_endpoints(sl.mesh, radius=0.045)
    merging = StreamlineMerging(sl.mesh, verbose=False)
    merged = [np.asarray(s, float) for s in merging.new_streamlines]

    for method, run in runs.items():
        wedge_arms = getattr(run["mesh"], "dropped_wedge_arms", []) or []
        showcase_simplification(
            sl.mesh, merged, wedge_arms,
            out_dir / f"step08_simplification_{part}_{method}.png")
        showcase_postprocessing(
            run, out_dir / f"step09_postprocessing_{part}_{method}.png")
        tfi = run.get("tfi")
        if tfi and tfi["grids"]:
            showcase_tfi(tfi, run["result"],
                         out_dir / f"step10_tfi_{part}_{method}.png")
            showcase_tiled(run["result"], tfi, float(run["mesh"].pitch_norm),
                           out_dir / f"step11_tiled_{part}_{method}.png")
