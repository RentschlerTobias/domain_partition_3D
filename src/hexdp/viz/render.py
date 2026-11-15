"""Static PNG views of geometries and block structures (matplotlib, offscreen).

Small and dependency-light on purpose: the MCP server returns these images to an agent, and the
tutorials / UI reuse the same drawing functions (interactive plotly views come with the UI step).

    render_geometry(geo)              boundary triangles coloured by patch role, feature curves in black
    render_blocks(bs)                 the 12 edges of every block, coloured by block, fine-grid edges if any
    render_quality(bs)                histogram of the fine-cell scaled Jacobian
"""

from __future__ import annotations

import io

import numpy as np

ROLE_COLORS = {"wall": "#8c8c8c", "inlet": "#2a9d8f", "outlet": "#e76f51", "periodic": "#457b9d",
               "interface": "#e9c46a", "symmetry": "#a8dadc", "unknown": "#c9b8e8"}


def _fig(size=(7, 6)):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=size, dpi=110)
    return fig, plt


def _png(fig, plt) -> bytes:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight")
    plt.close(fig)
    return buf.getvalue()


def _equal(ax, P):
    lo, hi = P.min(0), P.max(0)
    c, r = 0.5 * (lo + hi), 0.5 * float((hi - lo).max())
    ax.set_xlim(c[0] - r, c[0] + r); ax.set_ylim(c[1] - r, c[1] + r); ax.set_zlim(c[2] - r, c[2] + r)
    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("z")


def render_geometry(geo, elev=25, azim=-60, max_tris=20000, title=None) -> bytes:
    from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection

    fig, plt = _fig()
    ax = fig.add_subplot(projection="3d")
    allP = []
    rng = np.random.default_rng(0)
    for s in geo.surfaces.values():
        T = s.tris
        if len(T) > max_tris // max(1, len(geo.surfaces)):
            T = T[rng.choice(len(T), max_tris // max(1, len(geo.surfaces)), replace=False)]
        tri = s.nodes[T]
        allP.append(s.nodes)
        ax.add_collection3d(Poly3DCollection(tri, facecolor=ROLE_COLORS.get(s.role.value, "#cccccc"),
                                             edgecolor="none", alpha=0.35))
    lines = [np.c_[c.points[:-1], c.points[1:]].reshape(-1, 2, 3) for c in geo.feature_curves() if len(c.points) > 1]
    if lines:
        ax.add_collection3d(Line3DCollection(np.vstack(lines), colors="k", linewidths=0.6))
    _equal(ax, np.vstack(allP))
    ax.view_init(elev=elev, azim=azim)
    roles = sorted({s.role.value for s in geo.surfaces.values()})
    for r in roles:
        ax.plot([], [], "s", color=ROLE_COLORS.get(r, "#cccccc"), label=r)
    ax.legend(loc="upper left", fontsize=8)
    ax.set_title(title or f"{len(geo.surfaces)} surfaces, {len(geo.feature_curves())} feature curves", fontsize=9)
    return _png(fig, plt)


def _block_edge_lines(bs):
    """Per block its 12 edges as polylines: from the fine grid when present, else straight."""
    from hexdp.core.blockstructure import HEX_EDGES

    out = []
    for b, c8 in enumerate(bs.blocks):
        if bs.grids:
            G = bs.grids[b]
            X = bs.fine_points[G]
            for ax_ in range(3):
                o = [d for d in range(3) if d != ax_]
                for e1 in (0, -1):
                    for e2 in (0, -1):
                        idx = [None] * 3
                        idx[ax_] = slice(None); idx[o[0]] = e1; idx[o[1]] = e2
                        out.append((b, X[tuple(idx)]))
        else:
            for i, j in HEX_EDGES:
                out.append((b, bs.points[[c8[i], c8[j]]]))
    return out


def render_blocks(bs, elev=25, azim=-60, title=None, highlight=()) -> bytes:
    from mpl_toolkits.mplot3d.art3d import Line3DCollection

    fig, plt = _fig()
    ax = fig.add_subplot(projection="3d")
    cmap = plt.get_cmap("tab20")
    segs, cols, widths = [], [], []
    for b, L in _block_edge_lines(bs):
        s = np.stack([L[:-1], L[1:]], 1)
        segs.append(s)
        hl = b in highlight
        cols += [(0.85, 0.1, 0.1, 1.0) if hl else cmap(b % 20)] * len(s)
        widths += [2.0 if hl else 0.8] * len(s)
    ax.add_collection3d(Line3DCollection(np.vstack(segs), colors=cols, linewidths=widths))
    _equal(ax, bs.fine_points if bs.fine_points is not None else bs.points)
    ax.view_init(elev=elev, azim=azim)
    s = bs.summary()
    ax.set_title(title or f"{s['blocks']} blocks, min sJ {s['min_sj_fine'] if s['min_sj_fine'] is not None else s['min_sj_coarse']:.3f}, "
                          f"issues {s['issues'] or 'none'}", fontsize=9)
    return _png(fig, plt)


def render_quality(bs, bins=50) -> bytes:
    from hexdp.core.blockstructure import grid_cells, hex_corner_jacobians

    fig, plt = _fig((6, 3.5))
    ax = fig.add_subplot()
    if bs.grids:
        sj = np.concatenate([hex_corner_jacobians(bs.fine_points[grid_cells(G)]).min(1) for G in bs.grids])
    else:
        sj = bs.corner_jacobians().min(1)
    ax.hist(sj, bins=bins, color="#457b9d")
    ax.axvline(0, color="r", lw=1)
    ax.set_xlabel("scaled Jacobian (min over cell corners)"); ax.set_ylabel("cells")
    ax.set_title(f"{len(sj)} cells, min {sj.min():.3f}, p1 {np.percentile(sj, 1):.3f}, {int((sj <= 0).sum())} inverted",
                 fontsize=9)
    return _png(fig, plt)
