#!/usr/bin/env python3
"""Plot the integrated streamlines of ONE passage only: every separatrix is
folded back into the central pitch (s reduced modulo the pitch, relative to the
slanted seam), split where it crosses a seam. No neighbour copies, no tiling --
this shows the physical curves on the periodic surface.

Output: output/T1_9/hub_stage1/streamlines_central.png
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT.parent / "domain_partition_2D"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import unwrap_surface as us
import clean_separatrix as cs
import partition_surface as ps
from dp_adapter import build_dp_data
from tools import FrameField, StreamlineGenerator_v2
from tools.singularity_detector import detect_singularities

STL = str(REPO_ROOT / "T1_9_hub_raw.stl")
OUT = REPO_ROOT / "output" / "T1_9" / "hub_stage1"

# left seam runs (0,0)->(0.407,1): s_seam(t) = SEAM_SLOPE * t
SEAM_SLOPE = 0.407


def fold_curve(s, pitch):
    """Fold a cover-coordinate polyline into the fundamental pitch, splitting
    it into segments wherever it crosses a seam. The seam crossing point is
    interpolated exactly and appended to BOTH adjacent segments, so every
    winding is drawn seam-to-seam with no visual gap."""
    s = np.asarray(s, float)
    u = s[:, 0] - SEAM_SLOPE * s[:, 1]        # seam-aligned coordinate
    k = np.floor(u / pitch).astype(int)
    segs, cur = [], [s[0]]
    for i in range(1, len(s)):
        if k[i] != k[i - 1]:
            ub = pitch * max(k[i], k[i - 1])  # crossed seam level in u
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


def plot_state(mesh, title, out_png):
    """Fold-and-plot the current mesh.streamlines/separatrices state."""
    pitch = mesh.pitch_norm
    xy = mesh.x[:, 0:2].numpy()
    tris = mesh.faces.T.numpy()
    n_b = len(mesh.streamlines) - len(mesh.separatrices)

    # targets a curve may legitimately end on: boundaries + termination points
    boundary = [np.asarray(b, float) for b in mesh.streamlines[:n_b]]
    term = ps._termination_nodes(mesh)

    def _ends_ok(p):
        # judged modulo the pitch: a curve ending on a wrap image of a
        # termination node / boundary ended fine.
        for k in range(-ps.MAX_WRAPS, ps.MAX_WRAPS + 1):
            q = np.asarray(p, float) - [k * pitch, 0.0]
            if len(term) and np.min(np.linalg.norm(term - q, axis=1)) < 0.05:
                return True
            if ps._min_boundary_dist(q, boundary) < 0.05:
                return True
        return False

    fig, ax = plt.subplots(figsize=(11, 10))
    ax.triplot(xy[:, 0], xy[:, 1], tris, lw=0.1, color="0.92")
    for i, s in enumerate(mesh.streamlines):
        s = np.asarray(s, float)
        if s.ndim != 2 or len(s) < 2:
            continue
        if i < n_b:
            # boundary section: walls/blade stay as-is; closed orbit RINGS
            # live here too and are folded like curves
            if abs(s[-1, 0] - s[0, 0]) > 0.5 * pitch:   # ring (spans >=1 pitch)
                for seg in fold_curve(s, pitch):
                    ax.plot(seg[:, 0], seg[:, 1], "C0", lw=1.8)
            else:
                ax.plot(s[:, 0], s[:, 1], "0.45", lw=1.0)
        else:
            segs = fold_curve(s, pitch)
            # drop the trailing partial winding of a MAX_WRAPS-capped helix:
            # it ends mid-air (integration cap artifact, not physics). All
            # other segments end on a seam (interpolated) or a real target.
            if segs and not _ends_ok(s[-1]):
                segs = segs[:-1]
            for seg in segs:
                ax.plot(seg[:, 0], seg[:, 1], "C3", lw=1.3)
    sc = xy[tris[mesh.singularities.numpy() != 0]].mean(axis=1)
    ax.scatter(sc[:, 0], sc[:, 1], c="blue", s=70, zorder=6, label="singularity")
    cm = mesh.x[:, 2].numpy() == 0
    ax.scatter(xy[cm, 0], xy[cm, 1], c="green", marker="^", s=60, zorder=6,
               label="c0 corner")
    ax.set_aspect("equal")
    ax.legend()
    ax.set_title(title)
    fig.savefig(out_png, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_png}")


def main():
    us.set_blade_tip_corners(False)
    cs.set_emanate_outer_corners(True)
    mesh, tf = build_dp_data(STL)
    ff = FrameField(mesh)
    m = detect_singularities(ff.mesh)
    sl = StreamlineGenerator_v2(ff.mesh)
    ps._drop_degenerate_corner_seps(sl.mesh)
    pitch = sl.mesh.pitch_norm

    # integration is expensive -- run once, snapshot, then produce both states
    raw_sl = [np.array(s, float) for s in sl.mesh.streamlines]
    raw_sep = [dict(d) for d in sl.mesh.separatrices]

    # state A: NO merge (showcase: full spirals), snap only (dedups doubles)
    ps._snap_separatrix_endpoints(sl.mesh, radius=0.045)
    plot_state(sl.mesh,
               f"All streamlines, NO merge: {len(sl.mesh.separatrices)} "
               f"separatrices folded mod pitch ({pitch:.3f})",
               OUT / "streamlines_central_raw.png")

    # state B: convergence-based merge (helix -> prong + orbit ring), then snap
    sl.mesh.streamlines = raw_sl
    sl.mesh.separatrices = raw_sep
    ps._close_helical_streamlines(sl.mesh)
    ps._emit_dock_crossings(sl)
    ps._snap_separatrix_endpoints(sl.mesh, radius=0.045)
    plot_state(sl.mesh,
               f"Convergence-merged: {len(sl.mesh.separatrices)} separatrices "
               f"+ orbit rings, folded mod pitch ({pitch:.3f})",
               OUT / "streamlines_central_merged.png")


if __name__ == "__main__":
    main()
