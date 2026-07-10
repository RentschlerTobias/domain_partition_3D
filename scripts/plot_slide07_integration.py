#!/usr/bin/env python3
"""
Stage 1 / Slide 07: All streamlines emanating from singularities.

Single-panel figure (saved as a single PNG):
  - Boundary streamlines (walls/blade): black, lw=1.5
  - Separatrices (from singularities): light blue (#87CEEB), lw=1.3
  - No mesh background, no axes, no labels, no titles, no legends.
  - Transparent background (alpha channel).
  - Portrait orientation (figsize=(4, 6)).
  - DPI >= 200.

Reuses the fold-and-plot logic from scripts/stage1/plot_central_streamlines.py
without modifying it. The streamline integration is performed by
``StreamlineGenerator_v2`` (which sets ``mesh.streamlines`` and
``mesh.separatrices``). The fold-into-central-pitch logic is the same
``fold_curve`` helper used by ``plot_central_streamlines.plot_state``.

Output: plots/slide07_integration.png
"""

import sys
from pathlib import Path

# Make scripts/stage1/ importable so we can reuse plot_central_streamlines.py
# logic (StreamlineGenerator_v2, build_dp_data, FrameField, partition_surface,
# unwrap_surface, clean_separatrix, singularity_detector) as-is.
SCRIPT_DIR = Path(__file__).resolve().parent
STAGE1_DIR = SCRIPT_DIR / "stage1"
REPO_ROOT = SCRIPT_DIR.parent
DP2D_DIR = REPO_ROOT.parent / "domain_partition_2D"

if str(STAGE1_DIR) not in sys.path:
    sys.path.insert(0, str(STAGE1_DIR))
if str(DP2D_DIR) not in sys.path:
    sys.path.insert(0, str(DP2D_DIR))

STL_PATH = REPO_ROOT / "T1_9_hub_raw.stl"
OUT_PATH = REPO_ROOT / "plots" / "step07_integration.png"

import matplotlib

matplotlib.use("Agg")  # headless rendering
import matplotlib.pyplot as plt
import numpy as np

import unwrap_surface as us  # noqa: E402
import clean_separatrix as cs  # noqa: E402
import partition_surface as ps  # noqa: E402
from dp_adapter import build_dp_data  # noqa: E402
from tools import FrameField, StreamlineGenerator_v2  # noqa: E402
from tools.singularity_detector import detect_singularities  # noqa: E402

# Left seam runs (0,0)->(0.407,1): s_seam(t) = SEAM_SLOPE * t
SEAM_SLOPE = 0.407


def fold_curve(s, pitch):
    """Fold a cover-coordinate polyline into the fundamental pitch, splitting
    it into segments wherever it crosses a seam. The seam crossing point is
    interpolated exactly and appended to BOTH adjacent segments, so every
    winding is drawn seam-to-seam with no visual gap.

    Identical to plot_central_streamlines.fold_curve (reused verbatim).
    """
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


def plot_integration(mesh, out_png):
    """Fold-and-plot the mesh.streamlines/separatrices state with the
    slide07 color scheme: boundary streamlines black, separatrices light blue.
    No mesh background, no axes, no title, no legend.
    """
    pitch = mesh.pitch_norm
    n_b = len(mesh.streamlines) - len(mesh.separatrices)

    # Targets a curve may legitimately end on: boundaries + termination points.
    boundary = [np.asarray(b, float) for b in mesh.streamlines[:n_b]]
    term = ps._termination_nodes(mesh)

    def _ends_ok(p):
        # Judged modulo the pitch: a curve ending on a wrap image of a
        # termination node / boundary ended fine.
        for k in range(-ps.MAX_WRAPS, ps.MAX_WRAPS + 1):
            q = np.asarray(p, float) - [k * pitch, 0.0]
            if len(term) and np.min(np.linalg.norm(term - q, axis=1)) < 0.05:
                return True
            if ps._min_boundary_dist(q, boundary) < 0.05:
                return True
        return False

    fig, ax = plt.subplots(figsize=(4, 6))  # portrait

    # NO mesh background (clean) — per slide07 spec.

    for i, s in enumerate(mesh.streamlines):
        s = np.asarray(s, float)
        if s.ndim != 2 or len(s) < 2:
            continue
        if i < n_b:
            # Boundary section: walls/blade stay as-is; closed orbit RINGS
            # live here too and are folded like curves.
            if abs(s[-1, 0] - s[0, 0]) > 0.5 * pitch:   # ring (spans >=1 pitch)
                for seg in fold_curve(s, pitch):
                    ax.plot(seg[:, 0], seg[:, 1], color="black", lw=1.5)
            else:
                ax.plot(s[:, 0], s[:, 1], color="black", lw=1.5)
        else:
            segs = fold_curve(s, pitch)
            # Drop the trailing partial winding of a MAX_WRAPS-capped helix:
            # it ends mid-air (integration cap artifact, not physics). All
            # other segments end on a seam (interpolated) or a real target.
            if segs and not _ends_ok(s[-1]):
                segs = segs[:-1]
            for seg in segs:
                ax.plot(seg[:, 0], seg[:, 1], color="#87CEEB", lw=1.3)

    ax.set_aspect("equal")
    ax.set_axis_off()  # no axes, no labels, no titles

    fig.savefig(
        out_png,
        dpi=200,
        transparent=True,
        bbox_inches="tight",
    )
    plt.close(fig)


def main() -> None:
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)

    # --- 1. Build the mesh + cross-field (same as plot_central_streamlines.main) ---
    us.set_blade_tip_corners(False)
    cs.set_emanate_outer_corners(True)
    mesh, _transform = build_dp_data(STL_PATH)
    ff = FrameField(mesh)
    m = detect_singularities(ff.mesh)
    sl = StreamlineGenerator_v2(ff.mesh)
    ps._drop_degenerate_corner_seps(sl.mesh)
    pitch = sl.mesh.pitch_norm

    # --- 2. RAW integration state: NO merge, NO simplification/snap ---
    # Step 7 shows the streamlines directly as integrated from the
    # singularities (full spirals). The merge/simplification steps
    # (_close_helical_streamlines, _emit_dock_crossings,
    # _snap_separatrix_endpoints) come LATER (Step 8) and are omitted here.

    n_b = len(sl.mesh.streamlines) - len(sl.mesh.separatrices)
    print(f"  pitch:           {pitch:.3f}")
    print(f"  boundary curves: {n_b}")
    print(f"  separatrices:    {len(sl.mesh.separatrices)}")

    # --- 3. Plot ---
    plot_integration(sl.mesh, OUT_PATH)
    print(f"wrote {OUT_PATH}")


if __name__ == "__main__":
    main()
