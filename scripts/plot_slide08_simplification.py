#!/usr/bin/env python3
"""
Stage 1 / Step 08: Streamlines AFTER the Xiao streamline-simplification
(the 3 merging cases: Case 1 direct connect, Case 2 missed endpoints,
Case 3 pass-through — Xiao 2020, Alg. 2), still containing the streamlines
that the Ansatz T-a seam post-processing deletes later.

Single-panel figure (saved as a single PNG), style of slides 1-6:
  - Boundary streamlines (walls/blade): black, lw=1.5
  - Merged separatrices: light blue (#87CEEB), lw=1.3
  - Wedge arms that T-a deletes (``mesh.dropped_wedge_arms``): dashed
    dark orange, lw=2.0 — these survive Xiao but are removed by T-a.
  - No mesh background, no axes, no labels, no titles, no legends.
  - Transparent background, portrait (4, 6), dpi >= 200.

The merged state is the project's Xiao-equivalent streamline simplification
(``_close_helical_streamlines`` + ``_emit_dock_crossings`` +
``_snap_separatrix_endpoints``) — the SAME state the T-a pipeline consumes
before its seam post-processing. The deleted wedge arms are obtained by
running the T-a pipeline (``tmesh_partition.run_tmesh``) and reading
``mesh.dropped_wedge_arms``.

Output: plots/step08_streamline_simplification_merge.png
"""

import sys
import tempfile
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
STAGE1_DIR = SCRIPT_DIR / "stage1"
REPO_ROOT = SCRIPT_DIR.parent
DP2D_DIR = REPO_ROOT.parent / "domain_partition_2D"

if str(STAGE1_DIR) not in sys.path:
    sys.path.insert(0, str(STAGE1_DIR))
if str(DP2D_DIR) not in sys.path:
    sys.path.insert(0, str(DP2D_DIR))

STL_PATH = REPO_ROOT / "T1_9_hub_raw.stl"
OUT_PATH = REPO_ROOT / "plots" / "step08_streamline_simplification_merge.png"

import matplotlib

matplotlib.use("Agg")  # headless rendering
import matplotlib.pyplot as plt
import numpy as np

import unwrap_surface as us  # noqa: E402
import clean_separatrix as cs  # noqa: E402
import partition_surface as ps  # noqa: E402
import tmesh_partition as tp  # noqa: E402
from dp_adapter import build_dp_data  # noqa: E402
from tools import FrameField, StreamlineGenerator_v2  # noqa: E402
from tools.singularity_detector import detect_singularities  # noqa: E402
from tools.streamline_merging import StreamlineMerging  # noqa: E402

# Reuse the exact fold-into-central-pitch helper used by Step 07.
from plot_slide07_integration import fold_curve, SEAM_SLOPE  # noqa: E402,F401


def main() -> None:
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)

    us.set_blade_tip_corners(False)
    cs.set_emanate_outer_corners(True)
    # Same periodic cross-field configuration as the T-a pipeline.
    ps.set_periodic(True)
    ps.set_tile_periodic(False)

    # --- 1. Wedge arms deleted by T-a (run the T-a pipeline, no plots) ---
    # tmesh_partition.run_tmesh defaults to the hub STL; it sets
    # mesh.dropped_wedge_arms during collapse_seam_wedges.
    with tempfile.TemporaryDirectory() as tmp:
        ta = tp.run_tmesh(out_dir=Path(tmp), continue_seam_edges=False,
                          tag="ta", make_plots=False, verbose=False)
    wedge_arms = getattr(ta["mesh"], "dropped_wedge_arms", []) or []

    # --- 2. Merged streamline state after the 3 Xiao merging cases ---
    # This is the EXACT merge the partition pipeline uses (StreamlineMerging,
    # Xiao 2020 Alg. 2 Case 1/2/3), same as tmesh_partition.run_tmesh:1071 and
    # compare_xiao_vs_ta.run_xiao — NOT the plot_central helix-close variant.
    # We run it WITHOUT the T-a seam post-processing (collapse_seam_wedges), so
    # the wedge streamlines that T-a deletes are still present here.
    mesh, _tf = build_dp_data(STL_PATH)
    ff = FrameField(mesh)
    detect_singularities(ff.mesh)
    sl = StreamlineGenerator_v2(ff.mesh)
    ps._drop_degenerate_corner_seps(sl.mesh)
    ps._snap_separatrix_endpoints(sl.mesh, radius=0.045)
    pitch = sl.mesh.pitch_norm

    n_b = len(sl.mesh.streamlines) - len(sl.mesh.separatrices)
    boundary = [np.asarray(b, float) for b in sl.mesh.streamlines[:n_b]]

    merging = StreamlineMerging(sl.mesh, verbose=False)
    merged = [np.asarray(s, float) for s in merging.new_streamlines]

    # --- 3. Plot ---
    fig, ax = plt.subplots(figsize=(4, 6))  # portrait

    # Merged streamlines (boundary + separatrices) in light blue; the domain
    # boundary is overdrawn black on top so walls/blade read as boundary.
    for s in merged:
        if s.ndim != 2 or len(s) < 2:
            continue
        for seg in fold_curve(s, pitch):
            ax.plot(seg[:, 0], seg[:, 1], color="#87CEEB", lw=1.3, zorder=2)

    for b in boundary:
        b = np.asarray(b, float)
        if b.ndim != 2 or len(b) < 2:
            continue
        if abs(b[-1, 0] - b[0, 0]) > 0.5 * pitch:   # ring spanning >=1 pitch
            for seg in fold_curve(b, pitch):
                ax.plot(seg[:, 0], seg[:, 1], color="black", lw=1.5, zorder=4)
        else:
            ax.plot(b[:, 0], b[:, 1], color="black", lw=1.5, zorder=4)

    # Wedge arms deleted by T-a: dashed dark orange (folded like the rest).
    for arm in wedge_arms:
        arm = np.asarray(arm, float)
        if arm.ndim != 2 or len(arm) < 2:
            continue
        for seg in fold_curve(arm, pitch):
            ax.plot(seg[:, 0], seg[:, 1], "--", color="darkorange",
                    lw=2.0, zorder=7)

    ax.set_aspect("equal")
    ax.set_axis_off()

    fig.savefig(OUT_PATH, dpi=200, transparent=True, bbox_inches="tight")
    plt.close(fig)

    print(f"wrote {OUT_PATH}")
    print(f"  pitch:            {pitch:.3f}")
    print(f"  boundary curves:  {n_b}")
    print(f"  merged streamlines (Xiao 3-case): {len(merged)}")
    print(f"  T-a deleted arms: {len(wedge_arms)}")


if __name__ == "__main__":
    main()
