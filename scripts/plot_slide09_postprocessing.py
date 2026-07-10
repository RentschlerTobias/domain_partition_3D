#!/usr/bin/env python3
"""
Stage 1 / Step 09: Final T-a block structure (postprocessing), styled to match
slides 1-8 (pure image, transparent, portrait, no axes/title/legend).

Two marker families (both, per user request):
  - RED open circle rings: the irregular interior nodes (valence != 4) = the
    two blade field singularities (valence 5), the legitimate irregular block
    corners that remain.
  - ORANGE dashed arcs + rings: the wedge streamlines that T-a DELETES at the
    periodic seam (``mesh.dropped_wedge_arms``) — clarifies WHERE streamlines
    are removed (the same arms shown dashed in Step 8).

Block data + regularity come straight from the real partition pipeline
(``tmesh_partition.run_tmesh`` -> ``tmesh_faces.node_regularity``), so this is
the actual T-a result, not a re-drawing.

Output: plots/step09_postprocessing.png
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

OUT_PATH = REPO_ROOT / "plots" / "step09_postprocessing.png"

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import unwrap_surface as us  # noqa: E402
import clean_separatrix as cs  # noqa: E402
import partition_surface as ps  # noqa: E402
import tmesh_faces as tmf  # noqa: E402
import tmesh_partition as tp  # noqa: E402


def main() -> None:
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)

    us.set_blade_tip_corners(False)
    cs.set_emanate_outer_corners(True)

    # --- 1. Run the real T-a pipeline (no plots) ---
    with tempfile.TemporaryDirectory() as tmp:
        ta = tp.run_tmesh(out_dir=Path(tmp), continue_seam_edges=False,
                          tag="ta", make_plots=False, verbose=False)
    result = ta["result"]
    mesh = ta["mesh"]
    boundary_ref = ta["boundary_ref"]
    nodes = result["nodes"]
    wedge_arms = getattr(mesh, "dropped_wedge_arms", []) or []

    # Irregular interior nodes (valence != 4) = blade field singularities.
    def _bdist(p):
        return ps._min_boundary_dist(p, boundary_ref)

    _regular, _tnodes, irregular = tmf.node_regularity(result, _bdist)

    # --- 2. Plot (minimalist, slides 1-8 style) ---
    xy = mesh.x[:, 0:2].numpy()
    tris = mesh.faces.T.numpy()
    fig, ax = plt.subplots(figsize=(4, 6))  # portrait

    ax.triplot(xy[:, 0], xy[:, 1], tris, lw=0.1, color="0.93")

    # Block cells: subtle light-blue fill, blue edges, small black corners.
    for blk in result["blocks"]:
        ring = blk["ring"]
        ax.fill(ring[:, 0], ring[:, 1], alpha=0.12, color="#87CEEB")
        for s in blk["sides"]:
            ax.plot(s[:, 0], s[:, 1], color="#3aa6d0", lw=1.1, zorder=3)
        c = nodes[blk["corners"]]
        ax.scatter(c[:, 0], c[:, 1], c="black", s=10, zorder=5)

    # Domain boundary (walls / blade): black.
    for b in boundary_ref:
        b = np.asarray(b, float)
        if b.ndim == 2 and len(b) >= 2:
            ax.plot(b[:, 0], b[:, 1], color="black", lw=1.6, zorder=4)

    # ONLY the deletion story: each T-a-deleted wedge streamline drawn RED,
    # and a RED ring around the singularity it emanates from (the two
    # singularities where a streamline is removed). No other markers.
    sing_xy = nodes[irregular] if irregular else np.zeros((0, 2))
    ring_sings = set()
    for arm in wedge_arms:
        arm = np.asarray(arm, float)
        if arm.ndim != 2 or len(arm) < 2:
            continue
        ax.plot(arm[:, 0], arm[:, 1], color="red", lw=2.2, zorder=8)
        if len(sing_xy):
            # The singularity end of the arm: the endpoint nearest a singularity.
            d0 = np.min(np.linalg.norm(sing_xy - arm[0], axis=1))
            d1 = np.min(np.linalg.norm(sing_xy - arm[-1], axis=1))
            tip = arm[0] if d0 <= d1 else arm[-1]
            ring_sings.add(int(np.argmin(np.linalg.norm(sing_xy - tip, axis=1))))
    for si in ring_sings:
        ax.scatter(sing_xy[si, 0], sing_xy[si, 1], facecolors="none",
                   edgecolors="red", s=170, linewidths=2.2, zorder=9)

    ax.set_aspect("equal")
    ax.set_axis_off()
    fig.savefig(OUT_PATH, dpi=200, transparent=True, bbox_inches="tight")
    plt.close(fig)

    print(f"wrote {OUT_PATH}")
    print(f"  blocks:             {len(result['blocks'])}")
    print(f"  deleted streamlines (red): {len(wedge_arms)}")
    print(f"  ringed singularities:      {len(ring_sings)}")


if __name__ == "__main__":
    main()
