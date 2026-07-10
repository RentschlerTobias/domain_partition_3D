#!/usr/bin/env python3
"""
Stage 1 / Slide 06: Frame field with highlighted singularities.

Single-panel figure (saved as a single PNG):
  - Frame field: many small crosses (4 directions) covering the domain,
    drawn in light blue (#87CEEB). The cross directions are recovered
    from the representative vector via
    ``base = arctan2(u[:,1], u[:,0]) / 4`` and ``ang = base + k*π/2``
    for k=0,1,2,3 (the full 4-RoSy cross).
  - Singularities: colored dots overlaid on the frame field at the
    barycenter of each singular face.
      * Positive Poincaré index (+1): red dots
      * Negative Poincaré index (-1): blue dots (square marker)

Reuses the cross-field + singularity-detection logic from
scripts/stage1/field_diagnostic.py without modifying it. The cross-field
is computed by ``FrameField`` (which sets ``mesh.u`` as the (cos4θ, sin4θ)
representation field). The Poincaré index per face is computed by
``detect_singularities`` (which sets ``mesh.singularities``).

Visual recipe:
  - Mesh triangulation background (light gray, same as slide02/03/04/05).
  - Outer + inner contours (black, same as slide02/03/04/05).
  - Frame field: 4 small orthogonal crosses per node (light blue,
    scale=40, no arrow heads) — the full 4-RoSy cross at every node.
  - Singularity markers: red circles (positive index) + blue squares
    (negative index), s=60, zorder=5.
  - No axes, no labels, no titles, no legends.
  - Transparent background (alpha channel).
  - Portrait orientation (figsize=(4, 6)).
  - DPI >= 200.

Output: plots/slide06_singularities.png
"""

import sys
from pathlib import Path

# Make scripts/stage1/ importable so we can reuse field_diagnostic.py logic
# (FrameField, build_dp_data, unwrap_surface, detect_singularities) as-is.
SCRIPT_DIR = Path(__file__).resolve().parent
STAGE1_DIR = SCRIPT_DIR / "stage1"
if str(STAGE1_DIR) not in sys.path:
    sys.path.insert(0, str(STAGE1_DIR))

REPO_ROOT = SCRIPT_DIR.parent
STL_PATH = REPO_ROOT / "T1_9_hub_raw.stl"
OUT_PATH = REPO_ROOT / "plots" / "step06_singularities.png"
# Final partition singularities (the 4 prescribed ones the ta method keeps),
# as opposed to the ~13 raw frame-field detections.
HUB_MASTER_JSON = REPO_ROOT / "output" / "tmesh_hub" / "master" / "hub_master_ta.json"

import matplotlib

matplotlib.use("Agg")  # headless rendering
import matplotlib.pyplot as plt
import numpy as np
import torch

from dp_adapter import build_dp_data  # noqa: E402
from tools import FrameField  # noqa: E402
from tools.singularity_detector import detect_singularities  # noqa: E402
from unwrap_surface import unwrap  # noqa: E402


def main() -> None:
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)

    # --- 1. Build the mesh + cross-field (reuses field_diagnostic.py logic) ---
    mesh, _transform = build_dp_data(STL_PATH)
    ff = FrameField(mesh)

    # --- 2. Cross field (for the background). detect_singularities also sets
    # mesh.u (the (cos4θ, sin4θ) representation field) used below. Its raw
    # per-face Poincaré index (~13 detections) is NOT used for the markers —
    # we show the 4 FINAL partition singularities instead (loaded from JSON). ---
    m = detect_singularities(ff.mesh)

    # --- 3. Cross-field directions (same formula as field_diagnostic.py L37) ---
    u = m.u.numpy()  # (cos4θ, sin4θ) representation field
    base = np.arctan2(u[:, 1], u[:, 0]) / 4.0  # principal cross angle θ

    # --- 4. Unwrapped (s, t) coordinates for plotting ---
    # build_dp_data does NOT renumber nodes (see dp_adapter.py comment),
    # so node indices map directly into unwrap()["st"].
    st_data = unwrap(STL_PATH)
    st = st_data["st"]
    tris = st_data["tris"]
    loops = st_data["loops"]

    # --- 5. Final partition singularities (4 prescribed), loaded from the
    # hub-master export. These are the ones the ta method actually keeps —
    # unlike the ~13 raw frame-field detections. ---
    # ``position_st`` is stored in the pipeline's NORMALIZED [0,1]×[0,1] mesh
    # frame (= mesh.x[:,0:2]), not the raw unwrap (s,t) frame used for plotting.
    # Both share node ordering, so recover the per-axis linear map
    # normalized -> raw st by fitting mesh.x against st.
    import json
    if HUB_MASTER_JSON.exists():
        sing_json = json.loads(HUB_MASTER_JSON.read_text())["singularities"]
        pos_norm = np.array([s["position_st"] for s in sing_json], dtype=float)
        sing_idx = np.array([int(s["index"]) for s in sing_json])
        mx = mesh.x[:, 0:2].numpy()
        cx = np.polyfit(mx[:, 0], st[:, 0], 1)   # normalized x -> raw s
        cy = np.polyfit(mx[:, 1], st[:, 1], 1)   # normalized y -> raw t
        sc = np.column_stack([np.polyval(cx, pos_norm[:, 0]),
                              np.polyval(cy, pos_norm[:, 1])])
    else:
        print(f"[warn] {HUB_MASTER_JSON} missing — no singularity markers drawn. "
              f"Run the ta hub-master export first.")
        sc = np.zeros((0, 2))
        sing_idx = np.zeros(0, dtype=int)
    pos = sing_idx > 0

    # --- 6. Plot ---
    fig, ax = plt.subplots(figsize=(4, 6))  # portrait

    # Mesh triangulation background: very subtle, light gray (same as slide02/03/04/05).
    ax.triplot(st[:, 0], st[:, 1], tris, lw=0.1, color="0.93")

    # Outer contour: thick black (same as slide02/03/04/05).
    outer_loop = loops[0]
    ring = st[outer_loop + [outer_loop[0]]]
    ax.plot(ring[:, 0], ring[:, 1], lw=2.0, color="black")

    # Inner (blade) loops: also thick black so the blade-shaped hole is visible.
    for inner in loops[1:]:
        inner_ring = st[inner + [inner[0]]]
        ax.plot(inner_ring[:, 0], inner_ring[:, 1], lw=2.0, color="black")

    # Frame field: 4 small orthogonal crosses per node (the full 4-RoSy cross).
    # Back-mapping from the representative vector u to the 4 cross directions:
    #   base = arctan2(u[:,1], u[:,0]) / 4
    #   ang_k = base + k * π/2  for k = 0, 1, 2, 3
    for k in range(4):
        ang = base + k * (np.pi / 2)
        ax.quiver(
            st[:, 0], st[:, 1],
            np.cos(ang), np.sin(ang),
            color="#87CEEB",  # light blue
            scale=40,
            width=0.0015,
            headwidth=0,
            headlength=0,
            pivot="mid",
            alpha=0.9,
        )

    # Final singularity markers: red open circle rings (same style as the
    # non-4-valence markers in Step 9 / plot_blocks: facecolors none, red edge).
    ax.scatter(
        sc[:, 0], sc[:, 1],
        facecolors="none", edgecolors="red", s=150, linewidths=2.0,
        zorder=6,
    )

    ax.set_aspect("equal")
    ax.set_axis_off()

    fig.savefig(
        OUT_PATH,
        dpi=200,
        transparent=True,
        bbox_inches="tight",
    )
    plt.close(fig)

    print(f"wrote {OUT_PATH}")
    print(f"  total nodes:    {len(st)}")
    print(f"  triangles:      {len(tris)}")
    print(f"  singularities:  {len(sc)} (final, prescribed) "
          f"(+1: {int(pos.sum())}, -1: {int((~pos).sum())})")


if __name__ == "__main__":
    main()
