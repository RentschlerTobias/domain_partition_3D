#!/usr/bin/env python3
"""Combined overview PNG of the current Stage-1 hub partition result.

One figure, four panels:
  1. unwrapped mesh + corner types (outer vs blade-tip LE/TE)
  2. cross-field singularities + traced separatrices
  3. quad block partition, annotated (irregular nodes / high-aspect blocks)
  4. validation metrics text

Run:  MESH_SCRATCH=/tmp/dp_scratch /root/venv/bin/python scripts/stage1/plot_results.py
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import partition_surface as ps
import clean_separatrix as cs
from unwrap_surface import unwrap


def main(stl, out_png, annihilate=False):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cs.set_annihilate_pairs(annihilate)
    block_mesh, mesh, tf = ps.partition(stl)
    xy = mesh.x[:, 0:2].numpy()
    tris = mesh.faces.T.numpy()
    ct = mesh.corner_type.numpy()
    n_b = len(mesh.streamlines) - len(mesh.separatrices)
    sc = xy[tris[mesh.singularities.numpy() != 0]].mean(axis=1)

    bx = block_mesh.x.numpy()
    bf = block_mesh.faces.numpy().T
    irr, inv, asp = ps._block_annotations(bx, bf)
    v = ps.QuadPartitionValidator(block_mesh, mesh, frame_field=mesh.u, strict=True)
    valid = v.is_valid()
    soft = v.passes_soft_thresholds()
    qs = v.quality_score()

    fig, ax = plt.subplots(2, 2, figsize=(15, 13))

    # 1. unwrapped mesh + corner types
    a = ax[0, 0]
    a.triplot(xy[:, 0], xy[:, 1], tris, lw=0.2, color="0.8")
    a.scatter(xy[ct == 0, 0], xy[ct == 0, 1], c="red", s=55, zorder=5,
              label="outer corner")
    a.scatter(xy[ct == 1, 0], xy[ct == 1, 1], c="magenta", s=110, marker="*",
              zorder=6, label="blade-tip LE/TE")
    a.set_title("1. Unwrapped hub mesh + corner types")
    a.legend(loc="upper left", fontsize=9)
    a.set_aspect("equal")

    # 2. singularities + separatrices
    a = ax[0, 1]
    a.triplot(xy[:, 0], xy[:, 1], tris, lw=0.1, color="0.93")
    for i, s in enumerate(mesh.streamlines):
        s = np.asarray(s, float)
        if s.ndim != 2 or len(s) < 2:
            continue
        a.plot(s[:, 0], s[:, 1], "0.45" if i < n_b else "C3",
               lw=1.0 if i < n_b else 1.5)
    a.scatter(sc[:, 0], sc[:, 1], c="blue", s=70, zorder=6, label="singularity")
    a.scatter(xy[ct == 1, 0], xy[ct == 1, 1], c="magenta", s=110, marker="*",
              zorder=7, label="blade-tip")
    a.set_title(f"2. Cross field: {int((mesh.singularities.numpy()!=0).sum())} "
                f"singularities, {len(mesh.separatrices)} separatrices "
                f"(red), boundary (grey)")
    a.legend(loc="upper left", fontsize=9)
    a.set_aspect("equal")

    # 3. block partition annotated
    a = ax[1, 0]
    a.triplot(xy[:, 0], xy[:, 1], tris, lw=0.1, color="0.93")
    for qi, quad in enumerate(bf):
        ring = bx[list(quad) + [quad[0]]]
        if qi in inv:
            a.fill(ring[:, 0], ring[:, 1], color="red", alpha=0.4)
            a.plot(ring[:, 0], ring[:, 1], "red", lw=2)
        elif qi in asp:
            a.fill(ring[:, 0], ring[:, 1], color="orange", alpha=0.3)
            a.plot(ring[:, 0], ring[:, 1], "darkorange", lw=1.8)
        else:
            a.fill(ring[:, 0], ring[:, 1], alpha=0.25)
            a.plot(ring[:, 0], ring[:, 1], "C0", lw=1.3)
    a.scatter(bx[:, 0], bx[:, 1], c="k", s=10, zorder=5)
    if irr:
        a.scatter(bx[irr, 0], bx[irr, 1], facecolors="none", edgecolors="red",
                  s=180, linewidths=2, zorder=6)
    a.set_title(f"3. Quad block partition: {bf.shape[0]} blocks "
                f"(irregular={len(irr)}, inverted={len(inv)}, aspect>10={len(asp)})")
    a.set_aspect("equal")

    # 4. metrics text
    a = ax[1, 1]
    a.axis("off")
    lines = [
        "VALIDATION (QuadPartitionValidator, strict)",
        "",
        f"  is_valid (hard):        {valid}",
        f"  passes_soft_thresholds: {soft}",
        "",
        f"  blocks:                 {bf.shape[0]}",
        f"  block corner nodes:     {bx.shape[0]}",
        f"  irregular inner nodes:  {len(irr)}   (target 0)",
        f"  inverted quads:         {len(inv)}",
        f"  high-aspect (>10):      {len(asp)}",
        "",
        f"  scaled_jacobian_min:    {qs['scaled_jacobian_min']:.3f}   (>0.3 ok)",
        f"  scaled_jacobian_mean:   {qs['scaled_jacobian_mean']:.3f}   (>0.85 ok)",
        f"  interior angles:        [{qs['min_interior_angle']:.1f}, "
        f"{qs['max_interior_angle']:.1f}] deg",
        f"  edge-length ratio max:  {qs['edge_length_ratio_max']:.2f}   (<10 ok)",
        f"  singularity efficiency: {qs.get('singularity_efficiency')}",
        "",
        "Open: Euler -1 (1 net irregular node); interior cross",
        "singularities topologically capped by wall BC.",
        "Levers: pitchwise periodicity / accept & go 3D.",
    ]
    a.text(0.02, 0.98, "\n".join(lines), va="top", ha="left", fontsize=12,
           family="monospace", transform=a.transAxes)

    mode = "WITH singularity pair annihilation" if annihilate \
        else "no annihilation (default)"
    fig.suptitle(f"Stage-1 Hub Quad-Domain Partition - {mode}",
                 fontsize=16, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.985])
    fig.savefig(out_png, dpi=130, bbox_inches="tight")
    print(f"wrote {out_png}")


if __name__ == "__main__":
    stl = sys.argv[1] if len(sys.argv) > 1 else \
        "/root/repos/block_structured_meshing/T1_9_hub_raw.stl"
    out = Path("/root/repos/block_structured_meshing/output/T1_9/hub_stage1")
    out.mkdir(parents=True, exist_ok=True)
    main(stl, out / "results_overview.png", annihilate=False)
    main(stl, out / "results_overview_annihilated.png", annihilate=True)
