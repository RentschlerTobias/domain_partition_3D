#!/usr/bin/env python3
"""Xiao + custom clipping (pure) vs. Ansatz T-a (= same pipeline + seam
postprocessing).

User hypothesis to demonstrate: T-a is exactly the Xiao streamline
simplification with the project's own clipping (c0-corner separatrices, no
LE/TE tip corners) plus a seam postprocessing -- without it, the missing
cross-seam streamlines leave two 3-cornered wedge blocks at the periodic
walls; T-a deletes the two shorter wedge arms (and symmetrizes the junction
sets), which yields the clean block structure.

Both runs use the IDENTICAL periodic cross field (set_periodic(True),
set_tile_periodic(False)); they differ only in the block-stage
postprocessing:

  X   : ... -> snap (Xiao clipping)                        -> blocks
  T-a : ... -> snap -> kink fixes -> wedge collapse
            -> master-slave seam symmetrization            -> blocks + TFI

Writes to output/T1_9/hub_stage1/xiao_vs_ta/: xiao_blocks.png, ta_blocks.png,
xiao_vs_ta.png (side by side, dropped wedge arms dashed orange in the X
panel), compare.txt.
"""

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, "/home/t1dde/Duty/projects/domain_partition/domain_partition_2D")
sys.path.insert(0, str(Path(__file__).resolve().parent))

import unwrap_surface as us                                    # noqa: E402
import clean_separatrix as cs                                  # noqa: E402
import partition_surface as ps                                 # noqa: E402
import tmesh_faces as tmf                                      # noqa: E402
import tmesh_partition as tp                                   # noqa: E402
from dp_adapter import build_dp_data                           # noqa: E402
from tools import FrameField, StreamlineGenerator_v2           # noqa: E402
from tools.singularity_detector import detect_singularities    # noqa: E402
from tools.streamline_merging import StreamlineMerging         # noqa: E402
from tools.streamline_intersection_splitter import (           # noqa: E402
    StreamlineIntersectionSplitter)

STL = "/home/t1dde/Duty/projects/domain_partition/domain_partition_3D/T1_9_shroud_raw.stl"
OUT = Path("/home/t1dde/Duty/projects/domain_partition/domain_partition_3D/output/T1_9/shroud_stage1/"
           "xiao_vs_ta")


def run_xiao(stl=STL):
    """Pure Xiao + custom clipping: pipeline up to the snap, then blocks.
    No kink fixes, no wedge collapse, no seam symmetrization."""
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
        "approach": "X (Xiao + custom clipping, no seam postprocessing)",
        "singularities": n_sing,
        "blocks": len(result["blocks"]),
        "nonquad_regions": len(result["rejects"]),
        "irregular_interior_nodes": len(irregular),
        "runtime_s": round(time.time() - t0, 1),
    }
    return {"metrics": metrics, "result": result, "mesh": mesh}


def _draw_partition(ax, result, mesh, wedge_arms=None, seam_info=None,
                    title=""):
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
                    label="3-corner seam wedge" if ri == 0 else None)
            cen = r[:-1].mean(axis=0)
            ax.annotate(f"{rej['n_real']}-corner", cen, fontsize=8,
                        color="darkred", ha="center")
    if wedge_arms:
        for arm in wedge_arms:
            arm = np.asarray(arm, float)
            ax.plot(arm[:, 0], arm[:, 1], "--", color="darkorange", lw=2.0,
                    zorder=7)
        ax.plot([], [], "--", color="darkorange", lw=2.0,
                label=f"wedge arms deleted by T-a: {len(wedge_arms)}")
    if seam_info:
        for W in (seam_info["WL"], seam_info["WR"]):
            ax.plot(W[:, 0], W[:, 1], "green", lw=0.8, alpha=0.7)
    ax.set_aspect("equal")
    handles, labels = ax.get_legend_handles_labels()
    if labels:
        ax.legend(loc="upper left", fontsize=8)
    ax.set_title(title, fontsize=10)


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    us.set_blade_tip_corners(False)
    cs.set_emanate_outer_corners(True)
    OUT.mkdir(parents=True, exist_ok=True)

    print("=== X: Xiao + custom clipping (pure) ===")
    x = run_xiao()

    print("=== T-a: + kink fix, wedge collapse, seam symmetrization ===")
    ta = tp.run_tmesh(out_dir=OUT, continue_seam_edges=False, tag="ta")
    wedge_arms = getattr(ta["mesh"], "dropped_wedge_arms", [])

    # side-by-side
    fig, axs = plt.subplots(1, 2, figsize=(20, 8.5))
    _draw_partition(
        axs[0], x["result"], x["mesh"], wedge_arms=wedge_arms,
        title=f"X: Xiao + custom clipping -- "
              f"{x['metrics']['blocks']} blocks, "
              f"{x['metrics']['nonquad_regions']} non-quad seam wedges")
    _draw_partition(
        axs[1], ta["result"], ta["mesh"], seam_info=ta["seam_info"],
        title=f"T-a: + wedge collapse & master-slave seam -- "
              f"{ta['metrics']['blocks']} blocks, "
              f"{ta['metrics']['rejected_regions']} non-quad")
    out = OUT / "xiao_vs_ta.png"
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")

    fig, ax = plt.subplots(figsize=(10, 8))
    _draw_partition(ax, x["result"], x["mesh"], wedge_arms=wedge_arms,
                    title=f"Xiao + custom clipping: "
                          f"{x['metrics']['blocks']} blocks, "
                          f"{x['metrics']['nonquad_regions']} wedges")
    fig.savefig(OUT / "xiao_blocks.png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT}/xiao_blocks.png")

    tm = ta["metrics"]
    lines = [
        "Xiao + custom clipping vs. Ansatz T-a, T1_9 hub passage",
        "identical periodic cross field; difference = T-a seam postprocessing",
        "",
        f"{'metric':<44}{'X (Xiao+clipping)':<22}T-a",
        "-" * 90,
        f"{'blocks':<44}{x['metrics']['blocks']:<22}{tm['blocks']}",
        f"{'non-quad regions (seam wedges)':<44}"
        f"{x['metrics']['nonquad_regions']:<22}{tm['rejected_regions']}",
        f"{'irregular interior nodes':<44}"
        f"{x['metrics']['irregular_interior_nodes']:<22}"
        f"{tm['irregular_interior_nodes']}",
        f"{'hanging seam junctions':<44}{'n/a (no symmetrize)':<22}"
        f"{tm['hanging_seam_junctions']}",
        f"{'seam TFI conformity':<44}{'n/a':<22}{(tm['seam_tfi_dev'] or 0):.1e}",
        f"{'inverted TFI cells':<44}{'n/a':<22}"
        f"{tm['inverted_tfi_cells']}/{tm['total_tfi_cells']}",
        f"{'runtime [s]':<44}{x['metrics']['runtime_s']:<22}"
        f"{tm['runtime_s']}",
        "",
        f"wedge arms deleted by T-a: {len(wedge_arms)} "
        f"(dashed orange in xiao_blocks.png)",
    ]
    txt = "\n".join(lines)
    (OUT / "compare.txt").write_text(txt + "\n")
    (OUT / "compare_metrics.json").write_text(json.dumps(
        {"xiao": x["metrics"], "ta": tm}, indent=2))
    print(txt)
    print(f"wrote {OUT}/compare.txt")


if __name__ == "__main__":
    main()
