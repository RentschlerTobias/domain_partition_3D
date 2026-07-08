#!/usr/bin/env python3
"""Compare the two periodic-seam strategies on the T1_9 hub passage:

  Ansatz W (wrap):   streamlines integrate ACROSS the theta seam (universal
                     cover), helices merge onto limit-orbit rings, blocks are
                     built on a 3-pitch strip and the central pitch extracted.
  Ansatz T (t-mesh): streamlines END on the seam; junctions master-slave
                     symmetrized mod pitch; T-junction-tolerant block
                     extraction; TFI fill (tmesh_partition.py).

ONE run writes to output/T1_9/hub_stage1/periodic_compare/:
  wrap_blocks.png, tmesh_blocks.png, tmesh_tfi.png, compare.txt

The W path replicates partition_surface.partition() step by step (identical
calls/flags) but runs merging/splitter manually so the SAME splitter output can
be fed to (a) the upstream QuadFaceGenerator (= identical W blocks) and (b) the
region extractor, which counts the non-quad regions QuadFaceGenerator drops
silently (the "holes" metric).
"""

import json
import sys
import time
from pathlib import Path

import numpy as np
import torch
from torch_geometric.data import Data

sys.path.insert(0, "/root/repos/domain_partition")
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
from tools.streamlines_to_quad_faces import QuadFaceGenerator  # noqa: E402

STL = "/root/repos/block_structured_meshing/T1_9_hub_raw.stl"
OUT = Path("/root/repos/block_structured_meshing/output/T1_9/hub_stage1/"
           "periodic_compare")


def _seam_conformity_from_corners(bx, used, pitch):
    """Handoff sec 7 snippet: block corners on the two seams, right mapped by
    -pitch, matched by nearest neighbour. Returns (nL, nR, max_dev)."""
    def dl(p):
        return abs(p[0] - 0.407 * p[1])

    def dr(p):
        return abs(p[0] - (0.593 + 0.407 * p[1]))

    L = np.array([bx[i] for i in used if dl(bx[i]) < 0.02])
    R = np.array([bx[i] - [pitch, 0.0] for i in used if dr(bx[i]) < 0.02])
    if not len(L) or not len(R):
        return len(L), len(R), None
    d = np.linalg.norm(L[:, None, :] - R[None, :, :], axis=2)
    return len(L), len(R), float(max(d.min(axis=1).max(), d.min(axis=0).max()))


def run_wrap(stl=STL, out_dir=OUT):
    """Ansatz W, identical to partition_surface.partition() with periodicity on,
    plus the silent-drop (holes) count on the splitter output."""
    ps.set_periodicity(True)
    t0 = time.time()

    mesh, transform = build_dp_data(stl)
    pitch = float(mesh.pitch_norm)
    ff = FrameField(mesh)
    m = detect_singularities(ff.mesh)
    n_sing = int((m.singularities != 0).sum())

    sl = StreamlineGenerator_v2(ff.mesh)
    ps._drop_degenerate_corner_seps(sl.mesh)
    ps._close_helical_streamlines(sl.mesh)
    ps._emit_dock_crossings(sl)
    ps._snap_separatrix_endpoints(sl.mesh, radius=0.045)

    central_path = ps._central_parallelogram(sl.mesh)
    sl.mesh.streamlines = ps._tile_streamlines_periodic(
        sl.mesh.streamlines, pitch)

    merging = StreamlineMerging(sl.mesh, verbose=False)
    splitter = StreamlineIntersectionSplitter(offset_boundingBox=0.05,
                                              num_samples=5)
    updated = splitter.process_streamlines(merging.new_streamlines)

    # identical block path (upstream QuadFaceGenerator on the splitter output)
    fg = QuadFaceGenerator(updated)
    faces, edge_to_streamline, edge_index, nodes = fg.get_data()
    block_mesh = Data(x=nodes, edge_index=edge_index, faces=faces,
                      edge_to_streamline=edge_to_streamline)
    n_strip = faces.shape[1] if faces is not None else 0
    block_mesh = ps._extract_central_blocks(block_mesh, central_path)
    n_blocks = (block_mesh.faces.shape[1]
                if block_mesh.faces is not None else 0)
    print(f"[wrap] central extraction: {n_strip} strip -> {n_blocks} central")

    # holes: regions QuadFaceGenerator drops (len != 4) inside the central pitch
    from matplotlib.path import Path as MplPath
    n2, edges, e2s, n_par = tmf.build_connectivity(updated)
    regions = tmf.extract_regions(edges)
    holes = None
    if regions is not None:
        infos = []
        for face in regions:
            if len(face) < 3 or len(set(face)) != len(face):
                infos.append((list(face), None, 0.0))
                continue
            ring = tmf._ring_of(face, e2s)
            infos.append((list(face), ring, abs(tmf._shoelace(ring))))
        outer = int(np.argmax([i[2] for i in infos]))
        blade_paths = [MplPath(np.asarray(bl, float))
                       for bl in mesh.blade_loops]
        holes = 0
        for ri, (face, ring, _a) in enumerate(infos):
            if ri == outer or len(face) == 4:
                continue
            cen = (ring[:-1].mean(axis=0) if ring is not None
                   else n2[list(set(face))].mean(axis=0))
            if not central_path.contains_point(cen):
                continue
            if any(bp.contains_point(cen) for bp in blade_paths):
                continue
            holes += 1

    bx = block_mesh.x.numpy()
    bf = block_mesh.faces.numpy().T if block_mesh.faces is not None else \
        np.zeros((0, 4), int)
    irr, inv, hi = ps._block_annotations(bx, bf)
    used = np.unique(bf) if len(bf) else []
    nL, nR, dev = _seam_conformity_from_corners(bx, used, pitch)

    runtime = time.time() - t0
    metrics = {
        "approach": "W (wrap + orbit rings + tiling)",
        "singularities": n_sing,
        "blocks": int(len(bf)),
        "irregular_interior_nodes": len(irr),
        "inverted_blocks": len(inv),
        "dropped_nonquad_central_regions": holes,
        "seam_corner_lr": [int(nL), int(nR)],
        "seam_corner_dev": dev,
        "runtime_s": round(runtime, 1),
    }
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    ps._plot_blocks(block_mesh, mesh, out_dir / "wrap_blocks.png")
    print(f"[wrap] {metrics['blocks']} blocks, {len(irr)} irregular, "
          f"{len(inv)} inverted, holes={holes}, {runtime:.0f}s")
    return metrics


def _fmt(v):
    if v is None:
        return "n/a"
    if isinstance(v, float):
        return f"{v:.3e}" if (v != 0 and abs(v) < 1e-2) else f"{v:g}"
    return str(v)


def main(rerun_wrap=False):
    us.set_blade_tip_corners(False)
    cs.set_emanate_outer_corners(True)
    OUT.mkdir(parents=True, exist_ok=True)

    # Ansatz W is slow (~8 min) and unchanged since round 1 -> cache metrics
    cache = OUT / "compare_metrics.json"
    w = None
    if not rerun_wrap and cache.exists():
        w = json.loads(cache.read_text()).get("wrap")
        if w:
            print("[wrap] reusing cached metrics (pass --rerun-wrap to redo)")
    if w is None:
        print("=== Ansatz W (wrap) ===")
        w = run_wrap()

    print("=== Ansatz T-a (hanging seam T-nodes) ===")
    ta = tp.run_tmesh(out_dir=OUT, continue_seam_edges=False,
                      tag="ta")["metrics"]
    print("=== Ansatz T-b (seam edges continued) ===")
    tb = tp.run_tmesh(out_dir=OUT, continue_seam_edges=True,
                      tag="tb")["metrics"]

    def trow(name, wk, tk, wv=None):
        wval = w.get(wk, "n/a") if wk else (wv if wv is not None else "n/a")
        return (name, wval, ta.get(tk), tb.get(tk))

    rows = [
        trow("blocks", "blocks", "blocks"),
        trow("irregular interior nodes", "irregular_interior_nodes",
             "irregular_interior_nodes"),
        trow("inverted blocks", "inverted_blocks", "inverted_blocks"),
        ("inverted TFI cells", "n/a (no TFI)",
         f"{ta['inverted_tfi_cells']}/{ta['total_tfi_cells']}",
         f"{tb['inverted_tfi_cells']}/{tb['total_tfi_cells']}"),
        trow("non-quad regions (holes / rejects)",
             "dropped_nonquad_central_regions", "rejected_regions"),
        trow("hanging seam junctions", None, "hanging_seam_junctions", 0),
        trow("edge conformity max dev (grid pts)", None,
             "edge_conformity_dev", "n/a (no shared sampling)"),
        trow("wall clustering ratio (target 5)", None,
             "wall_cluster_ratio_mean", "n/a"),
        trow("seam corners L/R", "seam_corner_lr", "seam_corner_lr"),
        trow("seam corner max dev", "seam_corner_dev", "seam_corner_dev"),
        trow("seam TFI node max dev", None, "seam_tfi_dev", "n/a"),
        trow("singularities", "singularities", "singularities"),
        trow("runtime [s]", "runtime_s", "runtime_s"),
    ]
    lines = ["Periodic-seam strategy comparison, T1_9 hub passage",
             "config: blade_tip_corners=False, emanate_outer_corners=True",
             "",
             f"{'metric':<40}{'W (wrap)':<20}{'T-a (hanging)':<20}"
             f"T-b (continued)"]
    lines.append("-" * 100)
    for name, a, b, c in rows:
        lines.append(f"{name:<40}{_fmt(a):<20}{_fmt(b):<20}{_fmt(c)}")
    txt = "\n".join(lines)
    (OUT / "compare.txt").write_text(txt + "\n")
    (OUT / "compare_metrics.json").write_text(
        json.dumps({"wrap": w, "tmesh_ta": ta, "tmesh_tb": tb}, indent=2))
    print(txt)
    print(f"wrote {OUT}/compare.txt")


if __name__ == "__main__":
    main(rerun_wrap="--rerun-wrap" in sys.argv)
