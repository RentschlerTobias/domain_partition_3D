#!/usr/bin/env python3
"""Test continue_seam_edges option to reduce rejected regions."""

import json
import sys
import shutil
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent.parent / "domain_partition_2D"))

from tmesh_partition import run_tmesh


def test_continue_seam(stl_path, prescribed_path, flat_tol, verbose=False):
    """Test with continue_seam_edges=True."""
    with open(prescribed_path) as f:
        shroud_data = json.load(f)
    prescribed = shroud_data["singularities"]
    
    out_dir = Path("output/T1_9/shroud_stage1/tmesh_test")
    if out_dir.exists():
        shutil.rmtree(out_dir)
    
    try:
        result = run_tmesh(
            stl=stl_path,
            out_dir=out_dir,
            tag=f"cont_{flat_tol:.1f}",
            prescribed_singularities=prescribed,
            verbose=verbose,
            flat_tol_deg=flat_tol,
            continue_seam_edges=True,
            max_rounds=5,
        )
        metrics = result.get("metrics", {})
        return {
            "flat_tol": flat_tol,
            "blocks": metrics.get("blocks", 0),
            "rejected": metrics.get("rejected_regions", 0),
            "hanging": metrics.get("hanging_seam_junctions", 0),
            "irregular": metrics.get("irregular_interior_nodes", 0),
            "success": True
        }
    except Exception as e:
        return {
            "flat_tol": flat_tol,
            "success": False,
            "error": str(e)
        }


if __name__ == "__main__":
    stl = sys.argv[1] if len(sys.argv) > 1 else "T1_9_shroud_raw.stl"
    prescribed = sys.argv[2] if len(sys.argv) > 2 else "output/T1_9/shroud_stage1/shroud_prescribed_ta.json"
    
    print("Testing continue_seam_edges=True for Shroud partition...")
    print("=" * 60)
    
    results = []
    for flat_tol in [15.0, 25.0, 35.0, 45.0]:
        print(f"\n[test] flat_tol={flat_tol:.1f}, continue_seam=True")
        res = test_continue_seam(stl, prescribed, flat_tol, verbose=False)
        results.append(res)
        if res["success"]:
            print(f"  blocks={res['blocks']}, rejected={res['rejected']}, "
                  f"hanging={res['hanging']}, irregular={res['irregular']}")
        else:
            print(f"  FAILED: {res.get('error', 'unknown')}")
    
    print("\n" + "=" * 60)
    print("SUMMARY:")
    print(f"{'flat_tol':>10} {'blocks':>8} {'rejected':>10} {'hanging':>10} {'irregular':>10}")
    print("-" * 60)
    for r in results:
        if r["success"]:
            print(f"{r['flat_tol']:>10.1f} {r['blocks']:>8} {r['rejected']:>10} "
                  f"{r['hanging']:>10} {r['irregular']:>10}")
