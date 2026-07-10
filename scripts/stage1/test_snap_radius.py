#!/usr/bin/env python3
"""Test different snap radius values to reduce T-nodes without increasing rejects."""

import json
import sys
import shutil
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent.parent / "domain_partition_2D"))

from tmesh_partition import run_tmesh
import partition_surface as ps


def test_snap_radius(stl_path, prescribed_path, snap_radius, verbose=False):
    """Test a specific snap radius by temporarily modifying the global."""
    with open(prescribed_path) as f:
        shroud_data = json.load(f)
    prescribed = shroud_data["singularities"]
    
    out_dir = Path("output/T1_9/shroud_stage1/tmesh_test")
    if out_dir.exists():
        shutil.rmtree(out_dir)
    
    # Store original function
    original_snap = ps._snap_separatrix_endpoints
    
    # Create wrapper with custom radius
    def custom_snap(mesh, **kwargs):
        return original_snap(mesh, radius=snap_radius)
    
    # Monkey-patch
    ps._snap_separatrix_endpoints = custom_snap
    
    try:
        result = run_tmesh(
            stl=stl_path,
            out_dir=out_dir,
            tag=f"snap_{snap_radius:.3f}",
            prescribed_singularities=prescribed,
            verbose=verbose,
        )
        metrics = result.get("metrics", {})
        return {
            "snap_radius": snap_radius,
            "blocks": metrics.get("blocks", 0),
            "rejected": metrics.get("rejected_regions", 0),
            "hanging": metrics.get("hanging_seam_junctions", 0),
            "irregular": metrics.get("irregular_interior_nodes", 0),
            "success": True
        }
    except Exception as e:
        return {
            "snap_radius": snap_radius,
            "success": False,
            "error": str(e)
        }
    finally:
        # Restore original
        ps._snap_separatrix_endpoints = original_snap


if __name__ == "__main__":
    stl = sys.argv[1] if len(sys.argv) > 1 else "T1_9_shroud_raw.stl"
    prescribed = sys.argv[2] if len(sys.argv) > 2 else "output/T1_9/shroud_stage1/shroud_prescribed_ta.json"
    
    print("Testing various snap radius values for Shroud partition...")
    print("=" * 60)
    
    results = []
    for radius in [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.10]:
        print(f"\n[test] snap_radius={radius:.3f}")
        res = test_snap_radius(stl, prescribed, radius, verbose=False)
        results.append(res)
        if res["success"]:
            print(f"  blocks={res['blocks']}, rejected={res['rejected']}, "
                  f"hanging={res['hanging']}, irregular={res['irregular']}")
        else:
            print(f"  FAILED: {res.get('error', 'unknown')}")
    
    print("\n" + "=" * 60)
    print("SUMMARY:")
    print(f"{'radius':>10} {'blocks':>8} {'rejected':>10} {'hanging':>10} {'irregular':>10}")
    print("-" * 60)
    for r in results:
        if r["success"]:
            print(f"{r['snap_radius']:>10.3f} {r['blocks']:>8} {r['rejected']:>10} "
                  f"{r['hanging']:>10} {r['irregular']:>10}")
    
    best = min([r for r in results if r["success"]], key=lambda x: (x["rejected"], x["hanging"]))
    print(f"\nBEST: radius={best['snap_radius']:.3f} with {best['rejected']} rejected, {best['hanging']} hanging")
