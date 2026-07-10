#!/usr/bin/env python3
"""Test different singularity offsets to minimize rejected regions."""

import json
import sys
import subprocess
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent.parent / "domain_partition_2D"))

from create_shroud_prescribed import create_shroud_prescribed
from tmesh_partition import run_tmesh
import unwrap_surface as us


def test_offset(stl_path, offset, verbose=False):
    """Test a specific offset and return metrics."""
    mesh = us.unwrap(stl_path)
    r = float(mesh["r"])
    st = mesh["st"]
    blade_loops = mesh["blade_loops"]
    
    if not blade_loops:
        return None
    
    loop = np.asarray(blade_loops[0], float)
    le_idx = int(np.argmin(loop[:, 0]))
    te_idx = int(np.argmax(loop[:, 0]))
    le = loop[le_idx]
    te = loop[te_idx]
    
    if le_idx < te_idx:
        upper = loop[le_idx:te_idx+1]
        lower = np.vstack([loop[te_idx:], loop[:le_idx+1]])
    else:
        upper = np.vstack([loop[le_idx:], loop[:te_idx+1]])
        lower = loop[te_idx:le_idx+1]
    
    if np.mean(upper[:, 1]) > np.mean(lower[:, 1]):
        suction = upper
        pressure = lower
    else:
        suction = lower
        pressure = upper
    
    suction_mid = suction[len(suction)//2]
    pressure_mid = pressure[len(pressure)//2]
    
    tangent = te - le
    normal = np.array([-tangent[1], tangent[0]])
    normal = normal / (np.linalg.norm(normal) + 1e-8)
    
    prescribed = []
    for i, (pos, name) in enumerate([(le, "LE"), (te, "TE")]):
        sing_pos = pos - offset * normal
        prescribed.append({
            "id": f"sing_{i}",
            "index": -1,
            "separatrix_count": 5,
            "position_st": sing_pos.tolist(),
            "position_st_3d": [r * np.cos(sing_pos[0]/r), r * np.sin(sing_pos[0]/r), sing_pos[1]],
            "blade_parametric": {"location": name, "offset": offset}
        })
    
    for i, (pos, name) in enumerate([(suction_mid, "suction"), (pressure_mid, "pressure")]):
        sing_pos = pos + offset * normal
        prescribed.append({
            "id": f"sing_{i+2}",
            "index": 1,
            "separatrix_count": 3,
            "position_st": sing_pos.tolist(),
            "position_st_3d": [r * np.cos(sing_pos[0]/r), r * np.sin(sing_pos[0]/r), sing_pos[1]],
            "blade_parametric": {"location": name, "offset": offset}
        })
    
    data = {
        "schema_version": "1.0.0",
        "case_id": "T1_9",
        "surface": "shroud",
        "tag": "ta",
        "geometry": {
            "r": r,
            "pitch": mesh.get("pitch", 0),
            "blade_count": 8,
            "periodic": True,
        },
        "singularities": prescribed,
        "invariants": {
            "total_singularities": len(prescribed),
            "index_sum": sum(s["index"] for s in prescribed),
        }
    }
    
    # Write temp file
    tmp_path = Path("/tmp/shroud_prescribed_test.json")
    tmp_path.write_text(json.dumps(data, indent=2))
    
    import shutil
    out_dir = Path("output/T1_9/shroud_stage1/tmesh")
    if out_dir.exists():
        shutil.rmtree(out_dir)
    
    # Run tmesh
    try:
        result = run_tmesh(str(Path(stl_path).resolve()), str(tmp_path), tag=f"offset_{offset:.2f}", verbose=verbose)
        metrics = result.get("metrics", {})
        return {
            "offset": offset,
            "blocks": metrics.get("blocks", 0),
            "rejected": metrics.get("rejected_regions", 0),
            "hanging": metrics.get("hanging_seam_junctions", 0),
            "irregular": metrics.get("irregular_interior_nodes", 0),
            "inverted": metrics.get("inverted_blocks", 0),
            "success": True
        }
    except Exception as e:
        return {
            "offset": offset,
            "success": False,
            "error": str(e)
        }


if __name__ == "__main__":
    stl = sys.argv[1] if len(sys.argv) > 1 else "T1_9_shroud_raw.stl"
    
    print("Testing various offsets for Shroud prescribed singularities...")
    print("=" * 60)
    
    results = []
    for offset in [0.08, 0.10, 0.12, 0.15, 0.18, 0.20, 0.25, 0.30]:
        print(f"\n[test] offset={offset:.2f}")
        res = test_offset(stl, offset, verbose=False)
        results.append(res)
        if res["success"]:
            print(f"  blocks={res['blocks']}, rejected={res['rejected']}, "
                  f"hanging={res['hanging']}, irregular={res['irregular']}")
        else:
            print(f"  FAILED: {res.get('error', 'unknown')}")
    
    print("\n" + "=" * 60)
    print("SUMMARY:")
    print(f"{'offset':>8} {'blocks':>8} {'rejected':>10} {'hanging':>10} {'irregular':>10}")
    print("-" * 60)
    for r in results:
        if r["success"]:
            print(f"{r['offset']:>8.2f} {r['blocks']:>8} {r['rejected']:>10} "
                  f"{r['hanging']:>10} {r['irregular']:>10}")
    
    # Find best
    best = min([r for r in results if r["success"]], key=lambda x: (x["rejected"], x["hanging"]))
    print(f"\nBEST: offset={best['offset']:.2f} with {best['rejected']} rejected, {best['hanging']} hanging")
