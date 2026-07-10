#!/usr/bin/env python3
"""Verify S4 separatrix fix by inspecting post-merging streamlines directly.

The singularity analysis in run_xiao_only.py uses mesh.separatrices (pre-merging),
so it doesn't show the new streamlines added by _ensure_separatrix_count.
This script inspects merging.new_streamlines directly to verify the fix.
"""
import sys
from pathlib import Path

import numpy as np

_DP2D = Path(__file__).resolve().parent.parent.parent.parent / "domain_partition_2D"
sys.path.insert(0, str(_DP2D))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import partition_surface as ps
from dp_adapter import build_dp_data
from tools import FrameField, StreamlineGenerator_v2
from tools.singularity_detector import detect_singularities
from tools.streamline_merging import StreamlineMerging

ps.set_periodic(True)
ps.set_tile_periodic(False)

_ROOT = Path(__file__).resolve().parent.parent.parent
STL = str(_ROOT / "T1_9_shroud_raw.stl")

mesh, _tf = build_dp_data(STL)
ff = FrameField(mesh)
m = detect_singularities(ff.mesh)
sl = StreamlineGenerator_v2(ff.mesh)
ps._drop_degenerate_corner_seps(sl.mesh)
ps._snap_separatrix_endpoints(sl.mesh, radius=0.045)

merging = StreamlineMerging(sl.mesh, verbose=False)

# Singularity coordinates from mesh.singularities_coords
sing_coords = {fid: np.array(coord) for fid, coord in sl.mesh.singularities_coords.items()}
expected = sl.mesh.expected_separatrices

# Sort singularities by x-coordinate for S1, S2, ... labeling
sorted_sings = sorted(sing_coords.items(), key=lambda kv: kv[1][0])
sing_labels = {fid: f"S{i+1}" for i, (fid, _) in enumerate(sorted_sings)}

print("=" * 60)
print("POST-MERGING SINGULARITY ANALYSIS")
print("=" * 60)
print()

# For each singularity, count streamlines that start or end near it
results = {}
for fid, coords in sing_coords.items():
    label = sing_labels[fid]
    exp = expected.get(fid, "?")
    seps = []
    for i, sl_coords in enumerate(merging.new_streamlines):
        c = np.asarray(sl_coords, float)
        if len(c) < 2:
            continue
        start, end = c[0], c[-1]
        # Check if start or end is near this singularity
        d_start = np.linalg.norm(start - coords)
        d_end = np.linalg.norm(end - coords)
        if d_start < 0.02 or d_end < 0.02:
            # Determine the other end
            if d_start < d_end:
                other = end
                source = "start"
            else:
                other = start
                source = "end"
            # Check if other end is near another singularity
            target_sing = None
            for fid2, coords2 in sing_coords.items():
                if fid2 == fid:
                    continue
                if np.linalg.norm(other - coords2) < 0.02:
                    target_sing = sing_labels[fid2]
                    break
            seps.append((i, source, other, target_sing))
    results[label] = {"fid": fid, "coords": coords, "expected": exp, "seps": seps}

# Print results
for label in sorted(results.keys(), key=lambda l: int(l[1:])):
    r = results[label]
    c = r["coords"]
    exp = r["expected"]
    seps = r["seps"]
    print(f"{label}: ({c[0]:.4f}, {c[1]:.4f}) expected={exp}")
    print(f"  Separatrices ({len(seps)} total):")
    for idx, source, other, target in seps:
        if target:
            print(f"    sl {idx}: → {target} (direct sing↔sing)")
        else:
            print(f"    sl {idx}: → boundary ({other[0]:.4f},{other[1]:.4f})")
    print()

# Check S4 specifically
s4 = results.get("S4")
if s4:
    print("=" * 60)
    print("S4 VERIFICATION")
    print("=" * 60)
    print(f"S4 expected: {s4['expected']}")
    print(f"S4 actual: {len(s4['seps'])}")
    boundary_count = sum(1 for _, _, _, t in s4['seps'] if t is None)
    sing_count = sum(1 for _, _, _, t in s4['seps'] if t is not None)
    print(f"S4 to boundary: {boundary_count}")
    print(f"S4 to singularity: {sing_count}")
    print()
    if len(s4['seps']) >= s4['expected']:
        print(f"PASS: S4 has {len(s4['seps'])} separatrices (expected {s4['expected']})")
    else:
        print(f"FAIL: S4 has {len(s4['seps'])} separatrices (expected {s4['expected']})")
    print()
    # Check for S3 connection (should not be a new one)
    s3_connections = [s for s in s4['seps'] if s[3] == 'S3']
    if s3_connections:
        print(f"WARNING: S4 has {len(s3_connections)} connection(s) to S3")
    else:
        print("PASS: No incorrect S4→S3 connections")
