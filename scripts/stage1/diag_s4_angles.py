#!/usr/bin/env python3
"""Diagnostic: Analyze streamline angles at S4 after merge."""
import sys, time
from pathlib import Path
import numpy as np

_DP2D = Path(__file__).resolve().parent.parent.parent.parent / "domain_partition_2D"
sys.path.insert(0, str(_DP2D))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import unwrap_surface as us
import clean_separatrix as cs
import partition_surface as ps
import tmesh_faces as tmf
from dp_adapter import build_dp_data
from tools import FrameField, StreamlineGenerator_v2
from tools.singularity_detector import detect_singularities
from tools.streamline_merging import StreamlineMerging
from tools.streamline_intersection_splitter import StreamlineIntersectionSplitter

assert hasattr(StreamlineMerging.__init__, '__wrapped__') or \
       'find_missed_streamline_endpoints' in StreamlineMerging.__init__.__code__.co_names, \
       "StreamlineMerging monkeypatch missing! partition_surface must be imported."

_ROOT = Path(__file__).resolve().parent.parent.parent
STL = str(_ROOT / "T1_9_shroud_raw.stl")

def run_diagnostic():
    ps.set_periodic(True)
    ps.set_tile_periodic(False)

    print("=" * 70)
    print("DIAGNOSTIC: Xiao method streamline angles at S4")
    print("=" * 70)

    t0 = time.time()
    mesh, _tf = build_dp_data(STL)
    ff = FrameField(mesh)
    m = detect_singularities(ff.mesh)

    sl = StreamlineGenerator_v2(ff.mesh)
    ps._drop_degenerate_corner_seps(sl.mesh)
    ps._snap_separatrix_endpoints(sl.mesh, radius=0.045)

    # ── Run StreamlineMerging with verbose=True ─────────────────────────────────
    merging = StreamlineMerging(sl.mesh, verbose=True)

    # ── Collect termination-node info ─────────────────────────────────────────
    mask_c0 = sl.mesh.x[:, 2] == 0
    c0_nodes = sl.mesh.x[mask_c0, 0:2]
    n_c0 = len(c0_nodes)

    singularity_coords = np.array([
        np.asarray(sl.mesh.singularities_coords[s], float)
        for s in sl.mesh.singularities_coords
    ])
    n_sing = len(singularity_coords)

    print(f"\n{'='*70}")
    print(f"INTERNAL SINGULARITY POSITIONS AFTER SNAP")
    print(f"{'='*70}")
    print(f"  n_c0_boundary_nodes = {n_c0}")
    print(f"  n_internal_singularities = {n_sing}")

    # Map singularity face-id -> internal index in StreamlineMerging (offset n_c0)
    sing_fids = list(sl.mesh.singularities_coords.keys())

    for internal_idx in range(n_sing):
        fid = sing_fids[internal_idx]
        coords = singularity_coords[internal_idx]
        print(f"  internal_idx={internal_idx}  fid={fid}  coords=({coords[0]:.6f}, {coords[1]:.6f})")

    # ── Find S4: user coords (0.6589, 0.6647) in NORMALIZED [0,1]^2 ───────────
    user_s4 = np.array([0.6589, 0.6647], dtype=float)
    tol = 0.05

    s4_internal_idx = None
    for internal_idx in range(n_sing):
        if np.linalg.norm(singularity_coords[internal_idx] - user_s4) < tol:
            s4_internal_idx = internal_idx
            break

    print(f"\n{'='*70}")
    print(f"S4 IDENTIFICATION")
    print(f"{'='*70}")
    if s4_internal_idx is not None:
        print(f"  User S4 (0.6589, 0.6647) -> internal singularity index = {s4_internal_idx}")
        print(f"  S{s4_internal_idx} coords = ({singularity_coords[s4_internal_idx][0]:.6f}, "
              f"{singularity_coords[s4_internal_idx][1]:.6f})")
    else:
        print(f"  WARNING: No internal singularity within tol={tol} of (0.6589, 0.6647)")
        # Try larger tolerance
        best_dist = np.inf
        for internal_idx in range(n_sing):
            d = np.linalg.norm(singularity_coords[internal_idx] - user_s4)
            if d < best_dist:
                best_dist = d
                best_idx = internal_idx
        print(f"  Nearest: internal_idx={best_idx}  dist={best_dist:.6f}")
        print(f"  coords = ({singularity_coords[best_idx][0]:.6f}, {singularity_coords[best_idx][1]:.6f})")

    # ── Streamlines touching S4 after merge ───────────────────────────────────
    print(f"\n{'='*70}")
    print(f"STREAMLINES TOUCHING S4 (internal_idx={s4_internal_idx}) AFTER MERGE")
    print(f"{'='*70}")

    # Access internal Streamlines dict from merging
    SL = merging.Streamlines
    SG = merging.Singularities

    if s4_internal_idx is None:
        print("  (Cannot proceed without S4 internal index)")
        return

    s4_data = SG[s4_internal_idx]
    streamline_in_list = s4_data["streamline_in"]
    streamline_out_list = s4_data["streamline_out"]

    print(f"\n  S{s4_internal_idx} has {len(streamline_in_list)} incoming, {len(streamline_out_list)} outgoing")
    print(f"  (after merge, some may have been merged into new IDs)\n")

    # Collect all streamlines that have S4 as in or out
    touching = []
    for sl_id, sl_data in SL.items():
        if sl_data["singularity_in"] == s4_internal_idx or sl_data["singularity_out"] == s4_internal_idx:
            touching.append(sl_id)

    touching = sorted(set(touching))
    print(f"  Total streamlines touching S4: {len(touching)}")

    # For each touching streamline, determine which end is at S4
    for sl_id in touching:
        sl_data = SL[sl_id]
        coords = np.array(sl_data["coords"])

        # Determine which end is at S4 and which is the other end
        sin = sl_data["singularity_in"]
        sout = sl_data["singularity_out"]

        # Compute angle at S4 end
        angle_at_s4 = None
        other_end = None
        s4_is_out = False

        if sin == s4_internal_idx:
            # S4 is the IN end: angle_in is at the END of the streamline
            angle_at_s4 = float(sl_data["angle_in"])
            # Other end is the start (sout, or boundary)
            other_coords = coords[0]
            if sout is not None and not SG[sout]["is_boundary"]:
                other_type = f"S{sout}"
            elif sout is not None:
                other_type = "boundary"
            else:
                other_type = "open/dangling"
            s4_is_out = False
        elif sout == s4_internal_idx:
            # S4 is the OUT end: angle_out is at the START of the streamline
            angle_at_s4 = float(sl_data["angle_out"])
            # Other end is the end
            other_coords = coords[-1]
            if sin is not None and not SG[sin]["is_boundary"]:
                other_type = f"S{sin}"
            elif sin is not None:
                other_type = "boundary"
            else:
                other_type = "open/dangling"
            s4_is_out = True

        # Compute arclength
        if coords.ndim == 2 and len(coords) >= 2:
            seg_diffs = np.diff(coords, axis=0)
            arclen = float(np.linalg.norm(seg_diffs, axis=1).sum())
        else:
            arclen = 0.0

        # Non-degenerate: len >= 10 points
        non_degenerate = len(coords) >= 10

        print(f"  sl_id={sl_id:3d}  len={len(coords):4d}  arclen={arclen:.4f}  "
              f"angle_S4={angle_at_s4:.6f} rad ({np.degrees(angle_at_s4):.2f} deg)  "
              f"other_end=({other_coords[0]:.4f}, {other_coords[1]:.4f}) [{other_type}]  "
              f"non_deg={non_degenerate}")

    # ── Summary table ──────────────────────────────────────────────────────────
    print(f"\n{'='*70}")
    print(f"SUMMARY TABLE: Streamlines touching S4")
    print(f"{'='*70}")
    print(f"{'sl_id':>6} {'n_pts':>6} {'arclen':>8} {'angle_rad':>10} {'angle_deg':>9} "
          f"{'other_end_x':>10} {'other_end_y':>10} {'target':>8} {'non_deg':>7}")
    print("-" * 80)
    for sl_id in touching:
        sl_data = SL[sl_id]
        coords = np.array(sl_data["coords"])
        sin = sl_data["singularity_in"]
        sout = sl_data["singularity_out"]

        if sin == s4_internal_idx:
            angle_at_s4 = float(sl_data["angle_in"])
            other_coords = coords[0]
            target = f"S{sout}" if sout and not SG[sout]["is_boundary"] else ("bnd" if sout else "?")
            s4_is_out = False
        elif sout == s4_internal_idx:
            angle_at_s4 = float(sl_data["angle_out"])
            other_coords = coords[-1]
            target = f"S{sin}" if sin and not SG[sin]["is_boundary"] else ("bnd" if sin else "?")
            s4_is_out = True

        if coords.ndim == 2 and len(coords) >= 2:
            arclen = float(np.linalg.norm(np.diff(coords, axis=0), axis=1).sum())
        else:
            arclen = 0.0
        non_deg = len(coords) >= 10

        print(f"{sl_id:6d} {len(coords):6d} {arclen:8.4f} {angle_at_s4:10.6f} {np.degrees(angle_at_s4):9.2f} "
              f"{other_coords[0]:10.4f} {other_coords[1]:10.4f} {target:>8} {str(non_deg):>7}")

    print(f"\n  Total streamlines touching S4: {len(touching)}")
    print(f"  Non-degenerate count: {sum(1 for sl_id in touching if len(SL[sl_id]['coords']) >= 10)}")
    print(f"\n  (Angles computed at S4 end: angle_in for incoming, angle_out for outgoing)")
    print(f"  (Non-degenerate = 10+ points on streamline)")
    print(f"\n  S4 user coords (0.6589, 0.6647) -> internal singularity index = {s4_internal_idx}")

    elapsed = time.time() - t0
    print(f"\n  Diagnostic completed in {elapsed:.1f}s")


if __name__ == "__main__":
    run_diagnostic()
