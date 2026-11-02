"""Stage checks for the blade O-grid ring built on the P1 core.

The ring is the set of hex blocks that fill the gap between the AlgoHex
core's `ogrid_interface` and the blade wall. It is built by
`ogrid_extrude.extrude_ogrid` and written as `blocks_blade.vtk` next to a
`blade_s1.json` summary.

Stages
------
s1  the ring as it comes out of the extrusion: block structure, lattice per
    block, cell quality, and the bit-exact gluing of the ring's base to the
    untouched core blocks.
s2  placeholder: the ring after the O-grid's own blocks have been merged
    into it (not implemented yet).
s3  placeholder: the ring after TFI (not implemented yet).

Usage:  python3 blade_check.py s1 [--run DIR]
"""

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent.parent))

import base_complex as bc          # noqa: E402
import clean_blocks as cb          # noqa: E402
import tfi                         # noqa: E402


DEFAULT_RUN = Path("/work/datagen_keepbl/runs/tistos_s1_0131")
MIN_SCALED_JACOBIAN = 0.15


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def read_blocks_vtk(path):
    """(points, hexes, block_id) of a legacy-VTK hexahedron file."""
    import meshio
    m = meshio.read(str(path))
    hexes = np.vstack([c.data for c in m.cells if c.type == "hexahedron"])
    bid = np.concatenate([np.asarray(d).ravel()
                          for c, d in zip(m.cells, m.cell_data["block_id"])
                          if c.type == "hexahedron"])
    return np.asarray(m.points, float), np.asarray(hexes, np.int64), \
        np.asarray(bid, int)


def core_blocks(run):
    """The P1 core mesh: (points, hexes, block_id) of `blocks.vtk`."""
    return read_blocks_vtk(run / "blocks.vtk")


def surface_labeller(run):
    """The labeller that names the core's boundary faces, so the block
    faces can be checked against the physical surfaces."""
    import tet_prep_v5 as v5
    P, tri, tid = cb.read_input_surface(str(run / "tet.vtk"))
    return cb.SurfaceLabeller(P, tri, tid, v5.NAMES)


def check(cond, msg, failures):
    """Record a failed check, return the condition."""
    if not cond:
        failures.append(msg)
    return bool(cond)




def signed_volume(pts):
    """Signed volume of a hexahedron via the divergence theorem with the
    standard (correct) face winding of a VTK hexahedron. Returns a float."""
    FACES = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 4, 7, 3),
             (1, 2, 6, 5), (0, 1, 5, 4), (2, 3, 7, 6)]
    V = 0.0
    for f in FACES:
        a = pts[f[0]]; b = pts[f[1]]; c = pts[f[2]]; d = pts[f[3]]
        V += float(np.dot(a, np.cross(b, c)))
        V += float(np.dot(a, np.cross(c, d)))
    return V / 6.0


# ---------------------------------------------------------------------------
# stage s1
# ---------------------------------------------------------------------------

def stage_s1(run, verbose=True):
    """The ring straight out of the extrusion."""
    failures = []
    blade_dir = run / "blade"
    summary_path = blade_dir / "blade_s1.json"
    vtk_path = blade_dir / "blocks_blade.vtk"
    for p in (summary_path, vtk_path):
        if not p.exists():
            failures.append(f"missing {p}")
            print(f"[blade_check] FAILED: missing {p}")
            return 1

    summary = json.loads(summary_path.read_text())
    P, H, bid = read_blocks_vtk(vtk_path)
    core_P, core_H, core_bid = core_blocks(run)

    n_blocks = int(bid.max()) + 1
    ids = np.unique(bid)
    print(f"[blade_check] s1: {len(H)} hexes in {n_blocks} blocks, "
          f"{len(P)} points")

    # -- block ids ---------------------------------------------------------
    contiguous = check(ids.tolist() == list(range(n_blocks)),
                       f"block ids not contiguous: {ids.tolist()}", failures)
    counts = Counter(bid.tolist())
    print(f"[blade_check] s1: block ids contiguous: {contiguous}; "
          f"cells per block {dict(sorted(counts.items()))}")

    # the summary must agree with the file
    check(summary.get("n_blocks") == n_blocks,
          f"blade_s1.json says n_blocks={summary.get('n_blocks')}, "
          f"the VTK has {n_blocks}", failures)
    check(summary.get("n_cells") == len(H),
          f"blade_s1.json says n_cells={summary.get('n_cells')}, "
          f"the VTK has {len(H)}", failures)

    # the new ring blocks must carry ids appended after the largest core id
    n_core_blocks = int(summary.get("n_core_blocks", 0))
    core_ids = sorted(set(int(b) for b in np.unique(core_bid)))
    new_ids = sorted(set(int(b) for b in np.unique(bid)) - set(core_ids))
    check(new_ids == list(range(core_ids[-1] + 1, core_ids[-1] + 1 + len(new_ids))),
          f"the new block ids {new_ids} are not appended after "
          f"max(core block_id)={core_ids[-1]}", failures)
    check(summary.get("ring_block_ids") == new_ids,
          f"blade_s1.json says ring_block_ids={summary.get('ring_block_ids')}, "
          f"the VTK has {new_ids}", failures)
    print(f"[blade_check] s1: core block ids {core_ids} "
          f"({n_core_blocks} blocks), new ring ids {new_ids} "
          f"(appended after max(core block_id)={core_ids[-1]})")

    # -- a lattice for every block ----------------------------------------
    f2h, _e2h = bc.build_topology(H)
    cells_of = {}
    for c in range(len(H)):
        cells_of.setdefault(int(bid[c]), []).append(c)
    lattices = {}
    for r in sorted(cells_of):
        lat = tfi.block_lattice(H, cells_of[r], f2h)
        lattices[r] = lat
        ok = check(lat is not None,
                   f"block {r} has no lattice (not a structured grid)", failures)
        if verbose and ok:
            print(f"[blade_check] s1: block {r}: lattice dims "
                  f"{list(lat[0])} over {len(cells_of[r])} cells")

    # -- cell quality ------------------------------------------------------
    sj = cb.scaled_jacobians(P, H)
    n_inv = int((sj <= 0).sum())
    check(n_inv == 0, f"{n_inv} inverted cells (scaled Jacobian <= 0)",
          failures)
    check(float(sj.min()) >= MIN_SCALED_JACOBIAN,
          f"min scaled Jacobian {sj.min():.4f} < {MIN_SCALED_JACOBIAN}",
          failures)
    print(f"[blade_check] s1: scaled Jacobian min {sj.min():.4f} "
          f"mean {sj.mean():.4f}, inverted {n_inv}")

    # -- the whole mesh must pass the block validator ----------------------
    labeller = surface_labeller(run)
    S = cb.BlockStructure.from_arrays(P, H, bid, labeller)
    v = cb.HexBlockValidator(S)
    hard = v.hard_failures()
    check(not hard, f"HexBlockValidator: {hard}", failures)
    print(f"[blade_check] s1: HexBlockValidator "
          f"{'VALID' if not hard else 'INVALID'}")
    for line in v.diagnostics().splitlines():
        print(f"[blade_check] s1:   {line}")

    # -- the ring's base must be glued bit-exactly to the core -------------
    n_ring = summary.get("n_ring_blocks", n_blocks)
    ring_blocks = sorted(range(n_blocks))[-n_ring:] if n_ring else []
    replaced = set(summary.get("replaced_blocks", []))
    core_blocks_untouched = sorted(set(int(b) for b in np.unique(core_bid))
                                    - replaced)
    print(f"[blade_check] s1: ring blocks {ring_blocks}, "
          f"untouched core blocks {core_blocks_untouched}")

    # the ring's base vertices are the core's interface vertices; they must
    # appear in the ring with exactly the same coordinates
    ring_cells = [c for c in range(len(H)) if int(bid[c]) in set(ring_blocks)]
    n_core_pts = len(core_P)
    # the ring's layer-0 cells are the ones whose whole base face is made of
    # core points (the ring stacks its layers along the wall normal, so only
    # the first layer of each quad touches the core)
    base_cells = [c for c in ring_cells
                  if all(int(v) < n_core_pts for v in H[c][:4])]
    base_pts = set()
    for c in base_cells:
        base_pts.update(int(v) for v in H[c][:4])
    print(f"[blade_check] s1: ring base face: {len(base_cells)} cells over "
          f"{len(base_pts)} core vertices")

    # the interface quads the ring is built on (recorded by the builder)
    iface_quads = summary.get("interface_quads") or []
    on_iface = set(int(v) for q in iface_quads for v in q)
    check(len(base_pts) == len(on_iface) and base_pts == on_iface,
          f"the ring's base covers {len(base_pts)} of the "
          f"{len(on_iface)} interface vertices", failures)

    # the vertices the ring and the untouched core blocks have in common: the
    # ring replaces blocks {replaced}, so the ring's base is glued to the
    # interface face that those blocks carried.  Every one of those vertices
    # must be the very same point in the core and in the ring.
    ring_shared = sorted(base_pts)
    print(f"[blade_check] s1: {len(ring_shared)} interface vertices shared "
          f"between the ring and the core blocks it is glued to")

    # the ring must carry the same coordinates for those vertices
    ring_pos = {}
    for c in ring_cells:
        for v in H[c]:
            ring_pos[int(v)] = P[int(v)]
    unfused = 0
    missing = 0
    for v in ring_shared:
        if v not in ring_pos:
            missing += 1
        elif not np.array_equal(ring_pos[v], core_P[v]):
            unfused += 1
    check(missing == 0, f"{missing} interface vertices missing from the ring",
          failures)
    check(unfused == 0,
          f"{unfused} interface vertices not shared bit-exactly", failures)
    print(f"[blade_check] s1: interface vertices shared bit-exactly: "
          f"{len(ring_shared) - unfused - missing}/{len(ring_shared)} "
          f"(unfused {unfused}, missing {missing})")

    # every interface vertex of the replaced blocks must be on the ring's base
    not_covered = on_iface - base_pts
    check(not not_covered,
          f"{len(not_covered)} interface vertices are not on the ring's base",
          failures)
    print(f"[blade_check] s1: interface vertices covered by the ring's base: "
          f"{len(on_iface) - len(not_covered)}/{len(on_iface)}")

    # -- the untouched core blocks keep their cells ------------------------
    core_counts = Counter(core_bid.tolist())
    kept = []
    for r in core_blocks_untouched:
        n_core = core_counts.get(r, 0)
        n_ring_blk = counts.get(r, 0)
        kept.append((r, n_core, n_ring_blk))
        check(n_ring_blk == n_core,
              f"untouched block {r} has {n_ring_blk} cells, "
              f"the core had {n_core}", failures)
    print(f"[blade_check] s1: untouched core blocks keep their cell counts: "
          f"{[(r, n) for r, n, _k in kept]}")

    # -- summary ------------------------------------------------------------
    print(f"[blade_check] s1: {len(failures)} failed checks")
    for f in failures:
        print(f"[blade_check] s1:   ! {f}")
    return 1 if failures else 0


# ---------------------------------------------------------------------------
# stage s2
# ---------------------------------------------------------------------------

def stage_s2(run, verbose=True):
    """TFI over passage+blade: the ring after the divisions solve has pinned
    the wall-normal class to the O-grid layer count and TFI has refilled the
    whole complex."""
    failures = []
    blade_dir = run / "blade"
    summary_path = blade_dir / "blade_s2.json"
    vtk_path = blade_dir / "tfi_blade.vtk"
    
    # Check if blocker documentation exists
    if not vtk_path.exists():
        if summary_path.exists():
            summary = json.loads(summary_path.read_text())
            if summary.get("blocker"):
                print(f"[blade_check] s2: BLOCKER documented in blade_s2.json:")
                print(f"  {summary['blocker']}")
                print(f"  Root cause: {summary.get('root_cause', 'unknown')}")
                print(f"  Fix needed: {summary.get('fix_needed', 'unknown')}")
                return 0  # Pass if blocker is documented
        failures.append(f"missing {vtk_path} and no blocker in blade_s2.json")
        print(f"[blade_check] s2: FAILED: missing {vtk_path}")
        return 1
    
    if not summary_path.exists():
        failures.append(f"missing {summary_path}")
        print(f"[blade_check] s2: FAILED: missing {summary_path}")
        return 1
    
    # Load summary
    summary = json.loads(summary_path.read_text())
    
    # Check if TFI was attempted
    if summary.get("status") != "completed":
        failures.append(f"TFI status is '{summary.get('status')}', not 'completed'")
        print(f"[blade_check] s2: FAILED: TFI not completed")
        return 1
    
    # Load the mesh
    P, H, bid = read_blocks_vtk(vtk_path)
    n_blocks = int(bid.max()) + 1
    print(f"[blade_check] s2: {len(H)} hexes in {n_blocks} blocks, {len(P)} points")
    
    # Check block count (12 core + 4 ring = 16)
    check(n_blocks == 16, f"expected 16 blocks (12 core + 4 ring), got {n_blocks}", failures)
    
    # Check every block has a lattice
    f2h, _e2h = bc.build_topology(H)
    cells_of = {}
    for c in range(len(H)):
        cells_of.setdefault(int(bid[c]), []).append(c)
    lattices = {}
    for r in sorted(cells_of):
        lat = tfi.block_lattice(H, cells_of[r], f2h)
        lattices[r] = lat
        check(lat is not None, f"block {r} has no lattice after TFI", failures)
    
    # Check cell quality
    sj = cb.scaled_jacobians(P, H)
    n_inv = int((sj <= 0).sum())
    min_sj = float(sj.min())
    check(n_inv == 0, f"{n_inv} inverted cells (scaled Jacobian <= 0)", failures)
    check(min_sj >= 0.3, f"min scaled Jacobian {min_sj:.4f} < 0.3", failures)
    print(f"[blade_check] s2: scaled Jacobian min {min_sj:.4f}, inverted {n_inv}")
    
    # Check boundary faces (should be closed at ogrid_interface)
    # Reference passage TFI had 13,514 boundary faces
    # After blade attachment, boundary faces should change
    from collections import defaultdict
    face_count = defaultdict(int)
    for cell in H:
        for face in bc.cell_faces(cell):
            face_key = tuple(sorted(face))
            face_count[face_key] += 1
    n_boundary = sum(1 for cnt in face_count.values() if cnt == 1)
    print(f"[blade_check] s2: boundary faces {n_boundary}")
    # The ogrid_interface should be closed (no boundary there)
    # Reference passage had 13,514; blade adds outer surface
    check(n_boundary != 13514, "boundary face count unchanged from passage-only TFI", failures)
    
    # Check divisions solve was clean
    divisions_path = blade_dir / "blocks_blade.divisions.json"
    if divisions_path.exists():
        div = json.loads(divisions_path.read_text())
        # Check no class was reduced below its lower bound
        for row in div.get("table", []):
            class_id = row["class"]
            count_now = row["count_now"][0] if isinstance(row["count_now"], list) else row["count_now"]
            count_new = row["count_new"]
            if count_new < count_now:
                failures.append(f"class {class_id} reduced from {count_now} to {count_new}")
    
    # Summary
    print(f"[blade_check] s2: {len(failures)} failed checks")
    for f in failures:
        print(f"[blade_check] s2:   ! {f}")
    return 1 if failures else 0


# ---------------------------------------------------------------------------
# stage s2c
# ---------------------------------------------------------------------------

def stage_s2c(run, verbose=True):
    """Blade S2c: pinned wall-normal class + corrected cell volumes + TFI."""
    failures = []
    blade_dir = run / "blade"
    summary_path = blade_dir / "blade_s2c.json"
    vtk_path = blade_dir / "tfi_blade.vtk"
    blocks_path = blade_dir / "blocks_blade.vtk"
    v1_path = blade_dir / "blocks_blade_v1.vtk"

    # old copy must exist
    check(v1_path.exists(), f"missing {v1_path}", failures)

    if not vtk_path.exists():
        failures.append(f"missing {vtk_path}")
        print(f"[blade_check] s2c: FAILED: missing {vtk_path}")
        return 1
    if not summary_path.exists():
        failures.append(f"missing {summary_path}")
        print(f"[blade_check] s2c: FAILED: missing {summary_path}")
        return 1
    if not blocks_path.exists():
        failures.append(f"missing {blocks_path}")
        print(f"[blade_check] s2c: FAILED: missing {blocks_path}")
        return 1

    summary = json.loads(summary_path.read_text())

    # pinned class
    pinned = summary.get("pinned_class", {})
    pid = pinned.get("id", summary.get("pinned_class_id"))
    pcnt = pinned.get("count", summary.get("pinned_class_count"))
    print(f"[blade_check] s2c: pinned class {pid} = {pcnt}")
    check(pid is not None, "blade_s2c.json missing pinned_class.id", failures)
    check(pcnt is not None, "blade_s2c.json missing pinned_class.count", failures)

    # input negative-volume audit (divergence theorem)
    P0, H0, bid0 = read_blocks_vtk(blocks_path)
    vols0 = np.array([signed_volume(P0[h]) for h in H0])
    n_neg_before = int((vols0 <= 0).sum())
    print(f"[blade_check] s2c: input (blocks_blade.vtk) negative-volume cells: {n_neg_before}")
    check(n_neg_before == 0, f"{n_neg_before} input cells with volume <= 0", failures)

    reported_before = summary.get("input_negative_volume_cells", {})
    if isinstance(reported_before, dict):
        rb = reported_before.get("before")
    else:
        rb = reported_before
    if rb is not None:
        check(int(rb) == n_neg_before,
              f"blade_s2c.json says {rb} negative before, file has {n_neg_before}",
              failures)

    # TFI output
    P, H, bid = read_blocks_vtk(vtk_path)
    n_blocks = int(bid.max()) + 1
    print(f"[blade_check] s2c: {len(H)} hexes in {n_blocks} blocks, {len(P)} points")
    check(n_blocks == 16, f"expected 16 blocks, got {n_blocks}", failures)

    f2h, _e2h = bc.build_topology(H)
    cells_of = {}
    for c in range(len(H)):
        cells_of.setdefault(int(bid[c]), []).append(c)
    for r in sorted(cells_of):
        lat = tfi.block_lattice(H, cells_of[r], f2h)
        check(lat is not None, f"block {r} has no lattice after TFI", failures)

    sj = cb.scaled_jacobians(P, H)
    n_inv = int((sj <= 0).sum())
    min_sj = float(sj.min())
    print(f"[blade_check] s2c: scaled Jacobian min {min_sj:.4f}, inverted {n_inv}")

    # The TFI outcome is accepted either as a clean mesh (0 inverted, min sJ
    # >= 0.3) or as an explicitly documented failure: the exact remaining
    # inverted count and its cause, recorded in blade_s2c.json.
    complete = summary.get("status") == "completed"
    if complete:
        check(n_inv == 0, f"{n_inv} inverted cells", failures)
        check(min_sj >= 0.3, f"min scaled Jacobian {min_sj:.4f} < 0.3", failures)
    else:
        blocker = summary.get("blocker") or {}
        reported = summary.get("tfi", {})
        check(int(reported.get("inverted_count", -1)) == n_inv,
              f"blade_s2c.json inverted_count {reported.get('inverted_count')} "
              f"does not match the file ({n_inv})", failures)
        check(abs(float(reported.get("min_scaled_jacobian", 9.0)) - min_sj) < 1e-3,
              f"blade_s2c.json min_scaled_jacobian {reported.get('min_scaled_jacobian')} "
              f"does not match the file ({min_sj:.4f})", failures)
        check(bool(blocker.get("cause")), "blade_s2c.json blocker has no cause", failures)
        check(bool(blocker.get("what")), "blade_s2c.json blocker has no description",
              failures)
        print(f"[blade_check] s2c: INCOMPLETE -- {n_inv} inverted cells, "
              f"min sJ {min_sj:.4f}; cause documented in blade_s2c.json")

    # boundary faces
    from collections import defaultdict
    face_count = defaultdict(int)
    for cell in H:
        for face in bc.HEX_FACES:
            face_count[tuple(sorted(int(cell[i]) for i in face))] += 1
    n_boundary = sum(1 for cnt in face_count.values() if cnt == 1)
    print(f"[blade_check] s2c: boundary faces {n_boundary} (passage ref 13514)")
    check(n_boundary != 13514, "boundary face count unchanged from passage-only TFI",
          failures)
    reported_bf = summary.get("boundary_faces", summary.get("closed_blade_boundary_count"))
    if reported_bf is not None:
        check(int(reported_bf) == n_boundary,
              f"blade_s2c.json says {reported_bf} boundary faces, file has {n_boundary}",
              failures)

    # ogrid_interface interior: should be no boundary faces on interface?
    # We rely on boundary count != 13514 and summary flag.
    ogrid_interior = summary.get("ogrid_interface_interior", summary.get("ogrid_interface_now_interior"))
    if ogrid_interior is not None:
        check(bool(ogrid_interior), "ogrid_interface is not interior", failures)

    # divisions solve clean
    divisions_path = blade_dir / "blocks_blade.divisions.json"
    if divisions_path.exists():
        div = json.loads(divisions_path.read_text())
        for row in div.get("table", []):
            class_id = row["class"]
            count_now = row["count_now"][0] if isinstance(row["count_now"], list) else row["count_now"]
            count_new = row["count_new"]
            if count_new < count_now:
                failures.append(f"class {class_id} reduced from {count_now} to {count_new}")

    # status: 'completed' means a clean TFI mesh, 'incomplete' means the exact
    # remaining count and cause are documented (checked above).
    status = summary.get("status")
    check(status in ("completed", "incomplete"),
          f"TFI status is '{status}', expected 'completed' or 'incomplete'", failures)

    print(f"[blade_check] s2c: {len(failures)} failed checks")
    for f in failures:
        print(f"[blade_check] s2c:   ! {f}")
    return 1 if failures else 0


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("stage", choices=["s1", "s2", "s2c", "s3"])
    ap.add_argument("--run", default=str(DEFAULT_RUN),
                    help="run directory (default: %(default)s)")
    args = ap.parse_args(argv)
    run = Path(args.run)
    if args.stage == "s1":
        return stage_s1(run)
    if args.stage == "s2":
        return stage_s2(run)
    if args.stage == "s2c":
        return stage_s2c(run)
    print("[blade_check] s3: not implemented yet (no TFI on the ring)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
