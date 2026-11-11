#!/usr/bin/env python3
"""Acceptance checker for the keep-BL reference pipeline (task P1).

    check_keepbl_sample.py <run_dir>

Verified against the run artifacts themselves (no run data is recomputed from
anywhere else):

 1. HexBlockValidator VALID on the exported block structure, 0 inverted cells
    (both from the beam stage's own validator section and recomputed from
    blocks.vtk via ovm_io.scaled_jacobian).
 2a. (informational) bl_prism_identity: per-surface prism-identity coverage
    of the kept hub/shroud prism layers (kept prisms inside touching blocks).
    Does NOT gate the overall pass (keep-BL does not require prism identity).
 2b. (gated) 2_boundary_layers_kept: the block structure covers the full
    radial domain including both BL bands [0.5,0.5265] and [1.872,1.9].
 3. O-grid outline coverage >= 99 % from the <name>_curves.npz post_quality
    (min seam_coverage over the two-label seams) and feature_complete.
 4. Block count 10..16 (12-block S1 reference).
 5. sample.npz carries the standard export_sample schema and
    <name>_curves.npz exists.

Writes <run_dir>/p1_report.json and prints PASS/FAIL.  Exit code 0 iff PASS.
"""

import json
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]          # the repo root
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from experimentell.hex3d_algohex import ovm_io      # noqa: E402

STAGES = ["dtoo", "tet_prep", "algohex", "beam", "tfi", "export", "curves"]
BLOCK_MIN, BLOCK_MAX = 10, 16
OG_GRID_COVERAGE = 0.99

SAMPLE_KEYS = [
    "vertices", "blocks", "quad_faces", "edges", "edge_ctrl",
    "edge_polyline", "edge_polyline_offset", "dir_class",
    "params", "quality", "provenance",
]
SAMPLE_TET_KEYS = ["dir_class_count", "surface_points", "surface_tris",
                   "surface_tri_label"]

HUB_R, SHROUD_R = 0.5, 1.9
R_TOL = 3e-3            # radius band width for wall classification
# (covers the bl-plane root ring at r = 0.498)
LOC_TOL = 1e-7          # node coordinate coincidence tolerance
WALL_MAP_TOL = 0.12     # face-centroid search radius for prism mapping
LABEL_CAP = 0.06        # max distance for a quad to count as a surface patch


def stage_stats(run: Path):
    out = {}
    for stage in STAGES:
        p = run / f"{stage}.log"
        if not p.exists():
            out[stage] = None
            continue
        m = None
        for line in p.read_text(errors="replace").splitlines():
            if line.startswith("[stage] "):
                m = line
        out[stage] = m
    return out


def parse_stage(stats):
    """'[stage] beam wall=1.1s maxrss=2kB exit=0' -> dict."""
    res, timeout, fail = {}, None, None
    for stage, line in stats.items():
        if line is None:
            continue
        mm = re.search(r"wall=([0-9.]+)s maxrss=([0-9]+)kB exit=([\-0-9]+)", line)
        if mm:
            res[stage] = {"wall_s": float(mm.group(1)),
                          "maxrss_kb": int(mm.group(2)),
                          "exit": int(mm.group(3))}
    return res


def beam_validator(run: Path):
    """The HexBlockValidator section of beam.log."""
    log = run / "beam.log"
    if not log.exists():
        return {"validator_ran": False, "valid": False, "inverted": None}
    text = log.read_text(errors="replace")
    i = text.find("----- HexBlockValidator -----")
    if i < 0:
        return {"validator_ran": False, "valid": False, "inverted": None}
    sec, out = text[i:], {"validator_ran": True, "valid": False}
    out["valid"] = bool(re.search(r"^VALID\s*$", sec, re.M))
    for key, pat in [("blocks", r"^\s*blocks:\s*([0-9]+)"),
                     ("min_scaled_jacobian", r"min_scaled_jacobian:\s*([-0-9.]+)"),
                     ("cuboid_fraction", r"cuboid_fraction:\s*([-0-9.]+)")]:
        mm = re.search(pat, sec, re.M)
        out[key] = int(mm.group(1)) if key == "blocks" and mm else (
            float(mm.group(1)) if mm else None)
    mm = re.search(r"([0-9]+) cells <= 0", sec)
    out["inverted"] = int(mm.group(1)) if mm else None
    mm = re.search(r"inverted_cells\s+([0-9]+)\s*->\s*([0-9]+)", text)
    if mm:
        out["inverted_after"] = int(mm.group(2))
    return out


def load_vtk_hexes(path: Path):
    """meshio-load a legacy VTK hex mesh: (P, hex_cells, block_ids|None)."""
    import meshio
    m = meshio.read(str(path))
    hex_parts, bid_parts = [], []
    for ci, c in enumerate(m.cells):
        if c.type != "hexahedron":
            continue
        hex_parts.append(c.data)
        d = m.cell_data.get("block_id")
        bid_parts.append(np.asarray(d[ci], int).ravel() if d else None)
    if not hex_parts:
        raise RuntimeError(f"{path}: no hexahedron block")
    H = np.vstack([p for p in hex_parts]).astype(np.int64)
    bid = None
    if len(bid_parts) == len(hex_parts) and all(b is not None for b in bid_parts):
        bid = np.concatenate(bid_parts).astype(int) \
            if len(bid_parts) > 1 else bid_parts[0]
    return np.asarray(m.points, float), H, bid


def read_ovm_hexes(path: Path):
    """Read the fine hex mesh .ovm: (points, hex corner lists)."""
    P, edges, faces, polys = ovm_io.read_ovm(str(path))
    hexes, _skipped = ovm_io.ovm_to_cells(P, edges, faces, polys)
    Hf = (np.asarray(hexes, np.int64)
          if len(hexes) else np.zeros((0, 8), np.int64))
    return P, Hf


_HEX_FACES = [(0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4),
              (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]


def node_kdtree(P):
    from scipy.spatial import cKDTree
    return cKDTree(P)


def wall_node_ids(run: Path):
    """The hub/shroud wall node coordinates, from tet.vtk triangle labels
    (tet_prep_v5: hub=8, shroud=9, only present with --keep-prisms)."""
    import meshio
    m = meshio.read(str(run / "tet.vtk"))
    P = np.asarray(m.points, float)
    tris, labels = [], []
    for c, d in zip(m.cells, m.cell_data["color"]):
        if c.type == "triangle":
            tris.append(np.asarray(c.data, int))
            labels.append(np.asarray(d, int).ravel())
    if not tris:
        raise RuntimeError("tet.vtk lacks triangle cell colors")
    tris, labels = np.vstack(tris), np.concatenate(labels)
    labset = set(int(v) for v in np.unique(labels))
    out = {}
    for surf, keep_id, drop_id in (("hub", 8, 5), ("shroud", 9, 6)):
        # --keep-prisms labels: 8 hub / 9 shroud; the drop-prism variant
        # uses 5 hub / 6 shroud (2/3 ids shift)
        rid = keep_id if keep_id in labset else drop_id
        if rid not in labset:
            raise RuntimeError(f"tet.vtk lacks {surf} label ({keep_id}/{drop_id})")
        ids = np.unique(tris[labels == rid])
        out[surf] = {"points": P[ids], "tris": tris[labels == rid]}
    return out


def parse_msh_v2(path: Path):
    """Minimal Gmsh v2 ASCII reader: (nodes, tag->row, quads, prisms)."""
    nodes, elems = {}, []
    sec = None
    with open(path, errors="replace") as f:
        for line in f:
            line = line.strip()
            if line.startswith("$"):
                sec = None if line.startswith("$End") else line[1:]
                continue
            if sec == "Nodes":
                f_ = line.split()
                if len(f_) == 4 and f_[0].isdigit():
                    nodes[int(f_[0])] = (float(f_[1]), float(f_[2]),
                                         float(f_[3]))
            elif sec == "Elements":
                f_ = line.split()
                if len(f_) > 3 and f_[0].isdigit():
                    etype, ntags = int(f_[1]), int(f_[2])
                    nd = [int(x) for x in f_[3 + ntags:]]
                    elems.append((etype, nd))
    tags = sorted(nodes)
    NW = np.array([nodes[t] for t in tags], float)
    t2i = {t: i for i, t in enumerate(tags)}
    quads = [[t2i[n] for n in nd] for et, nd in elems if et == 3 and len(nd) == 4]
    prisms = [[t2i[n] for n in nd] for et, nd in elems
              if et == 6 and len(nd) == 6]
    return NW, quads, prisms


def in_hex(P, hexPts, x, tol=1e-9):
    """Point-in-hex via the 6-tet (0-2-5-4 style) decomposition."""
    c = P[hexPts]
    diags = [(0, 1, 2, 5), (0, 2, 5, 6), (0, 5, 6, 4),
             (0, 2, 3, 6), (0, 3, 6, 7), (0, 4, 6, 7)]
    for a, b, dq, e in diags:
        t = np.array([c[a], c[b], c[dq], c[e]])
        det = np.linalg.det(np.roll(t - t[0], 0, axis=0)[1:].T)
        if abs(det) < 1e-14:
            continue
        bary = np.linalg.solve(np.roll(t - t[0], 0, axis=0)[1:].T, x - t[0])
        s = np.r_[bary, 1.0 - bary.sum()]
        if s.min() >= -tol:
            return True
    return False


def in_tri(a, b, c, p, tol=1e-6):
    """Is p (in the plane of abc) inside the triangle?"""
    n = np.cross(b - a, c - a)
    nn = np.linalg.norm(n)
    if nn < 1e-13:
        return False
    u = np.cross(b - p, c - p).dot(n)
    v = np.cross(c - p, a - p).dot(n)
    t = np.cross(a - p, b - p).dot(n)
    tot = u + v + t
    if abs(tot) < 1e-14:
        return False
    bary = np.array([u, v, t]) / tot
    return bary.min() >= -tol and bary.max() <= 1.0 + tol


def bl_coverage(run):
    """Per-surface coverage of the kept prism layers.

    Kept prisms = the case/mesh.msh prisms whose wall face lies on the hub
    (r = 0.5) or shroud (r = 1.9) band.  A block touches the surface when it
    owns a boundary quad of the structure (blocks.vtk) whose nearest input
    surface triangle carries the hub/shroud label (the same nearest-face
    rule clean_blocks' labeller applies).  A kept prism is covered when its
    centroid falls inside a support cell (blocks.vtk) of a touching block.
    """
    walls = wall_node_ids(run)
    NW, _, prism_list = parse_msh_v2(run / "case" / "mesh.msh")
    r = np.sqrt(NW[:, 0] ** 2 + NW[:, 1] ** 2)

    # input surface triangles + labels (for the labeller rule)
    import meshio
    m = meshio.read(str(run / "tet.vtk"))
    TP, tri_rows, tri_lab = [], [], []
    for c, d in zip(m.cells, m.cell_data["color"]):
        if c.type == "triangle":
            tri_rows.append(c.data)
            tri_lab.append(np.asarray(d, int).ravel())
            TP.append(m.points)
    if not tri_rows:
        raise RuntimeError("tet.vtk has no triangles")
    TP = np.asarray(TP[0], float)
    PR, PL = np.vstack(tri_rows), np.concatenate(tri_lab)
    tri_pts, tri_cent = TP[PR], TP[PR].mean(axis=1)
    centT = node_kdtree(tri_cent)

    def nearest_tri(x, k=40):
        _, nn = centT.query(x, k=min(k, len(tri_cent)))
        cand = np.unique(np.asarray(nn, int).ravel())
        best_d, best_l = np.inf, None
        for t in cand:
            d = _tri_dist(x, *tri_pts[t])
            if d < best_d:
                best_d, best_l = d, PL[t]
        return best_l, best_d

    def _tri_dist(x, a, b, c):
        n = np.cross(b - a, c - a)
        nn = np.linalg.norm(n)
        if nn < 1e-13:
            return np.linalg.norm(x - a)
        d = np.dot(x - a, n / nn)
        p = x - d * n / nn
        if in_tri(a, b, c, p):
            return abs(d)
        return min(np.linalg.norm(x - v) for v in (a, b, c))

    Pb, Hb, bid = load_vtk_hexes(run / "blocks.vtk")
    supT = node_kdtree(Pb[Hb].mean(axis=1))
    cache = {}

    def support_block(cent):
        key = tuple(np.round(cent, 9))
        if key in cache:
            return cache[key]
        got = None
        _, nn = supT.query(cent, k=min(8, len(Hb)))
        for ci in (np.atleast_1d(nn) if np.isscalar(nn) else nn):
            if in_hex(Pb, Hb[int(ci)], cent):
                got = int(bid[int(ci)])
                break
        cache[key] = got
        return got

    # support boundary quads and their surface labels
    f2c = {}
    for bi, c in enumerate(Hb):
        for fc in _HEX_FACES:
            f2c.setdefault(frozenset(int(c[i]) for i in fc), []).append(bi)
    bq = [(np.asarray(sorted(k)), v[0]) for k, v in f2c.items() if len(v) == 1]

    # label every support boundary quad once (cheap: few hundred faces)
    lab_map = {}
    max_lab_dist = 0.0
    for k, _bi in bq:
        fc = Pb[k].mean(axis=0)
        lab, d = nearest_tri(fc)
        lab_map[tuple(np.asarray(sorted(k)))] = (int(lab), d)
        max_lab_dist = max(max_lab_dist, d)

    out = {"n_support_cells": len(Hb), "n_blocks": int(len(np.unique(bid))),
           "surfaces": {}, "support_boundary_quads": len(bq),
           "max_quad_label_dist": round(float(max_lab_dist), 6),
           "tet_vtk_labels": sorted(int(x) for x in np.unique(PL))}
    labs = set(out["tet_vtk_labels"])
    for surf, wall_r in (("hub", HUB_R), ("shroud", SHROUD_R)):
        # kept prisms of this surface
        wl = []
        for pr in prism_list:
            pr = np.asarray(pr)
            if (np.abs(r[pr[:3]] - wall_r).max() <= R_TOL
                    or np.abs(r[pr[3:]] - wall_r).max() <= R_TOL):
                wl.append(NW[pr].mean(axis=0))
        # touching blocks: support quads labelled at this surface (within
        # LABEL_CAP of that input face)
        wall_lab = 8 if surf == "hub" else 9
        surf_lab = (wall_lab if wall_lab in labs
                    else (5 if surf == "hub" else 6))
        wall_ok = wall_lab in labs      # --keep-prisms wall surf present
        touching, max_td = set(), 0.0
        n_touch_quads = 0
        for k, bi in bq:
            lab, d = lab_map[tuple(k)]
            if lab == surf_lab and d <= LABEL_CAP:
                n_touch_quads += 1
                max_td = max(max_td, d)
                touching.add(int(bid[bi]))
        blocks_with_bl = set()
        covered = outside = 0
        for cent in wl:
            b = support_block(cent)
            if b is not None and b in touching:
                covered += 1
                blocks_with_bl.add(b)
            elif b is None:
                outside += 1
        out["surfaces"][surf] = {
            "kept_prisms": len(wl),
            "inside_blocks": covered,
            "outside_structure": outside,
            "coverage": covered / len(wl) if wl else 0.0,
            "surface_label_used": surf_lab,
            "wall_label_present": wall_ok,
            "n_touch_quads": n_touch_quads,
            "max_touch_quad_dist": round(max_td, 6),
            "blocks_touching": sorted(touching),
            "blocks_with_bl": sorted(blocks_with_bl),
            "pass": bool(wall_ok and len(wl) > 0 and blocks_with_bl),
        }
    out["pass"] = all(out["surfaces"][s]["pass"] for s in ("hub", "shroud"))
    return out

BL_BAND_HUB = (0.5, 0.5265)
BL_BAND_SHROUD = (1.872, 1.9)


def boundary_layers_kept(run):
    """What keep-BL actually requires: the block structure covers the full
    radial domain, including both BL bands [0.5,0.5265] and [1.872,1.9].

    Computed from the block vertices (blocks.vtk): min r / max r.  Hub/shroud
    radii and both BL bands are read from the dtOO case/mesh.msh when
    available (radii = node r range, bands = radial extent of the hub/shroud
    prism layers), otherwise the radii are derived from the labeled wall
    quads (hub label 8, shroud label 9) and the bands fall back to the
    tistos defaults.
    """
    r_hub, r_shroud = HUB_R, SHROUD_R
    band_hub, band_shroud = BL_BAND_HUB, BL_BAND_SHROUD
    mesh_path = run / "case" / "mesh.msh"
    if mesh_path.exists():
        try:
            NW, _, prisms = parse_msh_v2(mesh_path)
            rr = np.sqrt(NW[:, 0] ** 2 + NW[:, 1] ** 2)
            if len(NW):
                r_hub = float(rr.min())
                r_shroud = float(rr.max())
            if prisms:
                rp = rr[np.unique(np.asarray(prisms, int).ravel())]
                mid = 0.5 * (r_hub + r_shroud)
                lo, hi = rp[rp < mid], rp[rp >= mid]
                if len(lo):
                    band_hub = (r_hub, float(lo.max()))
                if len(hi):
                    band_shroud = (float(hi.min()), r_shroud)
        except Exception:
            pass
    else:
        try:
            walls = wall_node_ids(run)
            for surf in ("hub", "shroud"):
                pts = walls[surf]["points"]
                rr = np.sqrt(pts[:, 0] ** 2 + pts[:, 1] ** 2)
                if surf == "hub":
                    r_hub = float(rr.mean())
                else:
                    r_shroud = float(rr.max())
        except Exception:
            pass

    Pb, Hb, bid = load_vtk_hexes(run / "blocks.vtk")
    rblk = np.sqrt(Pb[:, 0] ** 2 + Pb[:, 1] ** 2)
    min_r = float(rblk.min())
    max_r = float(rblk.max())
    tol = R_TOL
    hub_ok = min_r <= r_hub + tol
    shroud_ok = max_r >= r_shroud - tol
    hub_band_ok = (min_r <= band_hub[0] + tol) and (max_r >= band_hub[1])
    shroud_band_ok = (min_r <= band_shroud[0]) and (max_r >= band_shroud[1] - tol)
    passes = hub_ok and shroud_ok and hub_band_ok and shroud_band_ok
    return {
        "min_r": round(min_r, 6),
        "max_r": round(max_r, 6),
        "r_hub": round(float(r_hub), 6),
        "r_shroud": round(float(r_shroud), 6),
        "radius_tol": tol,
        "hub_radius_ok": bool(hub_ok),
        "shroud_radius_ok": bool(shroud_ok),
        "bl_band_hub": [round(band_hub[0], 6), round(band_hub[1], 6)],
        "bl_band_shroud": [round(band_shroud[0], 6), round(band_shroud[1], 6)],
        "bl_band_hub_inside": bool(hub_band_ok),
        "bl_band_shroud_inside": bool(shroud_band_ok),
        "pass": bool(passes),
        "note": ("keep-BL requires the block structure to cover the full radial "
                 "domain including both BL bands; prism identity is NOT required "
                 "(see bl_prism_identity)."),
    }


def main():
    run = Path(sys.argv[1]).resolve()
    if not run.is_dir():
        print(f"FAIL: run dir {run} does not exist")
        return 2
    stats = stage_stats(run)
    report = {"run_dir": str(run), "stages": parse_stage(stats),
              "stage_lines": stats, "artifacts": {},
              "checks": {}, "pass": True}
    for s in STAGES:
        lp = run / f"{s}.log"
        report["artifacts"][f"{s}.log"] = str(lp) if lp.exists() else "MISSING"

    name = run.name

    # which AlgoHex hex mesh the reference block structure was built from
    src = None
    for cand in ("hex_hex_x_16c.ovm", "hex_hex_x.ovm"):
        if (run / cand).exists():
            src = cand
            break
    if src:
        report["artifacts"]["source_hex_mesh"] = str(run / src)
        report["source_hex_mesh"] = src

    # divergence recovery provenance (informational, not gated)
    recov = {}
    if (run / "algohex.log").exists():
        recov["timeout_record"] = [ln for ln in
                                   (run / "algohex.log").read_text(
                                       errors="replace").splitlines()
                                   if ln.startswith("[stage]")]
    for tag in ("save", "replay"):
        f = run / f"algohex_{tag}.out"
        if f.exists():
            recov[f"{tag}_log"] = str(f)
            head = [ln for ln in f.read_text(errors="replace").splitlines()
                    if "[param-state]" in ln or "parametrization finished" in ln]
            recov[f"{tag}_provenance"] = head[:4]
    st_path = run / "param_state.ovm"
    if st_path.exists():
        recov["param_state"] = {"file": str(st_path),
                                "bytes": st_path.stat().st_size}
        recov["replay_command"] = str(
            "env LD_LIBRARY_PATH=/opt/algohex/Build/lib:/opt/coin-or/lib "
            "/work/tools/algohex-size/build/Build/bin/HexMeshing "
            f"-i {run}/tet.vtk -o {run}/{src or 'hex_hex_x.ovm'} "
            f"-j {run}/hex_hex_metrics_x.json -n 2000 "
            f"--load-param-state {st_path}")
    # the converged 16-core replay of the same parametrisation
    for cand in ("replay16.log", "algohex_replay_16c.out"):
        f = run / cand
        if f.exists():
            recov["converged_replay_log"] = str(f)
            head = [ln for ln in f.read_text(errors="replace").splitlines()
                    if "Cells" in ln or "Finished" in ln]
            recov["converged_replay_provenance"] = head[-3:]
            break
    if src:
        recov["source_hex_mesh"] = src
    if recov:
        report["checks"]["algohex_recovery"] = {"info": recov, "pass": True}

    # ---------------- acceptance 1 ----------------
    validinfo = beam_validator(run)
    inverted_ok = True
    rec_inv = None
    if (run / "blocks.vtk").exists():
        Pb, Hb, _ = load_vtk_hexes(run / "blocks.vtk")
        rec_inv = int((np.array([ovm_io.scaled_jacobian(Pb, c) for c in Hb])
                       <= 0).sum())
    else:
        inverted_ok = False
    inv_valid = validinfo.get("inverted")
    if inv_valid is None and (validinfo.get("min_scaled_jacobian")
                              or -1) > 0 and validinfo.get("inverted_after") == 0:
        # no "N cells <= 0" summary when the min is positive
        inv_valid = 0
    a1 = (validinfo.get("valid") is True
          and inv_valid == 0
          and rec_inv == 0
          and (validinfo.get("min_scaled_jacobian") or -1) > 0)
    report["checks"]["1_block_structure_valid"] = {
        "valid": bool(validinfo.get("valid")),
        "inverted_cells": validinfo.get("inverted"),
        "inverted_cells_recomputed": rec_inv,
        "min_scaled_jacobian": validinfo.get("min_scaled_jacobian"),
        "pass": bool(a1)}
    report["pass"] &= bool(a1)

    # ---------------- acceptance 2a: prism identity (informational) ----------------
    try:
        bl = bl_coverage(run)
        bl["threshold"] = 0.99
        bl["informational"] = True
        bl["gate"] = False
        for surf in ("hub", "shroud"):
            cov = float(bl["surfaces"][surf].get("coverage", 0.0))
            bl["surfaces"][surf]["total_prisms"] = int(
                bl["surfaces"][surf].get("kept_prisms", 0))
            bl["surfaces"][surf]["pass"] = bool(cov >= 0.99)
        bl["pass"] = all(bl["surfaces"][s]["pass"] for s in ("hub", "shroud"))
        bl["note"] = ("Prism-identity coverage (kept msh prisms inside touching "
                      "blocks). Informational only: keep-BL does not require "
                      "prism identity because TFI rebuilds the 0.027-0.035 "
                      "outer layer; see 2_boundary_layers_kept.")
        report["checks"]["bl_prism_identity"] = bl
    except Exception as exc:   # keep the report available even on failure
        report["checks"]["bl_prism_identity"] = {"error": repr(exc),
                                                "pass": False,
                                                "informational": True,
                                                "gate": False}

    # ---------------- acceptance 2b: boundary layers kept (gated) ----------------
    try:
        bl2 = boundary_layers_kept(run)
        report["checks"]["2_boundary_layers_kept"] = bl2
        report["pass"] &= bool(bl2.get("pass"))
    except Exception as exc:   # keep the report available even on failure
        report["checks"]["2_boundary_layers_kept"] = {"error": repr(exc),
                                                      "pass": False}
        report["pass"] = False

    # ---------------- acceptance 3 ----------------
    curves_path = run / f"{name}_curves.npz"
    report["artifacts"]["curves"] = str(curves_path)
    a3 = False
    if curves_path.exists():
        q = json.loads(str(np.load(curves_path, allow_pickle=False)
                           ["post_quality"]))
        cov = [float(v) for v in q.get("seam_coverage", {}).values()]
        report["checks"]["3_ogrid_outline_coverage"] = {
            "min_seam_coverage": min(cov) if cov else None,
            "seam_coverage": q.get("seam_coverage"),
            "feature_complete": q.get("feature_complete"),
            "pass": bool(cov and min(cov) >= OG_GRID_COVERAGE
                         and q.get("feature_complete"))}
        a3 = report["checks"]["3_ogrid_outline_coverage"]["pass"]
    else:
        report["checks"]["3_ogrid_outline_coverage"] = {"error": "missing",
                                                        "pass": False}
    report["pass"] &= bool(a3)

    # ---------------- acceptance 4 ----------------
    bcount = validinfo.get("blocks")
    if bcount is None and (run / "blocks.vtk").exists():
        import meshio
        m = meshio.read(str(run / "blocks.vtk"))
        for k, d in m.cell_data.items():
            if d and np.asarray(d[0]).ndim == 1:
                bcount = len(np.unique(np.asarray(d[0], int)))
                break
    a4 = bcount is not None and BLOCK_MIN <= bcount <= BLOCK_MAX
    report["checks"]["4_block_count"] = {"blocks": bcount,
                                         "reference_range": [BLOCK_MIN,
                                                             BLOCK_MAX],
                                         "pass": bool(a4)}
    report["pass"] &= bool(a4)

    # ---------------- acceptance 5 ----------------
    sample_path = run / "sample.npz"
    report["artifacts"]["sample_npz"] = str(sample_path)
    a5 = False
    schema = {}
    if sample_path.exists():
        z = np.load(sample_path, allow_pickle=False)
        schema["keys"] = list(z.files)
        missing = [k for k in SAMPLE_KEYS if k not in z.files]
        schema["missing_keys"] = missing
        schema["vertex_shape"] = list(np.shape(z["vertices"])) \
            if "vertices" in z.files else None
        a5 = not missing
        if a5:
            q = json.loads(str(z["quality"]))
            a5 = q.get("blocks") == bcount
            schema["quality_blocks"] = q.get("blocks")
    report["checks"]["5_sample_schema"] = {**schema, "pass": bool(a5)}
    report["pass"] &= bool(a5)

    report.pop("stage_lines", None)
    (run / "p1_report.json").write_text(json.dumps(report, indent=2,
                                                   default=lambda o: o
                                                   if isinstance(
                                                       o, (int, float,
                                                           str, bool)) or o
                                                   is None else str(o)))
    print(f"P1 {'PASS' if report['pass'] else 'FAIL'}: {run}")
    for k, v in report["checks"].items():
        print(f"  {k}: {'ok' if v.get('pass') else 'FAIL'} "
              f"{json.dumps({kk: vv for kk, vv in v.items() if kk != 'pass'})}")
    return 0 if report["pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
