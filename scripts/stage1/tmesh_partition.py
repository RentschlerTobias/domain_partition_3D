#!/usr/bin/env python3
"""Ansatz T: periodic seam as WALL + master-slave junctions + T-mesh blocks + TFI.

Alternative to the wrap approach (partition_surface.partition with periodicity
on): streamlines END on the theta seam instead of integrating across it. The
cross FIELD stays periodic (seam weld), so the separatrix geometry mirrors mod
pitch by itself; only the block stage treats the seam as a wall:

  1. field/separatrix pipeline with set_periodic(True) + set_tile_periodic(False)
  2. seam symmetrization: junction sets of both seams are unioned mod pitch
     (DLR Sauer/Morsbach 2023 sec 2.7 master-slave: only the left seam is
     parametrized, the right seam is the exact +pitch translate). Mirrored
     junctions without an interior curve are intentional hanging T-nodes.
  3. block extraction tolerating T-junctions (tmesh_faces: a block needs
     exactly 4 REAL corners; flat ~180deg nodes are allowed on its sides)
  4. TFI (Coons) fill per block; seam sides use the per-edge canonical
     sampling, so the right-seam discretization is the exact translate of the
     left one.

Outputs to output/T1_9/hub_stage1/tmesh/: tmesh_blocks.png, tmesh_tfi.png,
tmesh_metrics.json.
"""

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, "/root/repos/domain_partition")
sys.path.insert(0, str(Path(__file__).resolve().parent))

import unwrap_surface as us                                    # noqa: E402
import clean_separatrix as cs                                  # noqa: E402
import partition_surface as ps                                 # noqa: E402
import tmesh_faces as tmf                                      # noqa: E402
from dp_adapter import build_dp_data                           # noqa: E402
from tools import FrameField, StreamlineGenerator_v2           # noqa: E402
from tools.singularity_detector import detect_singularities    # noqa: E402
from tools.streamline_merging import StreamlineMerging         # noqa: E402
from tools.streamline_intersection_splitter import (           # noqa: E402
    StreamlineIntersectionSplitter)

STL = "/root/repos/block_structured_meshing/T1_9_hub_raw.stl"
OUT = Path("/root/repos/block_structured_meshing/output/T1_9/hub_stage1/tmesh")


# --------------------------------------------------------------------------
# seam identification + master-slave symmetrization
# --------------------------------------------------------------------------

def _dist_to_chain(p, chain):
    return ps._project_to_polyline(np.asarray(p, float), chain)[2]


def _periodic_chains(mesh):
    """(left, right) seam node chains from mesh.periodic_pairs, sorted by t."""
    pp = mesh.periodic_pairs
    pp = pp.numpy() if hasattr(pp, "numpy") else np.asarray(pp)
    X = mesh.x[:, 0:2].numpy()
    A, B = X[pp[:, 0]], X[pp[:, 1]]
    if A[:, 0].mean() > B[:, 0].mean():
        A, B = B, A
    return A[np.argsort(A[:, 1])], B[np.argsort(B[:, 1])]


def _chain_segments(segs):
    """Reassemble split wall segments into one polyline. Segment endpoints of a
    split wall coincide exactly, so exact-key chaining is safe."""
    from collections import Counter
    segs = [np.asarray(s, float) for s in segs]

    def key(p):
        return (round(float(p[0]), 9), round(float(p[1]), 9))

    cnt = Counter()
    for s in segs:
        cnt[key(s[0])] += 1
        cnt[key(s[-1])] += 1
    ends = [k for k, v in cnt.items() if v == 1]
    if len(ends) != 2:
        raise RuntimeError(f"wall segments do not chain (ends={len(ends)})")
    cur = min(ends, key=lambda k: k[1])          # start at smaller t
    remaining = list(range(len(segs)))
    parts = []
    while remaining:
        nxt = None
        for i in remaining:
            if key(segs[i][0]) == cur:
                nxt = (i, False)
                break
            if key(segs[i][-1]) == cur:
                nxt = (i, True)
                break
        if nxt is None:
            raise RuntimeError("wall segment chain broken")
        i, rev = nxt
        s = segs[i][::-1] if rev else segs[i]
        parts.append(s if not parts else s[1:])
        cur = key(s[-1])
        remaining.remove(i)
    W = np.vstack(parts)
    keep = np.concatenate([[True],
                           np.linalg.norm(np.diff(W, axis=0), axis=1) > 1e-12])
    return W[keep]


def _internal_junctions(segs, wall_ends, tol=5e-3):
    """Segment endpoints that are not the wall's own ends = T-junctions."""
    J = []
    for s in segs:
        for p in (np.asarray(s[0], float), np.asarray(s[-1], float)):
            if min(np.linalg.norm(p - e) for e in wall_ends) < tol:
                continue
            if not any(np.linalg.norm(p - q) < 1e-6 for q in J):
                J.append(p)
    return J


def _split_polyline_at(W, cuts):
    """Split polyline W at the given points (assumed on W). Returns segments."""
    marks = []
    for p in cuts:
        seg, t, dist, proj = ps._project_to_polyline(p, W)
        marks.append((seg + t, seg, np.asarray(p, float)))
    marks.sort(key=lambda m: m[0])
    out, cur, ptr = [], [W[0]], 0
    for i in range(len(W) - 1):
        while ptr < len(marks) and marks[ptr][1] == i:
            p = marks[ptr][2]
            if np.linalg.norm(cur[-1] - p) > 1e-12:
                cur.append(p)
            if len(cur) >= 2:
                out.append(np.array(cur))
            cur = [p]
            ptr += 1
        if np.linalg.norm(cur[-1] - W[i + 1]) > 1e-12:
            cur.append(W[i + 1])
    if len(cur) >= 2:
        out.append(np.array(cur))
    return out


def collapse_seam_wedges(mesh, gap=0.12, wall_tol=0.02):
    """Collapse 3-sided seam wedges (Xiao-style simplification at the wall).

    A singularity close to the theta seam sends two adjacent arms into the
    seam; with the seam as a WALL both dock a few hundredths apart and the
    region between them is a 3-corner triangle (in the wrap approach these
    arms integrate across the seam instead, so the wedge only exists here).
    Fix: drop the SHORTER arm of any same-singularity pair whose seam docks
    are adjacent on the wall (no other junction between) and closer than
    `gap`, and re-join the two wall segments at the freed junction."""
    A, B = _periodic_chains(mesh)
    n_b = len(mesh.streamlines) - len(mesh.separatrices)
    bnd = [np.asarray(s, float) for s in mesh.streamlines[:n_b]]
    seps = [np.asarray(s, float) for s in mesh.streamlines[n_b:]]
    dicts = list(mesh.separatrices)

    marks = []                       # (side, wall position, sep idx, end)
    for i, s in enumerate(seps):
        for side, C in (("L", A), ("R", B)):
            seg, t, dist, proj = ps._project_to_polyline(s[-1], C)
            if dist < wall_tol:
                marks.append((side, seg + t, i, s[-1]))
    from collections import defaultdict
    bysing = defaultdict(list)
    for side, pos, i, e in marks:
        oc = dicts[i].get("singularity_coords")
        if oc is None:
            continue
        oc = np.asarray(oc, float)
        bysing[(side, round(float(oc[0]), 6), round(float(oc[1]), 6))].append(
            (pos, i, e))
    all_pos = defaultdict(list)
    for side, pos, i, e in marks:
        all_pos[side].append(pos)

    drop, freed = set(), []
    for (side, _ocx, _ocy), lst in bysing.items():
        if len(lst) < 2:
            continue
        lst.sort(key=lambda x: x[0])
        for (p1, i1, e1), (p2, i2, e2) in zip(lst, lst[1:]):
            if np.linalg.norm(e2 - e1) > gap:
                continue
            if any(p1 + 1e-9 < q < p2 - 1e-9 for q in all_pos[side]):
                continue
            j = i1 if _arclen(seps[i1]) <= _arclen(seps[i2]) else i2
            drop.add(j)
            freed.append(e1 if j == i1 else e2)

    if not drop:
        return
    # re-join the two wall segments meeting at each freed junction
    for e in freed:
        hit = [k for k, b in enumerate(bnd)
               if min(np.linalg.norm(b[0] - e), np.linalg.norm(b[-1] - e))
               < 1e-9]
        if len(hit) != 2:
            continue
        k1, k2 = hit
        a, b = bnd[k1], bnd[k2]
        if np.linalg.norm(a[-1] - e) > 1e-9:
            a = a[::-1]
        if np.linalg.norm(b[0] - e) > 1e-9:
            b = b[::-1]
        bnd[k1] = np.vstack([a, b[1:]])
        bnd.pop(k2)
    mesh.streamlines = bnd + [s for k, s in enumerate(seps) if k not in drop]
    mesh.separatrices = [d for k, d in enumerate(dicts) if k not in drop]
    print(f"[wedge] collapsed {len(drop)} same-singularity seam wedge arm(s), "
          f"re-joined wall at {len(freed)} junction(s)")


def symmetrize_seam_junctions(mesh, tol_cls=0.012, tol_match=0.012):
    """Master-slave seam conformity (DLR sec 2.7). After the snap pass both
    seam walls are split at the T-junctions of the curves that ended there.
    This unions the junction sets mod pitch: the LEFT wall is the master; every
    junction (own or mirrored from the right) is projected onto it, the right
    wall is rebuilt as the EXACT +pitch translate (copy + shift, never
    re-projected), and the endpoints of the docking interior curves are moved
    onto the exact junction coordinates. Mirrored junctions without an interior
    curve stay hanging (T-nodes)."""
    pitch = float(mesh.pitch_norm)
    shift = np.array([pitch, 0.0])
    n_b = len(mesh.streamlines) - len(mesh.separatrices)
    bnd = [np.asarray(s, float) for s in mesh.streamlines[:n_b]]
    seps = [np.asarray(s, float) for s in mesh.streamlines[n_b:]]

    A, B = _periodic_chains(mesh)
    side = []
    for poly in bnd:
        dA = max(_dist_to_chain(q, A) for q in poly)
        dB = max(_dist_to_chain(q, B) for q in poly)
        side.append("L" if dA < tol_cls else "R" if dB < tol_cls else None)
    left = [bnd[i] for i in range(n_b) if side[i] == "L"]
    right = [bnd[i] for i in range(n_b) if side[i] == "R"]
    other = [bnd[i] for i in range(n_b) if side[i] is None]
    if not left or not right:
        print("[seam] walls not found -- symmetrization skipped")
        return None

    WL = _chain_segments(left)
    WR = _chain_segments(right)
    conf = max(np.linalg.norm(q - shift - ps._project_to_polyline(
        q - shift, WL)[3]) for q in WR)
    JL = _internal_junctions(left, [WL[0], WL[-1]])
    JR = _internal_junctions(right, [WR[0], WR[-1]])

    # canonical junction set on the master wall
    canon = []                                   # dicts: coord, srcL, srcR
    for p in JL:
        proj = ps._project_to_polyline(p, WL)[3]
        canon.append({"coord": np.asarray(proj, float), "srcL": p, "srcR": None})
    for p in JR:
        q = p - shift
        hit = next((c for c in canon
                    if np.linalg.norm(c["coord"] - q) < tol_match), None)
        if hit is None:
            proj = ps._project_to_polyline(q, WL)[3]
            canon.append({"coord": np.asarray(proj, float),
                          "srcL": None, "srcR": p})
        else:
            hit["srcR"] = p

    left_segs = _split_polyline_at(WL, [c["coord"] for c in canon])
    right_segs = [seg + shift for seg in left_segs]

    # re-terminate interior curves exactly on the (possibly moved) junctions
    moved, max_move = 0, 0.0
    for c in canon:
        for sd, old in (("L", c["srcL"]), ("R", c["srcR"])):
            if old is None:
                continue
            new = c["coord"] if sd == "L" else c["coord"] + shift
            d = float(np.linalg.norm(new - old))
            for k, s in enumerate(seps):
                if np.linalg.norm(s[-1] - old) < 1e-6:
                    seps[k] = np.vstack([s[:-1], new])
                elif np.linalg.norm(s[0] - old) < 1e-6:
                    seps[k] = np.vstack([new, s[1:]])
                else:
                    continue
                moved += 1
                max_move = max(max_move, d)

    mesh.streamlines = other + left_segs + right_segs + seps
    n_mirror = sum(1 for c in canon if c["srcL"] is None or c["srcR"] is None)
    print(f"[seam] junctions: L={len(JL)} R={len(JR)} union={len(canon)} "
          f"({n_mirror} hanging mirrors); wall conformity pre={conf:.2e} "
          f"post=0 (copy+shift); re-terminated {moved} curve ends "
          f"(max move {max_move:.2e})")
    return {"n_left": len(JL), "n_right": len(JR), "n_union": len(canon),
            "n_hanging": n_mirror, "wall_conformity_pre": conf,
            "max_endpoint_move": max_move, "WL": WL, "WR": WR}


# --------------------------------------------------------------------------
# TFI (Coons) with canonical per-edge sampling
# --------------------------------------------------------------------------

def _arclen(poly):
    return float(np.linalg.norm(np.diff(poly, axis=0), axis=1).sum())


def _resample(poly, n):
    poly = np.asarray(poly, float)
    seg = np.linalg.norm(np.diff(poly, axis=0), axis=1)
    d = np.concatenate([[0.0], np.cumsum(seg)])
    if d[-1] < 1e-15:
        return np.repeat(poly[:1], n, axis=0)
    ss = np.linspace(0.0, d[-1], n)
    return np.column_stack([np.interp(ss, d, poly[:, 0]),
                            np.interp(ss, d, poly[:, 1])])


def _coons(S, N, W, E):
    """S: c0->c1, N: c3->c2 (both len n_u); W: c0->c3, E: c1->c2 (len n_v).
    Returns grid (n_v, n_u, 2)."""
    n_u, n_v = len(S), len(W)
    u = np.linspace(0.0, 1.0, n_u)[None, :, None]
    v = np.linspace(0.0, 1.0, n_v)[:, None, None]
    c00, c10, c01, c11 = S[0], S[-1], N[0], N[-1]
    X = ((1 - v) * S[None, :, :] + v * N[None, :, :]
         + (1 - u) * W[:, None, :] + u * E[:, None, :]
         - ((1 - u) * (1 - v) * c00 + u * (1 - v) * c10
            + (1 - u) * v * c01 + u * v * c11))
    return X


def _inverted_cells(X):
    """Count TFI cells with non-positive signed area (shoelace)."""
    p00 = X[:-1, :-1]; p10 = X[:-1, 1:]; p11 = X[1:, 1:]; p01 = X[1:, :-1]
    area = 0.5 * ((p00[..., 0] * p10[..., 1] - p10[..., 0] * p00[..., 1])
                  + (p10[..., 0] * p11[..., 1] - p11[..., 0] * p10[..., 1])
                  + (p11[..., 0] * p01[..., 1] - p01[..., 0] * p11[..., 1])
                  + (p01[..., 0] * p00[..., 1] - p00[..., 0] * p01[..., 1]))
    return int((area <= 0).sum()), int(area.size)


def tfi_fill(result, WL, WR, pitch, h=0.04, seam_tol=0.01):
    """Coons patch per block. Every graph edge gets ONE canonical sample set
    (count from its arclength); a block side is the concatenation of its edge
    samples, so blocks sharing a full side share the discretization exactly and
    the right-seam sides are exact +pitch translates of the left ones (walls
    were rebuilt as copies). A side is only re-sampled when the opposite side
    of the same block dictates a different count (T-mesh hanging refinement);
    seam sides always keep the canonical sampling (master-slave)."""
    e2s = result["e2s"]
    edge_samples = {}
    for key, poly in e2s.items():
        n = max(2, int(np.ceil(_arclen(poly) / h)) + 1)
        edge_samples[key] = _resample(poly, n)

    def _edge_s(a, b):
        if (a, b) in edge_samples:
            return edge_samples[(a, b)]
        return edge_samples[(b, a)][::-1]

    def _side_samples(chain):
        parts = [_edge_s(chain[j], chain[j + 1])
                 for j in range(len(chain) - 1)]
        out = [parts[0]]
        out.extend(p[1:] for p in parts[1:])
        return np.vstack(out)

    def _seam_flag(poly):
        if max(_dist_to_chain(q, WL) for q in poly[:: max(1, len(poly) // 8)]) \
                < seam_tol and _dist_to_chain(poly[len(poly) // 2], WL) < seam_tol:
            return "L"
        if max(_dist_to_chain(q, WR) for q in poly[:: max(1, len(poly) // 8)]) \
                < seam_tol and _dist_to_chain(poly[len(poly) // 2], WR) < seam_tol:
            return "R"
        return None

    grids, n_inv, n_cells = [], 0, 0
    seam_pts = {"L": [], "R": []}
    for blk in result["blocks"]:
        samp = [_side_samples(ch) for ch in blk["side_chains"]]
        flag = [_seam_flag(s) for s in blk["sides"]]

        def _dim(i, j):
            if flag[i]:
                return len(samp[i])
            if flag[j]:
                return len(samp[j])
            return max(len(samp[i]), len(samp[j]))

        n_u, n_v = _dim(0, 2), _dim(1, 3)

        def _fit(arr, n):
            return arr if len(arr) == n else _resample(arr, n)

        S = _fit(samp[0], n_u)
        E = _fit(samp[1], n_v)
        N = _fit(samp[2], n_u)[::-1]      # side2: c2->c3, Coons wants c3->c2
        W = _fit(samp[3], n_v)[::-1]      # side3: c3->c0, Coons wants c0->c3
        X = _coons(S, N, W, E)
        inv, tot = _inverted_cells(X)
        n_inv += inv
        n_cells += tot
        grids.append(X)
        for i in range(4):
            if flag[i]:
                seam_pts[flag[i]].append(samp[i])

    # seam conformity of the TFI discretization
    seam_dev, n_l, n_r = None, 0, 0
    if seam_pts["L"] and seam_pts["R"]:
        Lp = np.unique(np.round(np.vstack(seam_pts["L"]), 12), axis=0)
        Rp = np.unique(np.round(np.vstack(seam_pts["R"]), 12), axis=0) \
            - np.array([pitch, 0.0])
        n_l, n_r = len(Lp), len(Rp)
        d = np.linalg.norm(Lp[:, None, :] - Rp[None, :, :], axis=2)
        seam_dev = float(max(d.min(axis=1).max(), d.min(axis=0).max()))
    return {"grids": grids, "inverted_cells": n_inv, "total_cells": n_cells,
            "seam_dev": seam_dev, "seam_nodes_lr": (n_l, n_r),
            "seam_pts": seam_pts}


# --------------------------------------------------------------------------
# pipeline
# --------------------------------------------------------------------------

def run_tmesh(stl=STL, out_dir=OUT, verbose=True, make_plots=True):
    ps.set_periodic(True)          # field seam weld stays ON
    ps.set_tile_periodic(False)    # block stage: seam = wall
    t0 = time.time()

    mesh, transform = build_dp_data(stl)
    pitch = float(mesh.pitch_norm)
    ff = FrameField(mesh)
    m = detect_singularities(ff.mesh)
    n_sing = int((m.singularities != 0).sum())

    sl = StreamlineGenerator_v2(ff.mesh)
    ps._drop_degenerate_corner_seps(sl.mesh)
    ps._close_helical_streamlines(sl.mesh)     # inert (TILE_PERIODIC off)
    ps._emit_dock_crossings(sl)                # inert (no docks)
    ps._snap_separatrix_endpoints(sl.mesh, radius=0.045)

    collapse_seam_wedges(sl.mesh)
    seam_info = symmetrize_seam_junctions(sl.mesh)
    n_boundary = len(sl.mesh.streamlines) - len(sl.mesh.separatrices)
    boundary_ref = [np.asarray(s, float)
                    for s in sl.mesh.streamlines[:n_boundary]]

    merging = StreamlineMerging(sl.mesh, verbose=False)
    splitter = StreamlineIntersectionSplitter(offset_boundingBox=0.05,
                                              num_samples=5)
    updated = splitter.process_streamlines(merging.new_streamlines)

    # flat_tol 15: hanging seam mirrors turn by ~0.2deg, while the flattest
    # genuine corner (a tip-cluster singularity that lost its wedge arm) still
    # turns by ~24deg -- 15 separates the two populations with margin.
    gen = tmf.TMeshFaceGenerator(updated, blade_loops=list(mesh.blade_loops),
                                 flat_tol_deg=15.0, verbose=verbose)
    result = gen.get_blocks()
    for rej in result["rejects"]:
        nds = result["nodes"][rej["cycle"]] if rej.get("cycle") else []
        turns = (tmf.corner_turns(rej["cycle"], result["e2s"])
                 if rej.get("cycle") and rej.get("ring") is not None else [])
        print(f"[tmesh] reject n_real={rej['n_real']} cycle={rej['cycle']} "
              f"nodes={np.round(np.asarray(nds), 3).tolist()} "
              f"turns={np.round(turns, 1).tolist()}")

    def _bdist(p):
        return ps._min_boundary_dist(p, boundary_ref)

    regular, tnodes, irregular = tmf.node_regularity(result, _bdist)

    WL, WR = (seam_info["WL"], seam_info["WR"]) if seam_info else (None, None)
    tfi = tfi_fill(result, WL, WR, pitch) if seam_info else None

    # block-corner seam conformity (block corners on the seams, matched mod pitch)
    corner_dev, corner_lr = None, (0, 0)
    if seam_info and result["blocks"]:
        cn = np.unique(np.concatenate([b["corners"] for b in result["blocks"]]))
        pts = result["nodes"][cn]
        Lc = np.array([p for p in pts if _dist_to_chain(p, WL) < 0.01])
        Rc = np.array([p for p in pts if _dist_to_chain(p, WR) < 0.01])
        corner_lr = (len(Lc), len(Rc))
        if len(Lc) and len(Rc):
            d = np.linalg.norm(Lc[:, None, :]
                               - (Rc[None, :, :] - np.array([pitch, 0.0])),
                               axis=2)
            corner_dev = float(max(d.min(axis=1).max(), d.min(axis=0).max()))

    runtime = time.time() - t0
    metrics = {
        "approach": "T (wall + T-mesh + TFI)",
        "singularities": n_sing,
        "blocks": len(result["blocks"]),
        "rejected_regions": len(result["rejects"]),
        "t_nodes": len(tnodes),
        "irregular_interior_nodes": len(irregular),
        "inverted_blocks": sum(1 for b in result["blocks"] if b["area"] <= 0),
        "inverted_tfi_cells": tfi["inverted_cells"] if tfi else None,
        "total_tfi_cells": tfi["total_cells"] if tfi else None,
        "seam_corner_lr": list(corner_lr),
        "seam_corner_dev": corner_dev,
        "seam_tfi_nodes_lr": list(tfi["seam_nodes_lr"]) if tfi else None,
        "seam_tfi_dev": tfi["seam_dev"] if tfi else None,
        "seam_junctions": {k: seam_info[k] for k in
                           ("n_left", "n_right", "n_union", "n_hanging")}
        if seam_info else None,
        "planar": result["planar"],
        "runtime_s": round(runtime, 1),
    }
    if verbose:
        print(f"[tmesh] GATE regions!=4 real corners: "
              f"{metrics['rejected_regions']} "
              f"({'PASS' if metrics['rejected_regions'] == 0 else 'FAIL'})")
        print(f"[tmesh] GATE seam conformity: corners {corner_lr} "
              f"max dev={corner_dev}  TFI dev={metrics['seam_tfi_dev']}")
        print(f"[tmesh] GATE inverted TFI cells: "
              f"{metrics['inverted_tfi_cells']}/{metrics['total_tfi_cells']}")
        print(f"[tmesh] {metrics['blocks']} blocks, {len(tnodes)} T-nodes, "
              f"{len(irregular)} irregular interior, {runtime:.0f}s")

    if make_plots:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        plot_blocks(result, mesh, tnodes, irregular, seam_info,
                    out_dir / "tmesh_blocks.png")
        if tfi:
            plot_tfi(tfi, result, pitch, out_dir / "tmesh_tfi.png")
        (out_dir / "tmesh_metrics.json").write_text(
            json.dumps(metrics, indent=2))
        print(f"wrote {out_dir}/tmesh_metrics.json")
    return {"metrics": metrics, "result": result, "tfi": tfi, "mesh": mesh,
            "seam_info": seam_info, "boundary_ref": boundary_ref}


# --------------------------------------------------------------------------
# plots
# --------------------------------------------------------------------------

def plot_blocks(result, mesh, tnodes, irregular, seam_info, out_png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    xy = mesh.x[:, 0:2].numpy()
    tris = mesh.faces.T.numpy()
    nodes = result["nodes"]
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.triplot(xy[:, 0], xy[:, 1], tris, lw=0.1, color="0.93")
    for blk in result["blocks"]:
        ring = blk["ring"]
        ax.fill(ring[:, 0], ring[:, 1], alpha=0.22, color="C0")
        for s in blk["sides"]:
            ax.plot(s[:, 0], s[:, 1], "C0", lw=1.2)
        c = nodes[blk["corners"]]
        ax.scatter(c[:, 0], c[:, 1], c="k", s=14, zorder=6)
    for rej in result["rejects"]:
        if rej.get("ring") is not None:
            r = rej["ring"]
            ax.plot(r[:, 0], r[:, 1], "red", lw=1.8)
            ax.fill(r[:, 0], r[:, 1], color="red", alpha=0.25)
    if tnodes:
        ax.scatter(nodes[tnodes, 0], nodes[tnodes, 1], marker="s",
                   facecolors="none", edgecolors="darkorange", s=70,
                   linewidths=1.6, zorder=7,
                   label=f"T-node (hanging): {len(tnodes)}")
    if irregular:
        ax.scatter(nodes[irregular, 0], nodes[irregular, 1],
                   facecolors="none", edgecolors="red", s=150, linewidths=2.0,
                   zorder=7, label=f"irregular interior: {len(irregular)}")
    if seam_info:
        for W in (seam_info["WL"], seam_info["WR"]):
            ax.plot(W[:, 0], W[:, 1], "green", lw=0.8, alpha=0.7)
    ax.set_aspect("equal")
    if tnodes or irregular:
        ax.legend(loc="upper left", fontsize=8)
    ax.set_title(f"Ansatz T: {len(result['blocks'])} T-mesh blocks "
                 f"({len(result['rejects'])} rejected, "
                 f"{len(tnodes)} T-nodes, {len(irregular)} irregular)")
    fig.savefig(out_png, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_png}")


def plot_tfi(tfi, result, pitch, out_png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(10, 8))
    for X in tfi["grids"]:
        for i in range(X.shape[0]):
            ax.plot(X[i, :, 0], X[i, :, 1], "0.4", lw=0.35)
        for j in range(X.shape[1]):
            ax.plot(X[:, j, 0], X[:, j, 1], "0.4", lw=0.35)
    for blk in result["blocks"]:
        for s in blk["sides"]:
            ax.plot(s[:, 0], s[:, 1], "C0", lw=1.0)
    for sd, col in (("L", "green"), ("R", "red")):
        if tfi["seam_pts"][sd]:
            P = np.vstack(tfi["seam_pts"][sd])
            ax.scatter(P[:, 0], P[:, 1], c=col, s=10, zorder=6,
                       label=f"seam nodes {sd}: {len(np.unique(np.round(P, 9), axis=0))}")
    ax.set_aspect("equal")
    ax.legend(loc="upper left", fontsize=8)
    dev = tfi["seam_dev"]
    ax.set_title(f"Ansatz T TFI: {len(tfi['grids'])} blocks, "
                 f"{tfi['inverted_cells']}/{tfi['total_cells']} inverted cells, "
                 f"seam dev={dev if dev is None else f'{dev:.1e}'}")
    fig.savefig(out_png, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_png}")


if __name__ == "__main__":
    us.set_blade_tip_corners(False)
    cs.set_emanate_outer_corners(True)
    run_tmesh()
