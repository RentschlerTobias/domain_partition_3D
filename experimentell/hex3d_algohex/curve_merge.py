#!/usr/bin/env python3
"""Whole-curve singular merge experiment for canada.

Implements collapse_curve_pair: collapse the hex sheet between two parallel
singular curves in a single pass, using a vertex correspondence built along
the two curves (not the wrong nearest-endpoint weld from T68). Boundary
vertices are kept fixed (excluded from merging).
"""

import argparse
import json
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO / "experimentell" / "hex3d_algohex"))

import base_complex as bc
import ovm_io
import singular_pair as sp


HEX_FACES = bc.HEX_FACES


def _boundary_verts(P, H):
    return set(sp._boundary_verts(P, H))


def curve_vertex_path(curve):
    """Return ordered list of vertex indices along a curve (list of edges)."""
    from collections import Counter, defaultdict as dd
    cnt = Counter()
    adj = dd(list)
    for ei, e in enumerate(curve):
        cnt[e[0]] += 1
        cnt[e[1]] += 1
        adj[e[0]].append((ei, e[1]))
        adj[e[1]].append((ei, e[0]))
    starts = [v for v, d in cnt.items() if d == 1]
    if not starts:
        starts = [curve[0][0]] if curve else []
    if not starts:
        return []
    start = starts[0]
    path = [start]
    used = set()
    cur = start
    while True:
        nxt = None
        for ei, nb in adj[cur]:
            if ei in used:
                continue
            nxt = (ei, nb)
            break
        if nxt is None:
            break
        ei, nb = nxt
        used.add(ei)
        path.append(nb)
        cur = nb
        if len(used) == len(curve):
            break
    return path


def _orient_curves(curve_a, curve_b, P):
    pa = curve_vertex_path(curve_a)
    pb = curve_vertex_path(curve_b)
    best = None
    for pa_o in (pa, list(reversed(pa))):
        for pb_o in (pb, list(reversed(pb))):
            d = (float(np.linalg.norm(P[pa_o[0]] - P[pb_o[0]]))
                 + float(np.linalg.norm(P[pa_o[-1]] - P[pb_o[-1]])))
            if best is None or d < best[0]:
                best = (d, pa_o, pb_o)
    return best[1], best[2]


def _arc_params(path, P):
    seg = [0.0]
    for i in range(1, len(path)):
        seg.append(seg[-1] + float(np.linalg.norm(P[path[i]] - P[path[i - 1]])))
    total = seg[-1]
    if total < 1e-12:
        return [0.0] * len(path), 0.0
    return [x / total for x in seg], total


def build_correspondence(curve_a, curve_b, P, boundary_verts=None):
    """Build vertex mapping from curve_b (interior verts) to curve_a."""
    if boundary_verts is None:
        boundary_verts = set()
    bnd = set(boundary_verts)
    path_a, path_b = _orient_curves(curve_a, curve_b, P)
    sa, _ = _arc_params(path_a, P)
    sb, _ = _arc_params(path_b, P)
    mapping = {}
    for idx_b, vb in enumerate(path_b):
        if vb in bnd:
            continue
        # nearest in a by normalized arc length, but target must not be boundary
        best_i = None
        best_cost = 1e18
        for i, va in enumerate(path_a):
            if va in bnd:
                continue
            cost = abs(sa[i] - sb[idx_b])
            if cost < best_cost:
                best_cost = cost
                best_i = i
        if best_i is None:
            # fallback: allow boundary target? skip to keep boundary fixed
            continue
        va = path_a[best_i]
        mapping[int(vb)] = int(va)
    return mapping, path_a, path_b


def _hex_volume(P, cell):
    pts = P[cell]
    vol = 0.0
    for f in HEX_FACES:
        v0 = pts[f[0]]; v1 = pts[f[1]]; v2 = pts[f[2]]; v3 = pts[f[3]]
        vol += float(np.dot(v0, np.cross(v1, v2))) / 6.0
        vol += float(np.dot(v0, np.cross(v2, v3))) / 6.0
    return vol


def collapse_curve_pair(P, H, curve_a, curve_b, strict_boundary=True,
                        boundary_tol=1e-3):
    P = np.asarray(P, float).copy()
    H = np.asarray(H, np.int64).copy()
    bnd_verts = _boundary_verts(P, H)

    mapping, path_a, path_b = build_correspondence(curve_a, curve_b, P, bnd_verts)

    parent = list(range(len(P)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for vb, va in mapping.items():
        union(vb, va)

    new_hexes = []
    for c in H:
        m = [find(int(v)) for v in c]
        if len(set(m)) == 8:
            new_hexes.append(m)
    new_hexes = (np.array(new_hexes, dtype=np.int64)
                 if new_hexes else np.zeros((0, 8), dtype=np.int64))

    Q = P.copy()
    groups = defaultdict(list)
    for v in range(len(P)):
        groups[find(v)].append(v)
    for r, vs_ in groups.items():
        if len(vs_) > 1:
            # keep boundary vertices fixed if present
            bvs = [v for v in vs_ if v in bnd_verts]
            if bvs:
                Q[vs_] = P[bvs[0]]
            else:
                Q[vs_] = P[vs_].mean(axis=0)

    # boundary check
    if strict_boundary and len(new_hexes) > 0:
        boundary_ok, bnd_reason, bnd_info = sp._check_boundary_geometric(
            P, H, Q, new_hexes, tol=boundary_tol)
    elif strict_boundary and len(new_hexes) == 0:
        boundary_ok = False
        bnd_reason = "no hexes remain"
        bnd_info = {}
    else:
        boundary_ok = True
        bnd_reason = None
        bnd_info = {}

    if len(new_hexes):
        sj_new = np.array([ovm_io.scaled_jacobian(Q, c) for c in new_hexes])
        vol_new = np.array([_hex_volume(Q, c) for c in new_hexes])
    else:
        sj_new = np.array([])
        vol_new = np.array([])

    if len(new_hexes):
        f2h_new, _ = bc.build_topology(new_hexes)
    else:
        f2h_new = {}
    watertight = (len(new_hexes) > 0) and all(
        len(hs) in (1, 2) for hs in f2h_new.values())

    # surviving original curve edges
    def surviving_edges(curve):
        tot = 0
        surv = 0
        for e in curve:
            tot += 1
            u, v = int(e[0]), int(e[1])
            if find(u) != find(v):
                surv += 1
        return tot, surv

    tot_a, surv_a = surviving_edges(curve_a)
    tot_b, surv_b = surviving_edges(curve_b)

    reason = None
    if len(new_hexes) == 0:
        reason = "no hexes remain after collapse"
    elif not watertight:
        reason = "mesh not watertight (a face is used 0 or >2 times)"
    elif strict_boundary and not boundary_ok:
        reason = f"boundary surface not preserved: {bnd_reason}"
    elif len(sj_new) and bool((sj_new <= 0).any()):
        reason = "new inverted cell introduced (sJ)"
    elif len(vol_new) and bool((vol_new <= 0).any()):
        reason = "new inverted cell introduced (volume)"

    info = {
        "dropped_hexes": int(len(H) - len(new_hexes)),
        "remaining_hexes": int(len(new_hexes)),
        "watertight": bool(watertight),
        "boundary_changed": bool(not boundary_ok),
        "boundary_reason": bnd_reason,
        "boundary_info": bnd_info if strict_boundary else {},
        "new_inverted_sj": int((sj_new <= 0).sum()) if len(sj_new) else 0,
        "new_inverted_volume": int((vol_new <= 0).sum()) if len(vol_new) else 0,
        "min_sj_new": float(sj_new.min()) if len(sj_new) else float("nan"),
        "min_volume_new": float(vol_new.min()) if len(vol_new) else float("nan"),
        "mapped_vertices": int(len(mapping)),
        "curve_a_edges_total": int(tot_a),
        "curve_a_edges_surviving": int(surv_a),
        "curve_b_edges_total": int(tot_b),
        "curve_b_edges_surviving": int(surv_b),
        "curves_disappear": bool(surv_a == 0 and surv_b == 0),
        "boundary_verts_fixed_count": int(len(bnd_verts)),
    }
    ok = reason is None
    return Q, new_hexes, ok, reason, info, mapping, path_a, path_b


def _load_mesh(mesh_path):
    mesh_path = Path(mesh_path)
    if mesh_path.suffix == ".ovm":
        pts, e, f, poly = ovm_io.read_ovm(str(mesh_path))
        P, H = pts, ovm_io.ovm_to_cells(pts, e, f, poly)[0]
        return P.astype(float), H.astype(np.int64)
    else:
        import meshio
        m = meshio.read(str(mesh_path))
        H = np.vstack([b.data for b in m.cells if b.type == "hexahedron"]).astype(np.int64)
        P = np.asarray(m.points, float)
        return P, H


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mesh", default="/work/canada/sweep/n100000/hex_n100000.ovm")
    ap.add_argument("--curve-analysis", default="/work/canada/sweep/pair_probe3/curve_analysis.json")
    ap.add_argument("--output-json", default="/work/canada/sweep/pair_probe3/curve_merge_test.json")
    ap.add_argument("--pairs", default=None, help="comma list of curve pairs like '1-5,0-2,0-5'")
    args = ap.parse_args()

    P, H = _load_mesh(args.mesh)
    print(f"Loaded mesh: {len(P)} pts, {len(H)} hexes", flush=True)

    f2h, e2h = bc.build_topology(H)
    sing = bc.singular_edges(H, P, f2h, e2h, interior_only=True)
    import viz_singular as vs
    curves, parallel = vs.analyze_curves(P, H, f2h, e2h, sing)
    orig_n_curves = len(curves)
    orig_n_sing = len(sing)

    # pick pairs
    if args.pairs:
        picks = []
        for tok in args.pairs.split(","):
            a, b = tok.split("-")
            picks.append((int(a), int(b)))
    else:
        # default: top 3 by (tangent, then length match, then dist)
        # use parallel list sorted
        # prefer matching lengths
        scored = []
        for pp in parallel:
            ci, cj = pp["i"], pp["j"]
            c1, c2 = curves[ci], curves[cj]
            lenmatch = abs(len(c1) - len(c2))
            score = (pp["ang_deg"], lenmatch, pp["min_dist"])
            scored.append((score, ci, cj))
        scored.sort(key=lambda x: x[0])
        picks = [(ci, cj) for _, ci, cj in scored[:3]]

    attempts = []
    for (ci, cj) in picks:
        print(f"\n=== Pair curves {ci}-{cj} ===", flush=True)
        curve_a = curves[ci]
        curve_b = curves[cj]
        t0 = time.time()
        Q, Hn, ok, reason, info, mapping, pa, pb = collapse_curve_pair(
            P.copy(), H.copy(), curve_a, curve_b, strict_boundary=True)
        dt = time.time() - t0
        # recompute singular curves on collapsed mesh
        try:
            f2h_n, e2h_n = bc.build_topology(Hn)
            sing_n = bc.singular_edges(Hn, Q, f2h_n, e2h_n, interior_only=True)
            curves_n, _ = vs.analyze_curves(Q, Hn, f2h_n, e2h_n, sing_n)
            new_n_curves = len(curves_n)
            new_n_sing = len(sing_n)
        except Exception as ex:
            new_n_curves = None
            new_n_sing = None
        # boundary vertex count new
        new_bnd = len(sp._boundary_verts(Q, Hn)) if len(Hn) else 0
        orig_bnd = len(sp._boundary_verts(P, H))
        # hausdorff from boundary_info
        bi = info.get("boundary_info", {}) or {}
        hausdorff = bi.get("max_dist", bi.get("hausdorff_max", 0.0))
        # validity
        valid = bool(ok and len(Hn) > 0)
        attempt = {
            "pair": f"{ci}-{cj}",
            "curve_i": ci,
            "curve_j": cj,
            "curve_i_edges": len(curve_a),
            "curve_j_edges": len(curve_b),
            "tangent_deviation_deg": None,
            "min_dist": None,
            "elapsed_s": float(dt),
            "ok": bool(ok),
            "reason": reason,
            "valid_hex_mesh": bool(valid),
            "boundary": {
                "preserved": bool(info.get("boundary_changed") is False),
                "orig_boundary_verts": int(orig_bnd),
                "new_boundary_verts": int(new_bnd),
                "boundary_vertex_count_change": int(new_bnd - orig_bnd),
                "hausdorff_max": (float(hausdorff) if hausdorff is not None else None),
                "tol": 0.001,
            },
            "inversion": {
                "new_inverted_sj": int(info["new_inverted_sj"]),
                "new_inverted_volume": int(info["new_inverted_volume"]),
                "min_sj_new": float(info["min_sj_new"]),
                "min_volume_new": float(info["min_volume_new"]),
            },
            "hexes": {
                "before": int(len(H)),
                "after": int(len(Hn)),
                "removed": int(info["dropped_hexes"]),
            },
            "watertight": bool(info["watertight"]),
            "curves": {
                "i_edges_total": int(info["curve_a_edges_total"]),
                "i_edges_surviving": int(info["curve_a_edges_surviving"]),
                "j_edges_total": int(info["curve_b_edges_total"]),
                "j_edges_surviving": int(info["curve_b_edges_surviving"]),
                "both_disappear": bool(info["curves_disappear"]),
            },
            "mapped_vertices": int(info["mapped_vertices"]),
            "singular_curves": {
                "before": int(orig_n_curves),
                "after": (int(new_n_curves) if new_n_curves is not None else None),
                "delta": (int(new_n_curves - orig_n_curves) if new_n_curves is not None else None),
                "both_curves_merged": bool(new_n_curves is not None and new_n_curves <= orig_n_curves - 2),
            },
            "singular_edges": {
                "before": int(orig_n_sing),
                "after": (int(new_n_sing) if new_n_sing is not None else None),
                "delta": (int(new_n_sing - orig_n_sing) if new_n_sing is not None else None),
            },
            "info": info,
        }
        # fill tangent/dist from parallel list
        for pp in parallel:
            if pp["i"] == ci and pp["j"] == cj:
                attempt["tangent_deviation_deg"] = float(pp["ang_deg"])
                attempt["min_dist"] = float(pp["min_dist"])
        attempts.append(attempt)
        print(json.dumps(attempt, indent=2), flush=True)
        # stop after first pass
        if ok:
            print("PASS - stopping", flush=True)
            break

    # verdict
    any_pass = any(a["ok"] for a in attempts)
    # determine wall
    if any_pass:
        verdict = "YES: whole-curve merge possible on canada (at least one pair passed)."
        wall = None
    else:
        # classify wall
        walls = set()
        for a in attempts:
            r = a["reason"] or ""
            if "boundary" in r:
                walls.add("boundary")
            if "inverted" in r:
                walls.add("inversion")
            if "watertight" in r:
                walls.add("lattice")
        if "boundary" in walls:
            wall = "boundary"
        elif "inversion" in walls:
            wall = "inversion"
        else:
            wall = "lattice"
        verdict = f"NO: whole-curve merge not possible on canada; primary wall = {wall}."
    result = {
        "mesh": str(args.mesh),
        "attempts": attempts,
        "num_attempts": len(attempts),
        "any_pass": bool(any_pass),
        "verdict": verdict,
        "wall": wall,
        "notes": "Whole-curve collapse: vertex correspondence along curves, boundary vertices excluded from merging (kept fixed). Inversion checked with both scaled Jacobian and signed hex volume.",
    }
    outp = Path(args.output_json)
    outp.parent.mkdir(parents=True, exist_ok=True)
    outp.write_text(json.dumps(result, indent=2))
    print(f"\nWrote {outp}", flush=True)
    print(verdict, flush=True)


if __name__ == "__main__":
    main()
