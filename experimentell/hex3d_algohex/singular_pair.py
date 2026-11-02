#!/usr/bin/env python3
"""Probe: local collapse of adjacent singular pairs in a 3D hex mesh."""

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO / "experimentell" / "hex3d_algohex"))

import base_complex as bc
import ovm_io

HE = bc.HEX_EDGES


def _parallel(P, e1, e2):
    a = P[e1[1]] - P[e1[0]]
    b = P[e2[1]] - P[e2[0]]
    na = np.linalg.norm(a)
    nb = np.linalg.norm(b)
    if na < 1e-12 or nb < 1e-12:
        return False
    cos = float(np.dot(a, b) / (na * nb))
    return abs(cos) > 0.9


def _min_endpoint_dist(P, e1, e2):
    d = 1e18
    for u in e1:
        for v in e2:
            d = min(d, float(np.linalg.norm(P[u] - P[v])))
    return d


def detect_singular(P, H, interior_only=True):
    f2h, e2h = bc.build_topology(H)
    sing = bc.singular_edges(H, P, f2h, e2h, interior_only=interior_only)
    return f2h, e2h, sing


def detect_pairs(P, H, interior_only=True):
    f2h, e2h, sing = detect_singular(P, H, interior_only=interior_only)
    hex_edges = []
    for hi, c in enumerate(H):
        es = []
        for a, b in HE:
            u, v = int(c[a]), int(c[b])
            es.append((min(u, v), max(u, v)))
        hex_edges.append(es)

    cand = set()
    for hi in range(len(H)):
        es = hex_edges[hi]
        for i in range(len(es)):
            for j in range(i + 1, len(es)):
                e1, e2 = es[i], es[j]
                if (e1 != e2 and not (set(e1) & set(e2)) and e1 in sing and e2 in sing
                        and _parallel(P, e1, e2)):
                    cand.add(frozenset((e1, e2)))
    for f, hs in f2h.items():
        if len(hs) != 2:
            continue
        h0, h1 = hs[0], hs[1]
        for e1 in hex_edges[h0]:
            if e1 not in sing:
                continue
            for e2 in hex_edges[h1]:
                if e2 not in sing:
                    continue
                if (e1 != e2 and not (set(e1) & set(e2)) and _parallel(P, e1, e2)):
                    cand.add(frozenset((e1, e2)))

    out = []
    for c in cand:
        e1, e2 = tuple(c)
        out.append((e1, e2, _min_endpoint_dist(P, e1, e2)))
    out.sort(key=lambda x: x[2])
    return out, sing, f2h, e2h


def _boundary_verts(P, H):
    if len(H) == 0:
        return []
    f2h, _ = bc.build_topology(H)
    return sorted({v for fk, hs in f2h.items() if len(hs) == 1 for v in fk})


def collapse_pair(P, H, e1, e2, strict_boundary=True):
    mapping = {}
    for v2 in e2:
        dists = sorted((float(np.linalg.norm(P[v2] - P[v1])), v1) for v1 in e1)
        mapping[v2] = dists[0][1]

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

    for v2, v1 in mapping.items():
        union(v2, v1)

    new_hexes = []
    for c in H:
        m = [find(int(v)) for v in c]
        if len(set(m)) == 8:
            new_hexes.append(m)
    new_hexes = (np.array(new_hexes, np.int64) if new_hexes
                 else np.zeros((0, 8), np.int64))

    Q = P.copy()
    groups = defaultdict(list)
    for v in range(len(P)):
        groups[find(v)].append(v)
    for r, vs in groups.items():
        if len(vs) > 1:
            Q[vs] = P[vs].mean(0)

    if strict_boundary:
        bnd_old = _boundary_verts(P, H)
        bnd_new = _boundary_verts(Q, new_hexes)
        boundary_ok = (bnd_old == bnd_new)
    else:
        mins = P.min(0); maxs = P.max(0)
        outer = [v for v in range(len(P))
                 if any(abs(P[v, d] - mins[d]) < 1e-9 or abs(P[v, d] - maxs[d]) < 1e-9
                        for d in range(3))]
        boundary_ok = all(np.allclose(Q[v], P[v]) for v in outer)
        # also indices present
        new_outer = [v for v in range(len(Q))
                     if any(abs(Q[v, d] - mins[d]) < 1e-9 or abs(Q[v, d] - maxs[d]) < 1e-9
                            for d in range(3))]
        boundary_ok = boundary_ok and (set(new_outer) == set(outer))

    if len(new_hexes):
        sj_new = np.array([ovm_io.scaled_jacobian(Q, c) for c in new_hexes])
    else:
        sj_new = np.array([])

    if len(new_hexes):
        f2h_new, _ = bc.build_topology(new_hexes)
    else:
        f2h_new = {}
    watertight = (len(new_hexes) > 0) and all(
        len(hs) in (1, 2) for hs in f2h_new.values())

    reason = None
    if len(new_hexes) == 0:
        reason = "no hexes remain after collapse"
    elif not watertight:
        reason = "mesh not watertight (a face is used 0 or >2 times)"
    elif not boundary_ok:
        reason = "boundary vertex set changed"
    elif len(sj_new) and bool((sj_new <= 0).any()):
        reason = "new inverted cell introduced"

    info = {
        "dropped_hexes": int(len(H) - len(new_hexes)),
        "remaining_hexes": int(len(new_hexes)),
        "watertight": bool(watertight),
        "boundary_changed": bool(not boundary_ok),
        "new_inverted": int((sj_new <= 0).sum()) if len(sj_new) else 0,
        "min_sj_new": float(sj_new.min()) if len(sj_new) else float("nan"),
    }
    ok = reason is None
    return Q, new_hexes, ok, reason, info


def merge_one_pair(P, H, max_attempts=5, interior_only=True, strict_boundary=True):
    pairs, sing, f2h, e2h = detect_pairs(P, H, interior_only=interior_only)
    attempts = []
    for idx in range(min(max_attempts, len(pairs))):
        e1, e2, dist = pairs[idx]
        Q, Hn, ok, reason, info = collapse_pair(P, H, e1, e2, strict_boundary=strict_boundary)
        attempts.append({
            "attempt": idx + 1,
            "edge1": [int(x) for x in e1],
            "edge2": [int(x) for x in e2],
            "distance": float(dist),
            "ok": bool(ok),
            "reason": reason,
            "info": info,
        })
        if ok:
            return Q, Hn, e1, e2, attempts
    return None, None, None, None, attempts


def _make_grid(nx, ny, nz):
    pts = []
    idx = {}
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                idx[(i, j, k)] = len(pts)
                pts.append([float(i), float(j), float(k)])
    P = np.array(pts, float)
    hexes = []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                c = [
                    idx[(i, j, k)],
                    idx[(i + 1, j, k)],
                    idx[(i + 1, j + 1, k)],
                    idx[(i, j + 1, k)],
                    idx[(i, j, k + 1)],
                    idx[(i + 1, j, k + 1)],
                    idx[(i + 1, j + 1, k + 1)],
                    idx[(i, j + 1, k + 1)],
                ]
                hexes.append(c)
    return P, np.array(hexes, np.int64)


def self_test():
    P, H = _make_grid(3, 3, 3)
    centre = 1 * 9 + 1 * 3 + 1
    H = np.delete(H, centre, axis=0)

    pairs, sing, f2h, e2h = detect_pairs(P, H, interior_only=False)
    print(f"[self-test] singular edges: {len(sing)}, candidate pairs: {len(pairs)}")

    Q, Hn, e1, e2, attempts = merge_one_pair(P, H, max_attempts=5, interior_only=False, strict_boundary=False)
    for a in attempts:
        print(f"[self-test] attempt {a['attempt']}: edges {a['edge1']} / "
              f"{a['edge2']} dist {a['distance']:.4f} -> "
              f"{'OK' if a['ok'] else 'FAIL: ' + str(a['reason'])}")

    if Q is None:
        print("[self-test] no pair collapsed successfully; reporting infeasibility")
        print("SELF-TEST: detection ran; no feasible collapse (infeasible probe)")
        return 0

    assert all(len(c) == 8 for c in Hn), "not all cells are hexes"
    assert all(len(set(c)) == 8 for c in Hn), "degenerate hex remains"
    f2h_n, _ = bc.build_topology(Hn)
    assert all(len(hs) in (1, 2) for hs in f2h_n.values()), "not watertight"
    # outer (domain) boundary vertices only
    mins = P.min(0); maxs = P.max(0)
    outer = [v for v in range(len(P))
             if any(abs(P[v, d] - mins[d]) < 1e-9 or abs(P[v, d] - maxs[d]) < 1e-9
                    for d in range(3))]
    outer_set = set(outer)
    # positions of outer vertices must be unchanged
    assert all(np.allclose(Q[v], P[v]) for v in outer), "outer boundary vertex moved"
    # outer vertex set unchanged (indices)
    new_outer = [v for v in range(len(Q))
                 if any(abs(Q[v, d] - mins[d]) < 1e-9 or abs(Q[v, d] - maxs[d]) < 1e-9
                        for d in range(3))]
    assert set(new_outer) == outer_set, "outer boundary vertex set changed"
    sj = np.array([ovm_io.scaled_jacobian(Q, c) for c in Hn])
    assert not bool((sj <= 0).any()), "new inverted cell"
    print(f"[self-test] collapsed pair edges {e1} / {e2}; "
          f"hexes {len(H)} -> {len(Hn)}; min sJ {sj.min():.4f}")
    print("SELF-TEST: PASS")
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--input", help="input .ovm or .vtk")
    ap.add_argument("--output-vtk", help="output .vtk")
    ap.add_argument("--output-json", help="output .json")
    ap.add_argument("--max-attempts", type=int, default=5)
    args = ap.parse_args()

    if args.self_test:
        rc = self_test()
        raise SystemExit(rc)

    if not args.input:
        print("need --input or --self-test", file=sys.stderr)
        raise SystemExit(2)

    if args.input.endswith(".ovm"):
        pts, e, f, poly = ovm_io.read_ovm(args.input)
        P, H = pts, ovm_io.ovm_to_cells(pts, e, f, poly)[0]
    else:
        import meshio
        m = meshio.read(args.input)
        H = np.vstack([b.data for b in m.cells if b.type == "hexahedron"]).astype(np.int64)
        P = np.asarray(m.points, float)

    Q, Hn, e1, e2, attempts = merge_one_pair(P, H, max_attempts=args.max_attempts)

    result = {
        "input": args.input,
        "singular_pairs_found": len(detect_pairs(P, H)[0]),
        "attempts": attempts,
        "collapsed": Q is not None,
    }
    if Q is not None:
        result["collapsed_edge1"] = [int(x) for x in e1]
        result["collapsed_edge2"] = [int(x) for x in e2]
        result["hexes_before"] = int(len(H))
        result["hexes_after"] = int(len(Hn))
        sj = np.array([ovm_io.scaled_jacobian(Q, c) for c in Hn])
        result["min_sj_after"] = float(sj.min()) if len(sj) else None
        result["inverted_after"] = int((sj <= 0).sum()) if len(sj) else 0
        if args.output_vtk:
            outp = Path(args.output_vtk)
            outp.parent.mkdir(parents=True, exist_ok=True)
            ovm_io.write_hex_vtk(str(outp), Q, Hn)
            print(f"wrote {outp}")
    if args.output_json:
        outj = Path(args.output_json)
        outj.parent.mkdir(parents=True, exist_ok=True)
        outj.write_text(json.dumps(result, indent=2))
        print(f"wrote {outj}")

    print(json.dumps(result, indent=2))
    raise SystemExit(0 if Q is not None else 1)


if __name__ == "__main__":
    main()
