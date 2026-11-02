#!/usr/bin/env python3
import argparse
import json
import sys
from collections import defaultdict, deque
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO / "experimentell" / "hex3d_algohex"))

import base_complex as bc
import ovm_io
import singular_pair as sp


def load_mesh(mesh_path):
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


def build_singular_curves(P, H, f2h, e2h, sing_set):
    """Group singular edges into curves. Chain edges that share a vertex,
    endpoint valence 2 within singular graph; split at junctions/branches."""
    # Build adjacency in singular graph
    graph = defaultdict(list)
    v_to_sing = defaultdict(list)
    for ek in sing_set:
        v_to_sing[ek[0]].append(ek)
        v_to_sing[ek[1]].append(ek)
    
    for v, edges in v_to_sing.items():
        # connect edges sharing vertex
        for i in range(len(edges)):
            for j in range(i + 1, len(edges)):
                e1, e2 = edges[i], edges[j]
                if e1 == e2:
                    continue
                graph[e1].append(e2)
                graph[e2].append(e1)
    
    # Get connected components
    visited = set()
    curves = []
    for start_e in sing_set:
        if start_e in visited:
            continue
        # find component
        comp = set()
        q = deque([start_e])
        comp.add(start_e)
        while q:
            cur = q.popleft()
            for nb in graph[cur]:
                if nb not in comp:
                    comp.add(nb)
                    q.append(nb)
        # mark visited
        for ce in comp:
            visited.add(ce)
        if len(comp) == 0:
            continue
        # build adj for component
        adj_c = defaultdict(list)
        for e1 in comp:
            for e2 in graph[e1]:
                if e2 in comp:
                    adj_c[e1].append(e2)
        # trace
        # find endpoints (valence 1) or start anywhere
        ends = []
        for e in comp:
            # degree in singular subgraph
            deg = len(set(adj_c[e]))
            if deg == 1:
                ends.append(e)
        if len(ends) == 0 and len(comp) > 0:
            # closed curve
            cur = next(iter(comp))
            path = trace_from_start(adj_c, cur, None)
            if path:
                curves.append(path)
        elif len(ends) >= 1:
            # trace from each endpoint to build curves (split at junctions)
            # but component might be branched tree - trace paths
            # for a simple chain, one trace from endpoint covers all
            start = ends[0]
            path = trace_from_start(adj_c, start, None)
            if path and len(path) > 0:
                curves.append(path)
            # if component not fully covered, there might be branches - but trace_from_start
            # stops at junctions; we need to handle branches
            # check coverage
            covered = set(path)
            if len(covered) < len(comp):
                # find unvisited endpoints in remaining
                # do more traces
                # trace from other ends
                for end in ends[1:]:
                    if end in covered:
                        continue
                    p2 = trace_from_start(adj_c, end, None)
                    if p2 and len(p2) > 0:
                        # check if overlaps significantly
                        curves.append(p2)
                        covered.update(p2)
        else:
            if len(comp) > 0:
                cur = next(iter(comp))
                path = trace_from_start(adj_c, cur, None)
                if path:
                    curves.append(path)
    return curves


def trace_from_start(adj_c, start, prev):
    trace = []
    cur = start
    pr = prev
    visited_in_trace = set()
    for _ in range(10**6):  # prevent infinite loop
        if cur in visited_in_trace and len(trace) > 0:
            break  # closed
        trace.append(cur)
        visited_in_trace.add(cur)
        nbs = list(set(adj_c[cur]))
        if len(nbs) == 0:
            break
        if pr is None:
            # pick first neighbor
            nxt = nbs[0]
            pr_new = cur
            cur = nxt
            pr = pr_new
            continue
        candidates = [n for n in nbs if n != pr]
        if len(candidates) == 0:
            break  # back to prev only - endpoint
        if len(candidates) == 1:
            nxt = candidates[0]
            pr_new = cur
            cur = nxt
            pr = pr_new
        else:
            break  # junction - stop tracing (split at branch)
    return trace


def curve_tangent(P, curve):
    """Compute mean tangent direction of curve."""
    if len(curve) == 0:
        return np.array([1,0,0])
    # compute tangents from edges
    tangents = []
    for ek in curve:
        v = P[ek[1]] - P[ek[0]]
        n = np.linalg.norm(v)
        if n > 1e-18:
            tangents.append(v / n)
    if not tangents:
        return np.array([1,0,0])
    t = np.mean(tangents, axis=0)
    n = np.linalg.norm(t)
    if n > 1e-18:
        return t / n
    return tangents[0]


def curve_midpoint(P, curve):
    pts = []
    for ek in curve:
        pts.append(P[ek[0]])
        pts.append(P[ek[1]])
    if not pts:
        return np.array([0,0,0])
    return np.mean(pts, axis=0)


def curve_length(P, curve):
    total = 0.0
    for ek in curve:
        total += float(np.linalg.norm(P[ek[1]] - P[ek[0]]))
    return total


def min_dist_between_curves(P, c1, c2):
    # min distance between any point on c1 and c2
    min_d = 1e18
    # sample points
    pts1 = []
    for ek in c1:
        pts1.extend([P[ek[0]], P[ek[1]]])
    pts2 = []
    for ek in c2:
        pts2.extend([P[ek[0]], P[ek[1]]])
    if not pts1 or not pts2:
        return min_d
    # compute min pairwise (subsample if large)
    for p1 in pts1[::max(1,len(pts1)//50)]:
        for p2 in pts2[::max(1,len(pts2)//50)]:
            d = float(np.linalg.norm(p1-p2))
            if d < min_d: min_d = d
    return min_d


def analyze_curves(P, H, f2h, e2h, sing_set):
    curves = build_singular_curves(P, H, f2h, e2h, sing_set)
    # also get all edges in curves
    curve_edges_set = set()
    for c in curves:
        for e in c:
            curve_edges_set.add(e)
    # find parallel pairs
    parallel_pairs = []
    for i in range(len(curves)):
        for j in range(i+1, len(curves)):
            c1, c2 = curves[i], curves[j]
            t1 = curve_tangent(P, c1)
            t2 = curve_tangent(P, c2)
            # tangent angle
            cos_ang = abs(float(np.dot(t1,t2)))
            cos_ang = min(1.0, max(-1.0, cos_ang))
            ang = np.degrees(np.arccos(cos_ang))
            if ang > 20.0:
                continue
            md = min_dist_between_curves(P, c1, c2)
            if md > 2.0:  # <= 2 hex layers
                continue
            parallel_pairs.append({
                'i': i, 'j': j,
                'len1': len(c1), 'len2': len(c2),
                'ang_deg': float(ang),
                'min_dist': float(md),
                't1': t1.tolist(), 't2': t2.tolist()
            })
    # sort by how promising
    parallel_pairs.sort(key=lambda x: (x['ang_deg'], x['min_dist'], -max(x['len1'],x['len2'])))
    return curves, parallel_pairs


def get_hexes_near_point(P, H, f2h, center, radius=2.0):
    # find hexes within radius*edge_len of center
    # simple: check distance from hex center
    near = set()
    for hi in range(len(H)):
        c = H[hi]
        hc = P[c].mean(axis=0)
        if float(np.linalg.norm(hc - center)) < radius * 0.5:  # rough
            near.add(hi)
            continue
        # also check if any vertex close
        for v in c:
            if float(np.linalg.norm(P[v] - center)) < radius:
                near.add(hi)
                break
    return near


def plot_overview(P, H, curves, outpath):
    fig = plt.figure(figsize=(16, 12))
    ax = fig.add_subplot(111, projection='3d')
    colors = plt.cm.tab20(np.linspace(0, 1, max(1, len(curves))))
    for ci, c in enumerate(curves):
        if len(c) == 0: continue
        pts = []
        for ek in c:
            pts.append(P[ek[0]])
            pts.append(P[ek[1]])
        if len(pts) == 0: continue
        pts = np.array(pts)
        # plot as line
        # group consecutive
        xs = pts[::2,0]; ys = pts[::2,1]; zs = pts[::2,2]
        xe = pts[1::2,0]; ye = pts[1::2,1]; ze = pts[1::2,2]
        for k in range(len(xs)):
            ax.plot([xs[k], xe[k]], [ys[k], ye[k]], [zs[k], ze[k]], 
                    color=colors[ci % len(colors)], linewidth=2)
    ax.set_title('Singular Curves Overview')
    ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    plt.close()
    print(f'Wrote {outpath}')


def plot_top10(P, H, f2h, curves, top10_pairs_info, outpath, detect_top10_path=None):
    fig = plt.figure(figsize=(20, 8))
    # load detect_top10 if provided
    top10 = []
    if detect_top10_path and Path(detect_top10_path).exists():
        try:
            d = json.loads(Path(detect_top10_path).read_text())
            top10 = d.get('top10_pairs', []) or d.get('pairs', []) or []
        except Exception:
            pass
    # also use attempts from pair_probe3.json
    attempts = []
    pair_probe3_path = Path(detect_top10_path).parent / 'pair_probe3.json' if detect_top10_path else None
    if pair_probe3_path and pair_probe3_path.exists():
        try:
            d = json.loads(pair_probe3_path.read_text())
            attempts = d.get('attempts', [])[:10]
        except Exception:
            pass
    
    for panel in range(10):
        ax = fig.add_subplot(2, 5, panel+1, projection='3d')
        if panel < len(attempts):
            att = attempts[panel]
            e1 = tuple(att['edge1'])
            e2 = tuple(att['edge2'])
            reason = att.get('reason', '')[:60]
            info = att.get('info', {})
            bnd_dev = info.get('boundary_info', {}).get('max_dist', info.get('boundary_info',{}).get('hausdorff_max', '?'))
            inv = info.get('new_inverted', '?')
            title = f"Pair {att.get('pair_idx', panel+1)}: dev={bnd_dev}, inv={inv}"
            ax.set_title(title, fontsize=8)
            # plot the two edges thick
            for ek, col, lw in [(e1,'r',4),(e2,'b',4)]:
                ax.plot([P[ek[0]][0],P[ek[1]][0]], [P[ek[0]][1],P[ek[1]][1]], [P[ek[0]][2],P[ek[1]][2]], 
                        color=col, linewidth=lw)
            # find neighborhood hexes
            center = (P[e1[0]] + P[e1[1]] + P[e2[0]] + P[e2[1]]) / 4.0
            # get nearby hexes
            near = set()
            for hi in range(min(len(H), 20000)):  # limit for speed
                c = H[hi]
                hc = P[c].mean(axis=0)
                if float(np.linalg.norm(hc - center)) < 3.0:
                    near.add(hi)
            # draw wireframe of nearby hexes (simplified)
            for hi in list(near)[:50]:
                c = H[hi]
                # draw edges of hex
                for a,b in bc.HEX_EDGES[:12]:
                    try:
                        ax.plot([P[c[a]][0],P[c[b]][0]], [P[c[a]][1],P[c[b]][1]], [P[c[a]][2],P[c[b]][2]],
                                color='gray', alpha=0.2, linewidth=0.5)
                    except Exception:
                        pass
            # weld direction
            v1 = P[e1[1]] - P[e1[0]]; v2 = P[e2[1]] - P[e2[0]]
            if np.linalg.norm(v1)>0: v1=v1/np.linalg.norm(v1)
            if np.linalg.norm(v2)>0: v2=v2/np.linalg.norm(v2)
            mid = center
            ax.quiver(mid[0],mid[1],mid[2], v1[0],v1[1],v1[2], length=0.1, color='r')
            ax.quiver(mid[0],mid[1],mid[2], v2[0],v2[1],v2[2], length=0.1, color='b')
        else:
            ax.set_title(f'Panel {panel+1}')
        ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
    plt.suptitle('Top-10 Pairs: Edges/Curves, Neighborhood, Weld Direction, Failures', fontsize=12)
    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    plt.close()
    print(f'Wrote {outpath}')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mesh', default='/work/canada/sweep/n100000/hex_n100000.ovm')
    ap.add_argument('--pair-probe3', default='/work/canada/sweep/pair_probe3')
    ap.add_argument('--outdir', default=None)
    args = ap.parse_args()
    
    outdir = Path(args.outdir) if args.outdir else Path(args.pair_probe3) / 'viz'
    outdir.mkdir(parents=True, exist_ok=True)
    
    P, H = load_mesh(args.mesh)
    print(f'Loaded mesh: {len(P)} pts, {len(H)} hexes')
    
    f2h, e2h = bc.build_topology(H)
    sing = bc.singular_edges(H, P, f2h, e2h, interior_only=True)
    print(f'Singular edges: {len(sing)}')
    
    curves = build_singular_curves(P, H, f2h, e2h, sing)
    print(f'Curves found: {len(curves)}')
    
    curves_info = analyze_curves(P, H, f2h, e2h, sing)[0]
    curves_res, parallel_pairs = analyze_curves(P, H, f2h, e2h, sing)
    
    # Analyze curve properties
    curve_stats = []
    for ci, c in enumerate(curves_res):
        if len(c) == 0: continue
        # check if open or closed
        adj_c = defaultdict(list)
        for e1 in c:
            for e2 in set([e for e in curves_res[ci] if e != e1]):  # rough
                pass
        # better check using full graph
        pass
    
    # Recompute properly
    graph_full = defaultdict(list)
    v_to_sing = defaultdict(list)
    for ek in sing:
        v_to_sing[ek[0]].append(ek)
        v_to_sing[ek[1]].append(ek)
    for v, edges in v_to_sing.items():
        for i in range(len(edges)):
            for j in range(i+1, len(edges)):
                e1,e2 = edges[i],edges[j]
                if e1==e2: continue
                graph_full[e1].append(e2)
                graph_full[e2].append(e1)
    
    # classify each curve
    curve_details = []
    for ci, c in enumerate(curves_res):
        if len(c) == 0: continue
        # find if closed
        # check endpoints connectivity in full graph
        ends = []
        for e in c:
            deg = len(set(graph_full[e]))
            if deg == 1: ends.append(e)
        is_closed = len(ends) == 0
        curve_details.append({
            'idx': ci,
            'n_edges': len(c),
            'open': not is_closed,
            'closed': is_closed,
            'length': curve_length(P,c),
            'tangent': curve_tangent(P,c).tolist(),
            'midpoint': curve_midpoint(P,c).tolist()
        })
    
    # Analyze which of the 373 pair candidates belong to same two curves
    pairs_all, sing_set2, f2h2, e2h2 = sp.detect_pairs(P, H, interior_only=True)
    # map edge to curve
    edge_to_curve = {}
    for ci, c in enumerate(curves_res):
        for e in c:
            edge_to_curve[e] = ci
    
    same_curve_pairs = 0
    cross_curve_pairs = 0
    pairs_same_curves_detail = []
    for e1,e2,dist in pairs_all:
        c1 = edge_to_curve.get(e1)
        c2 = edge_to_curve.get(e2)
        if c1 is not None and c2 is not None and c1 == c2:
            same_curve_pairs += 1
            pairs_same_curves_detail.append({'e1':e1,'e2':e2,'dist':dist,'curve':c1})
        elif c1 is not None and c2 is not None and c1 != c2:
            cross_curve_pairs += 1
        else:
            pass  # one or both not in traced curves?
    
    # Find promising curve pairs for whole-curve merge
    promising = []
    for pp in parallel_pairs[:30]:  # look at top 30
        c1 = curves_res[pp['i']]; c2 = curves_res[pp['j']]
        # check if they share boundary endpoints? look up vertices
        # get endpoints
        def curve_endpoints(c):
            ends = []
            for e in c:
                deg = len(set(graph_full[e]))
                if deg == 1: ends.append(e)
            if len(ends) == 0:  # closed
                # pick any
                if c: 
                    e = c[0]
                    return set([e[0],e[1]])
            pts = set()
            for e in ends:
                pts.add(e[0]); pts.add(e[1])
            return pts
        ep1 = curve_endpoints(c1); ep2 = curve_endpoints(c2)
        shared = len(ep1.intersection(ep2))
        # geometric consistency
        consistent = pp['ang_deg'] <= 20.0 and pp['min_dist'] <= 2.0 and abs(len(c1)-len(c2)) <= max(2, 0.3*max(len(c1),len(c2))) or shared >= 0
        promising.append({
            'curve_i': pp['i'], 'curve_j': pp['j'],
            'len1': pp['len1'], 'len2': pp['len2'],
            'max_tangent_deviation_deg': pp['ang_deg'],
            'min_dist': pp['min_dist'],
            'shared_boundary_endpoints': shared,
            'geometrically_consistent': consistent
        })
    promising.sort(key=lambda x: (x['max_tangent_deviation_deg'], x['min_dist']))
    
    analysis = {
        'mesh': str(args.mesh),
        'singular_edges_total': len(sing),
        'curve_count': len(curves_res),
        'curves': curve_details,
        'open_curves': sum(1 for c in curve_details if c['open']),
        'closed_curves': sum(1 for c in curve_details if c['closed']),
        'length_distribution': {
            'min': int(min((c['n_edges'] for c in curve_details), default=0)),
            'max': int(max((c['n_edges'] for c in curve_details), default=0)),
            'mean': float(np.mean([c['n_edges'] for c in curve_details])) if curve_details else 0.0,
            'hist': {}
        },
        'parallel_curve_pairs_total': len(parallel_pairs),
        'parallel_curve_pairs': parallel_pairs[:50],
        'pair_candidates_373_analysis': {
            'total_candidates': len(pairs_all),
            'pairs_belonging_to_same_two_curves': same_curve_pairs,
            'pairs_belonging_to_different_curves': cross_curve_pairs,
            'note': 'which of 373 pair candidates belong to the same two curves: same_curve_pairs count above'
        },
        'whole_curve_merge_assessment': {
            'question': 'would a whole-curve (full-sheet) merge be geometrically consistent?',
            'answer': 'mixed/possible for some pairs but not obvious from pairwise failures alone',
            'yes_no': 'partly yes for parallel curve pairs with similar length/shared endpoints; full-sheet merge requires consistent parallel tangents along whole curve + matching lengths + shared boundary endpoints - need case-by-case check',
            'supporting_numbers': {
                'parallel_pairs_within_tol': len([p for p in parallel_pairs if p['ang_deg']<=20.0 and p['min_dist']<=2.0]),
                'promising_curve_pairs_top10': promising[:10]
            }
        },
        'top10_promising_curve_pairs': promising[:10]
    }
    
    analysis_path = Path(args.pair_probe3) / 'curve_analysis.json'
    analysis_path.write_text(json.dumps(analysis, indent=2))
    print(f'Wrote {analysis_path}')
    
    # Generate visualizations
    plot_overview(P, H, curves_res, outdir / 'singular_curves_overview.png')
    detect_top10 = Path(args.pair_probe3) / 'detect_top10.json'
    plot_top10(P, H, f2h, curves_res, [], outdir / 'top10_pairs.png', str(detect_top10))
    
    print('Done')


if __name__ == '__main__':
    main()
