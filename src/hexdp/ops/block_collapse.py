"""Block-level sheet collapse with topological + geometric admissibility and beam search.

Input (generic):  --vtk  hex mesh with a per-cell block id (every block a structured patch, conforming)
                  --cad  CAD export of the domain (cad_export.py: face/edge samples, vertices, topology)
Steps
 1. coarse block mesh: one hex per block (block corners = coarse nodes), per-block fine grid;
    classification of every coarse node / edge / face on CAD vertex / edge / face / interior (distance to the
    CAD samples, majority vote of the fine boundary quads per block face). CAD-face surrogate = the fine
    boundary quads of the input, grouped by their CAD face.
 2. sheets = dual sheets of the coarse mesh (classes of parallel block edges).
 3. per sheet: topological collapse (weld the two ends of every sheet edge) + checks:
    sheet self-intersection / self-touching, placement (a block face on the boundary keeps its side),
    class compatibility of welded nodes / edges / faces (closure relation of the CAD topology), manifold /
    doublet checks, CAD-vertex and CAD-edge coverage.
 4. geometric realisation of the collapsed structure from the parent geometry: every welded node sits at
    fraction tau on its sheet edge (tau from the classes), merged block edges = old edge + sheet-edge part
    (union of the old curves) or the mid-curve of the sheet face, merged faces = the tau-surface of the
    sheet block, faces/blocks that absorb a half-sheet block are concatenated with it; everything is
    resampled to per-sheet interval counts (h), TFI-corrected to its new boundary and projected on the
    CAD surrogate (boundary nodes); interior by 3D TFI.  Fallback (pieces on two different axes): pure TFI.
 5. admissible = 0 inverted, min sJ >= floor, Hausdorff (CAD samples -> boundary, boundary -> CAD
    surrogate) <= baseline + tol.  Beam search (width/depth, best state anywhere).
usage: python block_collapse.py --vtk in.vtk --cad cad.npz --out DIR [--h 0.15 --sj-floor 0.05 ...]
"""
import sys, os, json, time, argparse, collections, math, itertools
import numpy as np
from multiprocessing import get_context
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bstruct import *


class Reject(Exception):
    def __init__(self, kind, msg):
        super().__init__(msg); self.kind = kind; self.msg = msg


# ================================================================================================ setup
def initial_state(P, grids, cad, tol):
    blocks = []
    for G in grids:
        n = np.array(G.shape) - 1
        c8 = tuple(int(G[tuple(o * n)]) for o in OFF)
        blocks.append(c8)
    X = [P[G] for G in grids]
    f2b = collections.defaultdict(list)
    for bi, c8 in enumerate(blocks):
        C = corner_cube(c8)
        for a, s in FACES:
            f2b[frozenset(np.take(C, s, axis=a).ravel().tolist())].append((bi, a, s))
    assert max(len(v) for v in f2b.values()) <= 2, "non-conforming input"
    # face classification: majority of the nearest CAD-face samples of the fine quad centres
    fcls, fstat, surf_tris = {}, {}, collections.defaultdict(list)
    for key, lst in f2b.items():
        if len(lst) == 2:
            fcls[key] = (3, -1); continue
        bi, a, s = lst[0]
        S = face_slice(X[bi], a, s)
        Q = np.stack([S[:-1, :-1], S[1:, :-1], S[1:, 1:], S[:-1, 1:]], 2).reshape(-1, 4, 3)
        lab, d = cad.label_points(Q.mean(1))
        cnt = collections.Counter(lab.tolist()); F, nF = cnt.most_common(1)[0]
        fcls[key] = (2, int(F)); fstat[key] = (nF / len(lab), float(d.max()))
        surf_tris[int(F)].append(Q)
    surf = {}
    for F, Qs in surf_tris.items():
        Q = np.vstack(Qs)
        surf[F] = TriSurface(np.r_[Q[:, 0], Q[:, 0]], np.r_[Q[:, 1], Q[:, 2]], np.r_[Q[:, 2], Q[:, 3]])
    if getattr(cad, "label_mode", False):
        return label_topology_state(blocks, X, fcls, fstat, surf, cad)
    # node classification
    pos = {}
    for bi, c8 in enumerate(blocks):
        n = np.array(X[bi].shape[:3]) - 1
        for k, o in enumerate(OFF):
            pos[c8[k]] = X[bi][tuple(o * n)]
    nlab = collections.defaultdict(set)
    for key, c in fcls.items():
        if c[0] == 2:
            for v in key:
                nlab[v].add(c[1])
    ncls, conflicts = {}, []
    for v, p in pos.items():
        if v not in nlab:
            ncls[v] = (3, -1); continue
        labs = nlab[v]
        dv, iv = cad.vtree.query(p)
        if dv < tol and all(cad.contains((0, int(iv)), (2, F)) for F in labs):
            ncls[v] = (0, int(iv)); continue
        R, de, ide = cad.edges.project(p[None], k=24)
        cand = None
        if de[0] < tol and labs <= cad.edge_faces[int(ide[0])]:
            cand = int(ide[0])
        if cand is not None:
            ncls[v] = (1, cand); continue
        if len(labs) == 1:
            ncls[v] = (2, next(iter(labs)))
        else:
            conflicts.append(("node", v, sorted(labs))); ncls[v] = (2, min(labs))
    # edge classification
    ecls, eline = {}, {}
    e2f = collections.defaultdict(set)
    for bi, c8 in enumerate(blocks):
        C = corner_cube(c8)
        for a, s in FACES:
            fk = frozenset(np.take(C, s, axis=a).ravel().tolist())
            if fcls[fk][0] == 2:
                F2 = np.take(C, s, axis=a)
                for e in ((F2[0, 0], F2[1, 0]), (F2[1, 0], F2[1, 1]), (F2[1, 1], F2[0, 1]), (F2[0, 1], F2[0, 0])):
                    e2f[frozenset(map(int, e))].add(fcls[fk][1])
        for ax in range(3):
            for p_ in (0, 1):
                for q_ in (0, 1):
                    L, u, v = block_edge_line(X[bi], C, ax, p_, q_)
                    eline.setdefault(frozenset((u, v)), L)
    for ek, L in eline.items():
        if ek not in e2f:
            ecls[ek] = (3, -1); continue
        labs = e2f[ek]; u, v = tuple(ek)
        done = False
        cu, cv = ncls[u], ncls[v]
        if cu[0] <= 1 and cv[0] <= 1:
            _, d, ids = cad.edges.project(L, k=24)
            E = collections.Counter(ids.tolist()).most_common(1)[0][0]
            _, dE, _ = cad.edges.project(L, only=E)
            if dE.max() < tol and labs <= cad.edge_faces[E] and cad.contains(cu, (1, E)) and cad.contains(cv, (1, E)):
                ecls[ek] = (1, int(E)); done = True
        if not done:
            if len(labs) == 1:
                ecls[ek] = (2, next(iter(labs)))
            else:
                conflicts.append(("edge", tuple(sorted(ek)), sorted(labs))); ecls[ek] = (2, min(labs))
    st = dict(c8s=blocks, X=X, ncls=ncls, ecls=ecls, fcls=fcls, path=[], origin=list(range(len(blocks))),
              next_id=int(max(pos) + 1))
    info = dict(conflicts=conflicts, face_majority_min=min(v[0] for v in fstat.values()) if fstat else 1.0,
                face_dist_max=max(v[1] for v in fstat.values()) if fstat else 0.0)
    return st, surf, info


class LabelCAD(CAD):
    """CAD given only as face samples (unconnected faces, e.g. a STEP without shared edges): the CAD edges and
    vertices are derived from the labels of the boundary block faces (curves where two CAD faces meet)."""
    def __init__(self, path):
        D = np.load(path, allow_pickle=False)
        self.fp, self.fid = D["face_pts"], D["face_id"]
        self.nf = int(self.fid.max()) + 1; self.ftree = cKDTree(self.fp); self.label_mode = True
        self.wall = set()


def label_topology_state(blocks, X, fcls, fstat, surf, cad):
    """virtual CAD edges = connected chains of block edges whose incident boundary block faces carry the same
    set of >= 2 labels; virtual vertices = chain ends, chain junctions, nodes with >= 3 labels"""
    pos, lines = {}, {}
    e2f = collections.defaultdict(set); nlab = collections.defaultdict(set)
    for bi, c8 in enumerate(blocks):
        C = corner_cube(c8); n = np.array(X[bi].shape[:3]) - 1
        for k, o in enumerate(OFF):
            pos[c8[k]] = X[bi][tuple(o * n)]
        for a, s_ in FACES:
            F2 = np.take(C, s_, axis=a); fk = frozenset(F2.ravel().tolist())
            if fcls[fk][0] == 2:
                for e in ((F2[0, 0], F2[1, 0]), (F2[1, 0], F2[1, 1]), (F2[1, 1], F2[0, 1]), (F2[0, 1], F2[0, 0])):
                    e2f[frozenset(map(int, e))].add(fcls[fk][1])
                for v in fk:
                    nlab[v].add(fcls[fk][1])
        for ax in range(3):
            for p_ in (0, 1):
                for q_ in (0, 1):
                    L, u, v = block_edge_line(X[bi], C, ax, p_, q_); lines.setdefault(frozenset((u, v)), L)
    feat = {e: frozenset(l) for e, l in e2f.items() if len(l) >= 2}
    # chains: connected components per label set
    chain = {}; nch = 0
    bylab = collections.defaultdict(list)
    for e, l in feat.items():
        bylab[l].append(e)
    edge_faces = []
    for l, es in sorted(bylab.items(), key=lambda x: sorted(x[0])):
        adj = collections.defaultdict(list)
        for e in es:
            a, b = tuple(e); adj[a].append(e); adj[b].append(e)
        for e0 in es:
            if e0 in chain:
                continue
            stack = [e0]; chain[e0] = nch
            while stack:
                e = stack.pop()
                for v in e:
                    for f in adj[v]:
                        if f not in chain:
                            chain[f] = nch; stack.append(f)
            edge_faces.append(set(l)); nch += 1
    deg = collections.defaultdict(collections.Counter)
    for e, c in chain.items():
        for v in e:
            deg[v][c] += 1
    verts = sorted(v for v, cnt in deg.items() if len(cnt) >= 2 or any(k != 2 for k in cnt.values()) or len(nlab[v]) >= 3)
    vid = {v: i for i, v in enumerate(verts)}
    cad.V = np.array([pos[v] for v in verts]) if verts else np.zeros((0, 3)); cad.vtree = cKDTree(cad.V) if verts else None
    cad.edge_faces = edge_faces
    cad.vert_edges = [set(deg[v].keys()) for v in verts]
    cad.edge_verts = collections.defaultdict(set)
    for i, es in enumerate(cad.vert_edges):
        for e in es:
            cad.edge_verts[e].add(i)
    Ps, ids = [], []
    for k, (e, c) in enumerate(chain.items()):
        L = lines[e]; Ps.append(L); ids.append(np.full(len(L), k))
    cad.edges = Polylines(np.vstack(Ps), np.concatenate(ids))
    remap = np.array([c for e, c in chain.items()]); cad.edges.id = remap[cad.edges.id]
    cad.ep = np.vstack(Ps); cad.eid = remap[np.concatenate(ids)]
    ncls, ecls, conflicts = {}, {}, []
    for v in pos:
        if v in vid:
            ncls[v] = (0, vid[v])
        elif v in deg:
            ncls[v] = (1, next(iter(deg[v])))
        elif v in nlab:
            if len(nlab[v]) > 1:
                conflicts.append(("node", v, sorted(nlab[v])))
            ncls[v] = (2, min(nlab[v]))
        else:
            ncls[v] = (3, -1)
    for e in lines:
        if e in chain:
            ecls[e] = (1, chain[e])
        elif e in e2f:
            ecls[e] = (2, next(iter(e2f[e])))
        else:
            ecls[e] = (3, -1)
    st = dict(c8s=blocks, X=X, ncls=ncls, ecls=ecls, fcls=fcls, path=[], origin=list(range(len(blocks))),
              next_id=int(max(pos) + 1))
    info = dict(conflicts=conflicts, face_majority_min=min(v[0] for v in fstat.values()) if fstat else 1.0,
                face_dist_max=max(v[1] for v in fstat.values()) if fstat else 0.0,
                virtual=dict(edges=nch, vertices=len(verts)))
    return st, surf, info


def wall_first_cells(Xs, c8s, fcls, wall):
    """first wall-normal spacing at every node of every block face on a wall CAD face"""
    out = []
    for X, c8 in zip(Xs, c8s):
        C = corner_cube(c8)
        for a, s_ in FACES:
            c = fcls[frozenset(np.take(C, s_, axis=a).ravel().tolist())]
            if c[0] == 2 and c[1] in wall:
                i0, i1 = (0, 1) if s_ == 0 else (-1, -2)
                out.append(np.linalg.norm(np.take(X, i1, axis=a) - np.take(X, i0, axis=a), axis=-1).ravel())
    return np.concatenate(out) if out else np.zeros(0)


def block_edge_line(X, C, ax, p, q):
    """fine polyline of the block edge along axis ax at the corners (p, q) of the other two axes; oriented
    from the smaller node id. Returns (line, u, v) with u < v."""
    o = [a for a in range(3) if a != ax]
    idx = [None] * 3; idx[ax] = slice(None); idx[o[0]] = 0 if p == 0 else -1; idx[o[1]] = 0 if q == 0 else -1
    ci = [0] * 3; ci[o[0]] = p; ci[o[1]] = q
    L = X[tuple(idx)]
    a_ = list(ci); a_[ax] = 0; b_ = list(ci); b_[ax] = 1
    u, v = int(C[tuple(a_)]), int(C[tuple(b_)])
    if u > v:
        return L[::-1].copy(), v, u
    return L.copy(), u, v


def coverage(st, cad):
    """CAD vertex ids represented by nodes, and per CAD edge the coarse edges classified on it"""
    verts = {c[1] for v, c in st["ncls"].items() if c[0] == 0 and v in used_nodes(st)}
    ce = collections.defaultdict(set)
    for ek in used_edges(st):
        c = st["ecls"][ek]
        if c[0] == 1:
            ce[c[1]].add(ek)
    return verts, ce


def used_nodes(st):
    return {v for c8 in st["c8s"] for v in c8}


def used_edges(st):
    out = set()
    for c8 in st["c8s"]:
        for grp in PARALLEL:
            for i, j in grp:
                out.add(frozenset((c8[i], c8[j])))
    return out


# ================================================================================================ topology
def sheet_orientation(c8s, mem):
    """consistent '-' side per sheet edge (dict edge -> node on the '-' side) or None if the sheet is one-sided"""
    bl = {b: ax for b, ax in mem}
    side = {}
    def dirs(b, flip):
        C = corner_cube(c8s[b]); ax = bl[b]
        A = np.take(C, 0, axis=ax).ravel(); B = np.take(C, 1, axis=ax).ravel()
        if flip:
            A, B = B, A
        return {frozenset((int(x), int(y))): int(x) for x, y in zip(A, B)}
    e2b = collections.defaultdict(list)
    for b in bl:
        for e in dirs(b, False):
            e2b[e].append(b)
    out = {}
    for b0 in bl:
        if b0 in side:
            continue
        side[b0] = False; stack = [b0]
        while stack:
            b = stack.pop()
            for e, lo in dirs(b, side[b]).items():
                if out.get(e, lo) != lo:
                    return None
                out[e] = lo
                for g in e2b[e]:
                    if g in side:
                        continue
                    dg = dirs(g, False)
                    side[g] = dg[e] != lo; stack.append(g)
    for b in bl:
        for e, lo in dirs(b, side[b]).items():
            if out[e] != lo:
                return None
    return out


def collapse_topo(st, sheet, cad, ref_cov, variant="mid"):
    c8s, ncls, ecls, fcls = st["c8s"], st["ncls"], st["ecls"], st["fcls"]
    E = sheet["edges"]; mem = sheet["members"]
    axes = collections.defaultdict(set)
    for b, ax in mem:
        axes[b].add(ax)
    if any(len(v) > 1 for v in axes.values()):
        raise Reject("topology", "self-intersecting sheet (a block is crossed twice)")
    cnt = collections.Counter(v for e in E for v in e)
    if cnt and max(cnt.values()) > 1:
        raise Reject("topology", "sheet touches itself (a node on two sheet edges)")
    stay = {}
    for b, ax in mem:
        C = corner_cube(c8s[b]); A = np.take(C, 0, axis=ax); B = np.take(C, 1, axis=ax)
        bA = fcls[frozenset(A.ravel().tolist())][0] == 2; bB = fcls[frozenset(B.ravel().tolist())][0] == 2
        wall = getattr(cad, "wall", set())
        for F_ in (A, B):
            c_ = fcls[frozenset(F_.ravel().tolist())]
            if c_[0] == 2 and c_[1] in wall:
                raise Reject("wall-layer", f"sheet is wall-normal at wall face {c_[1]} (would remove / regrade the wall layer)")
        if bA and bB:
            raise Reject("topology", "sheet block between two boundary faces (region would vanish)")
        if bA or bB:
            S = A if bA else B
            for x, y in zip(A.ravel(), B.ravel()):
                e = frozenset((int(x), int(y))); w = int(x) if bA else int(y)
                if stay.get(e, w) != w:
                    raise Reject("topology", "placement conflict (sheet edge forced to both ends)")
                stay[e] = w
    orient_ = sheet_orientation(c8s, mem) if variant != "mid" else None
    if variant != "mid" and orient_ is None:
        raise Reject("variant", "one-sided sheet: no lo/hi placement")
    nid = st["next_id"]; nodemap = {}; merged = {}; ncls2 = {}; nfree = 0
    for e in sorted(E, key=lambda e: tuple(sorted(e))):
        u, v = sorted(e); cu, cv = ncls[u], ncls[v]
        if e in stay:
            w = stay[e]; o = v if w == u else u
            c = ncls[w]
            if not cad.contains(c, ncls[o]):
                raise Reject("topology", f"boundary side forced, but welded node classes {ncls[w]}/{ncls[o]} incompatible")
            tau = 0.0 if w == u else 1.0
        else:
            c = merge_class(cad, cu, cv)
            if c is None:
                raise Reject("topology", f"welds nodes on different CAD entities {cu} / {cv}")
            tau = 0.0 if cu[0] < cv[0] else (1.0 if cv[0] < cu[0] else 0.5)
            if tau == 0.5 and variant != "mid":
                lo = orient_[e]; w = lo if variant == "lo" else (v if lo == u else u)
                tau = 0.0 if w == u else 1.0
            nfree += 1 if cu[0] == cv[0] else 0
        nodemap[u] = nodemap[v] = nid; merged[nid] = (u, v, tau); ncls2[nid] = c; nid += 1
    mp = lambda x: nodemap.get(x, x)
    sheet_blocks = set(axes)
    surv = [b for b in range(len(c8s)) if b not in sheet_blocks]
    if not surv:
        raise Reject("topology", "no block left")
    new_c8 = [tuple(mp(x) for x in c8s[b]) for b in surv]
    for c8 in new_c8:
        if len(set(c8)) != 8:
            raise Reject("topology", "degenerate block")
    if len({frozenset(c) for c in new_c8}) != len(new_c8):
        raise Reject("topology", "two blocks with the same corners")
    for v in used_nodes(dict(c8s=new_c8)):
        if v not in ncls2:
            ncls2[v] = ncls[v]
    # edges
    ecls2, epre = {}, collections.defaultdict(list)
    for ek, c in ecls.items():
        if ek in E:
            continue
        img = frozenset(mp(x) for x in ek)
        if len(img) < 2:
            continue
        epre[img].append(ek)
        if img in ecls2:
            m = merge_class(cad, ecls2[img], c)
            if m is None:
                raise Reject("topology", f"merges block edges on different CAD entities {ecls2[img]} / {c}")
            ecls2[img] = m
        else:
            ecls2[img] = c
    # faces
    fcls2, fpre = {}, collections.defaultdict(list)
    for fk, c in fcls.items():
        img = frozenset(mp(x) for x in fk)
        if len(img) < 4:
            continue
        fpre[img].append(fk)
        if img in fcls2:
            m = merge_class(cad, fcls2[img], c)
            if m is None:
                raise Reject("topology", f"merges block faces on different CAD faces {fcls2[img]} / {c}")
            fcls2[img] = m
        else:
            fcls2[img] = c
    # new structure checks
    f2b = collections.defaultdict(list); bfaces = []
    for bi, c8 in enumerate(new_c8):
        C = corner_cube(c8); fs = []
        for a, s in FACES:
            k = frozenset(np.take(C, s, axis=a).ravel().tolist()); f2b[k].append(bi); fs.append(k)
        bfaces.append(fs)
    for k, bl in f2b.items():
        if k not in fcls2:
            raise Reject("topology", "new face without preimage")
        if len(bl) > 2:
            raise Reject("topology", "non-manifold face (3 blocks)")
        if len(bl) == 2 and fcls2[k][0] == 2:
            raise Reject("topology", f"boundary faces glued together (CAD face {fcls2[k][1]} would become interior)")
        if len(bl) == 2 and bl[0] == bl[1]:
            raise Reject("topology", "block glued to itself")
        if len(bl) == 1 and fcls2[k][0] == 3:
            raise Reject("topology", "interior face exposed")
    pairs = collections.Counter()
    for k, bl in f2b.items():
        if len(bl) == 2:
            pairs[tuple(sorted(bl))] += 1
    if pairs and max(pairs.values()) > 1:
        raise Reject("topology", "doublet (two blocks share two faces)")
    be = bnd_edge_count(new_c8)
    be0 = bnd_edge_count(c8s)
    multi0 = {frozenset(mp(x) for x in e) for e, n in be0.items() if n > 2}
    if any(n % 2 for n in be.values()):
        raise Reject("topology", "boundary surface has an edge with an odd number of faces")
    if any(n > 2 and e not in multi0 for e, n in be.items()):
        raise Reject("topology", "new non-manifold boundary edge (boundary pinched / glued along an edge)")
    # edge fans: blocks around a block edge must stay one connected fan (unless already split in the parent,
    # e.g. a slit attached to a wall)
    split0 = {frozenset(mp(x) for x in e) for e in fan_split_edges(c8s)}
    bad = [e for e in fan_split_edges(new_c8) if e not in split0]
    if bad:
        raise Reject("topology", "non-manifold block edge (blocks around an edge no longer one fan)")
    child = dict(c8s=new_c8, ncls={v: ncls2[v] for v in used_nodes(dict(c8s=new_c8))},
                 ecls={k: ecls2[k] for k in used_edges(dict(c8s=new_c8))}, fcls={k: fcls2[k] for k in f2b},
                 origin=[st["origin"][b] for b in surv], next_id=nid)
    # CAD coverage
    verts0, ce0 = ref_cov
    verts1, ce1 = coverage(child, cad)
    if verts0 - verts1:
        raise Reject("topology", f"CAD vertex lost {sorted(verts0 - verts1)}")
    for Eid, es in ce0.items():
        es1 = ce1.get(Eid, set())
        if not es1:
            raise Reject("topology", f"CAD edge {Eid} no longer represented by block edges")
        deg = collections.Counter(v for e in es1 for v in e)
        adj = collections.defaultdict(set)
        for e in es1:
            a, b = tuple(e); adj[a].add(b); adj[b].add(a)
        start = next(iter(deg)); seen = {start}; stack = [start]
        while stack:
            x = stack.pop()
            for y in adj[x]:
                if y not in seen:
                    seen.add(y); stack.append(y)
        if len(seen) != len(deg) or max(deg.values()) > 2:
            raise Reject("topology", f"CAD edge {Eid}: block edges no longer a single chain")
    warn = 0
    for bi, c8 in enumerate(new_c8):
        C = corner_cube(c8)
        for (a1, s1), (a2, s2) in itertools.combinations(FACES, 2):
            if a1 == a2:
                continue
            k1 = frozenset(np.take(C, s1, axis=a1).ravel().tolist()); k2 = frozenset(np.take(C, s2, axis=a2).ravel().tolist())
            if child["fcls"][k1][0] == 2 and child["fcls"][k1] == child["fcls"][k2]:
                warn += 1
    if variant != "mid" and nfree == 0:
        raise Reject("variant", "no free sheet edge: variant identical to mid")
    geo = dict(nodemap=nodemap, merged=merged, surv=surv, mem=mem, E=E, epre=epre, fpre=fpre, flat_pairs=warn)
    return child, geo


def fan_split_edges(c8s):
    f2b = collections.defaultdict(list); e2b = collections.defaultdict(set)
    for bi, c8 in enumerate(c8s):
        C = corner_cube(c8)
        for a, s in FACES:
            f2b[frozenset(np.take(C, s, axis=a).ravel().tolist())].append(bi)
        for grp in PARALLEL:
            for i, j in grp:
                e2b[frozenset((c8[i], c8[j]))].add(bi)
    adj = collections.defaultdict(list)
    for k, bl in f2b.items():
        if len(bl) == 2:
            adj[bl[0]].append((bl[1], k)); adj[bl[1]].append((bl[0], k))
    out = []
    for ek, bl in e2b.items():
        if len(bl) < 2:
            continue
        start = next(iter(bl)); seen = {start}; stack = [start]
        while stack:
            x = stack.pop()
            for y, k in adj[x]:
                if ek <= k and y not in seen:
                    seen.add(y); stack.append(y)
        if len(seen) != len(bl):
            out.append(ek)
    return out


def bnd_edge_count(c8s):
    f2b = collections.Counter(); fc = {}
    for c8 in c8s:
        C = corner_cube(c8)
        for a, s in FACES:
            F2 = np.take(C, s, axis=a); k = frozenset(F2.ravel().tolist()); f2b[k] += 1; fc[k] = F2
    be = collections.Counter()
    for k, n in f2b.items():
        if n == 1:
            C = fc[k]
            for e in ((C[0, 0], C[1, 0]), (C[1, 0], C[1, 1]), (C[1, 1], C[0, 1]), (C[0, 1], C[0, 0])):
                be[frozenset(map(int, e))] += 1
    return be


# ================================================================================================ geometry
class ParentGeo:
    def __init__(self, st):
        self.st = st; self.X = st["X"]; self.C = [corner_cube(c) for c in st["c8s"]]
        self.lines = {}; self.faces = {}; self.pos = {}
        for bi in range(len(self.X)):
            X, C = self.X[bi], self.C[bi]
            n = np.array(X.shape[:3]) - 1
            for k, o in enumerate(OFF):
                self.pos[int(C[tuple(o)])] = X[tuple(o * n)]
            for ax in range(3):
                for p in (0, 1):
                    for q in (0, 1):
                        L, u, v = block_edge_line(X, C, ax, p, q)
                        self.lines.setdefault(frozenset((u, v)), (L, u))
            for a, s in FACES:
                k = frozenset(np.take(C, s, axis=a).ravel().tolist())
                self.faces.setdefault(k, (face_slice(X, a, s), np.take(C, s, axis=a)))

    def line(self, x, y):
        L, u = self.lines[frozenset((x, y))]
        return L if u == x else L[::-1]

    def face(self, target):
        G, C = self.faces[frozenset(int(t) for t in target.ravel() if t >= 0)] if (target >= 0).all() else (None, None)
        return orient(G, C, target)


def subpoly(L, f):
    """polyline from L[0] to fractional index f*(n) (f in [0,1]) -> (points, index parameter)"""
    n = len(L) - 1; x = f * n; i = int(math.floor(x))
    if i >= n:
        return L.copy(), np.arange(n + 1, dtype=float)
    w = x - i
    P = np.vstack([L[:i + 1], (1 - w) * L[i] + w * L[i + 1]]); t = np.r_[np.arange(i + 1, dtype=float), x]
    return P, t


def rev(Pt):
    P, t = Pt
    return P[::-1].copy(), t[-1] - t[::-1]


def concat(*parts):
    """concatenate (points, t) pieces sharing end points; t continues (index units of the old fine cells)"""
    P, t = parts[0]
    for Q, u in parts[1:]:
        if len(Q) > 1:
            P = np.vstack([P, Q[1:]]); t = np.r_[t, t[-1] + u[1:] - u[0]]
    keep = np.r_[True, np.diff(t) > 1e-9]
    if keep.sum() < 2:
        return P[[0, -1]], np.array([0.0, max(t[-1], 1e-9)])
    return P[keep], t[keep]


def resample_t(P, t, n):
    tt = np.linspace(t[0], t[-1], n + 1)
    return np.stack([np.interp(tt, t, P[:, k]) for k in range(3)], 1)


def realise(st, child, geo, cad, surf, h, nmin=2):
    """fine geometry of the collapsed structure `child` from the parent state `st`; geo = None -> identity
    (baseline re-fit). Returns list of block grids."""
    PG = ParentGeo(st)
    nodemap = geo["nodemap"] if geo else {}; merged = geo["merged"] if geo else {}
    mp = lambda x: nodemap.get(x, x)
    frac = {}                                         # old node -> fraction from it along its sheet edge
    sheet_edge_of = {}
    for m, (u, v, tau) in merged.items():
        frac[u] = tau; frac[v] = 1 - tau; sheet_edge_of[u] = v; sheet_edge_of[v] = u
    # ---- node positions
    pos = {}
    for v in used_nodes(child):
        if v in merged:
            u, w, tau = merged[v]; pos[v] = subpoly(PG.line(u, w), tau)[0][-1]
        else:
            pos[v] = PG.pos[v]
    for v, p in list(pos.items()):
        c = child["ncls"][v]
        if c[0] == 0:
            pos[v] = cad.V[c[1]].copy()
        elif c[0] == 1:
            pos[v] = cad.edges.project(p[None], only=c[1])[0][0]
        elif c[0] == 2:
            pos[v] = surf[c[1]].project(p[None])[0][0]
    # ---- edge curves (raw)
    if geo:
        epre = geo["epre"]
    else:
        epre = {k: [k] for k in child["ecls"]}
    raw = {}; hloc = {}
    for ek in used_edges(child):
        x_, y_ = sorted(ek); pre = epre[ek]
        hl = []
        for p_ in pre:
            Lp = PG.lines[frozenset(p_)][0]
            hl.append(np.linalg.norm(np.diff(Lp, axis=0), axis=1).sum() / max(len(Lp) - 1, 1))
        hloc[ek] = min(hl) if hl else h
        if len(pre) > 2:
            raise Reject("geometry", "block edge with > 2 preimages")
        L = None
        if len(pre) == 2:
            a, b = (tuple(sorted(p)) for p in pre)
            # find y1 (pre a) with mp == x_
            ya = a if mp(a[0]) == x_ else a[::-1]; yb = b if mp(b[0]) == x_ else b[::-1]
            T = np.array([[ya[0], yb[0]], [ya[1], yb[1]]], np.int64)
            k = frozenset(T.ravel().tolist())
            if k in PG.faces:
                G, _ = PG.face(T)
                f = np.linspace(frac[ya[0]], frac[ya[1]], G.shape[0])
                L = sample_surface_frac(G[:, None], f[:, None])[:, 0]; L = (L, np.arange(len(L), dtype=float))
            else:
                pre = [a]
        if L is None:
            a = tuple(sorted(pre[0])); a = a if mp(a[0]) == x_ else a[::-1]
            L = PG.line(a[0], a[1])
            parts = []
            if a[0] in frac:
                parts.append(rev(subpoly(PG.line(a[0], sheet_edge_of[a[0]]), frac[a[0]])))
            parts.append((L, np.arange(len(L), dtype=float)))
            if a[1] in frac:
                parts.append(subpoly(PG.line(a[1], sheet_edge_of[a[1]]), frac[a[1]]))
            L = concat(*parts)
        raw[ek] = L
    # ---- interval counts per sheet of the child
    sh = sheets_of(child["c8s"]); ecount = {}
    for s in sh:
        if h > 0:
            Lmax = max(np.linalg.norm(np.diff(raw[e][0], axis=0), axis=1).sum() for e in s["edges"])
            n = max(nmin, int(math.ceil(Lmax / h)))
        else:                                   # keep the parent's cells: count = old fine cells along the edge
            n = max(nmin, max(int(math.ceil(raw[e][1][-1] - 1e-6)) for e in s["edges"]))
        for e in s["edges"]:
            ecount[e] = n
    curves = {}
    for ek, L in raw.items():
        x_, y_ = sorted(ek)
        R = resample_axis(L[0], 0, ecount[ek]) if h > 0 else resample_t(L[0], L[1], ecount[ek])
        R[0] = pos[x_]; R[-1] = pos[y_]
        c = child["ecls"][ek]
        if c[0] == 1 and len(R) > 2:
            R[1:-1] = cad.edges.project(R[1:-1], only=c[1])[0]
        elif c[0] == 2 and len(R) > 2:
            R[1:-1] = surf[c[1]].project(R[1:-1])[0]
        curves[ek] = R

    def curve(x, y):
        R = curves[frozenset((x, y))]
        return R if x < y else R[::-1]
    # ---- faces
    fpre = geo["fpre"] if geo else {k: [k] for k in child["fcls"]}
    mem = dict(geo["mem"]) if geo else {}
    old_face_owner = {}                     # old face key -> (sheet block, axis, side)
    for b, ax in mem.items():
        C = PG.C[b]
        for s in (0, 1):
            old_face_owner[frozenset(np.take(C, s, axis=ax).ravel().tolist())] = (b, ax, s)

    def tau_field(b, ax, side, shape2):
        """fraction (from `side`) of the welded position over the cross-section of sheet block b"""
        C = PG.C[b]; A = np.take(C, side, axis=ax)
        Tc = np.array([[frac[int(A[i, j])] for j in (0, 1)] for i in (0, 1)])
        u = np.linspace(0, 1, shape2[0])[:, None]; v = np.linspace(0, 1, shape2[1])[None]
        return (1 - u) * (1 - v) * Tc[0, 0] + u * (1 - v) * Tc[1, 0] + (1 - u) * v * Tc[0, 1] + u * v * Tc[1, 1]

    def sheet_volume(b, ax, side):
        H = np.moveaxis(PG.X[b], ax, 2)
        if side == 1:
            H = H[:, :, ::-1]
        return H                                   # (nu+1, nv+1, nw+1, 3), w from `side`

    face_grid = {}
    for fk in child["fcls"]:
        pre = fpre[fk]
        prior, Cn = None, None
        if len(pre) == 2:
            b, ax, s = old_face_owner[pre[0]]
            H = sheet_volume(b, ax, s); T = tau_field(b, ax, s, H.shape[:2])
            prior = sample_surface_frac(H, T); Cn = np.vectorize(mp)(np.take(PG.C[b], s, axis=ax))
        else:
            G, C2 = PG.faces[pre[0]]
            strips = []
            for a in (0, 1):
                for s in (0, 1):
                    y = np.take(C2, s, axis=a)       # 2 ids along the other axis
                    if not (int(y[0]) in frac and int(y[1]) in frac):
                        continue
                    z = [sheet_edge_of[int(y[0])], sheet_edge_of[int(y[1])]]
                    T = np.array([[y[0], z[0]], [y[1], z[1]]], np.int64)
                    if frozenset(T.ravel().tolist()) not in PG.faces:
                        continue
                    fm = max(frac[int(y[0])], frac[int(y[1])])
                    if fm <= 0:
                        continue
                    Phi, _ = PG.face(T)                   # axis0 along y, axis1 along the sheet edge
                    nk = max(1, int(math.ceil((Phi.shape[1] - 1) * fm)))
                    f = np.linspace(frac[int(y[0])], frac[int(y[1])], Phi.shape[0])
                    Sx = np.stack([sample_surface_frac(Phi[:, None], (f * k / nk)[:, None])[:, 0] for k in range(nk + 1)], 1)
                    strips.append((a, s, Sx))
            if len({a for a, s, _ in strips}) > 1:
                prior = None
            else:
                prior = G.copy()
                for a, s, Sx in strips:
                    # Sx: (len along other axis, nk+1) ; resample to G's count along the other axis
                    o = 1 - a
                    Sx = resample_axis(Sx, 0, G.shape[o] - 1) if Sx.shape[0] != G.shape[o] else Sx
                    if a == 0:
                        Sx = np.swapaxes(Sx, 0, 1)        # (nk+1, n_other, 3)
                        prior = np.concatenate([prior, Sx[1:]], 0) if s == 1 else np.concatenate([Sx[::-1][:-1], prior], 0)
                    else:
                        prior = np.concatenate([prior, Sx[:, 1:]], 1) if s == 1 else np.concatenate([Sx[:, ::-1][:, :-1], prior], 1)
            Cn = np.vectorize(mp)(C2)
        Cn = np.asarray(Cn, np.int64)
        n0 = ecount[frozenset((int(Cn[0, 0]), int(Cn[1, 0])))]; n1 = ecount[frozenset((int(Cn[0, 0]), int(Cn[0, 1])))]
        Eb = {(0, 0): curve(Cn[0, 0], Cn[0, 1]), (0, 1): curve(Cn[1, 0], Cn[1, 1]),
              (1, 0): curve(Cn[0, 0], Cn[1, 0]), (1, 1): curve(Cn[0, 1], Cn[1, 1])}
        if prior is None:
            G2 = coons2(Eb[(0, 0)], Eb[(0, 1)], Eb[(1, 0)], Eb[(1, 1)])
        else:
            G2 = resample_face_param(prior, Eb)
        c = child["fcls"][fk]
        if c[0] == 2 and min(G2.shape[:2]) > 2:
            G2[1:-1, 1:-1] = surf[c[1]].project(G2[1:-1, 1:-1])[0]
        face_grid[fk] = canon2(G2, Cn)
    # ---- blocks
    Xn = []; nfallback = 0
    surv = geo["surv"] if geo else list(range(len(child["c8s"])))
    for nb, ob in enumerate(surv):
        Cold = PG.C[ob]; Cnew = corner_cube(child["c8s"][nb])
        prior = PG.X[ob].copy(); pieces = []
        for a, s in FACES:
            k = frozenset(np.take(Cold, s, axis=a).ravel().tolist())
            if k not in old_face_owner:
                continue
            b, bx, bs = old_face_owner[k]
            H = sheet_volume(b, bx, bs); T = tau_field(b, bx, bs, H.shape[:2])
            if T.max() <= 0:
                continue
            nk = max(1, int(math.ceil((H.shape[2] - 1) * T.max())))
            Pc = np.stack([sample_surface_frac(H, T * kk / nk) for kk in range(nk + 1)], 2)   # w from the face
            Cb = np.moveaxis(PG.C[b], bx, 2)
            if bs == 1:
                Cb = Cb[:, :, ::-1]
            Cp = Cb.copy(); Cp[:, :, 1] = -7                                              # far side: new ids
            tgt = -np.ones((2, 2, 2), np.int64)
            sl = [slice(None)] * 3; sl[a] = 0 if s == 1 else 1
            tgt[tuple(sl)] = np.take(Cold, s, axis=a)
            r = orient(Pc, Cp, tgt)
            if r is None:
                raise Reject("geometry", "piece orientation failed")
            pieces.append((a, s, r[0]))
        if len({a for a, s, _ in pieces}) > 1:
            prior = None; nfallback += 1
        else:
            for a, s, Pc in pieces:
                cnt = [prior.shape[i] for i in range(3)]
                for i in range(3):
                    if i != a and Pc.shape[i] != cnt[i]:
                        Pc = resample_axis(Pc, i, cnt[i] - 1)
                if s == 1:
                    prior = np.concatenate([prior, np.take(Pc, range(1, Pc.shape[a]), axis=a)], a)
                else:
                    prior = np.concatenate([np.take(Pc, range(0, Pc.shape[a] - 1), axis=a), prior], a)
        cnts = [ecount[frozenset((int(Cnew[0, 0, 0]), int(Cnew[tuple(np.eye(3, dtype=int)[i])])))] for i in range(3)]
        F = {}
        for a, s in FACES:
            tgt = np.take(Cnew, s, axis=a)
            G, Cf = face_grid[frozenset(tgt.ravel().tolist())]
            r = orient(G, Cf, tgt); F[(a, s)] = r[0]
        if prior is None:
            Xb = tfi_correct3(np.zeros(tuple(n + 1 for n in cnts) + (3,)), F)
        else:
            el = {}
            for ax in range(3):
                o = [q for q in range(3) if q != ax]; el[ax] = {}
                for p_ in (0, 1):
                    for q_ in (0, 1):
                        i0 = [0, 0, 0]; i0[o[0]] = p_; i0[o[1]] = q_; i1 = list(i0); i1[ax] = 1
                        el[ax][(p_, q_)] = curve(int(Cnew[tuple(i0)]), int(Cnew[tuple(i1)]))
            Xb = resample_block_param(prior, el, F)
        Xn.append(Xb)
    return Xn, dict(fallback_blocks=nfallback, face_grid=face_grid, cells=int(sum(np.prod(np.array(x.shape[:3]) - 1) for x in Xn)))


# ================================================================================================ evaluation
def evaluate(Xn, child, info, cad, surf_all, rev_pts, rev_base=None):
    sj = np.concatenate([hex_sj(X).ravel() for X in Xn])
    bmin = [float(hex_sj(X).min()) for X in Xn]
    Qs = []
    for fk, (G, C) in info["face_grid"].items():
        if child["fcls"][fk][0] == 2:
            Qs.append(np.stack([G[:-1, :-1], G[1:, :-1], G[1:, 1:], G[:-1, 1:]], 2).reshape(-1, 4, 3))
    Q = np.vstack(Qs)
    S = TriSurface(np.r_[Q[:, 0], Q[:, 0]], np.r_[Q[:, 1], Q[:, 2]], np.r_[Q[:, 2], Q[:, 3]], k=16)
    _, drev = S.project(rev_pts)
    _, dfwd = surf_all.project(Q.mean(1))
    extra = {}
    if rev_base is not None:
        extra["dh_rev"] = float((drev - rev_base).max())
    else:
        extra["_drev"] = drev
    wall = getattr(cad, "wall", set())
    if wall:
        fc = wall_first_cells(Xn, child["c8s"], child["fcls"], wall)
        extra.update(wall_h1_min=float(fc.min()), wall_h1_max=float(fc.max()), wall_nodes=int(fc.size))
    return dict(**extra, blocks=len(Xn), cells=int(sj.size), inverted=int((sj <= 0).sum()), min_sj=float(sj.min()),
                p1_sj=float(np.percentile(sj, 1)), median_sj=float(np.median(sj)),
                h_rev=float(drev.max()), h_rev_p99=float(np.percentile(drev, 99)), h_fwd=float(dfwd.max()),
                worst_blocks=sorted(zip(bmin, child["origin"]))[:4], fallback_blocks=info["fallback_blocks"],
                inverted_in=sorted(((int((hex_sj(X) <= 0).sum()), child["origin"][b]) for b, X in enumerate(Xn) if bmin[b] <= 0), reverse=True)[:6])


# ================================================================================================ beam
G_ = {}


def _eval_child(job):
    i, var = job
    st, sheets, cad, surf, surf_all, rev, h, ref_cov = (G_[k] for k in ("st", "sheets", "cad", "surf", "surf_all", "rev", "h", "ref_cov"))
    t0 = time.time()
    try:
        child, geo = collapse_topo(st, sheets[i], cad, ref_cov, var)
        Xn, inf = realise(st, child, geo, cad, surf, h)
        ev = evaluate(Xn, child, inf, cad, surf_all, rev, G_.get("rev_base"))
        ev["flat_pairs"] = geo["flat_pairs"]
        child["X"] = Xn
        return i, var, "ok", ev, child, time.time() - t0
    except Reject as e:
        return i, var, e.kind, e.msg, None, time.time() - t0


def describe_sheet(st, s, meta):
    """reporting only: block kinds / 2D ids of the members (from the VTK cell data) and the sheet edge
    direction relative to the z axis / azimuth"""
    orig = [st["origin"][b] for b, _ in s["members"]]
    kinds = collections.Counter(meta["kind"][o] for o in orig) if meta.get("kind") is not None else {}
    d_az = []
    PG = None
    for e in list(s["edges"])[:64]:
        u, v = sorted(e)
        pu, pv = st["_pos"][u], st["_pos"][v]
        m = (pu + pv) / 2; r = np.hypot(m[0], m[1])
        az = np.array([-m[1], m[0], 0]) / max(r, 1e-9); dv = pv - pu
        d_az.append(abs(np.dot(dv, az)) / max(np.linalg.norm(dv), 1e-12))
    return dict(n=len(s["members"]), kinds=dict(kinds), azimuthal_frac=float(np.mean(d_az)),
                blocks2d=sorted({int(meta["block2d"][o]) for o in orig}) if meta.get("block2d") is not None else [],
                origin=sorted(orig))


def node_pos(st):
    pos = {}
    for X, c8 in zip(st["X"], st["c8s"]):
        n = np.array(X.shape[:3]) - 1
        for k, o in enumerate(OFF):
            pos[c8[k]] = X[tuple(o * n)]
    return pos


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vtk", required=True); ap.add_argument("--cad", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--field", default="block")
    ap.add_argument("--h", type=float, default=0.0, help="target cell size of the TFI fill; <= 0: keep the input "
                    "cell size along every block edge (interval count per sheet = max over its edges)")
    ap.add_argument("--h-final", type=float, default=0.1)
    ap.add_argument("--tol", type=float, default=0.01, help="classification distance tolerance")
    ap.add_argument("--sj-floor", default="0.05", help="absolute floor, or 'input' = not below the input min sJ")
    ap.add_argument("--cad-topology", choices=("occ", "labels"), default="occ",
                    help="occ: CAD edges/vertices from cad_export.py; labels: --cad holds face samples only, edges and "
                    "vertices are derived from the block-face labels (CAD without shared topology)")
    ap.add_argument("--wall-faces", default="", help="CAD face ids carrying a graded wall layer: a sheet that is "
                    "wall-normal at such a face is not collapsed, and the first wall-normal cell must not change")
    ap.add_argument("--h1-tol", type=float, default=0.01, help="relative tolerance on the first wall cell")
    ap.add_argument("--rev-relative", action="store_true", help="reverse Hausdorff guard per CAD sample relative to "
                    "the baseline (for CAD samples the mesh never reaches, e.g. plane parts inside walls)")
    ap.add_argument("--htol", type=float, default=0.01, help="allowed Hausdorff increase over the baseline re-fit")
    ap.add_argument("--width", type=int, default=4); ap.add_argument("--depth", type=int, default=8)
    ap.add_argument("--workers", type=int, default=8); ap.add_argument("--nrev", type=int, default=40000)
    ap.add_argument("--meta-fields", default="kind,block2d", help="cell data used for reporting only")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    T0 = time.time()
    log = open(f"{a.out}/block_collapse.log", "w")

    def say(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); log.write(s + "\n"); log.flush()
    cad = LabelCAD(a.cad) if a.cad_topology == "labels" else CAD(a.cad)
    cad.wall = {int(x) for x in a.wall_faces.split(",") if x.strip()}
    P, grids = read_blocks_vtk(a.vtk, a.field)
    meta = {}
    try:
        import meshio
        m = meshio.read(a.vtk)
        bid = np.concatenate([d for d, c in zip(m.cell_data[a.field], m.cells) if c.type == "hexahedron"]).astype(int)
        ub = np.unique(bid)
        for f in a.meta_fields.split(","):
            if f and f in m.cell_data:
                v = np.concatenate([d for d, c in zip(m.cell_data[f], m.cells) if c.type == "hexahedron"])
                v = np.asarray(v).ravel()
                meta[f] = {i: int(v[np.where(bid == b)[0][0]]) for i, b in enumerate(ub)}
    except Exception as e:
        say("meta fields not read:", e)
    st0, surf, cinfo = initial_state(P, grids, cad, a.tol)
    say(f"[input] {len(st0['c8s'])} blocks, {len(used_nodes(st0))} coarse nodes, {len(used_edges(st0))} block edges, "
        f"{len(st0['fcls'])} block faces ({sum(c[0] == 2 for c in st0['fcls'].values())} on CAD faces)")
    say(f"[class] nodes by dim {dict(collections.Counter(c[0] for c in st0['ncls'].values()))}, "
        f"edges by dim {dict(collections.Counter(c[0] for c in st0['ecls'].values()))}; "
        f"face-label majority min {cinfo['face_majority_min']:.3f}, max quad-centre->CAD-sample dist {cinfo['face_dist_max']:.4f}; "
        f"conflicts {len(cinfo['conflicts'])} {cinfo['conflicts'][:6]}")
    ref_cov = coverage(st0, cad)
    be0 = bnd_edge_count(st0["c8s"])
    if cinfo.get("virtual"):
        say(f"[class] CAD topology from block-face labels: {cinfo['virtual']['edges']} edge chains, {cinfo['virtual']['vertices']} vertices")
    say(f"[input] boundary block edges by face count {dict(collections.Counter(be0.values()))}")
    allv = set(range(len(cad.V))); alle = set(cad.eid.tolist())
    say(f"[class] CAD vertices represented {len(ref_cov[0])}/{len(allv)} (missing {sorted(allv - ref_cov[0])}), "
        f"CAD edges represented {len(ref_cov[1])}/{len(alle)} (missing {sorted(alle - set(ref_cov[1]))})")
    allQ = []
    for F, S in surf.items():
        allQ.append((S.A, S.B, S.C))
    surf_all = TriSurface(np.vstack([x[0] for x in allQ]), np.vstack([x[1] for x in allQ]), np.vstack([x[2] for x in allQ]))
    rng = np.random.default_rng(0)
    rev = cad.fp[rng.choice(len(cad.fp), min(a.nrev, len(cad.fp)), replace=False)]
    # baseline re-fit (same pipeline, no collapse)
    st0["_pos"] = node_pos(st0)
    t = time.time()
    Xb, ib = realise(st0, dict(c8s=st0["c8s"], ncls=st0["ncls"], ecls=st0["ecls"], fcls=st0["fcls"]), None, cad, surf, a.h)
    base = evaluate(Xb, st0, ib, cad, surf_all, rev)
    rev_base = base.pop("_drev")
    if a.rev_relative:
        G_["rev_base"] = rev_base
    say(f"[baseline h={a.h}] {base}  ({time.time() - t:.1f}s)")
    st_in_q = evaluate([X for X in st0["X"]], st0, dict(fallback_blocks=0, face_grid={k: v for k, v in face_table(st0).items()}), cad, surf_all, rev)
    st_in_q.pop("_drev", None)
    say(f"[input as given] {st_in_q}")
    if a.sj_floor == "input":
        a.sj_floor = st_in_q["min_sj"] - 1e-9; floor_note = f"input min sJ {st_in_q['min_sj']:.4f}"
    else:
        a.sj_floor = float(a.sj_floor); floor_note = "absolute"
    h1lo = h1hi = None
    if cad.wall:
        h1lo = st_in_q["wall_h1_min"] * (1 - a.h1_tol); h1hi = st_in_q["wall_h1_max"] * (1 + a.h1_tol)
    sheets0 = sheets_of(st0["c8s"])
    inv = []
    for i, s in enumerate(sheets0):
        d = describe_sheet(st0, s, meta); d["id"] = i; inv.append(d)
        say(f"[sheet {i}] {d['n']} blocks, kinds {d['kinds']}, azimuthal {d['azimuthal_frac']:.2f}, 2D blocks {d['blocks2d']}")
    hrev_c = base["h_rev"] + a.htol; hfwd_c = base["h_fwd"] + a.htol
    say(f"[guard] inverted 0, min sJ >= {a.sj_floor:.4f} ({floor_note}), "
        + (f"per-sample CAD->mesh increase <= {a.htol}" if a.rev_relative else f"h_rev <= {hrev_c:.4f}")
        + f", h_fwd <= {hfwd_c:.4f}" + (f", wall first cell in [{h1lo:.4e}, {h1hi:.4e}]" if h1lo else ""))
    st1 = dict(st0); st1["X"] = Xb
    beam = [(base, st1)]; best = (base, st1); seen = {frozenset(frozenset(c) for c in st1["c8s"])}
    records = []
    for depth in range(a.depth):
        cand = []
        for ev0, st in beam:
            st["_pos"] = node_pos(st)
            sheets = sheets_of(st["c8s"])
            G_.update(st=st, sheets=sheets, cad=cad, surf=surf, surf_all=surf_all, rev=rev, h=a.h, ref_cov=ref_cov)
            with get_context("fork").Pool(a.workers) as pool:
                res = pool.map(_eval_child, [(i, v) for i in range(len(sheets)) for v in ("mid", "lo", "hi")], chunksize=1)
            bestvar = {}
            for i, var, kind, ev, child, dt in res:
                if kind == "variant":
                    continue
                d = describe_sheet(st, sheets[i], meta)
                rec = dict(depth=depth + 1, variant=var, parent_path=st["path"], parent_blocks=len(st["c8s"]), sheet_n=d["n"],
                           sheet_kinds=d["kinds"], azimuthal=round(d["azimuthal_frac"], 2), blocks2d=d["blocks2d"],
                           sheet_origin=d["origin"], time=round(dt, 1))
                if kind != "ok":
                    rec.update(result="reject", reason=kind, detail=ev)
                else:
                    why = []
                    if ev["inverted"] > 0:
                        why.append(f"inverted {ev['inverted']}")
                    if ev["min_sj"] < a.sj_floor:
                        why.append(f"min sJ {ev['min_sj']:.3f}")
                    if a.rev_relative and ev["dh_rev"] > a.htol:
                        why.append(f"CAD->mesh +{ev['dh_rev']:.4f}")
                    if not a.rev_relative and ev["h_rev"] > hrev_c:
                        why.append(f"h_rev {ev['h_rev']:.4f}")
                    if h1lo and not (h1lo <= ev["wall_h1_min"] and ev["wall_h1_max"] <= h1hi):
                        why.append(f"wall first cell {ev['wall_h1_min']:.3e}..{ev['wall_h1_max']:.3e}")
                    if ev["h_fwd"] > hfwd_c:
                        why.append(f"h_fwd {ev['h_fwd']:.4f}")
                    rec.update(result="admissible" if not why else "reject", reason="geometry" if why else "",
                               detail="; ".join(why), eval={k: v for k, v in ev.items()})
                    if not why:
                        sig = frozenset(frozenset(c) for c in child["c8s"])
                        if sig not in seen and (sig not in bestvar or ev["min_sj"] > bestvar[sig][0]["min_sj"]):
                            child["path"] = st["path"] + [dict(sheet_origin=d["origin"], variant=var)]
                            bestvar[sig] = (ev, child)
                records.append(rec)
                say(f"[d{depth + 1} parent {len(st['c8s'])}bl] sheet {i}/{var} ({d['n']} bl, {d['kinds']}, az {d['azimuthal_frac']:.2f}): "
                    f"{rec['result']} {rec.get('reason', '')} {rec.get('detail', '')} "
                    + (f"-> {ev['blocks']} bl, min sJ {ev['min_sj']:.3f}, inv {ev['inverted']}, h_rev {ev['h_rev']:.4f}, h_fwd {ev['h_fwd']:.4f}, fallback {ev['fallback_blocks']}" if kind == "ok" else "")
                    + f" ({dt:.1f}s)")
            for sig, c in bestvar.items():
                if sig not in seen:
                    seen.add(sig); cand.append(c)
        if not cand:
            say(f"[beam] depth {depth + 1}: no admissible collapse, stop"); break
        cand.sort(key=lambda c: (c[0]["blocks"], -c[0]["min_sj"]))
        beam = cand[:a.width]
        for ev, ch in beam:
            if (ev["blocks"], -ev["min_sj"]) < (best[0]["blocks"], -best[0]["min_sj"]):
                best = (ev, ch)
        say(f"[beam] depth {depth + 1}: {len(cand)} admissible, beam {[(e['blocks'], round(e['min_sj'], 3)) for e, _ in beam]}, best {best[0]['blocks']}")
    json.dump(dict(baseline=base, input=st_in_q, sheets=inv, records=records, best=best[0], best_path=best[1]["path"],
                   guard=dict(sj_floor=a.sj_floor, floor=floor_note, h_rev=hrev_c, h_fwd=hfwd_c, h=a.h, h1=[h1lo, h1hi],
                              rev_relative=a.rev_relative), args=vars(a)),
              open(f"{a.out}/records.json", "w"), indent=1, default=str)
    write_state(f"{a.out}/baseline_h{a.h}", st1, base, cad)
    if best[1] is not st1:
        write_state(f"{a.out}/best_h{a.h}", best[1], best[0], cad)
    say(f"[done] best {best[0]['blocks']} blocks via {len(best[1]['path'])} collapses; {time.time() - T0:.0f}s")


def face_table(st):
    out = {}
    for X, c8 in zip(st["X"], st["c8s"]):
        C = corner_cube(c8)
        for a, s in FACES:
            k = frozenset(np.take(C, s, axis=a).ravel().tolist())
            out.setdefault(k, (face_slice(X, a, s), np.take(C, s, axis=a)))
    return out


def write_state(prefix, st, ev, cad):
    pts, H, bl, sj = [], [], [], []
    off = 0
    for bi, X in enumerate(st["X"]):
        n1, n2, n3 = (s - 1 for s in X.shape[:3])
        idx = np.arange((n1 + 1) * (n2 + 1) * (n3 + 1)).reshape(n1 + 1, n2 + 1, n3 + 1) + off
        c = [idx[:-1, :-1, :-1], idx[1:, :-1, :-1], idx[1:, 1:, :-1], idx[:-1, 1:, :-1],
             idx[:-1, :-1, 1:], idx[1:, :-1, 1:], idx[1:, 1:, 1:], idx[:-1, 1:, 1:]]
        H.append(np.stack([x.ravel() for x in c], 1)); pts.append(X.reshape(-1, 3)); off += idx.size
        bl += [bi] * (n1 * n2 * n3); sj.append(hex_sj(X).ravel())
    write_hex_vtk(prefix + "_tfi.vtk", np.vstack(pts), np.vstack(H),
                  [("block", np.array(bl)), ("origin", np.array([st["origin"][b] for b in bl])), ("sJ", np.concatenate(sj))])
    lines, lid = [], []
    for bi, (X, c8) in enumerate(zip(st["X"], st["c8s"])):
        C = corner_cube(c8)
        for ax in range(3):
            for p in (0, 1):
                for q in (0, 1):
                    L, u, v = block_edge_line(X, C, ax, p, q); lines.append(L); lid.append(bi)
    write_lines_vtk(prefix + "_block_edges.vtk", lines, [("block", lid)])
    # coarse block mesh (one hex per block)
    pos = node_pos(st); ids = sorted(pos); im = {v: i for i, v in enumerate(ids)}
    Hc = np.array([[im[v] for v in c8] for c8 in st["c8s"]])
    write_hex_vtk(prefix + "_coarse.vtk", np.array([pos[v] for v in ids]), Hc, [("block", np.arange(len(Hc))), ("origin", np.array(st["origin"]))],
                  pdata=[("node_dim", np.array([st["ncls"][v][0] for v in ids])), ("node_cad", np.array([st["ncls"][v][1] for v in ids]))])
    json.dump(dict(eval=ev, c8s=[list(map(int, c)) for c in st["c8s"]], origin=st["origin"], path=st["path"]),
              open(prefix + ".json", "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
