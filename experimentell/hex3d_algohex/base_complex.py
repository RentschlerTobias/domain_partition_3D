"""Stage 4: base complex of a hex mesh -> coarse blocks for TFI.

The AlgoHex output is a FINE hex mesh; what the TFI pipeline needs is the
coarse block structure. That structure is the mesh's *base complex*: the
partition induced by the sheets that emanate from irregular (singular)
edges.

Algorithm (Gao et al. 2015 / 2017, "Hexahedral Mesh Reparameterization from
Aligned Base-Complex" / "Robust structure simplification for hex re-meshing",
as summarised in the hex-meshing survey): starting from all quad facets
incident to any singular edge, iteratively expand through *opposite* facets
across *regular* edges until termination. The resulting facet sets are the
base-complex separating surfaces ("sheets"); the hex components they cut the
mesh into are the coarse blocks.

Definitions used here:
  - interior edge: singular iff its incident-hex valence != 4
  - boundary edge: singular iff its incident-hex valence != 2, or the
    dihedral angle across the surface is sharp (a geometric feature must
    also break a block, otherwise a block would span a CAD corner)
"""

import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO / "experimentell" / "hex3d_algohex"))

import ovm_io                                                         # noqa: E402

# the 6 faces and 12 edges of a hex in VTK_HEXAHEDRON ordering
HEX_FACES = ((0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4),
             (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7))
HEX_EDGES = ((0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4),
             (0, 4), (1, 5), (2, 6), (3, 7))


def build_topology(hexes):
    """face -> incident hexes, edge -> incident hexes."""
    f2h, e2h = defaultdict(list), defaultdict(list)
    for hi, c in enumerate(hexes):
        for f in HEX_FACES:
            f2h[frozenset(int(c[i]) for i in f)].append(hi)
        for a, b in HEX_EDGES:
            u, v = int(c[a]), int(c[b])
            e2h[(min(u, v), max(u, v))].append(hi)
    return f2h, e2h


def singular_edges(hexes, points, f2h, e2h, sharp_deg=40.0,
                   interior_only=True):
    """Edges that must break the block structure.

    ``interior_only`` (default) uses only interior edges whose hex valence
    != 4 -- the actual topological singularities. Validated on AlgoHex's
    cylinder demo: it yields exactly the 5 blocks of the textbook O-grid
    decomposition (4 outer + core).

    Including boundary-valence and discrete sharp-angle criteria explodes
    the count (cylinder 5 -> 160 blocks, T1_9 21 -> 2485 with 991 single-cell
    blocks): on the boundary those tests fire on ordinary mesh irregularity
    rather than on real structure. Kept switchable for inspection.
    """
    # a boundary edge is one that belongs to at least one boundary face
    bnd_faces = {f for f, hs in f2h.items() if len(hs) == 1}
    bnd_edges = set()
    face_of_edge = defaultdict(list)
    for hi, c in enumerate(hexes):
        for f in HEX_FACES:
            key = frozenset(int(c[i]) for i in f)
            if key not in bnd_faces:
                continue
            loop = [int(c[i]) for i in f]
            for k in range(4):
                u, v = loop[k], loop[(k + 1) % 4]
                e = (min(u, v), max(u, v))
                bnd_edges.add(e)
                face_of_edge[e].append(loop)

    cos_thr = np.cos(np.radians(sharp_deg))

    def quad_normal(loop):
        P = points[loop]
        n = np.cross(P[2] - P[0], P[3] - P[1])
        ln = np.linalg.norm(n)
        return n / ln if ln > 0 else n

    sing = set()
    n_int = n_bnd = n_sharp = 0
    for e, hs in e2h.items():
        val = len(set(hs))
        if e in bnd_edges:
            if interior_only:
                continue
            if val != 2:
                sing.add(e)
                n_bnd += 1
            else:
                fs = face_of_edge.get(e, [])
                if len(fs) == 2:
                    d = float(np.dot(quad_normal(fs[0]), quad_normal(fs[1])))
                    if d < cos_thr:
                        sing.add(e)
                        n_sharp += 1
        else:
            if val != 4:
                sing.add(e)
                n_int += 1
    print(f"[base_complex] singular edges: {len(sing)} "
          f"({n_int} interior valence!=4, {n_bnd} boundary valence!=2, "
          f"{n_sharp} sharp feature)")
    return sing


def sheet_faces(hexes, f2h, sing):
    """Faces of the base-complex separating surfaces.

    Gao et al.: start from all facets incident to a singular edge and expand
    through *opposite facets across regular edges*. The expansion is across
    the EDGES of a face, which is what makes the result a 2D surface.

    Stepping to the opposite face *within* a hex instead walks a 1D chain of
    faces and never closes into a separating surface -- that mistake left the
    mesh in one giant component (9384 of 9792 hexes on the cylinder).

    "Opposite facet across edge g" = rotate twice around g: from face f enter
    the hex that owns both f and g, take that hex's other face containing g,
    then repeat once from there. For a regular interior edge (4 hexes) this
    lands on the face diametrically across g, i.e. straight ahead.
    """
    # face -> its 4 edges; edge -> incident faces; (face, hex) adjacency
    f_edges = {}
    e2f = defaultdict(set)
    for hi, c in enumerate(hexes):
        for f in HEX_FACES:
            loop = [int(c[i]) for i in f]
            fk = frozenset(loop)
            if fk not in f_edges:
                es = []
                for k in range(4):
                    u, v = loop[k], loop[(k + 1) % 4]
                    es.append((min(u, v), max(u, v)))
                f_edges[fk] = es
                for g in es:
                    e2f[g].add(fk)
            else:
                for g in f_edges[fk]:
                    e2f[g].add(fk)

    def other_face_in_hex(fk, g, hi):
        """the other face of hex hi that also contains edge g"""
        c = hexes[hi]
        for f in HEX_FACES:
            loop = [int(c[i]) for i in f]
            k2 = frozenset(loop)
            if k2 == fk:
                continue
            for k in range(4):
                u, v = loop[k], loop[(k + 1) % 4]
                if (min(u, v), max(u, v)) == g:
                    return k2
        return None

    def straight_across(fk, g):
        """rotate twice around g starting from fk"""
        hs = [h for h in f2h.get(fk, ()) ]
        for h0 in hs:
            f1 = other_face_in_hex(fk, g, h0)
            if f1 is None:
                continue
            h1 = [h for h in f2h.get(f1, ()) if h != h0]
            if not h1:
                continue
            f2 = other_face_in_hex(f1, g, h1[0])
            if f2 is not None and f2 != fk:
                return f2
        return None

    seeds = set()
    for fk, es in f_edges.items():
        if any(g in sing for g in es):
            seeds.add(fk)

    cut = set()
    stack = list(seeds)
    while stack:
        fk = stack.pop()
        if fk in cut:
            continue
        cut.add(fk)
        for g in f_edges[fk]:
            if g in sing:
                continue           # a singular edge terminates the surface
            nxt = straight_across(fk, g)
            if nxt is not None and nxt not in cut:
                stack.append(nxt)
    print(f"[base_complex] sheet faces: {len(cut)} of {len(f_edges)} "
          f"(from {len(seeds)} seed faces)")
    return cut


def blocks_from_cut(hexes, f2h, cut):
    """Connected hex components after removing the sheet faces = blocks."""
    adj = defaultdict(list)
    for fk, hs in f2h.items():
        if len(hs) == 2 and fk not in cut:
            a, b = hs
            adj[a].append(b)
            adj[b].append(a)
    seen, blocks = set(), []
    for s in range(len(hexes)):
        if s in seen:
            continue
        stack, comp = [s], []
        seen.add(s)
        while stack:
            u = stack.pop()
            comp.append(u)
            for w in adj[u]:
                if w not in seen:
                    seen.add(w)
                    stack.append(w)
        blocks.append(comp)
    return blocks


def analyse(ovm_path, sharp_deg=40.0):
    P, e, f, poly = ovm_io.read_ovm(ovm_path)
    hexes, _ = ovm_io.ovm_to_cells(P, e, f, poly)
    print(f"[base_complex] {len(hexes)} hexes, {len(P)} vertices")
    f2h, e2h = build_topology(hexes)
    sing = singular_edges(hexes, P, f2h, e2h, sharp_deg)
    cut = sheet_faces(hexes, f2h, sing)
    blocks = blocks_from_cut(hexes, f2h, cut)
    sizes = sorted((len(b) for b in blocks), reverse=True)
    print(f"[base_complex] BLOCKS: {len(blocks)}")
    print(f"[base_complex]   largest: {sizes[:10]}")
    print(f"[base_complex]   1-cell blocks: {sum(1 for s in sizes if s == 1)}")
    print(f"[base_complex]   mean cells/block: {np.mean(sizes):.1f}")
    return P, hexes, blocks, sing, cut


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("ovm", nargs="?",
                    default=str(REPO / "output" / "hex3d_algohex"
                               / "T1_9_hex_v4.ovm"))
    ap.add_argument("--sharp-deg", type=float, default=40.0)
    ap.add_argument("--out", default=str(REPO / "output" / "hex3d_algohex"
                                        / "vtk" / "07_base_complex.vtk"))
    a = ap.parse_args()
    P, hexes, blocks, sing, cut = analyse(a.ovm, a.sharp_deg)
    bid = np.zeros(len(hexes), int)
    for i, b in enumerate(blocks):
        for h in b:
            bid[h] = i
    import export_vtk as ev
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    ev.write_vtk(a.out, P, hexes, [12] * len(hexes), bid, "block_id",
                 f"base complex: {len(blocks)} blocks")
