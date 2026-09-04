"""Stage 4c: block topology for TFI, via labelled separating sheets.

This is the 3D analogue of the 2D separatrix tracing in
`experimentell/gmsh_pipeline/extract_separatrices.py`: there, one walks
edges out of a singular VERTEX and continues straight (`idx + n//2` in the
cyclic neighbour order) until hitting another singularity or the boundary;
the separatrix network then bounds the quad blocks.

In 3D the same idea applies one dimension up: separating SURFACES emanate
from singular EDGES, are propagated straight across regular edges, and stop
at another singular edge or at the boundary (`base_complex.sheet_faces`).
Each maximal such surface is one "sheet" -- the analogue of one separatrix.

A block face is then simply a maximal connected patch of the block's
boundary carrying ONE label, where the label is either
  * the id of the interior sheet the patch lies in, or
  * the physical surface (hub / shroud / blade / inlet / outlet / periodic)
    the patch lies on.
A block is a topological cuboid iff it has exactly 6 such faces.

Using the physical surfaces rather than a single "boundary" label matters:
lumping all domain-boundary faces together merges e.g. a cylinder's lateral
wall with its end caps across the rim, which made proper cuboids look like
4-faced blocks. It is also exactly what TFI needs -- every block face has to
lie on one surface.

NOTE: an earlier attempt assigned (i,j,k) grid coordinates per block by
walking cell to cell. That is unnecessary for TFI and was abandoned; block
faces/edges/corners come straight out of the sheet network, as in 2D.
"""

import sys
from collections import Counter, defaultdict, deque
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO / "experimentell" / "hex3d_algohex"))

import base_complex as bc                                             # noqa: E402
import ovm_io                                                         # noqa: E402
import tet_prep as tp                                                 # noqa: E402

HF = bc.HEX_FACES


def _face_edges_map(hexes):
    fe = {}
    for c in hexes:
        for fc in HF:
            lp = [int(c[i]) for i in fc]
            fk = frozenset(lp)
            if fk not in fe:
                fe[fk] = [(min(lp[k], lp[(k + 1) % 4]),
                           max(lp[k], lp[(k + 1) % 4])) for k in range(4)]
    return fe


def label_sheets(hexes, f2h, sing):
    """Label every sheet separately: the equivalence classes of the
    straight-across propagation. One sheet == one 3D separatrix."""
    fe = _face_edges_map(hexes)

    def other_face(fk, g, hi):
        c = hexes[hi]
        for fc in HF:
            lp = [int(c[i]) for i in fc]
            k2 = frozenset(lp)
            if k2 == fk:
                continue
            for k in range(4):
                u, v = lp[k], lp[(k + 1) % 4]
                if (min(u, v), max(u, v)) == g:
                    return k2
        return None

    def straight(fk, g):
        for h0 in f2h.get(fk, ()):
            f1 = other_face(fk, g, h0)
            if f1 is None:
                continue
            h1 = [h for h in f2h.get(f1, ()) if h != h0]
            if not h1:
                continue
            f2 = other_face(f1, g, h1[0])
            if f2 is not None and f2 != fk:
                return f2
        return None

    seeds = {fk for fk, es in fe.items() if any(g in sing for g in es)}
    lab, sid = {}, 0
    for s in seeds:
        if s in lab:
            continue
        dq = deque([s])
        lab[s] = sid
        while dq:
            fk = dq.popleft()
            for g in fe[fk]:
                if g in sing:
                    continue
                nxt = straight(fk, g)
                if nxt is not None and nxt not in lab:
                    lab[nxt] = sid
                    dq.append(nxt)
        sid += 1
    print(f"[block_faces] {sid} sheets over {len(lab)} faces")
    return lab, sid


def physical_of(P, loop):
    """Physical surface of a boundary quad (same criteria as
    tet_prep.classify_boundary, applied to the hex mesh)."""
    Q = P[loop]
    cen = Q.mean(0)
    n = np.cross(Q[2] - Q[0], Q[3] - Q[1])
    L = np.linalg.norm(n)
    n = n / L if L > 0 else n
    r = np.hypot(cen[0], cen[1])
    rad = np.array([cen[0] / max(r, 1e-12), cen[1] / max(r, 1e-12), 0.0])
    if abs(r - tp.R_HUB) < 0.03 and abs(np.dot(n, rad)) > 0.7:
        return tp.SURF_HUB
    if abs(r - tp.R_SHROUD) < 0.03 and abs(np.dot(n, rad)) > 0.7:
        return tp.SURF_SHROUD
    if abs(cen[2] - tp.Z_INLET) < 0.02 and abs(n[2]) > 0.7:
        return tp.SURF_INLET
    if abs(cen[2] - tp.Z_OUTLET) < 0.02 and abs(n[2]) > 0.7:
        return tp.SURF_OUTLET
    th = np.degrees(np.arctan2(cen[1], cen[0]))
    if th < -9.0:
        return tp.SURF_PER_A
    if th > 22.0:
        return tp.SURF_PER_B
    return tp.SURF_BLADE


def block_faces(hexes, P, f2h, blocks, blk_of, lab):
    """Per block: list of faces, each a (label, [quad loops]) patch."""
    out = []
    for bi, b in enumerate(blocks):
        bf = []
        for h in b:
            c = hexes[h]
            for fc in HF:
                lp = [int(c[i]) for i in fc]
                k = frozenset(lp)
                nb = [x for x in f2h[k] if x != h]
                if not nb:
                    bf.append((("P", physical_of(P, lp)), lp))
                elif blk_of[nb[0]] != bi:
                    bf.append((("S", lab.get(k, -1)), lp))
        e2f = defaultdict(list)
        for fi, (_lb, lp) in enumerate(bf):
            for a in range(4):
                u, v = lp[a], lp[(a + 1) % 4]
                e2f[(min(u, v), max(u, v))].append(fi)
        adj = defaultdict(list)
        for fl in e2f.values():
            for x in range(len(fl)):
                for y in range(x + 1, len(fl)):
                    if bf[fl[x]][0] == bf[fl[y]][0]:
                        adj[fl[x]].append(fl[y])
                        adj[fl[y]].append(fl[x])
        seen, patches = set(), []
        for s in range(len(bf)):
            if s in seen:
                continue
            comp, st = [], [s]
            seen.add(s)
            while st:
                u = st.pop()
                comp.append(u)
                for w in adj[u]:
                    if w not in seen:
                        seen.add(w)
                        st.append(w)
            patches.append((bf[s][0], [bf[i][1] for i in comp]))
        out.append(patches)
    return out


def analyse(ovm_path):
    P, e, f, poly = ovm_io.read_ovm(ovm_path)
    hexes, _ = ovm_io.ovm_to_cells(P, e, f, poly)
    f2h, e2h = bc.build_topology(hexes)
    sing = bc.singular_edges(hexes, P, f2h, e2h)
    lab, nsheet = label_sheets(hexes, f2h, sing)
    blocks = bc.blocks_from_cut(hexes, f2h, set(lab))
    blk_of = np.zeros(len(hexes), int)
    for i, b in enumerate(blocks):
        for h in b:
            blk_of[h] = i
    faces = block_faces(hexes, P, f2h, blocks, blk_of, lab)
    cnt = Counter(len(fs) for fs in faces)
    cub = cnt.get(6, 0)
    cells_ok = sum(len(b) for b, fs in zip(blocks, faces) if len(fs) == 6)
    print(f"[block_faces] {len(blocks)} blocks, faces per block: "
          f"{dict(sorted(cnt.items()))}")
    print(f"[block_faces] CUBOIDS (exactly 6 faces): {cub}/{len(blocks)} "
          f"({100 * cub / len(blocks):.0f}%), covering {cells_ok}/{len(hexes)} "
          f"cells ({100 * cells_ok / len(hexes):.0f}%)")
    return P, hexes, blocks, faces, blk_of


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("ovm", nargs="?",
                    default=str(REPO / "output" / "hex3d_algohex"
                               / "T1_9_hex_v5.ovm"))
    ap.add_argument("--out", default=str(REPO / "output" / "hex3d_algohex"
                                        / "deliverable" / "T1_9_blocks_v5.vtk"))
    a = ap.parse_args()
    P, hexes, blocks, faces, blk_of = analyse(a.ovm)
    import export_vtk as ev
    nf = np.zeros(len(hexes), int)
    for b, fs in zip(blocks, faces):
        for h in b:
            nf[h] = len(fs)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    ev.write_vtk(a.out, P, hexes, [12] * len(hexes), blk_of, "block_id",
                 f"{len(blocks)} blocks from sheet network")
    ev.write_vtk(a.out.replace(".vtk", "_nfaces.vtk"), P, hexes,
                 [12] * len(hexes), nf, "n_block_faces",
                 "faces per block (6 = cuboid, ready for TFI)")
