"""Evaluate an AlgoHex hex mesh of the slab test: size, inverted cells, boundary fit, coverage of the
(planar) plate mid surfaces by hex faces, raw base-complex block count."""
import sys, json, numpy as np, meshio
sys.path.insert(0, "/opt/stack/dp3d_keep_bl/experimentell/hex3d_algohex")
import ovm_io, base_complex as bc
from scipy.spatial import cKDTree
ovm, vtk = sys.argv[1], sys.argv[2]
P, e, f, poly = ovm_io.read_ovm(ovm); H, skipped = ovm_io.ovm_to_cells(P, e, f, poly)
sj = np.array([ovm_io.scaled_jacobian(P, c) for c in H])
m = meshio.read(vtk)
cd = {c.type: (c.data, np.asarray(d).ravel()) for c, d in zip(m.cells, m.cell_data["color"])}
T, TL = cd["triangle"]; Q = m.points
# boundary quads of the hex mesh
FACES = [(0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
from collections import Counter
fc = Counter(); rep = {}
for h in H:
    for fa in FACES:
        k = tuple(sorted(h[list(fa)])); fc[k] += 1; rep[k] = h[list(fa)]
bq = np.array([rep[k] for k, c in fc.items() if c == 1]); iq = np.array([rep[k] for k, c in fc.items() if c == 2])
tc = Q[T].mean(1); tree = cKDTree(tc)
d_b, _ = tree.query(P[bq].mean(1)) if len(bq) else (np.zeros(0), None)
h_avg = np.mean(np.linalg.norm(P[H[:, 0]] - P[H[:, 1]], axis=1))
res = dict(hexes=len(H), skipped=int(skipped), inverted=int((sj <= 0).sum()), min_sj=float(sj.min()),
           mean_edge=float(h_avg), boundary_quads=len(bq),
           boundary_quads_off_surface=int((d_b > 0.5 * h_avg).sum()), max_boundary_dist=float(d_b.max()))
# plate mid surfaces (labels 7, 8 and their crack copies 107, 108) -> fitted planes
for lab in (7, 8):
    sel = np.isin(TL, [lab, lab + 100]); 
    if not sel.any(): continue
    X = Q[np.unique(T[sel])]; c = X.mean(0); _, _, vt = np.linalg.svd(X - c); n = vt[2]
    flat = np.abs((X - c) @ n).max()
    area = 0.5 * np.linalg.norm(np.cross(Q[T[sel][:, 1]] - Q[T[sel][:, 0]], Q[T[sel][:, 2]] - Q[T[sel][:, 0]]), axis=1).sum()
    if flat > 1e-6:    # thick plate or envelope: not a mid surface
        continue
    if lab + 100 in TL: area /= 2
    allq = np.vstack([bq, iq]); V = P[allq]
    on = (np.abs((V - c) @ n) < 0.05 * h_avg).all(1)
    # inside the baffle footprint: nearest baffle triangle centroid close
    ctr = V[on].mean(1); dd, _ = cKDTree(tc[sel]).query(ctr); on_idx = np.where(on)[0][dd < 0.5 * h_avg]
    qa = 0.5 * (np.linalg.norm(np.cross(V[on_idx, 1] - V[on_idx, 0], V[on_idx, 2] - V[on_idx, 0]), axis=1)
                + np.linalg.norm(np.cross(V[on_idx, 2] - V[on_idx, 0], V[on_idx, 3] - V[on_idx, 0]), axis=1))
    # for a crack both sheets are boundary quads -> count unique positions once
    cov = qa.sum() / area
    if lab + 100 in TL: cov /= 2
    res[f"plate{lab}_coverage"] = float(cov)
f2h, e2h = bc.build_topology(H); sing = bc.singular_edges(H, P, f2h, e2h, 40.0)
cut = bc.sheet_faces(H, f2h, sing); blocks = bc.blocks_from_cut(H, f2h, cut)
res["singular_edges"] = int(len(sing)); res["raw_blocks"] = len(blocks)
print(json.dumps(res, indent=1))
json.dump(res, open(ovm.replace(".ovm", "_eval.json"), "w"), indent=1)
