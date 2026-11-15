"""Cell-size transition metrics of a hex mesh written by sweep3d*.py (legacy binary VTK with cell data block, kind, sJ).
For every interior face (2 hexes A,B):  volume ratio max/min, and per-direction size ratio = ratio of the two hex
thicknesses normal to that face (mean length of the 4 hex edges joining the face to its opposite face).
Groups: all faces; layer/outside interface (kind 3 vs 0); layer-internal wall-normal faces (both hexes kind 3 and the
face normal is the thinnest direction of both hexes = wall-normal); strip and butterfly (kinds 4, 1, 2) involvement.
usage: python3 eval_grading.py mesh.vtk out.json"""
import sys, json, numpy as np

def read_vtk(fn):
    b = open(fn, "rb").read(); pos = 0
    def line():
        nonlocal pos
        e = b.index(b"\n", pos); s = b[pos:e].decode(); pos = e + 1; return s
    data = {}
    while pos < len(b):
        s = line().strip()
        if s.startswith("POINTS"):
            n = int(s.split()[1]); P = np.frombuffer(b, ">f8", 3 * n, pos).reshape(n, 3); pos += 24 * n
        elif s.startswith("CELLS"):
            n, m = map(int, s.split()[1:]); C = np.frombuffer(b, ">i4", m, pos).reshape(n, 9)[:, 1:]; pos += 4 * m
        elif s.startswith("CELL_TYPES"):
            n = int(s.split()[1]); pos += 4 * n
        elif s.startswith("SCALARS"):
            name, typ = s.split()[1:3]; line()
            dt = ">i4" if typ == "int" else ">f8"; sz = 4 if typ == "int" else 8
            data[name] = np.frombuffer(b, dt, len(C), pos); pos += sz * len(C)
    return P, C, data

P, C, D = read_vtk(sys.argv[1]); kind = D["kind"]; sj = D["sJ"]
X = P[C]                                                    # (n, 8, 3)
def L(a, b): return np.linalg.norm(X[:, a] - X[:, b], axis=1)
thick = np.stack([(L(0, 4) + L(1, 5) + L(2, 6) + L(3, 7)) / 4,       # dir 0: faces 0 (0123) / 1 (4567)
                  (L(0, 3) + L(1, 2) + L(5, 6) + L(4, 7)) / 4,       # dir 1: faces 2 (0154) / 4 (2376)
                  (L(1, 0) + L(2, 3) + L(6, 7) + L(5, 4)) / 4], 1)   # dir 2: faces 3 (1265) / 5 (3047)
# volume: 6 tets of the hex (VTK ordering)
def tv(a, b, c, d): return np.einsum("ij,ij->i", np.cross(X[:, b] - X[:, a], X[:, c] - X[:, a]), X[:, d] - X[:, a]) / 6
V = np.abs(tv(0, 1, 3, 4) + tv(1, 2, 3, 6) + tv(1, 3, 4, 6) + tv(3, 4, 6, 7) + tv(1, 4, 5, 6) + tv(1, 3, 6, 4) * 0)
V = np.abs(tv(0, 1, 2, 5)) + np.abs(tv(0, 2, 3, 7)) + np.abs(tv(0, 5, 7, 4)) + np.abs(tv(2, 7, 5, 6)) + np.abs(tv(0, 2, 7, 5))
FACES = [(0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
FDIR = [0, 0, 1, 2, 1, 2]
keys = np.vstack([np.sort(C[:, f], axis=1) for f in FACES])
hexid = np.tile(np.arange(len(C)), 6); fdir = np.repeat(FDIR, len(C))
order = np.lexsort(keys.T[::-1]); ks = keys[order]
same = np.all(ks[1:] == ks[:-1], axis=1); i0 = order[:-1][same]; i1 = order[1:][same]
A, B = hexid[i0], hexid[i1]; dA, dB = fdir[i0], fdir[i1]
vr = np.maximum(V[A], V[B]) / np.maximum(np.minimum(V[A], V[B]), 1e-300)
tA = thick[A, dA]; tB = thick[B, dB]
dr = np.maximum(tA, tB) / np.maximum(np.minimum(tA, tB), 1e-300)
def st(x): return dict(n=int(len(x)), p50=float(np.percentile(x, 50)), p95=float(np.percentile(x, 95)),
                       p99=float(np.percentile(x, 99)), max=float(x.max())) if len(x) else dict(n=0)
kA, kB = kind[A], kind[B]
iface = ((kA == 3) & (kB == 0)) | ((kA == 0) & (kB == 3))
wn = (kA == 3) & (kB == 3) & (tA <= thick[A].min(1) * 1.0001) & (tB <= thick[B].min(1) * 1.0001)
strip = (kA == 4) | (kB == 4); bfly = np.isin(kA, (1, 2)) | np.isin(kB, (1, 2))
# first cell height on walls: layer hexes with a boundary face whose direction is the thinnest: thickness there
out = dict(hexes=int(len(C)), inverted=int((sj <= 0).sum()), min_sJ=float(sj.min()), p1_sJ=float(np.percentile(sj, 1)),
           median_sJ=float(np.median(sj)),
           volume_ratio=dict(all=st(vr), layer_outside=st(vr[iface]), strip=st(vr[strip]), butterfly=st(vr[bfly])),
           dir_ratio=dict(all=st(dr), layer_outside=st(dr[iface]), layer_wall_normal=st(dr[wn]), strip=st(dr[strip]),
                          butterfly=st(dr[bfly])),
           layer_thinnest_size=dict(min=float(thick[kind == 3].min(1).min()), median=float(np.median(thick[kind == 3].min(1)))))
json.dump(out, open(sys.argv[2], "w"), indent=1)
for k in ("volume_ratio", "dir_ratio"):
    for g, v in out[k].items():
        print(f"{k:12s} {g:18s} " + " ".join(f"{a}={b:.3g}" if isinstance(b, float) else f"{a}={b}" for a, b in v.items()))
print({k: out[k] for k in ("hexes", "inverted", "min_sJ", "p1_sJ", "median_sJ", "layer_thinnest_size")})
