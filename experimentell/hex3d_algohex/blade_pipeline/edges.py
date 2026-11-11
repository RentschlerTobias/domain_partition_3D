"""B4: every block edge once, as a curve on the exact geometry.

    edges.py <run_dir> [--blocks blade/blocks_b3.vtk] [--sub 8] [--out blade/edges_b4.npz]

A block edge is the lattice line between two block corners (tfi.lattices),
stored once however many blocks share it. Its nodes carry the reference
surfaces they were put on (b3_ring.json vertex_surfaces: B2 labels for the
core, B3 rows for the ring; never a distance test, which the facetted hub vs
the analytic blade root would fail); the cutout is interior after B3. The edge's
constraint is what ALL its nodes share:
  2 surfaces -> a curve: blade + hub / shroud = the iso line v = 0 / 1 of S,
                otherwise the dtOO 1D entity on both surfaces,
  1 surface  -> that surface,
  none       -> interior edge, free.
The curve is a centripetal Catmull-Rom spline through the lattice nodes,
sampled --sub points per lattice segment and then projected onto its
constraint (nodes stay where they are). Report: projection residual, how far
the curve leaves the straight chords, and a fold check (no sample may run
backwards along its segment).
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))

import tfi                                                            # noqa: E402
from geometry import Geometry                                         # noqa: E402
from ring import BladeSurface                                         # noqa: E402
from snap import curve_surfaces                                       # noqa: E402

TOL_FACET = 1e-7       # on a msh reference surface (nodes were projected there)
TOL_BLADE = 1e-4       # on the bilinear S(u, v) sample (B3 gate)


def block_edges(lat):
    """{key: node id array} for every lattice edge between block corners."""
    out = {}
    for blk, (dims, V) in lat.items():
        n = V.shape
        for axis in range(3):
            o = [a for a in range(3) if a != axis]
            for i in (0, n[o[0]] - 1):
                for j in (0, n[o[1]] - 1):
                    sl = [None] * 3
                    sl[axis] = slice(None)
                    sl[o[0]], sl[o[1]] = i, j
                    ids = V[tuple(sl)].astype(int)
                    k = (min(ids[0], ids[-1]), max(ids[0], ids[-1]), len(ids))
                    if k not in out:
                        out[k] = ids if ids[0] <= ids[-1] else ids[::-1]
    return out


def catmull_rom(X, sub, alpha=0.5):
    """Centripetal Catmull-Rom through X (n, 3); sub samples per segment, nodes exact."""
    n = len(X)
    if n == 2:
        t = np.linspace(0, 1, sub + 1)[:, None]
        return X[0] + t * (X[1] - X[0])
    P = np.vstack([2 * X[0] - X[1], X, 2 * X[-1] - X[-2]])   # reflected end points
    out = [X[:1]]
    for s in range(n - 1):
        p0, p1, p2, p3 = P[s:s + 4]
        t0 = 0.0
        t1 = t0 + max(np.linalg.norm(p1 - p0), 1e-12) ** alpha
        t2 = t1 + max(np.linalg.norm(p2 - p1), 1e-12) ** alpha
        t3 = t2 + max(np.linalg.norm(p3 - p2), 1e-12) ** alpha
        for t in np.linspace(t1, t2, sub + 1)[1:]:
            A1 = (t1 - t) / (t1 - t0) * p0 + (t - t0) / (t1 - t0) * p1
            A2 = (t2 - t) / (t2 - t1) * p1 + (t - t1) / (t2 - t1) * p2
            A3 = (t3 - t) / (t3 - t2) * p2 + (t - t2) / (t3 - t2) * p3
            B1 = (t2 - t) / (t2 - t0) * A1 + (t - t0) / (t2 - t0) * A2
            B2 = (t3 - t) / (t3 - t1) * A2 + (t - t1) / (t3 - t1) * A3
            out.append(((t2 - t) / (t2 - t1) * B1 + (t - t1) / (t2 - t1) * B2)[None])
    Y = np.vstack(out)
    Y[::sub] = X                                              # nodes stay exactly
    return Y


def node_surfaces(G, bs, P, ids):
    names = [n for n in G.surfaces if n not in ("cutout", "blade")]
    on = {int(v): set() for v in ids}
    X = P[ids]
    for n in names:
        d = np.linalg.norm(G.project(n, X) - X, axis=1)
        for v, dd in zip(ids, d):
            if dd < TOL_FACET:
                on[int(v)].add(n)
    for v, x in zip(ids, X):
        q, _ = bs.foot(x)
        if np.linalg.norm(q - x) < TOL_BLADE:
            on[int(v)].add("blade")
    return on


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run")
    ap.add_argument("--blocks", default=None)
    ap.add_argument("--sub", type=int, default=8)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    run = Path(a.run)
    blocks = Path(a.blocks) if a.blocks else run / "blade" / "blocks_b3.vtk"
    out = Path(a.out) if a.out else run / "blade" / "edges_b4.npz"
    P, H, B, f2h = tfi.load_blocks(str(blocks))
    lat, missing = tfi.lattices(P, H, B, f2h, verbose=False)
    if missing:
        raise SystemExit(f"[b4] blocks without lattice: {missing}")
    G = Geometry(run)
    bs = BladeSurface(run / "blade" / "machine_scratch" / "machine_blade.npz")
    csurf = curve_surfaces(G)
    u_te, u_le = bs.edges_le_te()
    per = bs.period
    E = block_edges(lat)
    vsurf = json.loads((run / "blade" / "b3_ring.json").read_text())["vertex_surfaces"]
    on = {}
    for v in np.unique(np.concatenate(list(E.values()))):
        on[int(v)] = set(vsurf.get(str(int(v)), []))

    polys, kinds, recs = [], [], []
    for k, ids in E.items():
        common = set.intersection(*[on[int(v)] for v in ids])
        X = P[ids]
        Y = catmull_rom(X, a.sub)
        kind = "interior"
        te_le_pts, extra, u_mono = [], [], True
        if "blade" in common:
            # spline in the blade's own parameters (u unwrapped along the edge), mapped
            # by S(u, v): exactly on the blade, and u cannot run backwards (0005 folded
            # when a space spline was projected back near the TE)
            uv = np.array([bs.foot(x)[1] for x in X])
            per = bs.period
            uv[:, 0] = uv[0, 0] + np.unwrap((uv[:, 0] - uv[0, 0]) * 2 * np.pi / per) * per / (2 * np.pi)
            other = [n for n in common if n != "blade"]
            if other:
                uv[:, 1] = bs.V[0] if "hub" in other else bs.V[-1]
                kind = f"curve:blade+{other[0]}"
            else:
                kind = "surface:blade"
            UV = catmull_rom(np.c_[uv, np.zeros(len(uv))], a.sub)[:, :2]
            # a sample exactly on the TE / LE line wherever the edge crosses it, so the
            # dense polyline does not cut the sharp edge (inserted, nodes keep their slots)
            extra = []
            for s_ in range(len(UV) - 1):
                for ue in (u_te, u_le):
                    lo, hi = sorted((UV[s_, 0], UV[s_ + 1, 0]))
                    ue = ue + per * np.ceil((lo - ue) / per)
                    if lo < ue < hi:
                        f = (ue - UV[s_, 0]) / (UV[s_ + 1, 0] - UV[s_, 0])
                        extra.append((s_, f, np.array([ue, UV[s_, 1] + f * (UV[s_ + 1, 1] - UV[s_, 1])])))
            # no fold = the parameter the edge runs along (u for chordwise, v for spanwise
            # edges) is strictly monotone
            dU = np.diff(UV, axis=0)
            ax = int(np.argmax(np.abs(UV[-1] - UV[0]) / np.r_[per, bs.V[-1] - bs.V[0]]))
            u_mono = bool(np.all(dU[:, ax] > 0) or np.all(dU[:, ax] < 0))
            Y = bs.at(UV[:, 0], np.clip(UV[:, 1], bs.V[0], bs.V[-1]))
            te_le_pts = [bs.at(np.array([e[2][0]]), np.array([np.clip(e[2][1], bs.V[0], bs.V[-1])]))[0] for e in extra]
        elif len(common) >= 2:
            cs = [t for t, ss in csurf.items() if common <= ss or len(common & ss) >= 2]
            cs = [t for t in cs if len(common & csurf[t]) >= 2]
            if cs:
                kind = "curve:" + "+".join(sorted(common))
                for i in range(len(Y)):
                    cand = np.array([G.curves[t].project(Y[i])[0] for t in cs])
                    Y[i] = cand[np.argmin(np.linalg.norm(cand - Y[i], axis=1))]
        elif len(common) == 1:
            n = next(iter(common))
            kind = "surface:" + n
            Y = G.project(n, Y)
        Y[::a.sub] = X                                         # lattice nodes untouched
        # residual of the samples on their constraint
        res, fac = 0.0, 0.0
        if kind.startswith("surface:") or kind.startswith("curve:"):
            for n in (kind.split(":")[1].split("+")):
                if n == "blade":
                    d = [np.linalg.norm(bs.foot(y)[0] - y) for y in Y]
                else:
                    d = np.linalg.norm(G.project(n, Y) - Y, axis=1)
                res = max(res, float(np.max(d)))
            # what part of it is only the facetting of the msh hub / shroud (the analytic
            # blade root lies on the true cylinder, the msh facets ~h^2/8r inside it)
            fac = 0.0
            for n in (kind.split(":")[1].split("+")):
                if n != "blade":
                    fac = max(fac, float(np.max(np.linalg.norm(G.project(n, X) - X, axis=1))))
        # chord deviation + fold check per lattice segment
        dev, fold = 0.0, False
        for s in range(len(ids) - 1):
            if kind.endswith("blade") or "blade+" in kind:
                fold = fold or not u_mono                 # around a sharp edge the chord test is meaningless
                break
            seg = Y[s * a.sub:(s + 1) * a.sub + 1]
            c0, c1 = X[s], X[s + 1]
            t = c1 - c0
            L = np.linalg.norm(t)
            proj = (seg - c0) @ (t / L)
            dev = max(dev, float(np.max(np.linalg.norm((seg - c0) - np.outer(proj, t / L), axis=1)) / L))
            fold |= bool(np.any(np.diff(proj) <= 0))
        if ("blade" in kind) and te_le_pts:
            Y = Y.copy()
            ins = sorted(extra, key=lambda e: e[0] + e[1])
            pieces, last = [], 0
            for (s_, f, _), q in zip(ins, [p_ for _, p_ in sorted(zip([e[0] + e[1] for e in extra], te_le_pts))]):
                pieces += [Y[last:s_ + 1], q[None]]
                last = s_ + 1
            Y = np.vstack(pieces + [Y[last:]])
        polys.append(Y)
        kinds.append(kind)
        recs.append({"nodes": [int(ids[0]), int(ids[-1])], "n": int(len(ids)), "kind": kind,
                     "residual": res, "node_offset_to_msh_facets": fac if kind != "interior" else 0.0,
                     "chord_dev_rel": dev, "fold": fold})
    kinds_count = {}
    for k in kinds:
        kinds_count[k] = kinds_count.get(k, 0) + 1
    res_max = {}
    for r in recs:
        res_max[r["kind"]] = max(res_max.get(r["kind"], 0.0), r["residual"])
    rep = {"blocks": str(blocks), "edges": len(recs), "sub": a.sub, "kinds": kinds_count,
           "residual_max_per_kind": res_max,
           "chord_dev_rel_max": max(r["chord_dev_rel"] for r in recs),
           "folded_edges": sum(r["fold"] for r in recs),
           "gate": {"no_fold": not any(r["fold"] for r in recs),
                    "residual_lt_1e-4": max(r["residual"] for r in recs) < 1e-4},
           "edges_detail": recs}
    off = np.cumsum([0] + [len(y) for y in polys])
    np.savez(out, polyline=np.vstack(polys), offset=off,
             nodes=np.array([np.r_[r["nodes"], r["n"]] for r in recs]),
             edge_ids=np.array([np.pad(E[k], (0, max(len(v) for v in E.values()) - len(E[k])), constant_values=-1)
                                for k in E]),
             kind=np.array(kinds))
    import meshio
    pts = np.vstack(polys)
    lines = np.vstack([np.c_[np.arange(off[i], off[i + 1] - 1), np.arange(off[i] + 1, off[i + 1])]
                       for i in range(len(polys))])
    kid = {k: i for i, k in enumerate(sorted(kinds_count))}
    cd = np.concatenate([[kid[kinds[i]]] * (off[i + 1] - off[i] - 1) for i in range(len(polys))]).astype(float)
    meshio.write(str(out.with_suffix(".vtk")), meshio.Mesh(pts, [("line", lines)], cell_data={"kind": [cd]}))
    (out.with_suffix(".json")).write_text(json.dumps(rep, indent=1))
    rep_s = {k: v for k, v in rep.items() if k != "edges_detail"}
    rep_s["kind_ids"] = kid
    print(json.dumps(rep_s, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
