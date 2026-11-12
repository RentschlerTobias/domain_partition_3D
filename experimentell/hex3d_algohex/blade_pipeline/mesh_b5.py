"""B5: fine mesh of the core blocks, edge -> face -> volume, every entity once.

    mesh_b5.py <run_dir> [--target-h 0.05] [--first 0.0009] [--growth 1.2]

The coarse lattice (blade/blocks_b3.vtk) only provides TOPOLOGY. Geometry comes
from the exact surfaces and curves:
  edges   the B4 curves (blade/edges_b4.npz: on the dtOO feature curves /
          surfaces), resampled by arc length; span edges that start or end on
          hub / shroud get the BL distribution with the absolute first cell
          --first at that end (per edge, so y1 is the same everywhere);
  faces   Coons patch of the 4 discretised edges in cylinder coordinates
          (r, theta, z), then put onto the face's surface: periodic faces by
          theta = theta_per(r, z) (circumferential, r and z untouched), the other
          surfaces by closest point (a small correction: hub / shroud / inlet /
          outlet are revolution surfaces / planes, exact in (r, theta, z) up to
          the Coons blend); interior faces stay the Coons patch;
  volume  Gordon-Hall TFI of the block's 6 finished faces.
Each edge and face is computed once and handed to every block that owns it
(in that block's orientation), so the mesh is conformal by construction.
Divisions: tfi.direction_classes + tfi.solve_block_divisions, BL classes get
the extra layers (core_tfi.class_cells). Ring blocks (ids >= core count) are B6;
the faces they share with the core (cutout) are stored in blade/b5_entities.npz.
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
import quality as qual                                                # noqa: E402
from core_tfi import graded, tfi_param                                # noqa: E402
from geometry import Geometry                                         # noqa: E402
from winslow import winslow, winslow2d                                # noqa: E402
from ring import BladeSurface, blend_segment                         # noqa: E402

WALLS = ("hub", "shroud", "blade")
# reference surface -> OpenFOAM patch, as dtOO names them (xml/ru_gridMesh.xml): the
# channel region R_0 has its faces in dtGmshModel order SOUTH, NORTH, FRONT, BACK,
# WEST, EAST = msh entities 1..6 = hub, shroud, inlet, outlet, periodic_B, periodic_A;
# R_0->WEST is RU_PERIOPS, R_0->EAST is RU_PERIOSS
PATCHES = {"hub": "RU_HUB", "shroud": "RU_SHROUD", "blade": "RU_BLADE", "inlet": "RU_INLET",
           "outlet": "RU_OUTLET", "periodic_B": "RU_PERIOPS", "periodic_A": "RU_PERIOSS"}


def cyl(X):
    return np.stack([np.hypot(X[..., 0], X[..., 1]), np.arctan2(X[..., 1], X[..., 0]), X[..., 2]], -1)


def cart(C, th0=0.0):
    return np.stack([C[..., 0] * np.cos(C[..., 1] + th0), C[..., 0] * np.sin(C[..., 1] + th0), C[..., 2]], -1)


def resample(poly, f):
    """Points at arc-length fractions f along a polyline."""
    L = np.r_[0, np.cumsum(np.linalg.norm(np.diff(poly, axis=0), axis=1))]
    t = f * L[-1]
    return np.stack([np.interp(t, L, poly[:, k]) for k in range(3)], -1), float(L[-1])


def resample_index(poly, nodes, n):
    """n cells along a polyline that passes through the coarse `nodes`: uniform in the
    coarse index, arc length inside every coarse segment."""
    L = np.r_[0, np.cumsum(np.linalg.norm(np.diff(poly, axis=0), axis=1))]
    at = [int(np.argmin(np.linalg.norm(poly - x, axis=1))) for x in nodes]
    at[0], at[-1] = 0, len(poly) - 1
    nc = len(nodes) - 1
    out = []
    for s in np.linspace(0, nc, n + 1):
        k = min(int(np.floor(s)), nc - 1)
        t = s - k
        a0, a1 = L[at[k]], L[at[k + 1]]
        out.append([np.interp(a0 + t * (a1 - a0), L, poly[:, q]) for q in range(3)])
    return np.array(out)


def coons(B, T, L, R):
    """Coons patch from edges B (s, j=0), T (s, j=m), L (i=0, t), R (i=n, t), all as
    arrays of points; blending parameters = mean arc-length fractions of the two
    opposite edges (the boundary is reproduced exactly for any choice)."""
    def frac(E):
        d = np.r_[0, np.cumsum(np.linalg.norm(np.diff(E, axis=0), axis=1))]
        return d / d[-1]
    u = 0.5 * (frac(B) + frac(T))[:, None, None]
    v = 0.5 * (frac(L) + frac(R))[None, :, None]
    P = ((1 - v) * B[:, None] + v * T[:, None] + (1 - u) * L[None] + u * R[None]
         - (1 - u) * (1 - v) * B[0] - u * (1 - v) * B[-1] - (1 - u) * v * T[0] - u * v * T[-1])
    return P


def class_cells(lat, classes, counts, vs, P, first, growth):
    """Cells per class incl. BL layers at hub / shroud ends (as in core_tfi)."""
    out = {}
    for ci, c in enumerate(classes):
        walls = set()
        L = []
        for r, ax in c:
            vert = lat[r][1]
            for side in (0, 1):
                g = tfi.face_grid(vert, ax, side)
                common = set.intersection(*[vs.get(int(v), set()) for v in g.ravel()])
                walls |= common & set(WALLS)
            L.append(tfi.axis_length(P, vert, ax))
        n = int(counts[ci])
        if walls:
            nbl = int(np.ceil(np.log(max(min(L) / n / first, 1.0)) / np.log(growth)))
            n += nbl * min(len(walls), 2)
        out[ci] = {"cells": n, "walls": sorted(walls)}
    return out


def strip_ok(cut, wall, t):
    """All quads (cut_i, cut_i+1, wall_i+1, wall_i) of one span row turn the same way
    about the local span direction t (n, 3): no crossing pairing lines, in 3D (the
    blade-to-blade projection with one mean radius is too coarse on the fine level)."""
    Q = np.stack([cut[:-1], cut[1:], wall[1:], wall[:-1]], 1)
    tm = 0.5 * (t[:-1] + t[1:])
    e1 = np.roll(Q, -1, 1) - Q
    e0 = Q - np.roll(Q, 1, 1)
    cr = np.einsum("qkc,qc->qk", np.cross(e0, e1), tm)
    sg = np.sign(np.median(cr))
    return bool((np.sign(cr) == sg).all())


def match(stored, wanted):
    """The dihedral transform (as a function on 2D arrays) taking the stored
    coarse id grid onto the wanted one."""
    for tr in (False, True):
        A = stored.T if tr else stored
        for fa in (False, True):
            for fb in (False, True):
                B = A[::-1] if fa else A
                B = B[:, ::-1] if fb else B
                if B.shape == wanted.shape and np.array_equal(B, wanted):
                    def f(X, tr=tr, fa=fa, fb=fb):
                        Y = np.swapaxes(X, 0, 1) if tr else X
                        Y = Y[::-1] if fa else Y
                        return Y[:, ::-1] if fb else Y
                    return f
    raise RuntimeError("face orientation not found")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run")
    ap.add_argument("--target-h", type=float, default=0.05)
    ap.add_argument("--first", type=float, default=0.0009)
    ap.add_argument("--growth", type=float, default=1.2)
    ap.add_argument("--winslow", type=int, default=300, help="Winslow iterations per block (0 = off)")
    ap.add_argument("--winslow2d", type=int, default=300,
                    help="2D Winslow iterations on boundary faces, on the surface (0 = off)")
    ap.add_argument("--ring", action="store_true", help="B6: mesh the blade ring blocks too")
    ap.add_argument("--min-spacing", type=float, default=0.3,
                    help="ring wall: smallest chord spacing as fraction of the row mean")
    ap.add_argument("--dump-block", type=int, nargs="*", default=None, help="diagnostic: save block arrays")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    run = Path(a.run)
    bd = run / "blade"
    out = Path(a.out) if a.out else bd / ("mesh_b6.vtk" if a.ring else "mesh_b5_core.vtk")
    b3 = json.loads((bd / "b3_ring.json").read_text())
    vs = {int(k): set(v) for k, v in b3["vertex_surfaces"].items()}
    ncore = int(b3["core_blocks"])
    P, H, B, f2h = tfi.load_blocks(str(bd / "blocks_b3.vtk"))
    lat, missing = tfi.lattices(P, H, B, f2h, verbose=False)
    classes = tfi.direction_classes(lat, f2h, H, B, verbose=False)
    classes = [list(c) for c in (classes.values() if isinstance(classes, dict) else classes)]
    counts = tfi.solve_block_divisions(lat, classes, P, a.target_h, verbose=False)
    cof = tfi.class_of_axis(classes)
    cc = class_cells(lat, classes, counts, vs, P, a.first, a.growth)
    G = Geometry(run)
    bs = BladeSurface(bd / "machine_scratch" / "machine_blade.npz")
    blade_pts = []
    # the cutout in its dtOO O-grid parameters (s around the blade, t hub -> shroud)
    Sc = G.cutout_grid()
    if Sc is None:
        raise SystemExit("[b5] cutout is not a structured quad grid")
    # parameters = mean arc length (the O-grid rows cluster at hub / shroud: by row
    # index, a straight line in (s, t) would bend flat onto hub and shroud)
    Uc = np.r_[0, np.cumsum(np.linalg.norm(np.diff(Sc, axis=0), axis=-1).mean(1))]
    Vc = np.r_[0, np.cumsum(np.linalg.norm(np.diff(Sc, axis=1), axis=-1).mean(0))]
    cs = BladeSurface(S=Sc, U=Uc, V=Vc / Vc[-1])

    def ring_sides(r):
        """(wall axis, blade-side grid, cutout-side grid) of ring block r."""
        vert = lat[r][1]
        for axl in range(3):
            sides = [tfi.face_grid(vert, axl, s_) for s_ in (0, 1)]
            onb = [all("blade" in vs.get(int(v), ()) for v in g.ravel()) for g in sides]
            if onb[0] != onb[1]:
                return (axl, sides[0], sides[1]) if onb[0] else (axl, sides[1], sides[0])
        raise RuntimeError(f"ring block {r}: no side on the blade")

    # the core's cutout faces = the outer faces of the ring blocks, and their edges
    cut_faces, cut_edges = set(), set()
    for r in lat:
        if r >= ncore:
            gc = ring_sides(r)[2]
            cut_faces.add(frozenset(int(x) for x in gc.ravel()))
            for ids_ in (gc[:, 0], gc[:, -1], gc[0, :], gc[-1, :]):
                cut_edges.add(tuple(sorted(int(x) for x in ids_)))

    def param_coons(surf, E4):
        """Coons patch of 4 boundary point rows in the (u, v) of a structured surface,
        mapped by it: exact on the surface, u unwrapped (periodic)."""
        uvs = [np.array([surf.foot(x)[1] for x in E]) for E in E4]
        per = surf.period
        u0 = uvs[0][0, 0]
        for uv in uvs:
            uv[:, 0] = u0 + np.mod(uv[:, 0] - u0 + per / 2, per) - per / 2
            uv[:, 0] = uv[0, 0] + np.unwrap((uv[:, 0] - uv[0, 0]) * 2 * np.pi / per) * per / (2 * np.pi)
        z = lambda uv: np.c_[uv, np.zeros(len(uv))]
        UV = coons(*(z(uv) for uv in uvs))[..., :2]
        return surf.at(UV[..., 0], np.clip(UV[..., 1], surf.V[0], surf.V[-1]))

    def cut_poly(ids0, per_seg=200):
        """Dense curve on the cutout through the coarse nodes ids0: monotone cubic in
        the cutout's (s, t) over the coarse index, mapped by S(s, t)."""
        from scipy.interpolate import PchipInterpolator
        st = np.array([cs.foot(P[v])[1] for v in ids0])
        st[:, 0] = st[0, 0] + np.unwrap((st[:, 0] - st[0, 0]) * 2 * np.pi / cs.period) * cs.period / (2 * np.pi)
        k = np.arange(len(ids0), dtype=float)
        tau = np.linspace(0, k[-1], per_seg * (len(ids0) - 1) + 1)
        if len(ids0) > 2:
            s_, t_ = PchipInterpolator(k, st[:, 0])(tau), PchipInterpolator(k, st[:, 1])(tau)
        else:
            s_, t_ = np.interp(tau, k, st[:, 0]), np.interp(tau, k, st[:, 1])
        X = cs.at(s_, np.clip(t_, cs.V[0], cs.V[-1]))
        X[::per_seg] = P[list(ids0)]
        return X

    # B4 curves by their coarse node sequence
    e4 = np.load(bd / "edges_b4.npz")
    b4 = {}
    for i in range(len(e4["nodes"])):
        ids = e4["edge_ids"][i]
        ids = tuple(int(x) for x in ids[ids >= 0])
        b4[tuple(sorted(ids))] = (ids, e4["polyline"][e4["offset"][i]:e4["offset"][i + 1]], str(e4["kind"][i]))

    # ---------------------------------------------------------------- edges
    edges = {}
    erep = {"max_proj": 0.0}
    # the core is built against the B4 cutout ("old": the coarse cutout, valid by
    # B2 / B3) and then morphed onto the exact cutout ("true"), which the ring uses
    mode = {"cut": "old"}

    def edge(ids, n):
        """Fine points of the edge with coarse ids `ids` (that orientation), n cells."""
        key0 = tuple(sorted(int(x) for x in ids))
        key = (key0, mode["cut"]) if key0 in cut_edges else key0
        if key not in edges:
            ids0, poly, kind = b4[key0]
            if key0 in cut_edges and mode["cut"] == "true":
                poly, kind = cut_poly(ids0), "surface:cutout"
            ends = [vs.get(ids0[0], set()), vs.get(ids0[-1], set())]
            common = set.intersection(*[vs.get(int(v), set()) for v in ids0])
            L_ = float(np.sum(np.linalg.norm(np.diff(poly, axis=0), axis=1)))
            lo = any(w in ends[0] and w not in common for w in WALLS)
            hi = any(w in ends[1] and w not in common for w in WALLS)
            if lo or hi:
                f = graded(n, a.first / L_ if lo else 0, a.first / L_ if hi else 0, a.growth)
                X, _ = resample(poly, f)
            else:
                # by coarse index, arc length only inside each coarse segment: keeps the
                # coarse pairing of opposite edges (wall vs cutout around the TE cap,
                # made valid in B3) instead of matching two unrelated arc lengths
                X = resample_index(poly, P[list(ids0)], n)
                d_ = np.r_[0, np.cumsum(np.linalg.norm(np.diff(X, axis=0), axis=1))]
                f = d_ / d_[-1]
            X[0], X[-1] = P[ids0[0]], P[ids0[-1]]
            # back onto the constraint (the polyline chords between B4 samples)
            names = [x for x in kind.split(":")[1].split("+")] if ":" in kind else []
            Y = X.copy()
            if names and "blade" not in names and "cutout" not in names:
                for _ in range(6 if len(names) == 2 else 1):
                    for nm in sorted(names, key=lambda s: s.startswith("periodic")):
                        Y[1:-1] = G.project_any(nm, Y[1:-1])
                erep["max_proj"] = max(erep["max_proj"], float(np.abs(Y - X).max()))
            edges[key] = (tuple(ids0), Y, f)
        ids0, Y, f = edges[key]
        return (Y, f) if tuple(int(x) for x in ids) == ids0 else (Y[::-1], 1 - f[::-1])

    # ---------------------------------------------------------------- faces
    faces = {}
    frep = {}

    def face_key(g):
        key0 = frozenset(int(x) for x in g.ravel())
        touch = key0 in cut_faces or any(tuple(sorted(int(x) for x in e_)) in cut_edges
                                         for e_ in (g[:, 0], g[:, -1], g[0, :], g[-1, :]))
        return (key0, mode["cut"]) if touch else key0

    def face(g, n0, n1):
        """Fine (n0+1, n1+1, 3) points of the face with coarse id grid g."""
        key = face_key(g)
        key0 = key[0] if isinstance(key, tuple) else key
        if key not in faces:
            Bt, _ = edge(g[:, 0], n0)
            Tp, _ = edge(g[:, -1], n0)
            Lf, _ = edge(g[0, :], n1)
            Rt, _ = edge(g[-1, :], n1)
            th0 = float(np.angle(np.exp(1j * np.arctan2(P[g.ravel()][:, 1], P[g.ravel()][:, 0])).mean()))
            rel = lambda X: (lambda C: np.stack([C[..., 0], np.angle(np.exp(1j * (C[..., 1] - th0))), C[..., 2]], -1))(cyl(X))
            F = cart(coons(rel(Bt), rel(Tp), rel(Lf), rel(Rt)), th0)
            common = set.intersection(*[vs.get(int(v), set()) for v in g.ravel()])
            # the cutout is a surface for the core (its outer face towards the ring)
            if key0 in cut_faces and mode["cut"] == "old":
                faces[key] = (g.copy(), F, "cutout")
                return match(g, g)(F)
            if key0 in cut_faces:
                Fc = param_coons(cs, (Bt, Tp, Lf, Rt))
                frep["cutout"] = max(frep.get("cutout", 0.0), float(np.linalg.norm(Fc - F, axis=-1).max()))
                Fc[0], Fc[-1], Fc[:, 0], Fc[:, -1] = F[0], F[-1], F[:, 0], F[:, -1]
                faces[key] = (g.copy(), Fc, "cutout")
                g0, F, name = faces[key]
                return match(g0, g)(F)
            name = sorted(common)[0] if len(common) == 1 else None
            if name == "blade":
                # on the blade: Coons in the dtOO parameters (u, v) of the 4 edges and
                # S(u, v) -- exact, and no closest-point jump across the sharp TE
                Fb = param_coons(bs, (Bt, Tp, Lf, Rt))
                frep["blade"] = max(frep.get("blade", 0.0), float(np.linalg.norm(Fb[1:-1, 1:-1] - F[1:-1, 1:-1], axis=-1).max()) if F.shape[0] > 2 and F.shape[1] > 2 else 0.0)
                Fb[0], Fb[-1], Fb[:, 0], Fb[:, -1] = F[0], F[-1], F[:, 0], F[:, -1]
                F = Fb
                blade_pts.append(F.reshape(-1, 3))
                faces[key] = (g.copy(), F, name)
                g0, F, name = faces[key]
                return match(g0, g)(F)
            inner = F[1:-1, 1:-1]
            if name is not None and inner.size:
                Y = G.project_any(name, inner.reshape(-1, 3)).reshape(inner.shape)
                frep[name] = max(frep.get(name, 0.0), float(np.linalg.norm(Y - inner, axis=-1).max()))
                F[1:-1, 1:-1] = Y
            if a.winslow2d > 0 and name is not None and inner.size:
                F, _ = winslow2d(F, project=lambda X, nm=name: G.project_any(nm, X), iters=a.winslow2d)
            faces[key] = (g.copy(), F, name)
        g0, F, name = faces[key]
        return match(g0, g)(F)

    # ---------------------------------------------------------------- ring wall (B6)
    te_le = bs.edges_le_te() if a.ring else ()
    rrep = {}

    def ring_wall(r, n):
        """Blade face of ring block r from its finished cutout face, as B3 on the coarse
        level: every cutout node of a span row goes to its foot point (u, v) on the
        blade; u is blended towards the cutout-arc-length map with the smallest weight
        that keeps the spacing and every (cutout, wall) quad uncrossed (strip_ok).
        Where the TE / LE lies inside the row, the cap apex (the
        cutout node where the foot-point sequence passes it) gets a wall node
        exactly on it, and the two parts are blended separately. The wall chord edges
        on hub / shroud and the blade face are stored as the shared entities."""
        axl, gw, gc = ring_sides(r)
        o = [x for x in range(3) if x != axl]
        rr = np.hypot(P[gw][..., 0], P[gw][..., 1])
        tr = np.ptp(rr, axis=0).mean() > np.ptp(rr, axis=1).mean()      # axis 0 = span
        gws, gcs = (gw.T, gc.T) if tr else (gw, gc)                     # (chord, span)
        nc, ns = (n[o[1]], n[o[0]]) if tr else (n[o[0]], n[o[1]])
        C = face(gcs, nc, ns)                                          # cutout, fixed
        EL, _ = edge(gws[0, :], ns)
        ER, _ = edge(gws[-1, :], ns)
        per = bs.period
        near = lambda u, ref: ref + np.mod(u - ref + per / 2, per) - per / 2
        UV = np.array([[bs.foot(x)[1] for x in C[:, j]] for j in range(ns + 1)])   # (span, chord, 2)
        uvL = np.array([bs.foot(x)[1] for x in EL])
        uvR = np.array([bs.foot(x)[1] for x in ER])
        hs = [all(w in vs.get(int(v), ()) for v in gws[:, k]) for w, k in (("hub", 0), ("shroud", -1))]
        T = np.gradient(C, axis=1)                                     # span direction
        T /= np.linalg.norm(T, axis=-1, keepdims=True)
        # pass 1, per span row: foot u (unwrapped, ends = the two span edges), v, the
        # cutout arc-length fractions and where the foot sequence passes TE / LE
        rows = []
        for j in range(ns + 1):
            uf = UV[j, :, 0].copy()
            uf[0] = uvL[j, 0]
            for i in range(1, nc + 1):
                uf[i] = near(uf[i], uf[i - 1])
            uf[-1] = near(uvR[j, 0], uf[-2])
            sgn = 1.0 if uf[-1] > uf[0] else -1.0
            seg = np.linalg.norm(np.diff(C[:, j], axis=0), axis=1)
            s = np.r_[0, np.cumsum(seg)] / seg.sum()
            v = (1 - s) * uvL[j, 1] + s * uvR[j, 1]
            if j == 0 and hs[0]:
                v[:] = bs.V[0]
            if j == ns and hs[1]:
                v[:] = bs.V[-1]
            cross = {}
            for e, ue in enumerate(te_le):
                ue = near(ue, 0.5 * (uf[0] + uf[-1]))
                if sgn * (ue - uf[0]) > 0 and sgn * (uf[-1] - ue) > 0:
                    i1 = int(np.argmax(sgn * (uf - ue) > 0))
                    cross[e] = (min(max(i1 - 1 if abs(uf[i1 - 1] - ue) < abs(uf[i1] - ue) else i1, 1), nc - 1), ue)
            rows.append((uf, sgn, s, v, cross))
        # one TE / LE chord index per block, only where the edge lies inside every
        # span row: a wall grid line along the whole edge. Where the ring split runs
        # past the edge over the span, a node line ending half way would kink the
        # wall grid; the blend alone keeps those rows uncrossed.
        kedge = {e: int(round(np.median([rw[4][e][0] for rw in rows])))
                 for e in range(len(te_le)) if all(e in rw[4] for rw in rows)}
        apex = sorted(kedge.values())

        def row_parts(j):
            uf, sgn, s, v, cross = rows[j]
            uf = uf.copy()
            cuts = [0]
            for e, (_, ue) in cross.items():
                if e not in kedge:
                    continue
                uf[kedge[e]] = ue
                cuts.append(kedge[e])
            cuts = sorted(set(cuts)) + [nc]
            return uf, sgn, s, v, list(zip(cuts[:-1], cuts[1:]))

        def blend(j, k0, k1, w=None):
            """u of row j on part (k0, k1): smallest valid weight, or the weight w."""
            uf, sgn, s, v, _ = row_parts(j)
            sub = s[k0:k1 + 1]
            up = uf[k0] + (sub - sub[0]) / (sub[-1] - sub[0]) * (uf[k1] - uf[k0])
            ok = lambda uu: strip_ok(C[k0:k1 + 1, j], bs.at(uu, np.clip(v[k0:k1 + 1], bs.V[0], bs.V[-1])),
                                     T[k0:k1 + 1, j])
            if w is None:
                w, uu = blend_segment(sgn * uf[k0:k1 + 1], sgn * up, a.min_spacing, valid=lambda x: ok(sgn * x))
                uu = sgn * uu
            else:
                uu = (1 - w) * uf[k0:k1 + 1] + w * up
            d = sgn * np.diff(uu)
            good = bool((d > 0).all() and d.min() >= a.min_spacing * d.mean() and ok(uu))
            return w, uu, good

        # pass 2: smallest valid weight per row and part; pass 3: one weight per part
        # (the largest), so u varies smoothly over the span; rows where that weight
        # is not valid keep their own
        wpart = {}
        for j in range(ns + 1):
            for k0, k1 in row_parts(j)[4]:
                w, _, good = blend(j, k0, k1)
                wpart[(k0, k1)] = max(wpart.get((k0, k1), 0.0), w if good else 0.0)
        Uw = np.zeros((ns + 1, nc + 1))
        Vw = np.zeros((ns + 1, nc + 1))
        own, bad = 0, 0
        for j in range(ns + 1):
            uf, sgn, s, v, parts = row_parts(j)
            u = uf.copy()
            for k0, k1 in parts:
                _, uu, good = blend(j, k0, k1, wpart[(k0, k1)])
                if not good:
                    _, uu, good = blend(j, k0, k1)
                    own += 1
                    bad += 0 if good else 1
                u[k0:k1 + 1] = uu
            Uw[j], Vw[j] = u, v
        wmax = max(wpart.values())
        W = bs.at(Uw, np.clip(Vw, bs.V[0], bs.V[-1])).transpose(1, 0, 2)  # (chord, span, 3)
        W[0], W[-1] = EL, ER
        for k, key_ids in ((0, gws[:, 0]), (-1, gws[:, -1])):
            row = W[:, k]
            d_ = np.r_[0, np.cumsum(np.linalg.norm(np.diff(row, axis=0), axis=1))]
            edges[tuple(sorted(int(x) for x in key_ids))] = (tuple(int(x) for x in key_ids), row, d_ / d_[-1])
        faces[frozenset(int(x) for x in gw.ravel())] = (gws.copy(), W, "blade")
        blade_pts.append(W.reshape(-1, 3))
        rrep[int(r)] = {"u_blend_max": float(wmax), "apex_index": sorted(set(apex)),
                        "chord_cells": int(nc), "span_cells": int(ns),
                        "weight_per_part": {f"{k0}-{k1}": float(w) for (k0, k1), w in wpart.items()},
                        "parts_with_own_weight": own, "parts_without_valid_blend": bad}

    # ---------------------------------------------------------------- volumes
    def build_block(r):
        """TFI of block r from its 6 finished faces; also (axis, side, grid) per side."""
        dims, vert = lat[r]
        n = [cc[cof[(r, ax)]]["cells"] for ax in range(3)]
        if r >= ncore:
            ring_wall(r, n)
        X = np.full((n[0] + 1, n[1] + 1, n[2] + 1, 3), np.nan)
        sides = []
        for ax in range(3):
            o = [x for x in range(3) if x != ax]
            for side in (0, 1):
                g = tfi.face_grid(vert, ax, side)
                F = face(g, n[o[0]], n[o[1]])
                sl = [slice(None)] * 3
                sl[ax] = 0 if side == 0 else -1
                X[tuple(sl)] = F
                sides.append((tuple(sl), g, faces[face_key(g)][2], (n[o[0]], n[o[1]])))
        # blending parameters per axis: mean fractions of the 4 parallel edges
        par = []
        for ax in range(3):
            fs = []
            for i in (0, -1):
                for j in (0, -1):
                    sl = [i, j]
                    sl.insert(ax, slice(None))
                    _, f = edge(vert[tuple(sl)], n[ax])
                    fs.append(f)
            par.append(np.mean(fs, axis=0))
        Xi = tfi_param(np.nan_to_num(X), *par)
        bnd = ~np.isnan(X[..., 0])
        Xi[bnd] = X[bnd]
        if a.dump_block is not None and int(r) in a.dump_block:
            np.save(bd / f"dump_block{int(r)}_tfi.npy", Xi)
        return Xi, (sides, par)

    def retfi(Xm, par):
        """Block interior again by TFI from its (morphed) 6 sides: the morph places
        the faces, TFI the interior (Laplace alone lets a cell row next to a strongly
        moved face cross)."""
        Xb = np.full_like(Xm, np.nan)
        for ax in range(3):
            for sd in (0, -1):
                sl = [slice(None)] * 3
                sl[ax] = sd
                Xb[tuple(sl)] = Xm[tuple(sl)]
        Xi = tfi_param(np.nan_to_num(Xb), *par)
        bnd = ~np.isnan(Xb[..., 0])
        Xi[bnd] = Xb[bnd]
        return Xi

    def morph_core(blk, sides):
        """Move the core onto the exact cutout: the displacement true - old on the cutout
        faces, extended over the welded core by a Laplace problem with edge weights
        1 / length^2 (physical, so a BL stack moves almost rigidly with its wall instead
        of being squeezed index by index). Components in the cylinder frame (dr, r dtheta,
        dz) with the constraint of each boundary surface: dr = 0 on hub / shroud, dz = 0
        on inlet / outlet, periodic nodes go back circumferentially onto their surface.
        Every node of a cutout face ends exactly on the ring's face."""
        from scipy.sparse import coo_matrix, diags
        from scipy.sparse.linalg import cg
        rs = sorted(blk)
        offs, pts_, hx, base_ = {}, [], [], 0
        for r in rs:
            ni, nj, nk = blk[r].shape[:3]
            offs[r] = np.arange(ni * nj * nk).reshape(ni, nj, nk) + base_
            pts_.append(blk[r].reshape(-1, 3))
            hx.append(tfi.block_cells(offs[r]))
            base_ += ni * nj * nk
        Pw_, remap_ = tfi.weld(np.vstack(pts_), 1e-9)
        Hw_ = remap_[np.vstack(hx)]
        nw = len(Pw_)
        C0 = cyl(Pw_)
        Dc = np.zeros((nw, 3))                                  # (dr, r dtheta, dz)
        fix = np.zeros((nw, 3), bool)
        surf = {}
        mode["cut"] = "true"
        for r in rs:
            for sl, g, name, (n0, n1) in sides[r][0]:
                ids_ = remap_[offs[r][sl]].ravel()
                if name == "cutout":
                    Ct = cyl(face(g, n0, n1).reshape(-1, 3))
                    Co = C0[ids_]
                    Dc[ids_] = np.c_[Ct[:, 0] - Co[:, 0], Co[:, 0] * np.angle(np.exp(1j * (Ct[:, 1] - Co[:, 1]))),
                                     Ct[:, 2] - Co[:, 2]]
                    fix[ids_] = True
                elif name is not None:
                    for v in ids_:
                        surf.setdefault(int(v), set()).add(name)
        mode["cut"] = "old"
        cons = {"hub": 0, "shroud": 0, "inlet": 2, "outlet": 2}
        for v, nm in surf.items():
            for x in nm:
                if x in cons:
                    fix[v, cons[x]] = True
            if sum(x.startswith("periodic") for x in nm) and len(nm) > 1:
                fix[v] = True                                   # periodic curves stay
        i_, j_ = (Hw_[:, [0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3]].ravel(),
                  Hw_[:, [1, 2, 3, 0, 5, 6, 7, 4, 4, 5, 6, 7]].ravel())
        w_ = 1.0 / np.maximum(((Pw_[i_] - Pw_[j_]) ** 2).sum(1), 1e-16)
        A = coo_matrix((w_, (i_, j_)), shape=(nw, nw)).tocsr()
        A = A + A.T
        A.data[:] = np.maximum(A.data, 0)
        L = (diags(np.asarray(A.sum(1)).ravel()) - A).tocsr()
        info = []
        for k in range(3):
            fr = np.nonzero(~fix[:, k])[0]
            fx = np.nonzero(fix[:, k])[0]
            Lff = L[fr][:, fr].tocsr()
            rhs = -(L[fr][:, fx] @ Dc[fx, k])
            x, it = cg(Lff, rhs, M=diags(1.0 / Lff.diagonal()), rtol=1e-9, maxiter=20000)
            Dc[fr, k] = x
            info.append(int(it))
        C1 = np.c_[C0[:, 0] + Dc[:, 0], C0[:, 1] + Dc[:, 1] / C0[:, 0], C0[:, 2] + Dc[:, 2]]
        X = cart(C1)
        per = {}
        for v, nm in surf.items():
            if not fix[v].all():
                for x in nm:
                    if x.startswith("periodic"):
                        per.setdefault(x, []).append(v)
        for nm, vv in per.items():
            vv = np.array(vv)
            X[vv] = G.project_any(nm, X[vv])
        q0 = qual.summary(Pw_, Hw_)
        q1 = qual.summary(X, Hw_)
        mrep = {"displacement_max": float(np.linalg.norm(X - Pw_, axis=1).max()), "cg_info": info,
                "inverted_before": int(q0["inverted"]), "inverted_after": int(q1["inverted"]),
                "min_sj_after": float(q1["min_sj"])}
        return {r: X[remap_[offs[r]]] for r in rs}, mrep

    core, sides = {}, {}
    for r in sorted(lat):
        if r < ncore:
            core[r], sides[r] = build_block(r)
    core, mrep = morph_core(core, sides)
    core = {r: retfi(core[r], sides[r][1]) for r in core}
    cq = [qual.summary(core[r].reshape(-1, 3), tfi.block_cells(np.arange(core[r][..., 0].size).reshape(core[r].shape[:3])))
          for r in sorted(core)]
    mrep["inverted_after_retfi"] = int(sum(c["inverted"] for c in cq))
    mode["cut"] = "true"
    blocks = dict(core)
    if a.ring:
        for r in sorted(lat):
            if r >= ncore:
                blocks[r], sides[r] = build_block(r)
    pts, hexes, bids = [], [], []
    bq, bpatch = [], []
    wrep = {}
    base = 0
    for r in sorted(blocks):
        Xi = blocks[r]
        if a.winslow > 0:
            Xi, wrep[int(r)] = winslow(Xi, iters=a.winslow)
        if a.dump_block is not None and int(r) in a.dump_block:
            np.save(bd / f"dump_block{int(r)}_win.npy", Xi)
        ni, nj, nk = Xi.shape[:3]
        idx = np.arange(ni * nj * nk).reshape(ni, nj, nk) + base
        for sl, _g, name, _n in sides[r][0]:
            if name in PATCHES:
                G2 = idx[sl]
                bq.append(np.c_[G2[:-1, :-1].ravel(), G2[1:, :-1].ravel(), G2[1:, 1:].ravel(), G2[:-1, 1:].ravel()])
                bpatch += [list(PATCHES).index(name)] * ((G2.shape[0] - 1) * (G2.shape[1] - 1))
        pts.append(Xi.reshape(-1, 3))
        hexes.append(tfi.block_cells(idx))
        bids += [r] * ((ni - 1) * (nj - 1) * (nk - 1))
        base += ni * nj * nk
    Pf, Hf = np.vstack(pts), np.vstack(hexes)
    Pw, remap = tfi.weld(Pf, 1e-9)
    Hw = remap[Hf]
    watertight, _ = tfi.check_watertight(Pw, Hw, verbose=False)
    qs = qual.summary(Pw, Hw)
    cl = qual.classify(Pw, Hw)
    fc = {}
    rr = np.hypot(Pw[:, 0], Pw[:, 1])
    for nm, r0 in (("hub", G.r_hub), ("shroud", G.r_shroud)):
        on = np.abs(rr - r0) < 1e-5
        cnt = on[Hw].sum(1)
        hs = [abs(rr[c][~on[c]].mean() - r0) for c in Hw[cnt == 4]]
        fc[nm] = {"cells": len(hs), "min": float(np.min(hs)), "median": float(np.median(hs)),
                  "max": float(np.max(hs))} if hs else None
    if blade_pts:
        from scipy.spatial import cKDTree
        onb = cKDTree(np.vstack(blade_pts)).query(Pw)[0] < 1e-9
        cnt = onb[Hw].sum(1)
        hs = []
        for c in Hw[cnt == 4]:
            w, o = Pw[c[onb[c]]].mean(0), Pw[c[~onb[c]]].mean(0)
            hs.append(float(np.linalg.norm(o - w)))
        fc["blade"] = {"cells": len(hs), "min": float(np.min(hs)), "median": float(np.median(hs)),
                       "max": float(np.max(hs))} if hs else None
    rep = {"out": str(out), "cells": int(len(Hw)), "points": int(len(Pw)), "watertight": bool(watertight),
           "edges": len(edges), "faces": len(faces), "edge_reprojection_max": erep["max_proj"],
           "face_correction_max": frep, "ring_wall": rrep, "cutout_morph": mrep, "winslow": wrep, "quality": qs, "first_cell_height": fc,
           "classes": [{"class": ci, "axes": len(classes[ci]), **cc[ci]} for ci in cc],
           "gate": {"watertight": bool(watertight), "inverted_0": qs["inverted"] == 0}}
    import meshio
    # boundary quads with their patch (index into PATCHES, names in field_data)
    Qw = remap[np.vstack(bq)]
    nq = len(Qw)
    zq = np.zeros(nq)
    meshio.write(str(out), meshio.Mesh(Pw, [("hexahedron", Hw), ("quad", Qw)],
                                       cell_data={"block_id": [np.array(bids), np.full(nq, -1)],
                                                  "sj": [cl["sj"], zq],
                                                  "inverted": [cl["inverted"].astype(float), zq],
                                                  "concave": [cl["concave"].astype(float), zq],
                                                  "patch": [np.full(len(Hw), -1), np.array(bpatch)]},
                                       field_data={v: np.array([i, 2]) for i, v in enumerate(PATCHES.values())}))
    out.with_suffix(".json").write_text(json.dumps(rep, indent=1))
    # entities for B6 (the ring meshes against the same cutout faces)
    fin = {(k[0] if isinstance(k, tuple) else k): v for k, v in faces.items()
           if not (isinstance(k, tuple) and k[1] == "old")}
    np.savez(bd / "b5_entities.npz",
             face_keys=np.array([np.array(sorted(k)) for k in fin], dtype=object),
             face_grids=np.array([v[0] for v in fin.values()], dtype=object),
             face_points=np.array([v[1] for v in fin.values()], dtype=object),
             face_surface=np.array([str(v[2]) for v in fin.values()]))
    print(json.dumps({k: rep[k] for k in ("cells", "watertight", "edge_reprojection_max", "face_correction_max",
                                          "quality", "first_cell_height", "gate")}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
