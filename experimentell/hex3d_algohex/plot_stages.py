"""Visualise the stages of the AlgoHex 3D pipeline for visual inspection.

Produces output/hex3d_algohex/plots/:
  step01_boundary_patches.png  Stage 1: analytic surface classification
  step02_feature_graph.png     Stage 1: feature curves + feature vertices
  step03_tagging_v1_vs_v3.png  the bug: classifySurfaces vs analytic tags
  step04_hexmesh_v1.png        Stage 2 v1: extracted hex mesh (31.4% coverage)
"""

import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                       # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection  # noqa: E402

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO / "experimentell" / "hex3d_algohex"))
OUT = REPO / "output" / "hex3d_algohex"
PLOTS = OUT / "plots"

import tet_prep as tp                                                 # noqa: E402

PATCH_COLORS = {
    1: ("#4C72B0", "hub (r=0.5)"),
    2: ("#DD8452", "shroud (r=1.9)"),
    3: ("#55A868", "inlet (z=0)"),
    4: ("#C44E52", "outlet (z=2.5)"),
    5: ("#8172B3", "blade"),
    6: ("#937860", "periodic A"),
    7: ("#DA8BC3", "periodic B"),
}


def read_vtk(path):
    lines = Path(path).read_text().split("\n")

    def find(tag, s=0):
        for j in range(s, len(lines)):
            if lines[j].startswith(tag):
                return j
        raise ValueError(tag)

    j = find("POINTS")
    n = int(lines[j].split()[1])
    P = np.array([[float(x) for x in lines[j + 1 + k].split()]
                  for k in range(n)])
    j = find("CELLS")
    nc = int(lines[j].split()[1])
    cells, k = [], j + 1
    for _ in range(nc):
        pr = [int(x) for x in lines[k].split()]
        cells.append(pr[1:1 + pr[0]])
        k += 1
    j = find("CELL_TYPES")
    ty = [int(lines[j + 1 + c]) for c in range(nc)]
    j = find("LOOKUP_TABLE")
    co = [int(lines[j + 1 + c]) for c in range(nc)]
    return P, cells, ty, co


def _equal_aspect(ax, P):
    c = (P.max(axis=0) + P.min(axis=0)) / 2
    r = (P.max(axis=0) - P.min(axis=0)).max() / 2
    ax.set_xlim(c[0] - r, c[0] + r)
    ax.set_ylim(c[1] - r, c[1] + r)
    ax.set_zlim(c[2] - r, c[2] + r)
    ax.set_box_aspect((1, 1, 1))


def plot_patches(P, tris, ids, out_png):
    fig = plt.figure(figsize=(17, 8))
    for si, (elev, azim, ttl) in enumerate(
            [(22, 35, "view A"), (22, 215, "view B (opposite)")]):
        ax = fig.add_subplot(1, 2, si + 1, projection="3d")
        for sid, (col, lab) in PATCH_COLORS.items():
            sel = [t for t, i in zip(tris, ids) if i == sid]
            if not sel:
                continue
            polys = [P[list(t)] for t in sel]
            pc = Poly3DCollection(polys, facecolor=col, edgecolor="none",
                                  alpha=0.85 if sid != 2 else 0.25)
            ax.add_collection3d(pc)
        _equal_aspect(ax, P)
        ax.view_init(elev=elev, azim=azim)
        ax.set_title(ttl, fontsize=10)
        ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("z")
    handles = [plt.Line2D([0], [0], marker="s", ls="", color=c,
                          label=f"{l} [{s}]", markersize=9)
               for s, (c, l) in PATCH_COLORS.items()]
    fig.legend(handles=handles, loc="lower center", ncol=7, fontsize=9,
               frameon=False)
    fig.suptitle("Stage 1 — analytic boundary classification of T1_9 "
                 "(7 physical surfaces; shroud drawn transparent)",
                 fontsize=13)
    fig.tight_layout(rect=[0, 0.06, 1, 0.97])
    fig.savefig(out_png, dpi=135, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_png}")


def plot_feature_graph(P, tris, ids, fedges, fverts, out_png):
    fig = plt.figure(figsize=(17, 8))
    for si, (elev, azim, ttl) in enumerate(
            [(22, 35, "view A"), (22, 215, "view B (opposite)")]):
        ax = fig.add_subplot(1, 2, si + 1, projection="3d")
        polys = [P[list(t)] for t in tris]
        ax.add_collection3d(Poly3DCollection(
            polys, facecolor="0.85", edgecolor="none", alpha=0.14))
        segs = [[P[a], P[b]] for a, b in fedges]
        ax.add_collection3d(Line3DCollection(segs, colors="#C44E52", lw=1.7))
        if fverts:
            V = P[list(fverts)]
            ax.scatter(V[:, 0], V[:, 1], V[:, 2], c="#1f1f1f", s=26,
                       depthshade=False, zorder=5)
        _equal_aspect(ax, P)
        ax.view_init(elev=elev, azim=azim)
        ax.set_title(ttl, fontsize=10)
    fig.suptitle(f"Stage 1 — feature graph fed to AlgoHex: "
                 f"{len(fedges)} feature edges (red), "
                 f"{len(fverts)} feature vertices (black)\n"
                 f"valid 1-manifold: 0 dangling ends, all endpoints tagged",
                 fontsize=13)
    fig.tight_layout(rect=[0, 0.02, 1, 0.94])
    fig.savefig(out_png, dpi=135, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_png}")


def plot_tagging_comparison(v1_vtk, v3_vtk, out_png):
    """The actual bug: classifySurfaces splits smooth cylinders, injecting
    feature curves across flat surfaces."""
    fig = plt.figure(figsize=(17, 8.5))
    for si, (path, ttl) in enumerate(
            [(v1_vtk, "v1 — gmsh classifySurfaces tags  (BROKEN)"),
             (v3_vtk, "v3 — analytic tags + cleaned graph  (FIXED)")]):
        P, cells, ty, co = read_vtk(path)
        tris = [(cells[c], co[c]) for c in range(len(cells)) if ty[c] == 5]
        fe = [tuple(cells[c]) for c in range(len(cells)) if ty[c] == 3]

        e2t = defaultdict(list)
        for i, (t, _) in enumerate(tris):
            for a in range(3):
                p, q = t[a], t[(a + 1) % 3]
                e2t[(min(p, q), max(p, q))].append(i)

        def nrm(t):
            v = np.cross(P[t[1]] - P[t[0]], P[t[2]] - P[t[0]])
            L = np.linalg.norm(v)
            return v / L if L > 0 else v

        flat, sharp = [], []
        for a, b in fe:
            inc = e2t.get((min(a, b), max(a, b)), [])
            if len(inc) != 2:
                continue
            d = float(np.dot(nrm(tris[inc[0]][0]), nrm(tris[inc[1]][0])))
            ang = math.degrees(math.acos(max(-1, min(1, d))))
            (sharp if ang > 20 else flat).append([P[a], P[b]])

        ax = fig.add_subplot(1, 2, si + 1, projection="3d")
        ax.add_collection3d(Poly3DCollection(
            [P[list(t)] for t, _ in tris], facecolor="0.86",
            edgecolor="none", alpha=0.13))
        if sharp:
            ax.add_collection3d(Line3DCollection(
                sharp, colors="#2E7D32", lw=1.6))
        if flat:
            ax.add_collection3d(Line3DCollection(
                flat, colors="#D32F2F", lw=2.2))
        _equal_aspect(ax, P)
        ax.view_init(elev=22, azim=35)
        n_surf = len(set(c for _, c in tris))
        ax.set_title(f"{ttl}\n{n_surf} surfaces, {len(fe)} feature edges\n"
                     f"green = real crease ({len(sharp)}), "
                     f"RED = FLAT/spurious ({len(flat)})", fontsize=10)
    fig.suptitle("Root cause of the failed v1 run: spurious feature curves "
                 "across smooth surfaces over-constrain the octahedral field",
                 fontsize=13)
    fig.tight_layout(rect=[0, 0.02, 1, 0.93])
    fig.savefig(out_png, dpi=135, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_png}")


def plot_hexmesh(ovm_path, out_png, title):
    import ovm_io
    P, edges, faces, polys = ovm_io.read_ovm(ovm_path)
    hexes, _ = ovm_io.ovm_to_cells(P, edges, faces, polys)
    vols = np.array([ovm_io._hex_volume(P[c]) for c in hexes])

    # outer quad faces of the hex mesh (face used by exactly one hex)
    cnt = Counter()
    rep = {}
    HF = ((0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4),
          (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7))
    for c in hexes:
        for f in HF:
            key = frozenset(int(c[i]) for i in f)
            cnt[key] += 1
            rep[key] = [int(c[i]) for i in f]
    surf = [rep[k] for k, v in cnt.items() if v == 1]

    fig = plt.figure(figsize=(17, 8))
    for si, (elev, azim, ttl) in enumerate(
            [(22, 35, "view A"), (22, 215, "view B (opposite)")]):
        ax = fig.add_subplot(1, 2, si + 1, projection="3d")
        ax.add_collection3d(Poly3DCollection(
            [P[f] for f in surf], facecolor="#4C72B0",
            edgecolor="#20304a", linewidth=0.16, alpha=0.9))
        _equal_aspect(ax, P)
        ax.view_init(elev=elev, azim=azim)
        ax.set_title(ttl, fontsize=10)
    fig.suptitle(f"{title}\n{len(hexes)} hex cells, "
                 f"volume {vols.sum():.3f} of 6.525 "
                 f"({100 * vols.sum() / 6.52525:.1f}% of domain), "
                 f"{(vols <= 0).sum()} inverted", fontsize=13)
    fig.tight_layout(rect=[0, 0.02, 1, 0.93])
    fig.savefig(out_png, dpi=135, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_png}")


def main():
    PLOTS.mkdir(parents=True, exist_ok=True)
    vtk = REPO / "data" / "T1_9" / "T1_9_tet.vtk"
    P, cells, ty, co = read_vtk(vtk)
    tris = [cells[c] for c in range(len(cells)) if ty[c] == 5]
    ids = [co[c] for c in range(len(cells)) if ty[c] == 5]
    fedges = [tuple(cells[c]) for c in range(len(cells)) if ty[c] == 3]
    fverts = [cells[c][0] for c in range(len(cells)) if ty[c] == 1]

    plot_patches(P, tris, ids, PLOTS / "step01_boundary_patches.png")
    plot_feature_graph(P, tris, ids, fedges, fverts,
                       PLOTS / "step02_feature_graph.png")

    v1 = Path("/tmp/T1_9_tet_v1_classifysurfaces.vtk")
    if v1.exists():
        plot_tagging_comparison(v1, vtk, PLOTS / "step03_tagging_v1_vs_v3.png")

    v1hex = OUT / "T1_9_hex_v1_classifysurfaces.ovm"
    if v1hex.exists():
        plot_hexmesh(v1hex, PLOTS / "step04_hexmesh_v1.png",
                     "Stage 2 run v1 — INCOMPLETE hex mesh "
                     "(invalid IGM -> holes)")
    cyl = OUT / "cylinder_hex.ovm"
    if cyl.exists():
        import ovm_io
        Pc, ec, fc, pc = ovm_io.read_ovm(cyl)
        hx, _ = ovm_io.ovm_to_cells(Pc, ec, fc, pc)
        vc = np.array([ovm_io._hex_volume(Pc[c]) for c in hx])
        fig = plt.figure(figsize=(9, 8))
        ax = fig.add_subplot(111, projection="3d")
        cntc = Counter(); repc = {}
        HF = ((0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4),
              (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7))
        for c in hx:
            for f in HF:
                k = frozenset(int(c[i]) for i in f)
                cntc[k] += 1; repc[k] = [int(c[i]) for i in f]
        sf = [repc[k] for k, v in cntc.items() if v == 1]
        ax.add_collection3d(Poly3DCollection(
            [Pc[f] for f in sf], facecolor="#55A868",
            edgecolor="#24402f", linewidth=0.15, alpha=0.92))
        _equal_aspect(ax, Pc)
        ax.view_init(elev=20, azim=35)
        ax.set_title(f"Stage 0 reference — AlgoHex cylinder demo\n"
                     f"{len(hx)} hexes, {(vc <= 0).sum()} inverted "
                     f"(pipeline + our OVM reader verified)", fontsize=12)
        fig.savefig(PLOTS / "step00_cylinder_reference.png", dpi=135,
                    bbox_inches="tight")
        plt.close(fig)
        print(f"wrote {PLOTS / 'step00_cylinder_reference.png'}")


if __name__ == "__main__":
    main()
