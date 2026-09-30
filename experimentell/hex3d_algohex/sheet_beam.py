"""Sheet-level block extraction: beam search (and greedy baseline) over the
SEPARATRIX CUT SET -- no mesh edits.

Premise (measured on base_a, 2026-09-29): the raw base complex over-segments.
Cutting ALL 24 labeled separatrix sheets gives 84 blocks; dropping the right
7 sheets yields the canonical 12-block structure. The drop decision is pure
subset selection on the cut set; the fine hex mesh is never modified, so all
mesh-quality guards of the cell-level collapse are irrelevant on this level.

State  = set of ACTIVE sheets (their faces still separate regions).
Child  = state minus one sheet (regions on both sides merge).
Search = beam over subsets; the partition depends ONLY on the final subset,
         so states are deduplicated by the active-sheet set itself.
Key    = fewer blocks first, then fewer (approximate) non-cuboid regions.

Inputs: perturb runs + smoke (raw AlgoHex .ovm) and dataset samples whose raw
mesh only survives as blocks.vtk (post-greedy fine mesh).

Outputs:
  /work/vtk/12_sheet_beam/<sample>/<sample>_final_regions.vtk
      the fine hex mesh re-coloured by the beam's final regions
  /work/vtk/12_sheet_beam/<sample>/drops.csv          drop sequence
  /work/analysis/sheet_beam_results.csv               summary over samples
"""
import contextlib
import csv
import io
import os
import sys
import time
from collections import defaultdict

import numpy as np

sys.path[:0] = ["/work/dp3d_datagen",
                "/work/dp3d_datagen/experimentell/hex3d_algohex"]
import base_complex as bc                                            # noqa: E402
import clean_blocks as cb                                            # noqa: E402
import ovm_io                                                        # noqa: E402
import meshio
import tet_prep_v5 as v5  # noqa: E402

OUT = "/work/vtk/12_sheet_beam"
CSV_OUT = "/work/analysis/sheet_beam_results.csv"
os.makedirs(OUT, exist_ok=True)


# ------------------------------------------------------------------ inputs
def load_mesh(sample_dir, mode):
    """Fine hex mesh + face->separatrix sheet labeling (bfx.label_sheets).

    This is the recorded-baseline variant of the cut-set analysis. A mesh
    classes variant (cb.mesh_sheets parallel-edge classes) was also tried: it
    can NOT be emulated by cut-set union (collapse welds vertices, which
    merges faces beyond dual adjacency), see the findings note."""
    if mode == "ovm":
        ovm = sample_dir / "hex_hex_x.ovm"
        if not ovm.exists():
            ovm = sample_dir / "hex_hex_smoke.ovm"
        P, e, f_, poly = ovm_io.read_ovm(str(ovm))
        H, _sk = ovm_io.ovm_to_cells(P, e, f_, poly)
    else:
        m = meshio.read(str(sample_dir / "blocks.vtk"))
        key = next(k for k in m.cells_dict
                   if str(k).lower().startswith("hexahedron"))
        H = np.asarray(m.cells_dict[key], np.int64)
        P = np.asarray(m.points, float)
    f2h, e2h = bc.build_topology(H)
    with contextlib.redirect_stdout(io.StringIO()):
        add = cb.fill_cavities(P, H, f2h, e2h)
    if len(add):
        H = np.vstack([H, add])
        f2h, e2h = bc.build_topology(H)
    with contextlib.redirect_stdout(io.StringIO()):
        sing = bc.singular_edges(H, P, f2h, e2h)
        sheet_of, _n = __import__("block_faces").label_sheets(H, f2h, sing)
    sheet_faces = defaultdict(list)
    for fc, s in sheet_of.items():
        sheet_faces[s].append(fc)
    return P, H, f2h, sheet_faces


# ----------------------------------------------------------------- searches
class SheetPartition:
    """Union-find over hexes across faces whose sheet is NOT active."""

    def __init__(self, H, f2h, sheet_of):
        self.n = len(H)
        self.f2h = f2h
        pairs, sheets = [], []
        for fk, hs in f2h.items():
            if len(hs) != 2:
                continue
            a, b = hs
            pairs.append((a, b))
            sheets.append(sheet_of.get(fk, -1))   # -1 = regular face, always joins
        self.pi = np.asarray([p[0] for p in pairs], np.int32)
        self.pj = np.asarray([p[1] for p in pairs], np.int32)
        self.ps = np.asarray(sheets, np.int32)

    def labels(self, active):
        """active: set of sheet ids still separating. Returns (lbl, n)."""
        allowed = ~np.isin(self.ps, list(active))
        parent = np.arange(self.n, dtype=np.int32)

        def find(x):
            root = x
            while parent[root] != root:
                root = parent[root]
            while parent[x] != root:
                parent[x], x = root, parent[x]
            return root

        for a, b in zip(self.pi[allowed], self.pj[allowed]):
            ra, rb = find(int(a)), find(int(b))
            if ra != rb:
                parent[rb] = ra
        roots = np.array([find(int(i)) for i in range(self.n)])
        uniq, lbl = np.unique(roots, return_inverse=True)
        return lbl.astype(np.int64), len(uniq)

    def noncuboid(self, lbl, n_regions):
        """Approximate regularity: distinct neighbour regions per region; a
        hexahedral block region borders exactly 6 neighbours (or fewer for
        boundary-only ones). Returns count of regions != 2..6."""
        nbr = defaultdict(set)
        for a, b in zip(self.pi, self.pj):
            ra, rb = lbl[a], lbl[b]
            if ra != rb:
                nbr[ra].add(rb)
                nbr[rb].add(ra)
        bad = 0
        for r in range(n_regions):
            k = len(nbr[r])
            if k > 6:
                bad += 1
        return bad


def census(P, H, lbl, n_regions, f2h, labeller):
    """Exact block-shape census, same semantics as the pipeline's
    _structure_stats: per region a face record (boundary label or neighbour
    region) -> _patch_split -> cuboid_status / _excess."""
    HF = cb.HF
    rec = [[] for _ in range(n_regions)]
    for fkey, hs in f2h.items():
        if len(hs) == 1:
            lp = None
            for fc in HF:
                c = H[hs[0]]
                if frozenset(int(c[i]) for i in fc) == fkey:
                    lp = [int(c[i]) for i in fc]
                    break
            ids = labeller.label(P[np.array(lp)].mean(0)[None, :])
            rec[lbl[hs[0]]].append((("P", int(ids[0])), lp))
        elif len(hs) == 2:
            a, b = hs
            if lbl[a] == lbl[b]:
                continue
            c = H[a]
            lp = None
            for fc in HF:
                if frozenset(int(c[i]) for i in fc) == fkey:
                    lp = [int(c[i]) for i in fc]
                    break
            rec[lbl[a]].append((("B", int(lbl[b])), lp))
            rec[lbl[b]].append((("B", int(lbl[a])), lp))
    from collections import defaultdict
    # neighbour faces merged to one patch per (kind, id) via _patch_split
    cub = exc = 0
    for r in range(n_regions):
        patches = cb._patch_split(rec[r])
        if cb.cuboid_status(patches)[0] == "cuboid":
            cub += 1
        exc += cb._excess(patches)
    return cub, exc


def search(P, H, f2h, sheet_faces, labeller, width, greedy=False,
           stop_patience=2, census_children=8, label="beam"):
    part = SheetPartition(H, f2h, _sheet_of_of(f2h, sheet_faces))
    sheet_ids = sorted(sheet_faces)
    allsheets = frozenset(sheet_ids)
    lbl0, n0 = part.labels(allsheets)
    cub0, exc0 = census(P, H, lbl0, n0, f2h, labeller)
    print(f"  [{label}] start: {n0} blocks, cuboids {cub0}, excess {exc0}",
          flush=True)

    def metrics(active):
        lbl, n_blocks = part.labels(active)
        cub, exc = census(P, H, lbl, n_blocks, f2h, labeller)
        return lbl, n_blocks, cub, exc

    best_lbl, best_blocks, best_cub, best_exc = metrics(allsheets)
    best_key = (n0 - cub0, exc0, n0)   # noncuboid count, excess, blocks
    best_state = dict(active=allsheets, blocks=n0, dropped=())
    beam = [best_state]
    drops_log = []
    patience, d = 0, 0
    t0 = time.time()
    while beam and d < len(sheet_ids) and patience < stop_patience:
        d += 1
        children = {}
        for st in beam:
            for s in st["active"]:
                act = st["active"] - {s}
                if act in children:
                    continue
                lbl, n_blocks = part.labels(act)
                children[act] = dict(active=act, blocks=n_blocks, lbl=lbl,
                                     dropped=st["dropped"] + (s,))
        if not children:
            break
        short = sorted(children.values(),
                       key=lambda c: c["blocks"])[:max(width, census_children)]
        for c in short:
            c["cub"], c["exc"] = census(P, H, c["lbl"], c["blocks"], f2h,
                                        labeller)
            c["key"] = (c["blocks"] - c["cub"], c["exc"], c["blocks"])
        cand = sorted(short, key=lambda c: c["key"])[:width]
        best_of_depth = cand[0]
        drops_log.append((d, best_of_depth["dropped"][-1],
                          best_state["blocks"], best_of_depth["blocks"]))
        if best_of_depth["key"] < best_key:
            best_state = best_of_depth
            best_key = best_of_depth["key"]
            best_lbl = best_of_depth["lbl"]
            patience = 0
        else:
            patience += 1
        if greedy:
            beam = [best_of_depth] if best_of_depth["blocks"] < \
                beam[0]["blocks"] else []
        else:
            beam = cand
    lbl, n_blocks = part.labels(best_state["active"])
    dropped = sorted(allsheets - best_state["active"])
    cub, exc = census(P, H, lbl, n_blocks, f2h, labeller)
    return dict(lbl=lbl, n_blocks=n_blocks, dropped=dropped,
                drops_log=drops_log, active=best_state["active"],
                seconds=time.time() - t0, bad=n_blocks - cub, exc=exc)


def _sheet_of_of(f2h, sheet_faces):
    """reconstruct face->sheet map from sheet_faces (face keys -> list)."""
    back = {}
    for s, faces in sheet_faces.items():
        for fc in faces:
            back[fc] = s
    return back


# ------------------------------------------------------------------ export
def write_regions_vtk(path, P, H, lbl):
    meshio.write(path, meshio.Mesh(np.asarray(P, float),
                                   [("hexahedron", np.asarray(H, np.int64))],
                                   cell_data={"block_id": [lbl]}),
                 file_format="vtk", binary=False)


def canonical_hash(npz_path):
    """Combination-free WL hash on the npz block structure: nodes = blocks +
    corner vertices, edges = incidence. Same convention as topology_stats'
    topo_comb minus corner surface labels."""
    import networkx as nx
    s = np.load(npz_path, allow_pickle=True)
    V, Bk = s["vertices"], s["blocks"]
    G = nx.Graph()
    for i in range(len(V)):
        G.add_node(("v", i), lab="v")
    for bi, blk in enumerate(Bk):
        G.add_node(("B", bi), lab="B")
        for v in blk:
            G.add_edge(("B", bi), ("v", int(v)))
    return nx.weisfeiler_lehman_graph_hash(G, node_attr="lab", iterations=6)


def comb_hash_regions(lbl, n_regions, part):
    """Same WL shape as canonical_hash, but nodes = regions + VERTEX
    positions? Regions don't carry corner vertices directly; approximate
    with the dual graph (region-region) -- documented, not mixed with the
    npz hash. Comparisons happen within this function only."""
    raise SystemExit("unused")


def make_labeller(sd, mode):
    """SurfaceLabeller from input_vtk (perturb tet.vtk) or npz surface."""
    tet = sd / "tet.vtk"
    if tet.exists():
        P0, tri, tid = cb.read_input_surface(str(tet))
        return cb.SurfaceLabeller(P0, tri, tid, v5.NAMES)
    npz = sd / "sample.npz"
    if npz.exists():
        s = np.load(str(npz), allow_pickle=True)
        return cb.SurfaceLabeller(s["surface_points"], s["surface_tris"],
                                  s["surface_tri_label"], v5.NAMES)
    raise FileNotFoundError(f"no label source in {sd}")


def main():
    samples = []
    for d in sorted(os.listdir("/work/perturb")):
        sd = __import__("pathlib").Path("/work/perturb") / d
        if (sd / "hex_hex_x.ovm").exists():
            samples.append((d, sd, "ovm"))
    smoke = __import__("pathlib").Path("/work/pipeline_smoke/sample")
    if (smoke / "hex_hex_smoke.ovm").exists():
        samples.append(("smoke", smoke, "ovm"))
    for geo in ("machine_0004_n2000", "machine_0034_n2000"):
        sd = __import__("pathlib").Path(
            f"/opt/stack/meshtron/data/hex3d_algohex/batch/{geo}")
        if (sd / "blocks.vtk").exists():
            samples.append((geo, sd, "blocks"))

    rows = []
    for name, sd, mode in samples:
        print(f"=== {name} ({mode}) ===", flush=True)
        try:
            P, H, f2h, sheet_faces = load_mesh(sd, mode)
        except Exception as e:
            print(f"  load FAILED: {e}", flush=True)
            continue
        part = SheetPartition(H, f2h, _sheet_of_of(f2h, sheet_faces))
        sheet_ids = sorted(sheet_faces)
        allsheets = frozenset(sheet_ids)
        raw_lbl, raw_blocks = part.labels(allsheets)

        outdir = __import__("pathlib").Path(OUT) / name
        outdir.mkdir(parents=True, exist_ok=True)

        for method in ("greedy", "beam"):
            labeller = make_labeller(sd, mode)
            res = search(P, H, f2h, sheet_faces, labeller, 16,
                         greedy=(method == "greedy"), label=method)
            tag = "beam" if method != "greedy" else "greedy"
            write_regions_vtk(outdir / f"{name}_{tag}_regions.vtk",
                              P, H, res["lbl"])
            with open(outdir / f"{tag}_drops.csv", "w", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(["depth", "sheet_dropped", "blocks_before",
                            "blocks_after"])
                w.writerows(res["drops_log"])
            row = dict(sample=name, mode=mode, n_hexes=len(H),
                       n_sheets=len(sheet_faces), raw_blocks=raw_blocks,
                       method=tag, final_blocks=res["n_blocks"],
                       noncuboid_apx=res["bad"],
                       sheets_dropped=len(res["dropped"]),
                       dropped_ids=";".join(map(str, res["dropped"])),
                       seconds=round(res["seconds"], 1))
            rows.append(row)
            print(f"  {tag:6s}: {raw_blocks} -> {res['n_blocks']} blocks "
                  f"(cuboids {len(H) and len(res['lbl']) and res['n_blocks'] - res['bad']}), "
                  f"excess {res['exc']}, dropped {len(res['dropped'])} "
                  f"sheets in {res['seconds']:.1f}s", flush=True)

    with open(CSV_OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("written", CSV_OUT, flush=True)


if __name__ == "__main__":
    main()
