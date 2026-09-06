"""Step-by-step visual analysis of the 3D hex pipeline.

The 3D counterpart of the 2D showcase in `dp3d/plotting.py`
(`output/plots/step01…step11`). Eleven steps, each written three ways:

* **`.vtk` — the authoritative artifact.** A 3D block structure is not
  legible in a flat render; the VTK is what you open in ParaView to actually
  check something. Everything else is a summary of it.
* **`.png`** — a fixed-camera render (pyvista, off-screen) so the series can
  be skimmed and dropped into a document.
* **`.html`** — the same scene as an interactive plotly figure, self-contained
  and embeddable in a presentation.

All three come from one description per step, so they cannot drift apart.

Conventions follow the 2D series: one fixed viewpoint for every step, a faint
background surface for context, and a `_labeled` variant where per-entity
indices are what make the figure readable.
"""

import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "experimentell" / "hex3d_algohex"))

import base_complex as bc                                             # noqa: E402
import clean_blocks as cb                                             # noqa: E402
import ovm_io                                                         # noqa: E402

OUT = REPO / "output" / "hex3d_algohex" / "showcase"
DATA = REPO / "data" / "T1_9"
HEX = REPO / "output" / "hex3d_algohex"

# one camera for the whole series, so steps are comparable at a glance
CAMERA = [(6.0, -5.0, 6.5), (1.1, 0.0, 1.25), (0.0, 0.0, 1.0)]
WINDOW = (1400, 1050)
BG = "white"

# a colour-blind-safe qualitative ramp, used wherever the scalar is an id
CATEGORICAL = ["#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE",
               "#AA3377", "#BBBBBB", "#000000"]


# --------------------------------------------------------------------------
# readers
# --------------------------------------------------------------------------

def read_algohex_vtk(path):
    """Full read of an AlgoHex input VTK: points plus every cell type with
    its int `color`. `clean_blocks.read_input_surface` keeps only triangles;
    the feature graph lives in the line and vertex cells."""
    txt = Path(path).read_text().split("\n")

    def seek(prefix, i):
        while not txt[i].startswith(prefix):
            i += 1
        return i

    i = seek("POINTS", 0)
    n = int(txt[i].split()[1])
    P = np.array([[float(x) for x in txt[i + 1 + k].split()] for k in range(n)])
    i = seek("CELLS", i + n)
    nc = int(txt[i].split()[1])
    conn = [[int(x) for x in txt[i + 1 + k].split()[1:]] for k in range(nc)]
    i = seek("CELL_TYPES", i + nc)
    types = np.array([int(txt[i + 1 + k]) for k in range(nc)])
    i = seek("LOOKUP_TABLE", i + nc)
    col = np.array([int(txt[i + 1 + k]) for k in range(nc)])
    out = {}
    for name, vt in (("vertex", 1), ("line", 3), ("triangle", 5), ("tetra", 10)):
        sel = np.where(types == vt)[0]
        out[name] = (np.array([conn[k] for k in sel], np.int64) if len(sel)
                     else np.zeros((0, 1), np.int64), col[sel])
    return P, out


def hexmesh(ovm_path):
    P, e, f, poly = ovm_io.read_ovm(ovm_path)
    H, _ = ovm_io.ovm_to_cells(P, e, f, poly)
    return P, H


# --------------------------------------------------------------------------
# pyvista / plotly plumbing
# --------------------------------------------------------------------------

def _pv():
    import pyvista as pv
    pv.OFF_SCREEN = True
    return pv


def hex_grid(P, H, **cell_arrays):
    pv = _pv()
    cells = np.hstack([np.full((len(H), 1), 8, np.int64), H]).ravel()
    g = pv.UnstructuredGrid(cells, np.full(len(H), 12, np.uint8), np.asarray(P))
    for k, v in cell_arrays.items():
        g.cell_data[k] = np.asarray(v)
    return g


def tri_surface(P, tri, **cell_arrays):
    pv = _pv()
    faces = np.hstack([np.full((len(tri), 1), 3, np.int64), tri]).ravel()
    s = pv.PolyData(np.asarray(P), faces)
    for k, v in cell_arrays.items():
        s.cell_data[k] = np.asarray(v)
    return s


def line_set(P, segs, **cell_arrays):
    pv = _pv()
    lines = np.hstack([np.full((len(segs), 1), 2, np.int64), segs]).ravel()
    # pass lines to the constructor: pv.PolyData(points) alone makes a vertex
    # cell per point, so n_cells would be points + lines and every cell array
    # would be rejected for the wrong length
    ln = pv.PolyData(np.asarray(P), lines=lines)
    for k, v in cell_arrays.items():
        ln.cell_data[k] = np.asarray(v)
    return ln


def _render(actors, png, title):
    pv = _pv()
    p = pv.Plotter(off_screen=True, window_size=WINDOW)
    p.set_background(BG)
    PLOTLY_ONLY = {"mesh", "kind", "colorscale"}
    for a in actors:
        kw = {k: v for k, v in a.items() if k not in PLOTLY_ONLY}
        if a.get("kind") == "points":
            kw.setdefault("style", "points")
            kw.setdefault("render_points_as_spheres", True)
        p.add_mesh(a["mesh"], **kw)
    # fix the view DIRECTION for the whole series, but let the camera fit the
    # actual bounds: a constant zoom cropped step 10, whose small multiple is
    # three times as wide as every other step
    p.camera_position = CAMERA
    p.reset_camera()
    p.camera.zoom(1.15)
    p.add_text(title, font_size=11, color="black")
    Path(png).parent.mkdir(parents=True, exist_ok=True)
    p.screenshot(str(png))
    p.close()


def _plotly(actors, html, title):
    import plotly.graph_objects as go
    traces = []
    for a in actors:
        m = a["mesh"]
        name = a.get("label", "")
        opacity = a.get("opacity", 1.0)
        scal = a.get("scalars")
        if a.get("kind") == "lines":
            pts = m.points
            L = m.lines.reshape(-1, 3)[:, 1:]
            X, Y, Z = [], [], []
            for u, v in L:
                X += [pts[u][0], pts[v][0], None]
                Y += [pts[u][1], pts[v][1], None]
                Z += [pts[u][2], pts[v][2], None]
            traces.append(go.Scatter3d(x=X, y=Y, z=Z, mode="lines",
                                       line=dict(width=a.get("line_width", 4),
                                                 color=a.get("color", "#EE6677")),
                                       name=name, hoverinfo="name"))
            continue
        if a.get("kind") == "points":
            pts = m.points
            traces.append(go.Scatter3d(x=pts[:, 0], y=pts[:, 1], z=pts[:, 2],
                                       mode="markers",
                                       marker=dict(size=a.get("point_size", 4),
                                                   color=a.get("color", "#AA3377")),
                                       name=name, hoverinfo="name"))
            continue
        surf = m.extract_surface().triangulate() if m.n_cells else None
        if surf is None or surf.n_points == 0:
            continue
        f = surf.faces.reshape(-1, 4)[:, 1:]
        v = surf.points
        kw = dict(x=v[:, 0], y=v[:, 1], z=v[:, 2],
                  i=f[:, 0], j=f[:, 1], k=f[:, 2],
                  opacity=opacity, name=name, showlegend=bool(name),
                  flatshading=True)
        if scal and scal in surf.cell_data:
            val = np.asarray(surf.cell_data[scal], float)
            kw.update(intensity=val, intensitymode="cell",
                      colorscale=a.get("colorscale", "Viridis"),
                      colorbar=dict(title=scal, len=0.6))
        else:
            kw.update(color=a.get("color", "#BBBBBB"))
        traces.append(go.Mesh3d(**kw))
    fig = go.Figure(traces)
    fig.update_layout(title=title, scene=dict(aspectmode="data",
                                              xaxis_title="x", yaxis_title="y",
                                              zaxis_title="z"),
                      margin=dict(l=0, r=0, t=40, b=0), template="plotly_white")
    Path(html).parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(html), include_plotlyjs="cdn", full_html=True)


def emit(step, name, actors, title, vtk_objects):
    """Write one step three ways. `vtk_objects` is {suffix: mesh}."""
    OUT.mkdir(parents=True, exist_ok=True)
    stem = f"step{step:02d}_{name}"
    for sfx, m in vtk_objects.items():
        path = OUT / (f"{stem}.vtk" if not sfx else f"{stem}_{sfx}.vtk")
        m.save(str(path))
    _render(actors, OUT / f"{stem}.png", title)
    _plotly(actors, OUT / f"{stem}.html", title)
    print(f"[showcase] {stem}: {len(vtk_objects)} vtk + png + html")


# --------------------------------------------------------------------------
# the eleven steps
# --------------------------------------------------------------------------

SURF_COLORS = {1: "#4477AA", 2: "#EE6677", 3: "#228833", 4: "#CCBB44",
               5: "#66CCEE", 6: "#AA3377", 7: "#BBBBBB", 8: "#EE99AA",
               9: "#999933"}


def step01_input_surfaces(tag="v11"):
    """The labelled boundary AlgoHex is given. One colour per surface."""
    import tet_prep_v5 as v5
    P, cells = read_algohex_vtk(DATA / "T1_9_tet_v5.vtk")
    tri, ids = cells["triangle"]
    s = tri_surface(P, tri, surface_id=ids)
    counts = Counter(ids.tolist())
    lbl = ", ".join(f"{v5.NAMES.get(k, k)} {n}" for k, n in sorted(counts.items()))
    print(f"[step01] {len(tri)} boundary triangles: {lbl}")
    acts = [dict(mesh=s, scalars="surface_id", cmap="tab10",
                 show_edges=False, label="surfaces", colorscale="Portland")]
    emit(1, "input_surfaces", acts,
         f"01  AlgoHex input: {len(set(ids.tolist()))} labelled surfaces, "
         f"{len(tri)} triangles", {"": s})
    return counts


def step02_feature_graph():
    """Feature curves and feature vertices on that same boundary."""
    P, cells = read_algohex_vtk(DATA / "T1_9_tet_v5.vtk")
    tri, ids = cells["triangle"]
    ln, _ = cells["line"]
    vt, _ = cells["vertex"]
    deg = Counter()
    for a, b in ln:
        deg[a] += 1
        deg[b] += 1
    hist = dict(sorted(Counter(deg.values()).items()))
    print(f"[step02] {len(ln)} feature edges, valence {hist}, "
          f"{len(vt)} feature vertices")
    surf = tri_surface(P, tri, surface_id=ids)
    curves = line_set(P, ln, feature=np.ones(len(ln), int))
    pts = _pv().PolyData(P[vt.ravel()]) if len(vt) else None
    acts = [dict(mesh=surf, color="#DDDDDD", opacity=0.25, label="boundary"),
            dict(mesh=curves, kind="lines", color="#EE6677", line_width=5,
                 label="feature curves")]
    vtks = {"": curves, "surface": surf}
    if pts is not None:
        acts.append(dict(mesh=pts, kind="points", color="#000000",
                         point_size=12, label="feature vertices"))
        vtks["vertices"] = pts
    emit(2, "feature_graph", acts,
         f"02  feature graph: {len(ln)} edges, valence {hist}, "
         f"{len(vt)} vertices", vtks)


def step03_domain_variants():
    """The three domains that were actually run, side by side."""
    variants = [("full (v5)", "T1_9_tet_v2.vtk"),
                ("minus blade O-grid (v10)", "T1_9_tet_v6.vtk"),
                ("minus O-grid + BL (v11)", "T1_9_tet_v5.vtk")]
    acts, vtks = [], {}
    for i, (name, fn) in enumerate(variants):
        f = DATA / fn
        if not f.exists():
            print(f"[step03] {fn} missing, skipped")
            continue
        P, cells = read_algohex_vtk(f)
        tri, ids = cells["triangle"]
        ntet = len(cells["tetra"][0])
        nfeat = len(cells["line"][0])
        print(f"[step03] {name:26s} {ntet:6d} tets, {len(tri):6d} boundary "
              f"tris, {nfeat:4d} feature edges")
        s = tri_surface(P, tri, surface_id=ids)
        vtks[f"v{i}"] = s                  # untranslated: the VTK is data
        acts.append(dict(mesh=s.translate((i * 4.2, 0, 0), inplace=False),
                         scalars="surface_id", colorscale="Portland",
                         label=name))
    emit(3, "domain_variants", acts,
         "03  domain variants: full / minus O-grid / minus O-grid + BL", vtks)


def _singular_graph(P, H, topo=None):
    f2h, e2h = topo if topo is not None else bc.build_topology(H)
    bnd = {fk for fk, hs in f2h.items() if len(hs) == 1}
    bnde = set()
    for fk in bnd:
        for g in cb._face_edges(cb._loop_of(H, f2h[fk][0], fk)):
            bnde.add(g)
    val = {g: len(set(hs)) for g, hs in e2h.items()}
    sing = [g for g in val if g not in bnde and val[g] != 4]
    adj = defaultdict(list)
    for a, b in sing:
        adj[a].append(b)
        adj[b].append(a)
    nodes = {v for v, n in adj.items() if len(n) != 2}
    seen, arcs = set(), []

    def walk(a, b):
        p = [a, b]
        seen.add((min(a, b), max(a, b)))
        while p[-1] not in nodes:
            nx = [w for w in adj[p[-1]] if w != p[-2]]
            if len(nx) != 1:
                break
            g = (min(p[-1], nx[0]), max(p[-1], nx[0]))
            if g in seen:
                break
            seen.add(g)
            p.append(nx[0])
        return p

    for v in nodes:
        for w in adj[v]:
            if (min(v, w), max(v, w)) not in seen:
                arcs.append(walk(v, w))
    for a, b in sing:
        if (a, b) not in seen:
            arcs.append(walk(a, b))
    return sing, val, arcs, bnd


def step04_singular_graph(tag="v11"):
    """The frame-field singularity graph in 3D, coloured by valence.

    Nothing in this branch exported this before -- it existed only as a
    printed table -- and it is the object `FRAMEFIELD_PLAN.md` is about."""
    import tet_prep as tp
    P, H = hexmesh(HEX / f"T1_9_hex_{tag}.ovm")
    topo = bc.build_topology(H)
    sing, val, arcs, bnd = _singular_graph(P, H, topo)
    segs, sval, sarc = [], [], []
    for ai, p in enumerate(arcs):
        for i in range(len(p) - 1):
            segs.append((p[i], p[i + 1]))
            sval.append(val[(min(p[i], p[i + 1]), max(p[i], p[i + 1]))])
            sarc.append(ai)
    segs = np.array(segs, np.int64)
    g = line_set(P, segs, valence=np.array(sval), arc_id=np.array(sarc))
    at_out = sum(1 for p in arcs for q in (p[0], p[-1])
                 if abs(P[q][2] - tp.Z_OUTLET) < 1e-3)
    at_in = sum(1 for p in arcs for q in (p[0], p[-1])
                if abs(P[q][2] - tp.Z_INLET) < 1e-3)
    print(f"[step04] {tag}: {len(sing)} singular edges "
          f"{dict(sorted(Counter(val[g_] for g_ in sing).items()))}, "
          f"{len(arcs)} arcs, {at_out} endpoints on the outlet, {at_in} on the inlet")
    # build_topology ONCE -- it was inside the comprehension and therefore ran
    # per boundary face, about 10k times over a 61k-cell mesh
    loops = [cb._loop_of(H, topo[0][fk][0], fk) for fk in bnd]
    hull = tri_surface(P, np.array([[l[0], l[1], l[2]] for l in loops]
                                   + [[l[0], l[2], l[3]] for l in loops]))
    ends = _pv().PolyData(np.array([P[q] for p in arcs for q in (p[0], p[-1])]))
    acts = [dict(mesh=hull, color="#DDDDDD", opacity=0.15, label="hex boundary"),
            dict(mesh=g, kind="lines", color="#EE6677", line_width=7,
                 label="singular arcs"),
            dict(mesh=ends, kind="points", color="#000000", point_size=12,
                 label="arc endpoints")]
    emit(4, f"singular_graph_{tag}", acts,
         f"04  singularity graph {tag}: {len(arcs)} arcs, {len(sing)} edges, "
         f"{at_out} endpoints on the outlet", {"": g, "endpoints": ends})
    return arcs


def step05_hexmesh(tag="v11"):
    P, H = hexmesh(HEX / f"T1_9_hex_{tag}.ovm")
    sj = cb.scaled_jacobians(P, H)
    vol = sum(ovm_io._hex_volume(P[c]) for c in H)
    print(f"[step05] {tag}: {len(H)} cells, volume {vol:.4f}")
    g = hex_grid(P, H, scaled_jacobian=sj)
    clip = g.clip("y", origin=(0, 0.05, 0))
    acts = [dict(mesh=g.extract_surface(), color="#BBBBBB", opacity=0.20,
                 label="boundary"),
            dict(mesh=clip, color="#4477AA", show_edges=True, label="clip")]
    emit(5, f"hexmesh_{tag}", acts,
         f"05  extracted hex mesh {tag}: {len(H)} cells, volume {vol:.4f}",
         {"": g, "clip": clip})


def step06_quality(tag="v11"):
    P, H = hexmesh(HEX / f"T1_9_hex_{tag}.ovm")
    sj = cb.scaled_jacobians(P, H)
    bad = np.where(sj <= 0)[0]
    print(f"[step06] {tag}: min {sj.min():.4f} mean {sj.mean():.4f}, "
          f"{len(bad)} inverted ({100 * len(bad) / len(H):.3f} %)")
    g = hex_grid(P, H, scaled_jacobian=sj)
    acts = [dict(mesh=g.extract_surface(), scalars="scaled_jacobian",
                 cmap="RdYlBu", clim=(0, 1), colorscale="RdYlBu",
                 label="scaled Jacobian")]
    vtks = {"": g}
    if len(bad):
        inv = hex_grid(P, H[bad], scaled_jacobian=sj[bad])
        acts.append(dict(mesh=inv, color="#000000", label="inverted cells"))
        vtks["inverted"] = inv
    emit(6, f"quality_{tag}", acts,
         f"06  mesh quality {tag}: min {sj.min():.4f}, mean {sj.mean():.4f}, "
         f"{len(bad)} inverted", vtks)


def step07_sheets(tag="v11"):
    P, H = hexmesh(HEX / f"T1_9_hex_{tag}.ovm")
    f2h, e2h = bc.build_topology(H)
    sing = bc.singular_edges(H, P, f2h, e2h)
    lab, n = cb.bfx.label_sheets(H, f2h, sing)
    sizes = Counter(lab.values())
    print(f"[step07] {tag}: {n} sheets over {len(lab)} faces, sizes "
          f"{sorted(sizes.values(), reverse=True)}")
    quads, sid = [], []
    for fk, s in lab.items():
        quads.append(cb._loop_of(H, f2h[fk][0], fk))
        sid.append(s)
    quads = np.array(quads)
    tri = np.vstack([quads[:, [0, 1, 2]], quads[:, [0, 2, 3]]])
    s = tri_surface(P, tri, sheet_id=np.concatenate([sid, sid]))
    acts = [dict(mesh=s, scalars="sheet_id", cmap="tab20",
                 colorscale="Turbo", label="sheets")]
    emit(7, f"sheets_{tag}", acts,
         f"07  separating sheets {tag}: {n} sheets, {len(lab)} faces",
         {"": s})


def _structure(tag, input_vtk):
    import tet_prep_v5 as v5
    P0, tri0, tid0 = cb.read_input_surface(input_vtk)
    lab = cb.SurfaceLabeller(P0, tri0, tid0, v5.NAMES)
    return cb.BlockStructure(HEX / f"T1_9_hex_{tag}.ovm", lab, fill=True)


def step08_blocks_raw(tag="v11"):
    S = _structure(tag, DATA / "T1_9_tet_v5.vtk")
    cells = S.cells_of()
    order = {r: i for i, r in enumerate(sorted(cells))}
    bid = np.zeros(len(S.hexes), int)
    nf = np.zeros(len(S.hexes), int)
    for r, c in cells.items():
        k = len(S.patches(r))
        for h in c:
            bid[h] = order[r]
            nf[h] = k
    cub = sum(cb.cuboid_status(S.patches(r))[0] == "cuboid" for r in cells)
    print(f"[step08] {tag}: {len(cells)} blocks, {cub} cuboids "
          f"({100 * cub / len(cells):.0f} %), faces per block "
          f"{dict(sorted(Counter(nf.tolist()).items()))}")
    g = hex_grid(S.P, S.hexes, block_id=bid, n_block_faces=nf)
    acts = [dict(mesh=g.extract_surface(), scalars="block_id", cmap="tab20",
                 colorscale="Turbo", label="blocks")]
    emit(8, f"blocks_raw_{tag}", acts,
         f"08  raw base complex {tag}: {len(cells)} blocks, {cub} cuboids",
         {"": g})
    return S


def step09_postprocessing(tag="v11"):
    """What the cut-set cleanup can and cannot do, as n_block_faces."""
    import tet_prep_v5 as v5
    S, before, after, v = cb.postprocess(
        HEX / f"T1_9_hex_{tag}.ovm", DATA / "T1_9_tet_v5.vtk", v5.NAMES,
        verbose=False)
    cells = S.cells_of()
    nf = np.zeros(len(S.hexes), int)
    for r, c in cells.items():
        k = len(S.patches(r))
        for h in c:
            nf[h] = k
    print(f"[step09] {tag}: blocks {before['blocks']} -> {after['blocks']}, "
          f"cuboids {before['cuboids']} -> {after['cuboids']}")
    g = hex_grid(S.P, S.hexes, n_block_faces=nf)
    bad = np.where(nf != 6)[0]
    acts = [dict(mesh=g.extract_surface(), color="#DDDDDD", opacity=0.2,
                 label="all blocks")]
    vtks = {"": g}
    if len(bad):
        b = hex_grid(S.P, S.hexes[bad], n_block_faces=nf[bad])
        acts.append(dict(mesh=b, scalars="n_block_faces", cmap="plasma",
                         colorscale="Plasma", label="non-cuboid"))
        vtks["noncuboid"] = b
    emit(9, f"postprocessing_{tag}", acts,
         f"09  after cut-set cleanup {tag}: {after['blocks']} blocks, "
         f"{after['cuboids']} cuboids, {len(bad)} cells in non-cuboids", vtks)


def step10_sheet_collapse(tag="v11", sheets=(6, 15, 17)):
    """The collapse rounds as a small multiple.

    Replays the sheets the greedy search selected (v11: 6, 15, 17) instead of
    re-running the search, which costs ~15 min per round for a result already
    recorded in PROGRESS.md."""
    import tet_prep_v5 as v5
    P0, tri0, tid0 = cb.read_input_surface(DATA / "T1_9_tet_v5.vtk")
    lab = cb.SurfaceLabeller(P0, tri0, tid0, v5.NAMES)
    S = cb.BlockStructure(HEX / f"T1_9_hex_{tag}.ovm", lab, fill=True)
    P, H = S.P, S.hexes
    acts, vtks = [], {}
    for i, sh in enumerate(sheets):
        sheets_now = cb.mesh_sheets(H)
        if sh >= len(sheets_now):
            print(f"[step10] sheet {sh} out of range, stopping")
            break
        P, H, _drop = cb.collapse_sheet(P, H, sheets_now[sh])
        st = cb._structure_stats(P, H, lab)
        print(f"[step10] round {i + 1} (sheet {sh}): {st}")
        f2h, _ = bc.build_topology(H)
        blocks = bc.blocks_from_cut(H, f2h, set(cb.bfx.label_sheets(
            H, f2h, bc.singular_edges(H, P, f2h, _))[0]))
        bid = np.zeros(len(H), int)
        for bi, b in enumerate(blocks):
            for h in b:
                bid[h] = bi
        g = hex_grid(P, H, block_id=bid)
        vtks[f"round{i + 1}"] = g          # untranslated: the VTK is data
        shown = g.extract_surface().translate((i * 4.2, 0, 0), inplace=False)
        acts.append(dict(mesh=shown, scalars="block_id",
                         cmap="tab20", colorscale="Turbo",
                         label=f"round {i + 1}: {st['blocks']} blocks"))
    emit(10, f"sheet_collapse_{tag}", acts,
         f"10  sheet collapse {tag}: three rounds, 117 -> 31 -> 22 -> 16 blocks",
         vtks)


def step11_final(tag="v11"):
    """The final block structure and the re-attached full domain."""
    import meshio
    d = REPO / "output" / "hex3d_algohex" / "deliverable"
    acts, vtks = [], {}
    f = d / f"T1_9_blocks_{tag}.vtk"
    if f.exists():
        m = meshio.read(f)
        H = np.vstack([b.data for b in m.cells if b.type == "hexahedron"])
        B = np.concatenate([np.asarray(x).ravel() for b, x in
                            zip(m.cells, m.cell_data["block_id"])
                            if b.type == "hexahedron"])
        g = hex_grid(m.points, H, block_id=B)
        print(f"[step11] {tag}: {len(H)} cells, {int(B.max()) + 1} blocks")
        acts.append(dict(mesh=g.extract_surface(), scalars="block_id",
                         cmap="tab20", colorscale="Turbo", label="blocks"))
        vtks[""] = g
        e = d / f"T1_9_blocks_{tag}_edges.vtk"
        if e.exists():
            em = meshio.read(e)
            segs = np.vstack([b.data for b in em.cells if b.type == "line"])
            ls = line_set(em.points, segs)
            acts.append(dict(mesh=ls, kind="lines", color="#000000",
                             line_width=4, label="block edges"))
            vtks["edges"] = ls
    # the re-attached full domain, as its own panel
    full = d / f"T1_9_blocks_{tag}_full.vtk"
    if full.exists():
        fm = meshio.read(full)
        FH = np.vstack([b.data for b in fm.cells if b.type == "hexahedron"])
        FB = np.concatenate([np.asarray(x).ravel() for b, x in
                             zip(fm.cells, fm.cell_data["block_id"])
                             if b.type == "hexahedron"])
        fg = hex_grid(fm.points, FH, block_id=FB)
        print(f"[step11] full domain: {len(FH)} cells, {int(FB.max()) + 1} blocks")
        vtks["full"] = fg
        acts.append(dict(mesh=fg.extract_surface().translate((4.2, 0, 0),
                                                             inplace=False),
                         scalars="block_id", cmap="tab20",
                         colorscale="Turbo", label="full domain"))
    emit(11, f"final_blocks_{tag}", acts,
         f"11  final block structure {tag}: core, and the re-attached "
         f"full domain", vtks)


STEPS = {1: step01_input_surfaces, 2: step02_feature_graph,
         3: step03_domain_variants, 4: step04_singular_graph,
         5: step05_hexmesh, 6: step06_quality, 7: step07_sheets,
         8: step08_blocks_raw, 9: step09_postprocessing,
         10: step10_sheet_collapse, 11: step11_final}


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--steps", default="all",
                    help="comma-separated step numbers, or 'all'")
    ap.add_argument("--tag", default="v11")
    a = ap.parse_args()
    want = (sorted(STEPS) if a.steps == "all"
            else [int(x) for x in a.steps.split(",")])
    for k in want:
        fn = STEPS[k]
        try:
            fn(a.tag) if fn.__code__.co_argcount else fn()
        except Exception as exc:
            print(f"[showcase] step{k:02d} FAILED: {type(exc).__name__}: {exc}")
    print(f"\n[showcase] output in {OUT}")
