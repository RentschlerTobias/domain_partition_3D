#!/usr/bin/env python3
"""Export the feature curves of the keep-BL tistos tetmesh input (tet.vtk) as VTK lines.

The source file (/work/datagen_keepbl/runs/tistos_s1_0131/tet.vtk, written by
tet_prep_v5 --keep-prisms --remall 0) is a legacy ASCII VTK UNSTRUCTURED_GRID:
  - CELL_TYPES 1 = vertex, 3 = edge (VTK_LINE), 5 = triangle, 10 = tetrahedron
  - CELL_DATA 'color': nonzero = feature element (feature verts/edges/tris + tets)

What the script does
  1. parses the source VTK, counts the feature edges (color != 0 among the edge
     cells) and groups them into connected curves (components sharing
     vertices), reporting the number of curves, their sizes and open/closed;
  2. counts the feature points and feature surfaces as well;
  3. exports the feature curves as explicit VTK_LINE cells plus the feature
     points as vertex cells to /work/datagen_keepbl/inspect/keepbl_features.vtk
     with a CELL_DATA 'feature_curve_id' per cell (vertices: -1), and prints
     the file size. All source points are copied so the exported lines keep
     the source vertex indices -- the edge sets then compare exactly;
  4. for context, parses the stripped reference
     /work/datagen/runs/tistos_s1_0131/tet.vtk (if present) and reports its
     feature-edge/curve counts next to the keep-BL ones;
  5. --check <export.vtk> re-parses the exported VTK and exits 0 only if the
     exported edge set and feature_curve_id assignment are exactly equal to
     the source ones (exact set equality, not just counts).

Read-only w.r.t. the run dirs; output goes to /work/datagen_keepbl/inspect/.

Usage:
  export_features_vtk.py [export.vtk]                # export + report
  export_features_vtk.py --check export.vtk          # verify export vs source
"""
import argparse
import collections
import os
import sys

DEFAULT_SOURCE = "/work/datagen_keepbl/runs/tistos_s1_0131/tet.vtk"
DEFAULT_EXPORT = "/work/datagen_keepbl/inspect/keepbl_features.vtk"
DEFAULT_REFERENCE = "/work/datagen/runs/tistos_s1_0131/tet.vtk"
DEFAULT_INSPECT_DIR = "/work/datagen_keepbl/inspect"

VTK_VERTEX = 1
VTK_LINE = 3
VTK_TRIANGLE = 5


def read_points(path, npoints_hint=None):
    """Return the source POINTS section as a list of string triples."""
    triples = []
    with open(path) as f:
        lines = f.read().split("\n")
    for i, l in enumerate(lines):
        if l.startswith("POINTS"):
            npoints = int(l.split()[1])
            k = i + 1
            seen = 0
            while seen < npoints:
                toks = lines[k].split()
                k += 1
                if not toks:
                    continue
                for q in range(0, len(toks), 3):
                    triples.append(toks[q:q + 3])
                    seen += 1
            return triples[:npoints]
    raise ValueError("%s: no POINTS section" % path)


def parse_vtk(path):
    """Parse a legacy ASCII VTK UNSTRUCTURED_GRID with per-cell int scalars.

    The first SCALARS block under CELL_DATA is read into 'colors'
    (one int per cell; None if the file has no COLS/SCALARS).
    """
    with open(path) as f:
        lines = f.read().split("\n")
    n_lines = len(lines)

    def next_line(i):
        while i < n_lines and lines[i].strip() == "" or (i < n_lines and lines[i].startswith("#")):
            i += 1
        return i

    def expect(i, prefix, what):
        if i >= n_lines or not lines[i].startswith(prefix):
            got = lines[i] if i < n_lines else "<eof>"
            raise ValueError("%s: expected %s at line %d, got %r"
                             % (path, what, i + 1, got))
        return lines[i]

    i = next_line(0)
    # legacy header = 2 lines (magic + title) before DATASET; scan to POINTS
    while i < n_lines and not lines[i].startswith("POINTS"):
        i = next_line(i + 1)
    expect(i, "POINTS", "POINTS header")
    npoints = int(lines[i].split()[1])
    point_start = i + 1
    point_end = point_start
    consumed = 0
    while consumed < 3 * npoints:
        if point_end >= n_lines:
            raise ValueError("%s: truncated POINTS" % path)
        toks = lines[point_end].split()
        if toks:
            consumed += len(toks)
        point_end += 1

    i = next_line(point_end)
    expect(i, "CELLS", "CELLS header")
    ncells, nints = (int(x) for x in lines[i].split()[1:3])
    cells = []
    j = i + 1
    for k in range(ncells):
        parts = lines[j + k].split()
        cells.append([int(x) for x in parts[1:1 + int(parts[0])]])

    j = next_line(j + ncells)
    expect(j, "CELL_TYPES", "CELL_TYPES header")
    types = [int(lines[j + 1 + k].split()[0]) for k in range(ncells)]

    colors = [None] * ncells
    k = j + 1 + ncells
    # scan forward for CELL_DATA (types section has exactly ncells int lines;
    # scan anyway so trailing whitespace or blanks do not break the lookup)
    while k < n_lines and not lines[k].startswith("CELL_DATA"):
        s = lines[k].strip()
        if s and not s.lstrip("+-").isdigit():
            break
        k += 1
    if k < n_lines and lines[k].startswith("CELL_DATA"):
        if int(lines[k].split()[1]) != ncells:
            raise ValueError("%s: CELL_DATA count mismatch" % path)
        k = next_line(k + 1)
        if k < n_lines and lines[k].startswith("SCALARS"):
            k = next_line(k + 1)
            tokens = []
            while k < n_lines:
                s = lines[k].strip()
                if not s or s.startswith(("POINT_DATA", "CELL_DATA", "VECTORS")):
                    break
                tokens.extend(s.split())
                k += 1
            # skip the LOOKUP_TABLE clauses that can be interleaved
            vals = []
            it = iter(tokens)
            it_skip = False
            for tok in it:
                if tok == "LOOKUP_TABLE":
                    it_skip = True
                    continue
                if it_skip and tok == "default":
                    it_skip = False
                    continue
                if tok == "default":
                    continue
                vals.append(int(tok))
            if len(vals) != ncells:
                raise ValueError("%s: got %d scalar values for %d cells"
                                 % (path, len(vals), ncells))
            colors = vals
    return {
        "npoints": npoints,
        "cells": cells,
        "types": types,
        "colors": colors,
    }


def feature_info(parsed):
    """Feature edges (as sorted (a,b) tuples), points, triangles, surfaces."""
    cells, types, colors = parsed["cells"], parsed["types"], parsed["colors"]
    edges = []
    points = []
    tris = 0
    surface_ids = set()
    all_ids = set()
    for cell, typ, col in zip(cells, types, colors):
        if col is None or col == 0:
            continue
        all_ids.add(col)
        if typ == VTK_LINE:
            edges.append(tuple(sorted(cell)))
        elif typ == VTK_VERTEX:
            points.append(cell[0])
        elif typ == VTK_TRIANGLE:
            tris += 1
            surface_ids.add(col)
    return {
        "edges": set(edges),
        "edge_count": len(edges),
        "points": sorted(set(points)),
        "tri_count": tris,
        "surface_ids": sorted(surface_ids),
        "all_ids": sorted(all_ids),
    }


def group_curves(edges):
    """Group a set of (a, b) edges into connected curves (components)."""
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in edges:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra
    member_edges = collections.defaultdict(list)
    deg = collections.defaultdict(collections.Counter)
    for e in sorted(edges):
        member_edges[find(e[0])].append(e)
        deg[find(e[0])][e[0]] += 1
        deg[find(e[0])][e[1]] += 1
    curves = []
    for root, cedges in member_edges.items():
        closed = all(d == 2 for d in deg[root].values())
        curves.append({"edges": sorted(cedges), "size": len(cedges), "closed": closed})
    curves.sort(key=lambda c: c["edges"])
    return curves


def report(tag, info, curves):
    print("[ %s ]" % tag)
    print("  feature edges      : %d" % info["edge_count"])
    print("  feature points (0d): %d" % len(info["points"]))
    print("  feature tris (2d)  : %d" % info["tri_count"])
    print("  feature surface ids: %s  (%d surfaces)" %
          (info["surface_ids"], len(info["surface_ids"])))
    print("  connected curves   : %d" % len(curves))
    for k, c in enumerate(curves, 1):
        lo = min(v for e in c["edges"] for v in e)
        hi = max(v for e in c["edges"] for v in e)
        print("    curve %3d: %4d edges, %-6s, vertex id range %d..%d" %
              (k, c["size"], "closed" if c["closed"] else "open", lo, hi))


def write_export(path, source_path, source_line, info, curves):
    """Write the export: all source points + feature vertices/lines + curve ids.

    Every source point is copied (same order), so the exported lines keep the
    source vertex indices and the check can compare the edge sets exactly.
    """
    lines = ["# vtk DataFile Version 2.0",
             source_line,
             "ASCII",
             "DATASET UNSTRUCTURED_GRID"]
    dp = read_points(source_path)
    n_points = len(dp)
    n_vcells = len(info["points"])
    n_edges = sum(c["size"] for c in curves)
    n_cells = n_vcells + n_edges
    conn = n_vcells + 2 * n_edges
    lines.append("POINTS %d float" % n_points)
    lines.extend(" ".join(t) for t in dp)
    lines.append("CELLS %d %d" % (n_cells, conn))
    for v in info["points"]:
        lines.append("1 %d" % v)
    for c in curves:
        for a, b in c["edges"]:
            lines.append("2 %d %d" % (a, b))
    lines.append("CELL_TYPES %d" % n_cells)
    lines.extend(str(VTK_VERTEX) for _ in range(n_vcells))
    lines.extend(str(VTK_LINE) for _ in range(n_edges))
    lines.append("CELL_DATA %d" % n_cells)
    lines.append("SCALARS feature_curve_id int 1")
    lines.append("LOOKUP_TABLE default")
    lines.extend(str(-1) for _ in info["points"])
    for c in curves:
        lines.extend(str(c["curve_id"]) for _ in c["edges"])
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    return n_cells


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("path", nargs="?", default=DEFAULT_EXPORT,
                    help="export target (default mode) or the exported file to check (--check);"
                         " default %(default)s")
    ap.add_argument("--source", default=DEFAULT_SOURCE,
                    help="source tet.vtk (default: %(default)s)")
    ap.add_argument("--reference", default=DEFAULT_REFERENCE,
                    help="stripped reference tet.vtk for context; pass '' to skip (default: %(default)s)")
    ap.add_argument("--check", action="store_true",
                    help="re-parse the exported file and verify the feature curve graph")
    args = ap.parse_args(argv)

    dp = read_points(args.source)
    info = feature_info(parse_vtk(args.source))
    curves = group_curves(info["edges"])
    for k, c in enumerate(curves, 1):
        c["curve_id"] = k
    report("SOURCE %s" % args.source, info, curves)

    if args.reference and os.path.exists(args.reference):
        rinfo = feature_info(parse_vtk(args.reference))
        rcurves = group_curves(rinfo["edges"])
        report("CONTEXT (stripped reference) %s" % args.reference, rinfo, rcurves)
        print("  keep-BL surface tags on feature tris: %s (%d) vs reference: %s (%d)"
              % (info["surface_ids"], len(info["surface_ids"]),
                 rinfo["surface_ids"], len(rinfo["surface_ids"])))
        extra = [s for s in info["surface_ids"] if s not in rinfo["surface_ids"]]
        gone = [s for s in rinfo["surface_ids"] if s not in info["surface_ids"]]
        print("  encoded surface tags (any feature cell): keep-BL %s vs reference %s"
              % (info["all_ids"], rinfo["all_ids"]))
        print("  extra keep-BL tags: %s%s -> those surfaces' edges/tris are part of"
              " the keep-BL feature graph, so the keep-BL edge/curve counts differ"
              " from the reference's; expected, not an export artefact."
              % (extra, "; tags %s of the reference have no feature tris in"
                        " keep-BL" % gone if gone else ""))

    edge_map = {}
    for c in curves:
        for e in c["edges"]:
            edge_map[e] = c["curve_id"]

    if args.check:
        exp_path = args.path
        ep = parse_vtk(exp_path)
        if ep["npoints"] != len(dp):
            print("MISMATCH: %d points in the export vs %d in the source"
                  % (ep["npoints"], len(dp)))
            print("CHECK FAILED: %s" % exp_path)
            return 1
        cells, types, colors = ep["cells"], ep["types"], ep["colors"]
        if any(c is None for c in colors):
            print("MISMATCH: no per-cell feature_curve_id in the export")
            print("CHECK FAILED: %s" % exp_path)
            return 1
        exp_points = sorted(c[0] for c, t in zip(cells, types) if t == VTK_VERTEX)
        exp_curves = collections.defaultdict(set)
        for c, t, cid in zip(cells, types, colors):
            if t == VTK_LINE:
                e = tuple(sorted(c))
                exp_curves[cid].add(e)
        exp_edges = set().union(*exp_curves.values()) if exp_curves else set()
        ok = True
        if exp_points != info["points"]:
            ok = False
            print("MISMATCH: feature points: %d exported vs %d in source"
                  % (len(exp_points), len(info["points"])))
        if exp_edges != info["edges"]:
            ok = False
            print("MISMATCH: edge set: %d exported vs %d in source;" % (len(exp_edges), len(info["edges"])))
            print("  only in export : %s" % sorted(exp_edges - info["edges"])[:10])
            print("  only in source : %s" % sorted(info["edges"] - exp_edges)[:10])
        else:
            print("EDGE SET: exact equality, %d feature edges" % len(exp_edges))
        src_part = {cid: tuple(sorted(curves[cid - 1]["edges"])) for cid in edge_map.values()}
        exp_part = {cid: tuple(sorted(es)) for cid, es in exp_curves.items()}
        if exp_part != src_part:
            ok = False
            bad = [cid for cid in exp_part if exp_part[cid] != src_part.get(cid)]
            print("MISMATCH: %d curve ids with different edge contents (ids: %s)"
                  % (len(bad), bad[:10]))
        else:
            print("CURVES: exact equality, %d curves by feature_curve_id" % len(exp_part))
        if not ok:
            print("CHECK FAILED: %s" % exp_path)
            return 1
        print("CHECK OK: exported feature curve graph == source feature curve graph")
        return 0

    # export mode
    source_line = "keep-BL tistos feature curves exported from %s" % os.path.basename(args.source)
    write_export(args.path, args.source, source_line, info, curves)
    size = os.path.getsize(args.path)
    print("EXPORTED %s (%d bytes)" % (args.path, size))
    print("  %d vertex cells + %d line cells, CELL_DATA feature_curve_id; %d curves"
          % (len(info["points"]), sum(c["size"] for c in curves), len(curves)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
