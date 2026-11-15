"""Regression through the plugins (not the command lines): same numbers as the legacy baselines.

Manual, long (~45 min, 8 cores): needs the baseline directories of /work/hexdp_reg (tistos_0131, canada),
created by running the legacy command lines on copies of the inputs. Not collected by pytest."""
import json, shutil, sys, time
from pathlib import Path
import numpy as np, meshio
from hexdp.contracts import Context, resolve_params
from hexdp.io.gmsh_model import load
from hexdp.io.blocks_vtk import read_blocks_vtk
from hexdp.strategies.algohex.plugin import AlgoHexStrategy
from hexdp.strategies.sweep.plugin import SweepStrategy, CGridStrategy
from hexdp.ops.plugins import BlockCollapseOp
from hexdp.fill.plugins import TfiFiller, B5Filler

R = Path("/work/hexdp_reg"); P = R / "plug"; T = R / "tistos_0131"; C = R / "canada"
which = sys.argv[1:] or ["canada", "tistos"]
out = []


def same_vtk(a, b):
    m, n = meshio.read(a), meshio.read(b)
    ok = m.points.shape == n.points.shape and np.array_equal(m.points, n.points) and \
        all(np.array_equal(x.data, y.data) for x, y in zip(m.cells, n.cells))
    return ok, f"{len(m.points)} pts"


def report(name, ok, info=""):
    out.append((name, ok)); print(f"{'OK  ' if ok else 'DIFF'} {name} {info}", flush=True)


def ctx(name):
    return Context(P / name, cores=8, cpu_list="0-7")


if "canada" in which:
    t = time.time()
    geo = load(C / "inputs" / "domain.brep")
    s = SweepStrategy(); c1 = ctx("sweep")
    print("sweep applicable:", s.applicable(geo, []), flush=True)
    bs = s.run(geo, [], resolve_params(s.params, {}), c1)
    secs = sorted((c1.workspace / "sections").glob("*.json"))
    report("sweep sections", len(secs) == 37 and all(json.load(open(f)) == json.load(open(C / "occ_base/sections" / f.name)) for f in secs), f"{len(secs)}")
    report("sweep layout_final", json.load(open(c1.workspace / "layout_final.json")) == json.load(open(C / "col_base/layout_final.json")))
    report("sweep blocks", *same_vtk(c1.workspace / "sweep_tfi_blocks.vtk", C / "sweep_phi/sweep_tfi_blocks.vtk"))
    print("  ", bs.summary(), flush=True)
    op = BlockCollapseOp(); c2 = ctx("block_collapse")
    print("block_collapse applicable:", op.applicable(bs, geo), flush=True)
    bc = op.apply(bs, geo, resolve_params(op.params, {}), c2)
    d0 = np.load(c2.workspace / "cad.npz"); d1 = np.load(C / "sweep_phi/block_collapse/cad.npz")
    report("cad_export cad.npz", d0.files == d1.files and all(np.array_equal(d0[k], d1[k]) for k in d0.files))
    report("block_collapse tfi", *same_vtk(c2.workspace / "block_collapse/best_h0.0_tfi.vtk", C / "sweep_phi/block_collapse/final/best_h0.0_tfi.vtk"))
    print("  ", bc.summary(), flush=True)
    cg = CGridStrategy(); c3 = ctx("cgrid")
    for f in ("sections", "layout_final.json"):          # same inputs as the sweep: reuse them
        src, dst = c1.workspace / f, c3.workspace / f
        shutil.copytree(src, dst) if src.is_dir() else shutil.copy(src, dst)
    bsc = cg.run(geo, [], resolve_params(cg.params, {"wall_cad": str(C / "inputs/cad_blk.stp")}), c3)
    report("cgrid walls.json", json.load(open(c3.workspace / "cgrid/walls.json")) == json.load(open(C / "occ_base/cgrid/walls.json")))
    report("cgrid blocks", *same_vtk(c3.workspace / "cgrid/sweep_tfi_blocks.vtk", C / "sweep_phi/cgrid/sweep_tfi_blocks.vtk"))
    print("  ", bsc.summary(), f"canada {time.time() - t:.0f}s", flush=True)

if "tistos" in which:
    t = time.time()
    tf = TfiFiller(); c4 = ctx("tfi")
    bs0 = read_blocks_vtk(T / "blocks.vtk")
    tf.fill(bs0, None, resolve_params(tf.params, {"input_vtk": str(T / "tet.vtk")}), c4)
    report("tfi filler", *same_vtk(c4.workspace / "tfi.vtk", T / "base/tfi.vtk"))
    b5 = B5Filler(); c5 = ctx("b5")
    b5.fill(bs0, None, resolve_params(b5.params, {"run_dir": str(T)}), c5)
    report("b5 filler", *same_vtk(c5.workspace / "mesh_b5_core.vtk", T / "base/mesh_b5_core.vtk"))
    geo = load(T / "case" / "mesh.msh")
    ah = AlgoHexStrategy(); c6 = ctx("algohex")
    print("algohex applicable:", ah.applicable(geo, []), flush=True)
    shutil.copy(T / "hex_hex_x.ovm", c6.workspace / "hex_hex_x.ovm")   # AlgoHex itself is not re-run (hours)
    bsa = ah.run(geo, [], resolve_params(ah.params, {}), c6)
    report("algohex tet_prep", *same_vtk(c6.workspace / "tet.vtk", T / "tet.vtk"))
    report("algohex beam", *same_vtk(c6.workspace / "blocks.vtk", T / "base/blocks_beam.vtk"))
    print("  ", bsa.summary(), f"tistos {time.time() - t:.0f}s", flush=True)
print("RESULT", "all identical" if all(ok for _, ok in out) else [n for n, ok in out if not ok])
