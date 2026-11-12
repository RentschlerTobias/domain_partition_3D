"""The dtOO mesh (case/mesh.msh) as an OpenFOAM case with the dtOO patch names.

    dtoo_foam.py <run_dir> --case <dir>

gmshToFoam -keepOrientation (without it gmshToFoam "inverts" valid but strongly
skewed O-grid hexes: 321 negative volumes on 0131; dtOO itself sets keepOrientation
in xml/tistos_ru_of.xml). The msh carries physical tag 0 on every element, so
gmshToFoam writes one patch; its faces are split here by the geometric entity of
the msh boundary element they come from, grouped as geometry.py does (hub, shroud,
blade, inlet, outlet, periodic_A/B) and named as mesh_b5.PATCHES. The cyclicAMI
pair is written as export_foam does.
"""
import argparse
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))

from geometry import Geometry                                         # noqa: E402
from mesh_b5 import PATCHES                                           # noqa: E402
from export_foam import write_boundary                                # noqa: E402

FOAM_BASHRC = "/usr/lib/openfoam/openfoam2606/etc/bashrc"


def read_list(path):
    """(header, entries) of an ascii OpenFOAM list file."""
    t = path.read_text()
    m = re.compile(r"\n(\d+)\s*\n\(\n").search(t, t.index("}", t.index("FoamFile")))
    return t[:m.start() + 1], [ln for ln in t[m.end():t.rindex(")")].split("\n") if ln.strip()]


def write_list(path, header, entries):
    path.write_text(header + f"{len(entries)}\n(\n" + "\n".join(entries) + "\n)\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run")
    ap.add_argument("--case", required=True)
    a = ap.parse_args()
    run, case = Path(a.run).resolve(), Path(a.case).resolve()
    (case / "system").mkdir(parents=True, exist_ok=True)
    (case / "system" / "controlDict").write_text(
        "FoamFile { version 2.0; format ascii; class dictionary; object controlDict; }\n"
        "application simpleFoam; startTime 0; endTime 1; deltaT 1; writeControl timeStep;\n"
        "writeInterval 1; writeFormat ascii;\n")
    subprocess.run(["bash", "-c", f"source {FOAM_BASHRC} >/dev/null 2>&1; cd '{case}'; "
                    f"gmshToFoam -keepOrientation '{run / 'case' / 'mesh.msh'}' > log.gmshToFoam 2>&1"], check=True)
    pm = case / "constant" / "polyMesh"
    bnd = (pm / "boundary").read_text()
    nf = int(re.search(r"nFaces\s+(\d+)", bnd).group(1))
    sf = int(re.search(r"startFace\s+(\d+)", bnd).group(1))
    fh, flines = read_list(pm / "faces")
    oh, olines = read_list(pm / "owner")
    faces = [list(map(int, re.match(r"\d+\(([\d ]+)\)", ln).group(1).split())) for ln in flines[sf:sf + nf]]
    # boundary element -> surface name by geometric entity
    import meshio
    m = meshio.read(str(run / "case" / "mesh.msh"))
    G = Geometry(run)
    ent2name = {t: n for n, ts in G.groups.items() for t in ts}
    key2name = {}
    for c, g in zip(m.cells, m.cell_data["gmsh:geometrical"]):
        if c.type in ("triangle", "quad"):
            for e, t in zip(c.data.tolist(), np.asarray(g).ravel().tolist()):
                if t in ent2name:
                    key2name[tuple(sorted(e))] = ent2name[t]
    names = list(PATCHES)
    pid = np.array([names.index(key2name[tuple(sorted(f))]) if tuple(sorted(f)) in key2name else -1
                    for f in faces])
    if (pid < 0).any():
        raise SystemExit(f"[dtoo_foam] {int((pid < 0).sum())} boundary faces without a surface")
    order = np.argsort(pid, kind="stable")
    flines[sf:sf + nf] = [flines[sf + i] for i in order]
    olines[sf:sf + nf] = [olines[sf + i] for i in order]
    write_list(pm / "faces", fh, flines)
    write_list(pm / "owner", oh, olines)
    patches = [(PATCHES[names[k]], int((pid == k).sum())) for k in sorted(set(pid.tolist()))]
    write_boundary(pm, patches, sf)
    print("[dtoo_foam]", dict(patches))


if __name__ == "__main__":
    main()
