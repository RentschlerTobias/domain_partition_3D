"""gmsh mesh of the CAD domain (domain_blk.brep) in the format prep_from_geometry.py reads:
  <out>_full.msh  all entities (1D curves = CAD edges, 2D, 3D), Mesh.SaveAll
  <out>_phys.msh  2D physical groups named by label (same entity tags)
Labels come from domain_blk_faces.json (an interior point per face), matched to the gmsh surfaces by
closest-point distance. Surface size: h_env on the envelopes, growing with distance; elsewhere curvature
based up to h_max.

usage: mesh_domain.py <brep> <faces.json> <out_prefix> [--h-env 0.03] [--h-max 0.15]
"""
import argparse, json
import numpy as np
import gmsh

ap = argparse.ArgumentParser()
ap.add_argument("brep"); ap.add_argument("faces"); ap.add_argument("out")
ap.add_argument("--h-env", type=float, default=0.03)
ap.add_argument("--h-max", type=float, default=0.15)
ap.add_argument("--threads", type=int, default=4)
a = ap.parse_args()

gmsh.initialize(); gmsh.option.setNumber("General.Verbosity", 2)
gmsh.option.setNumber("General.NumThreads", a.threads)
gmsh.model.occ.importShapes(a.brep); gmsh.model.occ.synchronize()
rec = json.load(open(a.faces))
surf = [t for _, t in gmsh.model.getEntities(2)]
lab = {}
for r in rec:
    d = {t: np.linalg.norm(np.array(gmsh.model.getClosestPoint(2, t, r["point"])[0]) - r["point"]) for t in surf}
    t = min(d, key=d.get)
    assert d[t] < 1e-5 and t not in lab, (r, d[t])
    lab[t] = r["label"]
assert len(lab) == len(surf), "unlabelled surfaces"
groups = {}
for t, l in lab.items():
    groups.setdefault(l, []).append(t)
for l, ts in sorted(groups.items()):
    gmsh.model.setPhysicalName(2, gmsh.model.addPhysicalGroup(2, ts), l)
print("surfaces per label:", {l: len(ts) for l, ts in sorted(groups.items())}, flush=True)

env = [t for t, l in lab.items() if l.startswith("ENV_")]
f_d = gmsh.model.mesh.field.add("Distance"); gmsh.model.mesh.field.setNumbers(f_d, "SurfacesList", env)
gmsh.model.mesh.field.setNumber(f_d, "Sampling", 200)
f_t = gmsh.model.mesh.field.add("Threshold")
gmsh.model.mesh.field.setNumber(f_t, "InField", f_d)
gmsh.model.mesh.field.setNumber(f_t, "SizeMin", a.h_env); gmsh.model.mesh.field.setNumber(f_t, "SizeMax", a.h_max)
gmsh.model.mesh.field.setNumber(f_t, "DistMin", 0.05); gmsh.model.mesh.field.setNumber(f_t, "DistMax", 1.0)
gmsh.model.mesh.field.setAsBackgroundMesh(f_t)
gmsh.option.setNumber("Mesh.MeshSizeMax", a.h_max)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 24)
gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 1)
gmsh.option.setNumber("Mesh.Algorithm", 6); gmsh.option.setNumber("Mesh.Algorithm3D", 10)
gmsh.model.mesh.generate(3)
ntri = sum(len(x) for x in gmsh.model.mesh.getElements(2)[1])
ntet = sum(len(x) for x in gmsh.model.mesh.getElements(3)[1])
print(f"mesh: {ntri} boundary triangles, {ntet} tets, {len(gmsh.model.getEntities(1))} curves", flush=True)
gmsh.option.setNumber("Mesh.SaveAll", 1); gmsh.write(a.out + "_full.msh")
gmsh.option.setNumber("Mesh.SaveAll", 0); gmsh.write(a.out + "_phys.msh")
gmsh.finalize()
