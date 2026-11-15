"""gmsh surface+volume meshing test of the simple ring envelopes and a VTK multiblock for inspection."""
import gmsh, meshio, numpy as np, sys
out = "/work/canada/vtk/cad_envelopes_simple"
blocks = []
def dump(tag, key, dim):
    nodes, X, _ = gmsh.model.mesh.getNodes(); idx = {int(n): i for i, n in enumerate(nodes)}
    X = X.reshape(-1, 3); cells, ent = [], []
    for d, t in gmsh.model.getEntities(dim):
        et, _, nt = gmsh.model.mesh.getElements(d, t)
        for e_, n_ in zip(et, nt):
            k = {1: "line", 2: "triangle", 4: "tetra"}.get(e_)
            if k is None: continue
            nn = {"line": 2, "triangle": 3, "tetra": 4}[k]
            c = np.array([idx[int(v)] for v in n_]).reshape(-1, nn); cells.append((k, c)); ent.append(np.full(len(c), t))
    meshio.write(f"{out}/parts/{tag}__{key}.vtu", meshio.Mesh(X, cells, cell_data={"cad_entity": ent}))
    blocks.append((tag, key))
gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0)
for ring in ("ring2", "ring3", "lip"):
    gmsh.clear(); gmsh.model.occ.importShapes(f"{ring}_simple_envelope.brep"); gmsh.model.occ.synchronize()
    gmsh.option.setNumber("Mesh.MeshSizeMax", 0.05)
    gmsh.model.mesh.generate(3)
    ntet = sum(len(t) for t in gmsh.model.mesh.getElements(3)[1])
    q = gmsh.model.mesh.getElementQualities(np.concatenate(gmsh.model.mesh.getElements(3)[1]), "minSICN")
    print(f"{ring}: {len(gmsh.model.getEntities(2))} faces, {len(gmsh.model.getEntities(1))} curves, "
          f"{ntet} tets, minSICN min {q.min():.3f}", flush=True)
    dump(ring, "envelope_faces", 2); dump(ring, "envelope_cad_edges", 1)
    gmsh.clear(); gmsh.model.occ.importShapes(f"{ring}_solid.stl" if False else "/work/canada/case/cad_blk.stp")
    gmsh.model.occ.synchronize()
    keep = {"ring2": [2, 5, 7, 8], "ring3": [3, 6, 10, 11], "lip": [1, 4, 9, 12, 13]}[ring]   # gmsh tags = export index + 1
    gmsh.model.occ.remove([(2, t) for _, t in gmsh.model.getEntities(2) if t not in keep], recursive=True)
    gmsh.model.occ.synchronize(); gmsh.option.setNumber("Mesh.MeshSizeMax", 0.05); gmsh.model.mesh.generate(2)
    dump(ring, "original_ring", 2)
gmsh.finalize()
with open(f"{out}/cad_envelopes_simple.vtm", "w") as f:
    f.write('<?xml version="1.0"?>\n<VTKFile type="vtkMultiBlockDataSet" version="1.0">\n<vtkMultiBlockDataSet>\n')
    for bi, ring in enumerate(("ring2", "ring3", "lip")):
        f.write(f'  <Block index="{bi}" name="{ring}">\n')
        for di, (_, k) in enumerate([b for b in blocks if b[0] == ring]):
            f.write(f'    <DataSet index="{di}" name="{k}" file="parts/{ring}__{k}.vtu"/>\n')
        f.write("  </Block>\n")
    f.write("</vtkMultiBlockDataSet>\n</VTKFile>\n")
