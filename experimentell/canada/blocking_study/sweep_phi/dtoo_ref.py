import gmsh, numpy as np, collections
gmsh.initialize(); gmsh.option.setNumber("General.Verbosity", 1)
gmsh.open("/work/canada/case/E1_12685_ingvrudtout_blkGridMesh.msh")
out = []
for d, t in gmsh.model.getEntities(3):
    ty, et, nt = gmsh.model.mesh.getElements(3, t)
    nodes = np.unique(np.concatenate(nt))
    xyz = np.array([gmsh.model.mesh.getNode(int(n))[0] for n in nodes[::max(1, len(nodes)//400)]])
    phi = np.degrees(np.arctan2(xyz[:,1], xyz[:,0])); r = np.hypot(xyz[:,0], xyz[:,1])
    kinds = {gmsh.model.mesh.getElementProperties(int(x))[0]: len(e) for x, e in zip(ty, et)}
    out.append((t, kinds, phi.min(), phi.max(), r.min(), r.max(), xyz[:,2].min(), xyz[:,2].max()))
c = collections.Counter((round(o[2]/5)*5, round(o[3]/5)*5) for o in out)
print("phi ranges (rounded):", c)
for o in out[:200]:
    print(o[0], o[1], "phi %.0f-%.0f r %.2f-%.2f z %.2f-%.2f" % o[2:])
