import sys, gmsh, numpy as np
from scipy.spatial import cKDTree
gmsh.initialize(); gmsh.option.setNumber('General.Verbosity',0)
gmsh.open(sys.argv[1])
groups={}
for dim,tag in gmsh.model.getPhysicalGroups(2):
    name=gmsh.model.getPhysicalName(dim,tag)
    pts=[]
    for e in gmsh.model.getEntitiesForPhysicalGroup(dim,tag):
        _,c,_=gmsh.model.mesh.getNodes(2,e,includeBoundary=True)
        pts.append(c.reshape(-1,3))
    if pts: groups[name]=np.vstack(pts)
print({k:len(v) for k,v in groups.items()})
want=sys.argv[2].split(',')
for i,a in enumerate(want):
    for b in want[i+1:]:
        if a in groups and b in groups:
            d,_=cKDTree(groups[b]).query(groups[a]); k=d.argmin()
            print(f'{a} - {b}: min {d.min():.4f} at {groups[a][k].round(3)}')
