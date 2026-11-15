import gmsh, numpy as np, sys
gmsh.initialize(); gmsh.option.setNumber("General.Verbosity",1)
gmsh.open(sys.argv[1])
for d,t in gmsh.model.getEntities(3):
    et,_,_=gmsh.model.mesh.getElements(d,t)
    up,_=gmsh.model.getAdjacencies(d,t)
    _,_,nt=gmsh.model.mesh.getElements(d,t)
    print(f"V{t}: types={[gmsh.model.mesh.getElementProperties(e)[0].split()[0] for e in et]} n={[len(n)//gmsh.model.mesh.getElementProperties(e)[3] for e,n in zip(et,nt)]} faces={sorted(gmsh.model.getBoundary([(3,t)],oriented=False,combined=False), key=lambda a:a[1])and sorted(abs(b) for _,b in gmsh.model.getBoundary([(3,t)],oriented=False,combined=False))}")
for d,t in gmsh.model.getEntities(2):
    tags,c,_=gmsh.model.mesh.getNodes(2,t,includeBoundary=True); c=c.reshape(-1,3)
    r=np.hypot(c[:,0],c[:,1]); th=np.degrees(np.arctan2(c[:,1],c[:,0]))
    et,_,nt=gmsh.model.mesh.getElements(2,t)
    up,_=gmsh.model.getAdjacencies(2,t)
    print(f"S{t:3d} vols={list(up)} n={sum(len(n)//gmsh.model.mesh.getElementProperties(e)[3] for e,n in zip(et,nt)):6d} {'+'.join(gmsh.model.mesh.getElementProperties(e)[0][:3] for e in et)} r[{r.min():.3f},{r.max():.3f}] z[{c[:,2].min():.3f},{c[:,2].max():.3f}] th[{th.min():.0f},{th.max():.0f}]")
