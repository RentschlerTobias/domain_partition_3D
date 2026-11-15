import meshio, numpy as np, sys
m=meshio.read(sys.argv[1]); p=m.points
for blk,data in zip(m.cells, m.cell_data[list(m.cell_data)[0]]):
    if blk.type!="triangle": continue
    t=blk.data; lab=np.asarray(data)
    for l in np.unique(lab):
        c=p[t[lab==l]].reshape(-1,3); r=np.hypot(c[:,0],c[:,1])
        a=0.5*np.linalg.norm(np.cross(p[t[lab==l]][:,1]-p[t[lab==l]][:,0],p[t[lab==l]][:,2]-p[t[lab==l]][:,0]),axis=1).sum()
        print(f"label {l}: {np.sum(lab==l):5d} tris area {a:.4f} r[{r.min():.3f},{r.max():.3f}] z[{c[:,2].min():.3f},{c[:,2].max():.3f}]")
print("cell_data keys", list(m.cell_data))
