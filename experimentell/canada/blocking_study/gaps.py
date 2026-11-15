"""Per-label-pair minimum distance between surface labels of an AlgoHex input tet.vtk
(triangles carry 'color'). Pairs sharing vertices are measured excluding a
neighbourhood of the shared curve (distance > 3*local edge) to get the true gap."""
import sys, numpy as np, meshio
from scipy.spatial import cKDTree
m = meshio.read(sys.argv[1])
pts = m.points
tris=[];cols=[];tets=[]
for cb, cd in zip(m.cells, m.cell_data['color']):
    if cb.type=='triangle': tris.append(cb.data); cols.append(np.asarray(cd).ravel())
    if cb.type=='tetra': tets.append(cb.data)
T=np.vstack(tris); C=np.concatenate(cols)
tet=np.vstack(tets)
v=pts[tet]; vol=np.abs(np.einsum('ij,ij->i',np.cross(v[:,1]-v[:,0],v[:,2]-v[:,0]),v[:,3]-v[:,0]))/6
V=vol.sum()
print(f'tets {len(tet)}  volume {V:.4f}  bbox {pts.min(0).round(3)} {pts.max(0).round(3)}')
for n in (2000,20000,100000):
    print(f'  n={n}: uniform hex edge (V/n)^(1/3) = {(V/n)**(1/3):.4f}')
labs=np.unique(C)
# sample points: triangle centroids
cent=pts[T].mean(1)
e=np.linalg.norm(pts[T[:,0]]-pts[T[:,1]],axis=1)
h=np.median(e)
vl={l:set(T[C==l].ravel()) for l in labs}
print('labels', {int(l):int((C==l).sum()) for l in labs}, 'median tri edge', round(h,4))
trees={l:cKDTree(cent[C==l]) for l in labs}
res=[]
for i,a in enumerate(labs):
    for b in labs[i+1:]:
        ca=cent[C==a]; d,j=trees[b].query(ca)
        shared = len(vl[a]&vl[b])>0
        if shared:
            # exclude points close to the shared curve
            sv=np.array(list(vl[a]&vl[b])); st=cKDTree(pts[sv])
            ds,_=st.query(ca); ok=ds>4*h
            # also require the closest b point be far from shared curve
            if ok.sum()==0: continue
            db,_=st.query(cent[C==b][j]); ok&=db>4*h
            if ok.sum()==0: continue
            d=d[ok]; ca=ca[ok]
        k=np.argmin(d)
        res.append((d[k],int(a),int(b),shared,ca[k]))
res.sort()
for d,a,b,s,p in res[:12]:
    print(f'  gap {d:.4f}  labels {a}-{b} {"(adjacent)" if s else ""} at {p.round(3)}')
