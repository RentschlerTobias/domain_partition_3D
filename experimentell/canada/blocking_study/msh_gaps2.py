import sys, meshio, numpy as np
from scipy.spatial import cKDTree
m=meshio.read(sys.argv[1])
names={v[0]:k for k,v in m.field_data.items() if v[1]==2}
groups={}
for cb,ph in zip(m.cells,m.cell_data['gmsh:physical']):
    if cb.type in('triangle','quad'):
        for t in np.unique(ph):
            groups.setdefault(names.get(t,str(t)),[]).append(cb.data[ph==t].ravel())
P=m.points
G={k:P[np.unique(np.concatenate(v))] for k,v in groups.items()}
def merge(prefix): return np.vstack([v for k,v in G.items() if k.startswith(prefix)])
pairs=[p.split(':') for p in sys.argv[2:]]
for a,b in pairs:
    A=merge(a);B=merge(b); d,_=cKDTree(B).query(A);k=d.argmin()
    print(f'{a} - {b}: min dist {d.min():.4f} at {A[k].round(3)}; extent A z [{A[:,2].min():.3f},{A[:,2].max():.3f}]')
