import sys; sys.path.insert(0,'/work/canada/cad')
import numpy as np
import domain_blk as db
from simple_envelope import faces, BRep_Tool
for i in (11,12,13,18,19,14,0,3):
    m = db.meridian(faces[i])
    for e in m:
        c,a,b = BRep_Tool.Curve(e)
        P = np.array([[*c.Value(t).Coord()] for t in np.linspace(a,b,5)])
        print(i, np.round(P[:,[0,2]],4).tolist())
