import json, glob, numpy as np
def ang(a, b):  # interior angle at the joint: incoming direction a, outgoing b (CCW loop, domain on the left)
    t = np.degrees(np.arctan2(a[0]*b[1]-a[1]*b[0], a@b))   # turn angle
    return 180 - t
for f in sorted(glob.glob("sections/sec_*.json")):
    S = json.load(open(f)); L = S["loops"][0]
    P = [np.array(d["pts"]) for d in L]
    # orientation: signed area
    allp = np.vstack([p[:-1] for p in P]); A = 0.5*np.sum(allp[:,0]*np.roll(allp[:,1],-1)-np.roll(allp[:,0],-1)*allp[:,1])
    out = []
    for k in range(len(P)):
        a = P[k][-1]-P[k][-4]; b = P[(k+1)%len(P)][3]-P[(k+1)%len(P)][0]
        th = ang(a/np.linalg.norm(a), b/np.linalg.norm(b))
        if A < 0: th = 360 - th
        out.append(f"{L[k]['face']}/{L[(k+1)%len(P)]['face']}:{th:.0f}")
    sl = [(d["face"], np.round(d["pts"][0],3).tolist(), np.round(d["pts"][-1],3).tolist()) for d in S["slits"]]
    print(S["phi"], "area", round(A,2), " ".join(out)); print("   slits", sl)
