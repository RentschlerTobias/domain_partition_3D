"""2D chord collapse of the swept section layout (greedy, validity-guarded).
A chord = the band of blocks between two neighbouring parallel coordinate lines of the parameter layout.
Collapsing it = merging the two coordinate groups (both get the mean physical coordinate).
Accepted if: parameter polygon valid, arrangement rectangular + T-junction free, block count drops, and the
Winslow re-fit of the sections at the guard phis has 0 inverted quads and min sJ >= GUARD.
usage: python3 collapse2d.py [guard_sJ] [h]  -> layout_final.json, collapse_log.json, layout_*.png"""
import sys, json, time, multiprocessing as mp
sys.path.insert(0, "/work/canada/blocking_study/sweep_phi")
from layout2d import *

GUARD = float(sys.argv[1]) if len(sys.argv) > 1 else 0.1
H = float(sys.argv[2]) if len(sys.argv) > 2 else 0.15
secs = load_sections(); ml = mean_loops(secs); ab = abstract(ml)
CHECK = [k for k, (phi, _) in enumerate(secs) if phi in (0.0, 10.0, 20.0, 30.0, 45.0, 60.0, 90.0)]


def evaluate(merges):
    lay = Layout(ml, ab, merges)
    ok, why = lay.valid()
    if not ok:
        return None, why
    try:
        nb = lay.arrangement()
    except RuntimeError as e:
        return None, str(e)
    if lay.tjunctions:
        return None, f"{lay.tjunctions} T-junctions"
    return lay, nb


def guard(merges):
    lay, nb = evaluate(merges)
    m = Mesh2D(lay, H)
    worst = (1e9, 1.0)
    for k in CHECK:
        X = m.solve(secs[k][1], iters=6000); sj = m.quality(X)
        worst = min(worst, (int((sj <= 0).sum()) * -1, float(sj.min())), key=lambda w: (w[0], w[1]))
        if (sj <= 0).any() or sj.min() < GUARD:
            return merges, nb, False, float(sj.min()), secs[k][0]
    return merges, nb, True, worst[1], None


def candidates(lay):
    out = []
    for ax, gname in ((0, "groups0"), (1, "groups1")):
        g = getattr(lay, gname)
        vals = sorted((lay.U[mem[0], ax], root) for root, mem in g.items())
        for (v0, r0), (v1, r1) in zip(vals[:-1], vals[1:]):
            out.append((abs(v1 - v0), (ax, r0, r1)))
    return [c for _, c in sorted(out)]


def neutral_pass(merges, log):
    # count-neutral pass: merge neighbouring coordinates closer than the mesh size (removes 1-cell slivers
    # from the fine grid without changing the block structure)
    lay, nb = evaluate(merges)
    for c in candidates(lay):
        ax, a, b = c
        if abs(lay.U[a, ax] - lay.U[b, ax]) >= H:
            continue
        l2, n2 = evaluate(merges + [c])
        if l2 is not None and n2 <= nb:
            r = guard(merges + [c])
            print(f"  neutral merge ax {'rz'[ax]} {lay.U[a, ax]:.4f}/{lay.U[b, ax]:.4f}: {'OK' if r[2] else 'REJECT'} min sJ {r[3]:.3f}", flush=True)
            if r[2]:
                merges = merges + [c]; lay, nb = evaluate(merges)
                log.append(dict(step=len(log), blocks=nb, merge=[int(x) for x in c], neutral=True, min_sJ=r[3]))
    return merges


if __name__ == "__main__":
    merges = []; log = []
    lay, nb = evaluate(merges)
    print(f"initial: {nb} blocks, {len(lay.rays)} rays", flush=True)
    log.append(dict(step=0, blocks=nb, merges=[]))
    pool = mp.Pool(8)
    merges = neutral_pass(merges, log)
    while True:
        lay, nb = evaluate(merges)
        tries = []
        for c in candidates(lay):
            l2, n2 = evaluate(merges + [c])
            if l2 is not None and n2 < nb:
                tries.append((n2, c))
        print(f"  {nb} blocks: {len(tries)} count-reducing candidates", flush=True)
        if not tries:
            break
        res = pool.map(guard, [merges + [c] for _, c in tries])
        good = [(r[1], -r[3], r[0]) for r in res if r[2]]
        for r in res:
            ax, a, b = r[0][-1]
            print(f"     merge ax {'rz'[ax]} {lay.U[a, ax]:.3f}/{lay.U[b, ax]:.3f}: blocks {r[1]} "
                  f"{'OK' if r[2] else 'REJECT'} min sJ {r[3]:.3f}" + (f" (phi {r[4]})" if not r[2] else ""), flush=True)
        if not good:
            break
        good.sort()
        n2, negq, mg = good[0]
        c = mg[-1]
        log.append(dict(step=len(log), blocks=n2, merge=[int(c[0]), int(c[1]), int(c[2])],
                        coords=[float(lay.U[c[1], c[0]]), float(lay.U[c[2], c[0]])], min_sJ=-negq))
        merges = mg
        print(f"  accepted: ax {'rz'[c[0]]} {lay.U[c[1], c[0]]:.3f} ~ {lay.U[c[2], c[0]]:.3f} -> {n2} blocks", flush=True)
    merges = neutral_pass(merges, log)
    lay, nb = evaluate(merges)
    json.dump(dict(guard=GUARD, h=H, merges=[[int(x) for x in m_] for m_ in merges], blocks=nb, log=log),
              open(f"{HERE}/layout_final.json", "w"), indent=1)
    print("final", nb, "blocks")
