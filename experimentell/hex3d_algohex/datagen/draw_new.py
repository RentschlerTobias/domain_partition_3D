"""Draw N NEW tistos design vectors (scrambled Sobol, own seed) -> <root>/<prefix>_<i>/params.json.

30-D tistos box, bounds from /dtOO/demo/tistos/build.py.
Reports the distance of every new point to the nearest existing dataset geometry
(params stored in $DATASET_BATCH/*/sample.npz, default /opt/stack/meshtron/data/hex3d_algohex/batch),
in box-normalised units.
usage: draw_new.py <root> <n> <seed> <prefix>
"""
import glob, json, os, re, sys
from pathlib import Path
import numpy as np
from scipy.stats import qmc

root, n, seed, prefix = Path(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
src = open('/dtOO/demo/tistos/build.py').read()
B = {m[0]: (float(m[1]), float(m[2])) for m in
     re.findall(r"'label': '([^']+)', 'min': ([-\d.e]+), 'max': ([-\d.e]+)", src)}
assert len(B) == 30, len(B)
labels = sorted(B)
lo, hi = np.array([B[k][0] for k in labels]), np.array([B[k][1] for k in labels])
unit = qmc.Sobol(d=30, scramble=True, seed=seed).random(n)
old = {}
batch = os.environ.get('DATASET_BATCH', '/opt/stack/meshtron/data/hex3d_algohex/batch')
for f in glob.glob(f'{batch}/machine_*/sample.npz'):
    p = json.loads(str(np.load(f, allow_pickle=True)['params']))
    old[f.split('/')[-2].rsplit('_n', 1)[0]] = [p[k] for k in labels]
O = (np.array(list(old.values())) - lo) / (hi - lo)
d = np.sqrt(((unit[:, None] - O[None]) ** 2).sum(-1)).min(1)
dO = np.sqrt(((O[:, None] - O[None]) ** 2).sum(-1)); np.fill_diagonal(dO, np.inf)
print(f"existing geometries {len(O)}; nearest-neighbour distance among them: median {np.median(dO.min(1)):.3f}")
print(f"new {n}: distance to nearest existing: min {d.min():.3f} median {np.median(d):.3f}")
for i, row in enumerate(qmc.scale(unit, lo, hi)):
    dd = root / f"{prefix}_{i:04d}"; dd.mkdir(parents=True, exist_ok=True)
    (dd / "params.json").write_text(json.dumps({k: round(float(v), 6) for k, v in zip(labels, row)}, indent=1))
print(n, "machines in", root)
