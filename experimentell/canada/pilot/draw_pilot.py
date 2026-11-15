"""Draw N canada design vectors (scrambled Sobol) in the box around state E1_12685.

Box: runner + guide-vane blade parameters (cV_ru_*, cV_gv_* without mesh-control ones),
value +-15 % of |value|, +-0.05 where the value is exactly 0. Channel constants stay fixed.
usage: draw_pilot.py <root> <n> <seed>
"""
import json, re, sys
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
from scipy.stats import qmc

root, n, seed = Path(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
r = ET.parse('/work/canada/case/E1_12685.xml').getroot()
s = [s for s in r.iter('state') if s.get('label') == 'E1_12685'][0]
v = {c.get('label'): float(c.get('value')) for c in s.iter() if c.get('label', '').startswith(('cV_ru_', 'cV_gv_'))}
mesh = re.compile(r'nElem|grid|Grading|firstElement|nBoundary|charC?L|lcInt|meshBlock|divide|stepN|mplane|nBlades'
                  r'|_in$|_out$|_le$|_te$|le>te|ruFine|gvFine')
labels = sorted(k for k in v if not mesh.search(k))
c = np.array([v[k] for k in labels])
w = np.where(c == 0, 0.05, 0.15 * np.abs(c))
lo, hi = c - w, c + w
(root).mkdir(parents=True, exist_ok=True)
json.dump({k: [float(a), float(b), float(x)] for k, a, b, x in zip(labels, lo, hi, c)},
          open(root / "box.json", "w"), indent=1)
for i, row in enumerate(qmc.scale(qmc.Sobol(d=len(labels), scramble=True, seed=seed).random(n), lo, hi)):
    d = root / "runs" / f"canada_p_{i:04d}"
    d.mkdir(parents=True, exist_ok=True)
    (d / "params.json").write_text(json.dumps({k: round(float(x), 7) for k, x in zip(labels, row)}, indent=1))
print(n, "samples,", len(labels), "parameters ->", root)
