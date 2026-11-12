"""Compare finished simpleFoam runs of the same machine on different meshes.

    foam_compare.py <case_a> <case_b> [...] [--last 100] [--omega 7.5398] [--rho 997]

Per case, from the dtOO function objects (postProcessing/) averaged over the last
--last iterations (and their spread, as a convergence measure):
  Q        volume flow through RU_INLET (m^3/s, of the simulated passage)
  dptot    total pressure drop inlet - outlet (kinematic, m^2/s^2, rho = 1 in ptot)
  Mz       torque about the machine axis on RU_BLADE (forces, rhoInf = --rho)
  P_hyd    rho * dptot * Q,  P_shaft = |Mz| * omega,  eta = P_shaft / P_hyd
  y+       per wall patch (min / max / mean at the last write)
plus the final initial residuals from log.simpleFoam and the cell count.
"""
import argparse
import json
import re
from pathlib import Path

import numpy as np


def series(case, fo, fname, col=1):
    f = sorted((case / "postProcessing" / fo).glob(f"*/{fname}"))[-1]
    rows = [ln.split() for ln in f.read_text().splitlines() if ln.strip() and not ln.startswith("#")]
    return np.array([[float(r[0]), float(r[col])] for r in rows])


def stat(s, last):
    v = s[-last:, 1]
    return float(v.mean()), float(np.ptp(v))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("cases", nargs="+")
    ap.add_argument("--last", type=int, default=100)
    ap.add_argument("--omega", type=float, default=7.5398223686155035)
    ap.add_argument("--rho", type=float, default=997.0)
    ap.add_argument("--json", default=None)
    a = ap.parse_args()
    out = {}
    for c in map(Path, a.cases):
        r = {}
        q, qs = stat(series(c, "Q_ru_in", "surfaceFieldValue.dat"), a.last)
        pi, pis = stat(series(c, "ptot_ru_in", "surfaceFieldValue.dat"), a.last)
        po, pos = stat(series(c, "ptot_ru_out", "surfaceFieldValue.dat"), a.last)
        mz, mzs = stat(series(c, "forces", "moment.dat", col=3), a.last)
        r.update({"Q": abs(q), "Q_spread": qs, "ptot_in": pi, "ptot_out": po, "dptot": pi - po,
                  "dptot_spread": pis + pos, "Mz": mz, "Mz_spread": mzs})
        r["P_hyd"] = a.rho * r["dptot"] * r["Q"]
        r["P_shaft"] = abs(mz) * a.omega
        r["eta"] = r["P_shaft"] / r["P_hyd"] if r["P_hyd"] else float("nan")
        yp = (sorted((c / "postProcessing" / "yPlus1").glob("*/yPlus.dat"))[-1]).read_text().splitlines()
        last_t = [ln.split() for ln in yp if ln.strip() and not ln.startswith("#")]
        tmax = last_t[-1][0]
        r["yplus"] = {p[1]: {"min": float(p[2]), "max": float(p[3]), "mean": float(p[4])}
                      for p in last_t if p[0] == tmax}
        log = (c / "log.simpleFoam").read_text()
        tail = log[log.rindex("\nTime = "):]
        r["final_initial_residuals"] = {m.group(1): float(m.group(2)) for m in
                                        re.finditer(r"Solving for (\w+), Initial residual = ([0-9.eE+-]+)", tail)}
        r["iterations"] = int(tmax)
        m = re.search(r"nCells:\s*(\d+)|cells:\s*(\d+)", (c / "log.checkMesh").read_text()) if (c / "log.checkMesh").exists() else None
        r["cells"] = int(m.group(1) or m.group(2)) if m else None
        out[c.name] = r
    keys = ["cells", "iterations", "Q", "Q_spread", "dptot", "dptot_spread", "Mz", "Mz_spread", "P_hyd", "P_shaft", "eta"]
    names = list(out)
    print(f"{'':16}" + "".join(f"{n:>16}" for n in names) + (f"{'diff b-a %':>14}" if len(names) == 2 else ""))
    for k in keys:
        vals = [out[n][k] for n in names]
        row = f"{k:16}" + "".join(f"{v:16.6g}" if isinstance(v, (int, float)) and v is not None else f"{str(v):>16}" for v in vals)
        if len(names) == 2 and all(isinstance(v, (int, float)) for v in vals) and vals[0]:
            row += f"{100 * (vals[1] - vals[0]) / abs(vals[0]):14.2f}"
        print(row)
    for n in names:
        print(n, "y+ mean", {p: round(v["mean"], 1) for p, v in out[n]["yplus"].items()},
              "residuals", {k: f"{v:.1e}" for k, v in out[n]["final_initial_residuals"].items()})
    if a.json:
        Path(a.json).write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
