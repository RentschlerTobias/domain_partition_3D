# LE/TE feature edges: measured and rejected (tistos family, 2026-10-03)

Question (from the singularity-graph VTKs): two singularity arcs run quasi
parallel to the blade LE/TE; LE/TE are not in the feature graph. Would
adding them simplify the graph or fix the historical full-domain
non-convergence?

## What the LE/TE actually are

The blade surface is built closed over the profile direction
(`closeU=true` in `xml/ru_bladeRunner2d.xml`); LE = u=0, TE = u=1 per the
`cV_ru_u_le/u_te` constants. The LE/TE are therefore PARAMETRIC boundary
curves, not geometric creases: the wall's internal dihedrals max out at
17 deg, so the 40 deg sharp-crease detector can never find them. In the
full-domain meshing the blade wall was already a feature SURFACE (all
boundary faces coloured); in the reduced dataset path the blade is cut
out entirely and only the two cut rings survive as features.

## Measurement

| run | setup | outcome |
|---|---|---|
| B1 | full domain, wall as feature surface, blade patch SPLIT at the projected LE/TE (5/8) -> LE/TE as feature edges | **non-converged**: 2 d 6 h, ~976k Newton rounds, f(x) stuck in the 1e2-4.9e3 range, 86 degenerate-cell warnings, quantization stage never resolves ("path constraints = 0" forever), no OVM |
| B2 | reduced domain, ogrid_interface split (8/9) along the projected LE/TE | parameterization + quantization completed, then silent crash at HexEx edge extraction ("Processing edge 0 of 210205"), no OVM, 2 h |
| base_a | reduced, no LE/TE | converged ~15 min (the dataset standard) |
| R1 | reduced + ring densify (--remesh-ogrid 0.030) | converged; 12 blocks / 12 cuboids / 0 inverted / no kink > 30 deg |

Full report with log evidence: `/work/analysis/le_te_report.md` (sandbox).

## Verdict

- The two LE/TE-parallel singularity arcs are the normal field topology,
  NOT a symptom of missing features.
- LE/TE feature edges are net-harmful in both paths. The reduced path
  without them stays the dataset standard; the full-domain route stays
  closed for now (its blocker is the blade-wall feature surface itself,
  not the missing edges).
- The ring-densify variant R1 remains the structural quality fix of this
  round (defect block healed).

Tooling kept: `/work/perturb/split_blade_le_te.py` (patch-split at
projected parametric curves; reusable for future constraint experiments),
`/work/perturb/le_te.npz` (the dtOO-sampled LE/TE curves).
