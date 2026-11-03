# Keep the 12-block structure, fix the misleading p1 check (2026-11-03)

Question: the P1 acceptance report on the keep-BL reference run
(`/work/datagen_keepbl/runs/tistos_s1_0131`) showed shroud prism-identity
coverage of 0.033634 against a hub value of 0.9987. Fix the boundary-layer
metric by re-blocking the structure along the radial direction, or keep the
12-block structure and fix the check?

## Decision

Keep the 12-block structure and fix the misleading p1 check. Do NOT split
the blocks along the radial (r-band) direction.

What keep-BL requires (user definition): the AlgoHex input covers the full
geometry, min r to max r, with no BL cut; the boundary layer is reproduced
by TFI; prism identity is NOT acceptance.

## The old check measured the wrong thing

Before commit `d43beef` the only gated key was `2_boundary_layers_kept`,
computed from `bl_coverage()`. Its per-surface pass was just
`bool(wall_ok and len(wl) > 0 and blocks_with_bl)` — a wall label, at least
one kept prism, at least one touching block, no coverage threshold — and
`report["pass"] &= bl["pass"]`. The report could therefore say PASS while
never quantifying how much of the shroud BL is reproduced.

## Measurement: 0.033634 is an artefact, not a missing layer

The cause is geometric, not a lost boundary layer: the 12-block structure's
outer faces are planar quads inscribed in r = 1.9 while the kept prism shell
sits at r = 1.897-1.8996, so the prism centroids fall a few 1e-3 outside the
flat faces.

* all 20,543 "outside" shroud prisms lie within 0.0033 of a support cell;
* a tolerance sweep (1e-9 -> 1e-2) moves coverage only 0.036 -> 0.102, so
  the deficit is the flat-face offset, not a missing layer;
* all 97 outer support quads reach r = 1.9, and 11 of the 12 blocks own
  shroud cells (only id 3 is missing);
* TFI rebuilds the 0.027-0.035 thick outer layer.

## The fix (commit d43beef)

* `bl_prism_identity` is now informational only: threshold 0.99, shroud
  coverage 0.033634 -> pass false, `gate` false, and it no longer feeds
  `report["pass"]`.
* `2_boundary_layers_kept` is the gate and measures radial coverage:
  `BL_BAND_HUB = (0.5, 0.5265)`, `BL_BAND_SHROUD = (1.872, 1.9)`; it passes
  iff `min_r <= r_hub + 3e-3`, `max_r >= r_shroud - 3e-3` and both bands lie
  inside `[min_r, max_r]`.
* Regenerated `p1_report.json`: min_r 0.489804, max_r 1.900002, both bands
  inside, `2_boundary_layers_kept` pass true, overall PASS.

## Rejected: the r-band split

Cut the blocks along cylinders r = 0.5265 and r = 1.872 so the kept prisms
sit in their own annular block layer. It touches all 12 block ids (all own
hub cells; the shroud `blocks_with_bl` misses only id 3) and adds at least
one block ring -> >= 24 blocks, outside `BLOCK_MIN, BLOCK_MAX = 10, 16`.
It also violates the user's "no re-blocking of untouched blocks".

## Why

The metric was measuring the wrong thing: prism identity is not acceptance,
and 0.0336 is an inscribed-quad artefact, not a missing BL. Splitting would
cost the block budget and re-block untouched regions for no acceptance gain.

## Consequences

* The P1 gate is now the honest keep-BL check (radial domain coverage);
  prism identity stays as a diagnostic under the informational key.
* Any future "the BL is missing" claim must be shown against
  `2_boundary_layers_kept`, not against `bl_prism_identity`.

## Sources

`d43beef` ("Fix misleading p1 boundary-layers check");
`experimentell/hex3d_algohex/datagen/check_keepbl_sample.py:12-19,378-440,533-564`;
`experimentell/hex3d_algohex/datagen/README.md:77-90,100-112`;
`/work/datagen_keepbl/runs/tistos_s1_0131/p1_report.json`; job tasks T78
(fix the check) and T88 (handoff doc) — the only places this decision is
recorded today.
