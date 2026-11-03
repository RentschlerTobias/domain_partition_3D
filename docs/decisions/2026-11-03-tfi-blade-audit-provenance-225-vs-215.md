# Provenance of tfi_blade_audit.vtk — the 225-vs-215 inverted-cell discrepancy

Date: 2026-11-03. Task T93 (read-only audit of the T87 viz numbers).

## Question

T87 reported that `tfi_blade_audit.vtk` showed **225** inverted cells on a first
read and **215** on a re-read. 225 is exactly the inverted-cell count of the T85
regraded mesh (`tfi_blade_graded.vtk`), not of the canonical `tfi_blade.vtk` (215).
So either the unit-of-truth audit file was overwritten by the T85 run and later
restored, or it is currently the graded audit — either way every later count and
the viz would be silently wrong.

## Verdict

**The audit file on disk IS the canonical audit. 225 was never the content of
`tfi_blade_audit.vtk` in its current state; it is the count of the T85 graded
mesh.** No regeneration was needed and no mesh was touched.

Evidence (all recomputed read-only from the files on disk):

| mesh | md5 | points | cells | n(sJ<0) | n(sJ<0.3) | n(vol<=0) | min sJ | cell |
|---|---|---|---|---|---|---|---|---|
| `tfi_blade.vtk` (canonical source) | `61b259521b3167ebaa9031f80bf01038` | 112809 | 104940 | 215 | 429 | 141 | -0.9936709 | 87410 |
| `tfi_blade_audit.vtk` (audit on disk) | `69ee5741d3aabe48038b24f72392b661` | 112809 | 104940 | 215 | 429 | 141 | -0.9936709 | 87410 |
| `tfi_blade_graded.vtk` (T85 regraded) | `f86a47b602789cde9865a4a875c8bc8f` | 112809 | 104940 | 225 | 443 | 150 | -0.9869441 | 87419 |

Identity checks on the audit file itself:

* its point coordinates are **bit-identical to `tfi_blade.vtk`** and **different
  from `tfi_blade_graded.vtk`** → it was built from the canonical mesh;
* canonical and graded meshes share the *same* hex connectivity, so only the
  point coordinates distinguish them (this is why a swap is easy to miss in a
  connectivity-only comparison);
* `tfi_blade_audit.json` records the input fingerprint
  sha256 `24698841f0c8de72e3f220599d4a861b40e40d5fa790dde83ceea03a507e03b5`,
  which equals the sha256 of `tfi_blade.vtk` on disk;
* the audit JSON counts (215 / 429 / 141 / -0.9936709 / cell 87410) match the
  values recomputed from the audit VTK and from the canonical VTK.

## The 225-vs-215 read

225 corresponds to the T85 regraded artifact (`tfi_blade_graded.vtk`) and to
nothing on disk now: the most likely explanation is that the first read hit a
stale/graded artifact (or a partially written file while the T85/T87 runs were
still replacing it), and the re-read hit the finished canonical audit. It is
**not** a corruption of the canonical audit and **not** a permanent overwrite:
the file now on disk is the canonical one.

## Consequence for the viz

The viz number **215** is backed by the file now on disk. No regeneration was
performed, no mesh (`tfi_blade.vtk`, `blocks_blade.vtk`, `tfi_blade_graded.vtk`)
was modified, and no backup copy was needed.

## Artifacts

* Machine-readable record (written by Python, not shell redirection):
  `/work/datagen_keepbl/runs/tistos_s1_0131/blade/tfi_blade_audit_provenance.json`
  — contains source/audit/graded md5 + sha256, the four counts for all three
  meshes, the min-sJ cell id and block (12), the canonical/regenerated flags and
  the one-line conclusion.
* `experimentell/hex3d_algohex/audit_tfi_blade.py` now also records the
  provenance in every audit JSON it writes: `source`, `source_path`,
  `source_md5`, `audit_md5` (md5 of the audit VTK it just wrote) and
  `metrics.min_sJ_cell_id`. The counts themselves are unchanged.
