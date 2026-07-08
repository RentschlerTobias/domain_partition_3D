# HANDOFF: Periodisches Streamline-Merging (konvergierende Helices)

## Kontext / Stand (2026-07-07, Branch `wip-backup`, letzter Commit `38d446d`)

Stage-1 Quad-Domain-Partition der T1_9-Hub-Passage. Pipeline: STL → Zylinder-Unwrap
(`unwrap_surface.py`) → periodisches Cross-Field (Seam-Weld, `partition_surface.py`
Monkeypatches) → Separatrix-Emanation (`clean_separatrix.py`) → RK-Integration
(`StreamlineGenerator_v2` aus `../domain_partition`, gepatcht) → Snap/Dedup
(`_snap_separatrix_endpoints`) → periodisches Tiling → `StreamlinePostProcessor` →
Quad-Blöcke. KEINE Edits an `../domain_partition` — alles Monkeypatch in
`scripts/stage1/partition_surface.py`.

**Was funktioniert (committed `38d446d`):**
- Wrap-Integration: Streamlines überqueren die Theta-Naht beliebig oft
  (`_wrap_shifts`, `MAX_WRAPS=4`; gepatchte `find_containing_face`,
  `get_best_cross_vector`, `check_termination_criteria`).
- Snap-Ziele/Boundaries werden um k Pitches getilt; Wand-Treffer auf zentrale
  Polyline zurückgemappt; Hausdorff-Dedup mod pitch.
- Akute Außenecken (0,0)/(1,1) emittieren nichts (physisch identisch mit den
  obtusen Naht-Partnern (0.593,0)/(0.407,1)).
- Ergebnis `38d446d`: 39 Blöcke, 1 irregulär (TE-Cluster), 0 invertiert.

**Offenes Problem:** Die 4-5 LE/TE-Prongs Richtung Inlet/Outlet erreichen NIE eine
nicht-periodische Boundary: sie laufen helikal ums Rad und **konvergieren gegen einen
geschlossenen Limit-Orbit** (Windungsabstand schrumpft pro Umlauf, sichtbar in
`output/T1_9/hub_stage1/streamlines_central.png` von `38d446d`). Bei `MAX_WRAPS`-Cap
enden sie mid-air.

**Fehlversuch (uncommitted, in `_close_helical_streamlines`):** Selbst-Treffer
mod pitch mit `tol=0.035` → Ring SOFORT beim ersten Wiedereintritt geschlossen.
User-Urteil: „es war viel besser vorher … jetzt wird sie direkt als Ring gemacht,
was die gesamte Partition bricht." Der Ring entsteht zu früh (Spirale noch nicht
konvergiert), Prong+Ring liegen ~0.03-0.05 auseinander → Klon-Ringe/Slivers.

## AUFGABE: Konvergenz-basiertes periodisches Merging

User-Spezifikation (wörtlich): „wir benötigen etwas, das die konvergierenden
Streamlines merged, aber erst wenn sie relativ nahe an der konvergierenden
Stromlinie ist und evtl. schon ein paar Umläufe gemacht hat, aber nicht direkt
am Anfang."

Konkret in `_close_helical_streamlines(mesh)` (partition_surface.py):
1. **Merge-Kriterium verschärfen** (statt `tol=0.035`, `align=0.7`):
   - `tol_merge ≈ 0.015` (eng: erst mergen wenn die Windung der vorherigen
     wirklich nahe ist),
   - `align ≥ 0.9`,
   - **Mindest-Umläufe**: Selbst-Treffer erst zulassen, wenn die Kurve vor dem
     Treffpunkt schon ≥ 2 Pitches netto in s zurückgelegt hat
     (`|s_j − s_0| ≥ 2·pitch_norm`) — nicht direkt am Anfang mergen.
2. **`MAX_WRAPS` auf 6 erhöhen** (Spirale braucht Raum zum Konvergieren, sonst
   Cap vor Konvergenz).
3. **Semantik beibehalten** (Prong+Ring-Modell, bereits implementiert):
   - Prong = `P[:i+1]` (enthält die früheren Spiral-Windungen! Die bleiben
     sichtbar/erhalten — das war das „vorher besser"),
   - Ring = `P[i:j+1]`, Drift linear verteilt so dass Ende == Start + k·pitch
     EXAKT; Ring bleibt am Orbit (KEINE Projektion auf den Quell-Sing — der
     liegt 0.24 entfernt, das hat die Partition zerstört),
   - Ring als BOUNDARY-Streamline einfügen (vor der Sep-Sektion), NICHT als
     Separatrix (keine Valenz-Störung; Snap-Pass-2 kann ihn splitten, Prongs
     docken als T-Junction wie am Blade),
   - Ring-Dedup: ein Ring pro Orbit (set-Hausdorff mod pitch, `ring_dup_tol`);
     nach Konvergenz sollten Klon-Ringe < 0.02 auseinander liegen → Dedup
     greift dann (beim Fehlversuch: 5 Ringe statt ~2, weil zu früh geschlossen).
4. Aufruf-Reihenfolge in `partition()`: `_drop_degenerate_corner_seps` →
   `_close_helical_streamlines` → `_snap_separatrix_endpoints` (Prong-Dock
   braucht den Snap).

## Verifikation / Deliverables

`plot_central_streamlines.py` so umbauen, dass EIN Lauf ZWEI PNGs schreibt
(Integration ist teuer, einmal integrieren, Streamline-Kopien ziehen):
1. `streamlines_central_raw.png` — ALLE Streamlines OHNE Merge (Showcase:
   volle Spiralen, gefaltet mod pitch, Naht-Kreuzungen interpoliert; Snap ohne
   Closing anwenden damit Doppel-Kurven dedupliziert sind).
2. `streamlines_central_merged.png` — MIT Merge (Closing + Snap).

Gate: Merged-Plot zeigt Spiral-Prongs die nach ≥2 Umläufen in einen sauberen
geschlossenen Ring münden; `plot_c0corner_seps.py` liefert wieder ≥ 36 Blöcke,
≤ 1 irregulär, 0 invertiert — Partition darf NICHT brechen (0 Blöcke = der
Fehlversuch). Ring-Anzahl im `[helix]`-Print: erwartet ~2 (ein Orbit oben, einer
unten), nicht 5.

## Dateien
- `scripts/stage1/partition_surface.py` — `_close_helical_streamlines` (Z. ~373),
  `MAX_WRAPS` (Z. ~104), `_wrap_shifts`, `_snap_separatrix_endpoints`.
- `scripts/stage1/plot_central_streamlines.py` — Fold-Plot (raw/merged Umbau).
- `scripts/stage1/plot_c0corner_seps.py` — Block-Gate.
- Python: `/root/venv/bin/python` (läuft 2-4 min; im Hintergrund starten).

## Nicht anfassen
- `../domain_partition` (nur Monkeypatch).
- CLAUDE.md nie committen/erwähnen. Keine Co-Author-Trailer in Commits.
