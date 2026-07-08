# HANDOFF: Periodik-Vergleich — Wrap-Ansatz vs. Wand + T-Mesh + TFI

Selbständig ausführbar, kein weiterer Kontext nötig. Ziel-Entscheidung des Users:
**"Beide vergleichen"** — den Standard-Ansatz (Streamlines enden an der periodischen
Boundary; T-Mesh-Blöcke; TFI; Master-Slave-Naht) bauen und gegen den bereits
implementierten Wrap-Ansatz (Streamlines laufen über die Naht weiter) vergleichen.

---

## 1. Mission

T1_9-Turbinen-Hub-Passage (Zylinderfläche, 1 von 8 Blade-Passagen) → 2D-Unwrap →
Cross-Field → Separatrizen → **grobe Quad-Block-Partition**. Die Passage ist
pitchweise periodisch (Theta-Nähte). Zwei Wege, die Periodik zu behandeln:

- **Ansatz W (Wrap, IMPLEMENTIERT, uncommitted):** Streamlines integrieren über
  die Naht hinweg (Universal-Cover), Helices werden auf Limit-Orbits gemergt,
  Orbit-Ringe + Dock-Quer-Kurven. Stand: 69 Blöcke, 1 irregulär, 0 invertiert.
- **Ansatz T (Wand + T-Mesh + TFI, ZU BAUEN):** Streamlines enden an der Naht.
  Naht-Junctions werden mod pitch auf die Gegen-Naht gespiegelt (Master-Slave,
  wie DLR Sauer/Morsbach 2023 §2.7: nur eine Naht-Seite parametrisiert, andere
  Seite = fixe Transformation). Block-Extraktion muss T-Junctions tolerieren
  (Quad mit unterteilter Kante statt Pentagon-Drop). TFI füllt jeden Block.

Dann beide mit identischen Metriken vergleichen (Abschnitt 7).

## 2. Umgebung

- Repo: `/root/repos/block_structured_meshing`, Branch `wip-backup`.
- Solver-Repo (READ-ONLY!): `/root/repos/domain_partition`.
- Python: `/root/venv/bin/python` (torch, torch_geometric, meshio, matplotlib,
  networkx, scipy, numba). Läufe dauern 2–10 min → **immer im Hintergrund
  starten** (`run_in_background`), Ausgabe tailen.
- Input-STL: `/root/repos/block_structured_meshing/T1_9_hub_raw.stl`
  (exakter Zylinder r=0.5, z∈[0,2.5]).
- Outputs: `output/T1_9/hub_stage1/`.

## 3. Harte Constraints

1. **KEINE Edits in `/root/repos/domain_partition`** — alle Solver-Änderungen als
   Monkeypatch in `scripts/stage1/partition_surface.py` (Muster dort vorhanden).
   Neue eigene Module in `scripts/stage1/` sind erlaubt.
2. CLAUDE.md niemals committen oder in Commits/PRs erwähnen.
3. Keine `Co-Authored-By`/KI-Attribution in Commit-Messages.
4. Vor Phasen-Wechseln committen (User pusht selbst oder `git push` versuchen;
   wenn Permission-Prompt scheitert, User bitten).
5. User bestätigt Zwischenstände anhand von Plots — nach jedem Meilenstein Plot
   erzeugen und Pfad nennen, TERSE berichten.

## 4. Pipeline-Karte (wo ist was)

Aufruf-Kette (`scripts/stage1/partition_surface.py::partition()`):
```
build_dp_data (dp_adapter.py)         # Unwrap + Normalisierung auf [0,1]^2
  → unwrap_surface.unwrap             # (s,t)=(r·θ,z); periodic_pairs; pitch
FrameField (../domain_partition/tools/frame_field.py, stark gepatcht)
  → periodischer Seam-Weld            # _compute_initial_frame_field_soft
detect_singularities                  # 4 Sings, alle idx=−1 (2 pro Blade-Tip!)
StreamlineGenerator_v2 (gepatcht)     # RK-Integration der Separatrizen
  Patches: _robust_find_containing_face (+_wrap_shifts, MAX_WRAPS=8),
           _periodic_get_best_cross_vector, _periodic_check_termination
CleanSeparatrixGenerator (clean_separatrix.py)  # Emanation; akute Ecken leer
_drop_degenerate_corner_seps
_close_helical_streamlines            # Ansatz W: Helix→Prong+Orbit-Ring
_emit_dock_crossings                  # Ansatz W: Dock-T → X (orthogonale Familie)
_snap_separatrix_endpoints            # Clipping/Dedup/Boundary-Split (T-Junctions)
_tile_streamlines_periodic            # Ansatz W: ±k-Pitch-Kopien (Dedupe)
StreamlinePostProcessor (../domain_partition/tools/streamline_post_processor.py)
  = StreamlineMerging → StreamlineIntersectionSplitter → QuadFaceGenerator
_extract_central_blocks               # Ansatz W: Faces mit Zentroid im Zentral-Pitch
```

Flags in `partition_surface.py`:
- `PERIODIC` / `set_periodic(f)` — Feld-Seam-Weld (Feld exakt periodisch).
- `TILE_PERIODIC` / `set_tile_periodic(f)` — GESAMTE Block-Stufen-Periodik
  (Wrap-Integration, Snap-Tiling, Helix-Merge, Docks, Streamline-Tiling,
  Zentral-Extraktion). Alle Funktionen prüfen dieses Flag selbst.
- `set_periodicity(f)` — Master-Schalter (beide zusammen).
- Für **Ansatz T**: `set_periodic(True)` + `set_tile_periodic(False)` — das Feld
  bleibt periodisch (Separatrizen-Geometrie spiegelt sich dann von selbst mod
  pitch), aber die Block-Stufe behandelt die Naht als Wand.

Kritischer Upstream-Code (READ-ONLY, Verhalten kennen):
- `../domain_partition/tools/streamlines_to_quad_faces.py::QuadFaceGenerator`
  - `build_connectivity` (Z.36): Graph-Knoten = NUR Streamline-Endpunkte
    (`s[0]`, `s[-1]`), Knoten-Merge-Toleranz `tol=1e-2`.
  - `get_quad_faces` (Z.74): planare Einbettung (`nx.check_planarity`),
    Regionen via `traverse_face`, **Filter `len(face)==4` droppt alles andere
    STILL** (Z.101) — Ursache aller "Löcher". Nicht-planar → 0 Faces.
- `StreamlineMerging.prepare_streamlines`: matcht Streamline-Enden gegen
  Terminations-Knoten (Sings + C0-Ecken), tol 1e-2; unmatched läuft unverändert
  durch.
- `StreamlineIntersectionSplitter`: splittet Kurven an transversalen Schnitten
  (Splines + minimize). Koinzidente/fast-parallele Kurven = pathologisch
  (deshalb überall Dedupe!). Kurven mit gemeinsamem Endpunkt werden übersprungen.

## 5. Aktueller Zustand (WICHTIG)

- Letzter Commit: `38d446d` (Wrap-Integration + Showcase-Plot). Davor `5dc596f`
  (Tiling), `b787876` (Corner-Subdivider), `9868bbf` (Xiao+Clipping).
- **UNCOMMITTED im Working Tree** (= Ansatz W komplett):
  `partition_surface.py` (Konvergenz-Merge `_close_helical_streamlines` mit
  tol=0.015/align=0.9/min_wraps=2, `_emit_dock_crossings`, `MAX_WRAPS=8`,
  `set_periodicity`), `plot_central_streamlines.py` (raw+merged Zwei-Plot-Lauf),
  Output-PNGs. **ERSTER SCHRITT: das als Ansatz-W-Baseline committen** (Message
  sachlich, Ergebniszeile 69/1/0, kein Co-Author).
- Ansatz-W-Reststand (bekannte offene Punkte, NICHT Blocker für den Vergleich):
  `[helix] closed 5 prongs onto 3 rings` — 1 Klon-Ring zu viel (`_ring_close`
  Dedup-Toleranz 0.06 knapp); wenige Rest-Lücken am unteren Ring-Band/unter TE.

## 6. Ansatz T — Bauplan (Schritt für Schritt)

Neues Skript/Modul, Vorschlag: `scripts/stage1/tmesh_partition.py` (+ Reuse von
`partition_surface`-Helpern via Import). Keine Verhaltensänderung am Ansatz-W-Pfad!

### 6.1 Lauf-Konfiguration
```python
import partition_surface as ps
ps.set_periodic(True)        # Feld periodisch (Seam-Weld)
ps.set_tile_periodic(False)  # Block-Stufe: Naht = Wand
```
Pipeline bis inkl. `_snap_separatrix_endpoints` wie in `partition()` (Code-Pfad
existiert; bei TILE_PERIODIC=False sind Wrap/Merge/Docks/Tiling automatisch inert).
Separatrizen, die die Naht treffen würden, enden dort (Integrator-Exit) und werden
vom Snap auf die Naht-Wand gesnappt + splitten sie (T-Junction) — existierendes
Verhalten.

### 6.2 Naht-Symmetrisierung (Master-Slave, DLR-Prinzip)
Nach dem Snap, vor dem PostProcessor:
- Naht-Polylines identifizieren: die Boundary-Streamlines, deren Endpunkte auf
  den Theta-Wänden liegen. Geometrie: linke Naht = Strecke (0,0)→(0.407,1),
  rechte = (0.593,0)→(1,1) in normalisierten Koordinaten;
  `pitch_norm = mesh.pitch_norm ≈ 0.5932` (reine s-Translation).
- Alle Junction-Punkte (Segment-Grenzen der gesplitteten Wand-Polylines) beider
  Nähte sammeln; Vereinigungsmenge mod pitch bilden; fehlende Punkte auf der
  jeweils anderen Naht einfügen (Wand-Polyline dort zusätzlich splitten).
  Ergebnis: rechte Naht-Segmentierung == linke + [pitch_norm, 0] exakt
  (Punktkoordinaten von der Master-Seite kopieren + verschieben, nicht neu
  projizieren — sonst Toleranz-Drift).
- Gespiegelte Junctions ohne eigene Innen-Kurve sind bewusst hängende Knoten
  (T-Junctions) — 6.3 behandelt sie.

### 6.3 T-Mesh-Block-Extraktion (Ersatz für QuadFaceGenerator-Filter)
Eigenes Modul `scripts/stage1/tmesh_faces.py`:
- Wie `QuadFaceGenerator.build_connectivity`/`get_quad_faces` den planaren
  Regionen-Graph bauen (Code-Logik übernehmen/kopieren, NICHT domain_partition
  editieren), aber statt `len(face)==4`-Filter:
  - pro Region die Innenwinkel an den Knoten messen (aus den anliegenden
    Polyline-Tangenten);
  - Knoten mit Winkel ≈ 180° (z.B. |180°−α| < 25°) = FLACH (hängender
    T-Knoten auf einer Kante), alle anderen = ECHTE Ecken;
  - Region akzeptieren wenn EXAKT 4 echte Ecken (beliebig viele Flachknoten);
    ihre 4 logischen Seiten = Kantenzüge zwischen echten Ecken.
  - Ausgabe pro Block: 4 Seiten als Polylines (inkl. Flachknoten), Ecken-Ids.
  - Außenregion/Blade-Loch wie gehabt verwerfen (größte Region / Orientierung).
- Validierung wie `_block_annotations` (irregulär/invertiert) sinngemäß auf
  echte Ecken bezogen.

### 6.4 TFI
Pro Block lineare transfinite Interpolation (Coons-Patch) aus den 4 Seiten:
`X(u,v) = (1−v)·S(u) + v·N(u) + (1−u)·W(v) + u·E(v) − Bilinear(Ecken)`.
- Seiten nach Bogenlänge parametrisieren; Auflösung z.B. 8×8 pro Block.
- **Master-Slave-Diskretisierung**: Naht-Seiten der rechten Naht bekommen exakt
  die Knotenverteilung der linken Naht + [pitch_norm, 0] (nicht unabhängig
  sampeln). Geteilte Innen-Kanten: beide anliegenden Blöcke nutzen dieselbe
  Seiten-Diskretisierung (Kante einmal sampeln, von beiden Blöcken referenzieren).
- Ausgabe: 2D-Mesh-Plot aller Blöcke (Zellkanten dünn), Naht-Knoten beider
  Seiten hervorheben.

### 6.5 Gate Ansatz T
- Jede zentrale Region hat exakt 4 echte Ecken (Report: Anzahl Regionen mit ≠4
  echten Ecken = 0), keine gedroppten Löcher.
- Naht-Konformität: Knotensatz rechte Naht == links + pitch (max. Abweichung
  < 1e-6 nach Konstruktion, prüfen + printen).
- TFI-Zellen: 0 invertierte Zellen (Shoelace pro Zelle).

## 7. Vergleich (Deliverable)

Skript `scripts/stage1/compare_periodic_approaches.py`, EIN Lauf, schreibt nach
`output/T1_9/hub_stage1/periodic_compare/`:
- `wrap_blocks.png` (Ansatz W: vorhandene `partition()` mit Periodik AN),
- `tmesh_blocks.png` + `tmesh_tfi.png` (Ansatz T),
- `compare.txt` mit identischen Metriken beider Ansätze:
  Blockzahl, irreguläre Innenknoten, invertierte Blöcke/Zellen, Anzahl
  gedroppter Nicht-Quad-Regionen im Zentral-Pitch (Löcher), Naht-Konformität
  (Korner-Sets links vs rechts−pitch: Anzahl + max Positions-Mismatch),
  Laufzeit.
- Naht-Konformität-Snippet (bewährt):
```python
used = np.unique(bf)
def dl(p): return abs(p[0] - 0.407*p[1])          # linke Naht
def dr(p): return abs(p[0] - (0.593 + 0.407*p[1]))  # rechte Naht
L = sorted([tuple(np.round(bx[i],3)) for i in used if dl(bx[i]) < 0.02], key=lambda p: p[1])
R = sorted([tuple(np.round(bx[i]-[pitch_norm,0],3)) for i in used if dr(bx[i]) < 0.02], key=lambda p: p[1])
```
- Danach: committen, User die zwei Block-Plots + compare.txt zeigen, WARTEN.

## 8. Bekannte Fallstricke (teuer erkauft — nicht wiederholen)

1. **Helix nicht früh schließen**: Selbst-Treffer-Toleranz 0.035 ohne
   Mindest-Umläufe → Ring zu früh → Partition bricht komplett (0 Blöcke).
   Konvergenz-Kriterium (tol 0.015, align 0.9, ≥2 Pitches) ist der Fix.
2. **Ring NIE auf den Quell-Sing projizieren** — Orbit liegt ~0.24 vom Sing;
   Verschieben deplatziert den Ring → Splitter-Explosion → nicht-planar.
3. **Koinzidente Kurven töten den Splitter** (Splines + minimize auf
   fast-identischen Kurven): überall geometrisch dedupen (set-Hausdorff mod
   pitch; Kurven mit gemeinsamen Endpunkten überspringt der Splitter selbst).
4. **`QuadFaceGenerator` droppt still**: jede Region ≠4 Knoten verschwindet
   (Löcher), nicht-planarer Graph → 0 Faces. 2-Punkt-Stub = Valenz-1-Knoten =
   Schlitz.
5. **Blade-Tip-Cluster**: jeder Tip trägt ZWEI idx=−1-Sings (~0.03–0.05
   auseinander, topologisch zwingend: Summe −1 fürs Blade-Loch im periodischen
   Band). Klon-Prongs beider Sings sind fast identisch → Sliver-Gefahr; der
   verbleibende 1 irreguläre Knoten sitzt am TE-Cluster. Nicht "wegoptimieren".
6. **dtype**: `pitch_norm` als python float halten (`float(...)`), sonst
   float64-Tensoren im torch-Integrator (dot-dtype-Crash).
7. Diagnose-Läufe brauchen die `partition_surface`-Patches: **immer
   `import partition_surface as ps` VOR FrameField-Nutzung**, sonst falsches
   Feld (12 statt 4 Singularitäten).
8. Plots im Hintergrund laufen lassen; `[helix]`/`[dock]`/`[tile]`/`[snap]`-
   Prints sind die schnellste Zustandskontrolle.

## 9. Referenzen

- DLR: Sauer & Morsbach 2023, IJNME 10.1002/nme.7308 — §2.7 periodische
  Master-Slave-Kontrollpunkte; Blöcke enden an der Naht; TFI/algebraisch als
  Initialisierung. (PDF im Netz frei.)
- Viertel et al., IMR 2019 (arXiv 1905.09097) — Limit-Zyklen: Tracing stoppen,
  T-Junctions akzeptieren, danach vereinfachen; manche T-Junctions ohne
  Zusatz-Singularitäten unvermeidbar.
- Myles/Pietroni/Zorin, SIGGRAPH 2014 — robustes Feld-Tracing, T-Mesh-Layout.
- Xiao 2020 (literature/1-s2.0-S0955799720300035-main.pdf) — Streamline-
  Simplification (Case 1/2, Zirkulär-Fall), Basis der Snap/Merge-Logik.
- Kowalski 2015 (../domain_partition/literature/quad_domain_partition.pdf) —
  3/5-Separatrix-Invariante, Vier-Seiter-Garantie.
