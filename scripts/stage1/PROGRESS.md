# PROGRESS — Stage-1 2D-Quad-Partition (lebende Datei)

Format je Eintrag: `[ ]/[x] Schritt — Status / letzter Lauf-Output / nächste Aktion`.
Vor Übergabe IMMER eintragen, wo du stehst (welche Datei halb-editiert, Lauf grün/rot).
Siehe `HANDOFF.md` für Kontext.

- [x] 0. `HANDOFF.md` + `PROGRESS.md` anlegen. — DONE.
- [x] 0b. `git rebase --abort` (pausierter Rebase). — DONE, Commit c997168 erhalten.
- [x] 1. `unwrap_surface.py`: `corner_type` (0=Außen,1=Tip,-1=keine) + `blade_loops` ausgeben.
      DONE. Lauf: 4 Außen-, 2 Blade-Tip-Ecken.
- [x] 2. `dp_adapter.py`: `mesh.corner_type` + `mesh.blade_loops` (normalisiert). DONE.
- [x] 3. `clean_separatrix.py`: `_emanate_from_blade_tips` ersetzt `_emanate_from_corners`.
      Außen: keine Emanation; LE/TE: 3 Cross-Arme via `_cross_dirs`, gefiltert mit
      `matplotlib.path.Path(blade_loop)`. DONE.
- [x] 4. `partition_surface.py`: `validate()` mit `QuadPartitionValidator(strict=True)` nach
      Postproc; Report nach `output/T1_9/hub_stage1/validation_report.txt`. DONE.
- [x] 5a. BC-Verfeinerung: **Mittelung im Representative-(4θ)-Raum** statt Tangenten-Raum
      (`r0=map(angle0); r1=map(angle1); ref=atan2(sin r0+sin r1, cos r0+cos r1)`).
      Entfernt die Eck-Streu-Singularitäten (roh 10→6). DONE.
      Hinweis: Einzel-Tangenten-Variante brach Blockbildung (0 Blöcke) → verworfen.
- [x] 5b. `bc_weight`-Knopf (weiche Penalty-BC, `set_bc_weight`/`BC_WEIGHT`,
      `_compute_initial_frame_field_soft`). IMPLEMENTIERT — aber **HILFT NICHT**.
      Sweep w∈{None,300,100,30,10,3,1}: `nsing=4` KONSTANT bei allen Gewichten;
      w≤30 bricht Blockbildung (0 Blöcke), w=1 Solver-Fehler. Befund: Singularitäten-
      ANZAHL ist topologisch fix durch die Wand-Alignment-BC (Poincaré-Hopf / Rand-
      Holonomie), NICHT durch Innen-Energiegewicht steuerbar. Default bleibt `None` (hart).
      Knopf bleibt drin (dokumentiert, harmlos).
- [x] 6. Annotierter Block-Plot (`_block_annotations`): irreguläre Innenknoten (rot),
      invertierte (rot gefüllt), High-Aspect>10 (orange). DONE.
- [ ] 7. Iterieren bis `is_valid()` & `passes_soft_thresholds()`. **OFFEN.** Aktueller Stand:
      19 Blöcke, 0 invertiert, SJ_min=0.79 SJ_mean=0.95, Winkel [52,121]. NOCH FAIL:
      Euler=−1 (statt 0), **1 irregulärer Innenknoten** (Valenz≠4, ~(0.46,0.64) über Blade),
      **2 Blöcke Aspect≈11**(>10), singularity_efficiency 0 (erwartet 0 Interior-Sing., ist 1).
      BOUNDARY-Check (197 Knoten) = Artefakt: vergleicht feine Tri-Randknoten gegen 36 grobe
      Block-Ecken (Knoten-Hausdorff), für grobe Blöcke prinzipbedingt — NICHT actionable.

## Nächste konkrete Aktion (für übernehmenden Agenten)
Schritt 7: die **4 verbleibenden Feld-Singularitäten am Blade** erzeugen Euler=−1 + 1 irregulären
Innenknoten. `bc_weight` (5b) erledigt das NICHT (Anzahl topologisch fix, s.o.). Verbleibende
echte Hebel:
1. **Pitchwise-Periodizität**: die θ-Seitenränder (s=links/rechts) werden aktuell als WÄNDE
   behandelt (MVP-Annahme, Plan). Echte Passage = periodisch. Wand-Alignment auf künstlichen
   θ-Wänden treibt vermutlich die Holonomie/Singularitäten. Periodische BC = physikalisch
   korrekt, wahrscheinlichste Ursache der Topologie-Defekte. War für Stufe 2 vorgesehen, ist
   aber jetzt der naheliegende Hebel.
2. **Singularitäten-Paar-Annihilation** (Kowalski): IMPLEMENTIERT als `_annihilate_pairs`
   (mutual-nearest opposite-index, `set_annihilate_pairs`/`ANNIHILATE_PAIRS`, default OFF).
   GETESTET: entfernt das einzige +1/−1-Paar (4→2 Sing.), Blöcke 19→11, ABER **bricht den
   Blade-Wrap** (großer unblockierter Bereich rechts, irregulärer Knoten schwebt im Blade-Loch)
   und Euler kippt −1→+1. Grund: das einzige annihilierbare Paar (+1@(0.61,0.49) / −1@(0.70,0.55))
   ist load-bearing für die rechte Blade-Umwicklung. Reine Löschung verliert die verbindende
   Separatrix. Richtige Kowalski-Variante müsste die Bridge-Streamline als reguläre Block-Kante
   ERHALTEN (noch nicht implementiert). Default OFF gelassen.
3. **Akzeptieren**: 19-Block-Ergebnis ist die beste Partition (0 invertiert, SJ 0.79/0.95, nur
   1 irregulärer Knoten + 2 leicht gestreckte Blöcke, wickelt Blade voll) → weiter zu
   3D-Mapping/Hexa-Extrusion.

## Synthese (beide Interior-Hebel topologisch gedeckelt)
`bc_weight` (Energie) UND reine Paar-Annihilation scheitern an derselben Ursache: bei
WAND-berandeter Passage ist die Singularitäten-Bilanz durch die Rand-Holonomie fixiert. Echte
Lösung = **Pitchwise-Periodizität** (ändert die Holonomie, Hebel 1) ODER 19-Block-Stand
akzeptieren. Alternativ: Annihilation-mit-Bridge-Erhalt implementieren.

## Iteration 3 (User-Feedback: 3 Separatrix-Fixes, bleiben in 2D)
- [x] Fix 3 — Blade-Tip-Separatrizen liefen ins Blade. Neu `_emanate_from_blade_tips`:
      primärer Arm = Cross-Arm max-aligned mit Fluid-Richtung; +2 EXAKT orthogonal (±90°);
      anti-primärer (ins-Blade) Arm gedroppt. (`clean_separatrix.py`)
- [x] Fix 1 — Snap-Dedup für 3/5-Valenz. `_snap_separatrix_endpoints` (`partition_surface.py`)
      umgebaut: trackt origin/target-Singularität je Separatrix; wenn Sep A→B snappt, wird am B
      der Prong mit kleinster Winkelabweichung zur Verbindungslinie (zeigt zurück zu A) gelöscht.
      Boundary-Splits erst nach Dedup (keine verwaisten T-Junctions).
      LAUF: deduped 2 Prongs, Sep 24→22, 17 Blöcke. ABER: Euler −1→+1, irreguläre Knoten 1→3.
      Grund vermutlich: Feld hat UNBALANCIERTE Indizes (1×+1, mehrere −1) → keine sauberen
      3/5-Paare möglich; überschüssige −1 hängen an Boundary, nicht aneinander → Dedup-Eingriff
      stört Tiling (gleiche topologische Wurzel wie bc_weight/Annihilation). VALENZ-CHECK offen
      (Bash-Klassifizierer war kurz down).
- [~] Fix 2 — Knicke. Befund: RK-Integrator `get_best_cross_vector` wählt SCHON max-dot
      (= min Winkelabweichung). Sichtbare Knicke v.a. von (a) Blade-Arme ins Blade (→ Fix 3)
      und (b) Snap-Geraden-Sprüngen. Nach Fix 1/3 neu bewerten.

## Stand vor dieser Iteration (Referenz)
17 Quad-Blöcke, Außenecken MIT Streu-Innen-Separatrix (Artefakt), LE/TE nur 1 Streamline,
keine Validierung. → jetzt behoben außer Interior-Singularitäten (Schritt 7).
