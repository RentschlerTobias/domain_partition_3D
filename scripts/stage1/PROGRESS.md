# PROGRESS — Stage-1 2D-Quad-Partition (lebende Datei)

Format je Eintrag: `[ ]/[x] Schritt — Status / letzter Lauf-Output / nächste Aktion`.
Vor Übergabe IMMER eintragen, wo du stehst (welche Datei halb-editiert, Lauf grün/rot).
Siehe `HANDOFF.md` für Kontext.

## Iteration 4 (2026-06-30): Separatrix-Count-Invariante + Pitchwise-Periodik
Env: `/root/venv/bin/python` (hat meshio+torch+torch_geometric; dp-env NICHT).
- [x] Phase 1 — Separatrix-Count-Invariante 3/5 erzwungen + validiert.
      * `clean_separatrix.py::_emanation_dirs`/`_all_local_minima`: nimmt jetzt EXAKT
        `expected_separatrices` (3 für idx+1, 5 für idx−1) tiefste Minima statt fixem thresh.
      * `partition_surface.py::_snap_separatrix_endpoints`: Dedup count-erhaltend
        (Prong-Drop nur wenn beide Endpunkte ≥ Soll bleiben); origin/target_sing auf
        Sep-Dicts gestempelt.
      * `partition_surface.py::_validate_separatrix_counts` + in `validate()` → Report.
      * BEFUND: Counts waren auf dem Hub schon 3/5 (Validator PASS), Dedup war schon sicher.
        Der Euler-Defekt kam NICHT vom Count, sondern von der EXISTENZ innerer Feld-
        Singularitäten (eine −1 IST ein Valenz-5-Innenknoten). User-Hypothese für diesen
        Case falsifiziert, Ursache lokalisiert.
- [x] Phase 2 — Pitchwise-Periodik (θ-Naht). WORKING + verifiziert.
      * `unwrap_surface.py::_detect_periodic_pair`: findet das θ-Seiten-Paar per Translat-
        Test. Hub: 36 conforme Knotenpaare, pitch=0.78540 (=π/4, 8 Blades), reines s-Translat.
      * `dp_adapter.py`: `mesh.periodic_pairs` (master_left,slave_right Knoten-IDs) + `mesh.pitch`.
      * `partition_surface.py`: Flag `PERIODIC=True`/`set_periodic`. `_compute_initial_frame_field_soft`
        weldet Slave-Steifigkeit auf Master + Constraint u_slave=u_master + KEINE Wand-BC auf
        Naht. `_generate_cross_field_seam_aware` nimmt Naht von Re-Imposition aus.
      * VERIFIZIERT: Naht-Stetigkeit exakt 0.0 (u UND frame_field). Periodik-Solve korrekt.
      * WIRKUNG: innere Blade-Umfeld-Singularitäten WEG; an Blade-Tips verschoben (4× −1/4,
        geom. Index-Summe −1 = χ). `detect_singularities`-Warnung ist Fehlalarm (nutzt
        nicht-periodische Euler-Zahl).
- [x] Phase 1d — Tip-Singularitäten nutzen statt löschen.
      * `clean_separatrix.py::_delete_corner_singularities`: nur an AUSSEN-Ecken (corner_type==0)
        löschen, Blade-Tips (==1) verschonen.
      * `_emanate_from_blade_tips`: Tip überspringen wenn dort schon eine Feld-Sing sitzt
        (radius 0.08) → keine Doppelbelegung.
      * LAUF (periodic): 4 Sing behalten, 18 Sep, **5 Quad-Blöcke**, Count-Invariante PASS,
        Element-Qualität GUT (0 invertiert, SJ_min 0.66, Winkel [52,139], Aspect 4.4).
- [ ] **OFFEN — Volle Passage-Tilung / periodischer Block-Schluss.** Aktueller FAIL:
      Euler=7, 5 irreguläre Innenknoten, BOUNDARY-Dist 0.46 — Partition deckt nur die
      Blade-Wrap-Region ab; die Parallelogramm-Außenbereiche (oben-links/unten-rechts)
      werden NICHT in Quads zerlegt. Ursache: Feld ist periodisch, aber `StreamlineGenerator_v2`/
      `StreamlinePostProcessor` behandeln die schrägen Naht-Seiten noch als WAND →
      Separatrizen/Block-Kanten laufen nicht über die Naht weiter, Blöcke schließen nicht
      periodisch. Separatrices-Plot zeigt sauberen O-Grid ums Blade + Wake-Linien zu den
      axialen Wänden, aber keine Tilung der Außenkanäle.
      NÄCHSTE AKTION: Periodik in die Streamline-/Block-Stufe ziehen (Naht-Continuation +
      periodischer Block-Schluss, Monkeypatch im Postprocessor). Siehe Plan-Datei Phase 2c-Folge.

## Iteration 5 (2026-07-01): Kowalski §4.3 Blend + Tip-Cluster-Befund
- [x] Blade-Tip-Ecken abschaltbar (`unwrap_surface.MARK_BLADE_TIP_CORNERS`/`set_blade_tip_corners`).
      TEST: Ecken AUS → schlechter (3 Blöcke, 8 irreg). Ecken AN (default) besser (5 Blöcke).
      Grund: die Ecken dienen als Boundary-Split-Andockpunkte; die Feld-Tip-Singularitäten
      machen die Emission ohnehin. Ecken AN lassen.
- [x] Kowalski §4.3 (Gl. 28) Separatrix-Blend statt Prong-Drop in `_snap_separatrix_endpoints`.
      Zwei zweimal-integrierte sing→sing-Linien (S0→S1 und S1→S0) werden linear verschmolzen
      `γ_b=(1−s)γ1+s·γ2rev` zu EINER genauen Kurve. LAUF: 2 Paare verschmolzen (die kurzen
      Intra-Tip-Verbinder 0↔1, 2↔3). Count-Invariante PASS.
- [BEFUND] Die 2 parallelen Ringe ums Blade sind KEINE Doppelung — es sind Ober-/Unterbogen
      des O-Grids (sep0: sing0(TE)→über-Oberseite→sing3(LE); sep12: sing2(LE)→Unterseite→sing1(TE)).
      Beide nötig. Kowalski Fig 20/§6.2 „thin part" = dasselbe Prinzip.
- [BEFUND] Eigentlicher Tiling-Blocker: an jedem Tip klumpen 3 fast-gleiche Ziele
      (2 Feld-Sing + 1 Blade-Tip-C0-Ecke). Wraps enden ~0.02–0.04 neben der richtigen Sing,
      snappen aber auf die Ecke (target=None) → O-Grid-Ring schließt nicht zu Blöcken.
- [~] Versuch „Snapping bevorzugt Feld-Sing" → BACKFIRE (über-Blend: 5 statt 2 Paare, ein
      Sing verliert Prong got=4, 4 Blöcke). REVERTIERT. Lokale Snap-Heuristiken rippeln
      unvorhersehbar in den Black-Box-`StreamlinePostProcessor`.
- [ ] OFFEN: Tiling-Blocker prinzipiell lösen, NICHT per weiterer Snap-Heuristik. Optionen:
      (a) verstehen wie `StreamlinePostProcessor` Zellen aus Streamlines baut (Quelle lesen),
      dann periodischen Block-Schluss + Tip-Cluster sauber einbauen;
      (b) Tip-Cluster als 1 Graph-Knoten je Tip zusammenfassen (2 Sing + Ecke), Wraps sing→sing
      schließen — aber User will beide Sing behalten;
      (c) Phase 3 Global-Optimum-Feld (evtl. 1 saubere Sing/Tip).
      Aktueller bester Stand: 5 Blöcke, 0 invertiert, SJ_min 0.66, Count PASS; Passage-Außen-
      kanäle + Blade-Ring ungetilt (Euler=7).

## Iteration 6 (2026-07-01): Exaktes Postproc-Paper gefunden (Xiao 2020)
- REFERENZ: `literature/1-s2.0-S0955799720300035-main.pdf` = Xiao, He, Xu, Chen, Wu (2020),
  Eng. Anal. Bound. Elem. 113. Ist die Referenz-Impl von domain_partition. Siehe Memory
  `xiao2020-streamline-simplification`. **§4.2.2 + Algorithm 2 + Fig 8/9 + Gl. 9/10** = das
  exakte Streamline-Postprocessing:
  * Flags f=in/f=out je Streamline an jeder Singularität.
  * Merge-Kriterium Gl. 9 (Winkel): 3-valent 120–240°, 5-valent 144–216°.
  * Fall 1 (beide erreichen) → Blend Gl. 10. Fall 2 (nur eine erreicht) → projizieren, splitten,
    mergen, Stack-Propagation. Zirkuläre Streamlines (Blade-Wrap!) → gleiche-Tri-zweimal
    erkennen, auf Quell-Sing projizieren, schließen.
- BEFUND: aktueller `_snap_separatrix_endpoints` macht NUR Fall 1 (Blend). Unsere Blade-Wraps
  sind Fall 2 / zirkulär → schließen nicht.
- [~] Versuche „prefer-singularity Snapping" + „Winkel-Kriterium Gl. 9" als Patch aufs bestehende
  Nearest-Snapping → BEIDE über-mergen (5 Paare, ein Sing got=4, 4 Blöcke). REVERTIERT.
  Grund: Nearest-Snapping truncatet Wraps an Singularitäten die sie nur passieren. Das Paper
  nutzt Flag-Reach + Projektion, NICHT Nearest. Patchen des alten Snappings ist der falsche Weg.
- [x] `StreamlinePostProcessor`-Kette gelesen: `StreamlineMerging` (=Alg 2) →
  `StreamlineIntersectionSplitter` → `QuadFaceGenerator` (baut 4-Seiter aus Streamline-Schnitten).
- BEFUND (wichtig): `StreamlineMerging` implementiert Alg 2 SCHON:
  * `search_duplicated_streamlines`+`merge_streamlines` = Fall 1 (start_i==end_j & start_j==end_i),
    Blend via `interpolate_streamlines` (Spline, Gl. 10). AKTIV.
  * `find_missed_streamline_endpoints` = Fall 2 (Winkel ~180° ±45°, schneidet Streamline die eine
    Sing passiert (min_dist<0.1) auf sie). **DEFINIERT ABER NIE AUFGERUFEN** (Dead Code in __init__).
  * prepare_streamlines matcht Endpunkte an Terminationsknoten (c0-Ecken + Sing) mit tol=0.01.
- [x] Fall-2-Monkeypatch in `partition_surface.py` (`_streamline_merging_init`): ruft
  find_missed zwischen prepare und merge. Mit MEINEM `_snap` davor: no-op (5 Blöcke, PASS) —
  weil mein Snap die Wraps schon auf die Ecke snappt, StreamlineMerging sieht sie als terminiert.
- TEST ohne mein `_snap` (nur domain_partition Merging+Fall2): 7 Blöcke, 9 irreg, 0 invertiert.
  Mehr Abdeckung, mehr Irreguläre. (mein Count-Validator misst hier nicht — Stamps fehlen.)
- BEFUND: mein `_snap` und domain_partitions StreamlineMerging überlappen/konfligieren. Zwei
  saubere Architektur-Optionen:
  (A) mein `_snap` behalten (macht Boundary-T-Junctions die StreamlineMerging NICHT macht) +
      StreamlineMerging nur als Fall-1-Ergänzung. Aktuell: 5 Blöcke PASS.
  (B) mein `_snap` fallenlassen, ganz auf domain_partition StreamlineMerging (Fall1+2) setzen +
      QuadFaceGenerator/Periodik dort fixen. Aktuell: 7 Blöcke, 9 irreg.
- [ ] NÄCHSTE AKTION: `QuadFaceGenerator` (streamlines_to_quad_faces.py) +
  `StreamlineIntersectionSplitter` lesen — verstehen wie/warum 4-Seiter NICHT über die ganze
  Passage schließen (Blade-Ring + Naht). Das ist der eigentliche Tiling-Blocker, nicht mehr das
  Separatrix-Postproc. Danach Architektur A oder B wählen + Periodik im Face-Builder.

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
