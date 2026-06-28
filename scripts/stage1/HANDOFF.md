# HANDOFF — Stage-1 2D-Quad-Partition (abgewickelter Hub)

Diese Datei + `PROGRESS.md` sind der **Übergabepunkt**. Jeder Agent ohne Vorkontext liest
zuerst beide, schaut dann in `PROGRESS.md` welcher Schritt offen ist, und macht dort weiter.
Bei jedem abgeschlossenen Teilschritt `PROGRESS.md` aktualisieren.

## Ziel dieser Iteration
Die 2D-Quad-Partition des abgewickelten Hub-Zylinders **qualitativ gut** machen, mit
objektiver Validierung. Konkret:
1. **Außenrand-Ecken**: emittieren KEINE Innen-Separatrix. Nur die 2 inzidenten Boundary-Edges
   bilden die Block-Ecke (Valenz 2). Wichtig: die Außenecken sind im abgewickelten 2D **schräg,
   nicht 90°** (nur in 3D 90°) — der schräge Bisektor-Dirichlet-BC erzeugt sonst eine
   Streu-Singularität direkt neben jeder Ecke.
2. **Leading/Trailing Edge (Blade-Tips)**: Knoten-Valenz **5** = 2 Runner-Profil-Splines
   (Blade-Boundary-Edges, laufen am Tip zusammen) + **3 Streamline-Separatrizen** ins Fluid,
   entlang dem lokalen Cross-Vektor, **keine** ragt ins Blade-Profil. (Aktuell nur 1 Streamline.)
3. **Validierung**: `QuadPartitionValidator` aus domain_partition verdrahten (Euler, Innenknoten
   Valenz 4, kein invertierter/entarteter Quad, Rand-Abdeckung, Singularity-Efficiency).
4. **Energie/BC-Tuning**: Frame-Field-BC verfeinern (Außenecke an 1 Wand-Tangente statt Bisektor;
   Blade-Tip an Blade-Tangente) + optionaler `bc_weight`-Knopf (weiche Dirichlet-BC).

## Pipeline-Reihenfolge (Validierung erst am Ende!)
```
unwrap_surface.unwrap        # 3D-Zylinder -> 2D (s,t), Loops, Eckentypen
  -> dp_adapter.build_dp_data # torch_geometric Data (x, faces, edges, streamlines, corner_type, blade_loops)
  -> FrameField               # harmonisches Cross-Field, harte Wand-Dirichlet-BC
  -> detect_singularities     # Poincaré-Index pro Face
  -> CleanSeparatrixGenerator # Emanation: Singularitäten + Blade-Tips (NICHT Außenecken)
  -> StreamlineGenerator_v2   # RK/Heun-Integration der Separatrizen
  -> _snap_separatrix_endpoints # Kowalski: Enden auf Singularität/Ecke/Rand snappen + Rand splitten
  -> StreamlinePostProcessor  # Merging + IntersectionSplitter + QuadFaceGenerator => block_mesh
  -> QuadPartitionValidator   # erst HIER validieren (auf block_mesh, nicht rohe Streamlines)
```
Postproc liefert erst die Quad-Faces; nur darauf ist Validierung definiert. Postproc-Fehler /
keine Faces = Validierungs-Fail (melden, nicht crashen).

## Datei-Map (`scripts/stage1/`)
- `unwrap_surface.py` — Zylinder-Abwicklung; gibt `corner_type` (0=Außen, 1=Blade-Tip, -1=keine)
  + `blade_loops` (innere Loop-Polygone in (s,t)) aus.
- `dp_adapter.py` — baut `Data`; hängt `mesh.corner_type` + `mesh.blade_loops` (normalisiert) an.
- `clean_separatrix.py` — `CleanSeparatrixGenerator`: Singularitäts-Emanation + typabhängige
  Ecken-Emanation (Außen: keine; Blade-Tip: 3 Cross-Arme, Blade-Polygon-gefiltert).
- `partition_surface.py` — Orchestrator + alle Monkeypatches; `QuadPartitionValidator`-Wiring;
  BC-Verfeinerung; annotierter Plot. Entry-Point.

## Run-Befehl & Env
```bash
MESH_SCRATCH=/tmp/dp_scratch /root/venv/bin/python \
  /root/repos/block_structured_meshing/scripts/stage1/partition_surface.py
```
Env `/root/venv` (torch CPU, torch_geometric, sklearn, numba, meshio). domain_partition liegt
unter `/root/repos/domain_partition` (per `sys.path.insert` eingebunden, NICHT editieren).
Artefakte nach `output/T1_9/hub_stage1/`.

## Monkeypatches (in `partition_surface.py`, kein domain_partition-Edit)
- `tools.streamline_generator_v2.SeparatrixGenerator` → `CleanSeparatrixGenerator`
  (domain_partition v1/v2-Separatrix-Gen ist buggy/ungetestet auf echten Singularitäten).
- `StreamlineGenerator_v2.find_containing_face` → `_robust_find_containing_face`
  (k-nächste-Knoten + Brute-Fallback statt 1-Knoten → keine mid-mesh-Abbrüche).
- `FrameField.add_cross_at_boundaries` → `_add_cross_at_boundaries_fixed`
  (speichert BC unter Knoten-ID statt Schleifenzähler; geweldetes STL hat gestreute IDs).

## Gotchas
- **Außenecken schräg ≠ 90°** im abgewickelten 2D (User-Korrektur). Keine 90°-Annahme.
- **domain_partition downstream buggy** auf echten Singularitäten → eigener Separatrix-Gen.
- **FrameField hat keine Stellschrauben** außer den BC-Werten. Einziger Hebel = BC (gepatcht)
  + optional weicher Alignment-Term (`bc_weight`).
- **CLAUDE.md** nie committen/erwähnen; keine Selbst-Referenz in Commit-Messages.
- **Validator `strict=False`** downgradet Hard-Fails zu Warnungen → für echte Prüfung `strict=True`.

## Referenzen
- Plan: `/root/.claude/plans/generell-claude-md-wird-nie-atomic-candy.md`
- domain_partition Validator: `/root/repos/domain_partition/tools/quad_partition_validator.py`
- Frame-Field: `/root/repos/domain_partition/tools/frame_field.py`
