# Open-source hex meshing libraries for hexdp (desk research, 2026-11-19)

Web research only (nothing downloaded or built). Licences checked first-hand: AlgoHex = AGPL-3.0 (local
/opt/algohex-src/LICENSE), HexDom = BSD-3-Clause and C4HexMeshing = GPLv3 (GitHub pages). The rest as stated on the
projects' pages; "unverified" where no licence file could be read.

## Most useful for the open hexdp problems

| problem in hexdp | candidate | licence | how | note |
|---|---|---|---|---|
| AlgoHex: invalid parametrisation / flipped tets in extraction (HexMe torus, canada gv) | **C4HexMeshing** (Brückler et al. 2023): motorcycle complex + collapse of zero-extent cells -> valid hex mesh from a seamless map | GPLv3 | CLI `c4hex_cli`, input .hexex (tet + seamless map), deps GMP, Clp, MC3D, QGP3D | external program like evocube; AlgoHex can write the seamless map (--igm-out-path) |
| block data model, base complex from a frame field | **gmds** (LIHPC-CEA): n-G-map blocking, SingGraphBuild, sheet ops | LGPL-2.1 | C++ + Python API (pygmds) | comparison / alternative to our beam collapse |
| polycube with holes (genus > 0) | **HexDom** (CMU): polycube -> hex-dominant, arbitrary genus | BSD-3 | multi-stage C++ pipeline, only 6 commits | permissive; maturity unclear |
| untangling / smoothing / quality metrics | **Mesquite** (Sandia) | LGPL-2.1 | C++ library | alternative to our untangle |
| sweep + boundary layers on CAD | **SALOME SMESH** | LGPL-2.1 | Python (smeshBuilder), OCC | large platform dependency |
| robust tet meshing of dirty input | **fTetWild** | MPL-2.0 | CLI | before AlgoHex for STL-like input (Charite) |
| quad layouts (2D) | QuadriFlow (MIT), Instant Meshes (BSD), QuadWild (GPL-3) | | CLI | comparison to layout2d |
| medial axis / skeleton decomposition | CGAL skeletonization / variational medial axis | GPL (package) | C++ | block topology seeding |
| hex-mesh structure simplification | Gao et al. "Robust Hexahedral Re-Meshing" | unverified (likely no licence) | C++ | reference for beam collapse |

## Licence flags
- AGPL: AlgoHex (relevant if hexdp/MCP is ever offered as a network service), TetGen (+ non-commercial terms).
- GPL: evocube, QGP3D, libHexEx, C4HexMeshing, QuadWild, TetWild, gmsh (GPL-2+ with linking exception), pygmsh, OpenFOAM,
  Graphite, libigl/copyleft (CoMISo, CGAL), most CGAL algorithm packages. Use as external programs only.
- LGPL: gmds, OpenVolumeMesh, SMESH, Mmg/ParMmg, Mesquite.
- Permissive: CinoLib (MIT), QuadriFlow (MIT), Instant Meshes (BSD), fTetWild (MPL-2), Geogram (BSD-3), HexDom (BSD-3),
  libigl core (MPL-2), HOHQMesh (MIT).
- Unverified: HexaLab, Q-MAT/SEG-MAT, Gao et al. re-meshing, HexBox code, DualCube code (repo 404), PolycubeNet /
  AI4polycube (May 2026, very new), Auto-Hex (RL + Cubit, no public code found).

Output of a GPL tool (our datasets) is generally not a derivative work; redistributing the tools themselves is.

## Sources
github.com/LIHPC-Computational-Geometry/gmds, github.com/HendrikBrueckler/C4HexMeshing, github.com/CMU-CBML/HexDom,
github.com/sandialabs/mesquite, docs.salome-platform.org, github.com/wildmeshing/fTetWild, github.com/hjwdzh/QuadriFlow,
github.com/wjakob/instant-meshes, github.com/nicopietroni/quadwild, github.com/mlivesu/cinolib, github.com/cgg-bern/AlgoHex,
cgal.org/license.html, github.com/gaoxifeng/Robust-Hexahedral-Re-Meshing, arXiv 2402.00652, 2410.18852, 2605.20274.
