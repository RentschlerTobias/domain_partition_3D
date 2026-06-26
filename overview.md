# Block-Structured Meshing for T1_9 — Project Overview

> **For AI Agents:** This file summarizes the project state, architecture, known bugs, and open questions. Read this before making changes.

---

## 1. Goal

Generate a **coarse block-structured quadrilateral mesh** for the **T1_9 turbine rotor passage**.
The primary targets are the **hub** and **shroud** surfaces extracted from the full 3D volume mesh.

The project follows the **quad-dominant block partitioning** pipeline:
1. Surface extraction
2. (Optional) Cross-field computation or Gmsh quasi-structured remeshing
3. Singularity detection
4. Separatrix tracing
5. Block partitioning (finding 4-sided cycles)

---

## 2. Directory Structure

```
/root/repos/block_structured_meshing/
├── README.md                          # Human-readable README
├── PLAN.md                            # T1_9 specific plan & current status
├── requirements.txt                   # Python dependencies (gmsh, meshio, scipy, networkx, ...)
├── main.py                            # General pipeline entry point (modified, uncommitted)
├── run_t1_9_hub.py                    # Monolithic T1_9 pipeline (older)
├── extract_hub_shroud_stl.py          # Extract hub/shroud from MSH reg_geom tags
├── build_blocks.py                    # Block structure builder from streamlines
├── T1_9_hub_raw.stl                   # Extracted hub (2,178 triangles)
├── T1_9_shroud_raw.stl              # Extracted shroud (7,912 triangles)
├── T1_9/                              # Source case files
│   ├── T1_9_ru_gridGmsh.msh           # Source mesh (130,198 nodes, 289,717 elements)
│   ├── stl_parts/                     # Per-part STLs
│   └── tistos_ru_of_*_T1_9/           # OpenFOAM case directories
├── scripts/                           # Pipeline scripts
│   ├── remesh_step1.py                # Gmsh Algorithm 11 quad remeshing
│   ├── extract_separatrices.py        # Extract separatrices from quad mesh
│   ├── coarsen_blocks.py             # Filter separatrices + detect 4-sided blocks
│   ├── phase1_surface_mesh.py         # Phase 1: Load STL, compute normals, boundary
│   ├── phase2_frame_field.py          # Phase 2: Solve cross-frame field on surface
│   ├── phase3_singularities.py      # Phase 3: Detect singularities + trace separatrices
│   └── phase4_blocks.py               # Phase 4: Build graph and find 4-sided faces
├── output/T1_9/
│   ├── hub/                           # Hub outputs
│   └── shroud/                        # Shroud outputs
├── literature/
│   └── Quasi-structured_quadrilateral_meshing_in_Gmsh_--_.pdf
└── html_files/                        # Showcase visualizations
```

---

## 3. Two Pipelines

### 3.1 Working Pipeline (Gmsh-based)
**Status: ✅ Functional**

**Scripts:**
1. `scripts/remesh_step1.py` — Runs Gmsh `Algorithm 11` (QuadQuasiStructured / Packing of Parallelograms). Requires `curveAngle=pi/2` to merge boundary edges, otherwise Gmsh infinite-loops.
2. `scripts/extract_separatrices.py` — Loads the quad mesh, detects singular vertices by **valence**, traces separatrices by following cyclic edge order.
3. `scripts/coarsen_blocks.py` — Filters separatrices by length, builds a graph with boundary edges, detects 4-sided cycles via `networkx`.

**Results:**
- **Hub**: 832 quads, 40 singular vertices, 118 separatrices → **18 coarse blocks**
- **Shroud**: 2873 quads, 48 singular vertices, 159 separatrices → **21 coarse blocks**

### 3.2 Experimental Pipeline (Own 3D Cross-Field)
**Status: ❌ Broken**

**Scripts:**
1. `scripts/phase1_surface_mesh.py` — Loads STL, computes facet normals, local tangent basis (`e1`, `e2`), identifies boundary edges/nodes.
2. `scripts/phase2_frame_field.py` — Solves cross-frame field on the 3D surface via Laplace equation with **Dirichlet boundary conditions**. Uses a **vertex-local basis** (projection onto tangent plane).
3. `scripts/phase3_singularities.py` — Detects singularities via **parallel transport** around faces. Traces separatrices using **RK2** with adaptive step size and **Kanten-Sampling** (sampling along triangle edges).
4. `scripts/phase4_blocks.py` — Builds a graph from separatrices + boundary rings, searches for 4-sided faces (`nx.simple_cycles`).

**Results:**
- **Phase 1**: ✅ Works. Outputs `hub_phase1_mesh.vtk`.
- **Phase 2**: ✅ Works. Outputs `hub_phase2_field.vtk` (1,191 points with frame vectors).
- **Phase 3**: ⚠️ **Broken**. Only **2 singularities** detected (both with index +1). **Expected: ~40 singularities with mixed indices (+1, -1 pairs).**
- **Phase 4**: ⚠️ **Broken**. Only **5 quads** produced. **Expected: ~18 blocks.**

---

## 4. Critical Technical Details

### 4.1 Source Data
- **Source Mesh**: `T1_9_ru_gridGmsh.msh` (130,198 nodes, 289,717 elements, MSH format 2.2)
- **Hub Surface**: `reg_geom=1` — 2,178 triangles, 11 lines, 1 point
- **Shroud Surface**: `reg_geom=2` — 7,912 triangles, 35 lines, 1 point
- **Hub Topology**: Annulus (χ=0). **2 boundary components**: outer ring (112 vertices), inner ring (92 vertices). **~90° corners** at the boundaries.

### 4.2 Cross-Field Representation
- The solver outputs a **representation vector** `(u_x, u_y)` where the field angle is `φ = 4θ`.
- The **Poincaré index in representation space** is **4× the geometric index**.
- A geometric index of `+1/4` corresponds to representation index `+1`.
- A geometric index of `-1/4` corresponds to representation index `-1`.

### 4.3 Singularity & Separatrix Rules
- **Singularities must appear in pairs** (Poincaré-Hopf theorem for a surface with boundary).
- **A singularity has either 3 or 5 separatrices**.
- The current code in `phase3_singularities.py` is **confused about the count formula**. The correct count depends on the geometric index and must be consistent with the representation index.
- **Boundary singularities**: If a singularity lies on the boundary (e.g., all 3 vertices of a face are on the boundary), some of its separatrices will run along the boundary and must be handled separately.

### 4.4 Boundary Conditions
- **Dirichlet boundary conditions** are applied on the boundary nodes.
- **Corner handling**: At C⁰ corners (e.g., 90°), the representation vectors of adjacent edges are averaged. For exact 90° corners, the vectors are identical and the field is **topologically singular** at the corner (this is a feature, not a bug).
- **Kowalski et al. 2015** approach: Dirichlet on ∂Ω with corner averaging. **Knöppel et al. 2013** uses natural Neumann boundary conditions.
- **Current approach**: No Dirichlet condition at the corners (free nodes). This causes singularities to move numerically into the interior, which is more stable.

### 4.5 Separatrix Tracing
- **Algorithm**: RK2 (2nd-order Runge-Kutta) with **adaptive step size**.
- **Step size**: Starts at `0.5 × avg_edge_length`. Halves when the predictor step goes outside the mesh.
- **Termination**: Tracing stops when the point leaves the mesh (barycentric coordinates fail) or lands on a boundary edge within distance `1e-6`.
- **Sampling**: Uses **Kanten-Sampling** (sampling along triangle edges) to find the initial separatrix directions from a singularity, instead of Kreis-Sampling.

### 4.6 Phase 4 Graph Problem
- **Current graph**: 7 nodes, 5 edges. **0 four-sided cycles.**
- **Root cause**: The 4 traced separatrices do **not cross each other**. They go directly from the singularity to the boundary.
- **To fix**: Need either more singularities (which will create more separatrices that cross), or boundary-corner singularities (to create additional separatrices from the 90° corners).

---

## 5. Known Bugs & Blockers

### 5.1 Phase 3 — Singularity Detection
- **Bug**: Only 2 singularities found (both with representation index +1). Expected: ~40 singularities with mixed +1 / -1 indices.
- **Location**: `scripts/phase3_singularities.py`.
- **Symptoms**: The parallel transport around faces is not detecting the correct topological defects. The sum of indices does not match the Euler characteristic.

### 5.2 Phase 3 — Separatrix Count
- **Bug**: The code is unsure whether the formula is `abs(2k - 1)` or `abs(index * 4 - 1)`.
- **Fix needed**: The number of separatrices must be **3 or 5** for each singularity, depending on its geometric index.

### 5.3 Phase 4 — Block Partitioning
- **Bug**: `phase4_blocks.py` produces only 5 quads instead of ~18 blocks.
- **Location**: `scripts/phase4_blocks.py`.
- **Symptoms**: Graph is too sparse (7 nodes, 5 edges). `nx.simple_cycles` finds no valid 4-cycles.
- **Fix needed**: Either get more separatrices (from more singularities), or implement a more sophisticated graph simplification (e.g., porting `streamline_simplificator.py` from `domain_partition`).

### 5.4 Visualization
- **Open**: The user requested exporting the **actual cross-field directions** (4 directions per vertex) instead of the representation vector. This is not yet implemented.

---

## 6. Next Steps / Open Questions

### Immediate (To Fix Phase 3)
1. **Debug singularity detection**: Verify the parallel transport and the barycentric interpolation for singularity locations against a known simple geometry (e.g., a flat disk with 4 corners).
2. **Fix the separatrix count formula**: Ensure it outputs **3 or 5** separatrices per singularity.
3. **Find the missing singularities**: The sum of Poincaré indices must match the Euler characteristic (`χ = 0` for the annulus). We are missing singularities with index **-1** (representation space).

### Short-term (To Fix Phase 4)
1. **Get separatrices to cross**: This requires finding the missing singularities (see above).
2. **Boundary corner singularities**: Implement detection of 90° corners as boundary singularities. These emit separatrices that can create the necessary crossings.
3. **Alternative graph algorithm**: If the graph remains too sparse, consider porting `streamline_simplificator.py` from the `domain_partition` reference.

### Long-term
1. **Visualization**: Export cross-field directions (4 vectors per vertex) to VTK.
2. **Shroud surface**: Run the fixed pipeline on the shroud (7,912 triangles).
3. **Full passage**: Extend the pipeline to handle the full 3D passage (not just hub/shroud).

---

## 7. Environment

- **OS**: Linux
- **Python**: 3.12
- **Branch**: `dtoo`
- **Key Packages**: `gmsh 4.15.2`, `meshio 5.3.5`, `networkx 3.6.1`, `scipy 1.17.1`, `numpy 2.4.4`
- **Working Directory**: `/root/repos/block_structured_meshing`
- **Output Directory**: `output/T1_9/hub/` and `output/T1_9/shroud/`

---

## 8. Literature

- **Kowalski et al. 2015**: `literature/quad_domain_partition.pdf` (boundary conditions, corner handling, separatrix count formula).
- **Gmsh Algorithm 11**: `literature/Quasi-structured_quadrilateral_meshing_in_Gmsh_--_.pdf` (QuadQuasiStructured / Packing of Parallelograms).
- **Knöppel et al. 2013**: "Stripe Patterns on Surfaces" (natural Neumann boundary conditions for cross fields).

---

## 9. Important Notes for AI Agents

- **Do NOT use the cross-field pipeline (`phase1-4`) for production.** It is experimental and currently broken. Use the **Gmsh-based pipeline** (`remesh_step1.py` → `extract_separatrices.py` → `coarsen_blocks.py`) if you need working blocks.
- **When debugging Phase 3**: Pay close attention to the **parallel transport angle**. The angle must be accumulated correctly around each face. The current code uses `e1/e2` from Phase 1 (do not use `np.cross` to reconstruct the basis).
- **When modifying Phase 4**: The graph is built from separatrices + boundary edges. If there are not enough crossings, the block detection will fail. The primary fix is in Phase 3 (finding more singularities).
- **Corner singularities**: The 90° corners on the hub boundary are **topologically required** to be singular (Poincaré-Hopf). The field is intentionally unsmooth there. Do not try to "smooth them out".
- **Git status**: There are uncommitted changes in `main.py`, `README.md`, `testing_s.py`, and `tistos_quasi_struc_msh.py`. Check `git status` before committing.
