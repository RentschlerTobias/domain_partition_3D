# T1_9 Coarse Block-Structured Meshing Plan

## Context
- **Case**: T1_9 turbine rotor passage
- **Source**: `/root/repos/block_structured_meshing/T1_9/T1_9_ru_gridGmsh.msh` (130,198 nodes, 289,717 elements)
- **OpenFOAM Case**: `/root/repos/block_structured_meshing/T1_9/tistos_ru_of_n_T1_9/`

## Mesh Structure (MSH Format 2.2)

### 2D Elementary Entity Tags (reg_geom)
- **reg_geom=1**: Hub → 2,178 triangles, 11 lines, 1 point
- **reg_geom=2**: Shroud → 7,912 triangles, 35 lines, 1 point
- **reg_geom=3**: Hub BL patch → 1,027 triangles, 612 quads, 11 lines, 1 point
- **reg_geom=4**: Shroud BL patch → 1,019 triangles, 612 quads, 35 lines, 1 point
- **reg_geom=5-6**: Periodic/Blade patches
- **reg_geom=7-31**: Structured quad patches (400-1600 quads each)

### Key Insight
**reg_geom=1 and 2 are the main hub/shroud surfaces** (without boundary layer).
All other reg_geom=3-31 are boundary layer or other patches.

## Current Status

### Completed
1. ✅ Analyzed MSH structure - identified reg_geom tags
2. ✅ Extracted hub and shroud surfaces as STL:
   - `/root/repos/block_structured_meshing/T1_9_hub_raw.stl` (107K, 2,178 tris)
   - `/root/repos/block_structured_meshing/T1_9_shroud_raw.stl` (387K, 7,912 tris)
3. ✅ Set up local Python environment with gmsh, numpy, meshio, networkx, scipy
4. ✅ Docker container: `atismer/dtoo-opensuse:stable` (dtOO + OpenFOAM 2406)
   - LD_LIBRARY_PATH: `/dtOO-install/lib64:/dtOO-install/lib:/usr/lib/openfoam/openfoam2406/platforms/linux64GccDPInt32Opt/lib/sys-openmpi:/usr/lib64/mpi/gcc/openmpi4/lib64:/usr/lib/hpc/gnu7/mpi/openmpi/4.1.6/lib64`
5. ✅ Installed gmsh in container at `/work/gmsh_pkg`

### Blocked
- **Algorithm 11 in container**: Missing `QUADMESHINGTOOLS` module
- **Solution**: Use local Python environment (gmsh 4.15.2 installed)

## Pipeline: Coarse Block Extraction

### Approach: Start from MSH file (not STL)
1. **Load 2D elements** from MSH file with reg_geom=1 and reg_geom=2
2. **Feed into Gmsh** as discrete geometry
3. **Apply Algorithm 11** (Quasi-Structured Quad) with background field
4. **Extract separatrices** from resulting quad mesh
5. **Detect 4-sided blocks** (cycles) in separatrix graph
6. **Export** coarse block structure as MSH

### Alternative: Start from STL
1. Load STL into Gmsh
2. `classifySurfaces()` → convert tri mesh to geometry
3. Algorithm 11 remeshing
4. Separatrix extraction
5. Block detection

## Background Field Strategy

To minimize singularities and control block count:
- Use `gmsh.model.mesh.field.add("MathEval")` 
- Set function to control element distribution
- Goal: ~20 blocks per surface
- Key: Background field affects singularity placement, not just element density

## File Locations
- **Source MSH**: `/root/repos/block_structured_meshing/T1_9/T1_9_ru_gridGmsh.msh`
- **Hub STL**: `/root/repos/block_structured_meshing/T1_9_hub_raw.stl`
- **Shroud STL**: `/root/repos/block_structured_meshing/T1_9_shroud_raw.stl`
- **Env**: `/root/venv` (Python 3.12, gmsh 4.15.2, numpy, meshio, networkx, scipy)
- **Scripts**: `/root/repos/block_structured_meshing/`

## Results

### Step 1: Quad Remeshing (Algorithm 11)
- **Hub**: 832 quads, 920 nodes (elem_size=0.1)
- **Shroud**: 2873 quads, 3016 nodes (elem_size=0.1)
- **Method**: STL → classifySurfaces(pi/2) → createGeometry → Algorithm 11
- **Critical**: curveAngle=pi/2 to merge 204 boundary edges into 2 curves
- **Script**: `remesh_step1.py` with configurable parameters

### Step 2: Separatrix Extraction
- **Hub**: 40 singular vertices, 118 separatrices
- **Shroud**: 48 singular vertices, 159 separatrices
- **Method**: Trace separatrices from singular vertices using cyclic edge order

### Step 3: Block Coarsening
- **Hub**: 18 coarse blocks (from 832 quads, using top 16 separatrices)
- **Shroud**: 21 coarse blocks (from 2873 quads, using top 14 separatrices)
- **Target**: ~20 blocks per surface ✓
- **Method**: Filter separatrices by length, keep only longest N
- **Files**: 
  - `T1_9_hub_blocks.vtk` (18 quads)
  - `T1_9_shroud_blocks.vtk` (21 quads)

## Key Scripts
- `remesh_step1.py`: Step 1 - Quad remeshing with all parameters
- `extract_separatrices.py`: Step 2 - Separatrix extraction
- `coarsen_blocks.py`: Step 3 - Block coarsening
- `plot_blocks.py`: Visualization of coarse blocks

## Key Code Snippets

```python
# Load MSH elements
import gmsh
gmsh.initialize()
# Option 1: Load MSH directly
gmsh.open("T1_9_ru_gridGmsh.msh")
# Then work with surface tags 1,2,3,4... (need identify which is hub/shroud)

# Option 2: Reconstruct from STL
gmsh.merge("T1_9_hub_raw.stl")
gmsh.model.mesh.classifySurfaces()

# Algorithm 11
gmsh.option.setNumber("Mesh.Algorithm", 11)

# Background field
field = gmsh.model.mesh.field.add("MathEval")
gmsh.model.mesh.field.setString(field, "F", "0.5")
gmsh.model.mesh.field.setAsBackgroundMesh(field)

# Generate
gmsh.model.mesh.generate(2)
```

## Docker Commands
```bash
# Run with libraries
export LD_LIBRARY_PATH=/dtOO-install/lib64:/dtOO-install/lib:/usr/lib/openfoam/openfoam2406/platforms/linux64GccDPInt32Opt/lib/sys-openmpi:/usr/lib64/mpi/gcc/openmpi4/lib64:/usr/lib/hpc/gnu7/mpi/openmpi/4.1.6/lib64
source /usr/lib/openfoam/openfoam2406/etc/bashrc

# Python with gmsh
python3.12 -c "import sys; sys.path.insert(0, '/work/gmsh_pkg'); import gmsh"
```

## OpenFOAM Info
- Version: 2406
- Boundary patches in T1_9:
  - RU_HUB: 3,298 faces
  - RU_SHROUD: 9,032 faces
  - RU_BLADE: 4,480 faces
  - RU_INLET: 1,639 faces
  - RU_OUTLET: 1,631 faces
  - RU_PERIOPS: 2,278 faces
  - RU_PERIOSS: 2,278 faces
