# Quadrilateral Block‑Structure Generator

A lightweight workflow that turns a 3‑D surface (STL or MSH) into a quadrilateral block‑structure.  


The pipeline is split into two parts:

* **GUI‑based** – easy to follow Gmsh steps for quick experimentation.  
* **Python‑based** – executed using 'python main.py', fully automated, extented functionallity.


---

## 1. What main.py does

1. **STL → Gmsh remeshing** (triangular mesh).  
2. **Triangular mesh → Quasi‑Structured Quad** (Gmsh algorithm).  
3. **Extract irregular vertices** – any vertex that is not shared by exactly four faces.  
4. **Generate “streamlines”** – chains of edges originating from the irregular vertices.  
5. **Split streamlines at intersection points** to obtain the final quadrilateral block structure.

The result can be visualised in Gmsh, exported as an MSH file or rendered interactively in a web browser (HTML).

Showcase HTML file of Blocking to open in any browser
-./html_files/Case09_post_LV_inner.html
-./html_files/Case09_post_LV_inner.html
Blocking applied to Case199_pre_all.stl (really slow rendering) with autmatically subdevided surface plotted individually

-./html_files/blocking_faces_independent_surfaces_s2.html
-./html_files/blocking_faces_independent_surfaces_s3.html
-./html_files/blocking_faces_independent_surfaces_s4.html


---

## 2. GUI Workflow

> *Prerequisite:* compile the C++ Gmsh code with **Quad‑Tools** enabled.

1. **Load geometry** and remshes to export stl file to msh file 
   ```bash
   gmsh remeshing.geo
   ```
   Inside `remeshing.geo` point to your STL:
   ```geo
   Merge "./stl_files/Case09_post/LV_outer.stl";
   ```

2. **Set meshing options**  
   * Right‑click → **View all options** → **Mesh → General → 2D → QuasiStructured‑Quad (experimental)**  
   * In the Gmsh main window left side, click **Register → Mesh → 2D** to generate a quasi‑structured quad mesh.

3. **Parameters**  
   `angle`, `curve_angle` and `element_size` are defined and described in `remeshing.geo`.  
   The remeshing splits the geometry into several surfaces depending on these parameters.

4. **Export** – Save the resulting MSH file (triangular or quadrilateral) for later use.

---

## 3. T1_9 Turbine Hub Pipeline

A dedicated pipeline for the T1_9 test case extracts the hub surface, remeshes it with the quasi-structured quad algorithm, and generates a coarse block structure.

```bash
python run_t1_9_hub.py
```

This produces:
- `T1_9_hub.stl` – extracted hub surface (z=0)
- `T1_9_hub_remeshed.msh` – quasi-structured quad mesh
- `T1_9_hub_blocks.vtk` – coarse block structure (4-sided blocks)

### 3.1 Automated Python Workflow

### 3.2 Installation

```bash
pip install -r requirements.txt
```

The repository ships with a `requirements.txt` containing (tested using python v3.12.3):
- `numpy`
- `gmsh` 
- `scipy` 
- `openmesh` 
- `networkx` 
- `meshio`
- `plotly` 

### 3.3 Usage

```python
from msh_extractor import MshExtractor
from streamline_extractor import StreamlineExtractor
from split_streamlines import split_streamlines
from plotting_tools import export_streamlines_to_msh, block_structure_to_html, faces_to_html

# ---- Step 1: Extract the mesh ------------------------------------
input_path   = './stl_files/Case09_post/LV_outer.stl'
output_path  = './msh_intermediate.msh'

extractor = MshExtractor(
    input_path=input_path,
    remesh=True,                # run Gmsh remeshing
    output_path=output_path,
    element_size=2.5,
    angle=50,
    curve_angle=180
)

vertices = extractor.vertices      # (N, 3) numpy array
surfaces = extractor.surfaces      # dict: surface_id -> (n_faces, 4)-array of quad connectivity

# ---- Step 2: Generate streamlines --------------------------------
streamline_extractor = StreamlineExtractor(vertices, surfaces)
streamlines = streamline_extractor.streamlines   # list of line segments

# ---- Step 3: Split at intersection points ------------------------
splitted_streamlines = split_streamlines(streamlines)

# ---- Step 4: Export / Visualise -----------------------------------
export_streamlines_to_msh(splitted_streamlines,
                          output_path='./streamlines.msh')

block_structure_to_html(
    splitted_streamlines,
    output_path='streamlines.html'
)

faces_to_html(
    vertices, surfaces,
    n_u=20, n_v=20,
    output_path='quadrilateral_blocks.html'
)
```

### 3.4 Class & Function Overview

| Component | Purpose |
|-----------|---------|
| `MshExtractor` | Automates the GUI workflow: STL → remeshing → triangular mesh → quasi‑structured quad mesh. Returns vertex and surface data. |
| `StreamlineExtractor` | Detects irregular vertices (valenz ≠ 4) and builds streamlines along edges of the quadrilateral mesh. |
| `split_streamlines` | Finds intersection points between streamlines and cuts them there, yielding the block boundaries. |
| `export_streamlines_to_msh` | Writes a new MSH file containing only the extracted streamlines (useful for debugging). |
| `block_structure_to_html` | Generates an interactive HTML page showing the streamlines with optional linear or bilinear interpolation. |
| `faces_to_html` | Renders the full quadrilateral block structure interactively; can interpolate faces (`n_u`, `n_v`) for smoother visualisation. |

---

