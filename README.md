# Quadrilateral Block‑Structure Generator

A lightweight workflow that turns a 2‑D surface (STL or MSH) into a quadrilateral block‑structure.  
[Blocking Case09_post/LV_outer](./html_files/quad_blocking.html)
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

## 3. Automated Python Workflow

### 3.1 Installation

```bash
pip install -r requirements.txt
```

The repository ships with a `requirements.txt` containing (tested using python v3.13.12):
- `numpy`
- `gmsh` 
- `scipy` 
- `openmesh` 
- `networkx` 
- `meshio`
- `plotly` 

### 3.2 Usage

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

### 3.3 Class & Function Overview

| Component | Purpose |
|-----------|---------|
| `MshExtractor` | Automates the GUI workflow: STL → remeshing → triangular mesh → quasi‑structured quad mesh. Returns vertex and surface data. |
| `StreamlineExtractor` | Detects irregular vertices (degree ≠ 4) and builds streamlines along edges of the quadrilateral mesh. |
| `split_streamlines` | Finds intersection points between streamlines and cuts them there, yielding clean block boundaries. |
| `export_streamlines_to_msh` | Writes a new MSH file containing only the extracted streamlines (useful for debugging). |
| `block_structure_to_html` | Generates an interactive HTML page showing the streamlines with optional linear or bilinear interpolation. |
| `faces_to_html` | Renders the full quadrilateral block structure interactively; can interpolate faces (`n_u`, `n_v`) for smoother visualisation. |

---

## 4. Example Output

| Step | Result |
|------|--------|
| **Remeshed triangular mesh** | ![triangular](docs/triangular.png) |
| **Quasi‑structured quad mesh** | ![quad](docs/quadrilateral.png) |
| **Streamlines only (MSH)** | `streamlines.msh` – can be opened in Gmsh. |
| **Interactive HTML** | <a href="streamlines.html">View streamlines</a> |
| **Full block structure** | <a href="quadrilateral_blocks.html">View blocks</a> |

---

## 5. Troubleshooting

| Symptom | Likely Cause | Fix |
|---------|--------------|-----|
| Gmsh crashes on `remeshing.geo` | Wrong path to STL or missing Quad‑Tools build | Verify the file exists and that you are using the compiled Gmsh with Quad‑Tools. |
| No streamlines produced | All vertices have degree 4 (mesh is already regular) | Check `angle`/`curve_angle`; try lowering `element_size`. |
| HTML page shows no geometry | `n_u`, `n_v` set to `None` but interpolation code expects numeric values | Pass integer values or modify the function call. |

---
