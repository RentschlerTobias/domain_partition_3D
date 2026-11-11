#!/usr/bin/env python3
"""Generate a VTK Multiblock (.vtm) with all tistos pipeline stages for debugging.

Each pipeline stage becomes a named component in the multiblock, allowing
isolation of the step that introduced a bug. Each component includes cell
fields: block_id, signed_volume, min_scaled_jacobian, and inverted.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import meshio

# Add parent directory to import helpers
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent.parent))

import ovm_io
import clean_blocks as cb


def hex_signed_volume_divtheorem(pts):
    """Signed volume by divergence theorem with standard VTK face winding.
    
    Uses the same formula as audit_tfi_blade.py for consistency.
    """
    VTK_HEX_FACES = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 4, 7, 3),
                     (1, 2, 6, 5), (0, 1, 5, 4), (2, 3, 7, 6)]
    V = 0.0
    for f in VTK_HEX_FACES:
        a = pts[f[0]]; b = pts[f[1]]; c = pts[f[2]]; d = pts[f[3]]
        V += float(np.dot(a, np.cross(b, c)))
        V += float(np.dot(a, np.cross(c, d)))
    return V / 6.0


def compute_mesh_metrics(points, cells, cell_type="hexahedron"):
    """Compute block_id, signed_volume, min_scaled_jacobian, and inverted flag.
    
    Args:
        points: (N, 3) array of point coordinates
        cells: (M, K) array of cell connectivity (K=8 for hex, K=4 for tet)
        cell_type: "hexahedron" or "tetra"
    
    Returns:
        dict with arrays: block_id, signed_volume, min_scaled_jacobian, inverted
    """
    n_cells = len(cells)
    
    # Initialize metrics
    block_id = np.zeros(n_cells, dtype=int)
    signed_volume = np.zeros(n_cells, dtype=float)
    min_scaled_jacobian = np.zeros(n_cells, dtype=float)
    inverted = np.zeros(n_cells, dtype=int)
    
    if n_cells == 0:
        return {
            "block_id": block_id,
            "signed_volume": signed_volume,
            "min_scaled_jacobian": min_scaled_jacobian,
            "inverted": inverted,
        }
    
    # Compute metrics based on cell type
    if cell_type == "hexahedron":
        # Signed volume for each hex
        for i in range(n_cells):
            pts = points[cells[i]]
            signed_volume[i] = hex_signed_volume_divtheorem(pts)
        
        # Scaled Jacobian (vectorized)
        min_scaled_jacobian = cb.scaled_jacobians(points, cells)
        
        # Inverted flag (sJ < 0)
        inverted = (min_scaled_jacobian < 0).astype(int)
    
    elif cell_type == "tetra":
        # For tetrahedra, compute volume and a simplified Jacobian
        for i in range(n_cells):
            pts = points[cells[i]]
            # Volume of tetrahedron
            a, b, c, d = pts
            vol = np.dot(b - a, np.cross(c - a, d - a)) / 6.0
            signed_volume[i] = vol
            
            # Simplified Jacobian: for tets, use volume sign as quality indicator
            # A negative volume indicates an inverted tet
            min_scaled_jacobian[i] = 1.0 if vol > 0 else -1.0
            inverted[i] = 1 if vol <= 0 else 0
    
    return {
        "block_id": block_id,
        "signed_volume": signed_volume,
        "min_scaled_jacobian": min_scaled_jacobian,
        "inverted": inverted,
    }


def read_mesh_with_metrics(filepath, cell_type="hexahedron"):
    """Read a mesh file and compute metrics.
    
    Args:
        filepath: Path to mesh file (.vtk, .vtu, .msh)
        cell_type: Expected cell type ("hexahedron" or "tetra")
    
    Returns:
        dict with points, cells, and metrics, or None if file not found
    """
    if not Path(filepath).exists():
        return None
    
    try:
        mesh = meshio.read(str(filepath))
    except Exception as e:
        print(f"ERROR reading {filepath}: {e}")
        return None
    
    # Extract cells of the requested type
    cells = []
    for cell_block in mesh.cells:
        if cell_block.type == cell_type:
            cells.append(cell_block.data)
    
    if not cells:
        print(f"WARNING: No {cell_type} cells found in {filepath}")
        return None
    
    cells_array = np.vstack(cells)
    points = np.asarray(mesh.points, dtype=float)
    
    # Extract block_id if present
    block_id = None
    if mesh.cell_data:
        for cell_block, data_dict in zip(mesh.cells, mesh.cell_data):
            if cell_block.type == cell_type and isinstance(data_dict, dict):
                for key in data_dict:
                    if key.lower() == 'block_id':
                        block_id = np.concatenate([np.asarray(x).ravel() 
                                                   for x in data_dict[key]])
                        break
                if block_id is not None:
                    break
    
    # Compute metrics
    metrics = compute_mesh_metrics(points, cells_array, cell_type)
    
    # Override block_id if found in file
    if block_id is not None:
        metrics["block_id"] = block_id
    
    return {
        "points": points,
        "cells": cells_array,
        "metrics": metrics,
        "cell_type": cell_type,
    }


def read_ovm_with_metrics(filepath):
    """Read an OVM file and convert to VTK with metrics.
    
    Args:
        filepath: Path to .ovm file
    
    Returns:
        dict with points, cells (hexahedra), and metrics, or None if file not found
    """
    if not Path(filepath).exists():
        return None
    
    try:
        points, edges, faces, polys = ovm_io.read_ovm(filepath)
        hexes, skipped = ovm_io.ovm_to_cells(points, edges, faces, polys)
        
        if len(hexes) == 0:
            print(f"WARNING: No hexahedra recovered from {filepath}")
            return None
        
        metrics = compute_mesh_metrics(points, hexes, cell_type="hexahedron")
        
        return {
            "points": points,
            "cells": hexes,
            "metrics": metrics,
            "cell_type": "hexahedron",
        }
    except Exception as e:
        print(f"ERROR reading {filepath}: {e}")
        return None


def create_vtm_block(mesh_data):
    """Create a meshio Mesh object from mesh_data.
    
    Args:
        mesh_data: dict with points, cells, metrics, cell_type
    
    Returns:
        meshio.Mesh object with cell_data
    """
    cell_type = mesh_data["cell_type"]
    cells = [(cell_type, mesh_data["cells"])]
    
    cell_data = {
        "block_id": [mesh_data["metrics"]["block_id"]],
        "signed_volume": [mesh_data["metrics"]["signed_volume"]],
        "min_scaled_jacobian": [mesh_data["metrics"]["min_scaled_jacobian"]],
        "inverted": [mesh_data["metrics"]["inverted"]],
    }
    
    return meshio.Mesh(
        points=mesh_data["points"],
        cells=cells,
        cell_data=cell_data,
    )


def extract_inverted_cells(mesh_data):
    """Extract only inverted cells (sJ < 0) as a separate mesh.
    
    Args:
        mesh_data: dict with points, cells, metrics
    
    Returns:
        meshio.Mesh with only inverted cells, or None if no inverted cells
    """
    inverted_mask = mesh_data["metrics"]["inverted"] > 0
    if not np.any(inverted_mask):
        return None
    
    inverted_cells = mesh_data["cells"][inverted_mask]
    if len(inverted_cells) == 0:
        return None
    
    # Create mesh with only inverted cells
    cells = [(mesh_data["cell_type"], inverted_cells)]
    cell_data = {
        "block_id": [mesh_data["metrics"]["block_id"][inverted_mask]],
        "signed_volume": [mesh_data["metrics"]["signed_volume"][inverted_mask]],
        "min_scaled_jacobian": [mesh_data["metrics"]["min_scaled_jacobian"][inverted_mask]],
        "inverted": [mesh_data["metrics"]["inverted"][inverted_mask]],
    }
    
    return meshio.Mesh(
        points=mesh_data["points"],
        cells=cells,
        cell_data=cell_data,
    )


def build_pipeline_stages(run_dir):
    """Build all pipeline stages from a run directory.
    
    Args:
        run_dir: Path to the run directory
    
    Returns:
        dict with stage data and metadata
    """
    run_dir = Path(run_dir)
    
    stages = {}
    
    # Stage 1: input (case/mesh.msh)
    input_path = run_dir / "case" / "mesh.msh"
    if input_path.exists():
        print(f"[pipeline] Reading input: {input_path}")
        mesh_data = read_mesh_with_metrics(input_path, cell_type="tetra")
        if mesh_data is not None:
            stages["input"] = {
                "source": str(input_path),
                "mesh_data": mesh_data,
            }
    
    # Stage 2: tet
    tet_path = run_dir / "tet.vtk"
    if tet_path.exists():
        print(f"[pipeline] Reading tet: {tet_path}")
        mesh_data = read_mesh_with_metrics(tet_path, cell_type="tetra")
        if mesh_data is not None:
            stages["tet"] = {
                "source": str(tet_path),
                "mesh_data": mesh_data,
            }
    
    # Stage 3: hexamesh (hex_hex_x_16c.ovm or fallback hex_hex.ovm)
    hexamesh_path = run_dir / "hex_hex_x_16c.ovm"
    if not hexamesh_path.exists():
        hexamesh_path = run_dir / "hex_hex.ovm"
    if hexamesh_path.exists():
        print(f"[pipeline] Reading hexamesh: {hexamesh_path}")
        mesh_data = read_ovm_with_metrics(hexamesh_path)
        if mesh_data is not None:
            stages["hexamesh"] = {
                "source": str(hexamesh_path),
                "mesh_data": mesh_data,
            }
    
    # Stage 4: sheet_collapse (blocks.vtk)
    blocks_path = run_dir / "blocks.vtk"
    if blocks_path.exists():
        print(f"[pipeline] Reading sheet_collapse: {blocks_path}")
        mesh_data = read_mesh_with_metrics(blocks_path, cell_type="hexahedron")
        if mesh_data is not None:
            stages["sheet_collapse"] = {
                "source": str(blocks_path),
                "mesh_data": mesh_data,
            }
    
    # Stage 5: blocks_linear (blade/blocks_blade.vtk)
    blocks_blade_path = run_dir / "blade" / "blocks_blade.vtk"
    if blocks_blade_path.exists():
        print(f"[pipeline] Reading blocks_linear: {blocks_blade_path}")
        mesh_data = read_mesh_with_metrics(blocks_blade_path, cell_type="hexahedron")
        if mesh_data is not None:
            stages["blocks_linear"] = {
                "source": str(blocks_blade_path),
                "mesh_data": mesh_data,
            }
    
    # Stage 6: blade_insert (blade/blocks_blade_smooth.vtk or blocks_blade_const.vtk)
    blade_insert_path = run_dir / "blade" / "blocks_blade_smooth.vtk"
    if not blade_insert_path.exists():
        blade_insert_path = run_dir / "blade" / "blocks_blade_const.vtk"
    if blade_insert_path.exists():
        print(f"[pipeline] Reading blade_insert: {blade_insert_path}")
        mesh_data = read_mesh_with_metrics(blade_insert_path, cell_type="hexahedron")
        if mesh_data is not None:
            stages["blade_insert"] = {
                "source": str(blade_insert_path),
                "mesh_data": mesh_data,
            }
    
    # Stage 7: block_edges_curved (blade/blocks_blade_curved.vtk)
    curved_path = run_dir / "blade" / "blocks_blade_curved.vtk"
    if curved_path.exists():
        print(f"[pipeline] Reading block_edges_curved: {curved_path}")
        mesh_data = read_mesh_with_metrics(curved_path, cell_type="hexahedron")
        if mesh_data is not None:
            stages["block_edges_curved"] = {
                "source": str(curved_path),
                "mesh_data": mesh_data,
            }
    
    # Stage 8: tfi (tfi.vtk and blade/tfi_blade.vtk)
    # Main tfi.vtk
    tfi_path = run_dir / "tfi.vtk"
    if tfi_path.exists():
        print(f"[pipeline] Reading tfi: {tfi_path}")
        mesh_data = read_mesh_with_metrics(tfi_path, cell_type="hexahedron")
        if mesh_data is not None:
            stages["tfi"] = {
                "source": str(tfi_path),
                "mesh_data": mesh_data,
            }
    
    # tfi_blade.vtk (blade-specific)
    tfi_blade_path = run_dir / "blade" / "tfi_blade.vtk"
    if tfi_blade_path.exists():
        print(f"[pipeline] Reading tfi_blade: {tfi_blade_path}")
        mesh_data = read_mesh_with_metrics(tfi_blade_path, cell_type="hexahedron")
        if mesh_data is not None:
            stages["tfi_blade"] = {
                "source": str(tfi_blade_path),
                "mesh_data": mesh_data,
            }
    
    return stages


def write_vtm_and_json(stages, output_dir):
    """Write VTM multiblock and JSON metadata.
    
    Args:
        stages: dict of stage data from build_pipeline_stages
        output_dir: Path to output directory
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Build VTM structure
    vtm_blocks = []
    json_records = []
    all_inverted = []
    
    for stage_name in ["input", "tet", "hexamesh", "sheet_collapse", 
                       "blocks_linear", "blade_insert", "block_edges_curved", 
                       "tfi", "tfi_blade"]:
        if stage_name not in stages:
            # Stage is missing
            json_records.append({
                "stage": stage_name,
                "source": None,
                "status": "missing",
                "n_cells": 0,
                "min_sJ": None,
                "n_inverted": 0,
                "n_volume_le0": 0,
            })
            continue
        
        stage = stages[stage_name]
        mesh_data = stage["mesh_data"]
        
        # Create mesh block
        mesh = create_vtm_block(mesh_data)
        vtm_blocks.append((stage_name, mesh))
        
        # Extract inverted cells
        inv_mesh = extract_inverted_cells(mesh_data)
        if inv_mesh is not None:
            vtm_blocks.append((f"{stage_name}_inverted", inv_mesh))
            all_inverted.append(inv_mesh)
        
        # Compute statistics
        metrics = mesh_data["metrics"]
        n_cells = len(mesh_data["cells"])
        min_sj = float(np.min(metrics["min_scaled_jacobian"])) if n_cells > 0 else None
        n_inverted = int(np.sum(metrics["inverted"] > 0))
        n_volume_le0 = int(np.sum(metrics["signed_volume"] <= 0))
        
        json_records.append({
            "stage": stage_name,
            "source": stage["source"],
            "status": "present",
            "n_cells": n_cells,
            "min_sJ": min_sj,
            "n_inverted": n_inverted,
            "n_volume_le0": n_volume_le0,
        })
    
    # Write VTM file using meshio's multiblock writer
    vtm_path = output_dir / "pipeline_stages.vtm"
    
    # Create a simple VTM XML structure manually
    vtm_xml = ['<?xml version="1.0" encoding="UTF-8"?>']
    vtm_xml.append('<VTKFile type="vtkMultiBlockDataSet" version="1.0">')
    vtm_xml.append('  <vtkMultiBlockDataSet>')
    
    for idx, (name, mesh) in enumerate(vtm_blocks):
        # Write each block as a separate .vtu file
        block_filename = f"{name}.vtu"
        block_path = output_dir / block_filename
        mesh.write(str(block_path))
        
        vtm_xml.append(f'    <Block name="{name}">')
        vtm_xml.append(f'      <DataSet index="0" file="{block_filename}"/>')
        vtm_xml.append(f'    </Block>')
    
    vtm_xml.append('  </vtkMultiBlockDataSet>')
    vtm_xml.append('</VTKFile>')
    
    vtm_path.write_text('\n'.join(vtm_xml))
    print(f"[pipeline] Wrote {vtm_path}")
    
    # Write JSON metadata
    json_path = output_dir / "pipeline_stages.json"
    json_path.write_text(json.dumps(json_records, indent=2))
    print(f"[pipeline] Wrote {json_path}")
    
    return vtm_path, json_path


def main():
    ap = argparse.ArgumentParser(description="Generate pipeline stages VTM for debugging")
    ap.add_argument("--run-dir", required=True,
                    help="Run directory with pipeline outputs")
    ap.add_argument("--output-dir", default=None,
                    help="Output directory (default: <run-dir>/blade/viz)")
    args = ap.parse_args()
    
    run_dir = Path(args.run_dir)
    output_dir = Path(args.output_dir) if args.output_dir else run_dir / "blade" / "viz"
    
    print(f"[pipeline] Building pipeline stages from {run_dir}")
    stages = build_pipeline_stages(run_dir)
    
    print(f"[pipeline] Found {len(stages)} stages")
    vtm_path, json_path = write_vtm_and_json(stages, output_dir)
    
    print(f"\n[pipeline] Done!")
    print(f"[pipeline] VTM: {vtm_path}")
    print(f"[pipeline] JSON: {json_path}")


if __name__ == "__main__":
    main()
