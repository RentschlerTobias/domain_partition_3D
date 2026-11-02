#!/usr/bin/env python3
"""Audit tfi_blade.vtk and make inversions visible."""
import argparse
import json
import hashlib
import os
from pathlib import Path
from collections import defaultdict, Counter

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import sys
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent.parent))

import base_complex as bc
import clean_blocks as cb
import tfi
import blade_check


DEFAULT_VTK = Path("/work/datagen_keepbl/runs/tistos_s1_0131/blade/tfi_blade.vtk")
DEFAULT_OUT_JSON = Path("/work/datagen_keepbl/runs/tistos_s1_0131/blade/tfi_blade_audit.json")
DEFAULT_OUT_VTK = Path("/work/datagen_keepbl/runs/tistos_s1_0131/blade/tfi_blade_audit.vtk")
DEFAULT_OUT_PNG = Path("/work/datagen_keepbl/runs/tistos_s1_0131/blade/tfi_blade_audit.png")
DEFAULT_S2C = Path("/work/datagen_keepbl/runs/tistos_s1_0131/blade/blade_s2c.json")


VTK_HEX_FACES = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 4, 7, 3),
                 (1, 2, 6, 5), (0, 1, 5, 4), (2, 3, 7, 6)]


def file_fingerprint(path):
    st = os.stat(path)
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        h.update(f.read())
    return {
        "path": str(path),
        "mtime": st.st_mtime,
        "size": st.st_size,
        "sha256": h.hexdigest(),
    }


def read_vtk_full(path):
    import meshio
    m = meshio.read(str(path))
    H = np.vstack([c.data for c in m.cells if c.type == "hexahedron"])
    bid = None
    # try to get block_id from cell_data
    for c, d in zip(m.cells, m.cell_data):
        if c.type == "hexahedron" and isinstance(d, dict):
            for k in d:
                if k == 'block_id' or k.lower() == 'block_id':
                    bid = np.concatenate([np.asarray(x).ravel() for x in d[k]])
                    break
            if bid is not None:
                break
    if bid is None and hasattr(m, 'cell_data'):
        for k in m.cell_data:
            if k == 'block_id' or k.lower() == 'block_id':
                try:
                    data = m.cell_data[k]
                    # data might be list of arrays
                    if isinstance(data, (list, tuple)):
                        bid = np.concatenate([np.asarray(x).ravel() for x in data])
                    else:
                        bid = np.asarray(data).ravel()
                    break
                except Exception:
                    pass
    P = np.asarray(m.points, float)
    return P, H, bid


def hex_signed_volume_divtheorem(pts):
    """Signed volume by divergence theorem with standard VTK face winding."""
    V = 0.0
    for f in VTK_HEX_FACES:
        a = pts[f[0]]; b = pts[f[1]]; c = pts[f[2]]; d = pts[f[3]]
        V += float(np.dot(a, np.cross(b, c)))
        V += float(np.dot(a, np.cross(c, d)))
    return V / 6.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vtk", default=str(DEFAULT_VTK))
    ap.add_argument("--out-json", default=str(DEFAULT_OUT_JSON))
    ap.add_argument("--out-vtk", default=str(DEFAULT_OUT_VTK))
    ap.add_argument("--out-png", default=str(DEFAULT_OUT_PNG))
    args = ap.parse_args()
    
    vtk_path = Path(args.vtk)
    print(f"[audit] Reading {vtk_path}")
    P, H, bid = read_vtk_full(vtk_path)
    n_points = len(P)
    n_cells = len(H)
    
    fp = file_fingerprint(vtk_path)
    
    # block info
    block_ids = []
    if bid is not None:
        block_ids = sorted(set(int(x) for x in bid))
    n_blocks = len(block_ids) if block_ids else (int(bid.max()) + 1 if bid is not None else 0)
    
    print(f"[audit] points={n_points}, cells={n_cells}, blocks={n_blocks}, block_ids={block_ids}")
    
    # Compute metrics for every hex
    vols = []
    for i in range(n_cells):
        pts = P[H[i]]
        vols.append(hex_signed_volume_divtheorem(pts))
    vols = np.array(vols)
    
    # Use project's scaled Jacobian
    sj_mins = cb.scaled_jacobians(P, H)
    
    n_vol_le0 = int((vols <= 0).sum())
    n_sj_neg = int((sj_mins < 0).sum())
    n_sj_lt03 = int((sj_mins < 0.3).sum())
    
    min_sj = float(np.min(sj_mins)) if len(sj_mins) > 0 else float('nan')
    max_sj = float(np.max(sj_mins)) if len(sj_mins) > 0 else float('nan')
    min_vol = float(np.min(vols)) if len(vols) > 0 else float('nan')
    max_vol = float(np.max(vols)) if len(vols) > 0 else float('nan')
    
    # find worst cells
    worst_idx = np.argsort(sj_mins)[:20]
    worst = []
    for idx in worst_idx:
        c = int(idx)
        centre = P[H[c]].mean(axis=0).tolist()
        b = int(bid[c]) if bid is not None else -1
        worst.append({
            "cell_id": c,
            "block_id": b,
            "centre": centre,
            "volume": float(vols[c]),
            "sJ": float(sj_mins[c]),
        })
    
    print(f"[audit] n(volume<=0)={n_vol_le0}, vol=[{min_vol:.6g}, {max_vol:.6g}]")
    print(f"[audit] n(sJ<0)={n_sj_neg}, min sJ={min_sj:.6g}, n(sJ<0.3)={n_sj_lt03}")
    print(f"[audit] worst: {[(w['cell_id'], w['block_id'], w['sJ']) for w in worst[:10]]}")
    
    # Cross-check with blade_s2c.json
    crosscheck = {}
    s2c = DEFAULT_S2C
    if s2c.exists():
        try:
            s2c_data = json.loads(s2c.read_text())
            rep_inv = s2c_data.get('tfi', {}).get('inverted_count')
            rep_min = s2c_data.get('tfi', {}).get('min_scaled_jacobian')
            crosscheck['blade_s2c_inverted_count_reported'] = rep_inv
            crosscheck['blade_s2c_min_sj_reported'] = rep_min
            crosscheck['our_n_sj_lt0'] = n_sj_neg
            crosscheck['our_min_sj'] = min_sj
            if rep_inv is not None and int(rep_inv) != n_sj_neg:
                crosscheck['difference'] = f"counts differ: s2c reports {rep_inv} inverted (sJ<0 per their measure), our audit finds {n_sj_neg} cells with sJ<0 using clean_blocks.scaled_jacobians (corner-tet min)."
            else:
                crosscheck['difference'] = "consistent"
        except Exception as e:
            crosscheck['error'] = str(e)
    else:
        crosscheck['note'] = 'blade_s2c.json not found'
    
    # Write audit VTK with CELL_DATA arrays sJ, volume, inverted, block_id
    out_vtk = Path(args.out_vtk)
    out_vtk.parent.mkdir(parents=True, exist_ok=True)
    try:
        import meshio
        cells = [("hexahedron", H)]
        cell_data = {}
        cell_data["sJ"] = [sj_mins.astype(float)]
        cell_data["volume"] = [vols.astype(float)]
        cell_data["inverted"] = [(sj_mins < 0).astype(int)]
        if bid is not None:
            cell_data["block_id"] = [bid.astype(int)]
        mesh_out = meshio.Mesh(P, cells, cell_data=cell_data)
        mesh_out.write(str(out_vtk))
        print(f"[audit] wrote {out_vtk}")
    except Exception as e:
        print(f"[audit] ERROR writing VTK: {e}")
        raise
    
    # Write audit JSON
    audit = {
        "file": fp,
        "mesh": {
            "n_points": n_points,
            "n_cells": n_cells,
            "n_blocks": n_blocks,
            "block_ids": block_ids,
        },
        "metrics": {
            "n_volume_le0": n_vol_le0,
            "min_volume": min_vol,
            "max_volume": max_vol,
            "n_sJ_lt0": n_sj_neg,
            "n_sJ_lt_0.3": n_sj_lt03,
            "min_sJ": min_sj,
            "max_sJ": max_sj,
            "convention": {
                "volume": "divergence theorem with VTK hex face winding (6 faces, two triangles per face), /6 (tfi._hex_signed_volume / blade_check.signed_volume style)",
                "scaled_jacobian": "clean_blocks.scaled_jacobians: 8 corner tetrahedra, per corner det([u,v,w])/(|u||v||w|) with u,v,w = edge vectors from corner; min over 8 corners"
            }
        },
        "worst_20_cells": worst,
        "crosscheck_blade_s2c": crosscheck,
        "claimed": {"inverted_cells": 215, "min_sJ": -0.9937},
    }
    if n_sj_neg == 215 and abs(min_sj + 0.9937) < 1e-3:
        audit["confirmed_statement"] = "215/-0.9937 is CONFIRMED."
    elif n_sj_neg == 215:
        audit["confirmed_statement"] = f"215 inverted cells (sJ<0) confirmed in count, but min sJ differs: our min_sJ={min_sj:.6g} vs claimed -0.9937 (small difference possible by definition details)."
    else:
        audit["confirmed_statement"] = f"215/-0.9937 is NOT confirmed by clean_blocks.scaled_jacobians; our audit gives n(sJ<0)={n_sj_neg}, min sJ={min_sj:.6g}."
    
    out_json = Path(args.out_json)
    out_json.write_text(json.dumps(audit, indent=2))
    print(f"[audit] wrote {out_json}")
    
    # Make PNG
    out_png = Path(args.out_png)
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle(f'tfi_blade.vtk Inversions - {n_sj_neg} inverted cells (sJ<0), min sJ={min_sj:.4f}', 
                 fontsize=16, fontweight='bold')
    
    # Left panel concept - put histogram in top-left spanning? Or do as specified: left panel histogram, right panel 2x2
    # "left panel a log-y histogram of sJ with lines at 0 and 0.3 and the inverted count annotated; right panel a 2x2 of the four blocks with the most inverted cells"
    plt.close()
    fig = plt.figure(figsize=(16, 12))
    gs = fig.add_gridspec(2, 2)
    ax_hist = fig.add_subplot(gs[:, 0])
    ax_hist.hist(sj_mins, bins=100, range=(-1.0, 1.0), log=True, alpha=0.7, edgecolor='black')
    ax_hist.axvline(0.0, color='red', linestyle='--', linewidth=2, label='sJ=0')
    ax_hist.axvline(0.3, color='orange', linestyle='--', linewidth=2, label='sJ=0.3')
    ax_hist.text(0.02, 0.98, f'Inverted (sJ<0): {n_sj_neg}\nCells with sJ<0.3: {n_sj_lt03}\nMin sJ: {min_sj:.4f}', 
                 transform=ax_hist.transAxes, verticalalignment='top',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax_hist.set_xlabel('Scaled Jacobian', fontsize=12)
    ax_hist.set_ylabel('Count (log scale)', fontsize=12)
    ax_hist.set_title('Histogram of sJ', fontsize=14, fontweight='bold')
    ax_hist.legend()
    ax_hist.grid(True, alpha=0.3)
    
    # Find 4 blocks with most inverted cells
    if bid is not None:
        inv_by_block = defaultdict(int)
        for i in range(n_cells):
            if sj_mins[i] < 0:
                inv_by_block[int(bid[i])] += 1
        sorted_blocks = sorted(inv_by_block.items(), key=lambda x: x[1], reverse=True)
        top4 = []
        for b, cnt in sorted_blocks:
            if len(top4) < 4:
                top4.append(b)
        if len(top4) < 4:
            cnt_by_block = Counter(int(bid[i]) for i in range(n_cells))
            for b, _ in cnt_by_block.most_common():
                if b not in top4:
                    top4.append(b)
                if len(top4) >= 4:
                    break
    else:
        top4 = []
    
    ax_list = [fig.add_subplot(gs[0,1]), fig.add_subplot(gs[1,1])]
    # put 2x2 in right side as 2 rows x 1 column of pairs? Or just do 2x2 on right - gs is 2x2 total, left takes col0 full height, right takes col1 split into 2x2? Wait gs(2,2) means 2 rows, 2 cols - so gs[:,0] is left column full height (1 subplot), gs[0,1] and gs[1,1] are right column split vertically - that's 2 subplots on right, not 2x2. But requirement says "right panel a 2x2 of the four blocks..." - need 2x2 = 4 panels. So better use different layout.
    plt.close()
    
    fig = plt.figure(figsize=(18, 12))
    gs = fig.add_gridspec(2, 3)
    ax_hist = fig.add_subplot(gs[:, 0])  # left panel
    ax_hist.hist(sj_mins, bins=100, range=(-1.0, 1.0), log=True, alpha=0.7, edgecolor='black')
    ax_hist.axvline(0.0, color='red', linestyle='--', linewidth=2, label='sJ=0')
    ax_hist.axvline(0.3, color='orange', linestyle='--', linewidth=2, label='sJ=0.3')
    ax_hist.text(0.02, 0.98, f'Inverted (sJ<0): {n_sj_neg}', 
                 transform=ax_hist.transAxes, verticalalignment='top',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax_hist.set_xlabel('Scaled Jacobian', fontsize=12)
    ax_hist.set_ylabel('Count (log scale)', fontsize=12)
    ax_hist.set_title('Histogram of sJ', fontsize=14, fontweight='bold')
    ax_hist.legend()
    ax_hist.grid(True, alpha=0.3)
    
    block_axes = [(0,1),(0,2),(1,1),(1,2)]
    for idx, b in enumerate(top4[:4]):
        ax = fig.add_subplot(gs[block_axes[idx][0], block_axes[idx][1]])
        if bid is None:
            ax.text(0.5, 0.5, 'No block_id', ha='center', va='center')
            ax.set_title(f'Block {b}')
            continue
        cells_in_block = np.where(bid == b)[0]
        inv_in_block = []
        for c in cells_in_block:
            if sj_mins[c] < 0:
                inv_in_block.append(P[H[c]].mean(axis=0))
        inv_in_block = np.array(inv_in_block) if inv_in_block else np.array([]).reshape(0,3)
        if len(cells_in_block) > 0:
            all_pts_in_block = P[H[cells_in_block]].reshape(-1, 3)
            ax.scatter(all_pts_in_block[:, 0], all_pts_in_block[:, 1], 
                      s=0.05, alpha=0.05, c='gray')
        if len(inv_in_block) > 0:
            ax.scatter(inv_in_block[:, 0], inv_in_block[:, 1], 
                      s=40, alpha=0.9, c='red', marker='o')
        ax.set_title(f'Block {b}: {len(inv_in_block)} inverted', fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.set_aspect('equal', adjustable='datalim')
    
    plt.tight_layout()
    plt.savefig(str(out_png), dpi=150, bbox_inches='tight')
    print(f"[audit] wrote {out_png}")
    
    print("\n=== AUDIT SUMMARY ===")
    print(f"File: {vtk_path}")
    print(f"Hash: {fp['sha256']}")
    print(f"Cells: {n_cells}, Points: {n_points}, Blocks: {n_blocks} ({block_ids})")
    print(f"n(sJ < 0): {n_sj_neg}")
    print(f"min sJ: {min_sj}")
    print(f"n(sJ < 0.3): {n_sj_lt03}")
    print(f"n(volume <= 0): {n_vol_le0}")
    print(f"\n{audit['confirmed_statement']}")


if __name__ == '__main__':
    main()
