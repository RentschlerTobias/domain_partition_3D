#!/usr/bin/env python3
"""
Plotte die Hexa-Block-Quad-Seiten auf Hub und Shroud abgewickelt.
Hintergrund: triplot des unwrapped Hub/STL (wie stage1 scripts).
Darüber: die dtOO Block-Quad-Grenzen als farbige Linien.
"""

import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import meshio
from pathlib import Path

# =============================================================================
# 1. Pfade
# =============================================================================
BASE = Path("/home/t1dde/Duty/projects/domain_partition/domain_partition_3D")
hub_stl = BASE / "T1_9_hub_raw.stl"
shroud_stl = BASE / "T1_9_shroud_raw.stl"
blocks_json = BASE / "block_instances_28.json"
out_png = BASE / "hexa_blocks_hub_shroud_unwrapped.png"

# =============================================================================
# 2. Hub STL laden und in (s, t) umrechnen
# =============================================================================
hub_mesh = meshio.read(str(hub_stl))
points_3d = hub_mesh.points
tris = hub_mesh.cells_dict.get("triangle", hub_mesh.cells_dict.get("triangle", []))

r_hub = 0.5
r_shroud = 1.9

theta = np.arctan2(points_3d[:, 1], points_3d[:, 0])
s_hub = r_hub * theta
t_hub = points_3d[:, 2]

# =============================================================================
# 3. dtOO Bloecke laden
# =============================================================================
with open(blocks_json, "r") as f:
    blocks_data = json.load(f)

instances = blocks_data["instances"]

# =============================================================================
# 4. Funktion: projiziere physikalische [x, y, z] auf (s, t)
# =============================================================================
def xyz_to_st(xyz, r):
    x, y, z = xyz
    theta = np.arctan2(y, x)
    return r * theta, z

# =============================================================================
# 5. Plot erstellen
# =============================================================================
fig, axes = plt.subplots(1, 2, figsize=(18, 8))

# --- Linker Plot: Hub (r = 0.5) ---
ax = axes[0]
ax.triplot(s_hub, t_hub, tris, lw=0.08, color="0.94", zorder=0)

# Block-Quad-Seiten auf Hub zeichnen
for inst in instances:
    corners = inst["corners"]
    
    # Finde die 4 Ecken mit z=0 (Hub-Seite)
    hub_corners = []
    for i, c in enumerate(corners):
        if abs(c[2]) < 0.01:  # z gleich 0
            hub_corners.append((i, c))
    
    if len(hub_corners) == 4:
        hub_indices = [c[0] for c in hub_corners]
        hub_indices_sorted = sorted(hub_indices)
        
        st_points = []
        for idx in hub_indices_sorted:
            s, t = xyz_to_st(corners[idx], r_hub)
            st_points.append([s, t])
        
        st_points = np.array(st_points)
        
        # Quad zeichnen (geschlossener Linienzug)
        ring = np.vstack([st_points, st_points[0]])
        ax.plot(ring[:, 0], ring[:, 1], 'C0-', lw=1.2, zorder=2)
        ax.fill(st_points[:, 0], st_points[:, 1], alpha=0.15, color="C0", zorder=1)

# Ecken markieren
for inst in instances:
    corners = inst["corners"]
    for i, c in enumerate(corners):
        if abs(c[2]) < 0.01:
            s, t = xyz_to_st(c, r_hub)
            ax.scatter(s, t, c="darkblue", s=30, zorder=3, edgecolors="white", linewidths=0.5)

ax.set_xlabel("s = r mal theta (circumferential)", fontsize=11)
ax.set_ylabel("t = z (axial)", fontsize=11)
ax.set_title(f"Hub (r gleich 0.5): {len(instances)} Hexa-Block-Quads", fontsize=12)
ax.set_aspect("equal")
ax.grid(True, alpha=0.3)

# --- Rechter Plot: Shroud (r = 1.9) ---
ax = axes[1]

# Shroud STL fuer Hintergrund
try:
    shroud_mesh = meshio.read(str(shroud_stl))
    points_shroud = shroud_mesh.points
    tris_shroud = shroud_mesh.cells_dict.get("triangle", [])
    theta_s = np.arctan2(points_shroud[:, 1], points_shroud[:, 0])
    s_shroud = r_shroud * theta_s
    t_shroud = points_shroud[:, 2]
    ax.triplot(s_shroud, t_shroud, tris_shroud, lw=0.08, color="0.94", zorder=0)
except Exception as e:
    print(f"Shroud STL nicht gefunden: {e}")
    ax.triplot(s_hub, t_hub, tris, lw=0.08, color="0.94", zorder=0)

# Block-Quad-Seiten auf Shroud zeichnen
for inst in instances:
    corners = inst["corners"]
    
    # Finde die 4 Ecken mit z=2.5 (Shroud-Seite)
    shroud_corners = []
    for i, c in enumerate(corners):
        if abs(c[2] - 2.5) < 0.01:  # z gleich 2.5
            shroud_corners.append((i, c))
    
    if len(shroud_corners) == 4:
        shroud_indices = [c[0] for c in shroud_corners]
        shroud_indices_sorted = sorted(shroud_indices)
        
        st_points = []
        for idx in shroud_indices_sorted:
            s, t = xyz_to_st(corners[idx], r_shroud)
            st_points.append([s, t])
        
        st_points = np.array(st_points)
        ring = np.vstack([st_points, st_points[0]])
        ax.plot(ring[:, 0], ring[:, 1], 'C1-', lw=1.2, zorder=2)
        ax.fill(st_points[:, 0], st_points[:, 1], alpha=0.15, color="C1", zorder=1)

# Ecken markieren
for inst in instances:
    corners = inst["corners"]
    for i, c in enumerate(corners):
        if abs(c[2] - 2.5) < 0.01:
            s, t = xyz_to_st(c, r_shroud)
            ax.scatter(s, t, c="darkorange", s=30, zorder=3, edgecolors="white", linewidths=0.5)

ax.set_xlabel("s = r mal theta (circumferential)", fontsize=11)
ax.set_ylabel("t = z (axial)", fontsize=11)
ax.set_title(f"Shroud (r gleich 1.9): {len(instances)} Hexa-Block-Quads", fontsize=12)
ax.set_aspect("equal")
ax.grid(True, alpha=0.3)

plt.tight_layout()
fig.savefig(str(out_png), dpi=150, bbox_inches="tight")
plt.close(fig)

print(f"Gespeichert: {out_png}")
