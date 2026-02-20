#!/usr/bin/env python3
"""
THE LOTUS EMBLEM: A three-petaled lotus for the last page.
Three petals = three Z_3 sectors. Colors = the spectral data.
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.path import Path
import matplotlib

matplotlib.rcParams['font.family'] = 'serif'

fig, ax = plt.subplots(figsize=(10, 12), facecolor='#0a0a1a')
ax.set_xlim(-5, 5)
ax.set_ylim(-4, 6.5)
ax.set_aspect('equal')
ax.axis('off')

# Petal shape: a teardrop using Bezier curves
def draw_petal(ax, angle_deg, color1, color2, alpha_val=0.85):
    """Draw a single lotus petal rotated by angle_deg."""
    angle = np.radians(angle_deg)
    
    # Petal outline points (teardrop shape)
    t = np.linspace(0, 2*np.pi, 200)
    # Cardioid-like shape for petal
    r = 2.2 * (1 - 0.3*np.cos(t)) * np.abs(np.sin(t/2))**0.8
    x = r * np.sin(t) * 0.6
    y = r * np.cos(t) * 1.0 + 0.5
    
    # Rotate
    xr = x*np.cos(angle) - y*np.sin(angle)
    yr = x*np.sin(angle) + y*np.cos(angle)
    
    # Gradient fill using multiple layers
    for i, a in enumerate(np.linspace(0.1, alpha_val, 8)):
        scale = 1.0 - i*0.02
        ax.fill(xr*scale, yr*scale, color=color1, alpha=a*0.15, zorder=2+i)
    
    # Outer glow
    ax.fill(xr*1.02, yr*1.02, color=color2, alpha=0.1, zorder=1)
    
    # Edge
    ax.plot(xr, yr, color=color2, linewidth=1.0, alpha=0.6, zorder=10)

# Three petals at 120-degree intervals (Z_3 symmetry!)
# Colors inspired by spectral data:
# Petal 1 (chi_0): electron sector — cool blue
# Petal 2 (chi_1): up-quark sector — warm gold  
# Petal 3 (chi_2): down-quark sector — deep rose

draw_petal(ax, 90, '#1a5276', '#5dade2', 0.9)      # chi_0: blue (top)
draw_petal(ax, 90+120, '#b7950b', '#f4d03f', 0.85)  # chi_1: gold (bottom-left)
draw_petal(ax, 90+240, '#922b21', '#ec7063', 0.85)  # chi_2: rose (bottom-right)

# Center glow (the bud — Z_3 fixed point S^3)
for r_glow in np.linspace(0.8, 0.05, 20):
    circle = plt.Circle((0, 0), r_glow, color='white', alpha=0.04, zorder=20)
    ax.add_patch(circle)
center = plt.Circle((0, 0), 0.15, color='white', alpha=0.9, zorder=25)
ax.add_patch(center)

# Fold walls (three radial lines at 0, 120, 240 degrees)
for angle_deg in [90, 210, 330]:
    angle = np.radians(angle_deg)
    x_end = 2.5 * np.cos(angle)
    y_end = 2.5 * np.sin(angle)
    ax.plot([0, x_end], [0, y_end], color='white', linewidth=0.5, 
            alpha=0.3, zorder=15, linestyle='--')

# The five spectral invariants around the lotus
specs = [
    (0, 4.0, r'$d_1 = 6$', '#5dade2'),
    (-3.2, -2.0, r'$\lambda_1 = 5$', '#f4d03f'),
    (3.2, -2.0, r'$K = 2/3$', '#ec7063'),
    (-2.5, 2.8, r'$\eta = 2/9$', '#a9cce3'),
    (2.5, 2.8, r'$p = 3$', '#f5cba7'),
]

for x, y, text, color in specs:
    ax.text(x, y, text, ha='center', va='center', fontsize=14,
            color=color, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#0a0a1a', 
                      edgecolor=color, alpha=0.8, linewidth=0.8))

# Title
ax.text(0, 5.8, r'$S^5 / \mathbb{Z}_3$',
        ha='center', fontsize=28, color='white', fontweight='bold')
ax.text(0, 5.2, 'The Spectral Geometry of Everything',
        ha='center', fontsize=14, color='#aab7c4', fontstyle='italic')

# Bottom text
ax.text(0, -3.0, 'One manifold.  Five invariants.  Zero free parameters.',
        ha='center', fontsize=12, color='#7f8c8d')
ax.text(0, -3.5, '87 predictions.  All Theorem.  All verified.',
        ha='center', fontsize=11, color='#5d6d7e')

out_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
plt.savefig(os.path.join(out_dir, 'lotus_emblem.png'), dpi=250, 
            bbox_inches='tight', facecolor='#0a0a1a')
print(f"Saved {os.path.join(out_dir, 'lotus_emblem.png')}")
plt.close()
