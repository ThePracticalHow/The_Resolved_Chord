#!/usr/bin/env python3
"""
FIGURE A — THE BARE ORBIFOLD
Cross-section of S^5/Z_3: three 120° sectors with finite-width fold walls.
The first time the reader sees the geometry.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Wedge, FancyArrowPatch
import matplotlib
import os, sys

matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['mathtext.fontset'] = 'cm'

# ── Palette (matching lotus_emblem.py) ──────────────────────────
BG       = '#fafafa'          # light background for paper figure
BLUE     = '#1a5276'          # sector 0 (χ₀)
BLUE_L   = '#5dade2'
GOLD     = '#b7950b'          # sector 1 (χ₁)
GOLD_L   = '#f4d03f'
ROSE     = '#922b21'          # sector 2 (χ₂)
ROSE_L   = '#ec7063'
WALL_C   = '#2c3e50'          # fold wall dark
TEXT_C   = '#1a1a2e'
GRAY     = '#7f8c8d'

fig, ax = plt.subplots(figsize=(8, 8), facecolor=BG)
ax.set_xlim(-5.5, 5.5)
ax.set_ylim(-5.5, 5.8)
ax.set_aspect('equal')
ax.axis('off')

R = 4.0          # main circle radius
WALL_W = 8       # wall angular half-width in degrees (for visual thickness)

# ── Three sectors as filled wedges ──────────────────────────────
sector_colors = [
    (BLUE,  BLUE_L,  r'Sector 0',  90 + 5),     # χ₀
    (GOLD,  GOLD_L,  r'Sector 1',  210 + 5),     # χ₁
    (ROSE,  ROSE_L,  r'Sector 2',  330 + 5),     # χ₂
]

# Draw three main sectors (120° each, starting from 90°)
for i, (dark, light, label, label_angle_deg) in enumerate(sector_colors):
    start_angle = 90 + i * 120 + WALL_W
    extent = 120 - 2 * WALL_W
    wedge = Wedge((0, 0), R, start_angle, start_angle + extent,
                  facecolor=light, edgecolor='none', alpha=0.35, zorder=2)
    ax.add_patch(wedge)
    # Slightly smaller for inner color gradient
    wedge2 = Wedge((0, 0), R * 0.92, start_angle + 3, start_angle + extent - 3,
                   facecolor=light, edgecolor='none', alpha=0.25, zorder=3)
    ax.add_patch(wedge2)

# ── Fold walls (finite-width bands at 90°, 210°, 330°) ─────────
for i in range(3):
    wall_center = 90 + i * 120  # degrees
    wedge_wall = Wedge((0, 0), R, wall_center - WALL_W, wall_center + WALL_W,
                       facecolor=WALL_C, edgecolor='none', alpha=0.18, zorder=4)
    ax.add_patch(wedge_wall)
    # Dashed center line
    angle_rad = np.radians(wall_center)
    x_end = R * np.cos(angle_rad)
    y_end = R * np.sin(angle_rad)
    ax.plot([0, x_end], [0, y_end], color=WALL_C, linewidth=1.2,
            alpha=0.5, zorder=5, linestyle='--')

# ── Outer circle boundary ──────────────────────────────────────
theta = np.linspace(0, 2*np.pi, 300)
ax.plot(R*np.cos(theta), R*np.sin(theta), color=TEXT_C, linewidth=1.5,
        alpha=0.7, zorder=6)

# ── Central singular point ─────────────────────────────────────
# Cone singularity at the origin
ax.plot(0, 0, 'o', color=TEXT_C, markersize=6, zorder=10)
ax.plot(0, 0, 'o', color='white', markersize=3, zorder=11)
ax.annotate('cone point', (0, 0), (0.6, -0.9),
            fontsize=10, color=GRAY, ha='left',
            arrowprops=dict(arrowstyle='->', color=GRAY, lw=0.8),
            zorder=12)

# ── Sector labels ──────────────────────────────────────────────
label_r = R * 0.55
for i, (dark, light, label, _) in enumerate(sector_colors):
    mid_angle = np.radians(90 + i*120 + 60)  # center of each sector
    lx = label_r * np.cos(mid_angle)
    ly = label_r * np.sin(mid_angle)
    # Sector number with χ character
    chi_labels = [r'$\chi_0$', r'$\chi_1$', r'$\chi_2$']
    ax.text(lx, ly, chi_labels[i], ha='center', va='center',
            fontsize=18, color=dark, fontweight='bold', zorder=8)

# ── ℤ₃ rotation arrow (curving around near the boundary) ──────
# Draw a curved arrow showing the 120° rotation
arrow_r = R + 0.6
arc_angles = np.linspace(np.radians(60), np.radians(-55), 60)
arc_x = arrow_r * np.cos(arc_angles)
arc_y = arrow_r * np.sin(arc_angles)
ax.plot(arc_x, arc_y, color=GRAY, linewidth=1.2, alpha=0.6, zorder=7)
# Arrowhead
ax.annotate('', xy=(arc_x[-1], arc_y[-1]), xytext=(arc_x[-3], arc_y[-3]),
            arrowprops=dict(arrowstyle='->', color=GRAY, lw=1.2), zorder=7)
# Label the rotation
ax.text(arrow_r + 0.3, (arc_y[0]+arc_y[-1])/2 + 0.5,
        r'$z \mapsto \omega z$',
        fontsize=13, color=TEXT_C, ha='left', va='center', zorder=8)
ax.text(arrow_r + 0.3, (arc_y[0]+arc_y[-1])/2 - 0.2,
        r'$\omega = e^{2\pi i/3}$',
        fontsize=11, color=GRAY, ha='left', va='center', zorder=8)

# ── Fold wall label ────────────────────────────────────────────
wall_label_angle = np.radians(90 + 15)
wx = (R + 0.5) * np.cos(wall_label_angle)
wy = (R + 0.5) * np.sin(wall_label_angle)
wall_mid_angle = np.radians(90)
wmx = R * 0.6 * np.cos(wall_mid_angle)
wmy = R * 0.6 * np.sin(wall_mid_angle)
ax.annotate('fold wall', (wmx, wmy), (-2.5, R + 1.0),
            fontsize=10, color=WALL_C, ha='center', fontstyle='italic',
            arrowprops=dict(arrowstyle='->', color=WALL_C, lw=0.8),
            zorder=12)

# ── S⁵ boundary label ─────────────────────────────────────────
bnd_angle = np.radians(15)
bx = R * np.cos(bnd_angle)
by = R * np.sin(bnd_angle)
ax.annotate(r'$S^5$', (bx, by), (R + 1.5, 1.5),
            fontsize=14, color=TEXT_C, ha='center', fontweight='bold',
            arrowprops=dict(arrowstyle='->', color=TEXT_C, lw=0.8),
            zorder=12)

# ── Title ──────────────────────────────────────────────────────
ax.text(0, 5.3, r'$S^5 / \mathbb{Z}_3$',
        ha='center', fontsize=24, color=TEXT_C, fontweight='bold')

# ── Save ───────────────────────────────────────────────────────
out_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
outpath = os.path.join(out_dir, 'bare_orbifold.png')
plt.savefig(outpath, dpi=300, bbox_inches='tight', facecolor=BG)
print(f"Saved {outpath}")
plt.close()
