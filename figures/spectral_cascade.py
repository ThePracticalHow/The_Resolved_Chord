#!/usr/bin/env python3
"""
FIGURE D — THE SPECTRAL CASCADE ON THE LOTUS
Four concentric rings on the orbifold, each representing a derivation level:
Level 0 (mₑ) → Level 1 (mₚ) → Level 2 (α) → Level 3 (v, mH).
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, FancyArrowPatch, Circle
import matplotlib
import os, sys

matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['mathtext.fontset'] = 'cm'

# ── Palette ─────────────────────────────────────────────────────
BG       = '#fafafa'
BLUE     = '#1a5276';  BLUE_L   = '#5dade2'
GOLD     = '#b7950b';  GOLD_L   = '#f4d03f'
ROSE     = '#922b21';  ROSE_L   = '#ec7063'
WALL_C   = '#2c3e50'
TEXT_C   = '#1a1a2e'
GRAY     = '#7f8c8d'

# Level colors (from innermost to outermost)
LV0_C = '#2c3e50'     # Core: dark (unit scale)
LV1_C = '#c0392b'     # Ring 1: QCD red
LV2_C = '#2980b9'     # Ring 2: EM blue
LV3_C = '#27ae60'     # Ring 3: EW green

fig, ax = plt.subplots(figsize=(9, 9), facecolor=BG)
ax.set_xlim(-6.5, 6.5)
ax.set_ylim(-6.0, 7.2)
ax.set_aspect('equal')
ax.axis('off')

R = 4.5
WALL_W = 6

# ── Ring radii for the four levels ─────────────────────────────
r0 = 0.35  # Level 0: mₑ (core dot)
r1 = 1.6   # Level 1: mₚ ring
r2 = 2.8   # Level 2: 1/α ring
r3 = R     # Level 3: v, mH ring (outer boundary)

# ── Draw four concentric annular zones ─────────────────────────
# Level 3 (outermost): EW scale
level3_theta = np.linspace(0, 2*np.pi, 300)
ax.fill(r3*np.cos(level3_theta), r3*np.sin(level3_theta),
        color=LV3_C, alpha=0.08, zorder=1)
ax.plot(r3*np.cos(level3_theta), r3*np.sin(level3_theta),
        color=LV3_C, linewidth=2.0, alpha=0.5, zorder=6)

# Level 2: EM scale
ax.fill(r2*np.cos(level3_theta), r2*np.sin(level3_theta),
        color=LV2_C, alpha=0.10, zorder=2)
ax.plot(r2*np.cos(level3_theta), r2*np.sin(level3_theta),
        color=LV2_C, linewidth=1.8, alpha=0.5, zorder=6)

# Level 1: QCD scale
ax.fill(r1*np.cos(level3_theta), r1*np.sin(level3_theta),
        color=LV1_C, alpha=0.12, zorder=3)
ax.plot(r1*np.cos(level3_theta), r1*np.sin(level3_theta),
        color=LV1_C, linewidth=1.8, alpha=0.5, zorder=6)

# Level 0: core
core = Circle((0, 0), r0, facecolor=LV0_C, edgecolor='white',
              linewidth=2, alpha=0.9, zorder=8)
ax.add_patch(core)

# ── Three sector boundaries (fold walls) visible across all rings
for i in range(3):
    wall_angle = np.radians(90 + i * 120)
    ax.plot([0, r3*np.cos(wall_angle)], [0, r3*np.sin(wall_angle)],
            color=WALL_C, linewidth=0.8, alpha=0.25, zorder=5, linestyle='--')
    # Small wedge shading for fold wall
    wc = 90 + i * 120
    w = Wedge((0, 0), r3, wc - WALL_W, wc + WALL_W,
              facecolor=WALL_C, edgecolor='none', alpha=0.06, zorder=4)
    ax.add_patch(w)

# ── Sector colors (very faint, visible in the outer ring) ──────
sector_colors_light = [BLUE_L, GOLD_L, ROSE_L]
for i, sc in enumerate(sector_colors_light):
    start = 90 + i * 120 + WALL_W
    extent = 120 - 2 * WALL_W
    w = Wedge((0, 0), r3, start, start + extent,
              facecolor=sc, edgecolor='none', alpha=0.10, zorder=1)
    ax.add_patch(w)

# ── Level labels (right side, with connecting arrows) ──────────
label_x = R + 1.5

# Level 0
ax.annotate(r'$m_e = 1$ (unit)', xy=(r0, 0), xytext=(label_x, 0.4),
            fontsize=12, color=LV0_C, fontweight='bold',
            arrowprops=dict(arrowstyle='->', color=LV0_C, lw=1.0),
            zorder=12)
ax.text(label_x, -0.1, r'Level 0: Koide ground state',
        fontsize=9, color=GRAY, fontstyle='italic')

# Level 1
ax.annotate(r'$m_p = 6\pi^5 \, m_e$', xy=(r1, 0.3), xytext=(label_x, 2.0),
            fontsize=12, color=LV1_C, fontweight='bold',
            arrowprops=dict(arrowstyle='->', color=LV1_C, lw=1.0),
            zorder=12)
ax.text(label_x, 1.5, r'Level 1: ghost fold energy',
        fontsize=9, color=GRAY, fontstyle='italic')

# Level 2
ax.annotate(r'$1/\alpha = 137.038$', xy=(r2, 0.3), xytext=(label_x, 3.6),
            fontsize=12, color=LV2_C, fontweight='bold',
            arrowprops=dict(arrowstyle='->', color=LV2_C, lw=1.0),
            zorder=12)
ax.text(label_x, 3.1, r'Level 2: APS lag ($G/p = 10/27$)',
        fontsize=9, color=GRAY, fontstyle='italic')

# Level 3
ax.annotate(r'$v,\, m_H$', xy=(r3 - 0.1, 0.3), xytext=(label_x, 5.2),
            fontsize=12, color=LV3_C, fontweight='bold',
            arrowprops=dict(arrowstyle='->', color=LV3_C, lw=1.0),
            zorder=12)
ax.text(label_x, 4.7, r'Level 3: EM budget ($2/\alpha - 35/3$)',
        fontsize=9, color=GRAY, fontstyle='italic')

# ── Cascade arrows (left side, showing derivation flow) ────────
arrow_x = -R - 0.7
arrow_levels = [(-0.2, 1.5, LV0_C), (1.5, 2.8, LV1_C),
                (2.8, 4.2, LV2_C)]
for y1, y2, c in arrow_levels:
    ax.annotate('', xy=(arrow_x, y2), xytext=(arrow_x, y1),
                arrowprops=dict(arrowstyle='->', color=c, lw=2.0,
                                connectionstyle='arc3,rad=0.0'),
                zorder=10)

# Cascade text labels (left)
ax.text(arrow_x - 0.2, 0.65, r'$\times d_1\pi^5$',
        fontsize=10, color=LV1_C, ha='right', rotation=90)
ax.text(arrow_x - 0.2, 2.15, r'$+ G/p$',
        fontsize=10, color=LV2_C, ha='right', rotation=90)
ax.text(arrow_x - 0.2, 3.5, r'$\times 2/\alpha$',
        fontsize=10, color=LV3_C, ha='right', rotation=90)

# ── "Cascade" label ────────────────────────────────────────────
ax.text(arrow_x - 1.2, 2.2, 'Spectral\nCascade',
        fontsize=12, color=TEXT_C, ha='center', va='center',
        fontweight='bold', rotation=90)

# ── Level numbers on the rings ─────────────────────────────────
ring_label_angle = np.radians(240)
for r, level_num, c in [(r1, '1', LV1_C), (r2, '2', LV2_C), (r3, '3', LV3_C)]:
    lx = r * 0.85 * np.cos(ring_label_angle)
    ly = r * 0.85 * np.sin(ring_label_angle)
    ax.text(lx, ly, level_num, ha='center', va='center',
            fontsize=16, color=c, fontweight='bold', alpha=0.4, zorder=7)

# ── Title ──────────────────────────────────────────────────────
ax.text(0, 6.6, r'The Spectral Cascade on $S^5/\mathbb{Z}_3$',
        ha='center', fontsize=18, color=TEXT_C, fontweight='bold')

# ── Save ───────────────────────────────────────────────────────
out_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
outpath = os.path.join(out_dir, 'spectral_cascade.png')
plt.savefig(outpath, dpi=300, bbox_inches='tight', facecolor=BG)
print(f"Saved {outpath}")
plt.close()
