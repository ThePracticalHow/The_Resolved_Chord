#!/usr/bin/env python3
"""
FIGURE B — THE FIVE INVARIANTS ON THE LOTUS
Same orbifold cross-section as Figure A, annotated with the five spectral
invariants {d₁, λ₁, K, η, p} mapped to specific geometric features.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, Circle, RegularPolygon
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
GREEN    = '#1e8449'    # for K invariant
PURPLE   = '#6c3483'    # for λ₁ ring

fig, ax = plt.subplots(figsize=(10, 9), facecolor=BG)
ax.set_xlim(-7.0, 7.0)
ax.set_ylim(-6.5, 7.0)
ax.set_aspect('equal')
ax.axis('off')

R = 4.0
WALL_W = 8

# ── Three sectors (same as Figure A) ───────────────────────────
sector_data = [
    (BLUE, BLUE_L), (GOLD, GOLD_L), (ROSE, ROSE_L)
]
for i, (dark, light) in enumerate(sector_data):
    start_angle = 90 + i * 120 + WALL_W
    extent = 120 - 2 * WALL_W
    wedge = Wedge((0, 0), R, start_angle, start_angle + extent,
                  facecolor=light, edgecolor='none', alpha=0.30, zorder=2)
    ax.add_patch(wedge)

# ── Fold walls ─────────────────────────────────────────────────
for i in range(3):
    wall_center = 90 + i * 120
    wedge_wall = Wedge((0, 0), R, wall_center - WALL_W, wall_center + WALL_W,
                       facecolor=WALL_C, edgecolor='none', alpha=0.15, zorder=4)
    ax.add_patch(wedge_wall)
    angle_rad = np.radians(wall_center)
    ax.plot([0, R*np.cos(angle_rad)], [0, R*np.sin(angle_rad)],
            color=WALL_C, linewidth=1.0, alpha=0.4, zorder=5, linestyle='--')

# ── Outer circle ───────────────────────────────────────────────
theta = np.linspace(0, 2*np.pi, 300)
ax.plot(R*np.cos(theta), R*np.sin(theta), color=TEXT_C, linewidth=1.3,
        alpha=0.6, zorder=6)

# ══════════════════════════════════════════════════════════════
# INVARIANT 1: p = 3 — the three sectors themselves
# ══════════════════════════════════════════════════════════════
for i, (dark, light) in enumerate(sector_data):
    mid_angle = np.radians(90 + i*120 + 60)
    lx = R * 0.5 * np.cos(mid_angle)
    ly = R * 0.5 * np.sin(mid_angle)
    ax.text(lx, ly, str(i+1), ha='center', va='center',
            fontsize=20, color=dark, fontweight='bold', alpha=0.5, zorder=8)

# p = 3 label
ax.text(R * 0.85, -R * 0.25,
        r'$p = 3$', ha='center', va='center',
        fontsize=14, color=ROSE,
        bbox=dict(boxstyle='round,pad=0.3', facecolor=BG,
                  edgecolor=ROSE, alpha=0.9, linewidth=1.2),
        zorder=15)

# ══════════════════════════════════════════════════════════════
# INVARIANT 2: d₁ = 6 — six ghost mode dots on fold walls
# ══════════════════════════════════════════════════════════════
ghost_positions = []
for i in range(3):
    wall_angle = np.radians(90 + i * 120)
    for frac in [0.35, 0.70]:
        gx = frac * R * np.cos(wall_angle)
        gy = frac * R * np.sin(wall_angle)
        ghost_positions.append((gx, gy))
        ax.plot(gx, gy, 'o', color=BLUE, markersize=10, zorder=12,
                markeredgecolor='white', markeredgewidth=1.5)

# d₁ = 6 label  
ax.text(-R - 1.6, R * 0.6,
        r'$d_1 = 6$', ha='center', va='center',
        fontsize=14, color=BLUE,
        bbox=dict(boxstyle='round,pad=0.3', facecolor=BG,
                  edgecolor=BLUE, alpha=0.9, linewidth=1.2),
        zorder=15)
ax.annotate('', xy=(ghost_positions[0][0]-0.15, ghost_positions[0][1]-0.15),
            xytext=(-R - 0.8, R * 0.5),
            arrowprops=dict(arrowstyle='->', color=BLUE, lw=1.0), zorder=14)
ax.text(-R - 1.6, R * 0.15, r'ghost modes', ha='center',
        fontsize=9, color=GRAY, fontstyle='italic')

# ══════════════════════════════════════════════════════════════
# INVARIANT 3: λ₁ = 5 — first harmonic ring
# ══════════════════════════════════════════════════════════════
lambda_r = R * 0.82
lambda_theta = np.linspace(0, 2*np.pi, 200)
ax.plot(lambda_r*np.cos(lambda_theta), lambda_r*np.sin(lambda_theta),
        color=PURPLE, linewidth=1.8, alpha=0.5, zorder=7, linestyle='-.')

ax.text(R + 1.6, R * 0.6,
        r'$\lambda_1 = 5$', ha='center', va='center',
        fontsize=14, color=PURPLE,
        bbox=dict(boxstyle='round,pad=0.3', facecolor=BG,
                  edgecolor=PURPLE, alpha=0.9, linewidth=1.2),
        zorder=15)
ax.annotate('', xy=(lambda_r*np.cos(np.radians(25)),
                    lambda_r*np.sin(np.radians(25))),
            xytext=(R + 0.8, R * 0.5),
            arrowprops=dict(arrowstyle='->', color=PURPLE, lw=1.0), zorder=14)
ax.text(R + 1.6, R * 0.15, r'first eigenvalue', ha='center',
        fontsize=9, color=GRAY, fontstyle='italic')

# ══════════════════════════════════════════════════════════════
# INVARIANT 4: K = 2/3 — inscribed equilateral triangle (Koide simplex)
# ══════════════════════════════════════════════════════════════
K_r = R * 0.42
tri_angles = [np.radians(90 + i*120) for i in range(3)]
tri_x = [K_r * np.cos(a) for a in tri_angles] + [K_r * np.cos(tri_angles[0])]
tri_y = [K_r * np.sin(a) for a in tri_angles] + [K_r * np.sin(tri_angles[0])]
ax.plot(tri_x, tri_y, color=GREEN, linewidth=2.0, alpha=0.7, zorder=9)
# Vertices
for a in tri_angles:
    ax.plot(K_r*np.cos(a), K_r*np.sin(a), 'D', color=GREEN,
            markersize=6, zorder=10, markeredgecolor='white', markeredgewidth=1)

ax.text(0, -R - 1.3,
        r'$K = 2/3$', ha='center', va='center',
        fontsize=14, color=GREEN,
        bbox=dict(boxstyle='round,pad=0.3', facecolor=BG,
                  edgecolor=GREEN, alpha=0.9, linewidth=1.2),
        zorder=15)
ax.annotate('', xy=(K_r*np.cos(np.radians(-30))*0.9,
                    K_r*np.sin(np.radians(-30))*0.9),
            xytext=(0.8, -R - 1.2),
            arrowprops=dict(arrowstyle='->', color=GREEN, lw=1.0), zorder=14)
ax.text(0, -R - 1.8, r'Koide simplex', ha='center',
        fontsize=9, color=GRAY, fontstyle='italic')

# ══════════════════════════════════════════════════════════════
# INVARIANT 5: η = 2/9 — fold wall width annotation
# ══════════════════════════════════════════════════════════════
# Annotate the fold wall width at the 210° wall
eta_wall_angle = 210
eta_r = R * 0.75
eta_a1 = np.radians(eta_wall_angle - WALL_W)
eta_a2 = np.radians(eta_wall_angle + WALL_W)
# Small double-headed arrow showing width
ax.annotate('', xy=(eta_r*np.cos(eta_a1), eta_r*np.sin(eta_a1)),
            xytext=(eta_r*np.cos(eta_a2), eta_r*np.sin(eta_a2)),
            arrowprops=dict(arrowstyle='<->', color=GOLD, lw=1.5), zorder=14)

ax.text(-R - 1.6, -R * 0.4,
        r'$\eta = 2/9$', ha='center', va='center',
        fontsize=14, color=GOLD,
        bbox=dict(boxstyle='round,pad=0.3', facecolor=BG,
                  edgecolor=GOLD, alpha=0.9, linewidth=1.2),
        zorder=15)
ax.annotate('', xy=(eta_r*np.cos(np.radians(eta_wall_angle)),
                    eta_r*np.sin(np.radians(eta_wall_angle))),
            xytext=(-R - 0.8, -R * 0.4),
            arrowprops=dict(arrowstyle='->', color=GOLD, lw=1.0), zorder=14)
ax.text(-R - 1.6, -R * 0.85, r'fold wall width', ha='center',
        fontsize=9, color=GRAY, fontstyle='italic')

# ── Center point ───────────────────────────────────────────────
ax.plot(0, 0, 'o', color=TEXT_C, markersize=5, zorder=11)

# ── Title ──────────────────────────────────────────────────────
ax.text(0, 6.2, r'The Five Spectral Invariants of $S^5/\mathbb{Z}_3$',
        ha='center', fontsize=18, color=TEXT_C, fontweight='bold')

# ── Save ───────────────────────────────────────────────────────
out_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
outpath = os.path.join(out_dir, 'five_invariants.png')
plt.savefig(outpath, dpi=300, bbox_inches='tight', facecolor=BG)
print(f"Saved {outpath}")
plt.close()
