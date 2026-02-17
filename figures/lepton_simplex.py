#!/usr/bin/env python3
"""
FIGURE C — THE LEPTON SIMPLEX (Koide moment map)
The standard 2-simplex with the three lepton masses as Z₃-orbit points.
Inset ghost overlay showing connection to the orbifold.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Arc, Circle, Wedge
import matplotlib
import os, sys

matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['mathtext.fontset'] = 'cm'

# ── Palette ─────────────────────────────────────────────────────
BG       = '#fafafa'
BLUE     = '#1a5276';  BLUE_L   = '#5dade2'
GOLD     = '#b7950b';  GOLD_L   = '#f4d03f'
ROSE     = '#922b21';  ROSE_L   = '#ec7063'
GREEN    = '#1e8449'
TEXT_C   = '#1a1a2e'
GRAY     = '#7f8c8d'

fig, ax = plt.subplots(figsize=(9, 8.5), facecolor=BG)
ax.set_xlim(-4.5, 4.5)
ax.set_ylim(-3.5, 5.5)
ax.set_aspect('equal')
ax.axis('off')

# ── The simplex triangle (equilateral, vertices at 120° intervals) ──
simplex_R = 3.2
verts = []
vert_labels = [r'$\sqrt{m_e}$', r'$\sqrt{m_\mu}$', r'$\sqrt{m_\tau}$']
vert_colors = [BLUE, GOLD, ROSE]
for i in range(3):
    angle = np.radians(90 + i * 120)
    verts.append((simplex_R * np.cos(angle), simplex_R * np.sin(angle)))

# Draw simplex edges
for i in range(3):
    j = (i + 1) % 3
    ax.plot([verts[i][0], verts[j][0]], [verts[i][1], verts[j][1]],
            color=TEXT_C, linewidth=1.5, alpha=0.3, zorder=2)

# Fill simplex lightly
tri_x = [v[0] for v in verts] + [verts[0][0]]
tri_y = [v[1] for v in verts] + [verts[0][1]]
ax.fill(tri_x[:-1], tri_y[:-1], color=GRAY, alpha=0.05, zorder=1)

# ── Inscribed circle (Koide circle, radius r = √2 scaled) ────
# In the simplex, the inscribed circle has radius = simplex_R / 2
# (for equilateral triangle inscribed circle = side / (2√3))
# But we use it as the Brannen circle
incircle_R = simplex_R * 0.5  # Inscribed circle proportional
inscribed_theta = np.linspace(0, 2*np.pi, 200)
ax.plot(incircle_R * np.cos(inscribed_theta),
        incircle_R * np.sin(inscribed_theta),
        color=GREEN, linewidth=1.8, alpha=0.5, zorder=4, linestyle='-')
ax.text(incircle_R + 0.4, 0.3, r'$r = \sqrt{2}$', fontsize=11,
        color=GREEN, ha='left', zorder=6)

# ── Three mass points on the inscribed circle ──────────────────
# The Koide phase δ = 2π/3 + 2/9 ≈ 2.342 radians
delta = 2*np.pi/3 + 2/9   # The exact Koide phase

# Brannen parametrization: sqrt(m_k) = μ(1 + √2 cos(δ + 2πk/3))
# We place points on the circle at angles δ + 2πk/3
mass_angles = [delta + 2*np.pi*k/3 for k in range(3)]
mass_labels = [r'$m_e$', r'$m_\mu$', r'$m_\tau$']
mass_colors = [BLUE, GOLD, ROSE]

for k in range(3):
    mx = incircle_R * np.cos(mass_angles[k])
    my = incircle_R * np.sin(mass_angles[k])
    ax.plot(mx, my, 'o', color=mass_colors[k], markersize=12, zorder=8,
            markeredgecolor='white', markeredgewidth=2)
    # Label offset
    label_r = incircle_R + 0.6
    lx = label_r * np.cos(mass_angles[k])
    ly = label_r * np.sin(mass_angles[k])
    ax.text(lx, ly, mass_labels[k], ha='center', va='center',
            fontsize=14, color=mass_colors[k], fontweight='bold', zorder=9)

# ── ℤ₃ rotation arrows between mass points ────────────────────
for k in range(3):
    k_next = (k + 1) % 3
    # Arc between consecutive mass points
    a1 = mass_angles[k]
    a2 = mass_angles[k_next]
    # Draw small curved arrow
    mid_a = (a1 + a2) / 2
    arc_r = incircle_R * 0.8
    n_pts = 30
    a_range = np.linspace(a1 + 0.15, a2 - 0.15, n_pts)
    arc_x = arc_r * np.cos(a_range)
    arc_y = arc_r * np.sin(a_range)
    ax.plot(arc_x, arc_y, color=GRAY, linewidth=0.8, alpha=0.5, zorder=5)
    ax.annotate('', xy=(arc_x[-1], arc_y[-1]),
                xytext=(arc_x[-3], arc_y[-3]),
                arrowprops=dict(arrowstyle='->', color=GRAY, lw=0.8), zorder=5)

# ── Phase angle annotation ─────────────────────────────────────
# Show δ as the angle from a reference direction to the first mass point
ref_angle = 0  # Reference: positive x-axis
# Draw reference line (thin, dashed)
ax.plot([0, incircle_R*1.1*np.cos(ref_angle)],
        [0, incircle_R*1.1*np.sin(ref_angle)],
        color=GRAY, linewidth=0.8, linestyle=':', alpha=0.4, zorder=3)
# Arc from reference to first mass point
arc_delta = np.linspace(ref_angle, mass_angles[0], 40)
arc_delta_r = incircle_R * 0.35
ax.plot(arc_delta_r * np.cos(arc_delta), arc_delta_r * np.sin(arc_delta),
        color=TEXT_C, linewidth=1.2, alpha=0.6, zorder=5)
# δ label
mid_delta = (ref_angle + mass_angles[0]) / 2
ax.text(arc_delta_r * 1.5 * np.cos(mid_delta),
        arc_delta_r * 1.5 * np.sin(mid_delta),
        r'$\delta = \frac{2\pi}{3} + \frac{2}{9}$',
        ha='center', va='center', fontsize=12, color=TEXT_C, zorder=6)

# ── Origin ─────────────────────────────────────────────────────
ax.plot(0, 0, '+', color=TEXT_C, markersize=8, markeredgewidth=1.5, zorder=7)

# ── Small orbifold inset (bottom-right) ────────────────────────
inset_cx, inset_cy = 3.0, -2.5
inset_r = 1.0
# Draw mini-orbifold
for i, (dark, light) in enumerate([(BLUE, BLUE_L), (GOLD, GOLD_L), (ROSE, ROSE_L)]):
    start = 90 + i * 120
    w = Wedge((inset_cx, inset_cy), inset_r, start + 4, start + 116,
              facecolor=light, edgecolor='none', alpha=0.3, zorder=2)
    ax.add_patch(w)
# Border
inset_th = np.linspace(0, 2*np.pi, 100)
ax.plot(inset_cx + inset_r*np.cos(inset_th),
        inset_cy + inset_r*np.sin(inset_th),
        color=TEXT_C, linewidth=1.0, alpha=0.4, zorder=3)
# Fold walls in inset
for i in range(3):
    a = np.radians(90 + i*120)
    ax.plot([inset_cx, inset_cx + inset_r*np.cos(a)],
            [inset_cy, inset_cy + inset_r*np.sin(a)],
            color=TEXT_C, linewidth=0.6, alpha=0.3, linestyle='--', zorder=4)
# Mini simplex inside inset
kr = inset_r * 0.45
for i in range(3):
    a1 = np.radians(90 + i*120)
    a2 = np.radians(90 + (i+1)*120)
    ax.plot([inset_cx + kr*np.cos(a1), inset_cx + kr*np.cos(a2)],
            [inset_cy + kr*np.sin(a1), inset_cy + kr*np.sin(a2)],
            color=GREEN, linewidth=1.2, alpha=0.6, zorder=5)
# Label
ax.text(inset_cx, inset_cy - inset_r - 0.4,
        r'$S^5/\mathbb{Z}_3$', ha='center', fontsize=10,
        color=GRAY, fontstyle='italic')

# Connecting dashed line from main simplex to inset
ax.plot([0.5, inset_cx - inset_r + 0.2], [-1.0, inset_cy + inset_r - 0.3],
        color=GRAY, linewidth=0.8, linestyle=':', alpha=0.4, zorder=1)

# ── Simplex vertex labels ─────────────────────────────────────
for i in range(3):
    offset = 0.5
    lx = verts[i][0] + offset * np.cos(np.radians(90 + i*120))
    ly = verts[i][1] + offset * np.sin(np.radians(90 + i*120))
    ax.text(lx, ly, vert_labels[i], ha='center', va='center',
            fontsize=11, color=GRAY, alpha=0.7, zorder=6)

# ── Title ──────────────────────────────────────────────────────
ax.text(0, 5.0, r'The Koide Simplex on $S^5/\mathbb{Z}_3$',
        ha='center', fontsize=18, color=TEXT_C, fontweight='bold')

# ── Save ───────────────────────────────────────────────────────
out_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
outpath = os.path.join(out_dir, 'lepton_simplex.png')
plt.savefig(outpath, dpi=300, bbox_inches='tight', facecolor=BG)
print(f"Saved {outpath}")
plt.close()
