#!/usr/bin/env python3
"""
THE MONEY SHOT: All predictions vs measurement on one plot.
Pulls live data from lotus.Universe() — never goes stale.
"""
import sys
import os

# Add parent so we can import lotus
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.lines import Line2D

matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.size'] = 10
matplotlib.rcParams['axes.linewidth'] = 0.8

# ── Pull live predictions from LOTUS ──
import lotus
u = lotus.Universe()
raw = u.predictions

# Compute fractional errors; skip exact matches and zero-measured
names = []
errors = []
for name, vals in raw.items():
    pred = vals['predicted']
    meas = vals['measured']
    if meas == 0 or pred == meas:
        continue
    err = (pred - meas) / meas
    # Skip structural identities (errors < 1e-12 are algebraic tautologies)
    if abs(err) < 1e-12:
        continue
    names.append(name)
    errors.append(err)

# Sort by absolute error (smallest at bottom → tightest predictions are most visible)
idx = np.argsort([abs(e) for e in errors])
names = [names[i] for i in idx]
errors = [errors[i] for i in idx]

# ── Figure ──
n_items = len(names)
fig_h = max(10, n_items * 0.30)
fig, ax = plt.subplots(figsize=(10, fig_h))

y_pos = np.arange(len(names))

# Color by error magnitude
colors = []
for err in errors:
    ae = abs(err) * 100
    if ae < 0.1:
        colors.append('#1a5276')   # dark blue — sub-0.1%
    elif ae < 1.0:
        colors.append('#2e86c1')   # medium blue — sub-1%
    else:
        colors.append('#c0392b')   # red — >1%

for i, (name, err, c) in enumerate(zip(names, errors, colors)):
    ax.scatter(err * 100, i, c=c, marker='D', s=50, zorder=3,
               edgecolors='white', linewidth=0.5)

# Zero line
ax.axvline(x=0, color='black', linewidth=0.5, zorder=1)

# Shade bands
ax.axvspan(-0.1, 0.1, alpha=0.08, color='green', zorder=0)
ax.axvspan(-1, -0.1, alpha=0.04, color='#f0c040', zorder=0)
ax.axvspan(0.1, 1, alpha=0.04, color='#f0c040', zorder=0)

ax.set_yticks(y_pos)
ax.set_yticklabels(names, fontsize=8)
ax.set_xlabel('Fractional error: (predicted − measured) / measured  (%)', fontsize=11)
ax.set_title(f'The Resolved Chord: All {len(raw)} Predictions vs Measurement\n'
             'One manifold.  Five spectral invariants.  Zero free parameters.',
             fontsize=13, fontweight='bold')

ax.set_xlim(-5, 5)
ax.set_ylim(-0.5, len(names) - 0.5)

# Legend
n_sub01 = sum(1 for e in errors if abs(e) * 100 < 0.1)
n_sub1 = sum(1 for e in errors if 0.1 <= abs(e) * 100 < 1.0)
n_over1 = sum(1 for e in errors if abs(e) * 100 >= 1.0)
n_exact = len(raw) - len(names)  # identities with zero error

legend_elements = [
    Line2D([0], [0], marker='D', color='w', markerfacecolor='#1a5276',
           markersize=8, label=f'< 0.1% error ({n_sub01})'),
    Line2D([0], [0], marker='D', color='w', markerfacecolor='#2e86c1',
           markersize=8, label=f'< 1% error ({n_sub1})'),
    Line2D([0], [0], marker='D', color='w', markerfacecolor='#c0392b',
           markersize=8, label=f'> 1% error ({n_over1})'),
    Line2D([0], [0], marker='s', color='w', markerfacecolor='#888',
           markersize=8, label=f'Exact identities ({n_exact})'),
]
ax.legend(handles=legend_elements, loc='upper left', fontsize=9,
          framealpha=0.95)

# Annotation
ax.text(0.02, 0.98, 'Green band: < 0.1% error\nYellow band: < 1% error',
        transform=ax.transAxes, fontsize=8, va='top', ha='left',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.9, pad=0.3))

plt.tight_layout(pad=1.2)
out_dir = os.path.dirname(os.path.abspath(__file__))
out_path = os.path.join(out_dir, 'prediction_scatter.png')
plt.savefig(out_path, dpi=200, bbox_inches='tight')
print(f"Saved {out_path}")
plt.close()
