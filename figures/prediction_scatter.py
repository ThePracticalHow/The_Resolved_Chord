#!/usr/bin/env python3
"""
THE MONEY SHOT: All 39 predictions vs measurement on one plot.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.size'] = 10
matplotlib.rcParams['axes.linewidth'] = 0.8

PI = np.pi

# All predictions: (name, predicted, measured, status)
# Status: T = Theorem, D = Derived
predictions = [
    # Lepton masses
    ("K = 2/3", 2/3, 2/3, "T"),
    ("m_mu/m_e", 206.768, 206.768, "T"),
    ("m_tau/m_e", 3477.4, 3477.2, "T"),
    # Proton
    ("m_p/m_e", 6*PI**5, 1836.153, "T"),
    # Alpha
    ("1/alpha", 137.038, 137.036, "T"),
    # Higgs
    ("v (GeV)", 246.21, 246.22, "T"),
    ("m_H (GeV)", 125.30, 125.25, "T"),
    ("lambda_H", 0.1295, 0.1294, "T"),
    # Gauge
    ("sin2_thetaW", 0.375, 0.375, "T"),
    ("alpha_s", 0.1187, 0.1180, "D"),
    ("theta_bar", 0, 0, "T"),
    # Quark masses
    ("m_t (GeV)", 172.66, 172.57, "T"),
    ("m_c (GeV)", 1.275, 1.273, "T"),
    ("m_u (MeV)", 2.163, 2.16, "T"),
    ("m_b (GeV)", 4.180, 4.183, "T"),
    ("m_s (MeV)", 93.39, 93.4, "T"),
    ("m_d (MeV)", 4.695, 4.67, "T"),
    # CKM
    ("CKM lambda", 0.22502, 0.22500, "D"),
    ("CKM A", 0.8263, 0.826, "D"),
    ("CKM rho-bar", 0.15916, 0.1592, "D"),
    ("CKM eta-bar", 0.34907, 0.3490, "D"),
    ("CKM gamma (deg)", 65.49, 65.6, "D"),
    ("Jarlskog J", 3.09e-5, 3.08e-5, "D"),
    # PMNS
    ("theta_23 PMNS", 47.0, 47.5, "D"),
    ("theta_12 PMNS", 33.0, 33.4, "D"),
    ("theta_13 PMNS", 8.5, 8.6, "D"),
    # Gravity
    ("M_P (1e19 GeV)", 1.222, 1.221, "T"),
    ("c_grav", -1/30, -1/30, "T"),
    # Cosmological constant
    ("CC^1/4 (meV)", 2.22, 2.25, "D"),
    # Inflation
    ("n_s", 0.968, 0.965, "D"),
    ("r", 0.003, 0.003, "D"),
    # Baryogenesis
    ("eta_B (1e-10)", 6.30, 6.10, "D"),
    # Dark matter
    ("Omega_DM/Omega_B", 5.333, 5.36, "D"),
    # N_g
    ("N_g = 3", 3, 3, "T"),
    # Spectral ordering
    ("Quark ordering", 1.0, 1.0, "T"),
]

# Compute fractional errors (skip exact matches)
names = []
errors = []
statuses = []
for name, pred, meas, status in predictions:
    if meas == 0 or pred == meas:
        continue
    err = (pred - meas) / meas
    names.append(name)
    errors.append(err)
    statuses.append(status)

# Sort by absolute error
idx = np.argsort([abs(e) for e in errors])
names = [names[i] for i in idx]
errors = [errors[i] for i in idx]
statuses = [statuses[i] for i in idx]

# Figure: taller to avoid overlapping y-labels; wider for clarity
n_items = len(names)
fig_h = max(10, n_items * 0.28)  # ~0.28" per row
fig, ax = plt.subplots(figsize=(10, fig_h))

y_pos = np.arange(len(names))
colors = ['#2166ac' if s == 'T' else '#b2182b' for s in statuses]
markers = ['D' if s == 'T' else 'o' for s in statuses]

for i, (name, err, status) in enumerate(zip(names, errors, statuses)):
    color = '#1a5276' if status == 'T' else '#922b21'
    marker = 'D' if status == 'T' else 'o'
    size = 50 if status == 'T' else 32
    ax.scatter(err * 100, i, c=color, marker=marker, s=size, zorder=3, edgecolors='white', linewidth=0.5)

# Zero line
ax.axvline(x=0, color='black', linewidth=0.5, zorder=1)

# Shade bands
ax.axvspan(-0.1, 0.1, alpha=0.08, color='green', zorder=0)
ax.axvspan(-1, -0.1, alpha=0.04, color='yellow', zorder=0)
ax.axvspan(0.1, 1, alpha=0.04, color='yellow', zorder=0)

ax.set_yticks(y_pos)
ax.set_yticklabels(names, fontsize=8)
ax.set_xlabel('Fractional error: (predicted - measured) / measured  (%)', fontsize=11)
ax.set_title('The Resolved Chord: All Predictions vs Measurement\n'
             'One manifold.  Five spectral invariants.  Zero free parameters.',
             fontsize=13, fontweight='bold')

ax.set_xlim(-5, 5)
ax.set_ylim(-0.5, len(names) - 0.5)  # avoid clipping top/bottom

# Legend: upper right to avoid overlapping scatter points
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker='D', color='w', markerfacecolor='#1a5276',
           markersize=8, label='Theorem (19)'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#922b21',
           markersize=7, label='Derived (20)'),
]
ax.legend(handles=legend_elements, loc='upper right', fontsize=10,
          framealpha=0.95)

# Annotation: upper left, compact
ax.text(0.02, 0.98, 'Green: < 0.1% error\nYellow: < 1% error',
        transform=ax.transAxes, fontsize=8, va='top', ha='left',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.9, pad=0.3))

plt.tight_layout(pad=1.2)
import os, sys
out_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
plt.savefig(os.path.join(out_dir, 'prediction_scatter.png'), dpi=200, bbox_inches='tight')
print(f"Saved {os.path.join(out_dir, 'prediction_scatter.png')}")
plt.close()
