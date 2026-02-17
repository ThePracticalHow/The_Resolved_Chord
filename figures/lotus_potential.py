#!/usr/bin/env python3
"""
THE LOTUS POTENTIAL: V(phi) with all key values labeled.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.size'] = 11
matplotlib.rcParams['axes.linewidth'] = 0.8

PI = np.pi
alpha = 1/137.036
d1 = 6; lam1 = 5; K = 2/3

phi_lotus = 1 - alpha*(d1 + lam1 + K)/2  # = 0.9574
phi_c = 0.60  # spectral phase transition
lambda_H = 0.1295
v_max = 2 * 0.938272 / alpha  # GeV (2*m_p/alpha)

phi = np.linspace(0, 1.05, 1000)
V = lambda_H / 4 * v_max**4 * (phi**2 - phi_lotus**2)**2
V_norm = V / V[0]  # normalize to V(0) = 1

fig, ax = plt.subplots(figsize=(10, 6))

ax.plot(phi, V_norm, 'k-', linewidth=2)

# Mark key points
ax.axvline(x=0, color='gray', linestyle=':', linewidth=0.5)
ax.axvline(x=phi_c, color='#e67e22', linestyle='--', linewidth=1.2, alpha=0.7)
ax.axvline(x=phi_lotus, color='#27ae60', linestyle='--', linewidth=1.2, alpha=0.7)
ax.axvline(x=1.0, color='#c0392b', linestyle='--', linewidth=1.2, alpha=0.7)

# Labels: positioned to avoid overlap
ax.annotate('$\\phi = 0$\nSmooth $S^5$\nNo physics',
            xy=(0, V_norm[0]), xytext=(0.08, 0.88),
            fontsize=9, ha='left', color='gray',
            arrowprops=dict(arrowstyle='->', color='gray', lw=0.8))

ax.annotate('$\\phi_c = 0.60$\nPhase transition\nInflation ends',
            xy=(phi_c, np.interp(phi_c, phi, V_norm)), xytext=(0.42, 0.62),
            fontsize=9, ha='center', color='#e67e22',
            arrowprops=dict(arrowstyle='->', color='#e67e22', lw=0.8))

ax.annotate('$\\phi_{\\mathrm{lotus}} = 0.9574$\nOur universe\n$V = 0$ (tree-level CC)',
            xy=(phi_lotus, 0), xytext=(0.78, 0.18),
            fontsize=9, ha='center', color='#27ae60',
            arrowprops=dict(arrowstyle='->', color='#27ae60', lw=0.8))

ax.annotate('$\\phi = 1$\nDried flower\n$v = 0$, all masses vanish',
            xy=(1.0, np.interp(1.0, phi, V_norm)), xytext=(0.92, 0.08),
            fontsize=9, ha='center', va='bottom', color='#c0392b',
            arrowprops=dict(arrowstyle='->', color='#c0392b', lw=0.8))

# Higgs mass curvature: moved down to avoid phi_c label
ax.annotate('$m_H^2 = V\'\'(\\phi_{\\mathrm{lotus}})$\n$= 125.3$ GeV',
            xy=(phi_lotus + 0.02, 0.002), xytext=(0.68, 0.35),
            fontsize=9, ha='center', color='#2c3e50',
            arrowprops=dict(arrowstyle='->', color='#2c3e50', lw=0.8))

# Barrier height: upper left, clear of phi=0
V_barrier = np.interp(0, phi, V_norm)
ax.annotate(f'Barrier $= V(0)^{{1/4}} \\approx 104$ GeV',
            xy=(0.02, V_barrier * 0.9), xytext=(0.18, 0.92),
            fontsize=9, ha='center', color='#7f8c8d',
            arrowprops=dict(arrowstyle='->', color='#7f8c8d', lw=0.6))

ax.set_xlabel('Fold depth $\\phi$', fontsize=13)
ax.set_ylabel('$V(\\phi) / V(0)$  (normalized)', fontsize=13)
ax.set_title('The LOTUS Potential: $V(\\phi) = \\frac{\\lambda_H}{4} v_{\\max}^4 (\\phi^2 - \\phi_{\\mathrm{lotus}}^2)^2$\n'
             'Lagrangian Of The Universe\'s Spectral State',
             fontsize=13, fontweight='bold')

ax.set_xlim(-0.05, 1.1)
ax.set_ylim(-0.05, 1.1)

# Add physics regions
ax.fill_betweenx([0, 1.1], 0, phi_c, alpha=0.03, color='blue', label='Substrate regime')
ax.fill_betweenx([0, 1.1], phi_c, 1.0, alpha=0.03, color='red', label='Information regime')

ax.legend(loc='upper right', fontsize=9, framealpha=0.9)  # upper right to avoid annotation overlap

plt.tight_layout(pad=1.0)
import os, sys
out_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
plt.savefig(os.path.join(out_dir, 'lotus_potential_figure.png'), dpi=200, bbox_inches='tight')
print(f"Saved {os.path.join(out_dir, 'lotus_potential_figure.png')}")
plt.close()
