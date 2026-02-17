#!/usr/bin/env python3
"""
THE ROADMAP: From Tr(f(D^2)) to all 39 predictions.
Polished flowchart with clean alignment.
"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib

matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.size'] = 9

fig, ax = plt.subplots(figsize=(12, 14))
ax.set_xlim(0, 12)
ax.set_ylim(0, 14)
ax.axis('off')

# Colors
c_axiom = '#1a3c5e'
c_theorem = '#d5efe3'
c_theorem_edge = '#1e8449'
c_derived = '#fef3e2'
c_derived_edge = '#d4880f'
c_bar = '#1a3c5e'
c_arrow = '#444444'

def rbox(x, y, w, h, text, fc, ec, fontsize=9, bold=False, text_color='#1a1a1a'):
    rect = mpatches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.2",
                                     facecolor=fc, edgecolor=ec,
                                     linewidth=1.2)
    ax.add_patch(rect)
    weight = 'bold' if bold else 'normal'
    ax.text(x + w/2, y + h/2, text, ha='center', va='center',
            fontsize=fontsize, fontweight=weight, color=text_color,
            linespacing=1.4)

def arrow_down(x, y1, y2, label=''):
    ax.annotate('', xy=(x, y2), xytext=(x, y1),
                arrowprops=dict(arrowstyle='->', color=c_arrow, lw=1.5))
    if label:
        ax.text(x + 0.15, (y1+y2)/2, label, fontsize=7.5, color='#666',
                fontstyle='italic', va='center')

def arrow_fan(x_from, y_from, x_to, y_to):
    ax.annotate('', xy=(x_to, y_to), xytext=(x_from, y_from),
                arrowprops=dict(arrowstyle='->', color=c_arrow, lw=1.0, alpha=0.7))

# ── Title ──
ax.text(6, 13.5, 'From the Spectral Action to All of Physics',
        ha='center', fontsize=18, fontweight='bold', color=c_axiom)
ax.text(6, 13.05, r'$\mathrm{Tr}(f(D^2/\Lambda^2))$ on $M^4 \times S^5/\mathbb{Z}_3$'
        r'  $\rightarrow$  39 predictions, zero free parameters',
        ha='center', fontsize=11, color='#555')

# ── AXIOM box ──
bw = 7.5  # box width
cx = 6 - bw/2  # center x
rbox(cx, 11.8, bw, 0.9,
     r'AXIOM:  $S^5 / \mathbb{Z}_3$' + '\n'
     r'$d_1\!=\!6 \quad \lambda_1\!=\!5 \quad K\!=\!2/3 \quad \eta\!=\!2/9 \quad p\!=\!3$',
     c_axiom, c_axiom, fontsize=11, bold=True, text_color='white')

# ── Level 0 ──
arrow_down(6, 11.8, 11.15, 'Koide moment map')
lw = 6.0
rbox(6-lw/2, 10.4, lw, 0.7,
     r'Level 0:  $m_e = 1$  (unit)' + '\n' + 'Koide ground state  —  Theorem',
     c_theorem, c_theorem_edge, fontsize=9.5)

# ── Level 1 ──
arrow_down(6, 10.4, 9.75, 'Ghost Parseval energy')
rbox(6-lw/2, 8.95, lw, 0.75,
     r'Level 1:  $m_p / m_e = 6\pi^5$' + '\n' + 'Parseval fold energy of ghost modes  —  Theorem',
     c_theorem, c_theorem_edge, fontsize=9.5)

# ── Level 2 ──
arrow_down(6, 8.95, 8.3, r'APS lag  $\eta\lambda_1/p = 10/27$')
rbox(6-lw/2, 7.5, lw, 0.75,
     r'Level 2:  $1/\alpha = 137.038$' + '\n' + r'APS spectral asymmetry correction  —  Theorem',
     c_theorem, c_theorem_edge, fontsize=9.5)

# ── Level 3 ──
arrow_down(6, 7.5, 6.85, 'EM budget  −  ghost cost')
rbox(6-lw/2-0.75, 5.95, lw+1.5, 0.85,
     r'Level 3:  $v = 246.2$ GeV    $m_H = 125.3$ GeV' + '\n'
     r'$v/m_p = 2/\alpha - 35/3$  $\quad$  $m_H/m_p = 1/\alpha - 7/2$  —  Theorem',
     c_theorem, c_theorem_edge, fontsize=9.5)

# ── Level 4: four branches ──
branch_y_top = 5.95
branch_y_arrow = 5.3
branch_y_box = 4.0
bh = 1.2
bw_branch = 2.4

positions = [
    (1.2, 'Masses & Mixing',
     r'$m_\mu,\; m_\tau$  (Thm)' + '\n'
     r'6 quark masses (Thm)' + '\n'
     r'CKM: 9 elements' + '\n' + r'PMNS: 3 angles',
     'spectral\nratios'),
    (4.0, 'Gauge Sector',
     r'$\alpha_s = 0.1187\;(0.6\%)$' + '\n'
     r'$\bar\theta = 0$  (Thm)' + '\n'
     r'$N_g = 3$  (Thm)' + '\n'
     r'$\sin^2\!\theta_W = 3/8$',
     'ghost\nsplitting'),
    (6.8, 'Gravity & CC',
     r'$M_P$  to  $0.10\%$  (Thm)' + '\n'
     r'$\Lambda^{1/4}\!=\!2.22$ meV' + '\n'
     r'$n_s = 0.968,\; r = 0.003$' + '\n'
     r'Inflation: $N \approx 63$',
     'KK +\ntunneling'),
    (9.6, 'DM & Cosmology',
     r'$\Omega_{\rm DM}/\Omega_B = 16/3$' + '\n'
     r'$\eta_B = \alpha^4\!\eta$' + '\n'
     r'BH: singularity resolved' + '\n'
     r'QG: dissolved',
     'phase\ntransition'),
]

for bx, title, content, arrow_label in positions:
    cx_b = bx
    # Fan arrow
    arrow_fan(6, branch_y_top, cx_b + bw_branch/2, branch_y_box + bh + 0.05)
    # Arrow label
    ax.text(cx_b + bw_branch/2, branch_y_arrow, arrow_label,
            ha='center', va='center', fontsize=7, color='#888', fontstyle='italic')
    # Title
    ax.text(cx_b + bw_branch/2, branch_y_box + bh + 0.15, title,
            ha='center', va='bottom', fontsize=8.5, fontweight='bold', color='#555')
    # Box
    rbox(cx_b, branch_y_box, bw_branch, bh, content,
         c_derived, c_derived_edge, fontsize=7.5)

# ── Scorecard bar ──
rbox(0.3, 2.4, 11.4, 0.8,
     '19 Theorem   |   20 Derived   |   0 Gaps   |   '
     '39 predictions from 5 invariants + π   |   Zero free parameters',
     c_bar, c_bar, fontsize=10.5, bold=True, text_color='white')

# ── Bottom note ──
ax.text(6, 1.9, 'Every arrow is an explicit derivation in Supplement X (The Math-to-Physics Map).',
        ha='center', fontsize=8.5, color='#888', fontstyle='italic')
ax.text(6, 1.5, 'Every prediction has a verification script.  120 scripts.  All runnable.',
        ha='center', fontsize=8.5, color='#888', fontstyle='italic')

# ── Legend ── (lower left to avoid overlapping axiom/title)
legend_t = mpatches.Patch(facecolor=c_theorem, edgecolor=c_theorem_edge, linewidth=1.2, label='Theorem level')
legend_d = mpatches.Patch(facecolor=c_derived, edgecolor=c_derived_edge, linewidth=1.2, label='Derived level')
legend_a = mpatches.Patch(facecolor=c_axiom, edgecolor=c_axiom, linewidth=1.2, label='Axiom / Summary')
ax.legend(handles=[legend_a, legend_t, legend_d], loc='lower left',
          fontsize=9, framealpha=0.95, edgecolor='#ccc')

fig.subplots_adjust(left=0.03, right=0.97, top=0.96, bottom=0.04)
import os, sys
out_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
plt.savefig(os.path.join(out_dir, 'roadmap_figure.png'), dpi=200, bbox_inches='tight')
print(f"Saved {os.path.join(out_dir, 'roadmap_figure.png')}")
plt.close()
