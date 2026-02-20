#!/usr/bin/env python3
"""
FIGURE E — THE SEVEN EQUATIONS FROM ONE TRACE
Shows how 7 classical physics equations all derive from Tr(f(D²/Λ²))
on M⁴ × S⁵/Z₃.  Visual companion to Supplement XIV.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib
import os, sys

matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.size'] = 10

# ── Palette ──
BG       = '#fafafa'
AXIOM_C  = '#1a3c5e'
GOLD     = '#b7950b'
TEXT_C   = '#1a1a2e'
GRAY     = '#7f8c8d'

# Equation colors — one per branch
EQ_COLORS = [
    '#c0392b',  # E=mc²  (red)
    '#2980b9',  # F=ma   (blue)
    '#27ae60',  # Maxwell (green)
    '#8e44ad',  # Dirac  (purple)
    '#d4880f',  # Schrödinger (orange)
    '#1a5276',  # Friedmann (dark blue)
    '#922b21',  # 2nd Law (dark red)
]

EQ_NAMES = [
    r'$E = mc^2$',
    r'$F = ma$',
    r'Maxwell',
    r'Dirac',
    r'Schrödinger',
    r'Friedmann',
    r'$dS \geq 0$',
]

EQ_LAWS = [
    r'$E^2 = p^2c^2 + m^2c^4$',
    r'$F = -m\nabla\Phi$',
    r'$\partial_\mu F^{\mu\nu} = J^\nu$',
    r'$(i\gamma^\mu\partial_\mu - m)\psi = 0$',
    r'$i\hbar\partial_t\psi = H\psi$',
    r'$H^2 = \frac{8\pi G}{3}\rho + \frac{\Lambda}{3}$',
    r'$dS \geq 0$  (time arrow)',
]

EQ_MECHANISMS = [
    r'$\eta_D = i/9 \Rightarrow$ Lorentz sig.',
    r'$a_2 \Rightarrow$ EH $\Rightarrow$ weak field',
    r'$D \to D+A \Rightarrow$ YM $\Rightarrow$ U(1)',
    r'$D_K\phi_n = m_n\phi_n$ on $S^5/\mathbb{Z}_3$',
    r'Dirac $\rightarrow$ NR limit $\rightarrow$ Schrödinger',
    r'$a_0, a_2 \Rightarrow \Lambda, G \Rightarrow$ FLRW',
    r'$\eta = 2/9 \neq 0 \Rightarrow$ T-breaking',
]

fig, ax = plt.subplots(figsize=(12, 10), facecolor=BG)
ax.set_xlim(0, 12)
ax.set_ylim(0, 10)
ax.axis('off')

def rbox(x, y, w, h, text, fc, ec, fontsize=9, bold=False, text_color='#1a1a1a'):
    rect = mpatches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.15",
                                     facecolor=fc, edgecolor=ec,
                                     linewidth=1.5)
    ax.add_patch(rect)
    weight = 'bold' if bold else 'normal'
    ax.text(x + w/2, y + h/2, text, ha='center', va='center',
            fontsize=fontsize, fontweight=weight, color=text_color,
            linespacing=1.3)

# ── Title ──
ax.text(6, 9.5, '7 Equations from One Trace',
        ha='center', fontsize=20, fontweight='bold', color=AXIOM_C)
ax.text(6, 9.05, r'Classical physics derived from $\mathrm{Tr}(f(D^2/\Lambda^2))$ '
        r'on $M^4 \times S^5/\mathbb{Z}_3$',
        ha='center', fontsize=12, color='#555')

# ── Central axiom ──
rbox(3.0, 7.8, 6.0, 0.9,
     r'$\mathrm{Tr}(f(D^2/\Lambda^2))$' + '\n'
     r'on $M^4 \times S^5/\mathbb{Z}_3$   ($d_1\!=\!6,\; \lambda_1\!=\!5,\; K\!=\!2/3,\; \eta\!=\!2/9,\; p\!=\!3$)',
     AXIOM_C, AXIOM_C, fontsize=11, bold=True, text_color='white')

# ── Seven equation boxes, arranged in two rows ──
# Row 1: 4 equations (top row)
row1_y = 5.0
row1_h = 2.0
row1_w = 2.6
row1_gap = 0.2
row1_start = 0.5

# Row 2: 3 equations (bottom row)
row2_y = 2.0
row2_h = 2.0
row2_w = 3.2
row2_gap = 0.6
row2_start = 1.0

# Place top row (4 equations)
for i in range(4):
    x = row1_start + i * (row1_w + row1_gap)
    c = EQ_COLORS[i]
    
    # Arrow from axiom to box
    ax.annotate('', xy=(x + row1_w/2, row1_y + row1_h),
                xytext=(6, 7.8),
                arrowprops=dict(arrowstyle='->', color=c, lw=1.2, alpha=0.6))
    
    # Mechanism label along arrow
    mid_y = (7.8 + row1_y + row1_h) / 2
    
    # Box
    rbox(x, row1_y, row1_w, row1_h,
         EQ_NAMES[i] + '\n\n' + EQ_LAWS[i] + '\n\n' + EQ_MECHANISMS[i],
         '#ffffff', c, fontsize=8)
    
    # Equation name header
    ax.text(x + row1_w/2, row1_y + row1_h - 0.15, EQ_NAMES[i],
            ha='center', va='top', fontsize=13, fontweight='bold', color=c)

# Place bottom row (3 equations)
for i in range(3):
    j = i + 4
    x = row2_start + i * (row2_w + row2_gap)
    c = EQ_COLORS[j]
    
    # Arrow from axiom to box
    ax.annotate('', xy=(x + row2_w/2, row2_y + row2_h),
                xytext=(6, 7.8),
                arrowprops=dict(arrowstyle='->', color=c, lw=1.2, alpha=0.6))
    
    # Box
    rbox(x, row2_y, row2_w, row2_h,
         EQ_NAMES[j] + '\n\n' + EQ_LAWS[j] + '\n\n' + EQ_MECHANISMS[j],
         '#ffffff', c, fontsize=8.5)
    
    # Equation name header
    ax.text(x + row2_w/2, row2_y + row2_h - 0.15, EQ_NAMES[j],
            ha='center', va='top', fontsize=13, fontweight='bold', color=c)

# ── Footer ──
rbox(1.5, 0.6, 9.0, 0.7,
     '7 equations  •  1 trace  •  0 assumptions beyond ' + r'$\mathrm{Tr}(f(D^2/\Lambda^2))$'
     + '  •  See Supplement XIV',
     AXIOM_C, AXIOM_C, fontsize=10, bold=True, text_color='white')

ax.text(6, 0.2, 'Born rule excepted — see Supplement XIV, Ch. 5',
        ha='center', fontsize=8, color=GRAY, fontstyle='italic')

fig.subplots_adjust(left=0.02, right=0.98, top=0.96, bottom=0.02)
out_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
outpath = os.path.join(out_dir, 'equations_derivation.png')
plt.savefig(outpath, dpi=200, bbox_inches='tight', facecolor=BG)
print(f"Saved {outpath}")
plt.close()
