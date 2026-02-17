#!/usr/bin/env python3
"""
QUARK KK CLOSURE: YUKAWA THRESHOLD CORRECTIONS
==============================================

Hypothesis: 
The discrepancy between the Geometric Koide (K=2/3) and the Observed SM Koide (K~0.9 at Mc)
is due to KK threshold corrections that ENHANCE the Top quark mass relative to others.

Mechanism:
The Top quark lives in the 'Heavy' sectors (Chi1/Chi2) which couple resonantly to the Higgs (Chi1).
The loop corrections from the KK tower of these sectors are stronger/positive, enriching the Top.
Light quarks live in the 'Light' sector (Chi0 - cone tip), with suppressed/negative corrections.

Assignments (from chiral analysis):
- Higgs: Chi1 (charge w)
- Top:   L(Chi2) -> R(Chi1). Coupling: w^-2 * w * w = 1. (Allowed, Resonant)
- Bottom: L(Chi2) -> R(Chi0). Coupling: w^-2 * 1 * w^2 = 1. (Allowed, Suppressed overlap)
- Charm/Up: involve Chi0, suppressed or forbidden tree-level.

Calculation:
Compute delta_y / y for the Top channel vs Bottom/Light channels.
delta_y ~ Sum_n (Overlap_n * LoopFactor_n)

Jixiang Leng & Claude, February 2026
"""

import numpy as np

print("="*60)
print("  QUARK KK CLOSURE: YUKAWA RENORMALIZATION")
print("="*60)

# Constants
R = 1.0  # Compactification radius (units of 1/Mc)
LAMBDA_CUTOFF = 10.0 # Cutoff in units of Mc (approximating finite string scale or similar)
zeta_3 = 1.202
PI = np.pi

def loop_factor(m_kk, m_higgs_bound):
    """
    Generic 1-loop vertex correction factor for a mode of mass m_kk.
    For Yukawa, structure is typically positive for scalar exchange?
    In 4D: ~ + y^2 / 16pi^2 * ln(...)
    In 5D: Sum over tower.
    
    If Higgs is a bound state of the geometry, coupling is geometric.
    
    Ansatz: The correction scales as 1/M_n^2 (decoupling) but summed over density.
    """
    if m_kk == 0: return 0
    return 1.0 / (m_kk * R)**2

def overlap_factor(l, sector_type):
    """
    Overlap integral of the mode l with the zero mode.
    Twisted sectors (Chi1/Chi2) have wavefunctions concentrated on the fold.
    Untwisted sector (Chi0) is concentrated at the cone tip.
    
    Model:
    - Twist (Chi1/Chi2): Overlap ~ 1 (Normal)
    - Tip (Chi0): Overlap ~ 1/sqrt(l) (Suppressed high modes?)
    
    Actually, let's use the degeneracy/multiplicity as a proxy for 'number of paths'.
    """
    # Degeneracy grows as l^4 for S^5 via dim_H
    deg = (l+1)*(l+2)**2*(l+3)/12.0
    
    if sector_type == 'Twist':
        # Full coupling
        return deg
    elif sector_type == 'Tip':
        # Suppressed coupling (reduced phase space at tip?)
        # Say suppressed by geometric factor 1/d1? 
        return deg * 0.1 
    return 0

def compute_correction(sector_L, sector_R):
    """
    Compute total correction delta_y/y summing over KK towers.
    L and R are sector names ('Chi0', 'Chi1', 'Chi2').
    """
    total_delta = 0.0
    
    # Tower summation loop
    for l in range(1, 100):
        # We need intermediate states. 
        # Loop: Q_L -> [Scalar_KK + Q_R_KK] -> Q_R
        # Simplified: Summing effective contributions from modes up to cutoff.
        
        m_kk = l # Mass approx linear in l
        
        # Determine overlap based on sectors involved
        # If Loop involves Chi1/Chi2 (Twist), use strong overlap.
        # If Loop involves Chi0 (Tip), use weak overlap.
        
        is_twist_L = (sector_L in ['Chi1', 'Chi2'])
        is_twist_R = (sector_R in ['Chi1', 'Chi2'])
        
        # Effective coupling strength of this mode
        strength = 1.0
        if is_twist_L: strength *= 1.0 
        else: strength *= 0.1
        
        if is_twist_R: strength *= 1.0
        else: strength *= 0.1
        
        # Density of states
        deg = (l+1)*(l+2)**2*(l+3)/12.0
        
        term = strength * deg * (1.0 / m_kk**2)
        total_delta += term
        
    # Normalized roughly by 16pi^2 factor for loop
    return total_delta * (1.0 / (16 * PI**2))

# ------------------------------------------------------------------
# 1. Compute corrections for Top vs Bottom
# ------------------------------------------------------------------

# Top: L(Chi2) -> R(Chi1)
dy_top = compute_correction('Chi2', 'Chi1')

# Bottom: L(Chi2) -> R(Chi0)
dy_bot = compute_correction('Chi2', 'Chi0')

# Charm (Hypothetical): L(Chi0) -> R(Chi1)?
dy_charm = compute_correction('Chi0', 'Chi1')

# Up (Hypothetical): L(Chi0) -> R(Chi0)?
dy_up = compute_correction('Chi0', 'Chi0')

print(f"Computed Yukawa Proportional Corrections (Delta y / y):")
print(f"  Top (Twist-Twist): {dy_top:.4e}")
print(f"  Bot (Twist-Tip):   {dy_bot:.4e}")
print(f"  Chm (Tip-Twist):   {dy_charm:.4e}")
print(f"  Up  (Tip-Tip):     {dy_up:.4e}")

# ------------------------------------------------------------------
# 2. Impact on Koide Ratio
# ------------------------------------------------------------------
# Assume starting with a 'Democratic' K=2/3 distribution.
# m1 = 1, m2 = 1, m3 = 1  => K = 1/3 (Too democratic)
# m1 = 0.001, m2 = 0.01, m3 = 1 => K ~ 1 (Too hierarchical)
# K=2/3 is intermediate.

# Let's take the UV masses implied by Yukawa Universality (Geometric):
# Scaled from Lepton masses.
m_u_geo = 0.511e-3
m_c_geo = 105.66e-3
m_t_geo = 1776.86e-3 
# (This gives K=2/3 exact, but the hierarchy is Lepton-like, not Top-like)

# Wait. m_t_geo is 1.77 GeV? 
# NO! m_t saturates the fold -> v/sqrt(2) = 174 GeV.
# But K depends on ratios.
# If m_t = 174, does m_c/m_t = m_mu/m_tau?
# Yes, that's Yukawa Universality.
# m_c = 174 * (0.105/1.77) = 10.3 GeV
# m_u = 174 * (0.0005/1.77) = 0.05 GeV
# Let's check K for this set.
# sqrt(174)=13.19, sqrt(10.3)=3.21, sqrt(0.05)=0.22. Sum=16.62. Sq=276.
# Sum m = 174 + 10.3 + 0.05 = 184.35
# K = 184.35 / 276 = 0.667. (2/3). 
# So the Geometric Starting Point is (0.05, 10.3, 174).

print("\nGeometric Starting Point (K=2/3):")
print(f"  m_u = 0.05, m_c = 10.3, m_t = 174.0")

# Apply Threshold Corrections
# y_eff = y_geo * (1 + delta)
# m_eff = m_geo * (1 + delta)

m_t_new = 174.0 * (1 + dy_top)
m_c_new = 10.3 * (1 + dy_charm)
m_u_new = 0.05 * (1 + dy_up)

print(f"\nCorrected Masses (with KK thresholds):")
print(f"  m_u = {m_u_new:.4f}, m_c = {m_c_new:.4f}, m_t = {m_t_new:.4f}")

# Compute new K
def calc_k(m1, m2, m3):
    s = m1+m2+m3
    sq = (np.sqrt(m1)+np.sqrt(m2)+np.sqrt(m3))**2
    return s/sq

K_new = calc_k(m_u_new, m_c_new, m_t_new)
print(f"\nResulting Koide Ratio: K = {K_new:.4f}")

gap_closed = (K_new - 0.6667)
target_gap = (0.90 - 0.6667)

print(f"Shift: {gap_closed:+.4f}")
print(f"Target Shift (to match SM running): +{target_gap:.4f}")

if gap_closed > 0:
    print("SUCCESS: KK corrections increase K (Hierarchy enhanced)!")
    print(f"Fraction of gap closed: {gap_closed/target_gap * 100:.1f}%")
else:
    print("FAILURE: KK corrections decrease K (Hierarchy flattened).")

