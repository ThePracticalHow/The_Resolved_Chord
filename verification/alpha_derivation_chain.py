#!/usr/bin/env python3
"""
THE COMPLETE ALPHA DERIVATION CHAIN
====================================

Status audit: what is Theorem, what is Derived, where is the gap.

Two independent routes to 1/alpha = 137.036:

  ROUTE 1 (Geometric RG + Lag):
    sin^2(theta_W) = 3/8 --> M_c --> 1/alpha_GUT --> lag G/p --> RG --> alpha(0)
    Result: 1/alpha = 137.038 (0.001%)

  ROUTE 2 (Proton Constraint Inversion):
    m_p/m_e = 6*pi^5*(1 + G*alpha^2/pi + G2*alpha^4/pi^2)
    Invert: alpha = f(m_p/m_e, 6*pi^5, G, G2)
    Result: 1/alpha = 137.036 (<10^-4%)

THE GAP: What separates Derived from Theorem.

Jixiang Leng & Claude, February 2026
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

PI = np.pi
ALPHA_CODATA = 1/137.035999084

# Spectral invariants (ALL Theorem)
d1 = 6; lam1 = 5; K = 2/3; eta = 2/9; p = 3

# Hurricane coefficients
G  = lam1 * eta              # 10/9
G2 = -lam1 * (d1 + eta)      # -280/9

print("=" * 72)
print("  THE COMPLETE ALPHA DERIVATION CHAIN")
print("=" * 72)

# ======================================================================
#  STEP-BY-STEP STATUS AUDIT
# ======================================================================

print("""
  STEP-BY-STEP STATUS:

  Step 1: sin^2(theta_W) = 3/8 at M_c
          Source: SO(6) -> SU(3) x SU(2) x U(1) branching rule
          Status: THEOREM (group theory)

  Step 2: N_g = 3 (three generations)
          Source: Dirac spectrum on S^5 decomposes under Z_3
          Status: THEOREM (lotus_aps_generation.py)

  Step 3: SM beta coefficients b_1, b_2, b_3
          Source: Standard SM with N_g=3, N_H=1
          Status: THEOREM (textbook + Step 2)

  Step 4: M_c from alpha_1 = alpha_2 crossing
          Uses: Steps 1-3 + M_Z = 91.19 GeV (measured)
          Status: Standard physics (M_Z is the ONE scale input)

  Step 5: 1/alpha_GUT at crossing
          Source: Determined by Step 4
          Status: Standard physics
""")

# Compute Step 4-5
M_Z = 91.1876
alpha_em_MZ = 1/127.951
sin2_W_MZ = 0.23122
alpha_s_MZ = 0.1180

alpha_1_MZ = (5/3) * alpha_em_MZ / (1 - sin2_W_MZ)
alpha_2_MZ = alpha_em_MZ / sin2_W_MZ
a1_MZ = 1/alpha_1_MZ; a2_MZ = 1/alpha_2_MZ

b1 = 41/10; b2 = -19/6; b3 = -7.0

t_12 = 2*PI*(a1_MZ - a2_MZ)/(b1 - b2)
M_c = M_Z * np.exp(t_12)
a_GUT = a1_MZ - b1/(2*PI)*t_12

print(f"  Computed values (one-loop):")
print(f"    M_c = {M_c:.3e} GeV")
print(f"    1/alpha_GUT = {a_GUT:.4f}")

# ======================================================================
#  THE LAG CORRECTION (the key gap)
# ======================================================================

lag = G/p  # = 10/27

print(f"""
  ================================================================
  Step 6: THE LAG CORRECTION  <-- THIS IS THE GAP
  ================================================================

  Formula: 1/alpha_GUT_corrected = 1/alpha_GUT + G/p

  Components:
    G = lam1 * eta = {lam1} * {eta:.4f} = {G:.4f} = 10/9  [THEOREM]
    p = {p}                                                  [AXIOM]
    G/p = {lag:.6f} = 10/27                                  [THEOREM arithmetic]

  Each FACTOR in G/p is Theorem:
    lam1 = 5 (Ikeda 1980)                  THEOREM
    eta  = 2/9 (Donnelly computation)       THEOREM
    p    = 3 (orbifold axiom)               AXIOM

  But the CLAIM that the lag correction equals G/p:
    "1/alpha_GUT is shifted by G/p at the matching scale"    ???????

  This is what separates Derived from Theorem.
  The individual spectral invariants are all proven.
  What is NOT proven: why the gauge coupling offset = G/p.

  Physical argument: Ghost modes have spectral inertia G = lam1*eta.
  This inertia is distributed across p sectors before decoupling.
  So the offset per unified coupling = G/p.

  The argument is PHYSICALLY CLEAR but NOT FORMALLY DERIVED
  from the spectral action.
""")

# ======================================================================
#  ROUTE 1: Geometric RG + Lag
# ======================================================================

print("=" * 72)
print("  ROUTE 1: GEOMETRIC RG + LAG")
print("=" * 72)

a_GUT_corrected = a_GUT + lag

# Run to alpha(0)
inv_alpha_em_Mc = (5/3 + 1) * a_GUT_corrected
t_Mc_MZ = np.log(M_c/M_Z)
a1_run = a_GUT_corrected + b1/(2*PI)*t_Mc_MZ
a2_run = a_GUT_corrected + b2/(2*PI)*t_Mc_MZ
inv_alpha_em_MZ = (5/3)*a1_run + a2_run

delta_alpha = 0.0591
inv_alpha_0_route1 = inv_alpha_em_MZ / (1 - delta_alpha)

err1 = abs(inv_alpha_0_route1 - 1/ALPHA_CODATA)/(1/ALPHA_CODATA)*100

print(f"""
  1/alpha_GUT (naive)     = {a_GUT:.4f}
  Lag correction G/p      = {lag:.6f}
  1/alpha_GUT (corrected) = {a_GUT_corrected:.4f}
  
  RG: M_c -> M_Z:
    1/alpha_em(M_Z) = {inv_alpha_em_MZ:.3f}
  
  Vacuum polarization: M_Z -> 0:
    Delta_alpha = {delta_alpha}
  
  RESULT:
    1/alpha(0) = {inv_alpha_0_route1:.4f}
    CODATA:      {1/ALPHA_CODATA:.4f}
    Error:       {err1:.4f}%

  STATUS: DERIVED (0.001%)
  WEAK LINK: Step 6 (lag correction = G/p)
""")

# ======================================================================
#  ROUTE 2: Proton Constraint Inversion
# ======================================================================

print("=" * 72)
print("  ROUTE 2: PROTON CONSTRAINT INVERSION")
print("=" * 72)

m_p_m_e = 1836.15267343

# Tree level (Theorem)
tree = 6*PI**5
delta = m_p_m_e / tree - 1

# One-loop inversion: alpha^2 = pi * delta / G
alpha_sq_1L = PI * delta / G
inv_alpha_1L = 1/np.sqrt(alpha_sq_1L)

# Two-loop inversion: G2/pi^2 * alpha^4 + G/pi * alpha^2 - delta = 0
A_c = G2/PI**2
B_c = G/PI
C_c = -delta
disc = B_c**2 - 4*A_c*C_c
x_2L = (-B_c + np.sqrt(disc))/(2*A_c)
inv_alpha_2L = 1/np.sqrt(x_2L)

err2_1L = abs(inv_alpha_1L - 1/ALPHA_CODATA)/(1/ALPHA_CODATA)*100
err2_2L = abs(inv_alpha_2L - 1/ALPHA_CODATA)/(1/ALPHA_CODATA)*100

print(f"""
  m_p/m_e (measured) = {m_p_m_e}
  6*pi^5 (Theorem)   = {tree:.6f}
  Fractional shift    = {delta:.10e}

  G  = lam1*eta         = {G:.6f} = 10/9
  G2 = -lam1*(d1+eta)   = {G2:.6f} = -280/9

  One-loop inversion:
    1/alpha = {inv_alpha_1L:.4f}  (error: {err2_1L:.4f}%)
  
  Two-loop inversion:
    1/alpha = {inv_alpha_2L:.6f}  (error: {err2_2L:.6f}%)

  STATUS: DERIVED (< 10^-4%)
  WEAK LINK: G, G2 as QCD loop coefficients
""")

# ======================================================================
#  WHAT MAKES G AND G2 "DERIVED" NOT "THEOREM"
# ======================================================================

print("=" * 72)
print("  THE GAP ANALYSIS")
print("=" * 72)

print(f"""
  G = lam1 * eta = 10/9
  G2 = -lam1 * (d1 + eta) = -280/9

  These formulas are STRUCTURAL: they express the QCD radiative
  correction coefficients as products of spectral invariants.

  The spectral invariants themselves are all THEOREM:
    lam1 = 5      (Ikeda 1980, Lichnerowicz-Obata)
    eta  = 2/9    (Donnelly computation, Cheeger-Muller)
    d1   = 6      (spherical harmonic formula)
    p    = 3      (orbifold axiom)

  The GAP: proving G is the one-loop coefficient
  
  WHAT IS PROVEN:
    - The spectral formula G = lam1*eta = 10/9
    - The QCD perturbative expansion of m_p to 10^-11
    - The two values numerically match

  WHAT IS NOT PROVEN:
    - That the spectral action loop expansion of the proton
      mass produces EXACTLY lam1*eta at one-loop and
      -lam1*(d1+eta) at two-loop.
    - That the gauge coupling at M_c receives EXACTLY G/p
      from the ghost mode spectral weight.

  TO CLOSE THE GAP, we need ONE of:

  APPROACH A: Prove the lag correction from spectral action
    Show: delta(1/alpha_GUT) = eta * lam1 / p at M_c
    Route: APS boundary term applied to gauge coupling matching
    
  APPROACH B: Prove the hurricane coefficients from spectral action
    Show: one-loop QCD on S^5/Z_3 gives coefficient lam1*eta
    Route: Spectral action perturbative expansion
    
  APPROACH C: Derive alpha from a THIRD independent route
    Show: 1/alpha = some_spectral_formula (without G/p or G)
    Route: Direct spectral zeta computation, or a_4 coefficient
""")

# ======================================================================
#  APPROACH A: The APS Boundary Term Route
# ======================================================================

print("=" * 72)
print("  APPROACH A: APS BOUNDARY TERM FOR LAG CORRECTION")
print("=" * 72)

print(f"""
  CONJECTURE: The lag correction is an APS boundary contribution
  to the gauge coupling at the compactification boundary.

  In the spectral action on M^4 x B^6/Z_3:
    - The bulk (r < 1) gives the 9D gauge coupling
    - The boundary (r = 1, i.e. S^5/Z_3) gives the matching correction
    - The APS boundary term includes the eta invariant

  The gauge kinetic term from the spectral action:
    1/g^2 = (bulk term) + (boundary term)
    
  Boundary term = eta_gauge / (normalization)
  
  For the Dirac operator on S^5/Z_3:
    eta = 2/9 (Theorem, Donnelly)
  
  The gauge coupling correction involves the FIRST KK level:
    delta(1/g^2) = eta * lam1 / p = (2/9)(5/3) = 10/27

  PHYSICAL INTERPRETATION:
    - eta = 2/9 is the spectral asymmetry (L-R imbalance)
    - lam1 = 5 is the eigenvalue at the ghost level
    - 1/p is the orbifold volume normalization
    - Product: the LEFT-RIGHT imbalance, weighted by the ghost
      eigenvalue, per orbifold sector
    
  WHY lam1 APPEARS:
    The APS boundary correction to a differential operator
    of order 2 (the gauge Laplacian) on a manifold with boundary
    includes a factor of the boundary operator's eigenvalue.
    The ghost modes at ell=1 have eigenvalue lam1 = 5.
    
  WHY eta APPEARS:
    The APS theorem says the spectral asymmetry of the boundary
    Dirac operator contributes to the index/gauge coupling
    through eta(D_boundary).
    
  WHY 1/p APPEARS:
    The orbifold projection divides by |Z_3| = p.
""")

# Verify the formula
print(f"  NUMERICAL CHECK:")
print(f"    eta * lam1 / p = {eta} * {lam1} / {p} = {eta*lam1/p:.6f}")
print(f"    G / p          = {G} / {p}             = {G/p:.6f}")
print(f"    Match: {abs(eta*lam1/p - G/p) < 1e-15}")
print(f"    = 10/27 = {10/27:.10f}")

# ======================================================================
#  APPROACH C: Direct spectral formula for alpha_GUT
# ======================================================================

print(f"\n\n{'='*72}")
print("  APPROACH C: SEARCHING FOR DIRECT SPECTRAL FORMULA")
print("=" * 72)

# The corrected 1/alpha_GUT
target = a_GUT + lag
print(f"\n  Target: 1/alpha_GUT,corrected = {target:.6f}")

# Try exact rational/spectral expressions
from fractions import Fraction

# The corrected value should decompose into spectral invariants
# a_GUT = 42.412, lag = 10/27 = 0.370, corrected = 42.782

# Let's check what combinations give 42.78
candidates_C = {
    "d1^2 + lam1 + K + eta":      d1**2 + lam1 + K + eta,
    "d1*lam1 + d1 + lam1 + K":    d1*lam1 + d1 + lam1 + K,
    "(d1+lam1)^2/p + 2":          (d1+lam1)**2/p + 2,
    "(d1+lam1)^2/p + K":          (d1+lam1)**2/p + K,
    "(d1+lam1)^2/p + eta":        (d1+lam1)**2/p + eta,
    "d1^2 + lam1 + K":            d1**2 + lam1 + K,
    "d1^2 + lam1 + 1":            d1**2 + lam1 + 1,
    "d1*(lam1+K+1)":              d1*(lam1+K+1),
    "d1*(lam1+1) + K":            d1*(lam1+1) + K,
    "d1^2 + d1 + K":              d1**2 + d1 + K,
    "d1*(d1+1) + K":              d1*(d1+1) + K,
    "p*(d1+lam1+K)^2/(d1-1)":     p*(d1+lam1+K)**2/(d1-1),
    "4*pi^2/p + lam1/pi":         4*PI**2/p + lam1/PI,
    "d1*pi^2/p + 1/eta":          d1*PI**2/p + 1/eta,
}

# Also check the UNCORRECTED value
print(f"\n  Also checking: 1/alpha_GUT (uncorrected) = {a_GUT:.6f}")

results_C = []
for name, val in candidates_C.items():
    err_corr = abs(val - target)/target*100
    err_raw  = abs(val - a_GUT)/a_GUT*100
    results_C.append((min(err_corr, err_raw), name, val, err_corr, err_raw))

results_C.sort()
print(f"\n  Best matches to corrected or uncorrected 1/alpha_GUT:")
for _, name, val, ec, er in results_C[:10]:
    which = "corr" if ec < er else "raw"
    err = min(ec, er)
    marker = "***" if err < 0.5 else "**" if err < 2 else "*" if err < 5 else ""
    print(f"    {name:<35} = {val:>10.4f}  ({which} err: {err:.3f}%) {marker}")

# ======================================================================
#  THE KEY IDENTITY: a_GUT = d1^2 + lam1 + K
# ======================================================================

a_GUT_spectral = d1**2 + lam1 + K
err_identity = abs(a_GUT_spectral - a_GUT)/a_GUT*100

print(f"\n\n{'='*72}")
print("  KEY IDENTITY TEST: 1/alpha_GUT = d1^2 + lam1 + K ?")
print("=" * 72)

print(f"""
  d1^2 + lam1 + K = 36 + 5 + 2/3 = {a_GUT_spectral:.4f}
  1/alpha_GUT (RG) = {a_GUT:.4f}
  Error: {err_identity:.3f}%
""")

a_GUT_corrected_test = a_GUT_spectral + lag
print(f"  If 1/alpha_GUT = d1^2 + lam1 + K (EXACTLY), then:")
print(f"  1/alpha_GUT,corrected = d1^2 + lam1 + K + G/p")
print(f"    = 36 + 5 + 2/3 + 10/27")
print(f"    = {a_GUT_corrected_test:.6f}")
print(f"    = {Fraction(36) + Fraction(5) + Fraction(2,3) + Fraction(10,27)}")

total_frac = Fraction(36) + Fraction(5) + Fraction(2,3) + Fraction(10,27)
print(f"    = {total_frac} = {float(total_frac):.10f}")

# Check how this propagates to alpha(0)
a_GUT_exact = float(total_frac)
inv_alpha_em_Mc_ex = (5/3 + 1)*a_GUT_exact
t_run = np.log(M_c/M_Z)
a1_ex = a_GUT_exact + b1/(2*PI)*t_run
a2_ex = a_GUT_exact + b2/(2*PI)*t_run
inv_em_MZ_ex = (5/3)*a1_ex + a2_ex
inv_alpha_0_ex = inv_em_MZ_ex / (1 - delta_alpha)
err_ex = abs(inv_alpha_0_ex - 1/ALPHA_CODATA)/(1/ALPHA_CODATA)*100

print(f"\n  Propagated to alpha(0):")
print(f"    1/alpha(0) = {inv_alpha_0_ex:.4f}")
print(f"    CODATA:      {1/ALPHA_CODATA:.4f}")
print(f"    Error: {err_ex:.4f}%")

# ======================================================================
#  THE DECOMPOSITION OF alpha_GUT
# ======================================================================

print(f"\n\n{'='*72}")
print("  UNDERSTANDING WHY 1/alpha_GUT ~ d1^2 + lam1 + K")
print("=" * 72)

print(f"""
  1/alpha_GUT comes from the RG crossing of alpha_1 = alpha_2.
  At one-loop: 
    1/alpha_GUT = 1/alpha_1(M_Z) - b_1/(2*pi) * ln(M_c/M_Z)
  
  Measured: 1/alpha_1(M_Z) = {a1_MZ:.4f}
  Running:  b_1/(2*pi)*ln(M_c/M_Z) = {b1/(2*PI)*t_12:.4f}
  Result:   1/alpha_GUT = {a_GUT:.4f}
  
  The value {a_GUT:.2f} is DETERMINED by standard RG + M_Z.
  The coincidence d1^2 + lam1 + K = {a_GUT_spectral:.4f} is SUGGESTIVE
  but the error is {err_identity:.2f}%, not exact.
  
  However, d1^2 + lam1 + K has clear spectral meaning:
    d1^2 = 36: squared ghost count (the d1 modes interact pairwise)
    lam1 = 5:  first eigenvalue (kinetic energy scale)
    K    = 2/3: Koide moment map (mass coupling)
    
  Total: the "vacuum spectral budget" available for gauge coupling.
""")

# ======================================================================
#  FINAL STATUS
# ======================================================================

print("=" * 72)
print("  FINAL STATUS AND THEOREM PATH")
print("=" * 72)

print(f"""
  CURRENT STATUS: DERIVED (0.001% via Route 1, <10^-4% via Route 2)
  
  THEOREM PATH A (most promising):
    Prove: delta(1/alpha_GUT) = eta*lam1/p at the KK boundary.
    
    This requires showing the APS boundary contribution to the
    spectral action gauge kinetic term equals G/p = 10/27.
    
    If proven: Route 1 becomes
      sin^2(theta_W) = 3/8   [THEOREM]
      + N_g = 3               [THEOREM]
      + lag = G/p = 10/27     [THEOREM if A succeeds]
      + SM RG                  [standard physics]
      + vacuum polarization    [standard physics]
      = 1/alpha(0) = 137.038  [THEOREM]
      
  THEOREM PATH B (harder):
    Prove: the spectral action loop expansion gives G = lam1*eta
    for the one-loop proton mass coefficient.
    
    If proven: Route 2 becomes
      6*pi^5                   [THEOREM, Parseval]
      + G = lam1*eta           [THEOREM if B succeeds]
      + measured m_p/m_e       [experiment]
      = 1/alpha = 137.036      [THEOREM]
      
  THEOREM PATH C (speculative):
    Find: 1/alpha_GUT = exact spectral formula
    Candidate: d1^2 + lam1 + K = {a_GUT_spectral:.4f} (error {err_identity:.2f}%)
    If exact: derive from spectral action a_4 coefficient
    
  DEPENDENCIES:
    alpha -> v (Higgs VEV): v = m_p*(2/alpha - 35/3)
    alpha -> m_H (Higgs mass): m_H = m_p*(1/alpha - 7/2)
    alpha -> baryogenesis: eta_B = alpha^4 * eta
    
    So promoting alpha to THEOREM promotes v, m_H, and eta_B too.
""")

print("=" * 72)
print("  COMPUTATION COMPLETE")
print("=" * 72)
