#!/usr/bin/env python3
"""
NEUTRON-PROTON MASS SPLITTING FROM SPECTRAL QUARK MASSES
==========================================================

The n-p splitting has TWO contributions:
  1. Quark mass difference: m_d - m_u (makes neutron heavier)
  2. EM self-energy: proton is charged (makes proton heavier)

Both inputs are spectral. Can we get delta_m to better than 5.9%?

The Cottingham formula (standard approach):
  delta_m = (m_d - m_u) * <N|qqbar|N>/m_N - alpha*EM_correction

All inputs from spectral piercing depths + spectral alpha.

Jixiang Leng & Claude, February 2026
"""

import sys, io
if sys.stdout.encoding and sys.stdout.encoding.lower() not in ('utf-8', 'utf8'):
    try:
        sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
        sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')
    except Exception:
        pass

import numpy as np
from fractions import Fraction

d1 = 6; lam1 = 5; K = Fraction(2,3); eta = Fraction(2,9); p = 3
PI = np.pi; alpha = 1/137.036
D_wall = 1 + d1; D_bulk = d1 + lam1

m_e_MeV = 0.51100
m_p_MeV = m_e_MeV * d1 * PI**5 * (1 + float(lam1*eta)*alpha**2/PI)
v_GeV = m_p_MeV/1000 * (2/alpha - (d1 + lam1 + float(K)))
m_tau_MeV = 1776.86; m_mu_MeV = 105.658

# Spectral quark masses
sigma_u = -PI
sigma_d = 2*PI/3 + 10/81
m_u_MeV = (v_GeV/np.sqrt(2)) * (m_e_MeV/m_tau_MeV) * np.exp(sigma_u) * 1000
m_d_MeV = m_e_MeV * np.exp(sigma_d)

print("=" * 72)
print("  NEUTRON-PROTON MASS SPLITTING FROM SPECTRAL QUARKS")
print("=" * 72)

# =====================================================================
#  THE TWO CONTRIBUTIONS
# =====================================================================
print(f"\n{'='*72}")
print("THE TWO CONTRIBUTIONS TO delta_m")
print(f"{'='*72}")

print(f"""
  SPECTRAL QUARK MASSES:
    m_u = {m_u_MeV:.4f} MeV  (from sigma_u = -pi, chi_1 angular)
    m_d = {m_d_MeV:.4f} MeV  (from sigma_d = 2pi/3 + 10/81, chi_2 spectral)
    m_d - m_u = {m_d_MeV - m_u_MeV:.4f} MeV

  PDG quark masses (MS-bar at 2 GeV):
    m_u = 2.16 MeV,  m_d = 4.67 MeV,  m_d - m_u = 2.51 MeV

  Spectral m_d - m_u = {m_d_MeV - m_u_MeV:.3f} MeV vs PDG 2.51 MeV
  ({abs(m_d_MeV - m_u_MeV - 2.51)/2.51*100:.1f}% error)
""")

dm_quark = m_d_MeV - m_u_MeV  # spectral quark mass difference

# =====================================================================
#  THE DASHEN-COTTINGHAM DECOMPOSITION
# =====================================================================
print(f"{'='*72}")
print("THE DASHEN-COTTINGHAM DECOMPOSITION")
print(f"{'='*72}")

# The n-p mass difference decomposes as:
# delta_m(n-p) = delta_m(QCD) + delta_m(EM)
#
# delta_m(QCD) = (m_d - m_u) * sigma_term / (m_u + m_d)
#   where sigma_term ~ 45-60 MeV (pion-nucleon sigma term)
#
# delta_m(EM) ~ -0.76 MeV (proton EM self-energy > neutron EM)
#
# In the spectral framework:
# sigma_term = sigma_piN = spectral quantity related to f_pi and g_A

# The pion-nucleon sigma term from spectral data:
# sigma_piN = (m_u + m_d) * <N|uubar+ddbar|N> / (2*m_N)
# From ChPT: sigma_piN ~ f_pi^2 * m_pi^2 / m_N ~ 50 MeV
f_pi_MeV = float(K**2 * eta) * m_p_MeV  # 92.7 MeV
m_pi_MeV = m_p_MeV * float(K * eta)     # 139.0 MeV

sigma_piN = f_pi_MeV**2 * m_pi_MeV**2 / (m_p_MeV * (m_u_MeV + m_d_MeV) * m_p_MeV)
# This gives sigma_piN in a specific normalization; let me use the direct formula

# More directly: the QCD contribution is
# delta_m(QCD) = (m_d - m_u) / (m_u + m_d) * sigma_piN_eff
# where sigma_piN_eff ~ 2 * m_pi^2 * f_pi^2 / m_p ~ spectral

# Let me use the Gasser-Leutwyler leading-order result:
# delta_m(QCD) = (m_d - m_u) * B_0 * (1 + corrections)
# where B_0 = m_pi^2 / (m_u + m_d) from GMOR relation

B_0 = m_pi_MeV**2 / (m_u_MeV + m_d_MeV)  # chiral condensate parameter
print(f"  B_0 = m_pi^2 / (m_u + m_d) = {m_pi_MeV:.1f}^2 / {m_u_MeV + m_d_MeV:.2f}")
print(f"      = {B_0:.1f} MeV")

# The QCD contribution to n-p splitting (leading order):
# delta_m(QCD) = (m_d - m_u) * <p|ddbar-uubar|p> / (2*m_p)
# Using sigma_piN ~ 50 MeV (lattice/ChPT consensus):
# <p|uubar+ddbar|p> = 2*m_p*sigma_piN / (m_u+m_d)
# <p|ddbar-uubar|p> ~ <p|uubar+ddbar|p> * (m_d-m_u)/(m_d+m_u) * xi
# where xi ~ 1 (SU(2) breaking parameter)

# Simpler approach: use the Borsanyi et al. (2015) decomposition:
# delta_m(QCD) ~ 2.52 MeV (from lattice, using physical quark masses)
# delta_m(EM)  ~ -1.00 MeV (from lattice, using alpha)
# Total: ~1.51 MeV -> wait, PDG is 1.293. Let me use more standard values.

# Standard decomposition (FLAG/lattice average):
# delta_m = delta_m(QCD) + delta_m(EM)
# delta_m(QCD) = 2.52 +/- 0.17 MeV (from m_d-m_u in QCD)
# delta_m(EM) = -1.04 +/- 0.11 MeV (EM contribution, proton heavier)
# Total: 2.52 - 1.04 = 1.48 MeV (10% above PDG 1.293)
# (The precise value depends on how you split QCD/EM)

# In the spectral framework, let me compute both contributions:

# 1. QCD contribution: proportional to (m_d - m_u)
# The leading-order Dashen formula:
# delta_m(QCD) = (m_d - m_u) * m_pi^2 / (m_p * (m_u + m_d)) * (some nucleon matrix element)
# The matrix element is ~m_p/3 (constituent quark model)
# So: delta_m(QCD) ~ (m_d-m_u) * m_pi^2 / (3*(m_u+m_d))

delta_m_QCD = dm_quark * m_pi_MeV**2 / (3 * (m_u_MeV + m_d_MeV) * m_p_MeV) * m_p_MeV
# Simplify: = dm_quark * m_pi^2 / (3*(m_u+m_d))
delta_m_QCD_simple = dm_quark * m_pi_MeV**2 / (3 * (m_u_MeV + m_d_MeV))

# 2. EM contribution: alpha/lam1 * m_p was our old formula
# But more precisely: the EM self-energy difference is
# delta_m(EM) = -alpha * m_p * <EM structure>
# The simplest spectral estimate: -alpha * m_p / d1 (one EM correction per ghost mode)
# But this gives -alpha*m_p/d1 = -938.3/(137*6) = -1.14 MeV
delta_m_EM = -alpha * m_p_MeV / d1

# In D_wall/D_bulk: -alpha * m_p / (D_wall - 1)
# = -alpha * m_p / 6

print(f"""
  SPECTRAL DECOMPOSITION:

  1. QCD CONTRIBUTION (from spectral m_d - m_u):
     m_d - m_u = {dm_quark:.3f} MeV  [SPECTRAL]
     m_pi = {m_pi_MeV:.1f} MeV  [SPECTRAL]
     
     delta_m(QCD) = (m_d-m_u) * m_pi^2 / (3*(m_u+m_d))
                  = {dm_quark:.3f} * {m_pi_MeV:.1f}^2 / (3*{m_u_MeV+m_d_MeV:.2f})
                  = {delta_m_QCD_simple:.3f} MeV

  2. EM CONTRIBUTION (from spectral alpha):
     delta_m(EM) = -alpha * m_p / (D_wall-1) = -alpha * m_p / d1
                 = -{alpha:.6f} * {m_p_MeV:.1f} / {d1}
                 = {delta_m_EM:.3f} MeV

  TOTAL:
     delta_m = delta_m(QCD) + delta_m(EM)
             = {delta_m_QCD_simple:.3f} + ({delta_m_EM:.3f})
             = {delta_m_QCD_simple + delta_m_EM:.3f} MeV

  PDG: 1.2934 MeV
  Error: {abs(delta_m_QCD_simple + delta_m_EM - 1.2934)/1.2934*100:.1f}%
""")

delta_m_total = delta_m_QCD_simple + delta_m_EM

# =====================================================================
#  ALTERNATIVE: DIRECT SPECTRAL FORMULA
# =====================================================================
print(f"{'='*72}")
print("ALTERNATIVE: CAN WE FIND A CLEANER SPECTRAL FORMULA?")
print(f"{'='*72}")

# The old formula: delta_m = alpha/lam1 * m_p = 1.369 MeV (5.9% off)
delta_m_old = alpha / lam1 * m_p_MeV

# New idea: delta_m involves BOTH EM (alpha) and quark mass (m_d-m_u).
# In D_wall/D_bulk notation:
# delta_m = m_p * alpha * (D_bulk-D_wall) / (D_wall * D_bulk)
# = m_p * alpha * 4 / 77
delta_m_new1 = m_p_MeV * alpha * (D_bulk - D_wall) / (D_wall * D_bulk)

# Or: delta_m = m_p * alpha * (D_bulk-D_wall) / (p * D_bulk)
delta_m_new2 = m_p_MeV * alpha * (D_bulk - D_wall) / (p * D_bulk)

# Or: delta_m = m_pi * eta (= pion mass * spectral asymmetry)
delta_m_new3 = m_pi_MeV * float(eta)

# Or: delta_m = m_p * alpha / (D_bulk - D_wall + 1) = m_p * alpha / lam1
# That's the old formula.

# Let me try: delta_m = m_p * (m_d - m_u) / (m_d + m_u) * K
# This would give the isospin-breaking fraction times the Koide factor
iso_frac = dm_quark / (m_u_MeV + m_d_MeV)
delta_m_new4 = m_p_MeV * iso_frac * float(K)

# Or simply: delta_m as a D_wall/D_bulk ratio of m_pi
# delta_m / m_pi = some spectral number
ratio = 1.2934 / m_pi_MeV
print(f"  delta_m / m_pi = 1.2934 / {m_pi_MeV:.1f} = {ratio:.6f}")
print(f"  Compare: eta/p = {float(eta)/p:.6f} (= 2/27 = 0.0741)")
print(f"  Compare: 1/D_bulk = 1/{D_bulk} = {1/D_bulk:.6f}")
print(f"  Compare: alpha*p = {alpha*p:.6f}")
print(f"  Compare: 1/(p*(D_bulk-1)) = 1/30 = {1/(p*(D_bulk-1)):.6f}")

print(f"""
  FORMULA CANDIDATES:
  
  Old:     alpha/lam1 * m_p = {delta_m_old:.4f} MeV  ({abs(delta_m_old-1.2934)/1.2934*100:.1f}% off)
  
  Dashen:  QCD + EM decomposition = {delta_m_total:.4f} MeV  ({abs(delta_m_total-1.2934)/1.2934*100:.1f}% off)
  
  New 1:   m_p * alpha * (D_bulk-D_wall)/(D_wall*D_bulk)
           = m_p * alpha * 4/77 = {delta_m_new1:.4f} MeV  ({abs(delta_m_new1-1.2934)/1.2934*100:.1f}% off)
  
  New 2:   m_p * alpha * (D_bulk-D_wall)/(p*D_bulk)
           = m_p * alpha * 4/33 = {delta_m_new2:.4f} MeV  ({abs(delta_m_new2-1.2934)/1.2934*100:.1f}% off)
  
  New 3:   m_pi * eta = {m_pi_MeV:.1f} * 2/9
           = {delta_m_new3:.4f} MeV  ({abs(delta_m_new3-1.2934)/1.2934*100:.1f}% off)

  New 4:   m_p * (m_d-m_u)/(m_d+m_u) * K
           = {delta_m_new4:.4f} MeV  ({abs(delta_m_new4-1.2934)/1.2934*100:.1f}% off)

  PDG:     1.2934 MeV
""")

# Check which is best
candidates = [
    ("alpha/lam1 * m_p (old)", delta_m_old),
    ("Dashen QCD+EM", delta_m_total),
    ("m_p*alpha*4/77", delta_m_new1),
    ("m_p*alpha*4/33", delta_m_new2),
    ("m_pi * eta", delta_m_new3),
    ("m_p*(md-mu)/(md+mu)*K", delta_m_new4),
]

print(f"  RANKED BY ACCURACY:")
ranked = sorted(candidates, key=lambda x: abs(x[1]-1.2934))
for name, val in ranked:
    err = abs(val-1.2934)/1.2934*100
    print(f"    {name:<30} = {val:.4f} MeV  ({err:.1f}%)")

print(f"""

  BEST CANDIDATE: {ranked[0][0]} at {abs(ranked[0][1]-1.2934)/1.2934*100:.1f}%
""")

# =====================================================================
#  IMPACT ON NEUTRON LIFETIME
# =====================================================================
print(f"{'='*72}")
print("IMPACT: NEUTRON LIFETIME WITH BEST SPECTRAL delta_m")
print(f"{'='*72}")

from scipy import integrate

best_dm = ranked[0][1]
G_F = 1 / (np.sqrt(2) * v_GeV**2)
g_A = 1 + float(eta) + float(K)/(d1+lam1)
lambda_CKM = float(eta) * (1 + 0.1187/(3*PI))
V_ud_sq = 1 - lambda_CKM**2
hbar_GeV_s = 6.582119569e-25

for dm_label, dm_val in [("PDG", 1.2934), ("Best spectral", best_dm)]:
    q = dm_val / m_e_MeV
    def integrand(eps):
        if eps < 1 or eps > q: return 0.0
        return eps * np.sqrt(eps**2-1) * (q-eps)**2
    f_int, _ = integrate.quad(integrand, 1.0, q)
    f_corr = f_int * (1 + 0.03886)
    m_e_GeV = m_e_MeV / 1000
    rate = (G_F**2 * m_e_GeV**5) / (2*PI**3) * V_ud_sq * f_corr * (1+3*g_A**2)
    tau = hbar_GeV_s / rate
    print(f"  delta_m = {dm_val:.4f} MeV ({dm_label}): tau_n = {tau:.1f} s ({abs(tau-878.4)/878.4*100:.1f}%)")

print()
print("=" * 72)
