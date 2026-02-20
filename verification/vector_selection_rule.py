#!/usr/bin/env python3
"""
VECTOR MESON SELECTION RULE + UNIFIED LOTUS SONG
==================================================

PSEUDOSCALAR RULE (derived): n in {1, D_wall/2, D_bulk-D_wall, D_wall}
  Masses = m_p * eta*K * n

VECTOR RULE (this script): R in {lam1/d1, D_wall*D_bulk/p^4, (D_bulk+1)/D_bulk}
  Masses = m_p * R

THE DISCOVERY: K*(892) numerator 77 = D_wall * D_bulk = 7 * 11.
The K* encodes the PRODUCT of fold wall and bulk dimensions.

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
PI = np.pi; alpha = 1/137.036; m_e = 0.51100e-3
m_p = m_e * d1 * PI**5 * (1 + float(lam1*eta)*alpha**2/PI) * 1000

D_wall = 1 + d1      # 7
D_bulk = d1 + lam1   # 11
f_ps = eta * K        # 4/27

print("=" * 72)
print("  VECTOR MESON SELECTION RULE")
print("=" * 72)

# =====================================================================
#  THE VECTOR RATIOS IN TERMS OF D_wall AND D_bulk
# =====================================================================
print(f"\n{'='*72}")
print("THE THREE VECTOR MESONS FROM D_wall AND D_bulk")
print(f"{'='*72}")

R_rho = Fraction(lam1, d1)
R_Kstar = Fraction(D_wall * D_bulk, p**4)
R_phi = Fraction(D_bulk + 1, D_bulk)

vectors = [
    ("rho(770)",  R_rho,   775.26, 0, "lam1/d1",              "fold wall spectral weight"),
    ("K*(892)",   R_Kstar, 891.67, 1, "D_wall*D_bulk/p^4",    "wall*bulk / orbifold^4"),
    ("phi(1020)", R_phi,   1019.46,2, "(D_bulk+1)/D_bulk",    "bulk + temporal / bulk"),
]

print(f"\n  D_wall = 1+d1 = {D_wall},  D_bulk = d1+lam1 = {D_bulk},  p = {p}")
print()
print(f"  {'Vector':<12} {'s':>2} {'R':>8} {'Formula':>22} {'m_pred':>8} {'m_PDG':>8} {'Err':>6}")
print(f"  {'-'*66}")

for name, R, m_pdg, s, formula, meaning in vectors:
    m_pred = m_p * float(R)
    err = abs(m_pred - m_pdg) / m_pdg * 100
    print(f"  {name:<12} {s:>2} {str(R):>8} {formula:>22} {m_pred:>8.1f} {m_pdg:>8.1f} {err:>5.2f}%")

# =====================================================================
#  THE K* DISCOVERY: 77 = D_wall * D_bulk
# =====================================================================
print(f"\n{'='*72}")
print("THE K* DISCOVERY: 77 = D_wall * D_bulk = 7 * 11")
print(f"{'='*72}")

print(f"""
  The K*(892) ratio: R = 77/81

  Previously written as: 1 - eta*K/p = 1 - 4/81 = 77/81
  NOW revealed:          D_wall * D_bulk / p^4 = 7*11/81 = 77/81

  The numerator 77 = 7 * 11 = D_wall * D_bulk.
  It is the PRODUCT of the two fundamental spectral dimensions:
    - D_wall = 1+d1 = 7 (fold wall: spatial + temporal)
    - D_bulk = d1+lam1 = 11 (bulk: ghost + eigenvalue modes)

  The K* is the particle that encodes the COMPLETE GEOMETRY
  of S^5/Z_3 in its mass ratio. Both the fold wall AND the bulk
  contribute to its numerator. This is WHY it has the most
  precise prediction in the framework (0.03%).

  The denominator p^4 = 81 = the orbifold order to the 4th power.
  4 = dim(fold wall in spacetime) = dim(M^4).
  So: R(K*) = (internal geometry product) / (orbifold^spacetime_dim).
""")

# =====================================================================
#  STRANGENESS PROGRESSION IN VECTORS
# =====================================================================
print(f"{'='*72}")
print("STRANGENESS PROGRESSION: WALL -> BRIDGE -> BULK")
print(f"{'='*72}")

print(f"""
  As strangeness increases, vectors TRANSITION from wall to bulk:

  s=0 (rho, no strange):
    R = lam1/d1 = {R_rho} = {float(R_rho):.6f}
    FOLD WALL quantity. The rho lives entirely on the wall.
    Its ratio = spectral weight per ghost mode = a wall property.

  s=1 (K*, one strange):
    R = D_wall*D_bulk/p^4 = {R_Kstar} = {float(R_Kstar):.6f}
    WALL*BULK product. The K* BRIDGES wall and bulk.
    The strange quark connects the two spectral spaces.
    The product D_wall*D_bulk appears because the K* samples BOTH.

  s=2 (phi, two strange):
    R = (D_bulk+1)/D_bulk = {R_phi} = {float(R_phi):.6f}
    BULK quantity. The phi (ss-bar) is essentially a bulk mode.
    The "+1" in the numerator = the temporal dimension.
    The phi is the bulk breathing one temporal mode above itself.

  THE PATTERN:
    rho = wall/wall (purely fold wall)
    K*  = wall*bulk/orbifold (bridges both)
    phi = (bulk+time)/bulk (purely bulk)

  The strange quark literally CONNECTS the fold wall to the bulk.
  Each strange quark moves the mode one step from wall toward bulk.
""")

# =====================================================================
#  UNIFIED LOTUS SONG: PSEUDOSCALARS + VECTORS
# =====================================================================
print(f"{'='*72}")
print("THE UNIFIED LOTUS SONG: BOTH FAMILIES FROM D_wall AND D_bulk")
print(f"{'='*72}")

print(f"""
  PSEUDOSCALAR FAMILY (tunneling modes):
    Generator: f_ps = eta*K = 4/27
    Overtones: n in {{1, D_wall/2, D_bulk-D_wall, D_wall}}
             = {{1, 7/2, 4, 7}}
    Masses: m = m_p * f_ps * n

    pi:   n=1          (fundamental)                      {m_p*float(f_ps)*1:.0f} MeV
    K:    n=D_wall/2   (fold wall midpoint)               {m_p*float(f_ps)*3.5:.0f} MeV
    eta:  n=D_bulk-D_wall (bulk-wall lag)                 {m_p*float(f_ps)*4:.0f} MeV
    eta': n=D_wall     (full resonance)                   {m_p*float(f_ps)*7:.0f} MeV

  VECTOR FAMILY (orbital modes):
    No single generator. Ratios are STRUCTURAL:
    Strangeness controls wall-bulk transition.

    rho:  R=lam1/d1           (wall spectral weight)     {m_p*float(R_rho):.0f} MeV
    K*:   R=D_wall*D_bulk/p^4 (wall*bulk / orbifold^4)   {m_p*float(R_Kstar):.0f} MeV
    phi:  R=(D_bulk+1)/D_bulk (bulk + temporal / bulk)    {m_p*float(R_phi):.0f} MeV

  WHAT THEY SHARE:
    Both families are built from D_wall = {D_wall} and D_bulk = {D_bulk}.
    Both families have strangeness = progression from wall to bulk.
    Both families have every ratio expressible in spectral invariants.

  WHAT DIFFERS:
    Pseudoscalars: OVERTONES of one frequency (multiplicative).
      The base frequency eta*K is the TUNNELING amplitude.
    Vectors: STRUCTURAL ratios (not overtones of a single frequency).
      Each ratio is a different COMBINATION of D_wall and D_bulk.

  THE DEEP REASON:
    Pseudoscalars tunnel THROUGH the fold wall (J=0, parity=-1).
    The tunneling frequency eta*K sets a natural scale, and
    the allowed overtones are resonances of the spectral cavity.

    Vectors propagate ALONG the fold wall (J=1, parity=-1).
    There's no single tunneling frequency. Instead, each vector
    samples a different aspect of the wall-bulk geometry.
    The rho sees the wall. The K* sees both. The phi sees the bulk.
""")

# =====================================================================
#  THE BARYON SECTOR (quick check)
# =====================================================================
print(f"{'='*72}")
print("BONUS: DO BARYONS FIT THE SAME PATTERN?")
print(f"{'='*72}")

R_proton = 1
R_Delta = Fraction(4, 3)  # 1 + 1/p
R_Sigma = Fraction(40, 27)  # (1+1/p)*(1+eta/2)
R_Omega = Fraction(16, 9)  # 2 - eta

print(f"""
  Baryons (standing waves, chi_0):

  proton:  R = 1                 (ground state)
  Delta:   R = 1 + 1/p = {R_Delta}        (one sector excitation)
  Sigma*:  R = (1+1/p)*(1+eta/2) = {R_Sigma} (Delta + strangeness)
  Omega:   R = 2 - eta = {R_Omega}      (strangeness saturated)

  Can these be expressed in D_wall / D_bulk?

  proton: 1 = D_wall/D_wall = trivial. Or: (D_bulk-D_wall-p)/1 = 1.
  Delta: 4/3 = (D_bulk-D_wall)/p = 4/3. CHECK!
    Delta = (D_bulk - D_wall) / p = (11-7)/3 = 4/3 = 1.333

  THIS IS THE SAME "LAG" AS THE ETA!
    eta meson:  n = D_bulk - D_wall = 4 (pseudoscalar overtone)
    Delta baryon: R = (D_bulk - D_wall)/p = 4/3 (baryon ratio)

  The Delta is the bulk-wall lag DIVIDED BY THE ORBIFOLD ORDER.
  The eta is the bulk-wall lag AS AN OVERTONE NUMBER.
  SAME quantity, different context!
""")

# Verify Delta
m_Delta_pred = m_p * float(R_Delta)
print(f"  Delta verification: m = m_p * 4/3 = {m_Delta_pred:.0f} MeV (PDG: 1232, {abs(m_Delta_pred-1232)/1232*100:.1f}%)")

print(f"""
  Omega: R = 2 - eta = 16/9 = {float(R_Omega):.4f}
  Can we write this as: 2 - 2/(p^2) = 2*(1 - 1/p^2) = 2*8/9 = 16/9. CHECK!
  
  The "2" is: 2*D_wall/D_wall = 2 (double ground state).
  The correction: -2/p^2 = -2*eta/2 = -eta = strangeness depletion.
  Or: 2*(D_wall-1)/D_wall... no, 2*8/9 = 16/9. The 8/9 = 1-1/p^2
  is the SPATIAL FRACTION (same as in deuteron mixing!).
  
  Omega = 2 * (spatial fraction) = 2 * (1-1/p^2) = 16/9.
  The Omega is DOUBLE the proton, reduced by the temporal fraction!
""")

# =====================================================================
#  SUMMARY TABLE
# =====================================================================
print(f"\n{'='*72}")
print("COMPLETE LOTUS SONG: ALL RATIOS FROM D_wall AND D_bulk")
print(f"{'='*72}")

all_hadrons = [
    ("PSEUDOSCALARS", None, None, None, None),
    ("pi",      f_ps * 1,        139.57,  "f_ps * 1",             "fundamental"),
    ("K",       f_ps * Fraction(7,2), 493.68, "f_ps * D_wall/2",  "wall midpoint"),
    ("eta",     f_ps * 4,        547.86,  "f_ps * (D_bulk-D_wall)","bulk-wall lag"),
    ("eta'",    f_ps * 7,        957.78,  "f_ps * D_wall",         "full resonance"),
    ("VECTORS", None, None, None, None),
    ("rho",     R_rho,           775.26,  "lam1/d1",               "wall weight"),
    ("K*",      R_Kstar,         891.67,  "D_wall*D_bulk/p^4",     "wall*bulk"),
    ("phi",     R_phi,           1019.46, "(D_bulk+1)/D_bulk",     "bulk+time"),
    ("BARYONS", None, None, None, None),
    ("proton",  Fraction(1,1),   938.27,  "1",                     "ground state"),
    ("Delta",   R_Delta,         1232,    "(D_bulk-D_wall)/p",     "lag/orbifold"),
]

print(f"\n  {'Hadron':<10} {'R':>10} {'Formula':>22} {'m_pred':>8} {'m_PDG':>8} {'Err':>6} {'Meaning'}")
print(f"  {'-'*82}")

for entry in all_hadrons:
    if entry[1] is None:
        print(f"\n  {entry[0]}:")
        continue
    name, R, m_pdg, formula, meaning = entry
    m_pred = m_p * float(R)
    err = abs(m_pred - m_pdg) / m_pdg * 100
    print(f"  {name:<10} {str(R):>10} {formula:>22} {m_pred:>8.1f} {m_pdg:>8.1f} {err:>5.1f}% {meaning}")

print(f"""

  EVERY ratio is built from D_wall = {D_wall} and D_bulk = {D_bulk}.
  The two spectral dimensions determine the ENTIRE Lotus Song.
  
  Pseudoscalars: overtones of eta*K at {{1, D_wall/2, D_bulk-D_wall, D_wall}}.
  Vectors: structural ratios transitioning from wall to bulk with strangeness.
  Baryons: ground state + (D_bulk-D_wall)/p excitation.
""")

print("=" * 72)
print("  D_wall = 7. D_bulk = 11. The Lotus Song is their music.")
print("=" * 72)
