#!/usr/bin/env python3
"""
Quark Koide with Full SM 1-Loop Running (QCD + EW + Yukawa)
============================================================

Tests whether sector-separated Koide ratios K(u,c,t) and K(d,s,b)
approach 2/3 at high scales when Yukawa running is included.

Key physics:
- 1-loop QCD PRESERVES K exactly (universal gamma_m) — cannot test the prediction
- Top Yukawa (y_t ~ 1) BREAKS universality: m_t runs faster than m_c, m_u
- This is the only SM effect that changes K between scales
- The LENG prediction also requires KK threshold corrections at M_c (not included)

Jixiang Leng & Claude, February 2026
"""

import numpy as np
from scipy.integrate import solve_ivp
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from leng_test_utils import koide_ratio, K_THEORY

# ─────────────────────────────────────────────────────────────────────────────
#  CONSTANTS
# ─────────────────────────────────────────────────────────────────────────────

V_HIGGS = 246.22   # GeV
M_Z = 91.1876      # GeV
PI = np.pi

# MS-bar quark masses at M_Z (GeV) — from PDG 2024 review
# Obtained from reference-scale values via standard QCD running
M_AT_MZ = {
    'u': 0.00127,   # from m_u(2 GeV) = 2.16 MeV
    'd': 0.00270,   # from m_d(2 GeV) = 4.67 MeV
    's': 0.0535,    # from m_s(2 GeV) = 93.4 MeV
    'c': 0.619,     # from m_c(m_c) = 1.27 GeV
    'b': 2.85,      # from m_b(m_b) = 4.18 GeV
    't': 171.7,     # from m_t(pole) = 172.76 GeV
}

# Gauge couplings at M_Z
ALPHA_S_MZ = 0.1180
ALPHA_EM_MZ = 1.0 / 127.95
SIN2_TW = 0.23122

G3_SQ = 4 * PI * ALPHA_S_MZ                              # ~ 1.48
G2_SQ = ALPHA_EM_MZ * 4 * PI / SIN2_TW                   # ~ 0.424
G1_SQ = ALPHA_EM_MZ * 4 * PI / (1 - SIN2_TW) * (5.0/3)  # ~ 0.170 (GUT norm)

# Tau Yukawa (for Yukawa trace)
Y_TAU_SQ = 2 * (1.777 / V_HIGGS)**2


# ─────────────────────────────────────────────────────────────────────────────
#  SM 1-LOOP RGEs
# ─────────────────────────────────────────────────────────────────────────────

def sm_rges(t, state):
    """
    Full SM 1-loop RGE system.

    State: [g1^2, g2^2, g3^2, ln(m_u), ln(m_d), ln(m_s), ln(m_c), ln(m_b), ln(m_t)]
    t = ln(mu / M_Z)

    Includes:
    - Gauge coupling running (1-loop beta functions, n_g=3)
    - Quark mass running with QCD + EW + Yukawa anomalous dimensions
    - Top Yukawa is the key universality-breaking effect (~y_t^2 ~ 0.5)
    """
    g1sq, g2sq, g3sq = state[0:3]
    ln_m = state[3:9]
    masses = np.exp(ln_m)

    fac = 1.0 / (16 * PI**2)

    # Yukawa couplings squared: y_q^2 = 2 m_q^2 / v^2
    yq_sq = 2 * masses**2 / V_HIGGS**2

    # Yukawa trace: T = 3(y_u^2 + y_c^2 + y_t^2) + 3(y_d^2 + y_s^2 + y_b^2) + y_tau^2
    # Indices: u=0, d=1, s=2, c=3, b=4, t=5
    T = 3 * (yq_sq[0] + yq_sq[3] + yq_sq[5]) + \
        3 * (yq_sq[1] + yq_sq[2] + yq_sq[4]) + Y_TAU_SQ

    # --- Gauge coupling beta functions (1-loop, SM with 3 generations) ---
    b1 = 41.0 / 6    # U(1)_Y
    b2 = -19.0 / 6   # SU(2)_L
    b3 = -7.0         # SU(3)_C

    dg1sq = fac * b1 * g1sq**2
    dg2sq = fac * b2 * g2sq**2
    dg3sq = fac * b3 * g3sq**2

    # --- Quark mass anomalous dimensions (1-loop SM) ---
    # Up-type:   16pi^2 d(ln m)/dt = -8 g3^2 - 9/4 g2^2 - 17/12 g1^2 + 3 y_q^2 + T
    # Down-type: 16pi^2 d(ln m)/dt = -8 g3^2 - 9/4 g2^2 -  5/12 g1^2 + 3 y_q^2 + T

    gamma_qcd = 8 * g3sq
    gamma_ew_up = 9.0/4 * g2sq + 17.0/12 * g1sq
    gamma_ew_down = 9.0/4 * g2sq + 5.0/12 * g1sq

    dln_m = np.zeros(6)

    # Up-type: u(0), c(3), t(5)
    for i in [0, 3, 5]:
        dln_m[i] = fac * (-gamma_qcd - gamma_ew_up + 3*yq_sq[i] + T)

    # Down-type: d(1), s(2), b(4)
    for i in [1, 2, 4]:
        dln_m[i] = fac * (-gamma_qcd - gamma_ew_down + 3*yq_sq[i] + T)

    return np.concatenate([[dg1sq, dg2sq, dg3sq], dln_m])


# ─────────────────────────────────────────────────────────────────────────────
#  MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main():
    print("=" * 80)
    print("  QUARK KOIDE WITH FULL SM 1-LOOP RUNNING (QCD + EW + YUKAWA)")
    print("=" * 80)

    # Compute Yukawa couplings at M_Z
    print(f"\nMS-bar masses at M_Z = {M_Z} GeV:")
    for q, m in M_AT_MZ.items():
        y_sq = 2 * m**2 / V_HIGGS**2
        print(f"  m_{q}(M_Z) = {m:.4e} GeV    y_{q}^2 = {y_sq:.4e}")

    y_t_sq = 2 * M_AT_MZ['t']**2 / V_HIGGS**2
    print(f"\n  Top Yukawa dominance: y_t^2 = {y_t_sq:.4f}")
    print(f"  QCD gamma_m = 8*alpha_s = {8*ALPHA_S_MZ:.4f}")
    print(f"  Yukawa contribution to gamma_m(top) = 3*y_t^2/(16pi^2) = {3*y_t_sq/(16*PI**2):.4f}")
    print(f"  Ratio (Yukawa/QCD for top): {3*y_t_sq/(16*PI**2) / (8*ALPHA_S_MZ/(4*PI)):.1%}")

    # Koide ratios at M_Z
    m_uct = np.array([M_AT_MZ['u'], M_AT_MZ['c'], M_AT_MZ['t']])
    m_dsb = np.array([M_AT_MZ['d'], M_AT_MZ['s'], M_AT_MZ['b']])
    m_cbt = np.array([M_AT_MZ['c'], M_AT_MZ['b'], M_AT_MZ['t']])

    print(f"\nKoide ratios at M_Z:")
    print(f"  K(u,c,t) = {koide_ratio(m_uct):.4f}  (target: 2/3 = {K_THEORY:.4f})")
    print(f"  K(d,s,b) = {koide_ratio(m_dsb):.4f}  (target: 2/3 = {K_THEORY:.4f})")
    print(f"  K(c,b,t) = {koide_ratio(m_cbt):.4f}  (Motl-Rivero triplet)")

    # Initial state
    masses_mz = np.array([
        M_AT_MZ['u'], M_AT_MZ['d'], M_AT_MZ['s'],
        M_AT_MZ['c'], M_AT_MZ['b'], M_AT_MZ['t']
    ])
    y0 = np.concatenate([[G1_SQ, G2_SQ, G3_SQ], np.log(masses_mz)])

    # Run to high scales
    t_max = np.log(1e16 / M_Z)  # up to 10^16 GeV
    t_eval = np.linspace(0, t_max, 200)

    sol = solve_ivp(sm_rges, [0, t_max], y0,
                    t_eval=t_eval, method='RK45',
                    rtol=1e-10, atol=1e-14, max_step=0.5)

    if not sol.success:
        print(f"\nERROR: Integration failed: {sol.message}")
        return

    # Print results at selected scales
    target_scales = [M_Z, 200, 500, 1e3, 1e4, 1e6, 1e8, 1e10, 1e12, 1e14, 1e16]

    print(f"\n{'mu (GeV)':>12} {'alpha_s':>8} {'m_t/m_c':>9} "
          f"{'K(u,c,t)':>9} {'K(d,s,b)':>9} {'K(c,b,t)':>9}")
    print("-" * 62)

    results_at_targets = []
    for target in target_scales:
        t_target = np.log(target / M_Z)
        idx = np.argmin(np.abs(sol.t - t_target))
        state = sol.y[:, idx]

        g3sq = state[2]
        alpha_s = g3sq / (4 * PI)
        masses = np.exp(state[3:9])
        m_u, m_d, m_s, m_c, m_b, m_t = masses

        K_uct = koide_ratio(np.array([m_u, m_c, m_t]))
        K_dsb = koide_ratio(np.array([m_d, m_s, m_b]))
        K_cbt = koide_ratio(np.array([m_c, m_b, m_t]))

        print(f"{target:>12.2e} {alpha_s:>8.4f} {m_t/m_c:>9.1f} "
              f"{K_uct:>9.4f} {K_dsb:>9.4f} {K_cbt:>9.4f}")

        results_at_targets.append({
            'mu': target, 'alpha_s': alpha_s,
            'K_uct': K_uct, 'K_dsb': K_dsb, 'K_cbt': K_cbt,
            'm_t_over_m_c': m_t/m_c,
            'm_u': m_u, 'm_c': m_c, 'm_t': m_t,
            'm_d': m_d, 'm_s': m_s, 'm_b': m_b,
        })

    # Summary
    first = results_at_targets[0]
    last = results_at_targets[-1]

    print(f"\n{'='*80}")
    print(f"  SUMMARY: K evolution from M_Z to {last['mu']:.0e} GeV")
    print(f"{'='*80}")

    for label, key in [("K(u,c,t)", "K_uct"), ("K(d,s,b)", "K_dsb"), ("K(c,b,t)", "K_cbt")]:
        k0, k1 = first[key], last[key]
        dk = k1 - k0
        dist0 = abs(k0 - K_THEORY)
        dist1 = abs(k1 - K_THEORY)
        direction = "TOWARD 2/3" if dist1 < dist0 else "AWAY from 2/3"
        print(f"\n  {label}: {k0:.4f} -> {k1:.4f}  (Delta = {dk:+.4f}, {direction})")
        print(f"    Distance from 2/3: {dist0:.4f} -> {dist1:.4f}")

    print(f"\n  Mass hierarchy change:")
    print(f"    m_t/m_c: {first['m_t_over_m_c']:.1f} -> {last['m_t_over_m_c']:.1f}")

    # Sector-separated mass ratios at highest scale
    r = last
    print(f"\n  Masses at {last['mu']:.0e} GeV:")
    print(f"    Up-type:   m_u={r['m_u']:.4e}  m_c={r['m_c']:.4e}  m_t={r['m_t']:.4e}")
    print(f"    Down-type: m_d={r['m_d']:.4e}  m_s={r['m_s']:.4e}  m_b={r['m_b']:.4e}")

    # What ratios would be needed for K = 2/3
    print(f"\n  For reference: lepton mass ratios (= UV quark ratios if Yukawa universality holds):")
    print(f"    m_e/m_tau = {0.511/1776.86:.4e}")
    print(f"    m_mu/m_tau = {105.66/1776.86:.4f}")
    print(f"    m_u/m_t at {last['mu']:.0e}: {r['m_u']/r['m_t']:.4e}")
    print(f"    m_d/m_b at {last['mu']:.0e}: {r['m_d']/r['m_b']:.4e}")

    print(f"""
{'='*80}
  PHYSICS INTERPRETATION
{'='*80}

  1-loop QCD alone preserves K exactly (universal anomalous dimension).
  The Yukawa coupling y_t ~ {np.sqrt(y_t_sq):.2f} is the ONLY SM effect
  that changes K between scales: it makes m_t run faster than m_c or m_u.

  The LENG prediction (K -> 2/3 at compactification scale M_c) also
  requires KK threshold corrections from the S^5/Z_3 spectrum, which
  are NOT included in this SM-only calculation. This script shows what
  standard SM running can contribute; the rest must come from the
  geometric tower.

  Key question: does K move TOWARD or AWAY from 2/3 with Yukawa running?
  This tells us whether the SM and LENG predictions are compatible.
""")


if __name__ == "__main__":
    main()
