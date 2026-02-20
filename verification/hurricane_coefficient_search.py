#!/usr/bin/env python3
"""
Hurricane Coefficient Search — Automated Loop Discovery
========================================================

Searches for hurricane correction coefficients that improve residuals.
Uses ONLY: Python stdlib (math, fractions, itertools) + lotus.

Grammar: coefficients as rational combinations of spectral invariants
  {1, d1, lam1, K, eta, p} — no external libraries.

For each observable with residual above threshold, tries:
  X_corr = X_bare * (1 + c * (alpha/pi)^n)  for n=1,2,3 and c from grammar.

Reports candidates that improve the residual.

Author: Jixiang Leng
Date: February 2026
"""

import sys
import os
import math

# Ensure project root on path
_here = os.path.dirname(os.path.abspath(__file__))
_root = os.path.dirname(_here)
if _root not in sys.path:
    sys.path.insert(0, _root)

from lotus import Universe
from lotus.core.geometry import S5Z3

# Measured values (from lotus.pdg_2024_reference — our project data)
PDG = {
    'proton_ratio': 1836.153,
    'inv_alpha': 137.036,
    'alpha_s': 0.1180,
    'sin2_weinberg': 0.23122,
    'higgs_vev': 246.22,
    'higgs_mass': 125.25,
    'ckm_lambda': 0.22500,
    'ckm_A': 0.826,
    'ckm_rho_bar': 0.1592,
    'ckm_eta_bar': 0.3490,
    'M_Planck_e19': 1.221,
    'CC_meV': 2.25,
    'eta_B': 6.1e-10,
}


def grammar_coefficients(M):
    """Rational coefficients from spectral invariants {1, d1, lam1, K, eta, p}."""
    d1 = M.d1
    lam1 = M.lambda1
    K = M.K
    eta = M.eta
    p = M.p
    G = lam1 * eta  # 10/9
    return [
        (1, "1"),
        (-1, "-1"),
        (1 / p, "1/p"),
        (-1 / p, "-1/p"),
        (1 / d1, "1/d1"),
        (-1 / d1, "-1/d1"),
        (1 / lam1, "1/lam1"),
        (eta, "eta"),
        (-eta, "-eta"),
        (K, "K"),
        (-K, "-K"),
        (G, "G"),
        (-G, "-G"),
        (G / p, "G/p"),
        (-G / p, "-G/p"),
        (1 / (d1 * lam1), "1/(d1*lam1)"),
        (-1 / (d1 * lam1), "-1/(d1*lam1)"),
        (lam1 / d1, "lam1/d1"),
        (-lam1 / d1, "-lam1/d1"),
        (d1 + lam1, "d1+lam1"),
        (-(d1 + lam1), "-(d1+lam1)"),
        (280 / 9, "280/9"),
        (-280 / 9, "-280/9"),
    ]


def residual(pred, meas):
    """Relative error in percent."""
    if meas == 0:
        return 0.0
    return abs(pred - meas) / abs(meas) * 100


def search_observable(name, X_bare, X_meas, alpha, M, threshold_pct=0.5):
    """Search for hurricane corrections that improve residual."""
    PI = math.pi
    a_over_pi = alpha / PI
    base_err = residual(X_bare, X_meas)
    if base_err < threshold_pct:
        return []
    coeffs = grammar_coefficients(M)
    improvements = []
    for c, c_str in coeffs:
        for n in (1, 2, 3):
            X_corr = X_bare * (1 + c * (a_over_pi ** n))
            err = residual(X_corr, X_meas)
            if err < base_err and err < base_err - 0.001:
                improvements.append({
                    'name': name,
                    'c': c,
                    'c_str': c_str,
                    'n': n,
                    'X_bare': X_bare,
                    'X_corr': X_corr,
                    'X_meas': X_meas,
                    'err_before': base_err,
                    'err_after': err,
                })
    return sorted(improvements, key=lambda x: x['err_after'])


def main():
    print("=" * 70)
    print("  HURRICANE COEFFICIENT SEARCH (stdlib + lotus only)")
    print("=" * 70)

    u = Universe()
    M = u.manifold
    alpha = u.alpha
    PI = math.pi

    observables = [
        ('proton_ratio', u.proton_ratio, PDG['proton_ratio'], M.d1 * PI ** 5),
        ('1/alpha', 1 / u.alpha, PDG['inv_alpha'], None),
        ('alpha_s', u.alpha_s, PDG['alpha_s'], None),
        ('sin2_theta_W', u.sin2_weinberg, PDG['sin2_weinberg'], None),
        ('higgs_vev', u.higgs_vev, PDG['higgs_vev'], None),
        ('higgs_mass', u.higgs_mass_val, PDG['higgs_mass'], None),
        ('CKM_lambda', u.ckm['lambda'], PDG['ckm_lambda'], M.eta),
        ('CKM_A', u.ckm['A'], PDG['ckm_A'], M.lambda1 / M.d1),
        ('M_Planck/1e19', u.M_Planck / 1e19, PDG['M_Planck_e19'], None),
    ]

    print(f"\n  Alpha = {alpha:.6f},  alpha/pi = {alpha/PI:.6f}")
    print(f"\n  Searching observables with residual > 0.5%...\n")

    all_improvements = []
    for name, pred, meas, bare in observables:
        if bare is None:
            continue
        imps = search_observable(name, bare, meas, alpha, M, threshold_pct=0.1)
        all_improvements.extend(imps)

    if not all_improvements:
        print("  No improvements found above threshold.")
        print("  (Current hurricane coefficients may already be optimal.)")
        return

    print("  Top candidates (correction improves residual):\n")
    for i, imp in enumerate(all_improvements[:15]):
        print(f"  [{i+1}] {imp['name']}")
        print(f"      c = {imp['c_str']},  n = {imp['n']}")
        print(f"      err: {imp['err_before']:.3f}% -> {imp['err_after']:.3f}%")
        print(f"      X_corr = {imp['X_corr']:.6f},  meas = {imp['X_meas']}")
        print()

    print("=" * 70)


if __name__ == '__main__':
    main()
