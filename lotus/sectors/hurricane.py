"""
lotus.sectors.hurricane — Hurricane Coefficient Registry
==========================================================

Extensible registry of radiative correction coefficients.
All coefficients are spectral invariants of S^5/Z_3.

Grammar: rational combinations of {1, d1, lam1, K, eta, p}.
No external libraries — pure lotus + stdlib.
"""

# Known hurricane coefficients: observable -> [(loop, formula, value)]
# loop = 1 or 2, formula = spectral expression string
# Actual computation lives in sector modules (quarks, mixing, gauge, gravity).

HURRICANE_COEFFICIENTS = {
    'proton_mass_ratio': [
        (1, 'G = lam1*eta', 10 / 9),
        (2, '-280/9', -280 / 9),
    ],
    'cabibbo': [
        (1, '1/p', 1 / 3),
    ],
    'wolfenstein_A': [
        (1, '-eta', -2 / 9),
    ],
    'alpha_lag': [
        (1, 'G/p', 10 / 27),
    ],
    'alpha_s_split': [
        (1, '-d1', -6),
    ],
    'gravity': [
        (1, '-1/(d1*lam1)', -1 / 30),
    ],
}
