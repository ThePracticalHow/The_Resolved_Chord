"""
lotus.gravity — Gravity and the Gauge Hierarchy
===================================================

Source: compile_universe.py lines 312-325, gravity_theorem_proof.py

Hurricane correction (1-loop):
  Tree:   X = X_bare = (d1+λ1)²/p
  1-loop: X = X_bare * (1 + c_grav) = X_bare * (1 − 1/30)
"""

import math
from lotus.constants import PI
from lotus.sectors.gauge import compactification_scale


def planck_mass(M, loops: int = 2) -> float:
    """Planck mass from the 5-lock overdetermined proof.

    X_bare = (d1+lambda1)^2 / p = 121/3
    X_corr = X_bare * (1 + c_grav) = X_bare * (1 - 1/30) = 3509/90
    M_P    = M_c * X_corr^(7/2) * sqrt(pi^3/3)

    Five locks:
      1. X_bare  = (d1+lambda1)^2/p  (spectral content)
      2. c_grav  = -1/(d1*lambda1)   (hurricane correction)
      3. Power   = 7/2               (dim S^5 + 2 KK levels)
      4. sqrt(pi^3/3)               (volume normalization)
      5. M_c                         (compactification scale)

    Args:
        M: S5Z3 manifold
        loops: 0 (tree — no c_grav), 1 or 2

    Returns:
        M_P in GeV
    """
    M_c = compactification_scale(M)
    X_bare = (M.d1 + M.lambda1) ** 2 / M.p
    if loops == 0:
        X = X_bare
    else:
        X = X_bare * (1 - 1 / (M.d1 * M.lambda1))
    return M_c * X ** (7 / 2) * math.sqrt(PI ** 3 / 3)


def gauge_hierarchy(M, loops: int = 2) -> float:
    """Gauge hierarchy ratio M_P / M_c.

    The hierarchy is a spectral invariant. No fine-tuning.

    Args:
        M: S5Z3 manifold
        loops: 0 (tree), 1, or 2

    Returns:
        M_P / M_c
    """
    M_c = compactification_scale(M)
    return planck_mass(M, loops=loops) / M_c

