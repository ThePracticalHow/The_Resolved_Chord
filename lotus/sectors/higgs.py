"""
lotus.sectors.higgs — Higgs VEV and Mass
============================================

Source: compile_universe.py lines 264-272

Hurricane corrections:
  Tree:   v = m_p * (2/α − (d1+λ1+K))
  1-loop: v *= (1 + c_grav * (α/π)²)   [gravity curvature feedback]

  Tree:   m_H = m_p * (1/α − 7/2)
  1-loop: m_H *= (1 + G * α² / π)      [EM self-energy]
"""

import math
from lotus.constants import PI


def higgs_vev(M, m_p: float, alpha: float, loops: int = 2) -> float:
    """Higgs VEV from spectral geometry.

    v = m_p * (2/alpha - (d1 + lambda1 + K))
    Bulk stiffness Theorem: exact spectral content determines v.

    Hurricane (1-loop): gravity curvature feedback c_grav*(α/π)².

    Args:
        M: S5Z3 manifold
        m_p: proton mass in GeV
        alpha: fine-structure constant
        loops: 0 (tree), 1, or 2

    Returns:
        v in GeV
    """
    v_bare = m_p * (2 / alpha - (M.d1 + M.lambda1 + M.K))
    if loops == 0:
        return v_bare
    # 1-loop: gravity curvature feedback
    c1 = M.c_grav * (alpha / PI) ** 2
    return v_bare * (1 + c1)


def higgs_mass(M, m_p: float, alpha: float, loops: int = 2) -> float:
    """Higgs boson mass from EM budget minus ghost spectral cost.

    m_H = m_p * (1/alpha - 7/2)
    7/2 = half the Dirac eigenvalue = ghost cost.

    Hurricane (1-loop): EM self-energy G*α²/π.

    Args:
        M: S5Z3 manifold
        m_p: proton mass in GeV
        alpha: fine-structure constant
        loops: 0 (tree), 1, or 2

    Returns:
        m_H in GeV
    """
    m_H_bare = m_p * (1 / alpha - 3.5)
    if loops == 0:
        return m_H_bare
    # 1-loop: EM self-energy correction
    c1 = M.G * alpha ** 2 / PI
    return m_H_bare * (1 + c1)


def quartic_coupling(m_H: float, v: float) -> float:
    """Higgs quartic coupling.

    lambda_H = m_H^2 / (2*v^2)

    Returns:
        lambda_H ~ 0.1295
    """
    return m_H ** 2 / (2 * v ** 2)

