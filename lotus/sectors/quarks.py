"""
lotus.sectors.quarks — Quark Masses and Proton Mass
======================================================

Source: compile_universe.py lines 200-221
"""

import math
from lotus.constants import PI, M_E_GEV


def proton_mass_ratio(M, alpha: float, loops: int = 2) -> float:
    """Proton-to-electron mass ratio from Parseval fold energy.

    Bare: m_p/m_e = d1 * pi^5 = 6*pi^5
    Hurricane corrections: geometry probing its own curvature.
      G1 = 10/9 (1-loop EM)
      G2 = -280/9 (2-loop)

    Args:
        M: manifold
        alpha: fine-structure constant
        loops: 0, 1, or 2

    Returns:
        m_p / m_e
    """
    bare = M.d1 * PI ** 5
    if loops == 0:
        return bare
    G1 = M.G                              # 10/9
    G2 = -(M.d1 + M.lambda1) * M.G      # Spec says -280/9 but let's compute exactly
    # Actually, from compile_universe.py: c1 = 10/9, c2 = 280/9
    # G2_coeff = -(280/9) but that doesn't factor cleanly from invariants alone.
    # The verified formula: bare * (1 + G*a^2/pi - (280/9)*a^4/pi^2)
    c1 = M.G * (alpha ** 2 / PI)
    if loops == 1:
        return bare * (1 + c1)
    c2 = (280 / 9) * (alpha ** 4 / PI ** 2)
    return bare * (1 + c1 - c2)


def proton_mass_GeV(M, alpha: float) -> float:
    """Proton mass in GeV.

    Args:
        M: manifold
        alpha: fine-structure constant

    Returns:
        m_p in GeV
    """
    return M_E_GEV * proton_mass_ratio(M, alpha, loops=0)


def quark_masses(M, v: float, m_mu: float, m_tau: float) -> dict:
    """Six quark masses from Yukawa universality + Z_3 piercing depths.

    Up-type (χ₁, angular steps = π/3 per half-sector):
      σ_t = -1/120     surface + hurricane correction (1/(d₁·λ₁·2²))
      σ_c = -2π/3      one Z₃ sector deep (angular representation ω)
      σ_u = -π          1.5 sectors deep (angular representation ω²)

    Down-type (χ₂, spectral steps = G/p² = 10/81):
      σ_b = 77/90  = A + 1/(p²λ₁) = λ₁/d₁ + 1/45  (Wolfenstein A + KK)
      σ_s = -10/81 = -G/p²  (one spectral step: ghost coupling per sector²)
      σ_d = 2π/3 + G/p²     (C1 constraint: angular + spectral offset)

    All exponents are rationals in {d₁, λ₁, K, η, p} or multiples of π/3.
    Derivation: verification/quark_piercing_rg.py

    Args:
        M: manifold
        v: Higgs VEV in GeV
        m_mu: muon mass in GeV
        m_tau: tau mass in GeV

    Returns:
        dict with 't', 'c', 'u', 'b', 's', 'd' in GeV
    """
    sqrt2 = math.sqrt(2)

    return {
        't': v / sqrt2 * math.exp(-1 / 120),
        'c': v / sqrt2 * (m_mu / m_tau) * math.exp(-2 * PI / 3),
        'u': v / sqrt2 * (M_E_GEV / m_tau) * math.exp(-PI),
        'b': m_tau * math.exp(77 / 90),
        's': m_mu * math.exp(-10 / 81),
        'd': M_E_GEV * math.exp(2 * PI / 3 + M.G / M.p ** 2),
    }
