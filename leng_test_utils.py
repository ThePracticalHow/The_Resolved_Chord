"""
LENG Test Utilities
===================
Shared constants, formulas, and helpers for the LENG testing framework.
All values and formulas sourced from LENG_Master.md and The_Resolved_Chord.
"""

import numpy as np
from math import comb
from typing import Tuple, List

# ─────────────────────────────────────────────────────────────────────────────
#  EXPERIMENTAL CONSTANTS (PDG 2024)
# ─────────────────────────────────────────────────────────────────────────────

# Charged lepton masses (MeV)
M_E_PDG = 0.51099895
M_MU_PDG = 105.6583755
M_TAU_PDG = 1776.86
M_TAU_ERR_PDG = 0.12  # MeV
M_TAU_BELLE2_ERR = 0.04  # Belle II target precision

# Proton/electron ratio
M_P_OVER_M_E_PDG = 1836.15267343
M_P_OVER_M_E_ERR = 0.00000011

# Fine structure constant
ALPHA = 1 / 137.035999084  # CODATA 2018

# Electroweak mixing angle at M_Z
SIN2_THETA_W_MZ_PDG = 0.23122

# Neutrino Oscillation Parameters (NuFIT 5.3, 2024, Normal Ordering)
SIN2_THETA_13_PDG = 0.02200
SIN2_THETA_13_ERR = 0.00069
SIN2_THETA_12_PDG = 0.307
SIN2_THETA_12_ERR = 0.013
SIN2_THETA_23_PDG = 0.546
SIN2_THETA_23_ERR = 0.021
DELTA_CP_PMNS_PDG = 195.0  # degrees
DELTA_CP_PMNS_ERR = 50.0   # degrees

# Mass Splittings (eV^2)
DM2_21_PDG = 7.53e-5
DM2_32_PDG = 2.453e-3
M3_MEAS_PDG = 0.05028  # eV (derived from splittings)

# ─────────────────────────────────────────────────────────────────────────────
#  THEORY CONSTANTS (from S⁵/Z₃ geometry)
# ─────────────────────────────────────────────────────────────────────────────

R = np.sqrt(2)  # Brannen amplitude, moment map simplex side
TWIST = 2 / 9
DELTA = 2 * np.pi / 3 + 2 / 9  # 2π/3 + 2/9 rad
K_THEORY = 2 / 3
D1 = 6  # l=1 mode degeneracy on S⁵
TAU_R = 1 / 27  # Reidemeister torsion L(3;1,1,1)
LAMBDA_1 = 5  # First eigenvalue on S⁵: l(l+4)|_{l=1}

# ─────────────────────────────────────────────────────────────────────────────
#  CORE FORMULAS
# ─────────────────────────────────────────────────────────────────────────────


def reidemeister_torsion(p: int, n: int) -> float:
    """τ_R for L(p;1,...,1) = S^{2n-1}/Z_p. Cheeger-Müller: τ_R = p^{-n}."""
    return p ** (-n)


def d1_degeneracy(n: int) -> int:
    """First excited mode degeneracy on S^{2n-1}."""
    return 2 * n


def twist_formula(n: int, p: int) -> float:
    """twist = d₁ × τ_R = 2n / p^n"""
    return 2 * n / (p ** n)


def koide_ratio(masses: np.ndarray) -> float:
    """K = (m_e + m_μ + m_τ) / (√m_e + √m_μ + √m_τ)²"""
    m = np.abs(np.asarray(masses))
    return float(np.sum(m) / np.sum(np.sqrt(m)) ** 2)


def koide_from_r(r: float) -> float:
    """K = (1 + r²/2)/3. For r=√2: K=2/3."""
    return (1 + r ** 2 / 2) / 3


def brannen_sqrt_masses(mu: float, r: float, delta: float, k_vals: List[int] = None) -> np.ndarray:
    """√m_k = μ(1 + r·cos(δ + 2πk/3)), k=0,1,2"""
    if k_vals is None:
        k_vals = [0, 1, 2]
    return np.array([mu * (1 + r * np.cos(delta + 2 * np.pi * k / 3)) for k in k_vals])


def brannen_masses(mu: float, r: float, delta: float, k_vals: List[int] = None) -> np.ndarray:
    """m_k = [μ(1 + r·cos(δ + 2πk/3))]²"""
    return brannen_sqrt_masses(mu, r, delta, k_vals) ** 2


def mu_from_electron(m_e: float, r: float, delta: float) -> float:
    """Solve μ from m_e = [μ(1 + r·cos(δ))]² → μ = √m_e / (1 + r·cos(δ))"""
    return np.sqrt(m_e) / (1 + r * np.cos(delta))


def predict_masses_from_e(m_e: float) -> Tuple[float, float, float]:
    """Predict m_μ, m_τ from m_e using LENG formula. Returns (m_e, m_μ, m_τ)."""
    mu = mu_from_electron(m_e, R, DELTA)
    masses = brannen_masses(mu, R, DELTA)
    return tuple(sorted(masses))


def eta_twisted_donnelly(p: int, n: int, m: int) -> complex:
    """η_D(χ_m) for L(p;1,...,1) via Donnelly (1978): (1/p) Σ ω^{mk}·(i·cot(πk/p))^n"""
    omega = np.exp(2j * np.pi / p)
    total = 0j
    for k in range(1, p):
        cot_k = np.cos(np.pi * k / p) / np.sin(np.pi * k / p)
        total += omega ** (m * k) * (1j * cot_k) ** n
    return total / p


def circulant_c0_c1(sqrt_masses: np.ndarray) -> Tuple[complex, complex]:
    """Fourier components of Z₃ circulant: c₀ = μ, c₁ = (1/3) Σ_k √m_k · ω^{-k}"""
    omega = np.exp(-2j * np.pi / 3)
    c0 = np.mean(sqrt_masses)
    c1 = (1/3) * sum(sqrt_masses[k] * (omega ** k) for k in range(3))
    return c0, c1


def harmonic_dim_s5(a: int, b: int) -> int:
    """Dimension of harmonic bidegree (a,b) on S⁵."""
    if a >= 1 and b >= 1:
        return comb(a + 2, 2) * comb(b + 2, 2) - comb(a + 1, 2) * comb(b + 1, 2)
    elif b == 0:
        return comb(a + 2, 2)
    else:
        return comb(b + 2, 2)


def proton_mass_leading(m_e: float) -> float:
    """m_p/m_e = 6π⁵ (leading)"""
    return 6 * np.pi ** 5


def proton_mass_corrected(m_e: float) -> float:
    """m_p/m_e = 6π⁵(1 + 10/9 · α²/π)"""
    return 6 * np.pi ** 5 * (1 + (10 / 9) * (ALPHA ** 2 / np.pi))


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    for i in range(2, int(n ** 0.5) + 1):
        if n % i == 0:
            return False
    return True


# ─────────────────────────────────────────────────────────────────────────────
#  UNIVERSE LANDSCAPE (UniverseLandscape.py compatibility)
# ─────────────────────────────────────────────────────────────────────────────

def koide_p(p: int) -> float:
    """K_p = 2/p for p generations with r=√2."""
    return 2.0 / p


def resonance_gap(n: int, p: int) -> float:
    """|p × twist - K_p|. Self-consistent when gap < 1e-9."""
    return abs(p * twist_formula(n, p) - koide_p(p))


def sqrt_masses_universe(n: int, p: int) -> np.ndarray:
    """√m_k for universe (n,p): 1 + r·cos(δ + 2πk/p), δ = 2π/p + twist."""
    delta = 2 * np.pi / p + twist_formula(n, p)
    return np.array([1.0 + R * np.cos(delta + 2 * np.pi * k / p) for k in range(p)])


def all_positive_universe(n: int, p: int) -> bool:
    """All √m_k > 0 for universe (n,p)."""
    return bool(np.all(sqrt_masses_universe(n, p) > 0))


def is_self_consistent(n: int, p: int, tol: float = 1e-9) -> bool:
    """Resonance lock: p × twist = K_p."""
    return resonance_gap(n, p) < tol
