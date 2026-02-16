"""
Negative Controls
=================
Deliberately wrong models that should fail loudly.
Blocks: "Your pipeline always returns something that looks good."
"""

import numpy as np
import pytest

from leng_test_utils import (
    twist_formula,
    koide_ratio,
    brannen_masses,
    mu_from_electron,
    eta_twisted_donnelly,
    d1_degeneracy,
    reidemeister_torsion,
    R,
    M_E_PDG,
    M_MU_PDG,
    M_TAU_PDG,
    K_THEORY,
)


# ── Control 1: Wrong cyclic order (Z_p with p ≠ 3) ─────────────────────────────

def brannen_masses_wrong_cyclic(mu: float, r: float, delta: float, p: int) -> np.ndarray:
    """Brannen with wrong p: √m_k = μ(1 + r·cos(δ + 2πk/p)), p ≠ 3."""
    return np.array([mu * (1 + r * np.cos(delta + 2 * np.pi * k / p)) for k in range(p)])


# ── Control 2: Randomized twist (phase scaling) ───────────────────────────────

def twist_formula_randomized(n: int, p: int, scale: float) -> float:
    """twist × scale. Scale ≠ 1 breaks resonance."""
    return twist_formula(n, p) * scale


def delta_with_wrong_twist(n: int, p: int, twist_scale: float) -> float:
    """δ = 2π/p + twist×scale."""
    return 2 * np.pi / p + twist_formula_randomized(n, p, twist_scale)


# ── Control 3: Wrong r (not √2) ──────────────────────────────────────────────

def koide_from_r(r: float) -> float:
    """K = (1 + r²/2)/3."""
    return (1 + r ** 2 / 2) / 3


class TestNegativeControlWrongCyclic:
    """Z_p with p ≠ 3 should not reproduce lepton masses."""

    def test_z2_gives_wrong_koide(self):
        delta = 2 * np.pi / 2 + twist_formula(4, 2)  # S⁷/Z₂
        mu = 1.0
        masses = brannen_masses_wrong_cyclic(mu, R, delta, 2)
        masses = masses ** 2
        K = koide_ratio(masses)
        assert abs(K - K_THEORY) > 0.1

    def test_z5_gives_wrong_masses(self):
        """Z₅ would give 5 masses, not 3; ratios wrong."""
        delta = 2 * np.pi / 5 + 0.1
        mu = mu_from_electron(M_E_PDG, R, delta)
        masses = brannen_masses_wrong_cyclic(mu, R, delta, 5)
        masses = np.sort(masses ** 2)[:3]  # Take smallest 3
        # These won't match (m_e, m_μ, m_τ) from PDG
        # Just check Koide is off
        K = koide_ratio(masses)
        assert abs(K - K_THEORY) > 0.05


class TestNegativeControlWrongTwist:
    """Twist scaling breaks resonance and mass predictions."""

    def test_twist_scale_1_1_breaks(self):
        """Scale 1.1: no longer resonance lock."""
        n, p = 3, 3
        twist_wrong = twist_formula(n, p) * 1.1
        K_p = 2 / p
        gap = abs(p * twist_wrong - K_p)
        assert gap > 0.05

    def test_twist_scale_prediction_wrong(self):
        """Scale 0.9: predicted masses diverge from PDG."""
        delta = 2 * np.pi / 3 + twist_formula(3, 3) * 0.9
        mu = mu_from_electron(M_E_PDG, R, delta)
        masses = brannen_masses(mu, R, delta)
        masses = np.sort(masses)
        m_mu_wrong = masses[1]
        rel_err = abs(m_mu_wrong - M_MU_PDG) / M_MU_PDG
        assert rel_err > 0.01  # > 1% error


class TestNegativeControlWrongManifold:
    """Wrong lens space L(p;1,1,1) breaks eta identity."""

    def test_L5_eta_sum_not_twist(self):
        """L(5;1,1,1) has Σ|η| ≠ 2/9 = twist for S⁵/Z₃."""
        eta_sum_p5 = sum(abs(eta_twisted_donnelly(5, 3, m)) for m in range(1, 5))
        assert abs(eta_sum_p5 - 2/9) > 0.01

    def test_L2_eta_sum_not_twist(self):
        """L(2;1,1,1) has Σ|η| ≠ 2/9."""
        eta_sum_p2 = sum(abs(eta_twisted_donnelly(2, 3, m)) for m in range(1, 2))
        assert abs(eta_sum_p2 - 2/9) > 0.01


class TestNegativeControlWrongR:
    """r ≠ √2 breaks K = 2/3."""

    def test_r_1_3_wrong_koide(self):
        K = koide_from_r(1.3)
        assert abs(K - K_THEORY) > 0.02

    def test_r_1_5_wrong_koide(self):
        K = koide_from_r(1.5)
        assert abs(K - K_THEORY) > 0.02


class TestNegativeControlSummary:
    """Control experiments table: pass/fail, typical deviations."""

    def test_control_table_values(self):
        """Document typical deviations for Supplement C table."""
        controls = []

        # Wrong cyclic Z₂
        delta = 2 * np.pi / 2 + twist_formula(4, 2)
        masses = brannen_masses_wrong_cyclic(1.0, R, delta, 2) ** 2
        k_dev = abs(koide_ratio(masses) - K_THEORY)
        controls.append(("Z₂ not Z₃", k_dev, k_dev > 0.1))

        # Twist scale 0.9
        delta = 2 * np.pi / 3 + twist_formula(3, 3) * 0.9
        mu = mu_from_electron(M_E_PDG, R, delta)
        masses = np.sort(brannen_masses(mu, R, delta))
        mu_err = abs(masses[1] - M_MU_PDG) / M_MU_PDG
        controls.append(("Twist ×0.9", mu_err, mu_err > 0.01))

        # Wrong r
        k_dev_r = abs(koide_from_r(1.4) - K_THEORY)
        controls.append(("r=1.4 not √2", k_dev_r, k_dev_r > 0.01))

        for name, dev, fails in controls:
            assert fails, f"Control {name} should fail: deviation={dev}"
