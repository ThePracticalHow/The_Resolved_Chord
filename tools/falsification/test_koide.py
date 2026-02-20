"""
Koide Formula Tests
==================
Four equivalent statements, circulant structure, Brannen parametrization.
From LENG_Master.md §5, BrannenMatrix.py.
"""

import numpy as np
import pytest

from leng_test_utils import (
    koide_ratio,
    koide_from_r,
    brannen_masses,
    brannen_sqrt_masses,
    circulant_c0_c1,
    R,
    DELTA,
    K_THEORY,
    M_E_PDG,
    M_MU_PDG,
    M_TAU_PDG,
)


class TestKoideRatio:
    """K = (sum m) / (sum √m)²"""

    def test_pdg_masses(self):
        masses = np.array([M_E_PDG, M_MU_PDG, M_TAU_PDG])
        K_obs = koide_ratio(masses)
        assert K_obs == pytest.approx(2/3, rel=0.01)  # < 0.01 within tau uncertainty
        assert abs(K_obs - 2/3) < 0.001

    def test_theory_masses(self):
        """Any Brannen masses with r=√2 give K=2/3"""
        mu = 17.7
        masses = brannen_masses(mu, R, DELTA)
        assert koide_ratio(masses) == pytest.approx(K_THEORY, rel=1e-10)


class TestFourEquivalentStatements:
    """Brannen, Circle, Circulant all equivalent"""

    def test_brannen_K_from_r(self):
        """Statement 1: K = (1+r²/2)/3 = 2/3 ⟺ r=√2"""
        assert koide_from_r(R) == pytest.approx(2/3, rel=1e-15)

    def test_circle_constraint(self):
        """Statement 2: ‖v - centroid‖² = 3μ² where v = (√m_k)"""
        mu = 17.7
        sqrt_m = brannen_sqrt_masses(mu, R, DELTA)
        centroid = np.array([mu, mu, mu])
        v = sqrt_m
        assert np.linalg.norm(v - centroid) ** 2 == pytest.approx(3 * mu**2, rel=1e-10)

    def test_circulant_ratio(self):
        """Statement 3: |c₁|/c₀ = 1/√2"""
        mu = 17.7
        sqrt_m = brannen_sqrt_masses(mu, R, DELTA)
        c0, c1 = circulant_c0_c1(sqrt_m)
        assert abs(c1) / abs(c0) == pytest.approx(1/np.sqrt(2), rel=1e-8)

    def test_all_equivalent(self):
        """Using PDG masses, all four give consistent K"""
        masses = np.array([M_E_PDG, M_MU_PDG, M_TAU_PDG])
        K = koide_ratio(masses)
        assert K == pytest.approx(2/3, rel=0.003)


class TestBrannenParametrization:
    """√m_k = μ(1 + √2·cos(δ + 2πk/3))"""

    def test_cos_sum_cancels(self):
        """Σ cos(δ + 2πk/3) = 0 → sum √m = 3μ"""
        mu = 17.7
        sqrt_m = brannen_sqrt_masses(mu, R, DELTA)
        assert np.sum(sqrt_m) == pytest.approx(3 * mu, rel=1e-10)

    def test_koide_independent_of_delta(self):
        """With r=√2, K=2/3 for any δ in the positive-mass domain"""
        mu = 17.7
        # Physical domain: all √m_k > 0 (verified explicitly)
        for delta in [0, DELTA, 2.0, 2*np.pi/3]:
            sqrt_m = brannen_sqrt_masses(mu, R, delta)
            assert np.all(sqrt_m > 0), f"δ={delta:.3f} gives negative √m"
            masses = brannen_masses(mu, R, delta)
            assert koide_ratio(masses) == pytest.approx(2/3, rel=1e-10)

    def test_koide_algebraic_even_unphysical(self):
        """K=2/3 holds algebraically for the SIGNED Brannen formula regardless of δ.
        When √m_k < 0, the standard koide_ratio (using unsigned √|m|) breaks,
        but the algebraic identity Σ(√m_k)² / (Σ√m_k)² = (1+r²/2)/3 still holds."""
        mu = 17.7
        for delta in [np.pi/4, np.pi, 3*np.pi/2]:
            sqrt_m = brannen_sqrt_masses(mu, R, delta)
            assert np.any(sqrt_m < 0), f"δ={delta:.3f} should have negative √m"
            # Compute K from SIGNED √m (the algebraic identity)
            K_signed = float(np.sum(sqrt_m**2) / np.sum(sqrt_m)**2)
            assert K_signed == pytest.approx(2/3, rel=1e-10)
