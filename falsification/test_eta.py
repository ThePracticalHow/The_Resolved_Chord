"""
Eta Invariant & Spectral Tests
==============================
Donnelly formula, twist = Σ|η|, landscape uniqueness.
From LENG_Master.md §7, EtaInvariant.py, eta_output.txt.
"""

import numpy as np
import pytest

from leng_test_utils import (
    eta_twisted_donnelly,
    reidemeister_torsion,
    d1_degeneracy,
    twist_formula,
    TWIST,
    TAU_R,
    D1,
)


class TestDonnellyEtaInvariant:
    """η_D(χ_m) for L(3;1,1,1)"""

    def test_eta1_magnitude(self):
        eta1 = eta_twisted_donnelly(3, 3, 1)
        assert abs(eta1) == pytest.approx(1/9, rel=1e-10)

    def test_eta2_conjugate(self):
        """η₂ = -η₁ in (i·cot)^n convention. Both real, opposite signs, |η|=1/9."""
        eta1 = eta_twisted_donnelly(3, 3, 1)
        eta2 = eta_twisted_donnelly(3, 3, 2)
        # Both essentially real in (i·cot)^n convention
        assert abs(eta1.imag) < 1e-14
        assert abs(eta2.imag) < 1e-14
        # η₁ = +1/9, η₂ = -1/9
        assert eta1.real == pytest.approx(1/9, rel=1e-10)
        assert eta2.real == pytest.approx(-1/9, rel=1e-10)
        # Equal magnitudes, opposite signs: η₁ + η₂ = 0
        assert abs(eta1) == pytest.approx(1/9, rel=1e-10)
        assert abs(eta2) == pytest.approx(1/9, rel=1e-10)
        assert abs(eta1 + eta2) < 1e-14

    def test_twist_sum(self):
        twist_computed = sum(abs(eta_twisted_donnelly(3, 3, m)) for m in range(1, 3))
        assert twist_computed == pytest.approx(2/9, rel=1e-10)
        assert twist_computed == pytest.approx(TWIST, rel=1e-15)


class TestEtaTorsionIdentity:
    """Σ|η_D| = d₁·τ_R for L(3;1,1,1)"""

    def test_s5_z3_identity(self):
        eta_sum = sum(abs(eta_twisted_donnelly(3, 3, m)) for m in range(1, 3))
        d1_tau = d1_degeneracy(3) * reidemeister_torsion(3, 3)
        assert eta_sum == pytest.approx(d1_tau, rel=1e-10)
        assert eta_sum == pytest.approx(2/9, rel=1e-15)


class TestLandscapeUniqueness:
    """Σ|η| = d₁·τ_R only for (n,p)=(3,3) among L(p;1,...,1)"""

    def test_s5_z3_only_match(self):
        matches = []
        for n in range(2, 6):
            for p in [2, 3, 5, 7, 11]:
                if p < 2:
                    continue
                d1_tau = d1_degeneracy(n) * reidemeister_torsion(p, n)
                eta_sum = sum(abs(eta_twisted_donnelly(p, n, m)) for m in range(1, p))
                if abs(eta_sum - d1_tau) < 1e-8:
                    matches.append((n, p))
        assert matches == [(3, 3)]

    @pytest.mark.parametrize("n,p", [(2, 2), (2, 3), (3, 2), (3, 5), (4, 3), (5, 3)])
    def test_others_dont_match(self, n, p):
        d1_tau = d1_degeneracy(n) * reidemeister_torsion(p, n)
        eta_sum = sum(abs(eta_twisted_donnelly(p, n, m)) for m in range(1, p))
        assert abs(eta_sum - d1_tau) > 1e-6
