"""
Geometry & Uniqueness Tests
===========================
Reidemeister torsion, uniqueness theorem, moment map, self-consistency.
From LENG_Master.md §2-5, §8, The_Resolved_Chord §2.
"""

import numpy as np
import pytest

from leng_test_utils import (
    reidemeister_torsion,
    d1_degeneracy,
    twist_formula,
    koide_from_r,
    R,
    TWIST,
    K_THEORY,
    D1,
    TAU_R,
    is_prime,
)


class TestReidemeisterTorsion:
    """τ_R = p^{-n} for L(p;1,...,1)"""

    def test_s5_z3_exact(self):
        assert reidemeister_torsion(3, 3) == pytest.approx(1/27, rel=1e-15)
        assert reidemeister_torsion(3, 3) == TAU_R

    @pytest.mark.parametrize("n,p", [(2, 2), (2, 3), (2, 5), (3, 2), (3, 4), (4, 2), (5, 3)])
    def test_formula_vs_product(self, n, p):
        """τ_R = 1 / Π|1-ω^k|^n via product formula"""
        omega = np.exp(2j * np.pi / p)
        product = 1.0
        for k in range(1, p):
            product *= abs(1 - omega ** k) ** n
        expected = 1.0 / product
        assert reidemeister_torsion(p, n) == pytest.approx(expected, rel=1e-10)

    def test_formula_p_negative_n(self):
        for n in range(2, 6):
            for p in range(2, 8):
                assert reidemeister_torsion(p, n) == pytest.approx(p ** (-n), rel=1e-12)


class TestUniquenessTheorem:
    """3n = p^{n-1} has unique prime solution n=p=3"""

    def test_n3_p3_satisfies(self):
        n, p = 3, 3
        assert 3 * n == p ** (n - 1)

    def test_unique_prime_solution(self):
        """Only (n,p)=(3,3) with p prime satisfies 3n = p^{n-1}"""
        solutions = []
        for n in range(1, 21):
            target = 3 * n
            for p in range(2, 101):
                if p ** (n - 1) == target and is_prime(p):
                    solutions.append((n, p))
        assert solutions == [(3, 3)]

    def test_n2_p6_not_prime(self):
        """n=2, p=6 satisfies but p is not prime"""
        assert 3 * 2 == 6 ** 1
        assert not is_prime(6)

    def test_n4_no_integer_p(self):
        """n=4: need p³=12, p not integer"""
        assert 12 ** (1/3) == pytest.approx(2.289, rel=0.01)


class TestMomentMap:
    """r = √2 from simplex geometry → K = 2/3"""

    def test_r_sqrt2(self):
        assert R == pytest.approx(np.sqrt(2), rel=1e-15)

    def test_K_from_r(self):
        assert koide_from_r(R) == pytest.approx(K_THEORY, rel=1e-15)
        assert koide_from_r(np.sqrt(2)) == pytest.approx(2/3, rel=1e-15)

    def test_K_equals_2_over_3(self):
        assert (1 + R**2 / 2) / 3 == pytest.approx(2/3, rel=1e-15)


class TestSelfConsistency:
    """p × twist = K when n=p=3"""

    def test_resonance_lock(self):
        p, twist, K = 3, TWIST, K_THEORY
        assert p * twist == pytest.approx(K, rel=1e-15)

    def test_twist_equals_d1_times_tau_R(self):
        assert D1 * TAU_R == pytest.approx(TWIST, rel=1e-15)
        assert twist_formula(3, 3) == pytest.approx(2/9, rel=1e-15)


class TestModeDegeneracy:
    """d₁ = 2n for S^{2n-1}"""

    def test_s5_d1_equals_6(self):
        assert d1_degeneracy(3) == 6
        assert d1_degeneracy(3) == D1

    @pytest.mark.parametrize("n", [2, 3, 4, 5])
    def test_d1_formula(self, n):
        assert d1_degeneracy(n) == 2 * n
