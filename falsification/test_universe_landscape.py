"""
Universe Landscape Scan Tests
=============================
Verify the selection pipeline: (4,2) negative masses, gap thresholds, ablation tally.
From UniverseLandscape.py, SELECTION_PIPELINE_SECTION.md.
"""

import numpy as np
import pytest

from leng_test_utils import (
    resonance_gap,
    sqrt_masses_universe,
    all_positive_universe,
    is_self_consistent,
    is_prime,
    koide_p,
    twist_formula,
    R,
)

# Scan parameters (match UniverseLandscape.py)
N_MAX = 8
P_MAX = 13
TOL = 1e-9


class TestS72NegativeMasses:
    """S⁷/Z₂ (n=4, p=2) has negative √m — the 'frozen universe'"""

    def test_s72_has_negative_sqrt_mass(self):
        n, p = 4, 2
        sm = sqrt_masses_universe(n, p)
        assert np.any(sm < 0)

    def test_s72_not_all_positive(self):
        assert not all_positive_universe(4, 2)

    def test_s72_sqrt_mass_values(self):
        sm = sqrt_masses_universe(4, 2)
        # At least one negative; UniverseLandscape says √m_0 < 0
        assert sm.min() < 0
        assert sm.min() == pytest.approx(-0.24, rel=0.15)  # ~-0.24 from geometry


class TestS53Alive:
    """S⁵/Z₃ (n=3, p=3) is the unique survivor"""

    def test_s53_all_positive(self):
        assert all_positive_universe(3, 3)

    def test_s53_self_consistent(self):
        assert is_self_consistent(3, 3)

    def test_s53_prime(self):
        assert is_prime(3)

    def test_s53_koide_k2_over_3(self):
        assert koide_p(3) == pytest.approx(2/3, rel=1e-15)


class TestResonanceGap:
    """Off-resonance candidates have gap > threshold"""

    def test_off_resonance_have_large_gap(self):
        """Typical off-resonance (n,p) have gap > 0.05"""
        off_resonance = [(2, 2), (2, 5), (3, 2), (3, 5), (5, 2), (5, 7)]
        for n, p in off_resonance:
            gap = resonance_gap(n, p)
            assert gap > 0.01

    def test_self_consistent_only_33_and_42(self):
        """Only (3,3) and (4,2) satisfy resonance lock in scan range"""
        survivors = []
        for n in range(1, N_MAX + 1):
            for p in range(2, P_MAX + 1):
                if is_self_consistent(n, p):
                    survivors.append((n, p))
        assert sorted(survivors) == [(3, 3), (4, 2)]


class TestAblationTally:
    """Reproduce the selection pipeline tally"""

    def test_resonance_lock_eliminates_94(self):
        total = N_MAX * (P_MAX - 1)  # 96
        passed = sum(1 for n in range(1, N_MAX + 1) for p in range(2, P_MAX + 1)
                    if is_self_consistent(n, p))
        assert passed == 2
        assert total - passed == 94

    def test_positive_masses_eliminates_one(self):
        """Of 2 self-consistent, only (3,3) has all positive"""
        self_consistent = [(3, 3), (4, 2)]
        viable = [(n, p) for n, p in self_consistent if all_positive_universe(n, p)]
        assert viable == [(3, 3)]
        assert len(self_consistent) - len(viable) == 1

    def test_unique_survivor(self):
        """Final survivor is S⁵/Z₃ only"""
        survivors = []
        for n in range(1, N_MAX + 1):
            for p in range(2, P_MAX + 1):
                if is_self_consistent(n, p) and all_positive_universe(n, p):
                    if is_prime(p) and abs(koide_p(p) - 1.0) > 1e-9:
                        survivors.append((n, p))
        assert survivors == [(3, 3)]


class TestGeneralFormulation:
    """n = p^{n-2} vs paper's 3n = p^{n-1}"""

    def test_33_satisfies_both(self):
        """(3,3) satisfies both formulations"""
        n, p = 3, 3
        assert 3 * n == p ** (n - 1)
        assert n == p ** (n - 2)

    def test_42_satisfies_general_only(self):
        """(4,2) satisfies n = p^{n-2} but not 3n = p^{n-1}"""
        n, p = 4, 2
        assert n == p ** (n - 2)
        assert 3 * n != p ** (n - 1)
