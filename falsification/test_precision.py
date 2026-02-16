"""
Precision / Numeric Stability Tests
====================================
Blocks: "Floating-point coincidence."
Run in double precision; optional high-precision check.
"""

import sys
import pytest

from leng_test_utils import (
    predict_masses_from_e,
    is_self_consistent,
    all_positive_universe,
    resonance_gap,
    koide_ratio,
    M_E_PDG,
)

try:
    from decimal import Decimal, getcontext
    HAS_DECIMAL = True
except ImportError:
    HAS_DECIMAL = False


class TestDoublePrecision:
    """Standard float64 (double) reproducibility."""

    def test_survivor_33_stable(self):
        """(3,3) consistently identified as survivor."""
        for _ in range(10):
            assert is_self_consistent(3, 3)
            assert all_positive_universe(3, 3)

    def test_masses_reproducible(self):
        """Mass predictions reproducible across calls."""
        m1 = predict_masses_from_e(M_E_PDG)
        m2 = predict_masses_from_e(M_E_PDG)
        assert m1 == m2

    def test_resonance_gap_exact_zero_for_33(self):
        """Resonance gap for (3,3) is numerically zero."""
        gap = resonance_gap(3, 3)
        assert gap < 1e-14  # Should be ~0 in double

    def test_koide_ratio_stable(self):
        """Koide ratio from predicted masses stable."""
        m_e, m_mu, m_tau = predict_masses_from_e(M_E_PDG)
        import numpy as np
        masses = np.array([m_e, m_mu, m_tau])
        K = koide_ratio(masses)
        assert abs(K - 2/3) < 1e-10


@pytest.mark.skipif(not HAS_DECIMAL, reason="decimal module required")
class TestHighPrecision:
    """Optional: higher precision gives same conclusions."""

    def test_resonance_exact_at_50_digits(self):
        """With Decimal(50 digits), resonance lock p×twist=K_p is exact for (3,3)."""
        getcontext().prec = 50
        n, p = 3, 3
        twist = Decimal(2) * n / (Decimal(p) ** n)
        K_p = Decimal(2) / p
        gap = abs(p * twist - K_p)
        assert gap < Decimal("1e-40")

    def test_no_spurious_survivor_at_high_prec(self):
        """(4,2) and (3,3) remain the only resonance-lock solutions."""
        getcontext().prec = 50
        solutions = []
        for n in range(2, 6):
            for p in range(2, 8):
                twist = Decimal(2) * n / (Decimal(p) ** n)
                K_p = Decimal(2) / p
                if abs(p * twist - K_p) < Decimal("1e-30"):
                    solutions.append((n, p))
        assert sorted(solutions) == [(3, 3), (4, 2)]
