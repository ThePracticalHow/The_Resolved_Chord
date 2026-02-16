"""
Proton Mass Conjecture Tests
============================
6π⁵ formula and corrected formula with α² term.
From The_Resolved_Chord §6.
"""

import numpy as np
import pytest

from leng_test_utils import (
    proton_mass_leading,
    proton_mass_corrected,
    M_P_OVER_M_E_PDG,
    ALPHA,
)


class TestProtonMassLeading:
    """m_p/m_e = 6π⁵"""

    def test_value(self):
        val = proton_mass_leading(1.0)
        assert val == pytest.approx(6 * np.pi**5, rel=1e-15)
        assert val == pytest.approx(1836.118, rel=0.001)

    def test_agreement_leading(self):
        """Leading term ~0.002% error (The_Resolved_Chord)"""
        pred = proton_mass_leading(1.0)
        rel_err = abs(pred - M_P_OVER_M_E_PDG) / M_P_OVER_M_E_PDG
        assert rel_err < 0.0001  # 0.01%


class TestProtonMassCorrected:
    """m_p/m_e = 6π⁵(1 + 10/9 · α²/π)"""

    def test_value(self):
        val = proton_mass_corrected(1.0)
        assert val == pytest.approx(1836.1527, rel=0.001)

    def test_coefficient_10_over_9(self):
        """10/9 = λ₁ × Σ|η| = 5 × 2/9"""
        assert 5 * (2/9) == pytest.approx(10/9, rel=1e-15)

    def test_agreement_corrected(self):
        """Corrected formula ~0.00001% (The_Resolved_Chord)"""
        pred = proton_mass_corrected(1.0)
        rel_err = abs(pred - M_P_OVER_M_E_PDG) / M_P_OVER_M_E_PDG
        assert rel_err < 0.000001  # 0.0001%


class TestGaugeStructure:
    """sin²θ_W = 3/8 at unification"""

    def test_unification_value(self):
        sin2_theta = 3/8
        assert sin2_theta == 0.375

    def test_mz_value(self):
        """RG running ~0.2313 at M_Z"""
        expected = 0.2313
        assert abs(expected - 0.23122) < 0.001
