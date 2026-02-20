"""
Neutrino Sector Falsification Tests
===================================
Predictions for PMNS mixing angles, CP phase, and absolute mass scale
derived from Supplement VII: The Neutrino Sector.
"""

import numpy as np
import pytest

from leng_test_utils import (
    M_E_PDG,
    M_P_OVER_M_E_PDG,
    TWIST,
    K_THEORY,
    D1,
    LAMBDA_1,
    SIN2_THETA_13_PDG,
    SIN2_THETA_12_PDG,
    SIN2_THETA_23_PDG,
    DELTA_CP_PMNS_PDG,
    DM2_21_PDG,
    DM2_32_PDG,
    M3_MEAS_PDG,
)

class TestPMNSMixing:
    """PMNS Mixing Angle Predictions (Supplement VII, §3-5)"""

    def test_theta13_reactor(self):
        """Reactor angle: sin²θ₁₃ = (ηK)² = (4/27)² = 16/729"""
        pred = (TWIST * K_THEORY)**2
        assert pred == pytest.approx(16/729, rel=1e-9)
        # Numerical: 0.02194 vs 0.02200 (NuFIT 5.3)
        assert pred == pytest.approx(SIN2_THETA_13_PDG, abs=0.001)

    def test_theta12_solar(self):
        """Solar angle: sin²θ₁₂ = p/π² = 3/π² (spectral impedance mismatch)"""
        p_orb = 3  # orbifold order
        pred = p_orb / np.pi**2
        assert pred == pytest.approx(0.30396, rel=1e-4)
        # Numerical: 0.3040 vs 0.307 (NuFIT 5.3), 1.0% error
        assert pred == pytest.approx(SIN2_THETA_12_PDG, abs=0.01)

    def test_theta23_atmospheric(self):
        """Atmospheric angle: sin²θ₂₃ = d₁ / (d₁ + λ₁) = 6/11"""
        pred = D1 / (D1 + LAMBDA_1)
        assert pred == pytest.approx(6/11, rel=1e-9)
        # Numerical: 0.5455 vs 0.546 (NuFIT 5.3)
        assert pred == pytest.approx(SIN2_THETA_23_PDG, abs=0.01)

class TestNeutrinoCP:
    """Leptonic CP Phase (Supplement VII, §6)"""

    def test_delta_cp_pmns(self):
        """δ_CP = 3 * arctan(2π²/9) ≈ 196.5°"""
        gamma_deg = np.degrees(np.arctan(2 * np.pi**2 / 9))
        pred = 3 * gamma_deg
        assert pred == pytest.approx(196.5, abs=0.1)
        # Numerical: 196.5 vs 195 (NuFIT 5.3)
        assert pred == pytest.approx(DELTA_CP_PMNS_PDG, abs=10.0)

class TestNeutrinoMassScale:
    """Neutrino Mass Scaling and Inversion (Supplement VII, §8-9)"""

    def test_mass_inversion_principle(self):
        """Inversion: p * m_p² * m_ν₃ = m_e³ -> m_ν₃ = m_e / (108π¹⁰)"""
        m_e_ev = M_E_PDG * 1e6
        # Formulation 1: p * m_p² * m_v = m_e³
        # Using ratios to avoid precision issues: m_v = m_e / (p * (m_p/m_e)²)
        m3_pred = m_e_ev / (3 * M_P_OVER_M_E_PDG**2)
        
        # Formulation 2: m_e / (108 * pi^10)
        m3_pred_pi = m_e_ev / (108 * np.pi**10)
        
        assert m3_pred == pytest.approx(m3_pred_pi, rel=0.001)
        assert m3_pred == pytest.approx(0.05052, abs=0.001) # 50.52 meV
        
        # Compare with measurement derived from splits (Normal Ordering)
        assert m3_pred == pytest.approx(M3_MEAS_PDG, rel=0.01)

    def test_mass_squared_ratio(self):
        """Splitting ratio: Δm²₃₂ / Δm²₂₁ = d₁² - p = 33"""
        ratio_pred = D1**2 - 3
        assert ratio_pred == 33
        
        ratio_meas = DM2_32_PDG / DM2_21_PDG
        # Numerical: 33 vs 32.58 (NuFIT 5.3)
        assert ratio_pred == pytest.approx(ratio_meas, rel=0.02)

    def test_normal_hierarchy_check(self):
        """Point/Side/Face hierarchy implies m₁ ≈ 0 (Normal Hierarchy)"""
        # Normal Ordering implies dm2_32 > 0
        assert DM2_32_PDG > 0
        
    def test_sum_of_masses(self):
        """Σ m_ν ≈ 59.2 meV (Normal Ordering, m₁=0)"""
        m3_pred = (M_E_PDG * 1e6) / (3 * M_P_OVER_M_E_PDG**2)
        # d_m2_32 / d_m2_21 = 33 -> m3² - m2² = 33(m2² - m1²)
        # With m1=0: m3² - m2² = 33m2² -> m3² = 34m2² -> m2 = m3 / sqrt(34)
        m2_pred = m3_pred / np.sqrt(34)
        m1_pred = 0.0
        sum_pred = m1_pred + m2_pred + m3_pred
        
        assert sum_pred == pytest.approx(59.2e-3, abs=0.1e-3) # 59.2 meV
        
        # Euclid/DESI window constraint mentioned in Supp VII: 50-70 meV
        assert 50e-3 < sum_pred < 70e-3
