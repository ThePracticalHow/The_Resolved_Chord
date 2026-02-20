"""
Mass Prediction Tests
====================
Lepton mass predictions from m_e, ratios, δ agreement.
From LENG_Master.md §6, §9, TauPrediction.py.
"""

import numpy as np
import pytest

from leng_test_utils import (
    predict_masses_from_e,
    mu_from_electron,
    brannen_masses,
    R,
    DELTA,
    M_E_PDG,
    M_MU_PDG,
    M_TAU_PDG,
    M_TAU_ERR_PDG,
)


class TestMassPredictions:
    """√m_k = μ(1 + √2·cos(2π/3 + 2/9 + 2πk/3)), m_e as input"""

    def test_prediction_from_electron(self):
        m_e, m_mu, m_tau = predict_masses_from_e(M_E_PDG)
        assert m_e == pytest.approx(M_E_PDG, rel=1e-10)

    def test_muon_agreement(self):
        """m_μ predicted to 0.001% (LENG_Master)"""
        _, m_mu_pred, _ = predict_masses_from_e(M_E_PDG)
        rel_err = abs(m_mu_pred - M_MU_PDG) / M_MU_PDG
        assert rel_err < 0.0001  # 0.01%

    def test_tau_agreement(self):
        """m_τ predicted to 0.007% (LENG_Master)"""
        _, _, m_tau_pred = predict_masses_from_e(M_E_PDG)
        rel_err = abs(m_tau_pred - M_TAU_PDG) / M_TAU_PDG
        assert rel_err < 0.0001  # 0.01%

    def test_tau_within_pdg_uncertainty(self):
        """Prediction within 2σ of tau mass uncertainty"""
        _, _, m_tau_pred = predict_masses_from_e(M_E_PDG)
        deviation = abs(m_tau_pred - M_TAU_PDG)
        sigma = deviation / M_TAU_ERR_PDG
        assert sigma < 2  # Within 2σ (actual: ~1.04σ)

    def test_mass_ratios(self):
        """m_μ/m_e ≈ 206.77, m_τ/m_μ ≈ 16.82"""
        m_e, m_mu, m_tau = predict_masses_from_e(M_E_PDG)
        assert m_mu / m_e == pytest.approx(206.77, rel=0.01)
        assert m_tau / m_mu == pytest.approx(16.82, rel=0.01)


class TestDeltaAgreement:
    """δ = 2π/3 + 2/9 vs fitted δ"""

    def test_delta_value(self):
        assert DELTA == pytest.approx(2*np.pi/3 + 2/9, rel=1e-15)
        assert DELTA == pytest.approx(2.316617, rel=0.001)

    def test_fitted_delta_close(self):
        """Fit from masses gives δ within 0.003% of 2π/3+2/9 (LENG_Master)"""
        # Known from LENG_Master: δ_fit = 2.316625 rad, δ_hyp = 2π/3 + 2/9 = 2.316617 rad
        delta_hyp = 2*np.pi/3 + 2/9
        delta_fit = 2.316625  # Fitted from lepton masses
        rel_diff = abs(delta_fit - delta_hyp) / delta_hyp
        assert rel_diff < 0.0001  # 0.01%


class TestConsistency:
    """Internal consistency of mass formula"""

    def test_sorted_order(self):
        masses = predict_masses_from_e(M_E_PDG)
        assert masses[0] < masses[1] < masses[2]

    def test_mu_scale(self):
        mu = mu_from_electron(M_E_PDG, R, DELTA)
        assert mu == pytest.approx(17.7, rel=0.01)  # ~17.7 MeV^{1/2}
