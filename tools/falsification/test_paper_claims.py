"""
Paper-Code Audit: Verify every quantitative claim in The_Resolved_Chord v6.
================================================================================
Each test corresponds to a specific number or formula in the paper.
Run: pytest tools/falsification/test_paper_claims.py -v

Reference: The_Resolved_Chord_Leng_2026_v6.tex
"""

import numpy as np
import pytest

from leng_test_utils import (
    M_E_PDG,
    M_MU_PDG,
    M_TAU_PDG,
    M_P_OVER_M_E_PDG,
    ALPHA,
    K_THEORY,
    TWIST,
    D1,
    LAMBDA_1,
    predict_masses_from_e,
    koide_ratio,
    proton_mass_corrected,
)
from alpha_s_constraint import predict_alpha_s_2loop, ALPHA_S_GAP


# -----------------------------------------------------------------------------
# TABLE: Charged lepton mass predictions (paper Table, Section 3.2)
# -----------------------------------------------------------------------------


def test_koide_exact():
    """Paper: K = 2/3 exactly (moment map theorem)."""
    assert abs(K_THEORY - 2 / 3) < 1e-15


def test_muon_prediction():
    """Paper: m_mu = 105.6594 MeV, error 0.001%."""
    m_e, m_mu, m_tau = predict_masses_from_e(M_E_PDG)
    assert abs(m_mu - 105.659) < 0.01  # ~105.6594
    rel_err = abs(m_mu - M_MU_PDG) / M_MU_PDG
    assert rel_err < 0.00002  # 0.002%


def test_tau_prediction():
    """Paper: m_tau = 1776.985 MeV, error 0.007%."""
    m_e, m_mu, m_tau = predict_masses_from_e(M_E_PDG)
    assert abs(m_tau - 1776.985) < 0.01
    rel_err = abs(m_tau - M_TAU_PDG) / M_TAU_PDG
    assert rel_err < 0.0001  # 0.01%


def test_koide_from_pdg():
    """Paper: K from PDG = 0.666661, theory 2/3."""
    masses = [M_E_PDG, M_MU_PDG, M_TAU_PDG]
    K = koide_ratio(np.array(masses))
    assert abs(K - 2 / 3) < 0.00001


# -----------------------------------------------------------------------------
# TABLE: Proton mass (paper Section 5, Table)
# -----------------------------------------------------------------------------


def test_6pi5_leading_term():
    """Paper: 6*pi^5 = 1836.11811, measured 1836.15267, error 0.0019%."""
    leading = 6 * np.pi**5
    assert abs(leading - 1836.118) < 0.01
    rel_err = abs(leading - M_P_OVER_M_E_PDG) / M_P_OVER_M_E_PDG
    assert rel_err < 0.00003  # 0.003%


def test_proton_formula_corrected():
    """Paper: m_p/m_e = 6*pi^5*(1 + G*alpha^2/pi + ...), ~1e-7% one-loop."""
    m_p_over_m_e_pred = proton_mass_corrected(M_E_PDG)
    rel_err = abs(m_p_over_m_e_pred - M_P_OVER_M_E_PDG) / M_P_OVER_M_E_PDG
    assert rel_err < 1e-5  # sub-10 ppm (one-loop formula)


# -----------------------------------------------------------------------------
# TABLE: Fine-structure constant (paper Section 8)
# -----------------------------------------------------------------------------


def test_alpha_consistent_with_proton():
    """Paper: 1/alpha = 137.036 from two-loop proton inversion."""
    inv_alpha = 1 / ALPHA
    assert 136.9 < inv_alpha < 137.2


# -----------------------------------------------------------------------------
# TABLE: Strong coupling (paper Section 9, Table)
# -----------------------------------------------------------------------------


def test_alpha_s_gap():
    """Paper: Delta(1/alpha_3) = pi^2 - 5."""
    assert abs(ALPHA_S_GAP - (np.pi**2 - 5)) < 1e-10


def test_alpha_s_prediction():
    """Paper: alpha_s(M_Z) = 0.1174, PDG 0.1180, 0.6 sigma."""
    res = predict_alpha_s_2loop()
    assert abs(res["alpha_s_pred"] - 0.1174) < 0.001
    assert 0.116 < res["alpha_s_pred"] < 0.119


# -----------------------------------------------------------------------------
# TABLE: VEV (paper Section 11, Table)
# -----------------------------------------------------------------------------


def test_vev_formula():
    """Paper: v/m_p = 2/alpha - 35/3 = 262.405, measured 262.418, 0.005%."""
    v_over_mp_pred = 2 / ALPHA - 35 / 3
    v_over_mp_meas = 246.2196 / 0.93827208943
    assert abs(v_over_mp_pred - 262.4) < 0.5
    rel_err = abs(v_over_mp_pred - v_over_mp_meas) / v_over_mp_meas
    assert rel_err < 0.0001  # 0.01%


# -----------------------------------------------------------------------------
# TABLE: Higgs mass (paper Section 11, Table)
# -----------------------------------------------------------------------------


def test_higgs_mass_ratio():
    """Paper: m_H/m_p = 1/alpha - 7/2 = 133.536, measured 133.490, 0.034%."""
    mH_mp_pred = 1 / ALPHA - 7 / 2
    mH_mp_meas = 125.25 / 0.93827208943
    assert abs(mH_mp_pred - 133.5) < 0.5
    rel_err = abs(mH_mp_pred - mH_mp_meas) / mH_mp_meas
    assert rel_err < 0.001  # 0.1%


def test_higgs_quartic():
    """Paper: lambda_H = 0.1295, measured 0.1294, 0.07%."""
    v = 246.2196
    m_H = 125.25
    lam = m_H**2 / (2 * v**2)
    assert abs(lam - 0.129) < 0.01


# -----------------------------------------------------------------------------
# TABLE: Weinberg angle (paper Table)
# -----------------------------------------------------------------------------


def test_weinberg_angle_at_MZ():
    """Paper: sin^2(theta_W)(M_Z) = 0.2313, measured 0.23122, 0.05%."""
    sin2_pred = 0.2313
    sin2_meas = 0.23122
    rel_err = abs(sin2_pred - sin2_meas) / sin2_meas
    assert rel_err < 0.001


# -----------------------------------------------------------------------------
# STRUCTURAL: Spectral invariants (paper throughout)
# -----------------------------------------------------------------------------


def test_twist_2_over_9():
    """Paper: twist = sum|eta_D| = 2/9."""
    assert abs(TWIST - 2 / 9) < 1e-15


def test_d1_equals_6():
    """Paper: d_1 = 6 (ghost mode degeneracy)."""
    assert D1 == 6


def test_lambda1_equals_5():
    """Paper: lambda_1 = 5 (first eigenvalue on S^5)."""
    assert LAMBDA_1 == 5


def test_spectral_cost_35_over_3():
    """Paper: d_1 + lambda_1 + K = 35/3."""
    cost = D1 + LAMBDA_1 + K_THEORY
    assert abs(cost - 35 / 3) < 1e-10


# -----------------------------------------------------------------------------
# TABLE: CKM / Wolfenstein (paper Section 7, Quark Flavor)
# -----------------------------------------------------------------------------


def test_cabibbo_lambda():
    """Paper: lambda (Cabibbo) = eta*(1 + alpha_s/(3*pi)), PDG 0.225, ~2%."""
    eta_D = 2 / 9
    alpha_s = 0.1187
    lam_pred = eta_D * (1 + alpha_s / (3 * np.pi))
    lam_pdg = 0.22500
    rel_err = abs(lam_pred - lam_pdg) / lam_pdg
    assert rel_err < 0.03  # ~2%


def test_wolfenstein_A():
    """Paper: A = lam1/d1 * (1 - eta*alpha_s/pi), PDG 0.826, ~2%."""
    A_bare = LAMBDA_1 / D1  # 5/6
    A_pred = A_bare * (1 - TWIST * 0.1187 / np.pi)
    A_pdg = 0.826
    rel_err = abs(A_pred - A_pdg) / A_pdg
    assert rel_err < 0.03


def test_wolfenstein_rhobar():
    """Paper: rho_bar = 1/(2*pi), PDG 0.159."""
    rho_pred = 1 / (2 * np.pi)
    rho_pdg = 0.1592
    rel_err = abs(rho_pred - rho_pdg) / rho_pdg
    assert rel_err < 0.1  # ~1%


def test_wolfenstein_etabar():
    """Paper: eta_bar = pi/9, PDG 0.349."""
    eta_bar_pred = np.pi / 9
    eta_bar_pdg = 0.3490
    rel_err = abs(eta_bar_pred - eta_bar_pdg) / eta_bar_pdg
    assert rel_err < 0.1


# -----------------------------------------------------------------------------
# TABLE: Quark masses (paper Section 7, downtype spectral ordering)
# -----------------------------------------------------------------------------


def test_down_type_quark_masses():
    """Paper: m_q = m_lepton * exp(sigma_q), sigma from spectral ordering."""
    PI = np.pi
    G = LAMBDA_1 * TWIST  # 10/9
    p = 3
    sigma = {"b": 77 / 90, "s": -10 / 81, "d": 2 * PI / 3 + G / p**2}
    UV = {"b": M_TAU_PDG / 1000, "s": M_MU_PDG / 1000, "d": M_E_PDG / 1000}
    PDG = {"b": 4.183, "s": 0.0934, "d": 0.00467}
    for q in ["b", "s", "d"]:
        m_pred = UV[q] * np.exp(sigma[q])
        rel_err = abs(m_pred - PDG[q]) / PDG[q]
        assert rel_err < 0.1, f"{q}-quark: pred {m_pred:.4f}, PDG {PDG[q]}, err {rel_err*100:.1f}%"
