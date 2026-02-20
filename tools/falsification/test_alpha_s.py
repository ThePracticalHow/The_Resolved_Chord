"""
Tests for alpha_s constraint derivation (Dirichlet at cone point of B^6/Z_3).
"""

import numpy as np
import pytest

from alpha_s_constraint import (
    ALPHA_S_GAP,
    predict_alpha_s_1loop,
    predict_alpha_s_2loop,
    proton_cross_check,
)


def test_gap_theory():
    """Geometric constraint: Delta(1/alpha_3) = pi^2 - 5."""
    assert abs(ALPHA_S_GAP - (np.pi**2 - 5)) < 1e-10
    assert 4.86 < ALPHA_S_GAP < 4.88


def test_alpha_s_1loop():
    """1-loop prediction returns sensible values."""
    res = predict_alpha_s_1loop()
    assert "alpha_s_pred" in res
    assert "M_c" in res
    assert "gap" in res
    assert 0.10 < res["alpha_s_pred"] < 0.14
    assert 1e12 < res["M_c"] < 1e19  # GUT scale ~10^13–10^16 GeV


def test_alpha_s_2loop():
    """2-loop prediction matches PDG to within ~1 sigma."""
    res = predict_alpha_s_2loop()
    alpha_s_pred = res["alpha_s_pred"]
    # PDG: 0.1180 +/- 0.0009
    assert 0.115 < alpha_s_pred < 0.121
    # Constraint prediction is ~0.1174
    assert abs(alpha_s_pred - 0.1174) < 0.01


def test_proton_cross_check():
    """Proton cross-check returns dict with predicted and actual m_p/m_e."""
    res_2L = predict_alpha_s_2loop()
    proton = proton_cross_check(res_2L["alpha_s_pred"])
    assert "m_p_over_m_e_pred" in proton
    assert "m_p_over_m_e_actual" in proton
    assert "rel_err" in proton
    # Should be within a few percent
    assert proton["rel_err"] < 0.05
