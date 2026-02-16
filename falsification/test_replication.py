"""
Independent Replication Tests
=============================
Verify leng_replication.py agrees with leng_test_utils.
Blocks: "Coding artifact / bug / hidden assumption."
"""

import pytest

# Import replication (standalone, no leng_test_utils)
import leng_replication as rep

from leng_test_utils import (
    predict_masses_from_e,
    is_self_consistent,
    all_positive_universe,
    resonance_gap,
    twist_formula,
    koide_ratio,
    eta_twisted_donnelly,
    reidemeister_torsion,
    d1_degeneracy,
    M_E_PDG,
    M_MU_PDG,
    M_TAU_PDG,
)


class TestSurvivorIdentity:
    """Both codebases agree on unique survivor."""

    def test_replication_survivors_match(self):
        surv = rep.run_scan(8, 13)
        assert surv == [(3, 3)]

    def test_replication_vs_utils_scan(self):
        survivors = []
        for n in range(1, 9):
            for p in range(2, 14):
                if is_self_consistent(n, p) and all_positive_universe(n, p):
                    if rep._is_prime(p) and abs(rep._K_p(p) - 1.0) > 1e-9:
                        survivors.append((n, p))
        assert survivors == rep.run_scan(8, 13)


class TestMassPredictions:
    """Both codebases agree on predicted masses."""

    def test_masses_match(self):
        u_e, u_mu, u_tau = predict_masses_from_e(M_E_PDG)
        r_e, r_mu, r_tau = rep.predict_masses(M_E_PDG)
        assert u_e == pytest.approx(r_e, rel=1e-10)
        assert u_mu == pytest.approx(r_mu, rel=1e-10)
        assert u_tau == pytest.approx(r_tau, rel=1e-10)

    def test_muon_agreement_both(self):
        _, u_mu, _ = predict_masses_from_e(M_E_PDG)
        _, r_mu, _ = rep.predict_masses(M_E_PDG)
        assert abs(u_mu - M_MU_PDG) / M_MU_PDG < 0.001
        assert abs(r_mu - M_MU_PDG) / M_MU_PDG < 0.001


class TestEtaIdentity:
    """Both codebases agree on eta sum = d₁·τ_R."""

    def test_eta_sum_match(self):
        rep_sum = rep.eta_sum_twist(3, 3)
        utils_sum = sum(abs(eta_twisted_donnelly(3, 3, m)) for m in range(1, 3))
        assert rep_sum == pytest.approx(utils_sum, rel=1e-10)
        assert rep_sum == pytest.approx(2/9, rel=1e-10)

    def test_d1_tau_match(self):
        rep_dt = rep.d1_tau(3, 3)
        utils_dt = d1_degeneracy(3) * reidemeister_torsion(3, 3)
        assert rep_dt == pytest.approx(utils_dt, rel=1e-10)


class TestTwistFormula:
    """Twist and resonance gap match."""

    def test_twist_match(self):
        for n, p in [(3, 3), (4, 2), (2, 5)]:
            assert rep._twist(n, p) == pytest.approx(twist_formula(n, p), rel=1e-12)

    def test_resonance_gap_match(self):
        for n, p in [(3, 3), (4, 2), (2, 2)]:
            assert rep._resonance_gap(n, p) == pytest.approx(resonance_gap(n, p), rel=1e-12)
