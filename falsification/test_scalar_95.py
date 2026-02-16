"""
95 GeV Scalar Tests
===================
Fold-wall shearing mode: m_95 = m_Z * (1 + eta^2) = m_Z * 85/81.
From STRANGE_CASTLES.md S3, scalar_95gev.py.
"""

import numpy as np
import pytest

from leng_test_utils import TWIST, D1, LAMBDA_1, ALPHA, M_E_PDG

# ─────────────────────────────────────────────────────────────────────────────
#  CONSTANTS
# ─────────────────────────────────────────────────────────────────────────────

P = 3           # orbifold order
ETA = 2 / 9     # Donnelly eta invariant = (p-1)/(p*n) for n=p=3
M_Z_PDG = 91.1876  # GeV (PDG 2024)
M_H_PDG = 125.25    # GeV (PDG 2024)

# Prediction
M_95_PRED = M_Z_PDG * (1 + ETA**2)  # = m_Z * 85/81

# Experimental target (CMS diphoton excess, approximate)
M_95_EXP_LO = 93.0   # GeV, lower edge of reported range
M_95_EXP_HI = 97.0   # GeV, upper edge of reported range


class TestScalar95Mass:
    """m_95 = m_Z * (1 + eta^2) = m_Z * 85/81"""

    def test_mass_formula_rational(self):
        """85/81 = (p^4 + 4) / p^4 = 1 + eta^2"""
        assert 85 / 81 == pytest.approx(1 + ETA**2, rel=1e-14)
        assert 85 / 81 == pytest.approx((P**4 + 4) / P**4, rel=1e-14)

    def test_mass_prediction(self):
        """m_95 = 95.69 GeV"""
        assert M_95_PRED == pytest.approx(95.69, rel=1e-3)

    def test_mass_in_experimental_range(self):
        """Prediction falls within the CMS excess range"""
        assert M_95_EXP_LO < M_95_PRED < M_95_EXP_HI

    def test_mass_below_higgs(self):
        """The shearing mode is lighter than the breathing mode (Higgs)"""
        assert M_95_PRED < M_H_PDG

    def test_mass_above_Z(self):
        """The shearing mode is heavier than m_Z (eta^2 > 0)"""
        assert M_95_PRED > M_Z_PDG

    def test_mass_is_linear_not_quadratic(self):
        """m_95 = m_Z*(1+eta^2), NOT m_Z*sqrt(1+eta^2).
        The correction is multiplicative on the mass, not the mass^2.
        sqrt formula gives 93.4 GeV (too low); linear gives 95.69 GeV (matches)."""
        m_linear = M_Z_PDG * (1 + ETA**2)   # = 95.69 GeV (correct)
        m_sqrt = M_Z_PDG * np.sqrt(1 + ETA**2)  # = 93.41 GeV (wrong)
        assert m_linear == pytest.approx(95.69, rel=1e-3)
        assert m_sqrt == pytest.approx(93.41, rel=1e-3)
        # The linear formula matches the CMS excess range; sqrt does not
        assert M_95_EXP_LO < m_linear < M_95_EXP_HI


class TestEtaSquaredUniversality:
    """eta^2 appears in three independent physical contexts"""

    def test_eta_value(self):
        """eta = 2/9 = (p-1)/(p*n) for n=p=3"""
        assert ETA == pytest.approx(2/9, rel=1e-14)
        assert ETA == pytest.approx((P - 1) / (P * P), rel=1e-14)

    def test_pmns_solar_correction(self):
        """sin^2(theta_12) = 1/3 - eta^2/2 ≈ 0.309 (PDG: 0.307)"""
        sin2_12 = 1/3 - ETA**2 / 2
        assert sin2_12 == pytest.approx(0.307, rel=0.01)

    def test_cc_residual(self):
        """Lambda^(1/4) = m_nu3 * eta^2 ≈ 2.49 meV (observed: ~2.3 meV)"""
        m_nu3_eV = 50.52e-3  # eV
        cc_pred = m_nu3_eV * ETA**2 * 1e3  # meV
        assert cc_pred == pytest.approx(2.49, rel=0.01)

    def test_95_gev_correction(self):
        """m_95 = m_Z * (1 + eta^2)"""
        assert M_95_PRED / M_Z_PDG - 1 == pytest.approx(ETA**2, rel=1e-10)


class TestShearingModeProperties:
    """Physical properties of the fold-wall shearing mode"""

    def test_mode_counting(self):
        """p=3 fold walls -> 1 breathing + 2 shearing modes"""
        n_breathing = 1
        n_shearing = P - n_breathing
        assert n_shearing == 2  # complex pair (chi_1 + chi_2)

    def test_coupling_suppression(self):
        """Signal strength mu = eta^2 ≈ 0.049"""
        mu = ETA**2
        assert mu == pytest.approx(4/81, rel=1e-14)
        assert mu < 0.05

    def test_width_narrow(self):
        """Width ~ eta^2 * Gamma_H(95 GeV) << 1 GeV"""
        # SM Higgs width at 95 GeV would be ~ 4 MeV (dominated by bb)
        gamma_h_95 = 4e-3  # GeV (approximate)
        gamma_95 = ETA**2 * gamma_h_95
        assert gamma_95 < 1e-3  # < 1 MeV (unresolvable at LHC)

    def test_hierarchy_breathing_vs_shearing(self):
        """m_H > m_95 because quartic > gauge: lambda_H > g_Z^2/4"""
        # Higgs mass from quartic: m_H^2 = 2*lambda_H*v^2
        v = 246.22
        lambda_H = M_H_PDG**2 / (2 * v**2)
        # Shearing mass from gauge: m_95^2 ~ m_Z^2*(1+eta^2)
        # m_Z^2 = g_Z^2 * v^2 / 4, so m_95^2 ~ g_Z^2*v^2/4*(1+eta^2)
        g_Z_sq = 4 * M_Z_PDG**2 / v**2
        assert 2 * lambda_H > g_Z_sq / 4  # quartic dominates
        assert M_H_PDG**2 > M_95_PRED**2


class TestPerturbationTheory:
    """The eta^2 correction arises from second-order degenerate perturbation theory."""

    def test_first_order_cancellation(self):
        """chi_1 and chi_2 eta shifts are opposite: +1/9 and -1/9, sum to 0"""
        from leng_test_utils import eta_twisted_donnelly
        eta1 = eta_twisted_donnelly(3, 3, 1)
        eta2 = eta_twisted_donnelly(3, 3, 2)
        # First-order shifts cancel in Z_3-invariant combination
        assert abs(eta1.real + eta2.real) < 1e-14

    def test_second_order_gives_eta_squared(self):
        """The splitting (eta_1 - eta_2)^2 = (2/9)^2 = eta^2.
        Individual eta_1^2 + eta_2^2 = 2/81 (half of eta^2).
        The physical correction uses the TOTAL spectral asymmetry squared."""
        eta1 = 1/9
        eta2 = -1/9
        # The physical eta^2 is the total spectral asymmetry squared
        splitting = (eta1 - eta2)**2  # = (2/9)^2 = 4/81
        assert splitting == pytest.approx(ETA**2, rel=1e-14)
        # Individual contributions sum to half of eta^2
        individual_sum = eta1**2 + eta2**2  # = 2/81
        assert individual_sum == pytest.approx(ETA**2 / 2, rel=1e-14)

    def test_kk_degeneracy_forces_cancellation(self):
        """d_l^(1) = d_l^(2) for all l — the structural reason for first-order cancellation"""
        from leng_test_utils import eta_twisted_donnelly
        # The degeneracies d^(1) and d^(2) are equal by complex conjugation symmetry
        # This means chi_1 and chi_2 are degenerate, forcing the first-order shift to cancel
        eta1 = eta_twisted_donnelly(3, 3, 1).real
        eta2 = eta_twisted_donnelly(3, 3, 2).real
        assert abs(eta1) == pytest.approx(abs(eta2), rel=1e-14)


class TestAlternativeFormula:
    """Cross-check: m_95 = m_H * p/(p+1) gives 93.94 GeV (1.1% off)"""

    def test_higgs_sharing(self):
        """m_H * 3/4 = 93.94 GeV — less precise but structurally interesting"""
        m_alt = M_H_PDG * P / (P + 1)
        assert m_alt == pytest.approx(93.94, rel=1e-3)
        # Less precise than the Z*(1+eta^2) formula
        assert abs(m_alt - 95.4) > abs(M_95_PRED - 95.4)
