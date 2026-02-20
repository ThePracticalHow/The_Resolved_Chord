"""
lotus.sectors.neutrinos — Neutrino Masses from Fold-Wall Tunneling
=====================================================================

Source: compile_universe.py lines 223-227
"""

from lotus.constants import M_E_GEV


def neutrino_mass(M, m_p: float) -> float:
    """Heaviest neutrino mass m_nu3 from the Geometric Seesaw.

    m_nu3 = m_e^3 / (p * m_p^2)

    Fold-wall tunneling: projection 1/p, round-trip penetration (m_e/m_p)^2.
    No fictitious right-handed neutrino required.

    Args:
        M: S5Z3 manifold
        m_p: proton mass in GeV

    Returns:
        m_nu3 in GeV
    """
    return M_E_GEV ** 3 / (M.p * m_p ** 2)


def neutrino_mass_meV(M, m_p: float) -> float:
    """Heaviest neutrino mass in meV.

    Args:
        M: S5Z3 manifold
        m_p: proton mass in GeV

    Returns:
        m_nu3 in meV
    """
    return neutrino_mass(M, m_p) * 1e12  # GeV -> meV


def mass_squared_ratio(M) -> float:
    """Neutrino mass-squared splitting ratio.

    Delta m^2_atm / Delta m^2_sol = d1^2 - p = 36 - 3 = 33

    Args:
        M: S5Z3 manifold

    Returns:
        33
    """
    return M.d1 ** 2 - M.p
