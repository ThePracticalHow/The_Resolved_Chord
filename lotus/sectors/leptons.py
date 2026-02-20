"""
lotus.sectors.leptons — Lepton Masses from Koide + Eta
=========================================================

Source: compile_universe.py lines 188-198, EtaInvariant.py
"""

import math
from lotus.constants import PI, M_E_GEV


def lepton_masses(M) -> dict:
    """Koide formula with eta phase correction.

    The Koide ratio K = 2/3 fixes the mass spectrum via the
    phase delta = 2*pi/3 + 2/9 (the eta invariant IS the
    Koide phase correction).

    Args:
        M: S5Z3 manifold

    Returns:
        dict with 'electron', 'muon', 'tau' (GeV) and
        'electron_ratio', 'muon_ratio', 'tau_ratio' (m/m_e)
    """
    delta_phase = 2 * PI / 3 + M.eta
    r_koide = math.sqrt(2)
    sqm = [1 + r_koide * math.cos(delta_phase + 2 * PI * k / 3)
           for k in range(3)]
    mu = math.sqrt(M_E_GEV) / sqm[0]
    m_e_p, m_mu_p, m_tau_p = [(mu * s) ** 2 for s in sqm]
    return {
        'electron': m_e_p,
        'muon': m_mu_p,
        'tau': m_tau_p,
        'electron_ratio': m_e_p / M_E_GEV,
        'muon_ratio': m_mu_p / M_E_GEV,
        'tau_ratio': m_tau_p / M_E_GEV,
    }
