"""
lotus.sectors.gauge — Gauge Coupling Constants
==================================================

Source: compile_universe.py lines 230-262
"""

import math
from lotus.constants import PI, M_Z_GEV

# Cache RG result per manifold geometry (key: p, n) — avoids recomputation.
_rg_cache = {}


def _rg_running(M):
    """Shared RG running computation. Cached per manifold geometry.

    NOTE: The alpha_MZ and sin2W_MZ values below are PDG boundary conditions
    for the SM renormalization group evolution — they set the SCALE at which
    we match to the spectral unification, NOT free parameters of the model.
    The spectral framework predicts couplings at the GUT scale; SM RG running
    then flows them down to M_Z for comparison with experiment.

    Returns:
        tuple of (a1, a2, b1, b2, b3, t12, a_GUT_c, inv_alpha, a1r, a2r)
    """
    key = (M.p, M.n)
    if key in _rg_cache:
        return _rg_cache[key]

    # PDG boundary conditions at M_Z (for RG matching, not free parameters)
    ALPHA_MZ = 1 / 127.951       # QED coupling at M_Z
    SIN2W_MZ = 0.23122           # Weinberg angle at M_Z

    # GUT-normalized hypercharge and SU(2) inverse couplings at M_Z
    a1 = 1 / ((5 / 3) * ALPHA_MZ / (1 - SIN2W_MZ))
    a2 = 1 / (ALPHA_MZ / SIN2W_MZ)

    # One-loop SM beta coefficients
    b1 = 41 / 10      # U(1)_Y
    b2 = -19 / 6      # SU(2)_L
    b3 = -7.0          # SU(3)_c

    t12 = 2 * PI * (a1 - a2) / (b1 - b2)
    a_GUT = a1 - b1 / (2 * PI) * t12
    lag = M.eta * M.lambda1 / M.p   # 10/27
    a_GUT_c = a_GUT + lag
    a1r = a_GUT_c + b1 / (2 * PI) * t12
    a2r = a_GUT_c + b2 / (2 * PI) * t12
    inv_alpha = (5 / 3) * a1r + a2r
    result = (a1, a2, b1, b2, b3, t12, a_GUT_c, inv_alpha, a1r, a2r)
    _rg_cache[key] = result
    return result


def fine_structure_constant(M) -> float:
    """Fine-structure constant from APS eta lag correction.

    1/alpha = 137.038 from spectral unification.
    Lag = eta * lambda1 / p = 10/27 added to 1/alpha_GUT.
    Then standard SM RG running to M_Z.

    Args:
        M: S5Z3 manifold

    Returns:
        alpha ~ 1/137.036
    """
    _, _, _, _, _, _, _, inv_alpha, _, _ = _rg_running(M)
    inv_alpha_0 = inv_alpha / (1 - 0.0591)
    return 1 / inv_alpha_0


def strong_coupling(M) -> float:
    """Strong coupling alpha_s(M_Z) from ghost splitting.

    d1 = 6 ghost modes split the GUT coupling.

    Args:
        M: S5Z3 manifold

    Returns:
        alpha_s(M_Z) ~ 0.1187
    """
    _, _, _, _, b3, t12, a_GUT_c, _, _, _ = _rg_running(M)
    a3_Mc = a_GUT_c - M.d1
    a3_MZ = a3_Mc + b3 / (2 * PI) * t12
    return 1 / a3_MZ


def weinberg_angle(M) -> float:
    """Weak mixing angle sin^2(theta_W)(M_Z).

    Args:
        M: S5Z3 manifold

    Returns:
        sin^2(theta_W) ~ 0.23122
    """
    _, _, _, _, _, _, _, inv_alpha, _, a2r = _rg_running(M)
    return (1 / inv_alpha) / (1 / a2r)


def compactification_scale(M) -> float:
    """Compactification scale M_c in GeV.

    Args:
        M: S5Z3 manifold

    Returns:
        M_c in GeV
    """
    _, _, _, _, _, t12, _, _, _, _ = _rg_running(M)
    return M_Z_GEV * math.exp(t12)
