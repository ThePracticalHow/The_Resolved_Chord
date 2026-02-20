"""
lotus.sectors.mixing — CKM and PMNS Matrices
================================================

Source: compile_universe.py lines 284-304

Hurricane corrections (1-loop):
  λ_CKM = η * (1 + α_s/(p·π))     [QCD dressing]
  A_CKM = (λ₁/d₁) * (1 − η·α_s/π) [QCD vertex]
  Tree:   λ = η,  A = λ₁/d₁  (pure spectral invariants)
"""

import math
from lotus.constants import PI


def ckm_matrix(M, alpha_s: float, loops: int = 2) -> dict:
    """CKM parameters from spectral invariants + hurricane corrections.

    Tree:
      lambda  = eta
      A       = lambda1 / d1

    1-loop (hurricane):
      lambda  = eta * (1 + alpha_s/(p*pi))
      A       = (lambda1/d1) * (1 - eta*alpha_s/pi)

    Always:
      rho_bar = 1/(2*pi)
      eta_bar = pi/9
      gamma   = arctan(eta_bar/rho_bar) = arctan(2*pi^2/9)

    Args:
        M: S5Z3 manifold
        alpha_s: strong coupling constant
        loops: 0 (tree), 1, or 2

    Returns:
        dict with 'lambda', 'A', 'rho_bar', 'eta_bar', 'gamma_deg', 'J'
    """
    if loops == 0:
        lam = M.eta
        A = M.lambda1 / M.d1
    else:
        lam = M.eta * (1 + alpha_s / (M.p * PI))
        A = (M.lambda1 / M.d1) * (1 - M.eta * alpha_s / PI)

    rho_bar = 1 / (2 * PI)
    eta_bar = PI / 9
    gamma = math.atan(eta_bar / rho_bar)
    J = A ** 2 * lam ** 6 * eta_bar

    return {
        'lambda': lam,
        'A': A,
        'rho_bar': rho_bar,
        'eta_bar': eta_bar,
        'gamma_deg': math.degrees(gamma),
        'J': J,
    }


def pmns_matrix(M) -> dict:
    """PMNS neutrino mixing angles from spectral impedance.

    sin^2(theta_23) = d1/(d1+lambda1) = 6/11   (ghost fraction)
    sin^2(theta_12) = p/pi^2 = 3/pi^2           (impedance mismatch)
    sin^2(theta_13) = (eta*K)^2 = (4/27)^2      (double crossing)

    Args:
        M: S5Z3 manifold

    Returns:
        dict with 'sin2_theta23', 'sin2_theta12', 'sin2_theta13'
    """
    return {
        'sin2_theta23': M.d1 / (M.d1 + M.lambda1),
        'sin2_theta12': M.p / PI ** 2,
        'sin2_theta13': (M.eta * M.K) ** 2,
    }


def cabibbo_angle(M, alpha_s: float, loops: int = 2) -> float:
    """Cabibbo mixing angle (sin theta_c).

    Tree:   sin theta_c = eta
    1-loop: sin theta_c = eta * (1 + alpha_s/(3*pi))

    Args:
        M: S5Z3 manifold
        alpha_s: strong coupling constant
        loops: 0 (tree), 1, or 2

    Returns:
        sin theta_c ~ 0.2250
    """
    if loops == 0:
        return M.eta
    return M.eta * (1 + alpha_s / (M.p * PI))

