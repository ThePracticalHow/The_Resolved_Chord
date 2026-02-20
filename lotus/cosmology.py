"""
lotus.cosmology — Cosmological Predictions
=============================================

Source: compile_universe.py lines 328-358

Hurricane corrections:
  CC:   1 + η²/π   (inside-outside inversion, 1-loop)
  η_B:  α⁴·η is already the base; no additional loop correction identified
"""

import math
import warnings
from lotus.constants import PI, M_E_GEV


def cosmological_constant(M, m_nu3: float, loops: int = 2) -> float:
    """Cosmological constant Lambda^{1/4}.

    Lambda^{1/4} = m_nu3 * eta^2 * (1 - K/d1) = m_nu3 * 32/729

    Heavy modes cancel by Z_3 equidistribution.  Only the lightest
    tunneling mode m_nu3 survives, suppressed by eta^2.

    Hurricane (1-loop): inside-outside inversion 1 + η²/π.

    Args:
        M: manifold
        m_nu3: heaviest neutrino mass in GeV
        loops: 0 (tree), 1, or 2

    Returns:
        Lambda^{1/4} in GeV
    """
    bare = m_nu3 * M.eta ** 2 * (1 - M.K / M.d1)
    if loops == 0:
        return bare
    # 1-loop: CC hurricane (inside-outside inversion)
    hurricane = 1 + M.eta ** 2 / PI
    return bare * hurricane


def cosmological_constant_meV(M, m_nu3: float, loops: int = 2) -> float:
    """Cosmological constant Lambda^{1/4} in meV.

    Args:
        M: manifold
        m_nu3: heaviest neutrino mass in GeV
        loops: 0 (tree), 1, or 2

    Returns:
        Lambda^{1/4} in meV
    """
    return cosmological_constant(M, m_nu3, loops=loops) * 1e12


def inflation_efolds(M) -> float:
    """Number of e-folds: N = (d1+lambda1)^2 * lambda1^2 / (p * 4^2) = 3025/48.

    Args:
        M: manifold

    Returns:
        N ~ 63
    """
    return (M.d1 + M.lambda1) ** 2 * M.lambda1 ** 2 / (M.p * 16)


def spectral_index(N: float) -> float:
    """Spectral tilt: n_s = 1 - 2/N (Starobinsky).

    Returns:
        n_s ~ 0.968
    """
    return 1 - 2 / N


def tensor_to_scalar(N: float) -> float:
    """Tensor-to-scalar ratio: r = 12/N^2 (Starobinsky).

    Returns:
        r ~ 0.003
    """
    return 12 / N ** 2


def baryogenesis(M, alpha: float, loops: int = 2) -> float:
    """Why there is something instead of nothing.

    eta_B = alpha^4 * eta = alpha^4 * (2/9)

    The Big Bang should have created exactly equal matter and antimatter.
    They should have annihilated completely, leaving only photons.
    But for every 10 billion antimatter particles, there were 10 billion
    AND ONE matter particles. That tiny excess (eta_B ~ 6.1e-10)
    is why galaxies, stars, and humans exist.

    The Donnelly eta invariant (eta = 2/9) IS the topological bias.
    A perfectly symmetric universe is forbidden by Z_3:
    the orbifold fold-walls are strictly asymmetric by construction.

    The exponent 4 = dim(fold wall): four EM vertices at the spectral
    phase transition. The coupling = alpha_em. The bias = eta.

    Predicted: 6.30e-10.  Measured: 6.10e-10.  Error: 3.3%.

    Args:
        M: S5Z3 manifold
        alpha: fine-structure constant
        loops: 0 (tree), 1, or 2 — currently no loop correction identified

    Returns:
        eta_B ~ 6.1e-10
    """
    return alpha ** 4 * M.eta


def dark_matter_ratio(M) -> float:
    """Dark matter to baryon ratio: Omega_DM/Omega_B = d1 - K = 16/3.

    Ghost mode counting minus Koide harmonic lock.

    Args:
        M: manifold

    Returns:
        16/3 = 5.333...
    """
    return M.d1 - M.K


def cosmic_snapshot(M) -> tuple:
    """Cosmic energy budget from spectral snapshot.

    Omega_Lambda/Omega_m = 2*pi^2/9
    Resolves the coincidence problem — only p=3 gives observed range.

    Args:
        M: manifold

    Returns:
        (ratio, Omega_Lambda, Omega_matter)
    """
    ratio = 2 * PI ** 2 / M.p ** 2
    OL = ratio / (1 + ratio)
    Om = 1 / (1 + ratio)
    return ratio, OL, Om


def hubble_constant(M, CC_meV: float, M_Planck: float,
                    omega_lambda: float) -> float:
    """Hubble constant from spectral Friedmann equation.

    H₀ = √(8π/3 · Λ / (M_P² · Ω_Λ))

    All inputs are spectral predictions:
      Λ^{1/4}  from CC (spectral tunneling)
      M_P      from 5-lock gravity theorem
      Ω_Λ      from cosmic snapshot

    Predicted: ≈67.7 km/s/Mpc.  Measured (Planck): 67.4 km/s/Mpc.

    Args:
        M: manifold
        CC_meV: Λ^{1/4} in meV
        M_Planck: Planck mass in GeV
        omega_lambda: Ω_Λ

    Returns:
        H₀ in km/s/Mpc
    """
    CC_eV = CC_meV * 1e-3
    CC_GeV = CC_eV * 1e-9
    H_0_natural = math.sqrt(8 * PI / 3 * CC_GeV ** 4 /
                            (M_Planck ** 2 * omega_lambda))
    # Convert H₀ from natural units (GeV) to km/s/Mpc:
    #   H₀ [s⁻¹] = H₀ [GeV] / ℏ [GeV·s]
    #   H₀ [km/s/Mpc] = H₀ [s⁻¹] / (3.2408e-20 s⁻¹ per km/s/Mpc)
    hbar_GeV_s = 6.582119569e-25
    km_s_Mpc_to_per_s = 3.2408e-20  # 1 km/s/Mpc in s⁻¹
    H_0_per_s = H_0_natural / hbar_GeV_s if H_0_natural > 0 else 0.0
    H_0_km_s_Mpc = H_0_per_s / km_s_Mpc_to_per_s
    return H_0_km_s_Mpc


def universe_age(H_0: float, omega_m: float = 0.3111,
                 omega_lambda: float = 0.6889) -> float:
    """Age of the universe from ΛCDM Friedmann integral.

    t_age = (1/H₀) × ∫₀^∞ dz / [(1+z) √(Ω_m(1+z)³ + Ω_Λ)]

    Uses spectral Ω_m and Ω_Λ (from cosmic snapshot) as defaults.

    Predicted: ≈13.7 Gyr.  Measured (Planck): 13.80 Gyr.

    Args:
        H_0: Hubble constant in km/s/Mpc
        omega_m: matter density parameter
        omega_lambda: dark energy density parameter

    Returns:
        t_age in Gyr
    """
    # Numerical integration via Simpson's rule (no scipy dependency)
    N = 10000
    z_max = 1000.0  # upper limit (effectively infinity)
    dz = z_max / N
    integral = 0.0
    for i in range(N):
        z = (i + 0.5) * dz
        integrand = 1.0 / ((1 + z) * math.sqrt(
            omega_m * (1 + z) ** 3 + omega_lambda))
        integral += integrand * dz

    # H₀ in km/s/Mpc → s⁻¹
    km_s_Mpc_to_per_s = 3.2408e-20
    H_0_per_s = H_0 * km_s_Mpc_to_per_s
    t_s = integral / H_0_per_s
    t_Gyr = t_s / (365.25 * 24 * 3600 * 1e9)
    return t_Gyr


def reheating_temperature(M, M_Planck: float, N: float,
                          m_H: float = None, v: float = None) -> float:
    """Reheating temperature after inflation.

    T_reheat = (90/(π² g*))^{1/4} · √(Γ_φ · M_P)

    In the LOTUS framework the inflaton IS the Higgs field (fold field phi).
    The Higgs-channel decay rate dominates:

      Γ_φ = m_H³ / (8π v²)

    which gives T_reheat ≈ 2.15×10⁹ GeV (all inputs spectral).

    If m_H and v are not supplied, falls back to the gravitational
    (Starobinsky) decay channel: Γ_φ = M_inflaton³ / (192π M_P²).

    Args:
        M: manifold
        M_Planck: Planck mass in GeV
        N: number of e-folds
        m_H: Higgs mass in GeV (optional; use for Higgs-channel decay)
        v: Higgs VEV in GeV (optional; use for Higgs-channel decay)

    Returns:
        T_reheat in GeV
    """
    g_star = 106.75
    if m_H is not None and v is not None:
        # Higgs-channel decay (LOTUS: inflaton = fold field = Higgs)
        # Gamma = m_H^3 / (8*pi*v^2)  [perturbative coupling to gauge bosons]
        Gamma_phi = m_H ** 3 / (8 * PI * v ** 2)
    else:
        # Fallback: gravitational (Starobinsky) decay channel
        M_inflaton = M_Planck / math.sqrt(6 * N)
        Gamma_phi = M_inflaton ** 3 / (192 * PI * M_Planck ** 2)
    T_reheat = (90 / (PI ** 2 * g_star)) ** 0.25 * \
               math.sqrt(Gamma_phi * M_Planck)
    return T_reheat


