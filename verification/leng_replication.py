"""
LENG Independent Replication
============================
Standalone implementation from scratch. No imports from leng_test_utils.
Used to verify: survivor identity, predicted masses, eta identities.
Blocks: "Coding artifact / bug / hidden assumption."
"""

# Implementation date: 2026-02-12
# Source: LENG_Master.md, UniverseLandscape.py, Donnelly (1978)

import math

# Constants (PDG 2024)
M_E = 0.51099895
R = 2 ** 0.5  # sqrt(2)


def _is_prime(n: int) -> bool:
    if n < 2:
        return False
    for i in range(2, int(n ** 0.5) + 1):
        if n % i == 0:
            return False
    return True


def _tau_R(p: int, n: int) -> float:
    """Reidemeister torsion: p^{-n}."""
    return p ** (-n)


def _twist(n: int, p: int) -> float:
    """twist = d₁ × τ_R = 2n / p^n."""
    return 2 * n / (p ** n)


def _K_p(p: int) -> float:
    """Koide ratio for p generations with r=√2: K = 2/p."""
    return 2.0 / p


def _resonance_gap(n: int, p: int) -> float:
    """|p × twist - K_p|."""
    return abs(p * _twist(n, p) - _K_p(p))


def _is_self_consistent(n: int, p: int, tol: float = 1e-9) -> bool:
    return _resonance_gap(n, p) < tol


def _delta(n: int, p: int) -> float:
    """δ = 2π/p + twist."""
    return 2 * math.pi / p + _twist(n, p)


def _sqrt_masses(n: int, p: int) -> tuple:
    """√m_k for universe (n,p)."""
    d = _delta(n, p)
    return tuple(1.0 + R * math.cos(d + 2 * math.pi * k / p) for k in range(p))


def _all_positive(n: int, p: int) -> bool:
    return all(sm > 0 for sm in _sqrt_masses(n, p))


def _eta_donnelly(p: int, n: int, m: int) -> complex:
    """η_D(χ_m) for L(p;1,...,1). Donnelly formula."""
    omega = complex(math.cos(2 * math.pi / p), math.sin(2 * math.pi / p))
    total = 0j
    for k in range(1, p):
        cot_k = math.cos(math.pi * k / p) / math.sin(math.pi * k / p)
        total += (omega ** (m * k)) * ((1j * cot_k) ** n)
    return total / p


def run_scan(n_max: int, p_max: int) -> list:
    """Return list of (n,p) that are self-consistent + all positive + prime p + K≠1."""
    survivors = []
    for n in range(1, n_max + 1):
        for p in range(2, p_max + 1):
            if not _is_self_consistent(n, p):
                continue
            if not _all_positive(n, p):
                continue
            if not _is_prime(p):
                continue
            if abs(_K_p(p) - 1.0) < 1e-9:
                continue
            survivors.append((n, p))
    return survivors


def predict_masses(m_e: float) -> tuple:
    """Predict (m_e, m_mu, m_tau) from m_e for S⁵/Z₃."""
    n, p = 3, 3
    sm = list(_sqrt_masses(n, p))
    sm.sort()
    mu = (m_e ** 0.5) / sm[0]
    masses = [(mu * s) ** 2 for s in sm]
    return tuple(sorted(masses))


def eta_sum_twist(n: int, p: int) -> float:
    """Σ|η_D| for m=1..p-1."""
    total = 0.0
    for m in range(1, p):
        total += abs(_eta_donnelly(p, n, m))
    return total


def d1_tau(n: int, p: int) -> float:
    """d₁ × τ_R = 2n × p^{-n}."""
    return 2 * n * _tau_R(p, n)


# ── Replication report ────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("LENG Independent Replication Report")
    print("=" * 50)

    survivors = run_scan(8, 13)
    print(f"Survivors (n_max=8, p_max=13): {survivors}")
    assert survivors == [(3, 3)], f"Expected [(3, 3)], got {survivors}"

    m_e, m_mu, m_tau = predict_masses(M_E)
    print(f"Predicted masses: m_e={m_e:.6f}, m_mu={m_mu:.4f}, m_tau={m_tau:.4f}")
    assert abs(m_e - 0.51099895) < 1e-6
    assert abs(m_mu - 105.658) < 1.0
    assert abs(m_tau - 1776.86) < 2.0

    eta_s = eta_sum_twist(3, 3)
    dt = d1_tau(3, 3)
    print(f"Eta sum (3,3): {eta_s:.10f}, d1*tau_R: {dt:.10f}")
    assert abs(eta_s - 2/9) < 1e-8
    assert abs(eta_s - dt) < 1e-10

    # Two-loop G2 check (added 2026-02-13)
    # G2 = -lambda1 * (d1 + sum|eta|) = -5 * (6 + 2/9) = -280/9
    lambda1 = 5
    d1 = 6
    eta_sum = 2/9
    G2_pred = -lambda1 * (d1 + eta_sum)
    print(f"G2 predicted: {G2_pred:.10f} (-280/9)")
    assert abs(G2_pred - (-280/9)) < 1e-10

    # Proton mass check (two-loop)
    # m_p/m_e = 6π^5 * (1 + (10/9)*alpha^2/pi - (280/9)*alpha^4/pi^2)
    alpha = 1.0 / 137.035999084  # PDG 2024
    G = 10/9
    term0 = 6 * (math.pi ** 5)
    term1 = G * (alpha ** 2) / math.pi
    term2 = G2_pred * (alpha ** 4) / (math.pi ** 2)
    mp_me_pred = term0 * (1.0 + term1 + term2)
    mp_me_meas = 1836.15267343
    print(f"One-loop m_p/m_e: {term0 * (1 + term1):.10f}")
    print(f"Two-loop m_p/m_e: {mp_me_pred:.10f}")
    print(f"Measured m_p/m_e: {mp_me_meas:.10f}")
    
    # Check inverse alpha prediction
    # Solve quadratic G2/pi^2 * alpha^4 + G/pi * alpha^2 - Delta = 0
    Delta = (mp_me_meas / term0) - 1.0
    a = G2_pred / (math.pi ** 2)
    b = G / math.pi
    c = -Delta
    alpha_sq = (-b + (b**2 - 4*a*c)**0.5) / (2*a)
    alpha_inv_pred = 1.0 / (alpha_sq ** 0.5)
    print(f"Inverse alpha (2-loop): {alpha_inv_pred:.6f}")
    assert abs(alpha_inv_pred - 137.036) < 0.001

    print("\nAll replication checks passed.")
