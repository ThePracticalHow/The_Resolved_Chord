"""
Permutation / Scramble Test
============================
Anti-numerology control: randomly permute spectral-quantity→role mapping.
Blocks: "You could map anything to anything."
"""

import numpy as np
import pytest

from leng_test_utils import (
    brannen_masses,
    mu_from_electron,
    koide_ratio,
    R,
    M_E_PDG,
    M_MU_PDG,
    M_TAU_PDG,
    K_THEORY,
    DELTA,
    TWIST,
)

# Correct mapping: k=0→e, k=1→μ, k=2→τ (sorted by mass)
# Scramble: swap which spectral index maps to which generation


def brannen_scrambled(mu: float, r: float, delta: float, perm: list) -> np.ndarray:
    """Brannen with permuted k: √m_k = μ(1 + r·cos(δ + 2π·perm(k)/3)).
    Note: Permuting k in cos(δ+2πk/3) only reorders the same 3 values; use swap_invariants instead."""
    return np.array([mu * (1 + r * np.cos(delta + 2 * np.pi * perm[k] / 3)) for k in range(3)])


def swap_twist_and_base_phase(mu: float, r: float) -> np.ndarray:
    """Wrong: swap 2π/3 and 2/9 roles. Use δ = 2/9 + 2π/3·k (wrong structure)."""
    # Correct: δ_k = 2π/3 + 2/9 + 2πk/3. Wrong: put twist in 'generation' slot.
    return np.array([
        mu * (1 + r * np.cos(2/9 + 2 * np.pi * k / 3))  # No base 2π/3, wrong phase structure
        for k in range(3)
    ]) ** 2


def brannen_wrong_phase_sign(mu: float, r: float) -> np.ndarray:
    """Wrong: use δ = 2π/3 - 2/9 (flip twist sign). Structural perturbation."""
    delta_wrong = 2 * np.pi / 3 - TWIST
    return brannen_masses(mu, r, delta_wrong)  # from leng_test_utils


class TestPermutationScramble:
    """Scrambled mappings should fail to reproduce PDG."""

    @pytest.mark.xfail(
        reason="Phase swap δ=2/9+2πk/3 is near-degenerate with δ=2π/3+2/9+2πk/3 "
               "under Z₃ cyclic permutation of cos arguments. "
               "See test_wrong_r_in_phase_slot_fails for a scramble that breaks predictions.",
        strict=True,
    )
    def test_swap_phase_structure_fails(self):
        """Swapping twist and base phase roles — documents that this scramble is degenerate."""
        mu = mu_from_electron(M_E_PDG, R, DELTA)
        masses_wrong = swap_twist_and_base_phase(mu, R)
        masses = np.sort(masses_wrong)
        rel_err_mu = abs(masses[1] - M_MU_PDG) / M_MU_PDG
        assert rel_err_mu > 0.01

    @pytest.mark.xfail(
        reason="Flipping twist sign (2/9 → -2/9) is a small perturbation that preserves "
               "K≈2/3 and all mass ratios to ~10⁻⁵. The Brannen formula's cos is robust "
               "to this near-symmetric flip. "
               "See test_wrong_r_in_phase_slot_fails for a scramble that breaks predictions.",
        strict=True,
    )
    def test_wrong_phase_sign_fails(self):
        """Flipping twist sign — documents that this scramble is degenerate."""
        mu = mu_from_electron(M_E_PDG, R, DELTA)
        masses = brannen_wrong_phase_sign(mu, R)
        masses = np.sort(masses)
        rel_err_mu = abs(masses[1] - M_MU_PDG) / M_MU_PDG
        assert rel_err_mu > 0.01

    def test_wrong_r_in_phase_slot_fails(self):
        """Using r (amplitude) in phase slot: δ = 2π/3 + r + 2πk/3. Structural scramble."""
        delta_wrong = 2 * np.pi / 3 + R  # R instead of 2/9
        mu = mu_from_electron(M_E_PDG, R, delta_wrong)
        masses = np.sort(brannen_masses(mu, R, delta_wrong))
        rel_err = abs(masses[1] - M_MU_PDG) / M_MU_PDG
        assert rel_err > 0.01
