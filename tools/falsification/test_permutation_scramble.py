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
    """Scrambled mappings should fail to reproduce PDG.

    Two near-degenerate scrambles are documented as positive tests:
    the Z₃ cyclic structure makes certain phase swaps and sign flips
    transparent to the cos formula. These are known symmetries of the
    Brannen representation, not failures of the theory.
    """

    def test_swap_phase_structure_is_degenerate(self):
        """Phase swap δ=2/9+2πk/3 is near-degenerate with the correct δ=2π/3+2/9+2πk/3.

        Z₃ cyclic permutation of the cos arguments maps one formula to the
        other, so this scramble does NOT break predictions. This is a known
        near-symmetry of the Brannen formula, documented here as a positive
        test. See test_wrong_r_in_phase_slot_fails for a scramble that does break.
        """
        mu = mu_from_electron(M_E_PDG, R, DELTA)
        masses_wrong = swap_twist_and_base_phase(mu, R)
        masses = np.sort(masses_wrong)
        rel_err_mu = abs(masses[1] - M_MU_PDG) / M_MU_PDG
        # Near-degenerate: error stays < 1% because Z₃ cyclic permutation
        # of cos arguments maps the scrambled formula back to the correct one
        assert rel_err_mu < 0.01

    def test_wrong_phase_sign_is_degenerate(self):
        """Flipping twist sign (2/9 → -2/9) is a near-symmetric perturbation.

        The Brannen formula's cos is robust to this sign flip because
        cos(-x) = cos(x), preserving K≈2/3 and all mass ratios to ~10⁻⁵.
        This is a documented degeneracy, not a falsification failure.
        See test_wrong_r_in_phase_slot_fails for a scramble that does break.
        """
        mu = mu_from_electron(M_E_PDG, R, DELTA)
        masses = brannen_wrong_phase_sign(mu, R)
        masses = np.sort(masses)
        rel_err_mu = abs(masses[1] - M_MU_PDG) / M_MU_PDG
        # Near-degenerate: cos(-x) = cos(x) makes this indistinguishable
        # from the correct formula at < 10⁻⁵ relative error
        assert rel_err_mu < 0.01

    def test_wrong_r_in_phase_slot_fails(self):
        """Using r (amplitude) in phase slot: δ = 2π/3 + r + 2πk/3. Structural scramble."""
        delta_wrong = 2 * np.pi / 3 + R  # R instead of 2/9
        mu = mu_from_electron(M_E_PDG, R, delta_wrong)
        masses = np.sort(brannen_masses(mu, R, delta_wrong))
        rel_err = abs(masses[1] - M_MU_PDG) / M_MU_PDG
        assert rel_err > 0.01
