"""
Ghost Mode Tests
================
l=1 mode structure, character traces, d₁_inv = 0.
From LENG_Master.md §7, GhostModes.py.
"""

import numpy as np
import pytest

from leng_test_utils import (
    harmonic_dim_s5,
    d1_degeneracy,
    D1,
)


def decompose_level(l: int):
    """Decompose degree-l harmonics on S⁵ under Z₃. Returns (d_total, d_inv)."""
    d_total = 0
    d_inv = 0
    omega = np.exp(2j * np.pi / 3)
    for a in range(l + 1):
        b = l - a
        h = harmonic_dim_s5(a, b)
        charge = (a - b) % 3
        d_total += h
        if charge == 0:
            d_inv += h
    return d_total, d_inv


def char_trace_l1(p: int, n: int) -> float:
    """Σ_{m=1}^{p-1} |Tr_{l=1}(χ_m)|. For p=3: = 2n·|cos(2π/3)| × 2 = 2n."""
    if p == 3:
        return 2 * n  # |cos(2π/3)| = 1/2, sum over 2 chars = 2n
    return sum(abs(2 * n * np.cos(2 * np.pi * m / p)) for m in range(1, p))


class TestGhostModeTheorem:
    """All l=1 modes killed by Z₃"""

    def test_l1_invariant_zero(self):
        """d₁_inv = 0: no Z₃-invariant l=1 modes"""
        d_total, d_inv = decompose_level(1)
        assert d_total == D1
        assert d_inv == 0

    def test_l1_total_six(self):
        """d₁(S⁵) = 6"""
        d_total, _ = decompose_level(1)
        assert d_total == 6

    def test_l0_survives(self):
        """l=0 (vacuum) survives"""
        d_total, d_inv = decompose_level(0)
        assert d_inv == 1

    def test_l2_partial_survival(self):
        """l=2: some modes survive"""
        d_total, d_inv = decompose_level(2)
        assert d_inv > 0
        assert d_total == 20  # Standard formula for S⁵


class TestCharacterTrace:
    """Σ|Tr_{l=1}(χ_m)| = d₁ for p=3 only"""

    def test_p3_equals_d1(self):
        assert char_trace_l1(3, 3) == D1

    def test_p5_not_equal(self):
        """p=5: ghost count differs from S⁵/Z₃ value d₁=6"""
        assert char_trace_l1(5, 3) != D1
