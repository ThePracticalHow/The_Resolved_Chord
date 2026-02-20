"""
Dictionary Validation Suite — LENG
=================================
Consolidated tests for Supplement G: Dictionary Validation Suite.

1. Negative Controls: Wrong parameters (alternate p, composite p) must change
   the survivor table or fail key claims.
2. Forking Paths: Conventions that shouldn't matter (basis order, tolerance)
   must leave results invariant.
3. Canonical Values: Frozen outputs must match Supplement G / branching_outputs.txt.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dictionary_spec import (
    decompose_level,
    survivor_rep_table,
    survivor_table_alternate_p,
    decompose_level_alternate_basis_order,
    burnside_check,
)


# ─────────────────────────────────────────────────────────────────────────────
#  NEGATIVE CONTROLS (Must Fail / Change Results)
# ─────────────────────────────────────────────────────────────────────────────

class TestDictionaryNegativeControls:
    """Alternate parameters must change survivor table or fail key claims."""

    def test_p3_ell1_all_killed(self):
        """Canonical: p=3, ell=1 has d_inv=0."""
        dec = decompose_level(1, p=3)
        assert dec.d_inv == 0, "ell=1 must have zero survivors for p=3"

    def test_p3_ell2_first_survivors(self):
        """Canonical: p=3, ell=2 has d_inv=8."""
        dec = decompose_level(2, p=3)
        assert dec.d_inv == 8, "ell=2 must have 8 survivors for p=3"

    def test_p2_vs_p3_ell2(self):
        """p=2 gives different survivor structure than p=3 at ell=2."""
        table_p3 = survivor_rep_table(ell_max=2, p=3)
        table_p2 = survivor_rep_table(ell_max=2, p=2)
        d_inv_p3_ell2 = table_p3[2][2]
        d_inv_p2_ell2 = table_p2[2][2]
        assert d_inv_p3_ell2 == 8
        assert d_inv_p2_ell2 == 20
        assert d_inv_p3_ell2 != d_inv_p2_ell2

    def test_p2_table_differs_from_p3(self):
        """p=2 survivor table differs from p=3 across multiple levels."""
        table_p2 = survivor_table_alternate_p(p=2, ell_max=6)
        table_p3 = survivor_rep_table(ell_max=6, p=3)
        p2_inv = {r[0]: r[1] for r in table_p2}
        p3_inv = {r[0]: r[2] for r in table_p3}
        assert p2_inv != p3_inv, "p=2 and p=3 must differ"

    def test_p5_ell3_differs_from_p3(self):
        """p=5 at ell=3: d_inv=0 vs p=3 d_inv=20."""
        dec_p5 = decompose_level(3, p=5)
        dec_p3 = decompose_level(3, p=3)
        assert dec_p3.d_inv == 20
        assert dec_p5.d_inv == 0
        assert dec_p3.d_inv != dec_p5.d_inv

    def test_p5_ell1_total_modes(self):
        """p=5: ell=1 has d_total=6 (same as p=3); resonance lock fails for p!=3."""
        dec = decompose_level(1, p=5)
        assert dec.d_total == 6, "ell=1 always has 6 modes on S^5"

    def test_p4_vs_p3_ell4(self):
        """p=4 (composite) gives different d_inv than p=3 at ell=4."""
        d_inv_p3 = decompose_level(4, p=3).d_inv
        d_inv_p4 = decompose_level(4, p=4).d_inv
        assert d_inv_p3 == 27
        assert d_inv_p4 == 57
        assert d_inv_p3 != d_inv_p4

    def test_p7_survivor_table_differs(self):
        """p=7 survivor table differs from p=3."""
        table_p7 = survivor_table_alternate_p(p=7, ell_max=8)
        table_p3 = survivor_rep_table(ell_max=8, p=3)
        p7_inv = {r[0]: r[1] for r in table_p7}
        p3_inv = {r[0]: r[2] for r in table_p3}
        assert p7_inv != p3_inv, "p=7 and p=3 survivor counts must differ"


# ─────────────────────────────────────────────────────────────────────────────
#  FORKING PATHS (Must Pass / Same Results)
# ─────────────────────────────────────────────────────────────────────────────

class TestDictionaryForkingPaths:
    """Basis order, tolerance, etc. must not change survivor multiplicities."""

    def test_alternate_basis_order_matches(self):
        """Iterating (b,a) instead of (a,b) yields same d_inv, d_total."""
        for ell in range(10):
            dec1 = decompose_level(ell)
            dec2 = decompose_level_alternate_basis_order(ell)
            assert dec1.d_inv == dec2.d_inv, f"ell={ell} d_inv must match"
            assert dec1.d_total == dec2.d_total, f"ell={ell} d_total must match"
            assert dec1.d_proj == dec2.d_proj, f"ell={ell} d_proj must match"

    def test_burnside_all_levels(self):
        """Burnside check holds for all ell <= 12."""
        for ell in range(13):
            assert burnside_check(ell), f"Burnside failed at ell={ell}"

    def test_burnside_tolerance_invariance(self):
        """Stricter tolerance still passes (within float precision)."""
        for ell in range(8):
            assert burnside_check(ell, tol=1e-12), f"Burnside tol=1e-12 at ell={ell}"

    def test_survivor_table_deterministic(self):
        """Survivor table is deterministic across calls."""
        t1 = survivor_rep_table(ell_max=8)
        t2 = survivor_rep_table(ell_max=8)
        assert t1 == t2, "Survivor table must be deterministic"


# ─────────────────────────────────────────────────────────────────────────────
#  CANONICAL VALUES (Frozen Outputs)
# ─────────────────────────────────────────────────────────────────────────────

class TestDictionaryCanonicalValues:
    """Canonical values must match Supplement G / branching_outputs.txt."""

    def test_canonical_values_frozen(self):
        """Canonical values match Supplement G / branching_outputs.txt."""
        table = survivor_rep_table(ell_max=12)
        assert table[0] == (0, 1, 1, 0, 0.0), "ell=0"
        assert table[1] == (1, 6, 0, 6, 5.0), "ell=1"
        assert table[2] == (2, 20, 8, 12, 12.0), "ell=2"
