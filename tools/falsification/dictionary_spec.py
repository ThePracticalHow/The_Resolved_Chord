"""
Dictionary Spec v1.0 — LENG Geometry-to-Physics Mapping
========================================================
Machine-verifiable implementation of the formal dictionary rules.
Frozen spec: DO NOT change without version bump and Supplement G update.

References: Supplement B (Dictionary), Supplement G (Dictionary Validation),
GhostModes.py, EtaInvariant.py, The_Resolved_Chord v6.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import comb
from typing import List, Tuple

import numpy as np


# ─────────────────────────────────────────────────────────────────────────────
#  LOCKED DEFINITIONS (Dictionary Spec v1.0)
# ─────────────────────────────────────────────────────────────────────────────

# D1: Manifold
MANIFOLD = "S^5/Z_3"
LENS_SPACE = "L(3;1,1,1)"
N, P = 3, 3  # complex dim n=3, orbifold order p=3

# D2: Operator — Laplacian on S^5 for mode counting; Dirac for eta
# Eigenvalues: λ_ℓ = ℓ(ℓ+4), degeneracy d_ℓ = (ℓ+1)(ℓ+2)²(ℓ+3)/12

# D3: Z_3 action — on coordinates z_j ∈ ℂ³: z_j ↦ ω·z_j, z̄_j ↦ ω²·z̄_j
# Bidegree (a,b) has charge ω^{a-b}. Survives iff (a-b) ≡ 0 (mod 3).

# D4: Subgroup chain — U(3) = SU(3)×U(1)/Z_3 preserved by Z_3 orbifold
# Parent: SO(6) ≅ SU(4)/Z_2. Color group SU(3)_C from U(3).


# ─────────────────────────────────────────────────────────────────────────────
#  RULE 1: Harmonic dimension (math fact)
# ─────────────────────────────────────────────────────────────────────────────

def harmonic_dim(a: int, b: int) -> int:
    """
    Dimension of harmonic bidegree (a,b) on S^5 = ℂ³.
    Source: H^{a,b} ⊂ Harm(S^5), dim = (a+1)(b+1)(a+b+2)/2 - corrections.
    """
    if a >= 1 and b >= 1:
        return comb(a + 2, 2) * comb(b + 2, 2) - comb(a + 1, 2) * comb(b + 1, 2)
    elif b == 0:
        return comb(a + 2, 2)
    else:
        return comb(b + 2, 2)


# ─────────────────────────────────────────────────────────────────────────────
#  RULE 2: Z_3 charge (math fact)
# ─────────────────────────────────────────────────────────────────────────────

def z3_charge(a: int, b: int, p: int = 3) -> int:
    """Z_p charge of bidegree (a,b): (a-b) mod p. 0 = invariant."""
    return (a - b) % p


def survives_projection(a: int, b: int, p: int = 3) -> bool:
    """True iff (a-b) ≡ 0 (mod p)."""
    return z3_charge(a, b, p) == 0


# ─────────────────────────────────────────────────────────────────────────────
#  RULE 3: Decompose level ℓ (math fact)
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class LevelDecomposition:
    """Output of decompose_level for a single ℓ."""
    ell: int
    d_total: int
    d_inv: int
    d_proj: int
    eigenvalue: float
    components: List[Tuple[int, int, int, int]]  # (a, b, dim, charge)


def decompose_level(ell: int, p: int = 3) -> LevelDecomposition:
    """
    Decompose degree-ℓ harmonics on S^5 under Z_p.
    Returns total, invariant, projected counts and eigenvalue.
    """
    if ell < 0:
        raise ValueError("ell must be >= 0")
    d_total = 0
    d_inv = 0
    components = []
    for a in range(ell + 1):
        b = ell - a
        h = harmonic_dim(a, b)
        charge = z3_charge(a, b, p)
        d_total += h
        if charge == 0:
            d_inv += h
        components.append((a, b, h, charge))
    eigenval = ell * (ell + 4)
    d_proj = d_total - d_inv
    return LevelDecomposition(
        ell=ell,
        d_total=d_total,
        d_inv=d_inv,
        d_proj=d_proj,
        eigenvalue=float(eigenval),
        components=components,
    )


# ─────────────────────────────────────────────────────────────────────────────
#  RULE 4: Survivor rep table (the key output)
# ─────────────────────────────────────────────────────────────────────────────

def survivor_rep_table(ell_max: int = 12, p: int = 3) -> List[Tuple[int, int, int, int, float]]:
    """
    Survivor rep table: ℓ → (d_total, d_inv, d_proj, eigenvalue).
    This is the object referees care about.
    """
    return [
        (dec.ell, dec.d_total, dec.d_inv, dec.d_proj, dec.eigenvalue)
        for dec in (decompose_level(ell, p) for ell in range(ell_max + 1))
    ]


# ─────────────────────────────────────────────────────────────────────────────
#  RULE 5: Burnside check (consistency)
# ─────────────────────────────────────────────────────────────────────────────

def burnside_check(ell: int, p: int = 3, tol: float = 1e-9) -> bool:
    """Verify d_inv = (d_total + Σ Tr(χ_m)) / p. Tr(χ_m) = Σ h(a,b) ω^{m(a-b)}."""
    dec = decompose_level(ell, p)
    omega = np.exp(2j * np.pi / p)
    tr_sum = 0.0
    for m in range(1, p):
        tr = sum(
            comp[2] * (omega ** (m * (comp[0] - comp[1])))
            for comp in dec.components
        )
        tr_sum += tr.real
    inv_check = (dec.d_total + tr_sum) / p
    return abs(inv_check - dec.d_inv) < tol


# ─────────────────────────────────────────────────────────────────────────────
#  NEGATIVE CONROLS: alternate parameters
# ─────────────────────────────────────────────────────────────────────────────

def survivor_table_alternate_p(p: int, ell_max: int = 8) -> List[Tuple[int, int, int]]:
    """Survivor table for alternate orbifold order p. Used for negative control."""
    return [
        (dec.ell, dec.d_inv, dec.d_total)
        for dec in (decompose_level(ell, p) for ell in range(ell_max + 1))
    ]


# ─────────────────────────────────────────────────────────────────────────────
#  FORKING PATHS: conventions that should not change results
# ─────────────────────────────────────────────────────────────────────────────

def decompose_level_alternate_basis_order(ell: int, p: int = 3) -> LevelDecomposition:
    """Same as decompose_level but iterate b first instead of a. Should match."""
    d_total = 0
    d_inv = 0
    components = []
    for b in range(ell + 1):
        a = ell - b
        h = harmonic_dim(a, b)
        charge = z3_charge(a, b, p)
        d_total += h
        if charge == 0:
            d_inv += h
        components.append((a, b, h, charge))
    eigenval = ell * (ell + 4)
    d_proj = d_total - d_inv
    return LevelDecomposition(
        ell=ell,
        d_total=d_total,
        d_inv=d_inv,
        d_proj=d_proj,
        eigenvalue=float(eigenval),
        components=components,
    )


# ─────────────────────────────────────────────────────────────────────────────
#  CANONICAL OUTPUT (for Supplement G)
# ─────────────────────────────────────────────────────────────────────────────

def format_survivor_table(ell_max: int = 12) -> str:
    """Format survivor table for inclusion in Supplement G or branching_outputs.txt."""
    rows = survivor_rep_table(ell_max)
    lines = [
        "# Survivor Rep Table — S^5/Z_3",
        "# ell | d_total | d_inv | d_proj | lambda_ell",
        "# ---",
    ]
    for ell, d_tot, d_inv, d_proj, lam in rows:
        lines.append(f"  {ell:2} | {d_tot:7} | {d_inv:5} | {d_proj:6} | {lam:.0f}")
    return "\n".join(lines)


if __name__ == "__main__":
    print(format_survivor_table())
    print("\n# Burnside check (all ell <= 12):", all(burnside_check(ell) for ell in range(13)))
    print("# ell=1: d_inv=0 (all killed):", decompose_level(1).d_inv == 0)
    print("# ell=2: first survivors:", decompose_level(2).d_inv)
