"""
lotus.landscape — The Universe Selection Pipeline
=====================================================

Scans all possible S^(2n-1)/Z_p universes and classifies them.
Only S^5/Z_3 survives: the rest are dead.

    >>> from lotus.landscape import survey, selection_pipeline
    >>> survey()                   # returns list of all candidate universes
    >>> selection_pipeline()       # prints the 96 → 2 → 1 ablation
    >>> death_certificate(4, 2)    # why S^7/Z_2 is dead

Ported from verification/UniverseLandscape.py.
"""

import math

PI = math.pi


# ── Core formulas for any (n, p) universe ──────────────────────

def _is_prime(n: int) -> bool:
    """Check if n is prime."""
    if n < 2:
        return False
    if n < 4:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True


R = math.sqrt(2)  # Brannen r, forced by moment map geometry


def koide_K(p: int) -> float:
    """Natural Koide ratio for p generations with r=√2: K = 2/p."""
    return 2.0 / p


def twist_formula(n: int, p: int) -> float:
    """Predicted Koide twist: d₁ × τ_R = 2n / p^n."""
    return 2.0 * n / p ** n


def delta_formula(n: int, p: int) -> float:
    """Full Brannen phase angle: δ = 2π/p + twist."""
    return 2 * PI / p + twist_formula(n, p)


def eta_invariant(n: int, p: int) -> float:
    """Donnelly eta invariant for S^(2n-1)/Z_p."""
    return n / p ** (n - 1)


def sqrt_masses(n: int, p: int) -> list:
    """Compute p square-root masses for universe (n, p).

    Uses the Brannen formula: √m_k = 1 + √2 · cos(δ + 2πk/p)
    """
    d = delta_formula(n, p)
    return [1.0 + R * math.cos(d + 2 * PI * k / p) for k in range(p)]


def all_positive(n: int, p: int) -> bool:
    """Check if all sqrt-mass modes are positive."""
    return all(m > 0 for m in sqrt_masses(n, p))


def resonance_gap(n: int, p: int) -> float:
    """Resonance gap: |n - p^(n-2)|. Zero means self-consistent."""
    if n < 2:
        return abs(n)
    return abs(n - p ** (n - 2))


def is_self_consistent(n: int, p: int) -> bool:
    """Check if (n, p) satisfies the resonance lock: n = p^(n-2)."""
    return resonance_gap(n, p) == 0


# ── Survey and selection pipeline ───────────────────────────────

def survey(n_max: int = 8, p_max: int = 13) -> list:
    """Survey all S^(2n-1)/Z_p candidate universes.

    Args:
        n_max: max sphere dimension parameter (default 8)
        p_max: max orbifold order (default 13)

    Returns:
        list of dicts, one per candidate universe, with classification
    """
    results = []
    for n in range(1, n_max + 1):
        for p in range(2, p_max + 1):
            eta = eta_invariant(n, p)
            K = koide_K(p)
            consistent = is_self_consistent(n, p)
            positive = all_positive(n, p) if consistent else False
            prime = _is_prime(p)

            # Classification
            if not consistent:
                status = 'DEAD: resonance lock fails'
            elif not positive:
                status = 'DEAD: negative mass modes'
            elif not prime:
                status = 'DEAD: non-prime orbifold'
            elif abs(K - 1.0) < 1e-9:
                status = 'DEAD: degenerate K=1'
            else:
                status = 'ALIVE'

            results.append({
                'n': n,
                'p': p,
                'manifold': f"S^{2*n-1}/Z_{p}",
                'eta': eta,
                'K': K,
                'self_consistent': consistent,
                'positive_masses': positive,
                'prime': prime,
                'status': status,
            })
    return results


def selection_pipeline(n_max: int = 8, p_max: int = 13, verbose: bool = True) -> dict:
    """Run the full selection pipeline, returning statistics.

    Stage 1: Resonance lock  n = p^(n-2)
    Stage 2: Positive masses  all √m_k ≥ 0
    Stage 3: Prime group      p is prime
    Stage 4: Non-degenerate   K ≠ 1

    Returns:
        dict with 'total', 'stage1', 'stage2', 'stage3', 'stage4', 'survivors'
    """
    total = n_max * (p_max - 1)

    stage1 = [(n, p)
              for n in range(1, n_max + 1)
              for p in range(2, p_max + 1)
              if is_self_consistent(n, p)]

    stage2 = [(n, p) for n, p in stage1 if all_positive(n, p)]
    stage3 = [(n, p) for n, p in stage2 if _is_prime(p)]
    stage4 = [(n, p) for n, p in stage3 if abs(koide_K(p) - 1.0) > 1e-9]

    result = {
        'total': total,
        'stage1_resonance': stage1,
        'stage2_positive': stage2,
        'stage3_prime': stage3,
        'stage4_survivors': stage4,
    }

    if verbose:
        print()
        print("  UNIVERSE SELECTION PIPELINE")
        print("  " + "=" * 50)
        print(f"  Total candidates:    {total}")
        print(f"  Stage 1 (resonance): {len(stage1)} survive  "
              f"({total - len(stage1)} eliminated)")
        for n, p in stage1:
            print(f"    -> S^{2*n-1}/Z_{p}  (n={n}, p={p})")
        print(f"  Stage 2 (positive):  {len(stage2)} survive  "
              f"({len(stage1) - len(stage2)} eliminated)")
        print(f"  Stage 3 (prime):     {len(stage3)} survive  "
              f"({len(stage2) - len(stage3)} eliminated)")
        print(f"  Stage 4 (K≠1):      {len(stage4)} survive  "
              f"({len(stage3) - len(stage4)} eliminated)")
        print()
        if stage4:
            n, p = stage4[0]
            print(f"  UNIQUE SURVIVOR: S^{2*n-1}/Z_{p}")
        else:
            print("  NO SURVIVORS")
        print()

    return result


def death_certificate(n: int, p: int) -> dict:
    """Why a specific universe (n, p) is dead.

    Args:
        n: sphere dimension parameter
        p: orbifold order

    Returns:
        dict with cause of death and details
    """
    eta = eta_invariant(n, p)
    manifold = f"S^{2*n-1}/Z_{p}"

    if not is_self_consistent(n, p):
        gap = resonance_gap(n, p)
        return {
            'manifold': manifold,
            'alive': False,
            'cause': 'Resonance lock failure',
            'detail': (f"n = p^(n-2) requires {n} = {p}^{n-2} = "
                       f"{p**(n-2) if n >= 2 else 'N/A'}. Gap = {gap}."),
            'eta': eta,
        }

    if not all_positive(n, p):
        masses = sqrt_masses(n, p)
        neg_modes = [(k, m) for k, m in enumerate(masses) if m < 0]
        return {
            'manifold': manifold,
            'alive': False,
            'cause': 'Negative mass modes',
            'detail': f"Modes with √m < 0: {neg_modes}",
            'eta': eta,
        }

    if not _is_prime(p):
        return {
            'manifold': manifold,
            'alive': False,
            'cause': 'Non-prime orbifold order',
            'detail': f"p = {p} is not prime.",
            'eta': eta,
        }

    return {
        'manifold': manifold,
        'alive': True,
        'cause': None,
        'detail': 'This universe is alive.',
        'eta': eta,
    }
