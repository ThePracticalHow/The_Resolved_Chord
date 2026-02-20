# Why Not Touch vev_overlap.py Directly

## What vev_overlap.py Does

The current `vev_overlap.py` is an **exploratory script** that tries many candidate formulas to derive v/m_p ≈ 262.4 from geometric overlap integrals on S⁵/Z₃. It:

1. **Tests multiple overlap definitions** (pairwise phase mismatch, triple amplitude, cross-mode triple, spectral zeta, volume×spectral combos)
2. **Searches for combinations** that might match v/m_p (e.g., 2/α − (d₁+λ₁+K), Vol(S⁴)×d₁λ₁/3, C(d₁,2)×λ₁×π)
3. **Runs Monte Carlo** on S⁵ to verify moments and explore sector-resolved overlap
4. **Contains the comment** at line 491: "Actually better: define by a Z₃-breaking quantity"

## The Paper's Actual Claim

The paper (P14, Supplement V) states:

- **v/m_p = 2/α − (d₁ + λ₁ + K)** — this is the *result*, not something to derive from overlap integrals
- The VEV is the **overlap amplitude** of three Z₃ sectors: ghost wavefunctions bleed through fold walls
- The Mexican hat is the **effective description** of overlap geometry in field-theory language

The paper does *not* give an explicit integral formula for the overlap. It gives the spectral-action derivation (higgs_vev_spectral_action.py) which yields v/m_p = 2/α − 35/3.

## Why Not Edit vev_overlap.py

1. **Exploratory value** — The script documents the search process. Deleting or rewriting it loses that history.
2. **The "Actually better" comment** — It suggests defining sectors by a Z₃-breaking quantity (e.g., arg(z₁) instead of z₁z₂z₃ phase). The current code uses `sector_phase = np.angle(z_mc[:, 0])` which *is* Z₃-breaking (z₁ shifts by ω under Z₃). So the "better" definition may already be partially implemented; the comment is a note about *why* that choice, not a missing implementation.
3. **No single "correct" implementation** — The paper gives the formula; it doesn't give the explicit overlap integral. vev_overlap explores candidates. A "paper-accurate" version would *start* from the formula and show how the overlap *interpretation* fits, not search for the formula.

## What vev_overlap_paper.py Will Do

A new script that:
- Starts from the paper's formula v/m_p = 2/α − (d₁ + λ₁ + K)
- Verifies S⁵ moments as the geometric basis
- Implements the overlap *interpretation*: sector labeling by Z₃-breaking arg(z₁), fold-wall overlap of ghost modes
- Shows the connection to the spectral action (higgs_vev_spectral_action.py)
- Does NOT search for formulas — it verifies and interprets the one the paper gives
