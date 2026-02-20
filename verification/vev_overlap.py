"""
Higgs VEV as the Overlap of Three Z₃ Sectors on S⁵
===================================================

Physical picture (Jixiang's insight):
- Three sectors of S⁵/Z₃ are three "spheres" that tile the space
- The ℓ=1 ghost modes DON'T respect sector boundaries — they bleed through
- The VEV measures the overlap amplitude: how much the three sectors
  fail to be independent due to the ghost modes
- "Not a field acquiring a VEV — the geometry itself has an overlap"

The question: what integral gives v/m_e ≈ 481,843?

Key reframing: v/m_e = (v/m_p) × (m_p/m_e) = (v/m_p) × 6π⁵
So the VEV question REDUCES to: what is v/m_p ≈ 262.4?
"""


import sys, io
if sys.stdout.encoding and sys.stdout.encoding.lower() not in ('utf-8', 'utf8'):
    try:
        sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
        sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')
    except Exception:
        pass

import numpy as np
from scipy.special import gamma as Gamma
import warnings

def safe_divide(a, b, default=0):
    """Safe division to avoid division by zero."""
    try:
        return a / b if b != 0 else default
    except (ZeroDivisionError, OverflowError):
        return default

def S5_moment(a1: int, a2: int, a3: int) -> float:
    """Normalized moment ⟨|z₁|^{2a₁}|z₂|^{2a₂}|z₃|^{2a₃}⟩ on S⁵.

    For S^{2n-1} with n complex coordinates, the formula is:
      (n-1)! × Γ(a₁+1)Γ(a₂+1)...Γ(aₙ+1) / Γ(a₁+a₂+...+aₙ+n)
    For S⁵ (n=3): prefactor is 2! = 2.
    """
    if any(x < 0 for x in [a1, a2, a3]):
        warnings.warn("Negative exponents not supported")
        return 0.0
    try:
        # n=3 complex coordinates on S^5, prefactor = (n-1)! = 2! = 2
        return float(2 * Gamma(a1+1) * Gamma(a2+1) * Gamma(a3+1) / Gamma(a1+a2+a3+3))
    except (OverflowError, ValueError):
        return 0.0

# ===================== PHYSICAL TARGETS =====================
v_higgs = 246.2196       # GeV (PDG)
m_e = 0.51099895e-3      # GeV
m_p = 0.93827208943      # GeV
alpha_em = 1/137.035999206
y_e = np.sqrt(2) * m_e / v_higgs

v_over_me = v_higgs / m_e     # 481,843
v_over_mp = v_higgs / m_p     # 262.4
mp_over_me = m_p / m_e        # 1836.15

print("=" * 70)
print("HIGGS VEV AS GHOST MODE OVERLAP ON S⁵/Z₃")
print("=" * 70)
print(f"\n  v = {v_higgs:.4f} GeV")
print(f"  m_e = {m_e*1e3:.6f} MeV")
print(f"  m_p = {m_p*1e3:.4f} MeV")
print(f"  v/m_e = {v_over_me:.1f}")
print(f"  v/m_p = {v_over_mp:.4f}")
print(f"  m_p/m_e = {mp_over_me:.5f}  (6π⁵ = {6*np.pi**5:.5f})")
print(f"  y_e = {y_e:.6e}")

# ===================== SPECTRAL DATA =====================
# S⁵ eigenvalues: λ_l = l(l+4), degeneracy d_l = (l+1)(l+2)²(l+3)/12 for scalar
# ℓ=1 ghost sector: λ₁=5, d₁=6 (killed by Z₃)
lambda1 = 5
d1 = 6
omega = np.exp(2j * np.pi / 3)

# Volumes
vol_S5 = np.pi**3                  # Vol(S⁵) = π³
vol_S4 = 8 * np.pi**2 / 3         # Vol(S⁴) = 8π²/3
vol_S5_Z3 = vol_S5 / 3            # Vol(S⁵/Z₃) = π³/3
vol_S3 = 2 * np.pi**2             # Vol(S³) = 2π²

print("\n" + "=" * 70)
print("PART 1: THE KEY RATIO v/m_p")
print("=" * 70)
print(f"\nSince m_p/m_e = 6π⁵, we have v/m_e = (v/m_p) × 6π⁵")
print(f"So the WHOLE VEV question reduces to: what is v/m_p = {v_over_mp:.4f}?")

# ===================== EXACT INTEGRALS ON S⁵ =====================
# Formula: ∫_{S^{2n-1}} |z₁|^{2a₁}|z₂|^{2a₂}|z₃|^{2a₃} dΩ / Vol(S^{2n-1})
#   = Γ(a₁+1)Γ(a₂+1)Γ(a₃+1) / Γ(a₁+a₂+a₃+3)
# For n=3 (S⁵)

def S5_moment_raw(a1, a2, a3):
    """Unnormalized moment Γ(a₁+1)Γ(a₂+1)Γ(a₃+1) / Γ(a₁+a₂+a₃+3) on S⁵"""
    return Gamma(a1+1)*Gamma(a2+1)*Gamma(a3+1) / Gamma(a1+a2+a3+3)

print("\n--- Exact moments on S⁵ ---")
print(f"⟨|z_j|²⟩ = {S5_moment(1,0,0):.6f}  (= 1/3)")
print(f"⟨|z_j|⁴⟩ = {S5_moment(2,0,0):.6f}  (= 1/6)")
print(f"⟨|z₁|²|z₂|²⟩ = {S5_moment(1,1,0):.6f}  (= 1/12)")
print(f"⟨|z₁|²|z₂|²|z₃|²⟩ = {S5_moment(1,1,1):.6f}  (= 1/60)")
print(f"⟨|z_j|⁶⟩ = {S5_moment(3,0,0):.6f}  (= 1/10)")

# ===================== OVERLAP COMPUTATIONS =====================
print("\n" + "=" * 70)
print("PART 2: OVERLAP INTEGRALS")
print("=" * 70)

# OVERLAP DEFINITION 1: Pairwise phase mismatch
# Ghost mode ψ = z_j acquires phase ω at each fold wall
# Pairwise overlap: Re(ψ₀*ψ₁)/|ψ|² = cos(2π/3) = -1/2
print("\n--- Def 1: Pairwise phase mismatch ---")
cos_mismatch = np.cos(2*np.pi/3)  # = -1/2
print(f"Phase mismatch per wall: cos(2π/3) = {cos_mismatch:.4f}")
print(f"Josephson energy per wall: 1-cos(2π/3) = {1-cos_mismatch:.4f} = 3/2")
print(f"Total (3 walls): {3*(1-cos_mismatch):.4f} = 9/2")
print(f"  → Too small for v/m_p")

# OVERLAP DEFINITION 2: Triple overlap integral
# ∫_{S⁵} (sector₀ wavefunction)(sector₁ wf)(sector₂ wf) 
# For ghost mode z_j: the "sector k piece" is z_j × indicator(sector k)
# But adjacent sectors DON'T overlap spatially → need wavefunction overlap
# The triple AMPLITUDE overlap: ψ(p)·ψ(Rp)·ψ(R²p) for ghost modes
# For ψ = z_j: z_j · ωz_j · ω²z_j = ω³ z_j³ = z_j³ (since ω³=1)
print("\n--- Def 2: Triple amplitude overlap ---")
print("⟨ψ·Rψ·R²ψ⟩ = ⟨z_j³⟩ = ∫z_j³ dΩ/Vol = 0 (by phase symmetry)")
print("  → Zero. Need |ψ|² version.")

# Triple intensity overlap: |ψ(p)|² × |ψ(Rp)|² × |ψ(R²p)|²
# For z_j: |z_j|² × |ωz_j|² × |ω²z_j|² = |z_j|⁶
triple_intensity = S5_moment(3, 0, 0)  # = 1/10
print(f"\n⟨|ψ|²|Rψ|²|R²ψ|²⟩ = ⟨|z_j|⁶⟩ = {triple_intensity:.6f} = 1/10")
print(f"Inverse: 1/triple = {1/triple_intensity:.1f}")

# OVERLAP DEFINITION 3: Cross-mode triple overlap
# What about the COMBINED overlap of all 6 ghost modes?
# This is ⟨|z₁|²|z₂|²|z₃|²⟩ = 1/60
cross_triple = S5_moment(1, 1, 1)  # = 1/60
print(f"\n--- Def 3: Cross-mode triple overlap ---")
print(f"⟨|z₁|²|z₂|²|z₃|²⟩ = {cross_triple:.6f} = 1/60")
print(f"Inverse: {1/cross_triple:.1f}")

# The 1/60 is interesting: 60 = 5! / 2 = |A₅| (alternating group)
# Also 60 = Vol(S⁵)/(π³/60) ... 

# OVERLAP DEFINITION 4: Spectral weight of ghost sector
# The ghost modes contribute to the heat kernel / spectral zeta
# ζ_ghost(s) = d₁ × λ₁^{-s} = 6/5^s
print(f"\n--- Def 4: Spectral zeta of ghost sector ---")
for s_val in [-3, -5/2, -2, -3/2, -1, -1/2, 0, 1/2, 1]:
    zeta = d1 * lambda1**(-s_val)
    print(f"  ζ_ghost({s_val:+.1f}) = 6 × 5^{-s_val:+.1f} = {zeta:.4f}")

# OVERLAP DEFINITION 5: Volume × spectral weight combinations
print(f"\n--- Def 5: Volume × spectral combinations ---")
combos = {}
combos["d₁ × λ₁"] = d1 * lambda1                    # 30
combos["d₁ × λ₁²"] = d1 * lambda1**2                # 150
combos["d₁ × λ₁² / (2π)"] = d1 * lambda1**2 / (2*np.pi)
combos["Vol(S⁵)/Vol(S⁴) × d₁λ₁"] = (vol_S5/vol_S4) * d1*lambda1
combos["(π/8) × d₁λ₁"] = (np.pi/8) * d1 * lambda1
combos["d₁² × λ₁²/6"] = d1**2 * lambda1**2 / 6

# The "overlap volume" interpretation
# Three S⁴ fold walls in S⁵, each with area Vol(S⁴) = 8π²/3
# Their pairwise intersection: S³ with volume 2π²
# Their triple intersection: S² with volume 4π
# Ratios:
combos["Vol(S⁵)/Vol(triple intersection S²)"] = vol_S5 / (4*np.pi)
combos["Vol(S⁵)/(Vol(S²)×Vol(S³))"] = vol_S5 / (4*np.pi * 2*np.pi**2)

# Ghost mode × fold wall combinations
combos["d₁ × Vol(S⁴)/Vol(S⁵/Z₃)"] = d1 * vol_S4 / vol_S5_Z3
combos["d₁λ₁ × Vol(S⁴)"] = d1 * lambda1 * vol_S4
combos["d₁λ₁ × Vol(S⁵)"] = d1 * lambda1 * vol_S5
combos["d₁ × Vol(S⁵)"] = d1 * vol_S5
combos["λ₁ × Vol(S⁵)²"] = lambda1 * vol_S5**2
combos["d₁λ₁ × π²"] = d1 * lambda1 * np.pi**2
combos["d₁λ₁ × (π²-5)"] = d1 * lambda1 * (np.pi**2 - 5)
combos["6 × 5 × π²"] = 6 * 5 * np.pi**2

# Key insight: v/m_p might decompose through α
combos["1/(2α)"] = 1/(2*alpha_em)           # ≈ 68.5
combos["1/α"] = 1/alpha_em                  # ≈ 137
combos["2/α"] = 2/alpha_em                  # ≈ 274 (close to 262!)
combos["(2/α) × (1 - 2π·α)"] = 2/alpha_em * (1 - 2*np.pi*alpha_em)
combos["π/6α²"] = np.pi/(6*alpha_em**2)

# Direct hit attempts for v/m_p ≈ 262.4
combos["2/α - 12"] = 2/alpha_em - 12        # ≈ 262.1
combos["2/α - 2·d₁"] = 2/alpha_em - 2*d1    # = 274.07 - 12 = 262.07
combos["2/α - 2d₁"] = 2/alpha_em - 2*d1
combos["(2/α)(1 - d₁α)"] = (2/alpha_em)*(1 - d1*alpha_em)

print(f"{'Quantity':<45} {'Value':>12} {'v/m_p ratio':>12} {'% error':>10}")
print("-" * 79)
for name, val in sorted(combos.items(), key=lambda x: abs(np.log(abs(x[1]/v_over_mp)) if x[1]>0 else 999)):
    if val > 0:
        ratio = v_over_mp / val
        pct = (ratio - 1) * 100
        marker = " <<<" if abs(pct) < 1 else (" **" if abs(pct) < 5 else "")
        print(f"{name:<45} {val:>12.4f} {ratio:>12.6f} {pct:>+9.4f}%{marker}")

# ===================== THE 2/α - 2d₁ HIT =====================
print("\n" + "=" * 70)
print("PART 3: THE 2/α − 2d₁ CANDIDATE")
print("=" * 70)

candidate = 2/alpha_em - 2*d1
print(f"\n  2/α = {2/alpha_em:.6f}")
print(f"  2d₁ = {2*d1}")
print(f"  2/α - 2d₁ = {candidate:.6f}")
print(f"  v/m_p (measured) = {v_over_mp:.6f}")
print(f"  Error: {(candidate/v_over_mp - 1)*100:+.4f}%")

# What this would mean:
# v/m_p = 2/α - 2d₁ = 2(1/α - d₁)
# v = m_p × 2(1/α - 6) = m_p × 2(137.036 - 6) = m_p × 262.072
# v/m_e = 6π⁵ × 2(1/α - 6) 
v_over_me_pred = 6*np.pi**5 * 2*(1/alpha_em - d1)
print(f"\n  v/m_e predicted = 6π⁵ × 2(1/α - d₁)")
print(f"                  = {v_over_me_pred:.1f}")
print(f"  v/m_e measured  = {v_over_me:.1f}")
print(f"  Error: {(v_over_me_pred/v_over_me - 1)*100:+.4f}%")

# What does this MEAN?
# 1/α is the electromagnetic coupling wall (137th prime wall)
# d₁ = 6 is the ghost mode degeneracy
# The VEV = 2 × (electromagnetic wall - ghost degeneracy) × m_p
# "How much of the EM wall is left after the ghosts eat their share"

# ===================== ALTERNATIVE: PURE SPECTRAL =====================
print("\n" + "=" * 70)
print("PART 4: PURE SPECTRAL APPROACHES (no α input)")
print("=" * 70)

# Can we get v/m_p from ONLY the S⁵/Z₃ spectral data?
# Available: d₁=6, λ₁=5, Vol(S⁵)=π³, Vol(S⁵/Z₃)=π³/3

spectral = {}
# The heat kernel at the ghost sector
# K_ghost(t) = d₁ × exp(-λ₁t) / Vol(S⁵/Z₃)
# At what "time" t does this equal 1/v_over_mp?

# Try: v/m_p = d₁ × λ₁^{n/2} for some dimension-dependent power
for n_half in np.arange(1, 8, 0.5):
    val = d1 * lambda1**(n_half)
    spectral[f"d₁ × 5^{n_half:.1f}"] = val

# Try: involving π from the bulk
for k in range(1, 6):
    spectral[f"d₁ × π^{k}"] = d1 * np.pi**k
    spectral[f"λ₁ × π^{k}"] = lambda1 * np.pi**k
    spectral[f"d₁λ₁ × π^{k}/6"] = d1*lambda1*np.pi**k/6

# Try: combinations with the overlap integrals
spectral["1/⟨|z|⁶⟩ × d₁"] = (1/triple_intensity) * d1   # 60
spectral["1/⟨|z₁z₂z₃|²⟩ × d₁/3"] = 60 * d1 / 3         # 120
spectral["5! × d₁/λ₁"] = 120 * 6/5                       # 144
spectral["d₁! × λ₁"] = 720 * 5                            # 3600
spectral["d₁! / (d₁-3)!"] = 720/6                         # 120

# Combinatorial: number of ways to arrange ghost modes
spectral["C(d₁²,2)"] = d1**2 * (d1**2-1)/2               # 630
spectral["d₁⁴/d₁"] = d1**3                                # 216
spectral["(d₁λ₁)^{3/2}"] = (d1*lambda1)**1.5             # 164.3

print(f"{'Quantity':<40} {'Value':>12} {'v/m_p ratio':>12}")
print("-" * 64)
for name, val in sorted(spectral.items(), key=lambda x: abs(np.log(x[1]/v_over_mp)) if x[1]>0 else 999):
    ratio = v_over_mp / val
    if 0.1 < ratio < 10:
        pct = (ratio - 1)*100
        marker = " <<<" if abs(pct) < 2 else ""
        print(f"{name:<40} {val:>12.4f} {ratio:>12.6f}{marker}")

# ===================== FULL v/m_e COMBINATIONS =====================
print("\n" + "=" * 70)
print("PART 5: DIRECT v/m_e COMBINATIONS")  
print("=" * 70)

# v/m_e = √2/y_e, so this is equivalent to finding y_e
# y_e = 2.935 × 10⁻⁶ — extremely small
# In framework: y_e should come from Koide structure + some spectral factor

# The Koide parameters: K = 2/3, δ = 2π/3 + 2/9
delta_koide = 2*np.pi/3 + 2/9
K_koide = 2/3

# Koide gives mass RATIOS (m_μ/m_e, m_τ/m_e) but not the overall scale
# The overall scale IS the electron Yukawa, which IS v/m_e

# What if the overall scale comes from the ghost overlap?
# The TOTAL ghost contribution to the spectral action:
# S_ghost = d₁ × λ₁ × Vol(S⁵/Z₃) = 6 × 5 × π³/3 = 10π³
S_ghost = d1 * lambda1 * vol_S5_Z3
print(f"\nTotal ghost spectral action: d₁λ₁ × Vol(S⁵/Z₃) = {S_ghost:.4f}")
print(f"= 10π³ = {10*np.pi**3:.4f}")

# exp(S_ghost)?
print(f"exp(10π³) = {np.exp(10*np.pi**3):.4e}")
print(f"v/m_e = {v_over_me:.4e}")

# Too big. What about exp(d₁) or exp(λ₁)?
print(f"\nexp(d₁) = exp(6) = {np.exp(6):.4f}")
print(f"exp(λ₁) = exp(5) = {np.exp(5):.4f}")
print(f"exp(d₁+λ₁) = exp(11) = {np.exp(11):.4f}")
print(f"exp(d₁×λ₁/3) = exp(10) = {np.exp(10):.4f}")

# Hmm, v/m_e ≈ 4.8 × 10⁵. Log of that:
print(f"\nln(v/m_e) = {np.log(v_over_me):.6f}")
print(f"This is close to: 4π+1 = {4*np.pi+1:.4f}")
print(f"Or: 2π² = {2*np.pi**2:.4f}? No, that's {2*np.pi**2:.4f}")
print(f"Or: d₁+λ₁+ln(d₁λ₁) = {d1+lambda1+np.log(d1*lambda1):.4f}")

# ln(v/m_e) = 13.085
# 13.085 ≈ ? Let's see: d₁+λ₁ = 11, π²+2 = 11.87, 
# 4π = 12.566, 13 = 13, λ₁+8 = 13, d₁+7 = 13, 2d₁+1 = 13!
print(f"\n2d₁ + 1 = {2*d1+1}")
print(f"exp(2d₁+1) = exp(13) = {np.exp(13):.1f}")
print(f"v/m_e = {v_over_me:.1f}")
print(f"Ratio: {v_over_me/np.exp(13):.6f}")
# 442413 vs 481843. Not great.

# More precisely:
print(f"\nSearching for ln(v/m_e) = {np.log(v_over_me):.6f}")
# Try: ln(v/m_e) = λ₁/2 × ln(d₁λ₁) + something?
for a in range(1, 10):
    for b in range(0, 10):
        for c in range(-3, 4):
            val = a * np.log(lambda1) + b * np.log(d1) + c * np.log(np.pi)
            if abs(val - np.log(v_over_me)) < 0.01:
                print(f"  {a}·ln(5) + {b}·ln(6) + {c}·ln(π) = {val:.6f}")
            val2 = a * np.log(lambda1) + b * np.log(d1) + c
            if abs(val2 - np.log(v_over_me)) < 0.01:
                print(f"  {a}·ln(5) + {b}·ln(6) + {c} = {val2:.6f}")

# ===================== THE ELECTRON AS COUPLING =====================
print("\n" + "=" * 70)
print("PART 6: 'THE PROTON IS THE UNIVERSE OF THE ELECTRON'")
print("=" * 70)
print("""
Jixiang's insight: "How many electrons fit in a proton?" → 6π⁵
The proton is a geometric container. The electron is the coupling to it.
Then: "How many protons fit in the VEV?" → v/m_p = ???

If the VEV is the "universe of the proton" the same way
the proton is the "universe of the electron":
""")

# The hierarchy: m_e → m_p → v
# m_p/m_e = 6π⁵ (geometric volume of proton cavity)
# v/m_p = ??? (geometric volume of... the fold?)

# What if v/m_p involves the SAME spectral data but at the NEXT level?
# Level 1: single mode (ℓ=1) → 6π⁵ (d₁ × π^dim)
# Level 2: pairwise overlap → involves pairs of modes

# Number of pairs from d₁=6 modes: C(6,2) = 15
# 15 × λ₁ × π = 15 × 5 × π = 235.6 (not bad! 10% off from 262)
val_pairs = 15 * lambda1 * np.pi
print(f"C(d₁,2) × λ₁ × π = 15 × 5 × π = {val_pairs:.4f}")
print(f"v/m_p = {v_over_mp:.4f}, ratio: {v_over_mp/val_pairs:.4f}")

# C(6,2) × λ₁ × something = 262.4
needed = v_over_mp / 15 / lambda1
print(f"\nC(d₁,2)×λ₁ × X = v/m_p → X = {needed:.6f}")
print(f"Compare: π = {np.pi:.6f}, e = {np.e:.6f}")
print(f"X/π = {needed/np.pi:.6f}")

# What about the TRIPLE overlap (all three sectors simultaneously)?
# C(6,3) = 20 ... the ℓ=2 degeneracy!
print(f"\nC(d₁,3) = {int(d1*(d1-1)*(d1-2)/6)} = d₂ (ℓ=2 degeneracy = 20)")
val_triples = 20 * lambda1 * np.pi / 6 * 5
print(f"This connects ghost sector to ℓ=2 sector!")

# Let's try: d₂ × λ₁ × ... = 262.4?
d2 = 20
lambda2 = 12
print(f"\nd₂ = {d2}, λ₂ = {lambda2}")
print(f"d₂ × λ₂ = {d2*lambda2}")
print(f"d₂ × λ₂ + d₁ × λ₁ = {d2*lambda2 + d1*lambda1}")
print(f"(d₂ × λ₂ + d₁ × λ₁)/π = {(d2*lambda2 + d1*lambda1)/np.pi:.4f}")

# Total spectral weight of first two levels
total_01 = d1*lambda1 + d2*lambda2  # = 30 + 240 = 270
print(f"\nTotal spectral weight: d₁λ₁ + d₂λ₂ = {total_01}")
print(f"v/m_p = {v_over_mp:.4f}")
print(f"Difference: {total_01 - v_over_mp:.4f}")
print(f"That's {total_01 - v_over_mp:.4f} ≈ {2*np.pi + 1/np.pi:.4f} = 2π + 1/π?")

# Close! 270 vs 262.4, difference ~7.6
# What if it's d₁λ₁ + d₂λ₂ - correction?
# correction = 7.6 ≈ d₁+1 = 7? Or λ₁+d₁/2 = 8?

print("\n" + "=" * 70)
print("PART 7: THE 2/α CONNECTION (DEEPER)")
print("=" * 70)

# 2/α = 274.072
# v/m_p = 262.37
# Difference = 11.7 ≈ d₁ + λ₁ = 11? Or 2d₁ = 12?
diff_2a = 2/alpha_em - v_over_mp
print(f"\n2/α = {2/alpha_em:.4f}")
print(f"v/m_p = {v_over_mp:.4f}")
print(f"Difference = {diff_2a:.4f}")
print(f"Compare: d₁ + λ₁ = {d1 + lambda1}")
print(f"         2d₁ = {2*d1}")
print(f"         2λ₁ + 1 = {2*lambda1+1}")
print(f"         d₁ + λ₁ + 2/3 = {d1 + lambda1 + 2/3:.4f}")

# v/m_p = 2/α - (d₁ + λ₁ + something small)?
for corr_name, corr_val in [
    ("d₁+λ₁", d1+lambda1),
    ("2d₁", 2*d1),
    ("2(d₁+λ₁)/2", d1+lambda1),
    ("d₁+λ₁+K", d1+lambda1+K_koide),
    ("d₁+λ₁+η(2/9)", d1+lambda1+2/9),
    ("d₁+λ₁+1/π", d1+lambda1+1/np.pi),
    ("2d₁-η", 2*d1-2/9),
    ("d₁+λ₁+α·(something)", d1+lambda1+0.7),
]:
    pred = 2/alpha_em - corr_val
    pct = (pred/v_over_mp - 1)*100
    print(f"  2/α - ({corr_name}) = 2/α - {corr_val:.4f} = {pred:.4f}  ({pct:+.3f}%)")

# ===================== THE FOLD STIFFNESS =====================
print("\n" + "=" * 70)
print("PART 8: FOLD STIFFNESS / CURVATURE INTERPRETATION")
print("=" * 70)

# The fold wall is S⁴ with curvature radius R
# The stiffness = curvature energy density × volume
# On S⁴ of radius R: Ricci scalar = 12/R² (for S⁴: R_scalar = n(n-1)/R² = 4×3/R² = 12/R²)
# Einstein-Hilbert action ∝ ∫ R_scalar × √g = 12/R² × Vol(S⁴) × R⁴ = 12 × Vol(S⁴) × R²
# For unit S⁴: action = 12 × 8π²/3 = 32π²

stiffness_S4 = 12 * vol_S4  # = 32π²
print(f"Einstein-Hilbert stiffness of fold wall S⁴: 32π² = {stiffness_S4:.4f}")
print(f"v/m_p / (32π²) = {v_over_mp/stiffness_S4:.6f}")
print(f"That's ≈ {v_over_mp/stiffness_S4:.4f}")

# Number of fold walls: 3
print(f"\n3 fold walls × 32π²/3 = 32π² = {stiffness_S4:.4f}")
print(f"Per fold wall: 32π²/3 = {stiffness_S4/3:.4f}")

# What if v/m_p = stiffness / (some normalization)?
print(f"\nstiffness / (d₁λ₁) = {stiffness_S4/(d1*lambda1):.4f}")
print(f"stiffness / d₁² = {stiffness_S4/d1**2:.4f}")
# 32π²/36 ≈ 8.77 ... not useful

# Actually: v/m_p ≈ 262.4 and 32π² ≈ 315.8
# 315.8/1.2 ≈ 263.2!
print(f"\n32π²/1.2 = {stiffness_S4/1.2:.4f}")
print(f"v/m_p = {v_over_mp:.4f}")
print(f"  Not clean...")

# But: 8π² × 10/3 = 263.2!
val_fold = 8*np.pi**2 * 10 / 3
print(f"\n8π²×10/3 = Vol(S⁴) × 10 = {val_fold:.4f}")
print(f"v/m_p = {v_over_mp:.4f}")  
print(f"Error: {(val_fold/v_over_mp - 1)*100:+.3f}%")

# Wait: 10 = d₁ × λ₁ / 3. So this is Vol(S⁴) × d₁λ₁/3 = 8π²/3 × 10 = 80π²/3
# Or: Vol(S⁴) × d₁λ₁/3 
val_fold2 = vol_S4 * d1 * lambda1 / 3
print(f"\nVol(S⁴) × d₁λ₁/3 = (8π²/3)(30/3) = 80π²/3 = {val_fold2:.4f}")
print(f"v/m_p = {v_over_mp:.4f}")
print(f"Error: {(val_fold2/v_over_mp - 1)*100:+.3f}%")

# ===================== MONTE CARLO VERIFICATION =====================
print("\n" + "=" * 70)
print("PART 9: MONTE CARLO OVERLAP ON S⁵")
print("=" * 70)

np.random.seed(42)
N_MC = 5_000_000

# Sample S⁵
x = np.random.randn(N_MC, 6)
norms = np.sqrt(np.sum(x**2, axis=1, keepdims=True))
x = x / norms
z_mc = np.zeros((N_MC, 3), dtype=complex)
z_mc[:, 0] = x[:, 0] + 1j*x[:, 1]
z_mc[:, 1] = x[:, 2] + 1j*x[:, 3]
z_mc[:, 2] = x[:, 4] + 1j*x[:, 5]

# Verify moments
print(f"MC ⟨|z₁|²⟩ = {np.mean(np.abs(z_mc[:,0])**2):.6f} (exact: 1/3)")
print(f"MC ⟨|z₁|²|z₂|²⟩ = {np.mean(np.abs(z_mc[:,0])**2 * np.abs(z_mc[:,1])**2):.6f} (exact: 1/12)")
print(f"MC ⟨|z₁z₂z₃|²⟩ = {np.mean(np.prod(np.abs(z_mc)**2, axis=1)):.6f} (exact: 1/60)")

# The Z₃ "transfer matrix" element
# T_{jk} = ∫ z̄_j (Rz_k) dΩ / Vol = ω × δ_{jk} × 1/3
# This is just group representation theory

# More interesting: sector-resolved overlap
# Define sectors by the phase of z₁z₂z₃ (Z₃ invariant quantity has phase mod 2π)
# Actually better: define by a Z₃-breaking quantity

# The key observable: how much does a ghost mode "leak" across sector boundaries?
# Define fold wall as locus where two sectors have equal "claim"
# This is where the Z₃ phase is exactly ±π/3 from a wall

# Product phase (Z₃ invariant → use its cube root as sector label)
prod = z_mc[:,0] * z_mc[:,1] * z_mc[:,2]
prod_phase = np.angle(prod)  # This is Z₃-invariant phase

# Sum of squares (another Z₃ invariant)
sum_sq = z_mc[:,0]**2 + z_mc[:,1]**2 + z_mc[:,2]**2
sum_sq_phase = np.angle(sum_sq)

# For a truly Z₃-breaking sector label, use the cubic root
# z₁z₂z₃ has charge 3 → its cube root has charge 1 → sectors
# But cube root is ambiguous. Instead use arg(z₁) which shifts by 2π/3 under Z₃

sector_phase = np.angle(z_mc[:, 0])
sector_label = np.floor(((sector_phase % (2*np.pi)) / (2*np.pi/3))).astype(int) % 3

# Ghost mode overlap between adjacent sectors:
# Integral of |z₁|² over fold wall region (thin strip near sector boundary)
wall_width = 0.1  # radians
wall_mask = np.abs(sector_phase % (2*np.pi/3) - np.pi/3) < wall_width/2
frac_wall = np.sum(wall_mask) / N_MC

print(f"\nFraction of S⁵ near fold walls (width={wall_width:.2f}): {frac_wall:.4f}")
print(f"Expected (3 walls × width / 2π): {3*wall_width/(2*np.pi):.4f}")

# Ghost mode intensity at wall vs bulk
ghost_at_wall = np.mean(np.abs(z_mc[wall_mask, 0])**2)
ghost_at_bulk = np.mean(np.abs(z_mc[~wall_mask, 0])**2)
print(f"⟨|z₁|²⟩ at wall: {ghost_at_wall:.6f}")
print(f"⟨|z₁|²⟩ in bulk: {ghost_at_bulk:.6f}")
print(f"Enhancement at wall: {ghost_at_wall/ghost_at_bulk:.6f}")

# ===================== SUMMARY =====================
print("\n" + "=" * 70)
print("SUMMARY: BEST CANDIDATES FOR v/m_p")
print("=" * 70)

results = [
    ("2/α − 2d₁ = 2(1/α − 6)", 2/alpha_em - 2*d1),
    ("2/α − (d₁+λ₁)", 2/alpha_em - (d1+lambda1)),
    ("Vol(S⁴) × d₁λ₁/3", vol_S4 * d1*lambda1/3),
    ("C(d₁,2) × λ₁ × π", 15 * lambda1 * np.pi),
    ("d₁λ₁ + d₂λ₂ = 270", d1*lambda1 + d2*lambda2),
    ("32π²/1.2", stiffness_S4/1.2),
]

print(f"\n{'Formula':<45} {'Value':>10} {'v/m_p':>10} {'Error':>10}")
print("-" * 75)
for name, val in sorted(results, key=lambda x: abs(x[1]/v_over_mp - 1)):
    pct = (val/v_over_mp - 1)*100
    print(f"{name:<45} {val:>10.4f} {v_over_mp:>10.4f} {pct:>+9.4f}%")

# Final: what v/m_e would be for best candidates
print(f"\n--- Implied v/m_e for best candidates ---")
for name, val in sorted(results, key=lambda x: abs(x[1]/v_over_mp - 1)):
    vme = val * 6 * np.pi**5
    pct = (vme/v_over_me - 1)*100
    print(f"  {name}: v/m_e = {vme:.1f} ({pct:+.4f}%)")

def test_S5_moments():
    """Test the S5 moment calculations."""
    # Test basic moments
    assert abs(S5_moment(1, 0, 0) - 1/3) < 1e-10, "⟨|z|²⟩ should be 1/3"
    assert abs(S5_moment(2, 0, 0) - 1/6) < 1e-10, "⟨|z|⁴⟩ should be 1/6"
    assert abs(S5_moment(1, 1, 0) - 1/12) < 1e-10, "⟨|z₁|²|z₂|²⟩ should be 1/12"
    assert abs(S5_moment(1, 1, 1) - 1/60) < 1e-10, "⟨|z₁|²|z₂|²|z₃|²⟩ should be 1/60"
    print("All S5 moment tests passed!")

def main():
    """Main function to run the VEV overlap analysis."""
    test_S5_moments()
    # All the computation is done above in global scope
    # This allows the script to be imported as a module
    pass

if __name__ == "__main__":
    main()
