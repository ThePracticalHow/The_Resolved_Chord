"""
CKM Matrix from Z₃ Fold Wall Tunneling
========================================

Physical picture:
- Sphere (high E): all sectors equivalent, democratic mixing (1/3 each)
- Triangle (low E): fold walls rigid, mixing hierarchical
- CKM = tunneling amplitudes through fold walls
- Adjacent walls: leading bleed = η = 2/9 (Cabibbo)
- Two walls: (2/9)² 
- Three walls (full loop): (2/9)³
- CP phase: from complex structure of Z₃ (ω = e^{2πi/3})
"""

import numpy as np

# ==================== PDG CKM VALUES ====================
# Wolfenstein parameters (PDG 2024)
lam_W = 0.22500   # ± 0.00067  (= sin θ_C)
A_W = 0.826       # ± 0.015
rho_bar = 0.159   # ± 0.010
eta_bar = 0.348   # ± 0.010

# CKM magnitudes (PDG 2024)
V_ud = 0.97435; V_us = 0.22500; V_ub = 0.00369
V_cd = 0.22486; V_cs = 0.97349; V_cb = 0.04182
V_td = 0.00857; V_ts = 0.04110; V_tb = 0.999118

# CP phase
delta_CP = 1.144  # radians (≈ 65.6°)

print("=" * 70)
print("CKM MATRIX: PDG VALUES")
print("=" * 70)
print(f"  |V_ud| = {V_ud:.5f}   |V_us| = {V_us:.5f}   |V_ub| = {V_ub:.5f}")
print(f"  |V_cd| = {V_cd:.5f}   |V_cs| = {V_cs:.5f}   |V_cb| = {V_cb:.5f}")
print(f"  |V_td| = {V_td:.5f}   |V_ts| = {V_ts:.5f}   |V_tb| = {V_tb:.5f}")
print(f"\n  CP phase δ = {delta_CP:.4f} rad = {np.degrees(delta_CP):.2f}°")

# ==================== SPECTRAL DATA ====================
eta = 2/9
d1 = 6; l1 = 5; K = 2/3
omega = np.exp(2j * np.pi / 3)

print("\n" + "=" * 70)
print("HYPOTHESIS: CKM FROM POWERS OF η = 2/9")  
print("=" * 70)

# The bleed hypothesis:
# sin θ₁₂ = η = 2/9
# sin θ₂₃ = η² or some spectral modification
# sin θ₁₃ = η³ or some spectral modification

s12 = eta              # = 2/9
s23_naive = eta**2     # = 4/81
s13_naive = eta**3     # = 8/729

print(f"\n  sin θ₁₂ = η = 2/9 = {s12:.6f}")
print(f"  PDG |V_us| = {V_us:.6f}")
print(f"  Error: {(s12/V_us - 1)*100:+.3f}%")

print(f"\n  sin θ₂₃ = η² = 4/81 = {s23_naive:.6f}")
print(f"  PDG |V_cb| = {V_cb:.6f}")
print(f"  Error: {(s23_naive/V_cb - 1)*100:+.3f}%")
# 0.0494 vs 0.0418 — 18% off. Need correction.

print(f"\n  sin θ₁₃ = η³ = 8/729 = {s13_naive:.6f}")
print(f"  PDG |V_ub| = {V_ub:.6f}")
print(f"  Error: {(s13_naive/V_ub - 1)*100:+.3f}%")
# 0.01097 vs 0.00369 — way off for pure η³

# ==================== WOLFENSTEIN CONNECTION ====================
print("\n" + "=" * 70)
print("WOLFENSTEIN PARAMETERIZATION")
print("=" * 70)

# Standard Wolfenstein:
# V_us ≈ λ
# V_cb ≈ Aλ²  
# V_ub ≈ Aλ³(ρ - iη)
# |V_ub| ≈ Aλ³√(ρ² + η²)

print(f"\n  If λ = 2/9:")
A_pred = V_cb / (2/9)**2
print(f"  A = V_cb/λ² = {V_cb}/{(2/9)**2:.6f} = {A_pred:.4f}")
print(f"  PDG A = {A_W:.4f}  (error: {(A_pred/A_W - 1)*100:+.3f}%)")

rho_eta_mag = V_ub / (A_pred * (2/9)**3)
print(f"  √(ρ²+η²) = |V_ub|/(Aλ³) = {rho_eta_mag:.4f}")
print(f"  PDG √(ρ̄²+η̄²) = {np.sqrt(rho_bar**2 + eta_bar**2):.4f}")

# ==================== GEOMETRIC CORRECTIONS ====================
print("\n" + "=" * 70)
print("GEOMETRIC CORRECTIONS TO POWER LAW")
print("=" * 70)

# The pure power law η^n is too naive. Each "wall crossing" sees
# the fold wall geometry, which has its own spectral weight.

# What if A (the Wolfenstein A parameter) is also spectral?
# A ≈ 0.826. What spectral quantity is that?

print("\n--- What is A? ---")
candidates_A = {
    "K + 1/d₁": K + 1/d1,                    # = 0.833
    "5/6 = λ₁/d₁": l1/d1,                    # = 0.833
    "1 - 1/d₁": 1 - 1/d1,                    # = 0.833
    "1 - K/4": 1 - K/4,                       # = 0.833
    "d₁/(d₁+K+1/d₁)": d1/(d1+K+1/d1),        # something
    "5·K = 10/3": 5*K,                        # 3.333 no
    "K/η": K/eta,                              # 3 no
    "√K": np.sqrt(K),                          # 0.8165
    "η·(d₁-λ₁+3)": eta*(d1-l1+3),            # 0.889
    "η·d₁/K²": eta*d1/K**2,                   # 3 no
    "(λ₁-1)/λ₁": (l1-1)/l1,                   # 0.8
    "π/(2+K+1)": np.pi/(2+K+1),               # 0.857
}

print(f"  PDG A = {A_W:.4f}")
for name, val in sorted(candidates_A.items(), key=lambda x: abs(x[1]/A_W - 1)):
    if 0.5 < val < 1.5:
        err = (val/A_W - 1)*100
        marker = " <<<" if abs(err) < 2 else ""
        print(f"    {name:<25} = {val:.6f}  ({err:+.3f}%){marker}")

# λ₁/d₁ = 5/6 ≈ 0.833 is 0.9% from PDG A = 0.826
print(f"\n*** BEST: A = λ₁/d₁ = 5/6 = {l1/d1:.6f} ***")
print(f"*** PDG: A = {A_W:.4f} ({(l1/d1/A_W-1)*100:+.3f}%) ***")

# ==================== FULL CKM WITH A = 5/6 ====================
print("\n" + "=" * 70)
print("CKM WITH λ = 2/9, A = 5/6")
print("=" * 70)

lam_pred = 2/9
A_pred2 = 5/6

# Wolfenstein to CKM (to order λ⁴)
s12_pred = lam_pred
s23_pred = A_pred2 * lam_pred**2
s13_pred_mag = A_pred2 * lam_pred**3  # × √(ρ²+η²) for |V_ub|

print(f"\n  sin θ₁₂ = λ = 2/9 = {s12_pred:.6f}")
print(f"  |V_us| PDG = {V_us:.6f}  ({(s12_pred/V_us-1)*100:+.3f}%)")

print(f"\n  sin θ₂₃ = Aλ² = (5/6)(2/9)² = {s23_pred:.6f}")
print(f"  |V_cb| PDG = {V_cb:.6f}  ({(s23_pred/V_cb-1)*100:+.3f}%)")

print(f"\n  |V_ub| base = Aλ³ = (5/6)(2/9)³ = {s13_pred_mag:.6f}")
print(f"  |V_ub| PDG = {V_ub:.6f}")
print(f"  Ratio |V_ub|/(Aλ³) = {V_ub/s13_pred_mag:.4f} = √(ρ̄²+η̄²)")

# So ρ̄² + η̄² = (V_ub/(Aλ³))²
rho_eta_sq = (V_ub / s13_pred_mag)**2
print(f"  ρ̄² + η̄² = {rho_eta_sq:.4f}")
print(f"  PDG: ρ̄² + η̄² = {rho_bar**2 + eta_bar**2:.4f}")

# ==================== THE CP PHASE ====================
print("\n" + "=" * 70)
print("CP PHASE: WHERE DOES δ ≈ 65.6° COME FROM?")
print("=" * 70)

print(f"\n  δ_CP = {delta_CP:.4f} rad = {np.degrees(delta_CP):.2f}°")
print(f"\n  Geometric candidates:")

cp_candidates = {
    "2π/9": 2*np.pi/9,
    "π/3": np.pi/3,
    "2π/3 - π/d₁": 2*np.pi/3 - np.pi/d1,
    "arctan(η̄/ρ̄)": np.arctan2(eta_bar, rho_bar),
    "π/3 + η": np.pi/3 + eta,
    "arcsin(2/3)": np.arcsin(2/3),
    "arccos(K)": np.arccos(K),
    "arctan(2)": np.arctan(2),
    "arctan(K×π)": np.arctan(K*np.pi),
    "1 + η": 1 + eta,
    "2η×π": 2*eta*np.pi,
    "π²/9": np.pi**2/9,
    "d₁η×π": d1*eta*np.pi,
    "arcsin(K×√K)": np.arcsin(K*np.sqrt(K)),
}

for name, val in sorted(cp_candidates.items(), key=lambda x: abs(x[1] - delta_CP)):
    deg = np.degrees(val)
    err = (val/delta_CP - 1)*100
    marker = " <<<" if abs(err) < 3 else (" **" if abs(err) < 10 else "")
    print(f"    {name:<25} = {val:.4f} rad = {deg:.2f}°  ({err:+.3f}%){marker}")

# ==================== JARLSKOG INVARIANT ====================
print("\n" + "=" * 70)
print("JARLSKOG INVARIANT")
print("=" * 70)

# J = Im(V_us V_cb V*_ub V*_cd) ≈ A²λ⁶η̄ (Wolfenstein)
# J ≈ 3.08 × 10⁻⁵ (PDG)
J_pdg = 3.08e-5

# With our parameters:
J_pred_base = A_pred2**2 * lam_pred**6
print(f"\n  J = A²λ⁶ × η̄ (needs η̄)")
print(f"  A²λ⁶ = (5/6)²(2/9)⁶ = {J_pred_base:.6e}")
print(f"  J_PDG = {J_pdg:.6e}")
print(f"  Needed η̄ = J_PDG/(A²λ⁶) = {J_pdg/J_pred_base:.4f}")
print(f"  PDG η̄ = {eta_bar:.4f}")

# What spectral quantity gives η̄ ≈ 0.348?
print(f"\n--- What is η̄ ≈ {eta_bar}? ---")
for name, val in [
    ("1/3", 1/3),
    ("η×π/2", eta*np.pi/2),
    ("K/2 + 1/90", K/2 + 1/90),
    ("1/π + 1/30", 1/np.pi + 1/30),
    ("sin(2π/9)", np.sin(2*np.pi/9)),
    ("2/(d₁-K/3)", 2/(d1-K/3)),
    ("η×√(π)", eta*np.sqrt(np.pi)),
    ("sin(δ_koide - π/2)", np.sin(2*np.pi/3 + 2/9 - np.pi/2)),
    ("7/(4π²)", 7/(4*np.pi**2)),
    ("λ₁/(4π+K)", l1/(4*np.pi + K)),
]:
    err = (val/eta_bar - 1)*100
    marker = " <<<" if abs(err) < 3 else (" **" if abs(err) < 10 else "")
    print(f"    {name:<30} = {val:.6f}  ({err:+.3f}%){marker}")

# ==================== FULL CKM RECONSTRUCTION ====================
print("\n" + "=" * 70)
print("SUMMARY: CKM FROM SPECTRAL DATA")
print("=" * 70)

print(f"""
  Wolfenstein parameters from S⁵/Z₃:
    λ = sin θ_C = 2/9 = η_Donnelly          [{(s12_pred/V_us-1)*100:+.3f}%]
    A = λ₁/d₁ = 5/6                         [{((5/6)/A_W-1)*100:+.3f}%]
    
  Derived CKM elements:
    |V_us| = 2/9                             = {2/9:.5f} (PDG: {V_us:.5f})
    |V_cb| = (5/6)(2/9)²  = 10/162 = 5/81   = {5/81:.5f} (PDG: {V_cb:.5f})
    |V_ub| = (5/6)(2/9)³ × √(ρ̄²+η̄²)       (needs ρ̄, η̄)
    
  What's derived: 2 of 4 Wolfenstein params (λ, A)
  What's needed: ρ̄, η̄ (or equivalently, the CP phase δ)
  
  Key insight: V_cb = λ₁/(d₁ × d₁²/K) ... no
  Simpler: V_cb = (λ₁/d₁)(η²) = (5/6)(4/81) = 20/486 = 10/243
  = {10/243:.6f}
  PDG: {V_cb:.6f}  ({(10/243/V_cb-1)*100:+.3f}%)
""")

# Clean ratio check
print("CLEAN RATIOS:")
print(f"  V_us = 2/9 = {2/9:.6f}")
print(f"  V_cb = 10/243 = {10/243:.6f}  (PDG: {V_cb:.6f}, {((10/243)/V_cb-1)*100:+.3f}%)")
print(f"  V_cb × 9/2 = {V_cb * 9/2:.6f} (= V_cb/V_us)")
print(f"  V_cb/V_us² = {V_cb/V_us**2:.6f} (= A)")

# The V_td, V_ts predictions
V_td_pred = A_pred2 * lam_pred**3  # ≈ Aλ³ (leading order)
V_ts_pred = A_pred2 * lam_pred**2  # ≈ V_cb (leading order, -Aλ²)
print(f"\n  |V_td| ≈ Aλ³ = {A_pred2 * lam_pred**3:.6f} (PDG: {V_td:.5f})")
print(f"  |V_ts| ≈ Aλ² = {A_pred2 * lam_pred**2:.6f} (PDG: {V_ts:.5f})")
