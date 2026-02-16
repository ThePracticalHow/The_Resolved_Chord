"""
Higgs Quartic λ: Quick Hunt
============================
m_H² = 2λv²  →  λ = m_H²/(2v²)
"""
import numpy as np

# PDG values
m_H = 125.25        # GeV (±0.17)
v = 246.2196         # GeV
m_p = 0.93827208943  # GeV
m_e = 0.51099895e-3  # GeV
alpha = 1/137.035999177

# The quartic
lam = m_H**2 / (2 * v**2)

# Key ratios
mH_over_v = m_H / v
mH_over_mp = m_H / m_p
mH_over_me = m_H / m_e

print("=" * 60)
print("HIGGS QUARTIC: THE NUMBERS")
print("=" * 60)
print(f"  m_H = {m_H} GeV")
print(f"  v   = {v} GeV")
print(f"  λ   = m_H²/(2v²) = {lam:.6f}")
print(f"  m_H/v   = {mH_over_v:.6f}")
print(f"  m_H/m_p = {mH_over_mp:.4f}")
print(f"  m_H/m_e = {mH_over_me:.1f}")
print(f"  (m_H/v)² = 2λ = {mH_over_v**2:.6f}")

# Spectral data
d1 = 6; l1 = 5; K = 2/3
d2 = 20; l2 = 12
eta = 2/9
pi2 = np.pi**2

print("\n" + "=" * 60)
print("RATIO HUNT")
print("=" * 60)

# m_H/v ≈ 0.5087 — close to 1/2!
print(f"\n--- m_H/v ≈ {mH_over_v:.6f} ---")
print(f"  1/2 = 0.500000  (error: {(mH_over_v/0.5-1)*100:+.3f}%)")
print(f"  1/2 + η/4 = {0.5 + eta/4:.6f}")
print(f"  1/2 + 1/18 = {0.5 + 1/18:.6f}  ({(mH_over_v/(0.5+1/18)-1)*100:+.3f}%)")
print(f"  λ₁/(π²) = {l1/pi2:.6f}")
print(f"  1/2 + α = {0.5 + alpha:.6f}")
print(f"  1/√(2π) = {1/np.sqrt(2*np.pi):.6f}")

# m_H/m_p ≈ 133.49 — close to 1/α - 3.5!
print(f"\n--- m_H/m_p ≈ {mH_over_mp:.4f} ---")
print(f"  1/α = {1/alpha:.4f}  (diff: {1/alpha - mH_over_mp:.4f})")
print(f"  1/α - d₁/2 = {1/alpha - 3:.4f}  ({(mH_over_mp/(1/alpha-3)-1)*100:+.3f}%)")
print(f"  1/α - λ₁+1 = {1/alpha - l1 + 1:.4f}")
print(f"  1/α - K×d₁ = {1/alpha - K*d1:.4f}  ({(mH_over_mp/(1/alpha-4)-1)*100:+.3f}%)")
print(f"  1/α - 35/9 = {1/alpha - 35/9:.4f}  ({(mH_over_mp/(1/alpha-35/9)-1)*100:+.4f}%)")

# λ ≈ 0.1293
print(f"\n--- λ ≈ {lam:.6f} ---")
print(f"  1/8 = {1/8:.6f}  ({(lam/0.125-1)*100:+.3f}%)")
print(f"  K/λ₁ = {K/l1:.6f}  ({(lam/(K/l1)-1)*100:+.3f}%)")
print(f"  2/(3λ₁) = {2/(3*l1):.6f}  ({(lam/(2/15)-1)*100:+.3f}%)")
print(f"  1/(d₁+K+1/d₁) = {1/(d1+K+1/d1):.6f}")
print(f"  η × λ₁/(d₁+K) = {eta*l1/(d1+K):.6f}")
print(f"  π²/(6×d₁²) = {pi2/(6*36):.6f}")
print(f"  (d₁-λ₁)²/(d₁²-1) = {1/35:.6f}")

# Actually: 2λ = (m_H/v)² ≈ 0.2588
two_lam = mH_over_v**2
print(f"\n--- 2λ = (m_H/v)² ≈ {two_lam:.6f} ---")
print(f"  1/4 = 0.250000  ({(two_lam/0.25-1)*100:+.3f}%)")
print(f"  K/K+K = {K/(K+K):.6f}")  # 1/2... no
print(f"  (π²-5)/(2π²) = {(pi2-5)/(2*pi2):.6f}  ({(two_lam/((pi2-5)/(2*pi2))-1)*100:+.3f}%)")
print(f"  λ₁/(d₁+λ₁+K+d₁²) = {l1/(d1+l1+K+d1**2):.6f}")

# The VEV formula gives v/m_p = 2/α - 35/3
# If m_H/m_p = 1/α - X, then m_H/v = (1/α - X)/(2/α - 35/3)
# For m_H/v ≈ 1/2, we need X ≈ 35/6 ≈ 5.833
print(f"\n--- IF m_H/v = (1/α - X)/(2/α - 35/3), solve for X ---")
X_needed = 1/alpha - mH_over_v * (2/alpha - 35/3)
print(f"  X = {X_needed:.6f}")
print(f"  35/6 = {35/6:.6f}  ({(X_needed/(35/6)-1)*100:+.4f}%)")
print(f"  d₁ - 1/6 = {d1-1/6:.6f}")
print(f"  d₁ - K/4 = {d1 - K/4:.6f}")
print(f"  (d₁+λ₁+K)/2 = {(d1+l1+K)/2:.6f}  ({(X_needed/((d1+l1+K)/2)-1)*100:+.4f}%)")

# WAIT: (d₁+λ₁+K)/2 = 35/6!
print(f"\n*** (d₁+λ₁+K)/2 = 35/6 = {35/6:.6f} ***")
print(f"*** X_needed = {X_needed:.6f} ***")
print(f"*** Match: {(X_needed/(35/6)-1)*100:+.4f}% ***")

# So: m_H/m_p = 1/α - (d₁+λ₁+K)/2 = 1/α - 35/6
mH_mp_pred = 1/alpha - 35/6
print(f"\n  m_H/m_p predicted = 1/α - 35/6 = {mH_mp_pred:.6f}")
print(f"  m_H/m_p measured  = {mH_over_mp:.6f}")
print(f"  Error: {(mH_mp_pred/mH_over_mp - 1)*100:+.4f}%")

# And: m_H/v = (1/α - 35/6) / (2/α - 35/3)
mH_v_pred = (1/alpha - 35/6) / (2/alpha - 35/3)
print(f"\n  m_H/v predicted = (1/α - 35/6)/(2/α - 35/3) = {mH_v_pred:.6f}")
print(f"  m_H/v measured  = {mH_over_v:.6f}")
print(f"  Error: {(mH_v_pred/mH_over_v - 1)*100:+.4f}%")

# Simplify: (1/α - 35/6)/(2/α - 35/3) = (1/α - 35/6)/(2(1/α - 35/6)) = 1/2 !!!!!
print(f"\n" + "!" * 60)
print("WAIT.")
print(f"  Numerator:   1/α - 35/6")
print(f"  Denominator: 2/α - 35/3 = 2(1/α - 35/6)")
print(f"  RATIO = 1/2 EXACTLY.")
print(f"  m_H/v = 1/2 EXACTLY if the same spectral cost appears")
print(f"  in both m_H and v with the SAME relative weight.")
print("!" * 60)

# So m_H = v/2 exactly in this framework?
print(f"\n  m_H = v/2 would give {v/2:.4f} GeV")
print(f"  m_H measured = {m_H:.4f} GeV")
print(f"  Error: {(v/2/m_H - 1)*100:+.3f}%")

# λ = m_H²/(2v²) = (v/2)²/(2v²) = 1/8
print(f"\n  λ = 1/8 = {1/8:.6f}")
print(f"  λ measured = {lam:.6f}")
print(f"  Error: {(0.125/lam - 1)*100:+.3f}%")

# But wait — m_H is measured as 125.25, v/2 = 123.11
# That's 1.7% off. Not great.
# The EXACT prediction is m_H/m_p = 1/α - 35/6, not m_H = v/2
# The 1/2 only holds if BOTH use the same spectral cost 35/3

# Let's check: does m_H have its OWN spectral cost?
print("\n" + "=" * 60)
print("REFINING: m_H/m_p with independent spectral cost")
print("=" * 60)

# v/m_p = 2/α - 35/3  (factor 2, cost 35/3)
# m_H/m_p = 1/α - C_H  (factor 1, cost C_H)
# Find C_H:
C_H = 1/alpha - mH_over_mp
print(f"  C_H = 1/α - m_H/m_p = {C_H:.6f}")
print(f"  35/6 = {35/6:.6f}  (half the VEV cost)")
print(f"  Diff from 35/6: {C_H - 35/6:.6f}")
print(f"  That diff ≈ {C_H - 35/6:.4f}")

# C_H ≈ 3.547. Compare to spectral things:
print(f"\n  Looking for C_H ≈ {C_H:.6f}:")
for name, val in [
    ("35/6", 35/6),
    ("d₁/2 + K/2", d1/2 + K/3),  
    ("d₁-λ₁/2", d1 - l1/2),
    ("K×λ₁", K*l1),
    ("d₁×K - 1/2", d1*K - 0.5),
    ("(d₁+λ₁)/π", (d1+l1)/np.pi),
    ("11/π", 11/np.pi),
    ("7/2", 7/2),
    ("(d₁²-1)/10", (d1**2-1)/10),
    ("K×d₁-η", K*d1 - eta),
    ("4-η×2", 4 - 2*eta),
    ("λ₁-√2", l1-np.sqrt(2)),
    ("π+η", np.pi + eta),
]:
    err = (val/C_H - 1)*100
    marker = " <<<" if abs(err) < 1 else (" **" if abs(err) < 5 else "")
    print(f"    {name:<20} = {val:.6f}  ({err:+.3f}%){marker}")

# Try: m_H/m_p = 1/α - K×d₁ + η  
# = 137.036 - 4 + 0.222 = 133.258? 
print(f"\n  1/α - (K×d₁ - η) = {1/alpha - (K*d1-eta):.4f} vs {mH_over_mp:.4f}")
print(f"  1/α - 7/2 = {1/alpha - 3.5:.4f} vs {mH_over_mp:.4f}  ({((1/alpha-3.5)/mH_over_mp-1)*100:+.3f}%)")
print(f"  1/α - (d₁²-1)/10 = {1/alpha - 3.5:.4f}")

# Quartic from various m_H/v predictions
print("\n" + "=" * 60)
print("QUARTIC λ FROM VARIOUS m_H PREDICTIONS")
print("=" * 60)
for name, mH_pred in [
    ("m_H = v/2 (exact half)", v/2),
    ("m_H/m_p = 1/α - 35/6", (1/alpha - 35/6) * m_p),
    ("m_H/m_p = 1/α - 7/2", (1/alpha - 3.5) * m_p),
    ("m_H/m_p = 1/α - (d₁²-1)/10", (1/alpha - 35/10) * m_p),
]:
    lam_pred = mH_pred**2 / (2*v**2)
    err = (mH_pred/m_H - 1)*100
    print(f"  {name:<35}: m_H = {mH_pred:.3f} GeV ({err:+.3f}%), λ = {lam_pred:.6f}")

print(f"\n  Measured: m_H = {m_H:.3f} GeV, λ = {lam:.6f}")
