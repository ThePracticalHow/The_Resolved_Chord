#!/usr/bin/env python3
"""
PIERCING DEPTH UNIQUENESS TEST
================================

For each quark, exhaustively search all admissible expressions
from the constraint grammar {rationals from spectral data, pi/3, ln 3}
and count how many match PDG within 1%.

If the paper's expressions are the ONLY matches (or the simplest),
the uniqueness claim holds. If there are many alternatives,
the claim needs refinement.
"""

import numpy as np

PI = np.pi
LN3 = np.log(3)

# UV masses
m_e = 0.51099895e-3; m_mu = 105.6584e-3; m_tau = 1776.86e-3
v = 246.22; alpha = 1/137.036

m_t_UV = v / np.sqrt(2)
m_c_UV = v / np.sqrt(2) * m_mu / m_tau
m_u_UV = v / np.sqrt(2) * m_e / m_tau
m_b_UV = m_tau
m_s_UV = m_mu
m_d_UV = m_e

pdg = {'t': 172.57, 'c': 1.273, 'u': 0.00216, 'b': 4.183, 's': 0.0934, 'd': 0.00467}
uv = {'t': m_t_UV, 'c': m_c_UV, 'u': m_u_UV, 'b': m_b_UV, 's': m_s_UV, 'd': m_d_UV}

sigma_paper = {
    't': -1/120,
    'c': -2*PI/3,
    'u': -PI,
    'b': 77/90,
    's': -10/81,
    'd': 4*PI/3 - 2*LN3 + 2/9,
}

print("=" * 70)
print("  PIERCING DEPTH UNIQUENESS TEST")
print("=" * 70)

# Verify paper values
print("\n  Paper's sigma values:")
for q in ['t', 'c', 'u', 'b', 's', 'd']:
    sigma = sigma_paper[q]
    m_pred = uv[q] * np.exp(sigma)
    err = abs(m_pred - pdg[q]) / pdg[q] * 100
    print(f"    {q}: sigma={sigma:>10.6f}  m_pred={m_pred:.4e}  PDG={pdg[q]:.4e}  err={err:.3f}%")

# Constraint grammar:
# sigma = a + b*(pi/3) + c*ln(3)
# a = rational from spectral data denominators
# b, c = small integers

# Denominators from spectral invariants and their products
denoms = [1, 2, 3, 4, 5, 6, 9, 10, 12, 15, 18, 27, 30, 36, 45, 54, 81, 90, 120]
rationals = set()
for d in denoms:
    for n in range(-15, 16):
        rationals.add(n / d)
rationals = sorted(rationals)

b_range = range(-4, 5)  # coefficient of pi/3
c_range = range(-4, 5)  # coefficient of ln(3)

total_candidates = len(rationals) * len(list(b_range)) * len(list(c_range))
print(f"\n  Search space: {len(rationals)} rationals x {len(list(b_range))} pi/3 x {len(list(c_range))} ln3")
print(f"  = {total_candidates} candidate expressions per quark")

# For each quark: find all matches within 1%, sorted by complexity
print(f"\n{'='*70}")

for q in ['t', 'c', 'u', 'b', 's', 'd']:
    target_sigma = np.log(pdg[q] / uv[q])
    matches = []
    
    for a in rationals:
        for b in b_range:
            for c in c_range:
                sigma_test = a + b * (PI/3) + c * LN3
                m_test = uv[q] * np.exp(sigma_test)
                err = abs(m_test - pdg[q]) / pdg[q] * 100
                if err < 1.0:
                    # Complexity = total coefficient magnitude
                    # Lower = simpler = better
                    denom_a = 1
                    for d in denoms:
                        if abs(a * d - round(a * d)) < 1e-9:
                            denom_a = d
                            break
                    numer_a = round(a * denom_a)
                    complexity = abs(numer_a) + abs(b) * 5 + abs(c) * 5 + (denom_a > 1) * 3
                    
                    matches.append((complexity, err, a, b, c, sigma_test))
    
    matches.sort()  # by complexity, then error
    
    # Check if paper's value is in the list
    paper_rank = -1
    for i, (cpx, err, a, b, c, sig) in enumerate(matches):
        if abs(sig - sigma_paper[q]) < 0.01:
            paper_rank = i + 1
            break
    
    print(f"\n  {q}-quark: {len(matches)} matches within 1% (paper's is rank #{paper_rank})")
    print(f"  Target sigma = {target_sigma:.6f}")
    print(f"  {'Rank':>4} {'Expression':<40} {'Error':>8} {'Complexity':>10}")
    print(f"  {'-'*65}")
    
    for i, (cpx, err, a, b, c, sig) in enumerate(matches[:8]):
        # Format expression
        parts = []
        if a != 0:
            # Find simplest fraction
            for d in denoms:
                n = round(a * d)
                if abs(a - n/d) < 1e-9:
                    if d == 1:
                        parts.append(f"{n}")
                    else:
                        parts.append(f"{n}/{d}")
                    break
            else:
                parts.append(f"{a:.4f}")
        if b != 0:
            if b == 1: parts.append("pi/3")
            elif b == -1: parts.append("-pi/3")
            else: parts.append(f"{b}*pi/3")
        if c != 0:
            if c == 1: parts.append("ln3")
            elif c == -1: parts.append("-ln3")
            else: parts.append(f"{c}*ln3")
        
        expr = " + ".join(parts) if parts else "0"
        expr = expr.replace("+ -", "- ")
        
        is_paper = " <<<PAPER" if abs(sig - sigma_paper[q]) < 0.01 else ""
        print(f"  {i+1:>4} {expr:<40} {err:>7.4f}% {cpx:>10.0f}{is_paper}")

print(f"\n{'='*70}")
print(f"  SUMMARY")
print(f"{'='*70}")

total_matches = 0
paper_ranks = {}
for q in ['t', 'c', 'u', 'b', 's', 'd']:
    target_sigma = np.log(pdg[q] / uv[q])
    count = 0
    rank = 0
    
    for a in rationals:
        for b in b_range:
            for c in c_range:
                sigma_test = a + b * (PI/3) + c * LN3
                m_test = uv[q] * np.exp(sigma_test)
                err = abs(m_test - pdg[q]) / pdg[q] * 100
                if err < 1.0:
                    count += 1
    
    total_matches += count
    print(f"  {q}: {count} matches within 1%")

print(f"\n  Average: {total_matches/6:.0f} matches per quark within 1%")
print(f"  from {total_candidates} candidates = {total_matches/6/total_candidates*100:.2f}% hit rate")
print(f"\n  Probability of ALL 6 matching by chance:")
hit_rate = total_matches / 6 / total_candidates
p_all_6 = hit_rate ** 6
print(f"  ({total_matches/6/total_candidates:.4f})^6 = {p_all_6:.2e}")
