"""
Anti-numerology Monte Carlo test for the Lotus Song hadron spectrum.
Tests whether random sets of small integers can match 23 hadron masses
at the same precision as the spectral invariants {d1=6, lam1=5, K=2/3, eta=2/9, p=3}.
"""
import random
from collections import Counter

# PDG hadron mass ratios (m/m_p)
HADRONS = {
    'pi': 0.1488, 'K': 0.5262, 'eta': 0.5839, 'eta_prime': 1.0208,
    'rho': 0.8263, 'K*': 0.9503, 'phi': 1.0865,
    'Delta': 1.3131, 'Sigma*': 1.4757, 'Xi*': 1.6326, 'Omega': 1.7825,
    'Lambda': 1.1891, 'Xi': 1.4047, 'Sigma': 1.2711,
    'D': 1.9927, 'D*': 2.1389, 'Ds': 2.0978,
    'B': 5.6267, 'Bs': 5.7200, 'Bc': 6.6877,
    'Jpsi': 3.3006, 'Upsilon': 10.0827, 'Upsilon2S': 10.6827,
}
TARGETS = list(HADRONS.values())
N_HADRONS = len(TARGETS)

def generate_pool(invariants):
    pool = set()
    vals = list(set(invariants))
    for v in vals:
        if 0.1 < v < 12:
            pool.add(round(v, 10))
    for a in vals:
        for b in vals:
            for op in [a+b, a-b, a*b]:
                if 0.1 < op < 12:
                    pool.add(round(op, 10))
            if b != 0:
                r = a/b
                if 0.1 < r < 12:
                    pool.add(round(r, 10))
    for a in vals:
        for b in vals:
            for c in vals:
                if c != 0:
                    for op in [a*b/c, (a+b)/c, a*b+c, a+b/c, a-b/c, a*b*c]:
                        if 0.1 < op < 12:
                            pool.add(round(op, 10))
                if b+c != 0:
                    r = a/(b+c)
                    if 0.1 < r < 12:
                        pool.add(round(r, 10))
                if b*c != 0:
                    r = a/(b*c)
                    if 0.1 < r < 12:
                        pool.add(round(r, 10))
    return pool

def count_matches(pool, targets, threshold=0.02):
    matched = 0
    for t in targets:
        for v in pool:
            if abs(v - t) / t < threshold:
                matched += 1
                break
    return matched

# Real invariants
real_inv = [6, 5, 2/3, 2/9, 3, 7, 11]
real_pool = generate_pool(real_inv)
real_matches = count_matches(real_pool, TARGETS)

print(f"REAL INVARIANTS: d1=6, lam1=5, K=2/3, eta=2/9, p=3, Dw=7, Db=11")
print(f"Pool size: {len(real_pool)}")
print(f"Hadrons matched (within 2%): {real_matches}/{N_HADRONS}")
print()

# Monte Carlo: random invariant sets
random.seed(42)
N_TRIALS = 100000
match_counts = []
beat_real = 0

for trial in range(N_TRIALS):
    inv = [random.randint(2, 12) for _ in range(5)]
    inv.append(1 + inv[0])
    inv.append(inv[0] + inv[1])
    # Add rational invariants like K and eta
    if inv[2] != 0:
        inv.append(inv[0] / inv[2])
    if inv[4] != 0:
        inv.append(inv[1] / inv[4])
    
    pool = generate_pool(inv)
    matches = count_matches(pool, TARGETS)
    match_counts.append(matches)
    if matches >= real_matches:
        beat_real += 1

print(f"MONTE CARLO: {N_TRIALS} trials")
print(f"Each: 5 random ints [2-12] + 2 derived + 2 rational combos")
avg = sum(match_counts) / len(match_counts)
mx = max(match_counts)
print(f"Average matched: {avg:.1f}/{N_HADRONS}")
print(f"Maximum matched: {mx}/{N_HADRONS}")
print(f"Trials >= {real_matches}: {beat_real}/{N_TRIALS}")
if beat_real > 0:
    print(f"p-value: {beat_real/N_TRIALS:.6f}")
else:
    print(f"p-value: < {1/N_TRIALS:.1e}")
print()

dist = Counter(match_counts)
print("Distribution:")
for k in sorted(dist.keys()):
    pct = dist[k] / N_TRIALS * 100
    bar = '#' * int(pct * 2)
    print(f"  {k:2d}: {dist[k]:6d} ({pct:5.2f}%) {bar}")
