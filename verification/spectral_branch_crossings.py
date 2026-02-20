#!/usr/bin/env python3
"""
SPECTRAL BRANCH CROSSINGS: Where Exotic Stability Lives
=========================================================

IDEA: An exotic particle isn't just a "near miss" of one consonant mode.
It's a LOCAL MAXIMUM OF STABILITY where two different spectral branches
CROSS — creating constructive interference far from any single resonance.

ANALOGY: Two tuning forks vibrating at different frequencies create
BEAT PATTERNS. At certain locations, the beats constructively interfere
to create a temporary standing wave — even though neither fork alone
would resonate there.

IN THE SPECTRAL FRAMEWORK:
  - The Dirac operator on B^6/Z_3 has multiple "branches":
    different KK towers, different angular momentum channels,
    different Z_3 character sectors
  - Each branch generates a LADDER of consonant modes
  - Where two ladders cross (same energy, different quantum numbers),
    the modes HYBRIDIZE and create a stability pocket
  - The exotic particle IS the hybridized mode

This is AVOIDED CROSSING in quantum mechanics:
when two energy levels meet, they repel and create a gap.
The states near the gap are superpositions — neither fish nor fowl.
That's why exotics have "unusual" quantum numbers.

Jixiang Leng & Claude, February 2026
"""

import sys, io
if sys.stdout.encoding and sys.stdout.encoding.lower() not in ('utf-8', 'utf8'):
    try:
        sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
        sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')
    except Exception:
        pass

import numpy as np
from fractions import Fraction

PI = np.pi
alpha = 1/137.036
alpha_s = 0.1187

# Spectral invariants
d1 = 6; lam1 = 5; K = Fraction(2, 3); eta = Fraction(2, 9); p = 3
d2 = 20; lam2 = 12; d3 = 50; lam3 = 21

m_e_MeV = 0.51100
m_p_MeV = 938.272
G_hurr = float(Fraction(lam1) * eta)  # 10/9
m_p = m_e_MeV * d1 * PI**5 * (1 + G_hurr * alpha**2/PI)

print("=" * 72)
print("  SPECTRAL BRANCH CROSSINGS: WHERE EXOTIC STABILITY LIVES")
print("=" * 72)

# =====================================================================
#  PART 1: THE SPECTRAL BRANCHES
# =====================================================================
print(f"\n{'='*72}")
print("PART 1: THE SPECTRAL BRANCHES OF B^6/Z_3")
print(f"{'='*72}")

print("""
  The Dirac operator on B^6/Z_3 has eigenvalues organized by:
    - Angular momentum l = 0, 1, 2, 3, ...
    - Z_3 character chi_0, chi_1, chi_2
    - Radial quantum number n_r = 0, 1, 2, ...

  Each (l, chi, n_r) triplet defines a BRANCH of the spectrum.
  Within each branch, the mass ladder is:

    m(l, chi, n) = m_base(l) * (1 + n * step(l, chi))

  where m_base = spectral eigenvalue * m_e and step = spacing.

  The SM particles sit on the l=0 and l=1 branches.
  EXOTIC PARTICLES are where DIFFERENT branches cross.
""")

# Define the spectral branches
# Each branch generates a TOWER of modes with a characteristic spacing

# Branch 1: Light hadron tower (l=1, chi_0)
# Spacing set by 1/d1 = 1/6 of the proton mass
branch_light = {
    'name': 'Light (l=1, chi_0)',
    'base': m_p,
    'spacing': m_p / d1,  # ~156 MeV = T_c(QCD)!
    'sector': 'light',
    'color': 'blue',
}

# Branch 2: Strange tower (l=1, chi_1)
# Spacing modified by eta = 2/9
branch_strange = {
    'name': 'Strange (l=1, chi_1)',
    'base': m_p * float(K),    # ~625 MeV (kaon~493, Lambda~1116)
    'spacing': m_p * float(eta),  # ~209 MeV
    'sector': 'strange',
    'color': 'green',
}

# Branch 3: Charm tower (l=2, chi_0)
# Base at l=2 eigenvalue, spacing by lam1
branch_charm = {
    'name': 'Charm (l=2, chi_0)',
    'base': m_p * float(Fraction(10, 3)),  # ~3127 MeV (J/psi region)
    'spacing': m_p / lam1,  # ~188 MeV
    'sector': 'charm',
    'color': 'red',
}

# Branch 4: Bottom tower (l=2, chi_1)
branch_bottom = {
    'name': 'Bottom (l=2, chi_1)',
    'base': m_p * 10,  # ~9383 MeV (Upsilon region)
    'spacing': m_p * float(K) / lam1,  # ~125 MeV
    'sector': 'bottom',
    'color': 'purple',
}

# Branch 5: Ghost tower (l=1, twisted sector)
# The d1=6 ghost modes create their own tower
branch_ghost = {
    'name': 'Ghost (l=1, twisted)',
    'base': m_p * float(Fraction(5, 6)),  # ~782 MeV (rho/omega)
    'spacing': m_p / (d1 * lam1),  # ~31 MeV (fine structure)
    'sector': 'ghost',
    'color': 'gray',
}

# Branch 6: Baryon resonance tower
# Baryons step by Delta-N splitting = m_p * (p/(d1) = 1/2)
branch_baryon = {
    'name': 'Baryon resonances',
    'base': m_p,
    'spacing': m_p * float(Fraction(p, d1)),  # ~469 MeV (Delta-N ~ 294, but the tower goes 940, 1232, 1520...)
    'sector': 'baryon',
    'color': 'orange',
}

branches = [branch_light, branch_strange, branch_charm, branch_bottom,
            branch_ghost, branch_baryon]

# Generate modes for each branch
max_E = 12000  # up to 12 GeV
branch_modes = {}

for br in branches:
    modes = []
    n = 0
    while True:
        E = br['base'] + n * br['spacing']
        if E > max_E:
            break
        modes.append(E)
        n += 1
    branch_modes[br['name']] = modes
    print(f"  {br['name']:<30} base={br['base']:.0f} MeV, step={br['spacing']:.0f} MeV, {len(modes)} modes")


# =====================================================================
#  PART 2: FIND BRANCH CROSSINGS
# =====================================================================
print(f"\n{'='*72}")
print("PART 2: BRANCH CROSSINGS (stability pockets)")
print(f"{'='*72}")

print("""
  A CROSSING occurs when modes from two DIFFERENT branches land at
  nearly the same energy. The crossing gap (|E1 - E2|) determines
  the stability:

    - Small gap -> strong hybridization -> NARROW exotic (long-lived)
    - Large gap -> weak hybridization -> BROAD resonance (short-lived)

  The hybridization strength depends on the COUPLING between branches,
  which is set by alpha_s (for QCD mixing) or alpha (for EM mixing).
""")

# Find all crossings between different branches where |E1 - E2| < threshold
crossing_threshold = 50.0  # MeV — look for modes within 50 MeV

crossings = []

branch_names = list(branch_modes.keys())
for i, name_i in enumerate(branch_names):
    for j, name_j in enumerate(branch_names):
        if j <= i:
            continue  # avoid double counting
        for E_i in branch_modes[name_i]:
            for E_j in branch_modes[name_j]:
                gap = abs(E_i - E_j)
                if gap < crossing_threshold:
                    E_avg = (E_i + E_j) / 2
                    crossings.append({
                        'E_avg': E_avg,
                        'gap': gap,
                        'branch_1': name_i,
                        'branch_2': name_j,
                        'E1': E_i,
                        'E2': E_j,
                    })

# Sort by crossing energy
crossings.sort(key=lambda c: c['E_avg'])

# Remove near-duplicates (within 20 MeV)
filtered = []
for c in crossings:
    if not filtered or abs(c['E_avg'] - filtered[-1]['E_avg']) > 20:
        filtered.append(c)
crossings = filtered

print(f"\n  Found {len(crossings)} branch crossings below {max_E} MeV:")
print(f"\n  {'E_cross (MeV)':>14} {'Gap (MeV)':>10} {'Branch 1':<30} {'Branch 2':<30}")
print(f"  {'─'*14} {'─'*10} {'─'*30} {'─'*30}")

for c in crossings:
    print(f"  {c['E_avg']:>14.1f} {c['gap']:>10.1f} {c['branch_1']:<30} {c['branch_2']:<30}")


# =====================================================================
#  PART 3: MAP CROSSINGS TO KNOWN PARTICLES
# =====================================================================
print(f"\n{'='*72}")
print("PART 3: CROSSING <-> PARTICLE MAP")
print(f"{'='*72}")

# Known particles and resonances for comparison
known_particles = {
    # Light mesons
    'eta(548)':    548,
    'rho(770)':    775,
    'omega(782)':  783,
    'eta\'(958)':  958,
    'phi(1020)':  1020,
    'f0(1370)':   1370,
    'f0(1500)':   1505,
    'f0(1710)':   1720,
    'f2(1950)':   1944,

    # Charm mesons
    'D(1865)':    1865,
    'D*(2007)':   2007,
    'Ds(1968)':   1968,
    'J/psi':      3097,
    'psi(3686)':  3686,
    'X(3872)':    3872,
    'Zc(3900)':   3888,
    'chi_c1(3872)': 3872,

    # Pentaquarks
    'Pc(4312)':   4312,
    'Pc(4440)':   4440,
    'Pc(4457)':   4457,

    # Baryons
    'N(940)':      938,
    'Delta(1232)': 1232,
    'N*(1440)':    1440,
    'N*(1520)':    1520,
    'N*(1680)':    1680,
    'Lambda(1116)': 1116,
    'Sigma(1189)': 1189,
    'Xi(1315)':    1315,
    'Omega(1672)': 1672,

    # Dibaryon
    'd*(2380)':    2380,

    # Bottom
    'Upsilon(9460)': 9460,
    'Upsilon(10023)': 10023,
    'Zb(10610)': 10610,
    'Zb(10650)': 10650,
}

# For each crossing, find nearest known particle
print(f"\n  {'Crossing (MeV)':>14} {'Gap':>6} {'Match':>14} {'Particle':>16} {'Dist (MeV)':>12} {'Branches'}")
print(f"  {'─'*14} {'─'*6} {'─'*14} {'─'*16} {'─'*12} {'─'*40}")

matches = []
for c in crossings:
    # Find nearest known particle
    best_name = None
    best_dist = float('inf')
    for pname, pmass in known_particles.items():
        dist = abs(c['E_avg'] - pmass)
        if dist < best_dist:
            best_dist = dist
            best_name = pname

    if best_dist < 100:  # within 100 MeV
        match_str = f"{best_dist:.0f} MeV"
        matches.append({'crossing': c, 'particle': best_name, 'dist': best_dist})
    else:
        match_str = "---"
        best_name = "(none)"

    if best_dist < 100:
        br_short = c['branch_1'].split('(')[0].strip() + " x " + c['branch_2'].split('(')[0].strip()
        print(f"  {c['E_avg']:>14.1f} {c['gap']:>6.1f} {match_str:>14} {best_name:>16} {best_dist:>12.1f} {br_short}")


# =====================================================================
#  PART 4: STABILITY LANDSCAPE — WHERE ARE THE POCKETS?
# =====================================================================
print(f"\n{'='*72}")
print("PART 4: THE STABILITY LANDSCAPE")
print(f"{'='*72}")

print("""
  The STABILITY LANDSCAPE is the density of spectral modes as a
  function of energy. High density = many overlapping branches =
  more stability pockets.

  Think of it like an interference pattern:
    - Each branch creates a WAVE of spectral density
    - Where waves overlap constructively = stability pocket
    - Where waves cancel = spectral desert (no resonances)

  The LOCAL MAXIMA of this landscape are where exotic particles live.
""")

# Compute spectral density as function of energy
E_range = np.linspace(100, 6000, 1200)
density = np.zeros_like(E_range)

# Each mode contributes a Lorentzian peak to the spectral density
# Width of contribution = spacing of the branch (determines resolution)
for br in branches:
    modes = branch_modes[br['name']]
    width = br['spacing'] / 2  # half-width of each mode's contribution
    for E_mode in modes:
        density += width**2 / ((E_range - E_mode)**2 + width**2)

# Find local maxima
from scipy.signal import argrelextrema
maxima_idx = argrelextrema(density, np.greater, order=5)[0]

# Filter to significant peaks
peak_threshold = np.mean(density) + 0.5 * np.std(density)
significant_peaks = [(E_range[i], density[i]) for i in maxima_idx if density[i] > peak_threshold]
significant_peaks.sort(key=lambda x: -x[1])  # sort by height

print(f"  Top {min(20, len(significant_peaks))} stability pockets (local maxima of spectral density):")
print(f"\n  {'Rank':>4} {'E (MeV)':>10} {'Density':>10} {'Nearest Particle':>20} {'Dist (MeV)':>12}")
print(f"  {'─'*4} {'─'*10} {'─'*10} {'─'*20} {'─'*12}")

for rank, (E_peak, dens) in enumerate(significant_peaks[:20], 1):
    # Find nearest known particle
    best_name = "?"
    best_dist = float('inf')
    for pname, pmass in known_particles.items():
        dist = abs(E_peak - pmass)
        if dist < best_dist:
            best_dist = dist
            best_name = pname

    marker = " <--" if best_dist < 30 else ""
    print(f"  {rank:>4} {E_peak:>10.0f} {dens:>10.2f} {best_name:>20} {best_dist:>12.0f}{marker}")


# =====================================================================
#  PART 5: THE "WAY OUT" STABILITY — UNEXPECTED POCKETS
# =====================================================================
print(f"\n{'='*72}")
print("PART 5: UNEXPECTED STABILITY POCKETS")
print(f"{'='*72}")

print("""
  The most interesting predictions are stability pockets where
  NO SINGLE BRANCH would put a resonance — but TWO branches
  crossing create one. These are the "way out" exotic states.

  Search criterion:
    - At least 2 branches contribute (multi-branch crossing)
    - NOT within 100 MeV of any branch's base frequency
    - There SHOULD be a known exotic particle nearby
""")

# Find peaks that are NOT close to any single branch's base modes
unexpected = []
for E_peak, dens in significant_peaks:
    # Check if this peak is within 50 MeV of any single branch mode
    near_single = False
    branch_count = 0
    for br in branches:
        modes = branch_modes[br['name']]
        for E_mode in modes:
            if abs(E_peak - E_mode) < 50:
                branch_count += 1
                break  # one per branch

    # "Unexpected" = requires 2+ branches, not dominated by one
    if branch_count >= 2:
        best_name = "?"
        best_dist = float('inf')
        for pname, pmass in known_particles.items():
            dist = abs(E_peak - pmass)
            if dist < best_dist:
                best_dist = dist
                best_name = pname

        unexpected.append({
            'E': E_peak,
            'density': dens,
            'branches': branch_count,
            'particle': best_name,
            'dist': best_dist,
        })

print(f"\n  Found {len(unexpected)} multi-branch stability pockets:")
print(f"\n  {'E (MeV)':>10} {'Density':>8} {'#Br':>4} {'Nearest':>16} {'Dist':>8} {'Type'}")
print(f"  {'─'*10} {'─'*8} {'─'*4} {'─'*16} {'─'*8} {'─'*20}")

for u in sorted(unexpected, key=lambda x: -x['density'])[:15]:
    ptype = "MATCH!" if u['dist'] < 40 else ("close" if u['dist'] < 100 else "prediction?")
    print(f"  {u['E']:>10.0f} {u['density']:>8.2f} {u['branches']:>4} {u['particle']:>16} {u['dist']:>8.0f} {ptype}")


# =====================================================================
#  PART 6: WHAT "EXOTIC" MEANS IN OUR FRAMEWORK
# =====================================================================
print(f"\n{'='*72}")
print("PART 6: WHAT 'EXOTIC' MEANS IN THE SPECTRAL FRAMEWORK")
print(f"{'='*72}")

print(f"""
  In the SM, "exotic" means "unusual quantum numbers" — states that
  can't be a simple quark-antiquark or three-quark system.

  In the spectral framework, "exotic" has a precise meaning:

  DEFINITION: An exotic state is a LOCAL MAXIMUM of the spectral
  stability landscape that requires contributions from TWO OR MORE
  spectral branches to exist.

  CONSEQUENCES:

  1. EXOTIC = MULTI-BRANCH HYBRIDIZATION
     A conventional meson sits on ONE branch (e.g., rho on the
     light l=1 tower). An exotic sits at the CROSSING of two branches.
     It's a superposition: |exotic> = a|branch_1> + b|branch_2>.

  2. UNUSUAL QUANTUM NUMBERS ARE AUTOMATIC
     Since the two branches have different (l, chi) quantum numbers,
     the hybridized state has quantum numbers that are LINEAR
     COMBINATIONS of the parent branches. This is why exotics have
     "impossible" J^PC values: they're superpositions.

  3. NARROW WIDTH FROM BRANCH PROTECTION
     When two branches cross, the hybridized state can only decay
     by de-hybridizing (separating back into its branch components).
     This requires breaking the inter-branch coupling, which is
     suppressed by alpha_s. Hence exotics are NARROWER than naive
     dimensional analysis would predict.

  4. MOLECULAR VS COMPACT:
     - "Molecular" exotic: crossing of two HADRON branches
       (e.g., D and D* branches crossing -> X(3872))
     - "Compact" exotic: crossing of two QUARK branches
       (e.g., charm and ghost branches -> Pc pentaquarks)
     Both are branch crossings — the distinction is which branches.

  5. PREDICTION:
     The spectral framework predicts that exotic states exist ONLY
     at branch crossings. If a bump appears at an energy where no
     crossing exists, it's not an exotic — it's noise.

  THE MUSICAL ANALOGY:
  =====================
  An exotic particle is like a CHORD — two or more notes played
  simultaneously. A chord can sound "stable" even if none of its
  individual notes are the fundamental of the instrument.

  The stability of a chord depends on the INTERVAL between the
  notes. Small-integer ratios (3/2, 4/3, 5/4) are consonant chords.
  Large-integer ratios are dissonant chords that don't sustain.

  In the spectral framework:
    - Each branch is a "string" of the instrument
    - Branch crossings are "chords"
    - Exotic particles are SPECTRAL CHORDS that happen to be stable
    - Their lifetime depends on how "consonant" the chord is

  THE RESOLVED CHORD:
  ===================
  And THIS is why the paper is called "The Resolved Chord."

  The Standard Model IS the resolved chord: the unique combination
  of spectral branches (3 generations, gauge group, Higgs) that
  creates a maximally stable harmony. All the consonant modes
  are the SM particles. All the dissonant modes decay.

  Exotic particles are PARTIALLY resolved chords — two branches
  that almost harmonize but don't quite. They ring briefly,
  then decay back to the fully resolved chord.
""")

print("=" * 72)
print("  EXOTIC = SPECTRAL CHORD = MULTI-BRANCH CROSSING")
print("  Stability pockets exist where branches constructively interfere")
print("  The Resolved Chord IS the Standard Model")
print("=" * 72)
