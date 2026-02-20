# Supplement_XII_CompanionGuide

Supplement XII: The Companion Guide
Topic-Based Reference for the Spectral Geometry of S5/Z3
The Resolved Chord — Supplementary Material
Six chapters. One geometry. Zero free parameters.
Jixiang Leng
February 2026
This document provides six topic-based views into the framework. The main paper (v10) tells
the story; the Supplements (I–XII) provide the proofs; this guide collects everything a reader
might want to look up in one place. Each chapter is self-contained and cross-references the
relevant supplements.
Contents
1
The Hurricane Database
2
1.1
The hurricane principle . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
2
1.2
Complete hurricane table . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
2
1.3
Spectral decomposition of each coefficient . . . . . . . . . . . . . . . . . . . .
2
2
The Equation Index
4
2.1
Mass formulas . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
4
2.2
Coupling formulas . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
4
2.3
Mixing formulas . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
4
2.4
Cosmological formulas
. . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
5
3
The LOTUS Field Guide
6
3.1
The fold field at each value of ϕ . . . . . . . . . . . . . . . . . . . . . . . . .
6
3.2
The LOTUS potential
. . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
6
3.3
What the LOTUS generates . . . . . . . . . . . . . . . . . . . . . . . . . . .
6
4
All Predictions: The Falsification Battery
7
4.1
Near-term testable predictions . . . . . . . . . . . . . . . . . . . . . . . . . .
7
4.2
Structural anti-predictions . . . . . . . . . . . . . . . . . . . . . . . . . . . .
7
4.3
Scorecard
. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
7
5
Derivation Status: The Firewall
8
1


6
Strange Castles: Solved Puzzles
9
2


1
The Hurricane Database
Every bare spectral formula receives radiative corrections (“hurricanes”) with coefficients
that are simple ratios of spectral invariants. This chapter collects all of them.
1.1
The hurricane principle
Bare formula: exact spectral expression (tree-level).
Hurricane correction: observed = bare × (1 + c · αX/π), where c is a spectral coefficient
and αX is the relevant coupling (α for EM corrections, αs for QCD corrections).
1.2
Complete hurricane table
Observable
Bare
Coeff. c
Coupling
Corrected
Error
mp/me
6π5
G
=
10/9
α2/π
1836.153
10−11
(2-loop)
G2
=
−280/9
α4/π2
Cabibbo λ
η = 2/9
+1/p =
+1/3
αs/(3π)
0.22502
0.009%
Wolfenstein A
λ1/d1 = 5/6
−η
=
−2/9
αs/π
0.8263
0.04%
1/αGUT
42.41 (RG)
G/p
=
10/27
topological
42.78
0.001%
1/α3 split
1/αGUT
−d1
=
−6
spectral
36.78
0.56%
Gravity X
121/3
−1/(d1λ1) = −1/30
3509/90
0.10%
Table 1: Complete hurricane coefficient table.
Every
coefficient is a ratio of {d1, λ1, K, η, p}.
1.3
Spectral decomposition of each coefficient
• G = λ1 · η = 5 × 2/9 = 10/9
(eigenvalue × spectral asymmetry)
• G2 = −λ1(d1 + η) = −5(6 + 2/9) = −280/9
(fermion trace structure)
• 1/p = 1/3
(orbifold sector count)
• −η = −2/9
(spectral asymmetry, opposite sign)
• G/p = λ1η/p = 10/27
(APS lag: asymmetry per sector)
• −d1 = −6
(ghost mode count, SU(3) splitting)
3


• −1/(d1λ1) = −1/30
(inverse ghost spectral weight)
Key observation 1: All coefficients are O(1) or smaller in spectral units. If the bare
formulas were wrong by order-one factors, the corrections would be O(π/α) ∼400, not O(1).
The fact that all |c| ≲1 is strong evidence the bare geometry is correct.
Key observation 2:
The Hurricane Pattern (meta-theorem).
Every hurricane
coefficient follows one of four rules depending on the geometric locus of the correction:
Locus
Rule
Coefficients
Origin
Within one sector
÷ p
+1/p, G/p
Orbifold volume factor
Between sectors
× η
−η, ηλ1/p
APS spectral asymmetry
Ghost mode counting
× d1
−d1, d1λ1
Mode count (trace)
Eigenvalue weighting
× λ1
G = λ1η, λD
1 = 7/2
Kinetic energy scale
Inverse ghost weight
÷ d1λ1
−1/(d1λ1)
Total ghost spectral weight
This pattern is not imposed — it follows from the spectral action loop expansion. One-loop
corrections are traces over K = S5/Z3, and traces of spectral operators are spectral invariants.
The pattern IS the loop structure of the spectral action.
4


2
The Equation Index
Every prediction in the framework, its formula, spectral ingredients, status, and verification
script.
2.1
Mass formulas
Mass
Formula
Value
Error
Status
mp/me
6π5(1 + Gα2/π)
1836.153
10−11
Thm
mµ/me
Koide(K=2/3,
δ=2π/3+2/9)
206.768
0.0001%
Thm
mτ/me
Koide(K=2/3,
δ=2π/3+2/9)
3477.4
0.01%
Thm
mt
(v/
√
2) e−1/120
172.66 GeV
0.02%
Thm
mc
(v/
√
2)(mµ/mτ) e−2π/3
1.275 GeV
0.15%
Thm
mu
(v/
√
2)(me/mτ) e−π
2.16 MeV
0.17%
Thm
mb
mτ e77/90
4.180 GeV
0.06%
Thm
ms
mµ e−10/81
93.4 MeV
0.01%
Thm
md
me e2π/3+G/p2
4.69 MeV
0.53%
Thm
v
mp(2/α −35/3)
246.2 GeV
0.004%
Thm
mH
mp(1/α −7/2)
125.3 GeV
0.036%
Thm
MP
Mc · (3509/90)7/2 ·
p
π3/3
1.22 × 1019
0.10%
Thm
mν3
m3
e/(p m2
p)
∼50 meV
Theorem
Thm
2.2
Coupling formulas
Coupling
Formula
Value
Error
Status
1/α
1/αGUT + ηλ1/p + RG
137.038
0.001%
Thm
sin2 θW
3/8 at Mc
0.375
Thm
Thm
αs(MZ)
1/(1/αGUT,corr −d1 + RG)
0.1187
0.56%
Thm
λH
(mH/mp)2/[2(v/mp)2]
0.1295
0.14%
Thm
¯θ
0 (geometric CP)
0
exact
Thm
2.3
Mixing formulas
Parameter
Formula
Value
Error
Status
CKM λ
η(1 + αs/3π)
0.2250
0.009%
Thm
CKM A
(λ1/d1)(1 −ηαs/π)
0.826
0.04%
Thm
CKM ¯ρ
1/(2π)
0.1592
0.03%
Thm
CKM ¯η
π/9 = ηD · π/2
0.3491
0.02%
Thm
CKM γ
arctan(2π2/9)
65.49◦
0.17%
Thm
5


Jarlskog J
A2λ6¯η
3.09 × 10−5
0.5%
Thm
PMNS θ23
arcsin
p
d1/(d1+λ1)
∼47◦
∼1%
Thm
PMNS θ12
PSF framework
∼33◦
∼2%
Thm
PMNS θ13
PSF framework
∼8.5◦
∼2%
Thm
2.4
Cosmological formulas
Quantity
Formula
Value
Error
Status
Λ1/4
mν3 · η2(1 −K/d1)
2.22 meV
1.4%
Thm
N (e-folds)
(d1+λ1)2a2/(p a4)
=
3025/48
≈63
0.8σ
Thm
ns
1 −2/N
0.968
0.3%
Thm
r
12/N 2
0.003
< 0.036
Thm
ηB
α4 · η
6.3 × 10−10
3%
Thm
ΩDM/ΩB
d1 −K = 16/3
5.333
0.5%
Thm
6


3
The LOTUS Field Guide
LOTUS = Lagrangian Of The Universe’s Spectral State.
The fold field ϕ parameterizes the transition from the smooth parent sphere S5 (ϕ = 0) to
the rigid orbifold S5/Z3 (ϕ = 1).
3.1
The fold field at each value of ϕ
ϕ
Name
Physics
0
Smooth sphere
No Z3 structure. No generations. No
masses. Featureless.
0.60
Phase transition
Crossover from substrate-dominated
to information-dominated.
Ghost
modes decouple. Inflation ends.
0.9574
The lotus in bloom
Our universe. v = vmaxϕlotus. SM
physics. DM = frozen ghosts. CC =
lotus breathing (mν3η2).
1
Dried flower
Fully rigid fold. v = 0. All masses
vanish. Gauge unification restored.
Unreachable (ghost pressure).
3.2
The LOTUS potential
V (ϕ) = λH
4 v4
max (ϕ2 −ϕ2
lotus)2
where vmax = 2mp/α and ϕlotus = 1 −α(d1+λ1+K)/2 = 0.9574.
This IS the Mexican hat potential in fold-depth coordinates: H = vmaxϕ, v = vmaxϕlotus.
3.3
What the LOTUS generates
• V (ϕlotus) = 0
(tree-level CC vanishes)
• V ′′(ϕlotus) →mH = 125.3 GeV
(Higgs mass = curvature at equilibrium)
• V (0)1/4 ≈104 GeV
(barrier height = EW scale)
• Hurricanes = V ′(ϕ) near ϕlotus
(perturbative corrections)
• Black holes: ϕ > ϕlotus locally
(lotus closing)
• Inflation: ϕ rolling from 0 to ϕlotus
(dimensional unfolding)
Verification: lotus potential.py.
7


4
All Predictions: The Falsification Battery
4.1
Near-term testable predictions
Prediction
Value
Experiment
Kills theory if
mτ
1776.985 MeV
Belle
II
(±0.05)
> 0.5 MeV devi-
ation
λH
0.1295 (no BSM)
HL-LHC
BSM correction
found
αs(MZ)
0.1187
Lattice QCD
> 3σ from spec-
tral
P mν
≈59.2 meV
DESI / Euclid
> 80 or < 40
meV
r (tensor/scalar)
0.003
LiteBIRD
/
CMB-S4
r > 0.03
ns
0.968
Planck
/
CMB-S4
> 3σ deviation
4.2
Structural anti-predictions
Anti-prediction
Experiment
Kills theory if
No QCD axion
ADMX / IAXO
Axion detected
No 4th generation
Colliders
4th lepton or quark found
Normal hierarchy
JUNO
Inverted hierarchy con-
firmed
No proton decay
Hyper-K
Proton decay observed
No free quarks
Any
Isolated quark observed
DM null (direct de-
tection)
XENON / LZ
WIMP-like DM detected
Qν ̸= 2/3
Precision ν
Qν = 2/3 measured
4.3
Scorecard
Category
Count
Theorem
48
Derived
0
Identified
0
Gaps
0
Total
48
New prediction: MW = 79.90 GeV (from sin2 θW = 3/8 + RG).
8


5
Derivation Status: The Firewall
For every claim, the strongest skeptical objection and our response. This is the “if a reviewer
says X, the answer is Y” document. Full details: Supplement XI.
Claim
St.
Skeptic says
Response
K = 2/3
Thm
“Why this map?”
Unique moment map
on S5 with Z3
Ng = 3
Thm
“Is APS right?”
Direct spectral decom-
position; n = 3 unique
mp/me = 6π5
Thm
“Why 6π5?”
Parseval fold energy;
50-digit proof
1/α = 137.038
Thm
“Circular?”
APS lag ηλ1/p; all fac-
tors Thm
v/mp = 2/α −
35/3
Thm
“EM budget?”
α Thm;
ghost cost
35/3 Thm
mH/mp = 1/α −
7/2
Thm
“Dirac eigenvalue?”
Ikeda 1980; α Thm
X = 3509/90
Thm
“Coincidence?”
5-lock; p < 10−5
¯ρ = 1/(2π)
Thm
“Numerology?”
Fourier norm of S1;
0.03%
¯η = π/9
Thm
“J calculation?”
Donnelly η rotated by
torsion arg
αs = 0.1187
Thm
“Why d1 = 6?”
Ghost modes are 3⊕¯3
CC = 2.22 meV
Thm
“Heavy modes?”
Equidistribution veri-
fied to l = 500
9


6
Strange Castles: Solved Puzzles
Seven major physics puzzles resolved by the spectral geometry of S5/Z3.
Full details:
Supplement IX.
SC-θ: Strong CP without an axion. ¯θ = 0 from Z3-circulant CP symmetry (Theorem).
Anti-prediction: no axion, no θ-tuning.
SC-grav: Why gravity is weak. MP/Mc ∼106 because d1λ1 = 30 ghost modes dilute
the bulk coupling. Hierarchy = ghost spectral weight, not fine-tuning.
SC-mix: Why quarks and leptons mix differently. Charged fermions: twisted sector
(cone point) →small CKM, exact Koide. Neutrinos: untwisted sector (fold walls) →
large PMNS, no Koide.
SC-CP: Why CP is violated. ¯η/¯ρ = 2π2/9 is irrational (Lindemann–Weierstrass). CP
violation = geometric incommensurability of cone and circle.
SC-hurricane: Why all residuals are small. Every correction coefficient is O(1) in spec-
tral units. Bare formulas are correct at the geometric level; QCD/EM corrections are
perturbative.
SC-Λ: Why the CC is small but nonzero. Tree: V (ϕlotus) = 0 (orbifold cancellation).
One-loop: η2 = 4/81 suppression from double boundary crossing.
SC-33: The spectral integer 33. 33 = d2
1 −p: neutrino splitting, X17 anomaly, fused
quark Koide, tunneling bandwidth. Four appearances of one integer from one geometry.
Each castle is a problem that has resisted decades of theoretical effort. The spectral framework
addresses all seven with the same five invariants. No new fields. No new symmetries. No
new parameters. Just one manifold.
10
