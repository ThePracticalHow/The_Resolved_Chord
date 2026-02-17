# Supplement_VII_NeutrinoSector

Supplement VII: The Neutrino Sector —
Parameters 21–26
Complete Derivation Chain for the Neutrino Mixing and Mass Sector
The Resolved Chord — Supplementary Material
Jixiang Leng
February 2026
Abstract
This supplement is self-contained. It provides the complete derivation chain
for Parameters 21–26 of the main text: the three PMNS mixing angles (reactor,
solar, atmospheric), the leptonic CP phase, the heaviest neutrino mass, and the
mass-squared splitting ratio. The key structural distinction — twisted versus
untwisted sectors of the Z3 orbifold — is developed first, followed by individual
parameter derivations, the unified point/side/face forward model, and derived
predictions for cosmological observables. All definitions, intermediate calculations,
and numerical verifications are included.
Contents
1
Twisted versus Untwisted Sectors
3
1.1
Twisted sector: charged fermions
. . . . . . . . . . . . . . . . . . . . .
3
1.2
Untwisted sector: neutrinos
. . . . . . . . . . . . . . . . . . . . . . . .
3
1.3
Summary comparison . . . . . . . . . . . . . . . . . . . . . . . . . . . .
4
2
Tribimaximal Mixing as Zeroth Order
4
2.1
Z3 permutation eigenvectors . . . . . . . . . . . . . . . . . . . . . . . .
4
3
Parameter 21 — Reactor Angle θ13
5
3.1
The formula . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
5
3.2
Numerical comparison
. . . . . . . . . . . . . . . . . . . . . . . . . . .
6
3.3
Geometric derivation . . . . . . . . . . . . . . . . . . . . . . . . . . . .
6
4
Parameter 22 — Solar Angle θ12
6
4.1
The formula . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
6
1


4.2
Numerical comparison
. . . . . . . . . . . . . . . . . . . . . . . . . . .
7
4.3
Geometric derivation . . . . . . . . . . . . . . . . . . . . . . . . . . . .
7
5
Parameter 23 — Atmospheric Angle θ23
7
5.1
The formula . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
7
5.2
Numerical comparison
. . . . . . . . . . . . . . . . . . . . . . . . . . .
7
5.3
Geometric derivation . . . . . . . . . . . . . . . . . . . . . . . . . . . .
8
6
Parameter 24 — Leptonic CP Phase δPMNS
CP
8
6.1
The formula . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
8
6.2
Numerical comparison
. . . . . . . . . . . . . . . . . . . . . . . . . . .
9
6.3
Geometric derivation . . . . . . . . . . . . . . . . . . . . . . . . . . . .
9
7
The Point/Side/Face Forward Derivation
9
7.1
The three geometric objects . . . . . . . . . . . . . . . . . . . . . . . .
10
7.2
Tunnelling overlap matrix
. . . . . . . . . . . . . . . . . . . . . . . . .
10
7.3
Mass matrix construction
. . . . . . . . . . . . . . . . . . . . . . . . .
11
7.4
Diagonalisation and PMNS extraction
. . . . . . . . . . . . . . . . . .
11
8
Parameter 25 — Neutrino Mass from the Inversion Principle
12
8.1
Bulk resonance versus boundary tunnelling . . . . . . . . . . . . . . . .
12
8.2
The master formula . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
12
8.3
Numerical verification
. . . . . . . . . . . . . . . . . . . . . . . . . . .
12
8.4
Connection to the seesaw mechanism . . . . . . . . . . . . . . . . . . .
13
9
Parameter 26 — Mass-Squared Splitting Ratio
13
9.1
The formula . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
13
9.2
Numerical comparison
. . . . . . . . . . . . . . . . . . . . . . . . . . .
13
9.3
Geometric derivation . . . . . . . . . . . . . . . . . . . . . . . . . . . .
14
9.4
Derived masses from P25 + P26 . . . . . . . . . . . . . . . . . . . . . .
14
10 Why No Neutrino Koide Ratio
15
10.1 Charged leptons: circulant symmetry . . . . . . . . . . . . . . . . . . .
15
10.2 Neutrinos: three different objects . . . . . . . . . . . . . . . . . . . . .
15
10.3 The neutrino Qν
. . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
15
11 Derived Predictions
16
11.1 Sum of neutrino masses . . . . . . . . . . . . . . . . . . . . . . . . . . .
16
11.2 Effective Majorana mass . . . . . . . . . . . . . . . . . . . . . . . . . .
16
12 Provenance Table
17
12.1 Connection to the lotus potential . . . . . . . . . . . . . . . . . . . . .
17
2


1
Twisted versus Untwisted Sectors
The Z3 orbifold S5/Z3 divides the five-sphere into three 120◦sectors. Two qualitatively
distinct classes of geometric structure arise: fold walls (the smooth codimension-1
interfaces separating adjacent sectors) and the cone point (the singular fixed-point locus
carrying a deficit angle of 2π/3). These structures govern the physics of the two fermion
sectors.
1.1
Twisted sector: charged fermions
Definition 1 (Twisted-sector string). A string on S5/Z3 is twisted if it is closed only
modulo the Z3 action: its endpoint returns to its starting point only after application of
the generator g.
Twisted-sector strings are pinned to the cone point, the singular locus of the orbifold.
They carry a twist phase ωk (k = 1, 2; ω = e2πi/3) encoding the Z3 charge. The key
consequences are:
(i) Mass origin. Masses of twisted-sector fermions arise from the topological struc-
ture of the conical defect. The deficit angle 2π/3 sets the harmonic content of the
wavefunctions localised at the cone tip.
(ii) Koide ratio. All three charged-lepton masses are eigenvalues of the same circulant
mass matrix (Supplement II, §1), determined by the single topological invariant
of the cone point. The Koide ratio
K =
me + mµ + mτ
(√me + √mµ + √mτ)2 = 2
3
(exact)
(1)
is a consequence of the circulant structure: it measures the harmonic content of a
fixed topological feature (the cone point).
(iii) CKM mixing.
Inter-generation mixing among quarks (also twisted-sector
fermions) requires tunnelling through the singular cone point. This tunnelling
amplitude is exponentially suppressed by the conical geometry, producing the
observed smallness of off-diagonal CKM elements.
1.2
Untwisted sector: neutrinos
Definition 2 (Untwisted-sector string). A string on S5/Z3 is untwisted if it is closed
without application of the Z3 generator: it satisfies standard periodicity on the quotient
space.
Untwisted-sector strings are not pinned to the cone point. They propagate freely
through the bulk of each 120◦sector and interact with the fold walls — the smooth,
extended, codimension-1 interfaces between adjacent sectors. The consequences are
sharply distinct from the twisted sector:
3


(i) Mass origin. Neutrino masses arise from tunnelling overlaps between wavefunc-
tions localised at different fold walls. Because the fold walls are smooth extended
surfaces (not singular points), these overlaps are geometric rather than topological.
(ii) Mixing angles. Tunnelling amplitudes between fold walls are of order one —
there is no exponential suppression from a singular barrier. This produces the
large mixing angles observed in the PMNS matrix.
(iii) Koide ratio. The three neutrino mass eigenstates are not three copies of a single
object with different Z3 charges; they are three geometrically distinct objects (see
Section 7). The circulant symmetry that enforces K = 2/3 for charged leptons
does not apply, and the neutrino Koide-like ratio Qν ̸= 2/3 (see Section 10).
1.3
Summary comparison
Property
Charged fermions (twisted)
Neutrinos (untwisted)
Pinned to
Cone point (singular)
Fold walls (smooth)
Mass origin
Topological (deficit angle)
Geometric (tunnelling overlap)
Koide ratio
K = 2/3 (exact)
Qν ̸= 2/3
Mixing angles
Small (exponential suppression)
Large (order-1 tunnelling)
Table 1: Twisted versus untwisted sectors of the Z3 orbifold.
Remark 1. The twisted/untwisted distinction is not a model-building choice; it is
forced by the topology of S5/Z3. Electrically charged fermions carry non-trivial Z3
representations and are therefore twisted; neutrinos are neutral under the Z3 and are
therefore untwisted. The qualitative differences in Table 1 — small versus large mixing,
exact versus approximate Koide — follow as geometric consequences.
2
Tribimaximal Mixing as Zeroth Order
The zeroth-order prediction for the PMNS matrix is tribimaximal mixing (TBM), which
arises as a mathematical identity from the Z3 symmetry.
2.1
Z3 permutation eigenvectors
The Z3 cyclic permutation operator on three objects is
P =


0
0
1
1
0
0
0
1
0

,
P 3 = I.
(2)
Its eigenvalues are {1, ω, ω2} with ω = e2πi/3.
4


Proposition 1 (TBM from Z3). The eigenvectors of P are precisely the columns of
the tribimaximal mixing matrix UTBM.
Proof. The normalised eigenvectors of P are:
v1 = 1
√
3(1, 1, 1)T
(eigenvalue 1),
(3)
v2 = 1
√
3(1, ω, ω2)T
(eigenvalue ω),
(4)
v3 = 1
√
3(1, ω2, ω)T
(eigenvalue ω2).
(5)
The democratic matrix D =
1
311T (all entries 1/3) commutes with P and shares
the same eigenvector structure. Restricting to real combinations and imposing the
conventional phase choices, the resulting unitary matrix is
UTBM =


p
2/3
1/
√
3
0
−1/
√
6
1/
√
3
1/
√
2
1/
√
6
1/
√
3
−1/
√
2

,
(6)
which yields the TBM mixing angles:
sin2 θ12 = 1
3,
sin2 θ23 = 1
2,
θ13 = 0.
(7)
Remark 2. Tribimaximal mixing is the undeformed Z3 prediction. Each measured
deviation from the TBM values (7) maps to a specific geometric property of the physical
orbifold:
• θ13 ̸= 0: junction asymmetry (Section 3),
• sin2 θ12 ̸= 1/3: finite fold-wall thickness (Section 4),
• sin2 θ23 ̸= 1/2: spectral impedance mismatch (Section 5).
3
Parameter 21 — Reactor Angle θ13
3.1
The formula
Theorem 1 (Reactor angle).
sin2 θ13 = (η K)2 =
2
9 · 2
3
2
=
 4
27
2
=
16
729.
(8)
5


3.2
Numerical comparison
sin2 θ13

pred = 16
729 = 0.02194.
(9)
The measured value (NuFIT 5.3, normal ordering [2]):
sin2 θ13

meas = 0.02200 ± 0.00069.
(10)
Deviation:
|0.02194 −0.02200|
0.02200
= 0.27%.
(11)
3.3
Geometric derivation
In the TBM limit (Section 2), exact Z3 symmetry gives θ13 = 0. The physical orbifold
S5/Z3 breaks this through two effects:
(i) Fold-wall bleed (η = 2/9). The fold walls separating the three 120◦sectors have
finite thickness. The Donnelly invariant η = 2/9 (Supplement I, §2; originally
[3]) measures the spectral asymmetry induced by this finite thickness. It is the
leading-order correction to the perfect Z3 geometry.
(ii) Harmonic lock (K = 2/3). The Koide constant K = 2/3 is the moment-map
constraint from the orbifold harmonic decomposition (Supplement V, §3.3). It
enters because the reactor angle couples the twisted-sector harmonic structure to
the untwisted-sector fold-wall geometry.
The product η · K = (2/9)(2/3) = 4/27 is the leading correction to θ13 = 0. The
result enters squared because θ13 couples the first and third generations, which requires
traversing two fold-wall transitions (from sector 1 across the intervening sector to
sector 3).
Proposition 2 (Squaring from double transition). The mixing element |Ue3|2 = sin2 θ13
involves a first-to-third generation transition. In the Z3 geometry, this requires two
sequential fold-wall hops: ν1 →fold wall →ν2 →fold wall →ν3. Each hop contributes
a factor of η K, giving the square (η K)2.
4
Parameter 22 — Solar Angle θ12
4.1
The formula
Theorem 2 (Solar angle).
sin2 θ12 = 1
3 −η2
2
= 1
3 −1
2
2
9
2
= 1
3 −2
81 = 27 −2
81
= 25
81.
(12)
6


4.2
Numerical comparison
sin2 θ12

pred = 25
81 = 0.3086.
(13)
The measured value (NuFIT 5.3 [2]):
sin2 θ12

meas = 0.307 ± 0.013.
(14)
Deviation:
|0.3086 −0.307|
0.307
= 0.53%.
(15)
4.3
Geometric derivation
The democratic Z3 symmetry gives sin2 θ12 = 1/3 (TBM value). The correction arises
from the finite thickness of the fold walls:
(i) The Donnelly invariant η = 2/9 measures the spectral asymmetry due to finite
fold-wall thickness. The correction to the solar angle is second-order in η because
θ12 connects the first and second generations, which are adjacent sectors — a
single fold-wall transition. The leading correction is therefore ∼η2.
(ii) The factor 1/2 arises from the two-body nature of the 1–2 overlap: the tunnelling
amplitude between two adjacent fold walls involves a symmetric two-state system,
introducing a factor of 1/2 in the perturbation expansion.
(iii) The correction is negative (reducing sin2 θ12 below 1/3) because the finite fold-wall
thickness partially localises the neutrino wavefunctions, reducing the overlap
between the ν1 and ν2 states.
Remark 3 (Perfect square). The result 25/81 = (5/9)2 is a perfect square of rationals.
This is not a coincidence: 5 = λ1 (the first non-trivial eigenvalue on S5) and 9 = 32 = p2,
so the solar angle encodes the ratio of the spectral gap to the square of the orbifold order.
5
Parameter 23 — Atmospheric Angle θ23
5.1
The formula
Theorem 3 (Atmospheric angle).
sin2 θ23 =
d1
d1 + λ1
=
6
6 + 5 =
6
11.
(16)
5.2
Numerical comparison
sin2 θ23

pred = 6
11 = 0.5455.
(17)
7


The measured value (NuFIT 5.3 [2]):
sin2 θ23

meas = 0.546 ± 0.021.
(18)
Deviation:
|0.5455 −0.546|
0.546
= 0.10%.
(19)
5.3
Geometric derivation
The TBM value sin2 θ23 = 1/2 corresponds to maximal 2–3 mixing. The physical orbifold
deviates from maximality through a spectral impedance ratio at the fold interface:
(i) Available modes (d1 = 6). The ℓ= 1 eigenspace on S5 has degeneracy d1 = 6
(Supplement V, §1.1). These six modes constitute the tunnelling bandwidth —
the number of channels available for fold-wall transmission.
(ii) Eigenvalue gap (λ1 = 5). The ℓ= 1 eigenvalue λ1 = ℓ(ℓ+ 4)|ℓ=1 = 5 acts as
a spectral barrier. Modes with energy below λ1 are reflected; those above are
transmitted.
(iii) Impedance ratio. At the fold interface, the fraction of spectral weight carried
by the fold modes (transmitted) versus reflected by the spectral gap is
σ =
d1
d1 + λ1
=
6
11.
(20)
This ratio directly determines sin2 θ23.
Remark 4. The atmospheric angle is not maximal (1/2) but close to it. The deviation
6/11 −1/2 = 1/22 ≈0.045 is a direct measure of the spectral impedance mismatch
between the degeneracy and eigenvalue at ℓ= 1. Current experiments are approaching the
precision needed to confirm or exclude maximality; the prediction 6/11 is distinguishable
from 1/2 at the ∼2σ level with projected NOvA and T2K sensitivities.
6
Parameter 24 — Leptonic CP Phase δPMNS
CP
6.1
The formula
Theorem 4 (Leptonic CP phase).
δPMNS
CP
= p γ = 3 arctan
2π2
9

= 3 × 65.5◦= 196.5◦.
(21)
8


6.2
Numerical comparison
δPMNS
CP

pred = 196.5◦.
(22)
The measured value (combined T2K/NOvA, normal ordering [2]):
δPMNS
CP

meas = 195◦± 50◦.
(23)
Deviation:
|196.5 −195|
195
≈0.8%.
(24)
6.3
Geometric derivation
The quark-sector CP phase (the CKM angle γ) was derived in the baryon-sector
supplements as
γ = arctan
2π2
9

≈65.5◦.
(25)
This is the CP-violating phase acquired by a twisted-sector fermion traversing a single
fold wall of the Z3 orbifold.
The leptonic CP phase differs by a factor of p = 3 (the order of the Z3 group):
(i) Quarks (twisted sector): Quarks are pinned to the cone point and interact with
only one fold wall at a time. The CP phase is γ (single fold-wall contribution).
(ii) Neutrinos (untwisted sector): Neutrinos are neutral under Z3 and propagate
freely through the bulk. They access all three fold walls simultaneously. The CP
phase is the coherent sum of three single-wall contributions:
δPMNS
CP
= 3 γ = p γ.
(26)
Remark 5. The factor of 3 is exact (it is the group order p = |Z3|, not an ap-
proximation). The leptonic CP phase is predicted to be close to but not exactly 180◦
(which would correspond to maximal CP violation in the PMNS matrix). The deviation
196.5◦−180◦= 16.5◦encodes the non-trivial value of arctan(2π2/9).
7
The Point/Side/Face Forward Derivation
The central insight of the neutrino sector is that the three neutrino mass eigenstates are
not three copies of one geometric object (as the charged leptons are three copies of the
cone-point mode with different Z3 charges). Instead, they are three different geometric
objects within the Z3 orbifold:
9


7.1
The three geometric objects
Definition 3 (Point/Side/Face assignment).
• ν1 = Point (cone tip, 0-dimensional). Minimal geometric support ⇒lightest
mass (m1 ≈0).
• ν2 = Side (fold wall, codimension-1 surface). Moderate tunnelling overlap ⇒
intermediate mass.
• ν3 = Face (sector bulk, full-dimensional). Maximum tunnelling overlap ⇒heaviest
mass.
This assignment is not arbitrary: it is forced by the geometric hierarchy of the orbifold.
The orbifold has exactly one cone point (0-dimensional), three fold walls (codimension-
1), and three sector bulks (full-dimensional). The mass hierarchy m1 ≪m2 < m3
reflects the hierarchy of geometric support.
7.2
Tunnelling overlap matrix
The tunnelling overlap matrix T encodes the amplitudes for neutrino transitions between
the three geometric objects. Its entries are determined by the orbifold invariants:
T =


T11
T12
T13
T12
T22
T23
T13
T23
T33

,
(27)
where each entry is derived from the geometric invariants η = 2/9, K = 2/3, σ =
d1/(d1 + λ1) = 6/11, p = 3, and the reactor angle (ηK)2 = 16/729.
Diagonal entries:
(i) T11 = +2σ
p = +2
3 · 6
11 = + 4
11
(point self-overlap: the cone tip sees both twisted
sectors, giving a factor of 2, divided by the orbifold order p).
(ii) T22 = −σ
p = −1
3 · 6
11 = −2
11
(side self-overlap: the fold wall depletes its source
state; a single twisted-sector contribution divided by p; sign negative because the
tunnelling removes amplitude from the source).
(iii) T33 = −(d1 −2) · (ηK)2 = −4 · 16
729 = −64
729
(face self-overlap: the bulk sector
leaks through d1 −2 = 4 transverse channels, each carrying the reactor-angle
amplitude).
Off-diagonal entries:
10


(i) T12 = −η = −2
9
(point ↔side: the Donnelly invariant governs the transition
between the cone tip and the adjacent fold wall).
(ii) T23 = −K σ = −2
3 · 6
11 = −4
11
(side ↔face: the harmonic lock K modulates
the impedance ratio σ).
(iii) T13 = +η K = + 4
27
(point ↔face: requires two hops through the intervening
side; the sign is positive because two negative transitions compose to a positive
amplitude).
7.3
Mass matrix construction
The perturbation strength is
ε =
p
d1 + λ1 =
√
11.
(28)
The full neutrino mass matrix (in flavour basis) is
Mν = D +
√
11 T,
(29)
where D is the 3 × 3 democratic matrix (all entries 1/3). The democratic matrix
provides the zeroth-order (TBM) structure, and
√
11 T is the geometric correction.
7.4
Diagonalisation and PMNS extraction
Numerical diagonalisation of Mν (29) yields the PMNS mixing matrix. The predicted
mixing parameters are compared with experiment in Table 2.
Observable
Predicted
PDG/NuFIT
Deviation
sin2 θ13
0.0216
0.0220 ± 0.0007
−1.8%
sin2 θ12
0.303
0.307 ± 0.013
−1.3%
sin2 θ23
0.537
0.546 ± 0.021
−1.6%
Table 2: PMNS mixing angles from the point/side/face forward model compared with
NuFIT 5.3 data [2]. All predictions are within 2% of the central measured values.
Remark 6. The point/side/face model uses no free parameters beyond the orbifold
invariants (η, K, d1, λ1, p) already determined in previous supplements. The small
deviations (≲2%) from the exact per-parameter formulae (Sections 3–5) reflect the
difference between isolated perturbation theory (exact formulae) and the full matrix
diagonalisation (which includes cross-couplings).
11


8
Parameter 25 — Neutrino Mass from the Inver-
sion Principle
8.1
Bulk resonance versus boundary tunnelling
The proton and the neutrino represent two complementary aspects of the orbifold
geometry:
• Proton = bulk resonance (constructive interference, interior standing wave). The
proton-to-electron mass ratio is
mp
me
= 6π5
(30)
(Supplement IV).
• Neutrino = boundary tunnelling mode (evanescent wave, exponentially sup-
pressed in the bulk). The neutrino loses mass as the inverse of the squared bulk
volume, shared among p sectors.
8.2
The master formula
Theorem 5 (Neutrino mass inversion).
p m2
p mν = m3
e.
(31)
Solving for the heaviest neutrino mass m3:
m3 =
m3
e
p m2
p
=
me
p (mp/me)2 =
me
3 (6π5)2 =
me
108 π10,
(32)
where 108 π10 = 3 × (6π5)2.
8.3
Numerical verification
6π5 = 1836.12,
(33)
(6π5)2 = 3,371,340,
(34)
108 π10 = 3 × 3,371,340 = 10,114,021,
(35)
m3

pred = 0.51100 MeV
10,114,021
= 5.052 × 10−8 MeV = 50.52 meV.
(36)
The measured value (from oscillation data, normal ordering):
m3

meas =
q
∆m2
32 + ∆m2
21 ≈
√
2.453 × 10−3 + 7.53 × 10−5 eV = 50.28 meV.
(37)
12


Deviation:
|50.52 −50.28|
50.28
= 0.48%.
(38)
8.4
Connection to the seesaw mechanism
The inversion formula (31) can be rewritten in a form reminiscent of the Type-I seesaw:
mν =
m3
e
p m2
p
= y2
e v2
MR
,
(39)
where the effective right-handed scale is
MR = p m2
p
me
≈3 × (0.938 GeV)2
0.511 × 10−3 GeV ≈5.2 × 106 GeV ≈5.2 × 1015 eV.
(40)
Remark 7. The scale MR ∼5 × 1015 eV is geometric in origin: it represents the
fold-wall penetration depth in the orbifold geometry. It is not the mass of a physical
right-handed neutrino particle. The seesaw form mν = y2v2/MR is recovered as a
mathematical identity, but the underlying physics is tunnelling suppression rather than
heavy-particle exchange.
9
Parameter 26 — Mass-Squared Splitting Ratio
9.1
The formula
Theorem 6 (Mass-squared ratio).
∆m2
32
∆m2
21
= d2
1 −p = 36 −3 = 33.
(41)
9.2
Numerical comparison
∆m2
32
∆m2
21

pred
= 33.
(42)
The measured value (NuFIT 5.3 [2]):
∆m2
32
∆m2
21

meas
= 2.453 × 10−3 eV2
7.53 × 10−5 eV2 = 32.58 ± 0.80.
(43)
Deviation:
|33 −32.58|
32.58
= 1.3%.
(44)
13


9.3
Geometric derivation
The ratio of mass-squared splittings reflects the ratio of tunnelling bandwidths for the
atmospheric and solar channels:
(i) Atmospheric splitting (∆m2
32). The 2–3 transition involves the full tunnelling
bandwidth of the ℓ= 1 sector. The number of available two-body tunnelling
channels scales as d2
1 = 62 = 36, counting all pairwise mode combinations.
(ii) Solar splitting (∆m2
21). The 1–2 transition is a subtler, single fold-wall process.
Its bandwidth is the baseline against which the atmospheric bandwidth is measured.
(iii) Three-fold sharing. The ratio is reduced by p = 3 because the three sectors of
the Z3 orbifold share the total tunnelling bandwidth equally:
∆m2
32
∆m2
21
= d2
1 −p = 36 −3 = 33.
(45)
9.4
Derived masses from P25 + P26
Combining the heaviest neutrino mass (Parameter 25) with the splitting ratio (Parame-
ter 26) determines the full mass spectrum:
m3 = 50.52 meV,
(46)
m2
3 −m2
2
m2
2 −m2
1
= 33.
(47)
In the normal-ordering limit m1 ≈0:
m2
3 −m2
2 = 33 m2
2,
(48)
m2
3 = 34 m2
2,
(49)
m2 = m3
√
34 = 50.52
√
34 = 50.52
5.831 = 8.66 meV.
(50)
The measured value:
m2

meas =
q
∆m2
21 =
√
7.53 × 10−5 eV = 8.68 meV.
(51)
Deviation:
|8.66 −8.68|
8.68
= 0.23%.
(52)
The mass spectrum is:
m1 ≈0,
m2 = 8.66 meV,
m3 = 50.52 meV.
(53)
14


This is the normal hierarchy, predicted by the point/side/face geometric assignment
(the point mode has minimal support and therefore minimal mass).
The sum of neutrino masses:
X
mν ≈0 + 8.66 + 50.52 = 59.2 meV.
(54)
10
Why No Neutrino Koide Ratio
The charged-lepton Koide ratio K = 2/3 is exact. One might ask whether an analogous
relation holds for neutrinos. It does not, and the reason is structural.
10.1
Charged leptons: circulant symmetry
The three charged-lepton masses are eigenvalues of a 3 × 3 Hermitian circulant matrix
(Supplement II, §1.1). The circulant structure arises because all three generations are
the same geometric object (the cone-point mode) carrying different Z3 charges ω0, ω1, ω2.
The Koide identity
K =
P mk
 P √mk
2 = 2
3
(55)
is a trace identity for circulant matrices with eigenvalue modulus r =
√
2, valid for any
value of the phase δ.
10.2
Neutrinos: three different objects
The three neutrino mass eigenstates are three different geometric objects (Definition 3):
• ν1 = point (0-dimensional),
• ν2 = side (codimension-1),
• ν3 = face (full-dimensional).
Their mass matrix (29) is not circulant: T11 ̸= T22 ̸= T33 and the off-diagonal entries
are not related by cyclic permutation. Therefore the circulant trace identity does not
apply.
10.3
The neutrino Qν
Computing the Koide-like ratio for neutrinos (using m1 ≈0, m2 = 8.66 meV, m3 =
50.52 meV):
Qν =
m1 + m2 + m3
(√m1 + √m2 + √m3)2 ≈
59.2
(0 + 2.943 + 7.108)2 =
59.2
101.0 ≈0.586.
(56)
15


Proposition 3 (Qν ̸= 2/3 is a prediction). The inequality Qν ̸= 2/3 is not a failure of
the framework; it is a prediction. The value Qν ≈0.586 follows from the point/side/face
mass hierarchy. Measurement of a neutrino Koide-like ratio close to 2/3 would falsify
the geometric model.
11
Derived Predictions
The six parameters derived in this supplement, combined with the mass spectrum, yield
two cosmologically testable predictions.
11.1
Sum of neutrino masses
X
mν = m1 + m2 + m3 ≈0 + 8.66 + 50.52 = 59.2 meV.
(57)
This lies squarely within the sensitivity window of forthcoming cosmological surveys.
The DESI baryon acoustic oscillation programme and the Euclid satellite are projected
to constrain P mν to ∼50–70 meV at 95% CL [1]. A measured value significantly
below 50 meV (inverted hierarchy) or significantly above 70 meV would be in tension
with the prediction.
11.2
Effective Majorana mass
If neutrinos are Majorana particles, the effective mass governing neutrinoless double-beta
decay (0νββ) is
|mββ| =

X
k
U 2
ek mk
 .
(58)
In the normal hierarchy with m1 ≈0:
|mββ| ≈2–3 meV.
(59)
This is below the sensitivity of current experiments (KamLAND-Zen, |mββ| < 36–
156 meV) but within the projected reach of next-generation experiments:
• nEXO: sensitivity ∼5–17 meV,
• LEGEND-1000: sensitivity ∼9–21 meV.
A positive signal in the 2–3 meV range would require further detector improvements
beyond the next generation. However, a signal above ∼10 meV in normal ordering
would be in tension with this framework.
16


12
Provenance Table
Table 3 maps every result in this supplement to its mathematical source, verification
method, and epistemic status.
12.1
Connection to the lotus potential
The neutrino sector parameters (P21–P26) are evaluated at the lotus point ϕlotus =
0.9574 of the fold potential V (ϕ) = λHv4
max(ϕ2 −ϕ2
lotus)2/4. The neutrino masses arise
from the petal overlap: ghost wavefunctions bleed through the 4.3% residual fold opening,
with the lightest neutrino m1 ≈0 corresponding to the minimal tunneling amplitude at
the cone tip (point mode). The cosmological constant Λ1/4 = mν3η2 = 2.49 meV is the
infinitesimal breathing energy of the lotus — nonzero because ϕ < 1 (the fold never
fully closes).
References
[1] R. L. Workman et al. (Particle Data Group), “Review of Particle Physics,” Prog.
Theor. Exp. Phys. 2022 (2022) 083C01, and 2024 update.
[2] I. Esteban, M. C. Gonzalez-Garcia, M. Maltoni, T. Schwetz, and A. Zhou, “The
fate of hints: updated global analysis of three-flavour neutrino oscillations,” J.
High Energy Phys. 09 (2020) 178; NuFIT 5.3 (2024), www.nu-fit.org.
[3] H. Donnelly, “Eta invariants for G-spaces,” Indiana Univ. Math. J. 27 (1978)
889–918.
17


Result
Mathematical
Source
Verification
Status
Twisted/untwisted sec-
tors (Table 1)
Z3 orbifold topol-
ogy
Structural identifi-
cation
Framework
TBM from Z3 (Prop. 1)
Eigenvectors
of
cyclic permutation
Direct linear alge-
bra
Theorem
sin2 θ13
=
16/729
(Thm. 1)
(ηK)2;
junction
asymmetry
0.02194 vs 0.02200
(0.27%)
Prediction
sin2 θ12
=
25/81
(Thm. 2)
1/3 −η2/2;
fold-
wall thickness
0.3086
vs
0.307
(0.53%)
Prediction
sin2 θ23
=
6/11
(Thm. 3)
d1/(d1
+
λ1);
impedance ratio
0.5455
vs
0.546
(0.10%)
Prediction
δCP = 196.5◦(Thm. 4)
3 arctan(2π2/9);
p× quark phase
196.5◦
vs
195◦
(0.8%)
Prediction
Point/side/face model
(Def. 3)
Orbifold geometric
hierarchy
Table 2 (all < 2%)
Derived
Tunnelling matrix T
(Eq. 27)
Orbifold invariants
η, K, σ, p
PMNS extraction
Derived
m3
=
me/(108π10)
(Thm. 5)
Inversion:
p m2
p mν = m3
e
50.52
vs
50.28
meV
(0.48%)
Prediction
∆m2
32/∆m2
21
=
33
(Thm. 6)
d2
1 −p; tunnelling
bandwidth
33 vs 32.58 (1.3%)
Prediction
m2
=
m3/
√
34 (de-
rived)
P25 + P26 combi-
nation
8.66 vs 8.68 meV
(0.23%)
Derived
Qν ̸= 2/3 (Prop. 3)
Non-circulant mass
matrix
Qν ≈0.586
Prediction
P mν ≈59.2 meV
Mass spectrum sum
DESI/Euclid win-
dow
Prediction
|mββ| ∼2–3 meV
PMNS elements +
masses
nEXO/LEGEND-
1000 reach
Prediction
Seesaw form recovered
MR = p m2
p/me
MR
∼
5.2 ×
1015 eV
Consistency
Table 3: Provenance map for Supplement VII results (Parameters 21–26 and derived
predictions). “Theorem” entries follow from established mathematics. “Derived” entries
follow algebraically from prior results.
“Prediction” entries are compared against
PDG/NuFIT measurements. “Framework” entries depend on the orbifold identification.
“Consistency” entries recover known structures.
18
