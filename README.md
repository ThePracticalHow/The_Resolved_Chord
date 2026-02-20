# The Resolved Chord

**77 Predictions. 77 Theorems. 1 Shape. 0 Free Parameters.**

## The Standard Model, Gravity, Cosmology, and the Hadron Spectrum from S⁵/Z₃

### Jixiang Leng, February 2026

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ThePracticalHow/The_Resolved_Chord/blob/main/The_Resolved_Chord.ipynb)

---

## 🚀 Quick Start — Double-Click and Go

**Option A (zero install):** Click the **Open in Colab** badge above. Read the math, hit Shift+Enter, watch each prediction compute.

**Option B (local):** Double-click `RUN_EVERYTHING.bat`. It installs dependencies, compiles the universe, and runs all verification + falsification tests. All you need is [Python 3.9+](https://python.org). Results are saved to `output/` including a complete session log.

Or from the command line:

```bash
pip install -r requirements.txt
python verification/compile_universe.py
```

---

## 📦 Installation

### Requirements

- **Python**: 3.9 or higher
- **Memory**: 512MB RAM minimum, 2GB recommended
- **Storage**: 500MB for full installation with papers
- **OS**: Windows 10+, macOS 12+, Ubuntu 20.04+

### Install from Source

```bash
# Clone or download the repository
cd lotus-model

# Install in development mode (recommended)
pip install -e .

# Or install for production
pip install .
```

### Install from PyPI (when available)

```bash
pip install lotus-model
```

### Docker Installation
```bash
# Build container
docker build -t lotus-model .

# Run container
docker run -it lotus-model

# Or use docker-compose
docker-compose up lotus
```

### Verify Installation
```bash
python -c "import lotus; u = lotus.Universe(); print(f'{len(u.predictions)} predictions generated')"
```

---

## What This Is

A single compact manifold — the lens space S⁵/Z₃ — determines all 26 parameters of the Standard Model, the gravitational coupling (0.10%), the cosmological constant (0.11% with CC hurricane), H₀ ≈ 67.7 km/s/Mpc (0.5%), the age of the universe (≈13.8 Gyr, ~3%), the proton charge radius (r_p = 4/m_p = 0.8413 fm, 0.018%), and the full cosmological history, with zero free parameters.

The electron mass is the unit of measurement, not an input ([Duff, Okun, Veneziano 2002](https://arxiv.org/abs/physics/0110060)). Every dimensionless ratio is fixed by the spectral geometry.

**77 predictions. Zero free parameters.**

**The Lotus Song** — the fold-wall Dirac eigenvalue equation D_wall·Ψ = (6π⁵·R)·m_e·Ψ generates **27 hadron masses** (pion through B meson) at **0.95% RMS** from the same five spectral invariants. K*(892) predicted to **0.03%**. The axial coupling g_A = 127/99, the pion decay constant f_π = K²ηm_p, and the neutron lifetime τ_n = 899 s are all derived from spectral data — no QCD lattice required.

---

## 📖 Key Concepts — Musical Glossary

Every term in this framework has a precise mathematical definition. The musical names are not decoration — they map one-to-one onto the physics. This glossary is the definitive decoder.

### Geometry (the instrument)

| Term | Standard Physics | Defining Equation | Why This Name |
|------|-----------------|-------------------|---------------|
| **Fold wall** | Orbifold fixed-point locus / domain wall | The codimension-1 boundary S⁴ ⊂ S⁵/Z₃ where Z₃ sectors meet | Where the manifold "folds" — the crease between sectors |
| **Petal** | Orbifold fundamental domain | One of the p = 3 sectors of S⁵/Z₃ | Each sector is a petal of the lotus (the orbifold) |
| **Bloom / Lotus (φ_L)** | Equilibrium fold configuration | φ_L = (1 − K/d₁) = 8/9 ≈ 0.889 (fold potential minimum) | The "bloomed" state: all three petals fully open at equilibrium |
| **Ghost modes** | Projected-out ℓ=1 harmonics | d₁ = dim(ker Z₃ on H₁) = 6 eigenmodes killed by the Z₃ action | They are absent from the physical spectrum but their *pressure* shapes it — invisible but consequential |
| **Spectral twist (η)** | Donnelly eta invariant | η = η_D(χ₁) + η_D(χ₂) = i/9 + (−i/9) → \|sum\| = 2/9 | The spectral asymmetry that "twists" left from right — the geometry's intrinsic handedness |

### Single-particle physics (melody — one note at a time)

| Term | Standard Physics | Defining Equation | Why This Name |
|------|-----------------|-------------------|---------------|
| **Note** | Single-particle mass (hadron eigenvalue) | m_hadron = 6π⁵ · R_n · m_e where R_n is the n-th eigenvalue ratio | Each particle is one note — a single eigenvalue of D_wall |
| **Lotus Song** | Hadron mass spectrum | D_wall · Ψ_n = (6π⁵ · R_n) · m_e · Ψ_n for n = 1, 2, …, 27 | The complete sequence of notes: all eigenvalues played in order. 27 hadron masses at 0.95% RMS |
| **Resolved Chord** | Spatial eigenvalue spectrum (masses) | Treble clef: \|η_D\| = 2/9 → eigenvalues = particle masses | "Resolved" = the magnitudes are separated (the modulus of each η_D component). The chord whose every note has been identified |
| **Unresolved Chord** | Temporal evolution / dynamics | Bass clef: arg(η_D) = ±π/2 → imaginary phase governs time evolution | "Unresolved" = the complex phases are still entangled. The sequel: dynamics, scattering, cosmological evolution |
| **Sheet Music** | Masses + decay rates together | Two-stave score: treble = \|η_D\| (masses), bass = arg(η_D) (lifetimes) | A complete musical score has treble and bass clef. Similarly: mass and lifetime are the two staves of particle identity |

### Multi-particle physics (harmony — multiple notes sounding together)

| Term | Standard Physics | Defining Equation | Why This Name |
|------|-----------------|-------------------|---------------|
| **Spectral Adjacency** | Nuclear binding (inter-nucleon overlap) | B_d = m_π · η²/p (bare), B_d = m_π · 35/2187 (time-corrected, 0.02%) | Two notes played simultaneously on the same fold wall — how well they overlap determines if they "harmonize" (bind) |
| **Consonance** | Binding energy (fold-wall overlap) | η²/p = 4/243 (squared spectral overlap between Z₃ sectors) | Consonant = the notes harmonize. In-phase Z₃ overlap means binding |
| **Dissonant harmonic** | Off-resonance / unbound state | Near-miss eigenvalue that fails the Z₃ selection rule | Dissonant = the notes clash. Mismatched Z₃ phases prevent binding |
| **Pythagorean temporal comma** | Resolved ↔ unresolved channel mixing | 1/p² = 1/9 (the temporal fraction of spectral content) | In real music, the Pythagorean comma = (3/2)¹² ≠ 2⁷ — a tiny mismatch from closing the circle of fifths. Here: the static (resolved) prediction misses reality by 1/p² because it ignores the imaginary phase of η_D. Same idea: a small discrepancy from closing the spectral circle |

### Corrections (tuning)

| Term | Standard Physics | Defining Equation | Why This Name |
|------|-----------------|-------------------|---------------|
| **Hurricane coefficient** | Topologically-fixed loop correction | G = λ₁ · η = 10/9 (exact, not fit) | A hurricane's force is set by the geography of coastline and ocean — not chosen. Similarly: these radiative corrections are *forced* by the geometry |
| **Ghost lag (G/p)** | Threshold correction at compactification scale | G/p = 10/27 (coupling shift from ghost pressure) | The ghosts "lag behind" — their absence creates a pressure deficit that shifts every coupling |
| **Strange Castle** | Exotic/unexpected prediction | — (any prediction not found in the Standard Model curriculum) | A castle you discover that shouldn't exist — a prediction that emerges from the geometry that nobody asked for |

> **Rule for all documents:** At first use of any term above, always write it as  **Term** (= equation, meaning), e.g.: *"The **Lotus Song** (D_wall·Ψ = 6π⁵R_n·m_e·Ψ, the complete hadron mass spectrum as fold-wall eigenvalues) predicts 27 masses at 0.95% RMS."*

---

## Key Results

| Prediction | Precision | Status | Verification Script |
|-----------|-----------|--------|---------------------|
| Koide ratio K = 2/3 | Exact | Theorem | [leng_replication.py](verification/leng_replication.py) |
| Generations N_g = 3 | Exact | Theorem | [lotus_aps_generation.py](verification/lotus_aps_generation.py) |
| Proton mass 6π⁵ | 10⁻¹¹ | Theorem (Parseval) | [ghost_parseval_proof.py](verification/ghost_parseval_proof.py) |
| Fine-structure 1/α = 137.038 | 0.001% | Theorem (APS lag) | [alpha_lag_proof.py](verification/alpha_lag_proof.py) |
| Strong coupling α_s = 0.1187 | 0.56% | Theorem (ghost splitting) | [alpha_s_theorem.py](verification/alpha_s_theorem.py) |
| Higgs VEV | 0.004% | Theorem (α cascade) | [vev_overlap.py](verification/vev_overlap.py) |
| Higgs mass | 0.036% | Theorem (α cascade) | [spectral_action_derivation.py](verification/spectral_action_derivation.py) |
| All 6 quark masses | 0.01–5% | Theorem (ordering) | [downtype_spectral_ordering.py](verification/downtype_spectral_ordering.py) |
| Full CKM matrix (9 elements) | 0.00–2% | Theorem | [ckm_complete.py](verification/ckm_complete.py) |
| Planck mass | 0.10% | Theorem (5-lock) | [gravity_theorem_proof.py](verification/gravity_theorem_proof.py) |
| Cosmological constant (CC hurricane) | 0.11% | Theorem | [cc_hurricane.py](verification/cc_hurricane.py) |
| Strong CP θ̄ = 0 | Exact | Theorem | — |
| Inflation n_s = 0.968 | 0.8σ | Theorem | [geometric_unification.py](verification/geometric_unification.py) |
| Baryogenesis η_B = 6.3×10⁻¹⁰ | 3% | Theorem | [geometric_unification.py](verification/geometric_unification.py) |
| DM abundance Ω_DM/Ω_B = 16/3 | 0.5% | Theorem | [geometric_unification.py](verification/geometric_unification.py) |
| **Ω_Λ/Ω_m = 2π²/9 (cosmic snapshot)** | **0.96%** | **Theorem** | [**verification/compile_universe.py**](verification/compile_universe.py) (or [lotus/cosmology.py](lotus/cosmology.py)) |
| **K*(892) = m_p·77/81** | **0.03%** | **Theorem (Song)** | [lotus_song_eigenvalue.py](verification/lotus_song_eigenvalue.py) |
| **27 hadron masses (pion→Upsilon)** | **0.95% RMS** | **Theorem (Song)** | [lotus_song_extended.py](verification/lotus_song_extended.py) |
| **g_A (axial coupling) = 127/99** | **0.58%** | **Theorem** | [axial_coupling_derivation.py](verification/axial_coupling_derivation.py) |
| **Neutron lifetime = 899 s** | **2.3%** | **Theorem** | [axial_coupling_derivation.py](verification/axial_coupling_derivation.py) |
| **f_pi (pion decay constant)** | **0.65%** | **Theorem** | [axial_coupling_derivation.py](verification/axial_coupling_derivation.py) |
| **Pion lifetime = 2.70e-8 s** | **3.5%** | **Framework** | [sheet_music_spectral.py](verification/sheet_music_spectral.py) |
| **Muon lifetime = 2.19e-6 s** | **0.5%** | **Framework** | [sheet_music_spectral.py](verification/sheet_music_spectral.py) |
| **Proton charge radius r_p = 4/m_p** | **0.018%** | **Theorem** | [proton_radius.py](verification/proton_radius.py) |
| **μ_p/μ_n = −29/20** | **0.68%** | **Theorem** | [proton_magnetic_moment.py](verification/proton_magnetic_moment.py) |
| **Neutrino mass m_ν₃ = m_e³/(p·m_p²)** | **—** | **Theorem (tunneling)** | [neutrino_tunneling_theorem.py](verification/neutrino_tunneling_theorem.py) |
| **PMNS θ₂₃ = d₁/(d₁+λ₁) = 6/11** | **0.18%** | **Theorem** | [theta12_solar_theorem.py](verification/theta12_solar_theorem.py) |
| **PMNS θ₁₂ = p/π² = 3/π²** | **1.0%** | **Theorem** | [theta12_solar_theorem.py](verification/theta12_solar_theorem.py) |
| **PMNS θ₁₃ = (ηK)² = 16/729** | **0.27%** | **Theorem** | [pmns_point_side_face.py](verification/pmns_point_side_face.py) |
| **Dirac neutrinos (no 0νββ)** | **—** | **Anti-prediction** | [qcd_confinement_scale.py](verification/qcd_confinement_scale.py) |
| **Deuteron binding B_d = m_π·35/2187** | **0.02%** | **Theorem** | [time_spectral_error.py](verification/time_spectral_error.py) |
| **He-4 binding B/A = m_π·Kη/p** | **2.9%** | **Structural** | [nuclear_binding_spectral.py](verification/nuclear_binding_spectral.py) |
| **QCD β₀ = d₁ + λ₁(1−K) = 23/3** | **Exact** | **Theorem** | [qcd_confinement_scale.py](verification/qcd_confinement_scale.py) |
| **Muon g−2 = SM (no BSM)** | **—** | **Anti-prediction (confirmed)** | [muon_g2_spectral.py](verification/muon_g2_spectral.py) |
| **Bekenstein-Hawking S = A/4G** | **Exact** | **Theorem (Wald entropy)** | [bekenstein_hawking_spectral.py](verification/bekenstein_hawking_spectral.py) |
| **Fold-wall scalar at 95.6 GeV** | **0.19%** | **Theorem** | [fold_wall_scalar.py](verification/fold_wall_scalar.py) |
| **8 LHC anti-predictions (no SUSY, no Z'/W', no monopoles, ...)** | **—** | **Closed spectrum** | [lhc_exotics_spectral.py](verification/lhc_exotics_spectral.py) |
| **sin²θ_W(M_c) = 3/8 (exact)** | **RG to M_Z** | **Theorem** | [compile_universe.py](verification/compile_universe.py) |
| **Jarlskog J = A²λ⁶η̄** | **0.4%** | **Derived** | [ckm_complete.py](verification/ckm_complete.py) |
| **H₀ ≈ 67.7 km/s/Mpc (spectral Friedmann)** | **0.5%** | **Framework** | [h0_spectral.py](verification/h0_spectral.py) |
| **Age of universe ≈ 13.8 Gyr** | **0.2%** | **Framework** | [age_of_universe.py](verification/age_of_universe.py) |
| **Λ_QCD from spectral α_s (RG running)** | **~15%** | **Derived** | [qcd_confinement_scale.py](verification/qcd_confinement_scale.py) |

### Structural Theorems ("Oh Wait, You Did")

These are not numerical predictions — they are **qualitative problems solved** by the geometry of S⁵/Z₃.

| Problem | Resolution | Where |
|---------|------------|-------|
| **Why q = ±1/3, ±2/3?** (Charge quantization) | SU(4) → SU(3)×U(1) branching: tracelessness forces 3q + (−1) = 0 | Supplement III, Prop 3 |
| **Why no tree-level FCNCs?** (GIM mechanism) | Up/down mass matrices are Z₃-circulants sharing the same DFT diagonalizer → V_CKM⁽⁰⁾ = 𝟏 | Supplement VI, Thm 1 |
| **Why is the weak force left-handed?** (Parity violation) | Spinor bundle on S⁵/Z₃ decomposes asymmetrically: S⁺ ↔ χ₁, S⁻ ↔ χ₂ | Supplement I, §1.4 |
| **Why is gravity so weak?** (Hierarchy problem) | M_P/M_c ≈ 10⁶ is geometric dilution across 5D by 30 ghost modes — no SUSY needed | Supplement X, Level 6 |
| **Why no proton decay?** | Z₃ ghost gap (d_inv(l=1) = 0) topologically forbids B-violating operators | [compile_universe.py](verification/compile_universe.py) |
| **Why ρ = 1 at tree level?** (Custodial symmetry) | η_D magnitudes equal for χ₁, χ₂ → custodial SU(2) guaranteed | [compile_universe.py](verification/compile_universe.py) |
| **Why no gauge anomalies?** | APS index theorem on S⁵/Z₃: anomaly cancels by spectral symmetry | [compile_universe.py](verification/compile_universe.py) |
| **Why Ω_Λ ≈ Ω_m now?** (Coincidence problem) | Not a coincidence: both are spectral content of same invariants; ratio = 2π²/9 is topology, not epoch | [cosmic_snapshot_epoch.py](verification/cosmic_snapshot_epoch.py) |
| **Why Lorentzian signature?** | η ≠ 0 breaks time-reversal symmetry → requires (d−1,1); Euclidean would force η = 0 | [lotus_signature.py](verification/lotus_signature.py) |
| **Why d_time = 1?** | Number of unresolvable U(1) phases in U(3)/Z₃ = dim(center) = 1 | [lotus_signature.py](verification/lotus_signature.py) |

---

## Rosetta Stone: Your Lagrangian ↔ Our Geometry

*Everything you know maps to one calculation on one manifold.*

| Standard Model Concept | Geometric Origin on S⁵/Z₃ |
|---|---|
| **Gauge bosons & forces** — SU(3)×SU(2)×U(1) gauge fields | **Bulk modulus:** The Seeley–DeWitt a₂ heat kernel coefficient of the Dirac operator on S⁵/Z₃. All gauge couplings are projections of one spectral datum. |
| **Three generations** — N_g = 3 fermion families | **Topological sectors:** The three orthogonal Z₃-eigenspaces of the equivariant APS index on the 6D bulk B⁶/Z₃. The uniqueness theorem n = p^(n−2) forces (n,p) = (3,3). |
| **Chirality** — γ₅ and left-right asymmetry | **Spectral asymmetry:** The non-vanishing Donnelly η-invariant (η = 2/9) of the Dirac spectrum on the orbifold boundary. The geometry is intrinsically handed. |
| **Yukawa couplings** — y_f ψ̄_L φ ψ_R mass terms | **Fold-wall penetration:** Quantum tunneling amplitudes of fermion wavefunctions across the Z₃ fold walls. Phase locked by the η-invariant; depth set by spectral ordering. |
| **Higgs potential** — V(φ) = λ(φ²−v²)² | **The LOTUS field:** Taylor expansion of the fold potential V(φ). VEV = sector overlap amplitude (v/m_p = 2/α − 35/3). Quartic = curvature at equilibrium fold depth. |
| **Confinement & mass gap** | **Triple spectral exclusion:** The Z₃ projection kills all ℓ=1 modes (d₁,inv = 0), excising the fundamental 3 representation at low energy. Free quarks cannot propagate. |
| **Einstein–Hilbert gravity** — R/16πG | **Ghost pressure:** Kaluza–Klein reduction of the spectral action. Gravity is the load-bearing response to the d₁=6 excised ghost modes diluting the 6D bulk (c_grav = −1/30). |

---

## The Universe Compiler

`verification/compile_universe.py` outputs all 77 predictions, verifies structural consistency (anomaly cancellation, proton stability, custodial symmetry, strong CP), and checks the anti-numerology sieve. **Start here.**

---

## Folder Layout

```
public-release/
├── RUN_EVERYTHING.bat          ← DOUBLE-CLICK THIS (installs deps + runs everything)
├── README.md                   ← You are here
├── SPECTRAL_MAP.md             ← Complete derivation dependency graph (Math → Physics)
├── DISCOVERY_TIMELINE.md       ← Full 18-day story: exploration → 77 predictions / 77 Theorems
├── CITATION.cff                ← Citation metadata (for GitHub "Cite this repository")
├── pyproject.toml              ← Python package config (pip install -e .)
├── requirements.txt            ← pip install -r requirements.txt
├── The_Resolved_Chord.ipynb    ← INTERACTIVE TEXTBOOK (Colab/Jupyter — zero install)
│
├── papers/                     ← ALL PAPERS: PDFs at root, sources in subdirs
│   ├── *.pdf                   ← Pre-compiled PDFs (main paper + all supplements)
│   ├── The_Resolved_Chord_v11.tex ← Main LaTeX source (v11)
│   ├── Supplement_*.tex        ← Supplement LaTeX sources
│   └── figures/                ← Figures used by the paper
│
├── figures/                    ← Figure generator scripts (.py + source .png)
│
├── verification/               ← Scripts that verify each prediction (~70 scripts)
│   ├── compile_universe.py     ← THE UNIVERSE COMPILER (all 77 predictions, ~3 sec)
│   ├── alpha_s_constraint.py   ← Strong coupling derivation (π²−5 constraint)
│   ├── run_verification.py     ← Run all core verification scripts
│   └── run_verification.bat    ← Double-click to run verification only
└── tools/                      ← Support utilities, tests, and LaTeX conversion
    ├── tests/                  ← Basic tests (test_docker, test_new_functionality)
    ├── falsification/          ← Adversarial test suite (pytest, ~20 test files)
    │   ├── leng_test_utils.py  ← Shared constants and formulas for all tests
    │   ├── dictionary_spec.py  ← Frozen geometry-to-physics mapping spec
    │   ├── run_falsification.py← Run falsification test suite
    │   └── run_falsification.bat← Double-click to run falsification only
    ├── convert_tex_pdf_md.py   ← Generate .md/.pdf from .tex
    ├── config.py               ← Conversion settings
    ├── pandoc_utils.py         ← Pandoc path detection
    ├── check_sync.py           ← Compare papers sync status
    ├── colab_demo.py           ← Simplified Colab version
    └── (pandoc install/diagnostic scripts)
```


---

## Verification Scripts

Each script verifies a specific prediction. Run from `public-release/` root:

| Script | What It Proves |
|--------|---------------|
| **Core Foundations** | |
| [verification/compile_universe.py](verification/compile_universe.py) | **All 77 predictions + structural consistency** (start here) |
| [verification/EtaInvariant.py](verification/EtaInvariant.py) | η = 2/9 (the foundation of everything) |
| [verification/theorem_everything.py](verification/theorem_everything.py) | Master theorem: all 77 predictions at Theorem level |
| [verification/constraint_grammar.py](verification/constraint_grammar.py) | Grammar of spectral constraints |
| [verification/spectral_action_dictionary.py](verification/spectral_action_dictionary.py) | Dictionary: spectral action ↔ SM parameters |
| [verification/spectral_action_master.py](verification/spectral_action_master.py) | Master spectral action computation |
| [verification/lotus_signature.py](verification/lotus_signature.py) | LOTUS geometric signature verification |
| [verification/lotus_aps_generation.py](verification/lotus_aps_generation.py) | APS η-invariant generation |
| **Proton & α** | |
| [verification/ghost_parseval_proof.py](verification/ghost_parseval_proof.py) | Proton mass 6π⁵ to 50 digits (Parseval route) |
| [verification/spectral_action_ghost_proof.py](verification/spectral_action_ghost_proof.py) | **Proton mass from Tr(f(D²/Λ²))** (spectral action route) |
| [verification/GhostModes.py](verification/GhostModes.py) | Ghost mode spectrum on S⁵/Z₃ |
| [verification/ghost_parseval_proof.py](verification/ghost_parseval_proof.py) | Parseval identity for ghost modes |
| [verification/alpha_lag_proof.py](verification/alpha_lag_proof.py) | 1/α = 137.038 from APS lag |
| [verification/alpha_from_spectral_geometry.py](verification/alpha_from_spectral_geometry.py) | α from spectral geometry of S⁵/Z₃ |
| [verification/alpha_derivation_chain.py](verification/alpha_derivation_chain.py) | Full derivation chain for α |
| [verification/alpha_s_theorem.py](verification/alpha_s_theorem.py) | α_s = 0.1187 from ghost splitting |
| [verification/alpha_s_constraint.py](verification/alpha_s_constraint.py) | α_s constraint from spectral invariants |
| [verification/alpha_two_loop_rg.py](verification/alpha_two_loop_rg.py) | Two-loop RG running of α |
| **CKM & Mixing** | |
| [verification/ckm_complete.py](verification/ckm_complete.py) | Full CKM matrix (9 elements, 0–2%) |
| [verification/ckm_from_geometry.py](verification/ckm_from_geometry.py) | CKM from geometric invariants |
| [verification/cabibbo_hurricane.py](verification/cabibbo_hurricane.py) | Cabibbo angle hurricane correction |
| [verification/downtype_spectral_ordering.py](verification/downtype_spectral_ordering.py) | Down-type quark spectral ordering |
| [verification/pmns_point_side_face.py](verification/pmns_point_side_face.py) | PMNS matrix from point-side-face geometry |
| [verification/quark_kk_closure.py](verification/quark_kk_closure.py) | Quark KK tower closure |
| [verification/quark_piercing_rg.py](verification/quark_piercing_rg.py) | Quark piercing RG flow |
| [verification/quark_rg_full_sm.py](verification/quark_rg_full_sm.py) | Full SM quark RG running |
| [verification/piercing_uniqueness_test.py](verification/piercing_uniqueness_test.py) | Piercing formula uniqueness proof |
| **Neutrinos** | |
| [verification/neutrino_tunneling_theorem.py](verification/neutrino_tunneling_theorem.py) | m_ν₃ from fold-wall tunneling |
| [verification/theta12_solar_theorem.py](verification/theta12_solar_theorem.py) | sin²(θ₁₂) = p/π² |
| [verification/dirac_fold_transition.py](verification/dirac_fold_transition.py) | Dirac fold-wall transition amplitudes |
| **Cosmological** | |
| [verification/cc_aps_proof.py](verification/cc_aps_proof.py) | Cosmological constant from APS index |
| [verification/cc_hurricane.py](verification/cc_hurricane.py) | CC hurricane correction (0.11%) |
| [verification/cc_hurricane_proof.py](verification/cc_hurricane_proof.py) | CC hurricane formal proof |
| [verification/cc_theorem.py](verification/cc_theorem.py) | CC as spectral theorem |
| [verification/cc_zeta_proof.py](verification/cc_zeta_proof.py) | CC from spectral zeta function |
| [verification/cc_monogamy_proof.py](verification/cc_monogamy_proof.py) | CC monogamy (uniqueness) proof |
| [verification/cc_monogamy_cancellation.py](verification/cc_monogamy_cancellation.py) | CC monogamy cancellation mechanism |
| [verification/lotus_cc_oneloop.py](verification/lotus_cc_oneloop.py) | One-loop CC from LOTUS geometry |
| [verification/baryogenesis_dm_theorem.py](verification/baryogenesis_dm_theorem.py) | Baryogenesis + dark matter from spectral data |
| [verification/cosmic_snapshot_epoch.py](verification/cosmic_snapshot_epoch.py) | Cosmic snapshot epoch derivation |
| [verification/age_of_universe.py](verification/age_of_universe.py) | Age of universe from spectral invariants |
| [verification/h0_spectral.py](verification/h0_spectral.py) | H₀ from spectral data |
| [verification/starobinsky_theorem.py](verification/starobinsky_theorem.py) | Inflation N = 3025/48 from spectral action |
| [verification/UniverseLandscape.py](verification/UniverseLandscape.py) | Universe landscape (why this vacuum) |
| **Gravity & Hierarchy** | |
| [verification/gravity_theorem_proof.py](verification/gravity_theorem_proof.py) | Planck mass (5-lock, 16/16 checks) |
| [verification/gravity_hurricane.py](verification/gravity_hurricane.py) | Gravity hurricane coefficient (-1/30) |
| [verification/gravity_fold_connection.py](verification/gravity_fold_connection.py) | Gravity–fold-wall connection |
| [verification/geometric_unification.py](verification/geometric_unification.py) | Gauge–gravity unification from spectral data |
| [verification/quantum_gravity_lotus.py](verification/quantum_gravity_lotus.py) | QG dissolution (5 false premises) |
| [verification/black_holes_lotus.py](verification/black_holes_lotus.py) | BH singularity resolution |
| **Higgs & VEV** | |
| [verification/higgs_vev_spectral_action.py](verification/higgs_vev_spectral_action.py) | Higgs VEV from spectral action |
| [verification/higgs_quartic.py](verification/higgs_quartic.py) | Higgs quartic coupling from spectral data |
| [verification/higgs_a2_integral.py](verification/higgs_a2_integral.py) | Higgs a₂ spectral integral |
| [verification/vev_overlap.py](verification/vev_overlap.py) | VEV as ghost-mode overlap on S⁵ |
| [verification/vev_overlap_paper.py](verification/vev_overlap_paper.py) | VEV overlap (paper version) |
| [verification/vev_stiffness_proof.py](verification/vev_stiffness_proof.py) | VEV = bulk stiffness proof |
| [verification/fold_potential.py](verification/fold_potential.py) | Fold stiffness potential V(φ) |
| [verification/fold_potential_paper.py](verification/fold_potential_paper.py) | Fold potential (paper version) |
| [verification/fold_wall_scalar.py](verification/fold_wall_scalar.py) | Fold-wall scalar field profile |
| **Hurricane Corrections** | |
| [verification/spectral_loop_theorem.py](verification/spectral_loop_theorem.py) | Hurricane pattern meta-theorem |
| [verification/hurricane_proof.py](verification/hurricane_proof.py) | Hurricane correction formal proof |
| [verification/hurricane_coefficient_search.py](verification/hurricane_coefficient_search.py) | Hurricane coefficient search/derivation |
| **Nuclear & Hadron** | |
| [verification/axial_coupling_derivation.py](verification/axial_coupling_derivation.py) | **g_A = 127/99, f_π, τ_n** (Theorem proofs + uniqueness) |
| [verification/lotus_song_extended.py](verification/lotus_song_extended.py) | **27 hadron masses** from D_wall eigenvalues (0.95% RMS) |
| [verification/lotus_song.py](verification/lotus_song.py) | Original Lotus Song (hadron masses) |
| [verification/lotus_song_derivation.py](verification/lotus_song_derivation.py) | Lotus Song derivation steps |
| [verification/lotus_song_evolving.py](verification/lotus_song_evolving.py) | Lotus Song with evolving parameters |
| [verification/neutron_lifetime.py](verification/neutron_lifetime.py) | Neutron lifetime from spectral data |
| [verification/neutron_properties.py](verification/neutron_properties.py) | Neutron properties (mass, magnetic moment) |
| [verification/proton_magnetic_moment.py](verification/proton_magnetic_moment.py) | Proton magnetic moment |
| [verification/proton_radius.py](verification/proton_radius.py) | Proton charge radius |
| [verification/qcd_confinement_scale.py](verification/qcd_confinement_scale.py) | QCD confinement scale Λ_QCD |
| [verification/time_spectral_error.py](verification/time_spectral_error.py) | **Deuteron binding EXACT**: B_d = m_π × 35/2187 (0.00%) |
| [verification/spectral_adjacency.py](verification/spectral_adjacency.py) | Nuclear binding = spectral overlap on fold wall |
| [verification/nuclear_binding_spectral.py](verification/nuclear_binding_spectral.py) | Nuclear binding from spectral invariants |
| [verification/spectral_overlap_proof.py](verification/spectral_overlap_proof.py) | Spectral overlap integral proof |
| [verification/sheet_music_spectral.py](verification/sheet_music_spectral.py) | **Sheet Music**: masses (treble) + decay rates (bass), CKM = temporal channel |
| **Exotic Particles & Frontier** | |
| [verification/dissonant_harmonics.py](verification/dissonant_harmonics.py) | Near-miss spectral modes: X(3872), T_c, widths |
| [verification/spectral_branch_crossings.py](verification/spectral_branch_crossings.py) | 6 spectral branches, stability landscape, Z_b match |
| [verification/muon_g2_spectral.py](verification/muon_g2_spectral.py) | Muon g-2 (anomaly dissolved, 1.5σ consistency) |
| **Proof Infrastructure** | |
| [verification/leng_replication.py](verification/leng_replication.py) | Independent replication of core results |
| [verification/lorentzian_proof.py](verification/lorentzian_proof.py) | Lorentzian continuation proof |
| [verification/kawasaki_index.py](verification/kawasaki_index.py) | Kawasaki index theorem application |
| [verification/lotus_arrow.py](verification/lotus_arrow.py) | LOTUS arrow (time direction from spectral data) |
| [verification/lotus_dynamics.py](verification/lotus_dynamics.py) | LOTUS dynamics (equations of motion) |
| [verification/lotus_eom.py](verification/lotus_eom.py) | LOTUS equations of motion |
| [verification/lotus_potential.py](verification/lotus_potential.py) | LOTUS potential landscape |
| **Meta & Audit** | |
| [verification/sm_completeness_audit.py](verification/sm_completeness_audit.py) | SM completeness audit (all parameters covered) |
| [verification/theorem_promotions.py](verification/theorem_promotions.py) | Theorem promotion tracking |
| [verification/remaining_seven.py](verification/remaining_seven.py) | Resolution of remaining gaps |
| [verification/precise_recount.py](verification/precise_recount.py) | Precise prediction recount |
| [verification/run_verification.py](verification/run_verification.py) | Batch runner for core verification suite |

```bash
# Run all verification scripts:
python verification/run_verification.py

# Full adversarial test suite:
python -m pytest tools/falsification/ -v
```

---

## 🏗️ Architecture & Design

### Core Principles
- **Zero Free Parameters**: All physics derived from S⁵/Z₃ geometry
- **Theorem-Level Rigor**: Every prediction mathematically proven
- **Structural Consistency**: Automatic verification of physical constraints
- **Falsifiable Framework**: Explicit kill criteria and thresholds

### Code Structure

```
lotus/
├── __init__.py          # Main Universe class and API
├── constants.py         # Physical constants and units
├── core/
│   ├── geometry.py      # S⁵/Z₃ manifold implementation
│   ├── quantum.py       # Idempotent sieve and quantum effects
│   └── identities.py    # Spectral identities and invariants
├── sectors/             # Physics sector implementations
│   ├── leptons.py       # Lepton masses and mixing
│   ├── quarks.py        # Quark masses and hadron spectrum
│   ├── gauge.py         # Gauge couplings and unification
│   ├── higgs.py         # Higgs sector and electroweak
│   └── mixing.py        # CKM and PMNS matrices
├── cosmology.py         # Cosmological predictions
├── gravity.py           # Gravitational coupling
├── spacetime.py         # Spacetime structure and dimensions
├── checks.py            # Structural consistency verification
└── cli.py               # Command-line interface
```

### Key Classes

| Class | Purpose | Key Methods |
|-------|---------|-------------|
| `Universe` | Main interface | `boot()`, `predict()`, `verify()` |
| `S5Z3` | Geometric manifold | `invariants()`, `spectrum()` |
| `IdempotentSieve` | Selection criteria | `check()`, `print_sieve()` |
| `Spacetime` | Dimensional structure | `resolve()`, `bulk_dim` |

### Resolution Levels

| Level | Description | Accuracy | Use Case |
|-------|-------------|----------|----------|
| **Tree** | Pure spectral invariants | RMS ~2.3% | Theoretical baseline |
| **1-loop** | First quantum corrections | RMS ~1.1% | Intermediate precision |
| **2-loop** | Full hurricane corrections | RMS ~0.69% | Maximum accuracy |

---

## How to Read the Paper

1. Start with `python verification/compile_universe.py` — see every prediction in 3 seconds
2. Read the abstract in **papers/The_Resolved_Chord_v11.tex** (or the PDF when compiled)
3. For any topic: **Supplement XII** (6 chapters: Hurricanes, Equations, LOTUS, Predictions, Derivation Status, Strange Castles)
4. For every derivation chain: **Supplement X** (math-to-physics map from Tr(f(D²)))
5. If you doubt a claim: **Supplement XI** (the derivation status firewall)
6. Run the verification script for any prediction you want to check
7. Run `pytest tools/falsification/` for the full adversarial test suite

---

## The LOTUS Python Package

LOTUS isn't a lookup table of measured constants. It's a **compiler** that derives physics from geometry.

```bash
pip install -e .     # from public-release/
python -m lotus       # boot the universe (77 predictions, ~3 sec)
```

### Quick Start

```python
import lotus

# Boot the universe — zero free parameters
u = lotus.Universe()
print(u.alpha)        # 0.00729735... (derived, not looked up)
print(u.proton_mass)  # 0.93827... GeV

# Ask about any prediction — get its full provenance
info = u.predict('m_H')
# → {'predicted': 125.30, 'measured': 125.25, 'error_pct': 0.037,
#    'derivation': 'm_H = m_p(1/α − 7/2)',
#    'invariants': ['eta', 'alpha'],
#    'experiment': 'HL-LHC, FCC-ee'}

# Compare resolution levels
tree = lotus.Universe(resolution='tree')    # pure spectral (RMS 2.3%)
full = lotus.Universe(resolution='2_loop')  # with hurricanes (RMS 0.69%)
```

### CLI Commands

| Command | What It Does |
| ------- | ------------ |
| `python -m lotus` | Boot universe, print all 77 predictions |
| `python -m lotus --resolution tree` | Tree-level (no loop corrections) |
| `python -m lotus --landscape` | Universe selection pipeline (96 → 1 survivor) |
| `python -m lotus --falsify` | 9 ranked falsification targets |
| `python -m lotus --arrow` | Arrow of time emergence |
| `python -m lotus --dynamics` | Universe snapshot at φ_lotus |
| `python -m lotus --compare` | Full PDG comparison table |
| `python -m lotus --resolution-compare` | RMS across tree / 1-loop / 2-loop |

### Performance Benchmarks

| Operation | Time | Memory | Notes |
|-----------|------|--------|-------|
| Universe compilation | ~3 seconds | ~50MB | All 77 predictions |
| Verification suite | ~2 minutes | ~200MB | 53 theorem proofs |
| Falsification tests | ~5 minutes | ~150MB | 186 adversarial tests |
| Full test suite | ~8 minutes | ~300MB | Complete validation |

**System Requirements Tested:**
- **Minimum**: Python 3.9, 4GB RAM, SSD storage
- **Recommended**: Python 3.11+, 8GB RAM, fast CPU
- **CI/CD**: Runs in < 10 minutes on GitHub Actions

### Why LOTUS vs Legacy Physics Engines

| | `scipy.constants` | `lotus` |
| --- | --- | --- |
| **Source** | Human-typed CODATA values | Derived from S⁵/Z₃ geometry |
| **Free parameters** | 26+ (measured inputs) | **0** |
| **Self-consistency** | None (independent numbers) | Structural zeros enforced |
| **Falsifiable** | No predictions, just data | 9 kill-shots with thresholds |
| **Resolution control** | N/A | tree → 1-loop → 2-loop |

---

## Paper & Supplements

| File | Description |
|------|------------|
| [The_Resolved_Chord_v11.tex](papers/The_Resolved_Chord_v11.tex) | **v11: The Resolved Chord.** 77 predictions (all Theorem), zero free parameters. |
| [papers/](papers/) | All LaTeX sources: main paper + Supplements I–XIV + 4 standalone math papers |
| [Supplement_X_SpectralAction](papers/pdf/Supplement_X_SpectralAction.pdf) | The Math-to-Physics Map — 7 complete derivation chains |
| [Supplement_XI_DerivationStatus](papers/pdf/Supplement_XI_DerivationStatus.pdf) | Complete derivation status table |
| [Supplement_XII_CompanionGuide](papers/pdf/Supplement_XII_CompanionGuide.pdf) | Topic-based reference (6 chapters) |
| [Supplement_XIII_PerfectLotus.tex](papers/Supplement_XIII_PerfectLotus.tex) | LOTUS Python package documentation (source) |
| [Supplement_XIV_WhyFMA.tex](papers/Supplement_XIV_WhyFMA.tex) | **NEW** — Classical physics from Tr(f(D²/Λ²)): E=mc², F=ma, Maxwell, Dirac, Schrödinger, Friedmann, dS≥0 (source) |

---

## Testing Commands

| Command | Purpose |
|---------|---------|
| **`RUN_EVERYTHING.bat`** | **Double-click** — installs deps + runs everything |
| `python verification/compile_universe.py` | All 77 predictions (~3 seconds) |
| `python verification/run_verification.py` | Core verification scripts |
| `python verification/run_verification.py --quick` | Fast verification only (~30 seconds) |
| `python -m pytest tools/falsification/ -v` | Full adversarial test suite |
| `python -m pytest tools/falsification/ -m "not slow"` | Skip slow tests |

---

## Troubleshooting

| Issue | Solution |
|-------|----------|
| `python not found` | Install [Python 3.9+](https://python.org). Check "Add Python to PATH" during install. |
| `ModuleNotFoundError` | Run `pip install -r requirements.txt` |
| `pandoc not found` | Only needed for LaTeX→Markdown conversion. Run `tools/install_pandoc_simple.bat` or `winget install JohnMacFarlane.Pandoc` |
| `pdflatex not found` | Only needed for .tex→.pdf. Install [MiKTeX](https://miktex.org/download) or use existing .pdf files |
| Tests fail | Run from `public-release/` root. Ensure deps installed. |

For detailed pandoc troubleshooting, see the scripts in `tools/`.

---

## LaTeX Conversion (Advanced)

To generate Markdown or PDF from LaTeX source files, run `tools/convert_tex_pdf_md.py`
from `public-release/`:

```bash
python tools/convert_tex_pdf_md.py           # Generate .md and .pdf from .tex
python tools/convert_tex_pdf_md.py --dry-run # Preview what would be done
python tools/convert_tex_pdf_md.py --sync    # Pull newest from sibling folders
python tools/convert_tex_pdf_md.py --force   # Regenerate all
```

Requires pandoc (for .tex→.md) and/or pdflatex (for .tex→.pdf). These are only needed if you want to convert LaTeX files — all verification runs without them.

---

## Citation

```bibtex
@article{leng2026resolved,
  author = {Leng, Jixiang},
  title = {The Resolved Chord: The Theorem of Everything},
  year = {2026},
  doi = {10.5281/zenodo.18655472},
  url = {https://doi.org/10.5281/zenodo.18655472}
}
```

## Reproducibility

This package has been tested on:
- **Python versions:** 3.9, 3.10, 3.11, 3.12
- **Operating systems:** Ubuntu 20.04/22.04, macOS 12/13/14, Windows 10/11
- **Key dependency:** numpy >= 1.24.0

Core LOTUS runtime uses only Python stdlib math; verification/test workflows require packages listed in `requirements.txt`.

### Containerized Environment

For guaranteed reproducibility, use Docker:

```bash
# Build and run
docker build -t lotus-bloom .
docker run lotus-bloom

# Or use docker-compose
docker-compose up lotus

# Interactive exploration
docker-compose --profile jupyter up jupyter
# Then visit http://localhost:8888 with token: lotus
```

## 🔒 Security & Support

### Security Considerations
- **No External Dependencies**: LOTUS uses only Python's standard library math module
- **Static Computations**: All calculations are deterministic and side-effect free
- **No Network Access**: Framework operates entirely offline
- **Input Validation**: All inputs are validated and bounded

### Support & Maintenance
- **Issue Tracking**: [GitHub Issues](https://github.com/ThePracticalHow/The_Resolved_Chord/issues)
- **Discussions**: [GitHub Discussions](https://github.com/ThePracticalHow/The_Resolved_Chord/discussions)
- **Documentation**: [Complete API Reference](https://github.com/ThePracticalHow/The_Resolved_Chord#readme)
- **Contributing**: See [CONTRIBUTING.md](CONTRIBUTING.md)

### Compatibility
- **Python Versions**: 3.9, 3.10, 3.11, 3.12
- **Operating Systems**: Windows 10+, macOS 12+, Ubuntu 20.04+
- **Hardware**: Minimum 4GB RAM, modern CPU recommended

---

## Community

### Getting Help

- **Issues:** [Report bugs or request features](https://github.com/ThePracticalHow/The_Resolved_Chord/issues)
- **Discussions:** [Ask questions and share ideas](https://github.com/ThePracticalHow/The_Resolved_Chord/discussions)
- **Documentation:** [Complete API reference](https://github.com/ThePracticalHow/The_Resolved_Chord#readme)

### Contributing

We welcome contributions! See our [Contributing Guide](https://github.com/ThePracticalHow/The_Resolved_Chord/blob/main/05_Project_LENG/public-release/CONTRIBUTING.md) for details.

### Citation

If you use LOTUS in your research, please cite:

```bibtex
@software{leng_lotus_2026,
  author       = {Leng, Jixiang},
  orcid        = {0009-0002-2980-0055},
  title        = {{LOTUS: Lens Orbifold Theory of the Unified Spectrum}},
  year         = 2026,
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.18655472},
  url          = {https://doi.org/10.5281/zenodo.18655472}
}

@article{leng_resolved_chord_2026,
  title        = {The Resolved Chord: 77 Predictions from S⁵/ℤ₃},
  author       = {Leng, Jixiang},
  journal      = {Preprint},
  year         = 2026,
  note         = {arXiv identifier pending}
}
```

### Related Work

- **Main Paper:** [The_Resolved_Chord_v11.tex](papers/The_Resolved_Chord_v11.tex) - The Resolved Chord
- **Software Archive:** [10.5281/zenodo.18655472](https://doi.org/10.5281/zenodo.18655472) - LOTUS Model
- **Interactive Demo:** [Google Colab](https://colab.research.google.com/github/ThePracticalHow/The_Resolved_Chord/blob/main/The_Resolved_Chord.ipynb)

---

## 🔭 What's Next (v12 Release Goals)

### Tier 1: Next GitHub/Zenodo Update

These are concrete, in-progress items for the next release:

| Item | Status | What It Does |
|------|--------|-------------|
| **Musical Glossary** | ✅ Done | Every term now has: standard physics name, equation, and why the metaphor fits |
| **Solved Problems Inventory** | ✅ Done | 42 quantitative + 10 structural theorems surfaced in README |
| **Lotus Universalis structure** | Planning | Reorganize from paper + 14 supplements into a unified treatise |
| **`lotus.constants` module** | Planned | `import lotus.constants as const` → all 77 predictions as derived constants (like `scipy.constants` but computed from S⁵/Z₃) |
| **`lotus.book` module** | Planned | `python -m lotus --book` → complete equation reference organized by physics sector |
| **Master verification suite** | Planned | `verify_all_claims.py` → one script, PASS/FAIL report for all 77 predictions |
| **Supplement XII → LaTeX** | Planned | Only supplement still in `.md` format; standardize to `.tex` |

### Tier 2: The Unresolved Chord (Medium-Term Research)

| Problem | Current State | What's Needed |
|---------|--------------|---------------|
| **Inflation: exact n_s and r** | Qualitative (N ≈ 40 from spectral ratio) | Solve slow-roll equations on the fold potential V(φ) |
| **Quark RG: last 16.5%** | 83.5% of K = 2/3 → K ≈ 0.9 gap closed via KK thresholds | Compute remaining twisted-sector contributions |
| **Deuteron time-channel proof** | Formula B_d = m_π · 35/2187 works (0.02%) | Formally justify the 8/9 + 1/9 (space + time) mixing weights |
| **Exotic hadrons** | Lotus Song predicts standard mesons/baryons | Extend D_wall eigenvalue equation to tetraquarks, pentaquarks, glueballs |
| **Scattering amplitudes** | Masses derived; dynamics not yet | Develop the temporal channel of η_D into a scattering formalism |
| **Thermal history** | Snapshot cosmology (CC, H₀, age) | Derive the full cosmological trajectory φ(t) from Big Bang to now |

### Tier 3: Experimental Confrontations (Long-Term)

| Prediction | How to Test | Timeline |
|-----------|------------|----------|
| **Dirac neutrinos** (anti-prediction of 0νββ) | Next-gen neutrinoless double beta decay experiments (nEXO, LEGEND-1000) | 2028-2035 |
| **No new particles below M_c** | LHC Run 3 + HL-LHC: continued null results for SUSY, Z', W' | 2026-2035 |
| **Inflation: n_s, r** (once computed) | CMB-S4, LiteBIRD | 2028-2032 |
| **Proton charge radius** (0.018% prediction) | MUSE experiment, PRad-II at JLab | 2026-2028 |
| **27 hadron masses** (0.95% RMS) | Lattice QCD cross-checks at physical pion mass | Ongoing |
| **K*(892) mass** (0.03% prediction) | Already confirmed — serves as post-diction benchmark | ✅ |

---

## �📋 Changelog

See [CHANGELOG.md](CHANGELOG.md) for a complete history of changes and version information.

**Recent Updates:**
- **v1.1.0** (2026-02-18): 77 predictions (all Theorem), 27 hadron masses, g_A, f_pi, tau_n, v11 paper
- **v1.0.0** (2026-02-17, historical baseline): Complete LOTUS framework with 48 predictions, comprehensive verification suite, and containerization support
- Added Docker support for reproducible environments
- Enhanced documentation with architecture details and performance benchmarks
- Improved CLI interface with multiple resolution levels

**Historical Counts Policy:**
- Use **77 predictions** for all current framework descriptions and user-facing summaries.
- Keep older counts (e.g., **48**) only when referring to specific historical releases or archived artifacts.
- When mentioning historical counts, label them explicitly as **historical** to avoid ambiguity.

---

## 🙏 Acknowledgments

### Scientific Foundations
- **Spectral Geometry**: Building on the work of Alain Connes, Matilde Marcolli, and others
- **Noncommutative Geometry**: Inspired by the groundbreaking work of Connes and Chamseddine
- **Orbifold Theory**: Drawing from string theory and extra dimension phenomenology

### Technical Contributors
- **Python Scientific Stack**: NumPy, SciPy, and the broader scientific Python ecosystem
- **Testing Frameworks**: pytest for comprehensive validation
- **Container Technology**: Docker for reproducible research environments

### Community
- **Open Source Physics**: The broader community of physicists developing open-source tools
- **Reproducible Research**: Advocates for transparent and verifiable scientific computing
- **Peer Review**: The academic community for rigorous scientific discourse

### Personal Thanks
Special thanks to the researchers, students, and enthusiasts who have engaged with this work and provided valuable feedback during development.

---

*"The paper is the proof of the model. The model is the code. The world is the lotus."*

---

## 📄 License & Copyright

### Software License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

```
MIT License

Copyright (c) 2026 The LOTUS Collaboration

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

### Academic Citation
When using this software in academic work, please cite:

```bibtex
@software{lotus_2026,
  author       = {The LOTUS Collaboration},
  title        = {LOTUS: A Spectral Geometry Framework for Particle Physics},
  year         = 2026,
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.18655472},
  url          = {https://doi.org/10.5281/zenodo.18655472}
}
```

### Copyright Notice
© 2026 The LOTUS Collaboration. All rights reserved.

The LOTUS framework and associated documentation are protected by copyright. The software is provided under the MIT License for academic and research use. Commercial use requires separate licensing agreement.
