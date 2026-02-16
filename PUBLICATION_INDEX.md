# Publication Materials Index

This document lists all files in the `resolved-chord` public release and their purpose.

## Paper

| File | Description |
|------|--------------|
| `paper/The_Resolved_Chord_v9.tex` | Main paper — 18 sections, 1615 lines. The spectral geometry of S⁵/Z₃ determines all 26 SM parameters, gravity, and the cosmological constant. |

## Supplements (Proof Chains)

| File | Description |
|------|--------------|
| `supplements/Supplement_I_Geometry.tex` | Manifold, spectral invariants, Donnelly eta, uniqueness theorem |
| `supplements/Supplement_II_LeptonSector.tex` | Koide, Yukawa phase, N=1 bridge |
| `supplements/Supplement_III_GaugeConfinement.tex` | Triple spectral exclusion, confinement, Weinberg angle |
| `supplements/Supplement_IV_BaryonSector.tex` | Proton mass, ghost modes, spectral action |
| `supplements/Supplement_V_HiggsSector.tex` | VEV, Higgs mass, quartic coupling |
| `supplements/Supplement_VI_QuarkFlavor.tex` | CKM, spectral ordering, quark masses |
| `supplements/Supplement_VII_NeutrinoSector.tex` | PMNS, neutrino masses |
| `supplements/Supplement_VIII_AdversarialDefense.tex` | Negative controls, look-elsewhere |
| `supplements/Supplement_IX_StrangeCastles.tex` | Beyond-SM, gravity, CC |
| `supplements/Supplement_X_SpectralAction.tex` | Spectral dictionary derivation |
| `supplements/Supplement_XI_DerivationStatus.tex` | Complete derivation status table |

## Math Papers (Standalone)

| File | Description |
|------|--------------|
| `math-papers/N1_Bridge_Theorem.tex` | N=1 Yukawa bridge from spectral monogamy |
| `math-papers/Eta_Ghost_Identity.tex` | η = d₁/pⁿ identity |
| `math-papers/APS_Three_Generations.tex` | Equivariant APS index → N_g = 3 |
| `math-papers/Spectral_Exclusion.tex` | Triple spectral exclusion theorem |

## Verification Scripts

Each script verifies a specific prediction. Run from `public-release/` root:

| Script | Verifies |
|--------|----------|
| `verification/EtaInvariant.py` | η = 2/9 |
| `verification/GhostModes.py` | d₁ = 6, λ₁ = 5 |
| `verification/leng_replication.py` | Koide masses |
| `verification/ghost_parseval_proof.py` | Proton mass 6π⁵ |
| `verification/gravity_theorem_proof.py` | Planck mass |
| `verification/gravity_hurricane.py` | c_grav = -1/30 |
| `verification/cc_aps_proof.py` | Cosmological constant |
| `verification/spectral_action_dictionary.py` | π² = λ₁ + α_s |
| `verification/alpha_from_spectral_geometry.py` | 1/α |
| `verification/vev_overlap.py` | v/m_p |
| `verification/higgs_quartic.py` | λ_H |
| `verification/downtype_spectral_ordering.py` | Quark masses |
| ... | (see README.md for full list) |

## Falsification Tests

Run `python -m pytest falsification/ -v` for the full adversarial test suite (186 tests).

## Root Dependencies

- `leng_test_utils.py` — Shared constants and formulas for tests
- `alpha_s_constraint.py` — Strong coupling derivation (test dependency)
- `dictionary_spec.py` — Geometry-to-physics mapping (test dependency)
