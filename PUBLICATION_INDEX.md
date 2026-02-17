# Publication Materials Index

This document lists all files in the `public-release` folder and their purpose.

## Paper

| File | Description |
|------|--------------|
| `paper/The_Resolved_Chord_v10.tex` | **v10: The Theorem of Everything.** 43 predictions, all Theorem level, zero free parameters. Includes alpha Theorem (APS lag), alpha_s (ghost splitting), complete CKM, neutrino tunneling, Starobinsky inflation, QG dissolution, W mass prediction, 8 figures. |

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
| `supplements/Supplement_X_SpectralAction.tex` | **EXPANDED:** The Math-to-Physics Map — roadmap + 7 complete derivation chains from Tr(f(D^2)) to all predictions (proton, alpha, Higgs VEV, alpha_s, gravity, identity chain, CC) |
| `supplements/Supplement_XI_DerivationStatus.tex` | Complete derivation status table |
| `supplements/Supplement_XII_CompanionGuide.tex` | **Topic-based reference:** 6 chapters — (1) Hurricane Database, (2) Equation Index, (3) LOTUS Field Guide, (4) All Predictions / Falsification Battery, (5) Derivation Status Firewall, (6) Strange Castles (7 solved puzzles) |

## Math Papers (Standalone)

| File | Description |
|------|--------------|
| `math-papers/N1_Bridge_Theorem.tex` | **EXPANDED:** Equivariant spectral decomposition with coefficient one — universal tool for orbifold spectral theory. 5 applications (spectral action, heat kernel, Casimir, APS index, crystals), counterexample, obstructions, K-theory connection. |
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
| `verification/alpha_from_spectral_geometry.py` | 1/α (one-loop baseline) |
| `verification/alpha_lag_proof.py` | **NEW:** 1/α Theorem via APS lag correction η·λ₁/p = 10/27 |
| `verification/alpha_s_theorem.py` | **NEW:** α_s(M_Z) = 0.1187 from ghost splitting d₁ = 6 |
| `verification/alpha_derivation_chain.py` | **NEW:** Complete α derivation chain audit |
| `verification/ckm_complete.py` | **NEW:** Full CKM matrix from spectral data (9 elements, 0.00–2%) |
| `verification/vev_stiffness_proof.py` | **NEW:** VEV from bulk stiffness theorem (v/m_p = 2/α - 35/3, 0.005%) |
| `verification/fold_potential_paper.py` | **NEW:** V(φ) Mexican hat from P14–P16 (paper-accurate) |
| `verification/vev_overlap_paper.py` | **NEW:** v/m_p overlap interpretation (paper-accurate) |
| `verification/black_holes_lotus.py` | **NEW:** Black holes as the closing lotus |
| `verification/quantum_gravity_lotus.py` | **NEW:** QG dissolution: graviton=KK, UV finite, topology protected, perturbative, BH bounce. ToE 10-point checklist |
| `verification/sm_completeness_audit.py` | **NEW:** Full SM completeness audit (43 predictions) |
| `verification/vev_overlap.py` | v/m_p |
| `verification/higgs_quartic.py` | λ_H |
| `verification/downtype_spectral_ordering.py` | Quark masses |
| `verification/quark_kk_closure.py` | Quark mass closure |
| ... | (see README.md for full list) |

## Run All Verification

Run `python run_verification.py` to execute all core verification scripts. Pytest also runs them: `python -m pytest falsification/test_verification_scripts.py -v`.

## Falsification Tests

Run `python -m pytest falsification/ -v` for the full adversarial test suite (includes verification script tests).

## Root Dependencies

- `leng_test_utils.py` — Shared constants and formulas for tests
- `alpha_s_constraint.py` — Strong coupling derivation (test dependency)
- `dictionary_spec.py` — Geometry-to-physics mapping (test dependency)

## Conversion & Setup Tools

| File | Purpose |
|------|---------|
| `run_verification.py` | Run all core verification scripts. Use `--list` to list, `-q` for quiet. |
| `run_verification.bat` | Double-click to run verification; writes `verification_results.txt` |
| `run_falsification.py` | Run falsification test suite; writes `falsification_results.txt` |
| `run_falsification.bat` | Double-click to run falsification; writes `falsification_results.txt` |
| `convert_tex_pdf_md.py` | Generate .md and .pdf from .tex. Use `--sync` to pull from sibling folders, `--force` to regenerate all. |
| `check_sync.py` | Compare tex/pdf/md modification times; report out-of-sync files. |
| `requirements.txt` | Python deps: `pip install -r requirements.txt` |
| `GETTING_STARTED.md` | **Start here** — requirements, layout, quick start |
| `TROUBLESHOOTING.md` | Pandoc, pdflatex, PATH issues |
| `test_pandoc_install.bat` | Check if pandoc is in PATH |
| `install_pandoc_simple.bat` | Install pandoc (winget/Chocolatey) |
| `install_pandoc_download.bat` | Direct download (when winget/choco fail) |
| `find_pandoc.bat` | Search for pandoc on system |
