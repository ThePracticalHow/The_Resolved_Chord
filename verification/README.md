# Verification Scripts

Each script verifies a specific prediction or theorem. Run from `public-release/` root.

## Entry Points

| Script | Purpose |
|--------|---------|
| **compile_universe.py** | **Start here.** Compiles all 87 predictions in ~3 seconds |
| **run_verification.py** | Runs all core verification scripts |
| **run_verification.bat** | Double-click to run verification only |

## Scripts by Topic

### Alpha / Couplings
- `alpha_derivation_chain.py`, `alpha_from_spectral_geometry.py`, `alpha_lag_proof.py`
- `alpha_s_theorem.py`, `alpha_s_constraint.py`, `alpha_two_loop_rg.py`
- `quark_rg_full_sm.py`, `quark_piercing_rg.py`, `quark_kk_closure.py`

### Cosmological Constant
- `cc_theorem.py`, `cc_hurricane.py`, `cc_hurricane_proof.py`
- `cc_aps_proof.py`, `cc_zeta_proof.py`, `cc_monogamy_proof.py`, `cc_monogamy_cancellation.py`
- `lotus_cc_oneloop.py`

### Gravity
- `gravity_theorem_proof.py`, `gravity_fold_connection.py`, `gravity_hurricane.py`
- `quantum_gravity_lotus.py`, `black_holes_lotus.py`

### Cosmology
- `cosmic_snapshot_epoch.py`, `age_of_universe.py`, `h0_spectral.py`
- `starobinsky_theorem.py`, `baryogenesis_dm_theorem.py`, `kawasaki_index.py`

### Higgs / VEV
- `higgs_vev_spectral_action.py`, `higgs_a2_integral.py`, `higgs_quartic.py`
- `vev_overlap.py`, `vev_overlap_paper.py`, `vev_stiffness_proof.py`
- `fold_potential.py`, `fold_potential_paper.py`, `fold_wall_scalar.py`

### Hadrons / Baryons
- `neutron_lifetime.py`, `neutron_properties.py`, `proton_magnetic_moment.py`
- `proton_radius.py`, `axial_coupling_derivation.py`
- `qcd_confinement_scale.py`

### Mixing (CKM / PMNS)
- `ckm_from_geometry.py`, `ckm_complete.py`, `cabibbo_hurricane.py`
- `pmns_point_side_face.py`, `downtype_spectral_ordering.py`

### Neutrinos
- `neutrino_tunneling_theorem.py`, `theta12_solar_theorem.py`

### Spectral Action / Geometry
- `EtaInvariant.py`, `GhostModes.py` — Foundation (η = 2/9)
- `spectral_action_master.py`, `spectral_action_derivation.py`, `spectral_action_dictionary.py`
- `spectral_loop_theorem.py`, `hurricane_proof.py`, `hurricane_coefficient_search.py`
- `constraint_grammar.py`, `geometric_unification.py`

### LOTUS Song / Hadron Spectrum
- `lotus_song.py`, `lotus_song_derivation.py`, `lotus_song_eigenvalue.py`
- `lotus_song_evolving.py`, `lotus_song_extended.py`
- `remaining_seven.py`

### Fold / Dynamics
- `dirac_fold_transition.py`, `lorentzian_proof.py`
- `lotus_eom.py`, `lotus_dynamics.py`, `lotus_potential.py`
- `lotus_arrow.py`, `lotus_signature.py`, `lotus_aps_generation.py`

### Meta / Audits
- `precise_recount.py`, `theorem_promotions.py`, `theorem_everything.py`
- `sm_completeness_audit.py`, `leng_replication.py`
- `UniverseLandscape.py`, `piercing_uniqueness_test.py`
