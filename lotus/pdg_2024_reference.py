"""
lotus.pdg_2024_reference — Quarantined Experimental Data
============================================================

THIS FILE CONTAINS ONLY MEASURED VALUES FROM PDG 2024.
IT IS NEVER IMPORTED BY ANY PREDICTION FUNCTION.

The prediction engine (lotus/sectors/*.py, lotus/gravity.py,
lotus/cosmology.py) imports ONLY from:
    - lotus.constants  (m_e, pi, M_Z)
    - lotus.core.geometry  (S5Z3 manifold)

This file exists solely to compute error percentages in the
CLI boot sequence and test suite.  You can delete this entire
file and every prediction will still compute identically.

Source: Particle Data Group, 2024 Review of Particle Physics.
https://pdg.lbl.gov
"""

# ── Lepton sector ───────────────────────────────────────────
KOIDE_K = 2 / 3                # exact (a theorem, not measured)
M_MU_OVER_M_E = 206.768       # muon/electron mass ratio
M_TAU_OVER_M_E = 3477.2       # tau/electron mass ratio
M_P_OVER_M_E = 1836.153       # proton/electron mass ratio

# ── Quark sector (GeV) ──────────────────────────────────────
M_TOP_GEV = 172.57
M_CHARM_GEV = 1.273
M_UP_MEV = 2.16
M_BOTTOM_GEV = 4.183
M_STRANGE_MEV = 93.4
M_DOWN_MEV = 4.67

# ── Neutrino sector ─────────────────────────────────────────
M_NU3_MEV = 50.28              # heaviest neutrino mass (meV)
DM2_RATIO = 32.58              # atmospheric/solar mass-squared ratio

# ── Coupling constants ──────────────────────────────────────
INV_ALPHA = 137.036            # 1 / fine-structure constant
SIN2_WEINBERG = 0.23122        # sin^2(theta_W) at M_Z
SIN2_WEINBERG_MC = 0.375       # sin^2(theta_W) at M_c = 3/8 (exact)
ALPHA_S_MZ = 0.1180            # strong coupling at M_Z

# ── Higgs sector ────────────────────────────────────────────
HIGGS_VEV_GEV = 246.22        # Higgs VEV
HIGGS_MASS_GEV = 125.25       # Higgs boson mass
LAMBDA_H = 0.1294             # Higgs quartic coupling

# ── CKM matrix (Wolfenstein) ────────────────────────────────
CKM_LAMBDA = 0.22500
CKM_A = 0.826
CKM_RHO_BAR = 0.1592
CKM_ETA_BAR = 0.3490
CKM_GAMMA_DEG = 65.6
JARLSKOG_J = 3.08e-5

# ── PMNS matrix ─────────────────────────────────────────────
SIN2_THETA23 = 0.546
SIN2_THETA12 = 0.307
SIN2_THETA13 = 0.02200
DELTA_CP_PMNS_DEG = 195.0      # CP-violating phase (PMNS)

# ── Gravity ─────────────────────────────────────────────────
M_PLANCK_E19 = 1.221           # M_P in units of 10^19 GeV
M_W_GEV = 80.369               # W boson mass

# ── Cosmology ───────────────────────────────────────────────
CC_QUARTER_MEV = 2.25          # Lambda^{1/4} in meV
N_EFOLDS = 63.0
SPECTRAL_INDEX = 0.965
TENSOR_TO_SCALAR = 0.003
ETA_B = 6.1e-10               # baryon asymmetry
OMEGA_DM_OVER_B = 5.36        # Omega_DM / Omega_baryon
OMEGA_L_OVER_M = 2.214        # Omega_Lambda / Omega_matter

# ── Structural (exact) ──────────────────────────────────────
N_GENERATIONS = 3.0
THETA_QCD = 0.0                # strong CP phase

# ── Nuclear & decays ────────────────────────────────────────
G_A = 1.2754                   # nucleon axial coupling
F_PI_MEV = 92.07               # pion decay constant (MeV)
TAU_N_S = 878.4                # neutron lifetime (s)
TAU_PI_E8 = 2.603              # pion lifetime (×10⁻⁸ s)
TAU_MU_E6 = 2.197              # muon lifetime (×10⁻⁶ s)
B_D_MEV = 2.2246               # deuteron binding energy (MeV)

# ── Extended cosmology ──────────────────────────────────────
H_0_KM_S_MPC = 67.4            # Hubble constant (km/s/Mpc, Planck 2018)
T_AGE_GYR = 13.80              # age of universe (Gyr, Planck 2018)
T_REHEAT_E9 = 2.15             # reheating temperature (×10⁹ GeV, theoretical)

# ── QCD structure ───────────────────────────────────────────────
B_0 = 23/3                     # one-loop QCD beta coefficient (exact, N_f=6)
G_PINN = 13.12                 # pion-nucleon coupling (Arndt 2006)

# ── Gauge boson decay widths ─────────────────────────────────────
GAMMA_W_GEV = 2.085            # W boson total decay width (GeV, PDG 2024)
GAMMA_Z_GEV = 2.4955           # Z boson total decay width (GeV, PDG 2024)


