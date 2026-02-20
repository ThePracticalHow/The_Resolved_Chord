"""
Hadron PDG reference values (2024) in GeV.

Measured masses for cross-checking the Lotus Song predictions.
These values are QUARANTINED from the prediction engine.
"""

# All masses in GeV
HADRON_PDG = {
    # Pseudoscalar mesons (plucked strings)
    'pi':         0.13957,
    'K':          0.49368,
    'eta':        0.54786,
    'eta_prime':  0.95778,
    # Vector mesons (bowed strings)
    'rho':        0.77526,
    'omega':      0.78266,
    'K_star':     0.89167,
    'phi':        1.01946,
    # Baryon decuplet + ground state (drums)
    'proton':     0.93827,
    'neutron':    0.93957,
    'Delta':      1.23200,
    'Sigma_star': 1.38400,
    'Xi_star':    1.53180,
    'Omega':      1.67250,
    # Baryon octet — strangeness ladder
    'Lambda':     1.11568,
    'Sigma':      1.19264,
    'Xi':         1.31800,
    # Charm mesons (half charm-loop)
    'D':          1.86723,
    'D_star':     2.00685,
    'Ds':         1.96835,
    # Bottom mesons (eigenvalue + Koide)
    'B':          5.27934,
    'Bs':         5.36688,
    'Bc':         6.27447,
    # Heavy quarkonia (organ pipes)
    'J_psi':      3.09690,
    'psi_2S':     3.68610,
    'Upsilon':    9.46030,
    'Upsilon_2S': 10.02326,
}
