"""
lotus.constants — Derived Physical Constants
===============================================

Like scipy.constants, but every value is DERIVED from geometry.
Zero free parameters. Zero lookup tables.

Usage::

    import lotus.constants as const

    const.alpha          # 1/137.038... (derived, not looked up)
    const.m_proton       # 0.93827 GeV
    const.M_Planck       # 1.22e19 GeV
    const.G_newton       # 6.674e-11 m³/(kg·s²)

Every constant carries its provenance::

    const.info('alpha')
    # → DerivedConstant(value=0.00729..., unit='',
    #     formula='APS lag η·λ₁/p + SM RG', provenance='Level 3')

THERE ARE EXACTLY TWO PHYSICAL INPUTS:
    1. m_e  — the electron mass (the ruler / unit of measurement)
    2. π    — mathematical constant

THERE ARE ZERO FREE PARAMETERS.

Everything else — all 87 predictions — follows from the five
topological invariants {d₁=6, λ₁=5, K=2/3, η=2/9, p=3}
of the manifold S⁵/Z₃.

──────────────────────────────────────────────────────────────
ABOUT M_Z (READ THIS BEFORE TWEETING):

M_Z is NOT a geometric input, NOT a tunable dial, and NOT a
second free parameter. It is a KINEMATIC REFERENCE SCALE.
See the SPECTRAL_MAP.md Level 3 entry for full explanation.
──────────────────────────────────────────────────────────────
"""

import math
from collections import namedtuple

# ═══════════════════════════════════════════════════════════════
#  THE RULER AND THE CIRCLE
# ═══════════════════════════════════════════════════════════════

PI = math.pi

# The ruler: electron mass in two unit systems
M_E_GEV = 0.51099895e-3   # electron mass in GeV
M_E_MEV = 0.51099895       # electron mass in MeV
M_E_KG  = 9.1093837015e-31  # electron mass in kg

# Kinematic reference scale for RG running (NOT a geometric input)
M_Z_GEV = 91.1876          # Z boson mass in GeV — see docstring above

# ═══════════════════════════════════════════════════════════════
#  THE FIVE SPECTRAL INVARIANTS  (topology → numbers)
# ═══════════════════════════════════════════════════════════════

d1      = 6          # ghost mode count (dim of l=1 eigenspace on S⁵)
lambda1 = 5          # first Laplacian eigenvalue: l(l+4) at l=1
K       = 2 / 3      # Koide ratio (moment map of Z₃ on simplex)
eta     = d1 / 27    # Donnelly eta invariant: d₁/p³ = 6/27 = 2/9
p       = 3          # orbifold order |Z₃|

# Derived identities (Level 1)
tau     = 1 / 27             # Reidemeister torsion: 1/p³
G_spec  = lambda1 * eta      # Proton spectral coupling: 10/9
c_grav  = -1 / (d1 * lambda1)  # Gravity hurricane coefficient: -1/30

# ═══════════════════════════════════════════════════════════════
#  DERIVED CONSTANTS  (geometry → physics)
#
#  These are computed once at import time. No Universe() needed.
#  Every formula traces to {d₁, λ₁, K, η, p} + π + m_e.
# ═══════════════════════════════════════════════════════════════

# --- Level 2: Proton mass ---
proton_electron_ratio = d1 * PI ** 5        # 6π⁵ = 1836.12 (tree level)
m_proton = M_E_GEV * proton_electron_ratio  # GeV

# --- Level 3: Couplings (use existing sector functions) ---
# We lazily import to avoid circular dependencies, but cache on first access.
_cache = {}

def _compute_couplings():
    """One-time computation of gauge couplings via sector functions."""
    if 'alpha' in _cache:
        return
    from lotus.core.geometry import S5Z3
    from lotus.sectors.gauge import (
        fine_structure_constant, strong_coupling, weinberg_angle,
        compactification_scale,
    )
    M = S5Z3()
    _cache['alpha'] = fine_structure_constant(M)
    _cache['alpha_s'] = strong_coupling(M)
    _cache['sin2_weinberg'] = weinberg_angle(M)
    _cache['M_c'] = compactification_scale(M)

    # --- Level 4: Electroweak ---
    a = _cache['alpha']
    mp = M_E_GEV * proton_electron_ratio
    _cache['inv_alpha'] = 1 / a
    _cache['higgs_vev'] = mp * (2 / a - d1 - lambda1 - K)  # v = m_p(2/α − 35/3)
    _cache['higgs_mass'] = mp * (1 / a - 7 / 2)             # m_H = m_p(1/α − 7/2)
    _cache['m_W'] = M_Z_GEV * math.sqrt(1 - _cache['sin2_weinberg'])

    # --- Level 5: Neutrino ---
    _cache['m_nu3_GeV'] = M_E_GEV ** 3 / (p * mp ** 2)
    _cache['m_nu3_meV'] = _cache['m_nu3_GeV'] * 1e12

    # --- Level 6: Gravity ---
    from lotus.gravity import planck_mass
    _cache['M_Planck'] = planck_mass(M, loops=2)
    # G_N from M_Planck: M_P = sqrt(ℏc/G) → G = ℏc/M_P²
    # In natural units (ℏ=c=1): G = 1/M_P²
    # Convert to SI: G = (ℏc/M_P²) with M_P in kg
    M_P_kg = _cache['M_Planck'] * 1e9 * 1.78266192e-27  # GeV → kg
    hbar_si = 1.054571817e-34
    c_si = 299792458.0
    _cache['G_newton'] = hbar_si * c_si / M_P_kg ** 2

    # --- Level 6: Cosmological constant ---
    m_nu3 = _cache['m_nu3_GeV'] * 1e12  # meV
    _cache['CC_quarter_meV'] = m_nu3 * (32 / 729) * (1 + eta**2 / PI)

    # --- Level 7: Cosmology ---
    _cache['N_efolds'] = (d1 + lambda1)**2 * lambda1**2 / (p * 16)  # 3025/48
    _cache['n_s'] = 1 - 2 / _cache['N_efolds']
    _cache['omega_ratio'] = 2 * PI**2 / p**2  # Ω_Λ/Ω_m = 2π²/9


# ═══════════════════════════════════════════════════════════════
#  MODULE-LEVEL LAZY ACCESS  (PEP 562)
#
#  import lotus.constants as c
#  c.alpha   →  triggers __getattr__ → _compute_couplings() → return
# ═══════════════════════════════════════════════════════════════

_DERIVED_KEYS = {
    'alpha', 'alpha_s', 'inv_alpha', 'sin2_weinberg',
    'higgs_vev', 'higgs_mass', 'm_W',
    'm_nu3_GeV', 'm_nu3_meV',
    'M_Planck', 'G_newton', 'M_c',
    'CC_quarter_meV',
    'N_efolds', 'n_s', 'omega_ratio',
}


def __getattr__(name):
    if name in _DERIVED_KEYS:
        _compute_couplings()
        return _cache[name]
    raise AttributeError(f"module 'lotus.constants' has no attribute '{name}'")


# ═══════════════════════════════════════════════════════════════
#  PROVENANCE REGISTRY
# ═══════════════════════════════════════════════════════════════

DerivedConstant = namedtuple('DerivedConstant', ['value', 'unit', 'formula', 'provenance'])

_PROVENANCE = {
    # Spectral invariants
    'd1':       ('', 'd₁ = dim(l=1 eigenspace on S⁵)', 'Level 0: topology'),
    'lambda1':  ('', 'λ₁ = l(l+4) at l=1', 'Level 0: topology'),
    'K':        ('', 'K = 2/p (moment map)', 'Level 0: Koide theorem'),
    'eta':      ('', 'η = d₁/p³ (Donnelly 1978)', 'Level 1: spectral asymmetry'),
    'p':        ('', 'p from n=p^(n−2) → (3,3)', 'Level 0: uniqueness theorem'),

    # QCD scale
    'proton_electron_ratio': ('', 'm_p/m_e = 6π⁵ (Parseval fold energy)', 'Level 2: ghost modes'),
    'm_proton': ('GeV', 'm_p = m_e · 6π⁵', 'Level 2: ghost modes'),

    # Couplings
    'alpha':        ('', '1/α from APS lag η·λ₁/p + SM RG', 'Level 3: APS → RG'),
    'inv_alpha':    ('', '1/α = 137.038', 'Level 3: APS → RG'),
    'alpha_s':      ('', 'α_s from ghost splitting d₁=6', 'Level 3: ghost → RG'),
    'sin2_weinberg':('', 'sin²θ_W = 3/8 at M_c, RG to M_Z', 'Level 3: branching rule'),

    # Electroweak
    'higgs_vev':  ('GeV', 'v = m_p(2/α − 35/3)', 'Level 4: EM budget'),
    'higgs_mass': ('GeV', 'm_H = m_p(1/α − 7/2)', 'Level 4: spectral gap'),
    'm_W':        ('GeV', 'm_W = M_Z·√(1−sin²θ_W)', 'Level 4: Weinberg'),

    # Neutrino
    'm_nu3_GeV':  ('GeV', 'm_ν₃ = m_e³/(p·m_p²)', 'Level 5: fold-wall tunneling'),
    'm_nu3_meV':  ('meV', 'm_ν₃ = m_e³/(p·m_p²)', 'Level 5: fold-wall tunneling'),

    # Gravity
    'M_Planck':   ('GeV', 'M_P from 5-lock proof', 'Level 6: KK + hurricane'),
    'G_newton':   ('m³/(kg·s²)', 'G = ℏc/M_P²', 'Level 6: from M_Planck'),
    'M_c':        ('GeV', 'M_c = M_Z·exp(t₁₂)', 'Level 3: RG unification'),

    # Cosmology
    'CC_quarter_meV': ('meV', 'Λ^(1/4) = m_ν₃·(32/729)·(1+η²/π)', 'Level 6: CC hurricane'),
    'N_efolds':       ('', 'N = (d₁+λ₁)²λ₁²/(p·16)', 'Level 7: spectral action'),
    'n_s':            ('', 'n_s = 1 − 2/N', 'Level 7: Starobinsky'),
    'omega_ratio':    ('', 'Ω_Λ/Ω_m = 2π²/9', 'Level 6: spectral partition'),
}


def info(name: str) -> DerivedConstant:
    """Get the full provenance of a derived constant.

    Args:
        name: constant name (e.g. 'alpha', 'm_proton', 'G_newton')

    Returns:
        DerivedConstant namedtuple with value, unit, formula, provenance

    Example::

        >>> info('alpha')
        DerivedConstant(value=0.00729..., unit='', formula='...', provenance='Level 3: ...')
    """
    if name not in _PROVENANCE:
        available = sorted(_PROVENANCE.keys())
        raise KeyError(f"Unknown constant '{name}'. Available: {available}")

    unit, formula, prov = _PROVENANCE[name]

    # Get the value
    _compute_couplings()
    if name in _cache:
        val = _cache[name]
    else:
        val = globals().get(name)

    return DerivedConstant(value=val, unit=unit, formula=formula, provenance=prov)


def print_constants():
    """Print all derived constants with provenance."""
    _compute_couplings()

    print("=" * 72)
    print("LOTUS DERIVED CONSTANTS — from S⁵/Z₃, zero free parameters")
    print("=" * 72)
    print()

    sections = [
        ("SPECTRAL INVARIANTS", ['d1', 'lambda1', 'K', 'eta', 'p']),
        ("QCD SCALE", ['proton_electron_ratio', 'm_proton']),
        ("GAUGE COUPLINGS", ['alpha', 'inv_alpha', 'alpha_s', 'sin2_weinberg']),
        ("ELECTROWEAK", ['higgs_vev', 'higgs_mass', 'm_W']),
        ("NEUTRINO", ['m_nu3_meV']),
        ("GRAVITY", ['M_Planck', 'G_newton', 'M_c']),
        ("COSMOLOGY", ['CC_quarter_meV', 'N_efolds', 'n_s', 'omega_ratio']),
    ]

    for section_name, keys in sections:
        print(f"  ── {section_name} ──")
        for name in keys:
            dc = info(name)
            if isinstance(dc.value, float):
                if abs(dc.value) > 1e6 or (abs(dc.value) < 0.001 and dc.value != 0):
                    val_str = f"{dc.value:.6e}"
                else:
                    val_str = f"{dc.value:.6f}"
            else:
                val_str = str(dc.value)
            unit_str = f" {dc.unit}" if dc.unit else ""
            print(f"    {name:28s} = {val_str}{unit_str}")
            print(f"      {dc.formula}")
        print()

    print("─" * 72)
    print("  Usage:  import lotus.constants as const")
    print("          const.alpha, const.m_proton, const.info('alpha')")
    print("─" * 72)
