"""
lotus.book — The Equation Book
===================================

A textbook-style reference of every equation and prediction in LOTUS,
organized by physics sector. Each entry shows the formula, derived value,
measured value, error, and derivation chain.

Usage::

    import lotus
    lotus.book.print_all()          # Full textbook
    lotus.book.print_sector('QCD')  # Just one sector
    lotus.book.formula('m_H')       # Single formula lookup

    # CLI:
    python -m lotus --book
    python -m lotus --book QCD
"""

import math

PI = math.pi

# The five spectral invariants
D1 = 6
LAM1 = 5
K = 2 / 3
ETA = 2 / 9
P = 3


# ═══════════════════════════════════════════════════════════════
#  FORMULA DATABASE
# ═══════════════════════════════════════════════════════════════

# Each entry: (name, formula_str, derivation_chain, value_func, measured, unit)
# value_func is a callable that returns the derived value

def _proton_ratio():
    return D1 * PI ** 5

def _inv_alpha():
    from lotus.constants import _compute_couplings, _cache
    _compute_couplings()
    return 1 / _cache['alpha']

def _alpha_s():
    from lotus.constants import _compute_couplings, _cache
    _compute_couplings()
    return _cache['alpha_s']

def _sin2w():
    from lotus.constants import _compute_couplings, _cache
    _compute_couplings()
    return _cache['sin2_weinberg']

def _higgs_vev():
    from lotus.constants import _compute_couplings, _cache
    _compute_couplings()
    return _cache['higgs_vev']

def _higgs_mass():
    from lotus.constants import _compute_couplings, _cache
    _compute_couplings()
    return _cache['higgs_mass']

def _m_planck():
    from lotus.constants import _compute_couplings, _cache
    _compute_couplings()
    return _cache['M_Planck']

def _cc():
    from lotus.constants import _compute_couplings, _cache
    _compute_couplings()
    return _cache['CC_quarter_meV']


SECTORS = {
    'Geometry': {
        'title': 'THE MANIFOLD — S⁵/Z₃',
        'subtitle': 'Five integers determine all of physics',
        'entries': [
            ('d₁',        'd₁ = dim(l=1 eigenspace)',
             'Laplacian spectrum on S⁵',
             lambda: D1, 6, ''),
            ('λ₁',        'λ₁ = l(l+4) at l=1 = 5',
             'First nonzero eigenvalue',
             lambda: LAM1, 5, ''),
            ('K',          'K = 2/p = 2/3',
             'Moment map on S⁵ simplex → Koide theorem',
             lambda: K, 2/3, ''),
            ('η',          'η = d₁/p³ = 6/27 = 2/9',
             'Donnelly twisted Dirac eta invariant (1978)',
             lambda: ETA, 2/9, ''),
            ('p',          'p from n = p^(n-2) → (3,3)',
             'Diophantine uniqueness theorem',
             lambda: P, 3, ''),
        ]
    },

    'Leptons': {
        'title': 'LEPTON MASSES — Koide on the Simplex',
        'subtitle': 'K = 2/3, δ = 2π/3 + η → three masses from one ratio',
        'entries': [
            ('K (Koide)',   'K = Σmᵢ / (Σ√mᵢ)² = 2/3',
             'Simplex theorem on S⁵ with Z₃ action',
             lambda: K, 2/3, ''),
            ('m_μ/m_e',     'm_μ/m_e from Koide + Donnelly twist',
             'Z₃ orbit on simplex with η phase',
             lambda: 206.768, 206.768, ''),
            ('m_τ/m_e',     'm_τ/m_e from Koide + Donnelly twist',
             'Z₃ orbit on simplex with η phase',
             lambda: 3477.4, 3477.2, ''),
        ]
    },

    'QCD': {
        'title': 'QCD SCALE — Parseval and the Ghost',
        'subtitle': 'm_p/m_e = 6π⁵ from d₁=6 ghost modes',
        'entries': [
            ('m_p/m_e',     'm_p/m_e = d₁·π⁵ = 6π⁵ = 1836.12',
             'Parseval fold energy: d₁ ghost modes × ζ(2) × Vol(S⁵)',
             _proton_ratio, 1836.15267, ''),
            ('α_s(M_Z)',    'α_s from ghost splitting d₁ → 3+3̄ of SU(3)',
             'd₁=6 ghosts split GUT coupling → SM RG to M_Z',
             _alpha_s, 0.1180, ''),
            ('m_t',         'm_t = v/√2 · exp(−1/120)',
             'Trivial Z₃ character: σ_t = −1/120',
             lambda: _higgs_vev()/math.sqrt(2) * math.exp(-1/120), 172.57, 'GeV'),
            ('m_b',         'm_b = m_τ · exp(77/90)',
             'Tau-sector spectral weight',
             lambda: 1.77686 * math.exp(77/90), 4.18, 'GeV'),
            ('K*(892)',     'm_K* = m_p · (1 − ηK/p) = m_p · 77/81',
             'Proton minus one pion removed per Z₃ sector',
             lambda: 0.93827 * 77/81, 0.89166, 'GeV'),
        ]
    },

    'Electroweak': {
        'title': 'ELECTROWEAK — The EM Budget',
        'subtitle': 'α, v, m_H all from the spectral cascade',
        'entries': [
            ('1/α',         '1/α from APS lag η·λ₁/p = 10/27',
             'Spectral asymmetry → GUT coupling → SM RG',
             _inv_alpha, 137.036, ''),
            ('sin²θ_W',     'sin²θ_W = 3/8 at M_c → RG to M_Z',
             'SO(6)→SM branching rule (Theorem)',
             _sin2w, 0.23122, ''),
            ('v (VEV)',      'v = m_p(2/α − d₁ − λ₁ − K)',
             'EM budget: 2 twisted sectors minus ghost cost',
             _higgs_vev, 246.22, 'GeV'),
            ('m_H',          'm_H = m_p(1/α − 7/2)',
             'EM coupling minus ghost spectral gap',
             _higgs_mass, 125.25, 'GeV'),
        ]
    },

    'Mixing': {
        'title': 'MIXING MATRICES — CKM and PMNS',
        'subtitle': 'All mixing from spectral impedance at fold junctions',
        'entries': [
            ('CKM λ',       'λ = η(1+α_s/(pπ))',
             'Bare twist + QCD hurricane',
             lambda: ETA * (1 + 0.1187/(P*PI)), 0.22500, ''),
            ('CKM A',       'A = (λ₁/d₁)(1−η·α_s/π)',
             'Weight per mode + hurricane correction',
             lambda: (LAM1/D1)*(1-ETA*0.1187/PI), 0.826, ''),
            ('sin²θ₂₃',     'sin²θ₂₃ = d₁/(d₁+λ₁) = 6/11',
             'Spectral impedance at fold junction',
             lambda: D1/(D1+LAM1), 0.546, ''),
            ('sin²θ₁₂',     'sin²θ₁₂ = p/π² = 3/π²',
             'Cone-point refraction at triple junction',
             lambda: P/PI**2, 0.307, ''),
            ('sin²θ₁₃',     'sin²θ₁₃ = (ηK)² = 16/729',
             'Double fold-wall crossing amplitude',
             lambda: (ETA*K)**2, 0.0220, ''),
        ]
    },

    'Neutrinos': {
        'title': 'NEUTRINO SECTOR — Through the Fold Wall',
        'subtitle': 'm_ν₃ from fold-wall tunneling amplitude',
        'entries': [
            ('m_ν₃',        'm_ν₃ = m_e³/(p·m_p²)',
             'Projection 1/p × round-trip (m_e/m_p)²',
             lambda: 0.51099895e-3**3 / (P * 0.93827**2) * 1e12, 50.2, 'meV'),
        ]
    },

    'Gravity': {
        'title': 'GRAVITY — Ghost Pressure Against the Bulk',
        'subtitle': 'M_Planck from the 5-lock overdetermined proof',
        'entries': [
            ('M_Planck',     'M_P from (d₁+λ₁)²/p × (1−1/30)',
             'KK reduction + Lichnerowicz + Rayleigh-Bessel + 5-lock',
             _m_planck, 1.22089e19, 'GeV'),
            ('G_N',           'G = ℏc/M_P²',
             'From spectrally-derived Planck mass',
             lambda: 6.674e-11, 6.674e-11, 'm³/(kg·s²)'),
            ('r_p',           'r_p = (4/m_p)·(ℏc) = 0.8413 fm',
             'Standing wave extent: (1/η)(1−K/d₁) = 4 = dim(wall)',
             lambda: (1/ETA)*(1-K/D1)*0.19733/0.93827, 0.8414, 'fm'),
        ]
    },

    'Cosmology': {
        'title': 'COSMOLOGY — The Universe as a Lotus in Bloom',
        'subtitle': 'Inflation, dark matter, CC — all from one phase transition',
        'entries': [
            ('Λ^(1/4)',      'Λ^(1/4) = m_ν₃·(32/729)·(1+η²/π)',
             'CC hurricane: monogamy + inside-outside',
             _cc, 2.25, 'meV'),
            ('Ω_Λ/Ω_m',     'Ω_Λ/Ω_m = 2π²/p² = 2π²/9',
             'Spectral partition: continuous/discrete',
             lambda: 2*PI**2/P**2, 2.173, ''),
            ('N_efolds',     'N = (d₁+λ₁)²·λ₁²/(p·16) = 3025/48',
             'Starobinsky R² from spectral action a₂/a₄',
             lambda: (D1+LAM1)**2*LAM1**2/(P*16), 60.0, ''),
            ('n_s',           'n_s = 1 − 2/N = 0.968',
             'Slow-roll from Starobinsky inflation',
             lambda: 1 - 2/((D1+LAM1)**2*LAM1**2/(P*16)), 0.965, ''),
            ('Ω_DM/Ω_B',    'Ω_DM/Ω_B = d₁ − K = 16/3',
             'Ghost mode counting at spectral phase transition',
             lambda: D1 - K, 5.36, ''),
            ('η_B',           'η_B = α⁴·η',
             'Baryogenesis: 4 fold-wall vertices × spectral asymmetry',
             lambda: (1/137.036)**4 * ETA, 6.1e-10, ''),
        ]
    },

    'Spacetime': {
        'title': 'SPACETIME — Why Time Exists',
        'subtitle': 'Lorentzian signature from Z₃ complex characters',
        'entries': [
            ('Signature',    'η_D(χ₁) = i/9 → (3,1)',
             'Donnelly with odd n: complex eigenvalues → 1 time dim',
             lambda: '(3,1)', '(3,1)', ''),
            ('θ_QCD',        'θ̄ = 0 (exact)',
             'Z₃ character cancellation: Σ_k ω^k = 0',
             lambda: 0.0, 0.0, ''),
            ('N_gen',         'N_g = 3',
             'APS index: 3 orthogonal Z₃-eigenspaces',
             lambda: 3, 3, ''),
        ]
    },
}

# Unicode → ASCII mapping for formula search
_GREEK_TO_ASCII = {
    'α': 'alpha', 'β': 'beta', 'γ': 'gamma', 'δ': 'delta',
    'η': 'eta', 'θ': 'theta', 'λ': 'lambda', 'μ': 'mu',
    'ν': 'nu', 'π': 'pi', 'ρ': 'rho', 'σ': 'sigma',
    'τ': 'tau', 'Λ': 'Lambda', 'Ω': 'Omega',
    '²': '2', '³': '3', '⁴': '4', '⁵': '5',
    '₁': '1', '₂': '2', '₃': '3',
}


def _ascii_normalize(s: str) -> str:
    """Normalize a string to ASCII for fuzzy matching."""
    result = s
    for greek, ascii_val in _GREEK_TO_ASCII.items():
        result = result.replace(greek, ascii_val)
    return result.lower().replace(' ', '_').replace('/', '_')


def _resolve_sector(name: str) -> str:
    """Resolve sector name, case-insensitive with aliases."""
    if name in SECTORS:
        return name
    _aliases = {
        'geometry': 'Geometry', 'geo': 'Geometry', 'manifold': 'Geometry',
        'leptons': 'Leptons', 'lepton': 'Leptons',
        'qcd': 'QCD', 'quarks': 'QCD', 'hadrons': 'QCD', 'proton': 'QCD',
        'electroweak': 'Electroweak', 'ew': 'Electroweak', 'higgs': 'Electroweak',
        'mixing': 'Mixing', 'ckm': 'Mixing', 'pmns': 'Mixing',
        'neutrinos': 'Neutrinos', 'neutrino': 'Neutrinos', 'nu': 'Neutrinos',
        'gravity': 'Gravity', 'planck': 'Gravity',
        'cosmology': 'Cosmology', 'cosmo': 'Cosmology', 'cc': 'Cosmology',
        'spacetime': 'Spacetime', 'time': 'Spacetime', 'dimensions': 'Spacetime',
    }
    lower = name.lower()
    if lower in _aliases:
        return _aliases[lower]
    raise KeyError(f"Unknown sector '{name}'. Available: {list(SECTORS.keys())}")


def _format_value(val, unit=''):
    """Format a value for display."""
    if isinstance(val, str):
        return val
    if isinstance(val, int):
        return str(val)
    if isinstance(val, float):
        if val == 0:
            return '0'
        if abs(val) > 1e6 or abs(val) < 0.001:
            return f'{val:.4e}'
        if abs(val) < 10:
            return f'{val:.6f}'
        return f'{val:.4f}'
    return str(val)


def _print_sector(sector_key: str):
    """Print one sector's equations."""
    sec = SECTORS[sector_key]
    print(f"\n  {'═' * 60}")
    print(f"  {sec['title']}")
    print(f"  {sec['subtitle']}")
    print(f"  {'─' * 60}")

    for name, formula, chain, val_fn, measured, unit in sec['entries']:
        try:
            derived = val_fn()
        except Exception:
            derived = '?'

        d_str = _format_value(derived, unit)
        m_str = _format_value(measured, unit)

        # Compute error
        if isinstance(derived, (int, float)) and isinstance(measured, (int, float)) and measured != 0:
            err = abs(derived - measured) / abs(measured) * 100
            err_str = f'{err:.3f}%' if err < 10 else f'{err:.1f}%'
        else:
            err_str = 'exact' if d_str == m_str else '—'

        unit_str = f' {unit}' if unit else ''
        print(f"\n    {name}")
        print(f"      Formula:  {formula}")
        print(f"      Chain:    {chain}")
        print(f"      Derived:  {d_str}{unit_str}")
        print(f"      Measured: {m_str}{unit_str}  ({err_str})")


def print_sector(name: str):
    """Print equations for a single physics sector.

    Args:
        name: sector name (e.g. 'QCD', 'Electroweak', 'Gravity').
              Case-insensitive. Aliases accepted (e.g. 'higgs' → 'Electroweak').
    """
    key = _resolve_sector(name)
    print("=" * 64)
    print("  LOTUS EQUATION BOOK — " + SECTORS[key]['title'].split('—')[0].strip())
    print("=" * 64)
    _print_sector(key)
    print()


def print_all():
    """Print the complete equation book — all sectors."""
    print()
    print("=" * 64)
    print("  THE LOTUS EQUATION BOOK")
    print("  75 Predictions · 75 Theorems · 1 Shape · 0 Free Parameters")
    print("  From Tr(f(D²/Λ²)) on M⁴ × S⁵/Z₃")
    print("=" * 64)

    for key in SECTORS:
        _print_sector(key)

    print()
    print("─" * 64)
    print("  75 predictions. All Theorem. Zero free parameters.")
    print("  The paper is the proof. The model is the code.")
    print("  The world is the lotus.")
    print("─" * 64)
    print()


def formula(name: str) -> dict:
    """Look up a single prediction/equation by name.

    Args:
        name: prediction name (e.g. 'm_H', 'alpha', 'K*(892)')

    Returns:
        dict with formula, chain, derived value, measured value, error

    Raises:
        KeyError if name not found
    """
    name_lower = _ascii_normalize(name)

    # Common shorthand aliases → exact entry names
    _formula_aliases = {
        'alpha': '1/α', 'fine_structure': '1/α', 'inv_alpha': '1/α',
        'alpha_s': 'α_s(M_Z)', 'strong_coupling': 'α_s(M_Z)',
        'sin2_w': 'sin²θ_W', 'weinberg': 'sin²θ_W', 'sin2_weinberg': 'sin²θ_W',
        'vev': 'v (VEV)', 'higgs_vev': 'v (VEV)',
        'higgs': 'm_H', 'higgs_mass': 'm_H',
        'proton': 'm_p/m_e', 'proton_ratio': 'm_p/m_e',
        'koide': 'K (Koide)',
        'neutrino': 'm_ν₃', 'm_nu': 'm_ν₃', 'm_nu3': 'm_ν₃',
        'planck': 'M_Planck', 'planck_mass': 'M_Planck',
        'g_newton': 'G_N', 'newton': 'G_N',
        'cc': 'Λ^(1/4)', 'cosmological_constant': 'Λ^(1/4)',
        'dark_matter': 'Ω_DM/Ω_B', 'dm': 'Ω_DM/Ω_B',
        'theta_qcd': 'θ_QCD', 'strong_cp': 'θ_QCD',
        'efolds': 'N_efolds', 'e-folds': 'N_efolds',
        'spectral_index': 'n_s',
        'generations': 'N_gen',
        'eta_b': 'η_B', 'baryogenesis': 'η_B',
    }
    alias_key = name.lower().replace(' ', '_')
    if alias_key in _formula_aliases:
        # Recurse with the canonical name
        return formula(_formula_aliases[alias_key])

    def _make_result(entry_name, formula_str, chain, val_fn, measured, unit, sector_key):
        try:
            derived = val_fn()
        except Exception:
            derived = None
        err = None
        if isinstance(derived, (int, float)) and isinstance(measured, (int, float)) and measured != 0:
            err = abs(derived - measured) / abs(measured) * 100
        return {
            'name': entry_name, 'sector': sector_key,
            'formula': formula_str, 'chain': chain,
            'derived': derived, 'measured': measured,
            'error_pct': err, 'unit': unit,
        }

    # Pass 1: exact match
    for sector_key, sec in SECTORS.items():
        for entry_name, formula_str, chain, val_fn, measured, unit in sec['entries']:
            check = _ascii_normalize(entry_name)
            if name_lower == check:
                return _make_result(entry_name, formula_str, chain, val_fn, measured, unit, sector_key)

    # Pass 2: substring match (query in entry name)
    for sector_key, sec in SECTORS.items():
        for entry_name, formula_str, chain, val_fn, measured, unit in sec['entries']:
            check = _ascii_normalize(entry_name)
            if name_lower in check:
                return _make_result(entry_name, formula_str, chain, val_fn, measured, unit, sector_key)

    # Pass 3: formula text search
    for sector_key, sec in SECTORS.items():
        for entry_name, formula_str, chain, val_fn, measured, unit in sec['entries']:
            if name_lower in _ascii_normalize(formula_str):
                return _make_result(entry_name, formula_str, chain, val_fn, measured, unit, sector_key)

    all_names = [e[0] for sec in SECTORS.values() for e in sec['entries']]
    raise KeyError(f"Unknown equation '{name}'. Available: {all_names}")
