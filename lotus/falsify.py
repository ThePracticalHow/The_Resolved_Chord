"""
lotus.falsify — Falsification Targets and PDG Comparison
============================================================

Which predictions are most vulnerable to near-term experiment?
What exactly kills the theory?

    >>> from lotus.falsify import compare_pdg, falsification_targets
    >>> compare_pdg()           # full prediction vs PDG table
    >>> falsification_targets() # ranked by testability

Design principle: a theory that can't be killed isn't physics.
"""

import math

PI = math.pi


def compare_pdg(resolution: str = '2_loop') -> list:
    """Full comparison table: LOTUS prediction vs PDG measurement.

    Args:
        resolution: 'tree', '1_loop', or '2_loop' (default)

    Returns:
        list of dicts with name, predicted, measured, error_pct, unit
    """
    # Import here to avoid circular imports
    from lotus import Universe
    u = Universe(resolution=resolution)
    preds = u.predictions

    table = []
    for name, vals in preds.items():
        table.append({
            'name': name,
            'predicted': vals['predicted'],
            'measured': vals['measured'],
            'error_pct': vals['error_pct'],
            'unit': vals['unit'],
        })

    return sorted(table, key=lambda x: x['error_pct'], reverse=True)


def print_comparison(resolution: str = '2_loop'):
    """Print a formatted comparison table to stdout.

    Args:
        resolution: 'tree', '1_loop', or '2_loop'
    """
    table = compare_pdg(resolution)

    print()
    print(f"  LOTUS vs PDG  [resolution={resolution}]")
    print("  " + "=" * 70)
    print(f"  {'Name':<22} {'Predicted':>12}  {'Measured':>12}  {'Error':>8}  {'Unit'}")
    print("  " + "-" * 70)

    for row in sorted(table, key=lambda x: x['error_pct']):
        pred = row['predicted']
        meas = row['measured']
        err = row['error_pct']
        unit = row['unit']

        # Format prediction
        if abs(pred) > 100:
            p_str = f"{pred:.2f}"
        elif abs(pred) > 1:
            p_str = f"{pred:.4f}"
        elif abs(pred) > 0.001:
            p_str = f"{pred:.5f}"
        else:
            p_str = f"{pred:.3e}"

        # Format measurement
        if abs(meas) > 100:
            m_str = f"{meas:.2f}"
        elif abs(meas) > 0.001:
            m_str = f"{meas:.4f}"
        else:
            m_str = f"{meas:.3e}"

        # Format error
        if err == 0:
            e_str = "exact"
        elif err < 0.01:
            e_str = f"{err:.4f}%"
        elif err < 1:
            e_str = f"{err:.3f}%"
        else:
            e_str = f"{err:.1f}%"

        print(f"  {row['name']:<22} {p_str:>12}  {m_str:>12}  {e_str:>8}  {unit}")

    # Summary stats
    errs = [r['error_pct'] for r in table if r['error_pct'] > 0]
    rms = math.sqrt(sum(e ** 2 for e in errs) / len(errs)) if errs else 0
    median = sorted(errs)[len(errs) // 2] if errs else 0
    print("  " + "-" * 70)
    print(f"  RMS error: {rms:.3f}%  |  Median: {median:.3f}%  |  "
          f"Max: {max(errs):.1f}%  |  N = {len(table)}")
    print()


def falsification_targets() -> list:
    """Predictions most vulnerable to near-term experimental falsification.

    Returns a ranked list of kill-shots: if any of these disagree
    with experiment beyond their stated threshold, LOTUS is dead.

    Each entry has:
        - name: prediction name
        - predicted: LOTUS value
        - measured: current PDG value
        - error_pct: current discrepancy
        - kill_threshold: how much error would kill the theory (%)
        - experiment: which experiment tests this
        - timeframe: when the decisive measurement is expected
    """
    from lotus import Universe
    u = Universe()
    M = u.manifold

    targets = [
        {
            'name': 'Proton radius r_p',
            'predicted': (1/M.eta) * (1 - M.K/M.d1) * 0.19733 / u.proton_mass,
            'measured': 0.8414,
            'unit': 'fm',
            'kill_threshold': 1.0,
            'experiment': 'PRad-II at JLab, MUSE at PSI',
            'timeframe': '2025-2027',
            'why_fatal': 'r_p is overdetermined by 3 invariants (η, K, d₁). No tuning possible.',
        },
        {
            'name': 'Neutrino mass m_ν₃',
            'predicted': u.m_nu3_meV,
            'measured': 50.28,
            'unit': 'meV',
            'kill_threshold': 5.0,
            'experiment': 'KATRIN, Project 8, JUNO',
            'timeframe': '2025-2030',
            'why_fatal': 'Geometric seesaw m_ν = m_e³/(p·m_p²). No free parameter to adjust.',
        },
        {
            'name': 'Higgs mass m_H',
            'predicted': u.higgs_mass_val,
            'measured': 125.25,
            'unit': 'GeV',
            'kill_threshold': 0.5,
            'experiment': 'HL-LHC, FCC-ee',
            'timeframe': '2025-2035',
            'why_fatal': 'm_H = m_p(1/α − 7/2). Ghost spectral cost fixed by geometry.',
        },
        {
            'name': 'sin²θ₂₃ (PMNS)',
            'predicted': u.pmns['sin2_theta23'],
            'measured': 0.546,
            'unit': '',
            'kill_threshold': 2.0,
            'experiment': 'T2K, NOvA, DUNE, Hyper-K',
            'timeframe': '2025-2030',
            'why_fatal': 'sin²θ₂₃ = d₁/(d₁+λ₁) = 6/11. Exact rational prediction.',
        },
        {
            'name': 'sin²θ₁₃ (PMNS)',
            'predicted': u.pmns['sin2_theta13'],
            'measured': 0.02200,
            'unit': '',
            'kill_threshold': 3.0,
            'experiment': 'JUNO, Daya Bay successor',
            'timeframe': '2025-2030',
            'why_fatal': 'sin²θ₁₃ = (ηK)² = (4/27)². Double-crossing amplitude locked.',
        },
        {
            'name': 'm_95 (fold-wall scalar)',
            'predicted': 91.1876 * math.sqrt(1 + 2 * M.eta ** 2),
            'measured': 95.4,
            'unit': 'GeV',
            'kill_threshold': 2.0,
            'experiment': 'CMS/ATLAS di-photon, LEP reanalysis',
            'timeframe': '2025-2028',
            'why_fatal': 'Unique BSM prediction. If no signal at ~95 GeV: LOTUS is wrong.',
        },
        {
            'name': 'Jarlskog J',
            'predicted': u.ckm['J'],
            'measured': 3.08e-5,
            'unit': '',
            'kill_threshold': 5.0,
            'experiment': 'LHCb, Belle II',
            'timeframe': '2025-2030',
            'why_fatal': 'J = A²λ⁶η̄. All ingredients from spectral invariants.',
        },
        {
            'name': 'CKM γ angle',
            'predicted': u.ckm['gamma_deg'],
            'measured': 65.6,
            'unit': 'deg',
            'kill_threshold': 3.0,
            'experiment': 'LHCb, Belle II',
            'timeframe': '2025-2030',
            'why_fatal': 'γ = arctan(2π²/9). Rational function of π only.',
        },
        {
            'name': 'Baryogenesis η_B',
            'predicted': u.eta_B,
            'measured': 6.1e-10,
            'unit': '',
            'kill_threshold': 10.0,
            'experiment': 'CMB-S4, Simons Observatory',
            'timeframe': '2025-2030',
            'why_fatal': 'η_B = α⁴η. Currently the largest error; most fragile prediction.',
        },
    ]

    # Compute current error
    for t in targets:
        if t['measured'] != 0:
            t['error_pct'] = abs(t['predicted'] - t['measured']) / abs(t['measured']) * 100
        else:
            t['error_pct'] = 0.0

    return sorted(targets, key=lambda x: x['kill_threshold'])


def print_falsification_targets():
    """Print the falsification target table to stdout."""
    targets = falsification_targets()

    print()
    print("  LOTUS FALSIFICATION TARGETS")
    print("  " + "=" * 72)
    print(f"  {'#':<3} {'Prediction':<22} {'Error':>7}  {'Kill @':>7}  {'Experiment':<24}  {'When'}")
    print("  " + "-" * 72)

    for i, t in enumerate(targets, 1):
        print(f"  {i:<3} {t['name']:<22} {t['error_pct']:>6.2f}%  {t['kill_threshold']:>6.1f}%  "
              f"{t['experiment']:<24}  {t['timeframe']}")

    print("  " + "-" * 72)
    print()
    print("  KILL CONDITIONS:")
    for t in targets:
        print(f"    • {t['name']}: {t['why_fatal']}")
    print()
    print("  If ANY of these disagree beyond threshold: LOTUS is dead.")
    print("  No patches. No epicycles. The geometry is exact or it is wrong.")
    print()


def resolution_comparison():
    """Compare predictions at all three resolution levels.

    Shows how hurricane corrections improve accuracy.

    Returns:
        dict mapping resolution -> {name: error_pct}
    """
    from lotus import Universe

    result = {}
    for res in ('tree', '1_loop', '2_loop'):
        u = Universe(resolution=res)
        preds = u.predictions
        result[res] = {}
        for name, vals in preds.items():
            result[res][name] = vals['error_pct']

    return result


def print_resolution_comparison():
    """Print a side-by-side RMS comparison across resolution levels."""
    data = resolution_comparison()

    print()
    print("  RESOLUTION COMPARISON")
    print("  " + "=" * 50)

    for res in ('tree', '1_loop', '2_loop'):
        errs = [v for v in data[res].values() if v > 0]
        rms = math.sqrt(sum(e ** 2 for e in errs) / len(errs)) if errs else 0
        median = sorted(errs)[len(errs) // 2] if errs else 0
        mx = max(errs) if errs else 0
        print(f"  {res:<8}  RMS: {rms:>6.3f}%  "
              f"Median: {median:>6.3f}%  Max: {mx:>5.1f}%")

    print()
