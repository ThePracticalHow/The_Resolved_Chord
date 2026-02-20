"""
python -m lotus -- Run the LOTUS engine from the command line.

Usage:
    python -m lotus                          Boot the universe (default: 2_loop)
    python -m lotus --resolution tree        Boot at tree level
    python -m lotus --predict m_H            Show full provenance for one prediction
    python -m lotus --sieve                  Print the idempotent sieve
    python -m lotus --json                   Print all predictions as JSON
    python -m lotus --compare                Print prediction vs PDG table
    python -m lotus --falsify                Print falsification targets
    python -m lotus --landscape              Run the universe selection pipeline
    python -m lotus --dynamics               Print universe snapshot at phi_lotus
    python -m lotus --arrow                  Print arrow-of-time analysis
    python -m lotus --resolution-compare     Compare all resolution levels
    python -m lotus --book                   Print the equation book (all sectors)
    python -m lotus --book QCD               Print equations for one sector
    python -m lotus --constants              Print all derived constants
    python -m lotus --equations              Print all derived classical equations
"""
import sys
import argparse


def main():
    parser = argparse.ArgumentParser(
        prog='python -m lotus',
        description='LOTUS — Lens Orbifold Theory of the Unified Spectrum',
    )
    parser.add_argument(
        '--resolution', choices=['tree', '1_loop', '2_loop'],
        default='2_loop',
        help='Resolution level: tree, 1_loop, or 2_loop (default: 2_loop)',
    )
    parser.add_argument('--json', action='store_true',
                        help='Print all predictions as JSON')
    parser.add_argument('--sieve', action='store_true',
                        help='Print the idempotent sieve')
    parser.add_argument('--compare', action='store_true',
                        help='Print prediction vs PDG comparison table')
    parser.add_argument('--falsify', action='store_true',
                        help='Print falsification targets with kill thresholds')
    parser.add_argument('--landscape', action='store_true',
                        help='Run the universe selection pipeline')
    parser.add_argument('--dynamics', action='store_true',
                        help='Print universe snapshot at phi_lotus')
    parser.add_argument('--arrow', action='store_true',
                        help='Print arrow-of-time analysis')
    parser.add_argument('--resolution-compare', action='store_true',
                        dest='resolution_compare',
                        help='Compare predictions across all resolution levels')
    parser.add_argument('--book', nargs='?', const='ALL', default=None,
                        metavar='SECTOR',
                        help='Print the equation book (optionally for a single sector)')
    parser.add_argument('--constants', action='store_true',
                        help='Print all derived constants with provenance')
    parser.add_argument('--predict', metavar='NAME', default=None,
                        help='Show full provenance for one prediction (e.g. m_H, 1/alpha)')
    parser.add_argument('--equations', action='store_true',
                        help='Print all classical equations derived from Tr(f(D²/Λ²))')

    args = parser.parse_args()

    if args.predict is not None:
        import lotus
        u = lotus.Universe(resolution=args.resolution)
        try:
            result = u.predict(args.predict)
        except KeyError as e:
            print(f"\nError: {e.args[0]}")
            sys.exit(1)
        unit = f" {result['unit']}" if result['unit'] else ''
        pred = result['predicted']
        meas = result['measured']
        err  = result['error_pct']
        print()
        print(f"  {result['name']}")
        print(f"  {'─' * max(len(result['name']), 40)}")
        print(f"  Predicted:   {pred:.6g}{unit}")
        print(f"  Measured:    {meas:.6g}{unit}  ({err:.4f}%)")
        print(f"  Formula:     {result['derivation']}")
        print(f"  Invariants:  {', '.join(result['invariants'])}")
        print(f"  Experiment:  {result['experiment']}")
        print(f"  Resolution:  {result['resolution']}")
        print()

    elif args.json:
        import json
        import lotus
        u = lotus.Universe(resolution=args.resolution)
        print(json.dumps(u.predictions, indent=2, ensure_ascii=False))

    elif args.sieve:
        import lotus
        u = lotus.Universe(resolution=args.resolution)
        u.sieve.print_sieve()

    elif args.compare:
        from lotus.falsify import print_comparison
        print_comparison(resolution=args.resolution)

    elif args.falsify:
        from lotus.falsify import print_falsification_targets
        print_falsification_targets()

    elif args.landscape:
        from lotus.landscape import selection_pipeline
        selection_pipeline()

    elif args.dynamics:
        from lotus.dynamics import universe_at_phi, phi_lotus
        snap = universe_at_phi(phi_lotus)
        print()
        print(f"  UNIVERSE SNAPSHOT at φ = {phi_lotus:.4f}")
        print("  " + "=" * 50)
        for k, v in snap.items():
            if isinstance(v, float):
                if abs(v) > 100:
                    print(f"    {k:<20} {v:.2f}")
                elif abs(v) > 0.001:
                    print(f"    {k:<20} {v:.6f}")
                else:
                    print(f"    {k:<20} {v:.3e}")
            else:
                print(f"    {k:<20} {v}")
        print()

    elif args.arrow:
        from lotus.dynamics import arrow_of_time, phi_lotus, eta_of_phi
        print()
        print("  THE ARROW OF TIME")
        print("  " + "=" * 50)
        for phi_val in [0.0, 0.25, 0.5, 0.6, 0.75, phi_lotus, 1.0]:
            a = arrow_of_time(phi_val)
            label = " ← OUR UNIVERSE" if abs(phi_val - phi_lotus) < 0.01 else ""
            print(f"    φ={phi_val:.3f}  η={a['eta']:.6f}  "
                  f"θ={a['theta_deg']:>7.2f}°  "
                  f"T-broken={a['T_broken']}{label}")
        print()
        a_us = arrow_of_time(phi_lotus)
        print(f"  {a_us['description']}")
        print()

    elif args.resolution_compare:
        from lotus.falsify import print_resolution_comparison
        print_resolution_comparison()

    elif args.book is not None:
        from lotus import book
        if args.book == 'ALL':
            book.print_all()
        else:
            book.print_sector(args.book)

    elif args.constants:
        from lotus.constants import print_constants
        print_constants()

    elif args.equations:
        from lotus.equations import print_equations
        print_equations()

    else:
        import lotus
        u = lotus.Universe(resolution=args.resolution)
        u.boot()


if __name__ == "__main__":
    main()

