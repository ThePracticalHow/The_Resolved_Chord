"""
lotus.cli — Console script entry point for 'lotus-boot'.
"""
import sys


def main():
    """Boot the universe. Called by the 'lotus-boot' console script."""
    import lotus

    args = sys.argv[1:]

    # --equations: print all derived classical equations
    if '--equations' in args:
        from lotus.equations import print_equations
        print_equations()
        return

    # --derive NAME: show provenance for a specific equation
    if '--derive' in args:
        idx = args.index('--derive')
        if idx + 1 < len(args):
            name = args[idx + 1]
        else:
            print("Usage: lotus-boot --derive <equation>")
            print("  e.g.: lotus-boot --derive F=ma")
            sys.exit(1)
        from lotus.equations import derive
        import json
        result = derive(name)
        # Pretty-print the chain prominently
        print(f"\n  {result.get('name', name)}")
        print(f"  {'=' * len(result.get('name', name))}")
        if 'law' in result:
            print(f"  Law:   {result['law']}")
        elif 'equation' in result:
            print(f"  Law:   {result['equation']}")
        print(f"  Chain: {result.get('chain', '?')}")
        print()
        origin = result.get('spectral_origin', {})
        if origin:
            print("  Spectral origin:")
            for k, v in origin.items():
                print(f"    {k}: {v}")
        print()
        return

    # --predict NAME: show full provenance for one prediction
    if '--predict' in args:
        idx = args.index('--predict')
        if idx + 1 < len(args):
            name = args[idx + 1]
        else:
            print("Usage: lotus-boot --predict <name>")
            print("  e.g.: lotus-boot --predict m_H")
            print("  e.g.: lotus-boot --predict alpha_s")
            sys.exit(1)
        u = lotus.Universe()
        try:
            result = u.predict(name)
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
        return

    # --book [SECTOR]: print the equation book
    if '--book' in args:
        from lotus.book import print_all, print_sector
        idx = args.index('--book')
        if idx + 1 < len(args) and not args[idx + 1].startswith('--'):
            print_sector(args[idx + 1])
        else:
            print_all()
        return

    # Default: boot the universe
    # Check for --resolution flag
    resolution = '2_loop'
    if '--resolution' in args:
        idx = args.index('--resolution')
        if idx + 1 < len(args):
            resolution = args[idx + 1]

    u = lotus.Universe(resolution=resolution)
    u.boot()



if __name__ == "__main__":
    main()
