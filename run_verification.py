#!/usr/bin/env python3
"""
Run all core verification scripts.
Ensures documentation-listed scripts execute successfully.

Usage:
  python run_verification.py           # Run all, report pass/fail
  python run_verification.py --list   # List scripts only
  python run_verification.py -q       # Quiet: only show failures
  python run_verification.py --results verification_results.txt  # Write summary to file

Run from public-release/ root.
"""

import argparse
from datetime import datetime
import subprocess
import sys
from pathlib import Path

# Core verification scripts (matches README + GETTING_STARTED + new paper-accurate)
VERIFICATION_SCRIPTS = [
    "verification/EtaInvariant.py",
    "verification/ghost_parseval_proof.py",
    "verification/gravity_theorem_proof.py",
    "verification/alpha_lag_proof.py",
    "verification/alpha_s_theorem.py",
    "verification/ckm_complete.py",
    "verification/cc_aps_proof.py",
    "verification/vev_stiffness_proof.py",
    "verification/quantum_gravity_lotus.py",
    "verification/fold_potential_paper.py",
    "verification/vev_overlap_paper.py",
]

# Fast scripts only (for --quick: ~30s total)
QUICK_SCRIPTS = [
    "verification/EtaInvariant.py",
    "verification/vev_stiffness_proof.py",
    "verification/fold_potential_paper.py",
    "verification/vev_overlap_paper.py",
    "verification/alpha_lag_proof.py",
]


def run_script(script_path: Path, quiet: bool = False) -> bool:
    """Run a verification script; return True if exit code 0."""
    try:
        r = subprocess.run(
            [sys.executable, str(script_path)],
            cwd=script_path.parent.parent,
            capture_output=quiet,
            text=True,
            timeout=120,
        )
        return r.returncode == 0
    except subprocess.TimeoutExpired:
        if not quiet:
            print(f"  TIMEOUT: {script_path.name}")
        return False
    except Exception as e:
        if not quiet:
            print(f"  ERROR: {script_path.name}: {e}")
        return False


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--list", "-l", action="store_true", help="List scripts only")
    ap.add_argument("-q", "--quiet", action="store_true", help="Only show failures")
    ap.add_argument("--quick", action="store_true", help="Run only fast scripts (~30s)")
    ap.add_argument("--results", metavar="FILE", help="Write pass/fail summary to file (e.g. verification_results.txt)")
    args = ap.parse_args()

    root = Path(__file__).resolve().parent
    script_list = QUICK_SCRIPTS if args.quick else VERIFICATION_SCRIPTS
    scripts = [root / s for s in script_list]

    if args.list:
        for s in (QUICK_SCRIPTS if args.quick else VERIFICATION_SCRIPTS):
            print(s)
        return

    if not args.quiet:
        print("Running verification scripts...")
        print()

    ok = 0
    fail = 0
    for script in scripts:
        if not script.exists():
            if not args.quiet:
                print(f"  SKIP (missing): {script.name}")
            fail += 1
            continue
        success = run_script(script, quiet=args.quiet)
        if success:
            ok += 1
            if not args.quiet:
                print(f"  OK: {script.name}")
        else:
            fail += 1
            if not args.quiet:
                print(f"  FAIL: {script.name}")

    if not args.quiet:
        print()
        print(f"Result: {ok} ok, {fail} failed")

    exit_code = 1 if fail else 0

    if args.results:
        try:
            ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            with open(args.results, "w", encoding="utf-8") as f:
                f.write(f"{ts} | Verification: {ok} ok, {fail} failed | exit {exit_code}\n")
        except OSError as e:
            if not args.quiet:
                print(f"Warning: could not write results to {args.results}: {e}")

    sys.exit(exit_code)


if __name__ == "__main__":
    main()
