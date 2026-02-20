#!/usr/bin/env python3
"""
Run the falsification test suite (pytest) and write results to a file.

Usage:
  python tools/falsification/run_falsification.py              # Run full suite
  python tools/falsification/run_falsification.py -m "not slow" # Skip slow tests
  python tools/falsification/run_falsification.py --results falsification_results.txt

Run from public-release/ root.
"""

import argparse
import subprocess
import sys
from datetime import datetime
from pathlib import Path


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-m", "--markers", help="Pytest markers (e.g. 'not slow')")
    ap.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    ap.add_argument("--results", metavar="FILE", default="falsification_results.txt",
                    help="Write pass/fail summary to file (default: falsification_results.txt)")
    args = ap.parse_args()

    project_root = Path(__file__).resolve().parent.parent.parent
    cmd = [sys.executable, "-m", "pytest", "tools/falsification/", "-v" if args.verbose else "-q"]
    if args.markers:
        cmd.extend(["-m", args.markers])

    print("=== Falsification Test Suite ===\n")
    r = subprocess.run(cmd, cwd=project_root)
    exit_code = r.returncode

    try:
        ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        status = "PASS" if exit_code == 0 else "FAIL"
        with open(project_root / args.results, "w", encoding="utf-8") as f:
            f.write(f"{ts} | Falsification: {status} (exit {exit_code})\n")
    except OSError as e:
        print(f"Warning: could not write results: {e}")

    sys.exit(exit_code)


if __name__ == "__main__":
    main()
