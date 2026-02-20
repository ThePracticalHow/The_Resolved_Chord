#!/usr/bin/env python3
"""
Run all core verification scripts.
Ensures documentation-listed scripts execute successfully.

Usage:
  python verification/run_verification.py           # Run all, report pass/fail
  python verification/run_verification.py --list   # List scripts only
  python verification/run_verification.py -q       # Quiet: only show failures
  python verification/run_verification.py --results verification_results.txt

Run from public-release/ root.
"""

import argparse
import hashlib
import json
from datetime import datetime
import os
import re
import subprocess
import sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import py_compile

# Force UTF-8 in subprocesses (Windows cp1252 chokes on Greek/math chars)
_UTF8_ENV = {**os.environ, "PYTHONIOENCODING": "utf-8"}

# Core verification scripts (matches README + GETTING_STARTED + new paper-accurate)
VERIFICATION_SCRIPTS = [
    "EtaInvariant.py",
    "ghost_parseval_proof.py",
    "gravity_theorem_proof.py",
    "alpha_lag_proof.py",
    "alpha_s_theorem.py",
    "ckm_complete.py",
    "cc_aps_proof.py",
    "vev_stiffness_proof.py",
    "quantum_gravity_lotus.py",
    "fold_potential_paper.py",
    "vev_overlap_paper.py",
    "hurricane_proof.py",
    "lotus_aps_generation.py",
]

# Fast scripts only (for --quick: ~30s total)
QUICK_SCRIPTS = [
    "EtaInvariant.py",
    "vev_stiffness_proof.py",
    "fold_potential_paper.py",
    "vev_overlap_paper.py",
    "alpha_lag_proof.py",
]

CACHE_FILE = ".verification_cache.json"
HYGIENE_TARGET_DIRS = [
    "public-release/verification",
    "scripts",
    "tests",
]
HYGIENE_EXCLUDE_PARTS = {
    "__pycache__",
    ".pytest_cache",
    "_archive",
}
HYGIENE_EXCLUDE_FILES = {
    "__init__.py",
}


def _file_hash(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def _load_cache(cache_path: Path) -> dict:
    if not cache_path.exists():
        return {}
    try:
        return json.loads(cache_path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return {}


def _save_cache(cache_path: Path, cache: dict) -> None:
    cache_path.write_text(json.dumps(cache, indent=2, sort_keys=True), encoding="utf-8")


def discover_scripts(verification_dir: Path) -> list[Path]:
    ignored = {"run_verification.py", "__init__.py"}
    scripts = [p for p in verification_dir.glob("*.py") if p.name not in ignored]
    return sorted(scripts)


def filter_recent_scripts(scripts: list[Path], minutes: int) -> list[Path]:
    if minutes <= 0:
        return scripts
    now = datetime.now().timestamp()
    threshold_seconds = minutes * 60
    return [p for p in scripts if (now - p.stat().st_mtime) <= threshold_seconds]


def collect_hygiene_targets(project_root: Path, selected_scripts: list[Path], hygiene_all: bool) -> list[Path]:
    if not hygiene_all:
        return [p for p in selected_scripts if p.exists()]

    targets: list[Path] = []
    for rel_dir in HYGIENE_TARGET_DIRS:
        base = project_root / rel_dir
        if not base.exists():
            continue
        for path in base.rglob("*.py"):
            if any(part in HYGIENE_EXCLUDE_PARTS for part in path.parts):
                continue
            if path.name in HYGIENE_EXCLUDE_FILES:
                continue
            targets.append(path)
    return sorted(set(targets))


def hygiene_check(scripts: list[Path]) -> list[str]:
    issues = []
    bare_except_re = re.compile(r"^\s*except\s*:\s*$")
    wildcard_import_re = re.compile(r"^\s*from\s+\S+\s+import\s+\*\s*$")

    for script in scripts:
        try:
            py_compile.compile(str(script), doraise=True)
        except py_compile.PyCompileError as e:
            issues.append(f"{script.name}: syntax error ({e.msg})")
            continue

        try:
            text = script.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            issues.append(f"{script.name}: non-UTF-8 file encoding")
            continue

        lines = text.splitlines()
        for idx, line in enumerate(lines, 1):
            if line.rstrip() != line:
                issues.append(f"{script.name}:{idx}: trailing whitespace")
            if bare_except_re.match(line):
                issues.append(f"{script.name}:{idx}: bare except")
            if wildcard_import_re.match(line):
                issues.append(f"{script.name}:{idx}: wildcard import")

        if text and not text.endswith("\n"):
            issues.append(f"{script.name}: missing newline at EOF")

    return issues


def run_script(script_path: Path, quiet: bool = False, timeout: int = 120) -> tuple[bool, str]:
    """Run a verification script; return (ok, reason)."""
    try:
        # Ensure the script exists and is executable
        if not script_path.exists():
            if not quiet:
                print(f"  MISSING: {script_path.name}")
            return False, "missing"

        # Robust environment setup
        env = dict(os.environ)
        env.update(_UTF8_ENV)

        # Try to run the script
        r = subprocess.run(
            [sys.executable, str(script_path)],
            cwd=script_path.parent.parent,  # public-release root
            capture_output=quiet,
            text=True,
            encoding='utf-8',
            timeout=timeout,
            env=env,
        )
        return r.returncode == 0, f"exit={r.returncode}"
    except subprocess.TimeoutExpired:
        if not quiet:
            print(f"  TIMEOUT: {script_path.name}")
        return False, "timeout"
    except FileNotFoundError:
        if not quiet:
            print(f"  NO PYTHON: {script_path.name}")
        return False, "no-python"
    except PermissionError:
        if not quiet:
            print(f"  PERMISSION DENIED: {script_path.name}")
        return False, "permission"
    except OSError as e:
        if not quiet:
            print(f"  OS ERROR: {script_path.name} - {e}")
        return False, f"os-error: {e}"
    except Exception as e:
        if not quiet:
            print(f"  UNEXPECTED ERROR: {script_path.name} - {e}")
        return False, f"unexpected: {e}"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--list", "-l", action="store_true", help="List scripts only")
    ap.add_argument("-q", "--quiet", action="store_true", help="Only show failures")
    ap.add_argument("--quick", action="store_true", help="Run only fast scripts (~30s)")
    ap.add_argument("--results", metavar="FILE", help="Write pass/fail summary to file (e.g. verification_results.txt)")
    ap.add_argument("--discover-all", action="store_true", help="Auto-discover all *.py scripts under verification/")
    ap.add_argument("--new-minutes", type=int, default=0, help="Only run scripts modified in the last N minutes")
    ap.add_argument("--changed-only", action="store_true", help="Skip unchanged scripts that previously passed")
    ap.add_argument("--workers", type=int, default=1, help="Parallel workers for script execution")
    ap.add_argument("--timeout", type=int, default=120, help="Per-script timeout in seconds")
    ap.add_argument("--hygiene", action="store_true", help="Run coding hygiene checks on selected scripts")
    ap.add_argument("--strict-hygiene", action="store_true", help="Fail if hygiene issues are found")
    ap.add_argument(
        "--hygiene-all",
        action="store_true",
        help="Run hygiene checks across public-release/verification, scripts, and tests",
    )
    args = ap.parse_args()

    project_root = Path(__file__).resolve().parent.parent
    verification_dir = Path(__file__).resolve().parent
    if args.discover_all:
        scripts = discover_scripts(verification_dir)
    else:
        script_list = QUICK_SCRIPTS if args.quick else VERIFICATION_SCRIPTS
        scripts = [verification_dir / s for s in script_list]

    if args.new_minutes > 0:
        scripts = filter_recent_scripts(scripts, args.new_minutes)

    cache_path = project_root / CACHE_FILE
    cache = _load_cache(cache_path)

    if args.changed_only:
        filtered = []
        for script in scripts:
            if not script.exists():
                filtered.append(script)
                continue
            digest = _file_hash(script)
            cached = cache.get(str(script.resolve()), {})
            if cached.get("hash") == digest and cached.get("ok") is True:
                continue
            filtered.append(script)
        scripts = filtered

    if args.list:
        for s in scripts:
            print(s.name)
        return

    if not args.quiet:
        print("Running verification scripts...")
        print()

    if args.hygiene:
        hygiene_targets = collect_hygiene_targets(project_root, scripts, args.hygiene_all)
        hygiene_issues = hygiene_check(hygiene_targets)
        if hygiene_issues and not args.quiet:
            print("Hygiene issues found:")
            for issue in hygiene_issues:
                print(f"  - {issue}")
            print()
        if args.strict_hygiene and hygiene_issues:
            if not args.quiet:
                print("Result: hygiene failed")
            sys.exit(2)

    ok = 0
    fail = 0
    results = []
    if args.workers > 1:
        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            future_map = {
                executor.submit(run_script, script, args.quiet, args.timeout): script
                for script in scripts
            }
            for future in as_completed(future_map):
                script = future_map[future]
                success, reason = future.result()
                results.append((script, success, reason))
    else:
        for script in scripts:
            success, reason = run_script(script, quiet=args.quiet, timeout=args.timeout)
            results.append((script, success, reason))

    for script, success, reason in sorted(results, key=lambda x: x[0].name.lower()):
        if success:
            ok += 1
            if not args.quiet:
                print(f"  OK: {script.name}")
        else:
            fail += 1
            if not args.quiet:
                print(f"  FAIL: {script.name} ({reason})")

        if script.exists():
            cache[str(script.resolve())] = {
                "hash": _file_hash(script),
                "ok": success,
                "ts": datetime.now().isoformat(timespec="seconds"),
            }

    if not args.quiet:
        print()
        print(f"Result: {ok} ok, {fail} failed")

    exit_code = 1 if fail else 0

    try:
        _save_cache(cache_path, cache)
    except OSError:
        pass

    if args.results:
        try:
            ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            with open(args.results, "w", encoding="utf-8") as f:
                hygiene_suffix = ""
                if args.hygiene:
                    hygiene_suffix = f" | hygiene_issues={len(hygiene_issues) if 'hygiene_issues' in locals() else 0}"
                f.write(f"{ts} | Verification: {ok} ok, {fail} failed{hygiene_suffix} | exit {exit_code}\n")
        except OSError as e:
            if not args.quiet:
                print(f"Warning: could not write results to {args.results}: {e}")

    sys.exit(exit_code)


if __name__ == "__main__":
    main()
