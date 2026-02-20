#!/usr/bin/env python3
"""Auto-clean utility for public-release.

Default mode is dry-run. Use --apply to delete.
"""

from __future__ import annotations

import argparse
import fnmatch
import shutil
from pathlib import Path

DIR_NAMES = {
    "__pycache__",
    ".pytest_cache",
    ".mypy_cache",
    "build",
    "dist",
    "results",
    "output",
}
FILE_PATTERNS = ["*.pyc", "*.pyo"]


def collect_targets(root: Path, include_results: bool) -> tuple[list[Path], list[Path]]:
    dirs: list[Path] = []
    files: list[Path] = []

    for p in root.rglob("*"):
        if ".git" in p.parts or "_archive" in p.parts:
            continue
        if p.is_dir() and p.name in DIR_NAMES:
            dirs.append(p)
        elif p.is_file() and any(fnmatch.fnmatch(p.name, pat) for pat in FILE_PATTERNS):
            files.append(p)

    for egg in root.glob("*.egg-info"):
        if egg.is_dir():
            dirs.append(egg)

    coverage_file = root / ".coverage"
    if coverage_file.exists():
        files.append(coverage_file)

    if include_results:
        results_file = root / "results.txt"
        if results_file.exists():
            files.append(results_file)

    return sorted(dirs), sorted(files)


def main() -> int:
    parser = argparse.ArgumentParser(description="Clean generated artifacts from public-release.")
    parser.add_argument("--apply", action="store_true", help="Actually delete matched files/directories.")
    parser.add_argument("--include-results", action="store_true", help="Also remove results.txt in project root.")
    args = parser.parse_args()

    root = Path(__file__).resolve().parent.parent
    dirs, files = collect_targets(root, include_results=args.include_results)

    mode = "APPLY" if args.apply else "DRY-RUN"
    print("=" * 64)
    print(f"LOTUS AUTO CLEAN ({mode})")
    print("=" * 64)
    print(f"Root: {root}")
    print(f"Directories matched: {len(dirs)}")
    print(f"Files matched: {len(files)}")

    for d in dirs[:20]:
        print(f"  [dir]  {d.relative_to(root)}")
    if len(dirs) > 20:
        print(f"  ... and {len(dirs) - 20} more directories")

    for f in files[:20]:
        print(f"  [file] {f.relative_to(root)}")
    if len(files) > 20:
        print(f"  ... and {len(files) - 20} more files")

    if not args.apply:
        print("\nNo changes made. Re-run with --apply to delete.")
        return 0

    for d in dirs:
        shutil.rmtree(d, ignore_errors=True)
    for f in files:
        try:
            f.unlink(missing_ok=True)
        except TypeError:
            if f.exists():
                f.unlink()

    print("\n✅ Cleanup complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
