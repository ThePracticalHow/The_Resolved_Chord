#!/usr/bin/env python3
"""Release gate checks for public-release robustness.

Run from public-release root:
  python tools/release_gate.py
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

REQUIRED_FILES = [
    "README.md",
    "CONTRIBUTING.md",
    "CHANGELOG.md",
    "verification/compile_universe.py",
    "verification/run_verification.py",
    "tools/falsification/run_falsification.py",
    "papers/The_Resolved_Chord_v11.tex",
]

FORBIDDEN_PATTERNS = [
    (re.compile(r"file:///", re.IGNORECASE), "absolute file URI"),
    (re.compile(r"public-release/falsification/", re.IGNORECASE), "legacy falsification path"),
    (re.compile(r"papers/tex/", re.IGNORECASE), "legacy papers/tex path"),
]

MARKDOWN_LINK_RE = re.compile(r"\[[^\]]+\]\(([^)]+)\)")

DOC_SURFACE_FILES = [
    "README.md",
    "CONTRIBUTING.md",
    "CHANGELOG.md",
    "tools/README.md",
]

POLICY_FILE_PATTERNS = [
    ("README.md", re.compile(r"Python\s+3\.8\+", re.IGNORECASE), "stale Python minimum (3.8+)", True),
    ("RUN_EVERYTHING.bat", re.compile(r"Python\s+3\.8\+", re.IGNORECASE), "stale Python minimum (3.8+)", True),
    ("README.md", re.compile(r"https://github\.com/(?:ThePracticalHow|jixiang-leng)/lotus-bloom", re.IGNORECASE), "stale repository URL", True),
    ("CONTRIBUTING.md", re.compile(r"https://github\.com/(?:ThePracticalHow|jixiang-leng)/lotus-bloom", re.IGNORECASE), "stale repository URL", True),
    ("pyproject.toml", re.compile(r"https://github\.com/jixiang-leng/lotus-bloom", re.IGNORECASE), "stale repository URL in package metadata", True),
    ("CITATION.cff", re.compile(r"https://github\.com/(?:ThePracticalHow|jixiang-leng)/lotus-bloom", re.IGNORECASE), "stale repository URL in citation metadata", True),
]

POLICY_REQUIRED_PATTERNS = [
    ("pyproject.toml", re.compile(r"https://github\.com/ThePracticalHow/The_Resolved_Chord", re.IGNORECASE), "canonical repository URL missing in pyproject"),
    ("CITATION.cff", re.compile(r"https://github\.com/ThePracticalHow/The_Resolved_Chord", re.IGNORECASE), "canonical repository URL missing in citation"),
]


def collect_markdown_files(root: Path) -> list[Path]:
    files: list[Path] = []
    for rel in DOC_SURFACE_FILES:
        p = root / rel
        if p.exists():
            files.append(p)
    return files


def is_external_link(url: str) -> bool:
    return url.startswith(("http://", "https://", "mailto:", "#"))


def normalize_local_link(url: str) -> str:
    return url.split("#", 1)[0].strip()


def check_required_files(root: Path) -> list[str]:
    issues: list[str] = []
    for rel in REQUIRED_FILES:
        if not (root / rel).exists():
            issues.append(f"Missing required file: {rel}")
    return issues


def check_forbidden_patterns(root: Path, markdown_files: list[Path]) -> list[str]:
    issues: list[str] = []
    for md in markdown_files:
        text = md.read_text(encoding="utf-8", errors="replace")
        rel = md.relative_to(root)
        for pattern, label in FORBIDDEN_PATTERNS:
            if pattern.search(text):
                issues.append(f"{rel}: found {label}")
    return issues


def check_broken_local_links(root: Path, markdown_files: list[Path]) -> list[str]:
    issues: list[str] = []
    for md in markdown_files:
        text = md.read_text(encoding="utf-8", errors="replace")
        rel = md.relative_to(root)
        for raw in MARKDOWN_LINK_RE.findall(text):
            url = raw.strip()
            if is_external_link(url):
                continue
            link = normalize_local_link(url)
            if not link:
                continue
            target = (md.parent / link).resolve() if not link.startswith("/") else (root / link.lstrip("/")).resolve()
            if not target.exists():
                issues.append(f"{rel}: broken local link -> {url}")
    return issues


def check_policy_patterns(root: Path) -> list[str]:
    issues: list[str] = []

    for rel_path, pattern, label, is_forbidden in POLICY_FILE_PATTERNS:
        file_path = root / rel_path
        if not file_path.exists():
            continue
        text = file_path.read_text(encoding="utf-8", errors="replace")
        matched = bool(pattern.search(text))
        if is_forbidden and matched:
            issues.append(f"{rel_path}: found {label}")

    for rel_path, pattern, label in POLICY_REQUIRED_PATTERNS:
        file_path = root / rel_path
        if not file_path.exists():
            issues.append(f"{rel_path}: missing file for policy check")
            continue
        text = file_path.read_text(encoding="utf-8", errors="replace")
        if not pattern.search(text):
            issues.append(f"{rel_path}: {label}")

    return issues


def main() -> int:
    parser = argparse.ArgumentParser(description="Run robustness release gate checks.")
    parser.add_argument("--strict", action="store_true", help="Fail on warnings (currently informational only).")
    args = parser.parse_args()

    root = Path(__file__).resolve().parent.parent
    markdown_files = collect_markdown_files(root)

    print("=" * 64)
    print("LOTUS RELEASE GATE")
    print("=" * 64)
    print(f"Root: {root}")
    print(f"Markdown files scanned: {len(markdown_files)}")

    errors: list[str] = []
    errors.extend(check_required_files(root))
    errors.extend(check_forbidden_patterns(root, markdown_files))
    errors.extend(check_broken_local_links(root, markdown_files))
    errors.extend(check_policy_patterns(root))

    if errors:
        print("\n❌ FAIL: release gate found issues")
        for issue in errors:
            print(f"  - {issue}")
        return 1

    if args.strict:
        print("\n✅ PASS (strict mode): no issues found")
    else:
        print("\n✅ PASS: no issues found")
    return 0


if __name__ == "__main__":
    sys.exit(main())
