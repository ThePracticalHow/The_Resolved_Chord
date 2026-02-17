#!/usr/bin/env python3
"""
Shared pandoc utilities for convert_tex_pdf_md.py.
Single source of truth for finding the pandoc executable.
"""

import os
import shutil
import subprocess
from pathlib import Path


def _pandoc_search_paths() -> list[tuple[str, Path]]:
    """Paths to check for pandoc.exe (name, path) for display/debug."""
    return [
        ("LOCALAPPDATA", Path(os.environ.get("LOCALAPPDATA", "")) / "Pandoc" / "pandoc.exe"),
        ("ProgramFiles", Path(os.environ.get("ProgramFiles", "C:\\Program Files")) / "Pandoc" / "pandoc.exe"),
        ("ProgramFiles(x86)", Path(os.environ.get("ProgramFiles(x86)", "C:\\Program Files (x86)")) / "Pandoc" / "pandoc.exe"),
        ("C:\\pandoc", Path("C:/pandoc/pandoc.exe")),
    ]


def find_pandoc() -> str | None:
    """
    Find pandoc: PATH first, then common Windows locations, then cmd where.
    Returns path to pandoc.exe or None.
    """
    exe = shutil.which("pandoc")
    if exe:
        return exe

    for _name, base in _pandoc_search_paths():
        if base and base.exists():
            return str(base)

    try:
        r = subprocess.run(
            ["cmd", "/c", "where", "pandoc"],
            capture_output=True,
            text=True,
            timeout=5,
        )
        if r.returncode == 0 and r.stdout.strip():
            return r.stdout.strip().splitlines()[0].strip()
    except Exception:
        pass

    return None


def get_pandoc_status() -> dict:
    """Return dict with pandoc path, PATH snippet, and search-path results for debug output."""
    path_val = os.environ.get("PATH", "")
    return {
        "path_snippet": path_val[:500],
        "which": shutil.which("pandoc"),
        "found": find_pandoc(),
        "search_paths": _pandoc_search_paths(),
    }
