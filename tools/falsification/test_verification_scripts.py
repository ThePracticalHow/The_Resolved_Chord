"""
Pytest: Verify that all core verification scripts run successfully (exit 0).

These scripts are listed in README.md and GETTING_STARTED.md.
Run: pytest tools/falsification/test_verification_scripts.py -v
"""

import os
import subprocess
import sys
from pathlib import Path

import pytest

# Must match run_verification.py VERIFICATION_SCRIPTS
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
    "verification/hurricane_proof.py",
    "verification/lotus_aps_generation.py",
]

# Scripts that take >10s; marked slow for pytest -m "not slow"
SLOW_SCRIPTS = {
    "verification/ghost_parseval_proof.py",
    "verification/gravity_theorem_proof.py",
    "verification/quantum_gravity_lotus.py",
}

_script_params = [
    pytest.param(s, marks=pytest.mark.slow) if s in SLOW_SCRIPTS else s
    for s in VERIFICATION_SCRIPTS
]

# Force UTF-8 encoding for subprocess stdout/stderr on Windows
# (CP1252 can't handle Unicode chars like α, φ, ⁵, box-drawing
#  used in verification script output).
_UTF8_ENV = {**os.environ, "PYTHONIOENCODING": "utf-8"}


@pytest.fixture(scope="module")
def project_root():
    return Path(__file__).resolve().parent.parent


@pytest.mark.parametrize("script", _script_params)
def test_verification_script_runs(project_root, script):
    """Each verification script must run and exit 0."""
    path = project_root / script
    if not path.exists():
        pytest.skip(f"Script not found: {script}")
    r = subprocess.run(
        [sys.executable, str(path)],
        cwd=project_root,
        capture_output=True,
        text=True,
        timeout=120,
        env=_UTF8_ENV,
    )
    assert r.returncode == 0, (
        f"{script} failed (exit {r.returncode})\n"
        f"stderr: {r.stderr[:500] if r.stderr else '(none)'}"
    )
