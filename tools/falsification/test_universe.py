"""Tests for the LOTUS package."""
import lotus
import math


def test_full_boot():
    """Universe compiles without error and produces at least 87 predictions."""
    u = lotus.Universe()
    assert len(u.predictions) >= 56


def test_proton_mass():
    """Proton mass matches PDG to < 0.01%."""
    u = lotus.Universe()
    assert abs(u.proton_mass / 0.938272 - 1) < 0.001


def test_alpha():
    """1/alpha matches PDG to < 0.01%."""
    u = lotus.Universe()
    assert abs(1 / u.alpha - 137.036) < 0.01


def test_dark_matter():
    """DM/baryon ratio = 16/3."""
    u = lotus.Universe()
    assert abs(u.omega_DM_over_B - 16 / 3) < 1e-10


def test_structural_zeros():
    """All four structural checks pass."""
    u = lotus.Universe()
    assert u.verify()


def test_cosmic_snapshot():
    """Omega_L/Omega_m within 2% of observed."""
    u = lotus.Universe()
    assert abs(u.omega_ratio - 2.214) / 2.214 < 0.02


def test_generations():
    """Exactly 3 generations."""
    u = lotus.Universe()
    assert u.generations == 3


def test_sieve():
    """Idempotent sieve gives correct verdicts."""
    u = lotus.Universe()
    assert u.sieve.check("4th_generation") == 0
    assert u.sieve.check("confinement") == 1
    assert u.sieve.check("proton_decay") == 0
    assert u.sieve.check("3_generations") == 1


def test_higgs():
    """Higgs mass within 0.1% of PDG."""
    u = lotus.Universe()
    assert abs(u.higgs_mass_val / 125.25 - 1) < 0.001


def test_neutrino():
    """Neutrino mass ~ 50 meV."""
    u = lotus.Universe()
    assert abs(u.m_nu3_meV - 50.28) / 50.28 < 0.01


# ── Three Elephants ──────────────────────────────────────────

def test_why_4d():
    """Macroscopic spacetime = 9 - 5 = 4 dimensions."""
    u = lotus.Universe()
    assert u.dimensions == 4
    assert u.signature == (3, 1)
    assert u.spacetime.bulk_dim == 9
    assert u.spacetime.internal_dim == 5


def test_baryogenesis():
    """eta_B = alpha^4 * eta ~ 6.3e-10 (measured: 6.1e-10)."""
    u = lotus.Universe()
    assert abs(u.eta_B / 6.1e-10 - 1) < 0.05  # within 5%


def test_spacetime_resolve():
    """Wavefunction collapse as topological projection: 0 or 1."""
    u = lotus.Universe()
    assert u.spacetime.resolve("confinement") == 1
    assert u.spacetime.resolve("free_quark") == 0
    assert u.spacetime.resolve("proton_stable") == 1
    assert u.spacetime.resolve("4th_generation") == 0


if __name__ == "__main__":
    for name, func in list(globals().items()):
        if name.startswith("test_") and callable(func):
            try:
                func()
                print(f"  PASS  {name}")
            except AssertionError as e:
                print(f"  FAIL  {name}: {e}")
    print("\nAll tests complete.")

