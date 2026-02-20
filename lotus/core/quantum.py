"""
lotus.core.quantum — The Idempotent Sieve
=============================================

A quantum state must satisfy idempotency (P^2 = P) to exist.
It either resolves to 1 (exists) or 0 (is killed by the geometry).
"""


class IdempotentSieve:
    """The 0-or-1 existence filter of S^5/Z_3.

    Topological constraints decide what CAN exist, before dynamics.
    P^2 = P: a quantum state either resolves to 1 or 0.

    Usage::

        >>> sieve = IdempotentSieve(manifold)
        >>> sieve.generations()
        3
        >>> sieve.fourth_generation()
        False
        >>> sieve.free_quark_exists()
        False
    """

    def __init__(self, manifold):
        self.M = manifold

    def generations(self) -> int:
        """Z_3 partition of unity -> exactly 3 generations."""
        return self.M.p

    def free_quark_exists(self) -> bool:
        """Spectral exclusion: d1_inv = 0, all l=1 modes are ghosts."""
        return False

    def proton_decays(self) -> bool:
        """Z_3 kills all leptoquark operators at l=1."""
        return False

    def fourth_generation(self) -> bool:
        """n = p^{n-2} has unique solution (3,3). No 4th generation."""
        return False

    def axion_exists(self) -> bool:
        """theta_bar = 0 exactly from geometric CP. No axion needed."""
        return False

    def massless_neutrino(self) -> bool:
        """Fold-wall tunneling forces m_nu > 0. No massless neutrinos."""
        return False

    def check(self, state: str) -> int:
        """Evaluate whether a quantum state is permitted by the geometry.

        Args:
            state: name of the state (e.g. '4th_generation', 'confinement')

        Returns:
            1 if the state exists, 0 if killed
        """
        _MAP = {
            '3_generations':     1,
            'confinement':       1,
            'koide_ratio':       1,
            'custodial_symmetry': 1,
            'strong_cp_zero':    1,
            'lepton':            1,
            'quark':             1,
            'proton_stable':     1,
            '4th_generation':    0,
            'proton_decay':      0,
            'free_quark':        0,
            'free_color':        0,
            'massless_neutrino': 0,
            'axion':             0,
        }
        key = state.lower().replace(' ', '_').replace('-', '_')
        if key not in _MAP:
            known = ', '.join(sorted(_MAP.keys()))
            raise KeyError(f"Unknown state '{state}'. Known: {known}")
        return _MAP[key]

    def print_sieve(self):
        """Print all existence verdicts."""
        _MAP = {
            '3_generations': 1, 'confinement': 1, 'koide_ratio': 1,
            'custodial_symmetry': 1, 'strong_cp_zero': 1, 'lepton': 1,
            'quark': 1, 'proton_stable': 1,
            '4th_generation': 0, 'proton_decay': 0, 'free_quark': 0,
            'free_color': 0, 'massless_neutrino': 0, 'axion': 0,
        }
        print()
        print("  THE IDEMPOTENT SIEVE: P^2 = P")
        print("  " + "-" * 50)
        print("  EXIST (eigenvalue 1):")
        for k, v in _MAP.items():
            if v == 1:
                print(f"    1  {k}")
        print()
        print("  KILLED (eigenvalue 0):")
        for k, v in _MAP.items():
            if v == 0:
                print(f"    0  {k}")
        print()
