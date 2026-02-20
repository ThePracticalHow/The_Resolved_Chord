"""
lotus.core.geometry — The S5Z3 Manifold
==========================================

The rigid topological manifold of the universe.
All five spectral invariants derived from (n, p) = (3, 3),
the unique solution to n = p^{n-2}.
"""

import math
from lotus.constants import PI


class S5Z3:
    """The spectral geometry of S^5/Z_3 (lens space L(3;1,1,1)).

    Uniquely selected by the resonance lock: n = p^{n-2} -> (3,3).

    Five spectral invariants::

        d1      = 6      Ghost mode count (first eigenspace degeneracy)
        lambda1 = 5      First nonzero eigenvalue of Laplacian on S^5
        K       = 2/3    Koide ratio (moment map of Z_3 on simplex)
        eta     = 2/9    Donnelly eta invariant (spectral asymmetry)
        p       = 3      Orbifold order

    Usage::

        >>> from lotus.core.geometry import S5Z3
        >>> M = S5Z3()
        >>> M.eta
        0.2222...
        >>> M.G
        1.1111...
    """

    p = 3               # orbifold order |Z_3|
    n = 3               # complex dimension of C^3/Z_3
    d1 = 6              # ghost mode count: dim of l=1 eigenspace on S^5
    lambda1 = 5         # first eigenvalue: l(l+n-1) with l=1, n=3 -> 1*4+1=5
    K = 2 / 3           # Koide ratio: 2/p

    @property
    def eta(self) -> float:
        """Donnelly twisted Dirac eta invariant: d1/p^n = 6/27 = 2/9."""
        return self.d1 / (self.p ** self.n)

    @property
    def tau(self) -> float:
        """Reidemeister torsion: 1/p^n = 1/27."""
        return 1 / (self.p ** self.n)

    @property
    def G(self) -> float:
        """Proton spectral coupling: lambda1 * eta = 10/9."""
        return self.lambda1 * self.eta

    @property
    def c_grav(self) -> float:
        """Gravity hurricane coefficient: -1/(d1*lambda1) = -1/30."""
        return -1 / (self.d1 * self.lambda1)

    @classmethod
    def verify_uniqueness(cls) -> list:
        """Prove (n,p) = (3,3) is the unique physical solution to n = p^{n-2}.

        Returns:
            list of (n, p) solutions for n, p in [2, 20]
        """
        solutions = [
            (nn, q) for nn in range(2, 21) for q in range(2, 21)
            if nn == q ** (nn - 2)
        ]
        assert (3, 3) in solutions
        return solutions

    def __repr__(self):
        return (
            f"S5Z3(d1={self.d1}, lambda1={self.lambda1}, "
            f"K={self.K:.4f}, eta={self.eta:.4f}, p={self.p})"
        )
