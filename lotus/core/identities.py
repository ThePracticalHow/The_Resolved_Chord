"""
lotus.core.identities — The Identity Chain
===============================================

The master identity chain connecting all sectors through p^n = 27.
"""


def identity_chain(M) -> dict:
    """Compute the full identity chain from the manifold.

    The chain:
      tau   = 1/p^n          = 1/27
      eta   = d1/p^n         = 2/9
      G     = lambda1 * eta  = 10/9
      c_grav = -tau/G        = -1/(d1*lambda1) = -1/30
      eta^2 = (p-1)*tau*K    [unique to (3,3)]

    Returns:
        dict of identity values
    """
    return {
        'tau': M.tau,
        'eta': M.eta,
        'G': M.G,
        'c_grav': M.c_grav,
        'eta_sq': M.eta ** 2,
        'eta_sq_identity': (M.p - 1) * M.tau * M.K,
        'pn': M.p ** M.n,
    }
