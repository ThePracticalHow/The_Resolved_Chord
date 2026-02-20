#!/usr/bin/env python3
"""
PMNS from Point/Side/Face Self-Overlap
=======================================

Key insight (J. Leng, Feb 2026): neutrinos are NOT three copies of one thing.
They are three DIFFERENT geometric measurements of the Z_3 fold:
  nu_1 = POINT (cone tip, 0-dimensional)
  nu_2 = SIDE  (fold wall, codimension-1)
  nu_3 = FACE  (sector bulk, full-dimensional)

This explains:
  - Why PMNS angles are large (different objects overlap strongly)
  - Why no Koide ratio (not three copies of same thing)
  - Why mass hierarchy m1 << m2 << m3 (dimensional hierarchy)

The tunneling overlap matrix T has ALL entries derived from spectral data:

  Diagonal (self-overlap):
    T_11 = +2*sigma/p     = +4/11   (point: 2 polarizations, p sectors)
    T_22 = -sigma/p       = -2/11   (side: wall depletion, 1/p walls)
    T_33 = -(d1-2)*r      = -64/729 (face: bulk leaks via 2(n-1) channels)

  Off-diagonal (coupling between objects):
    T_12 = -eta           = -2/9    (point<->side: Donnelly fold bleed)
    T_23 = -K*sigma       = -4/11   (side<->face: harmonic lock x impedance)
    T_13 = +eta*K         = +4/27   (point<->face: full fold-wall amplitude)

  Perturbation strength:
    epsilon = sqrt(d1 + lam1) = sqrt(11)  (total ell=1 spectral content)

Result: M = D + sqrt(11)*T gives all three PMNS angles within 2% of PDG
with ZERO fitted parameters.

Jixiang Leng & Claude, February 2026
"""

import numpy as np
from typing import Tuple

# ====================================================================
# SPECTRAL INVARIANTS (fixed input from S^5/Z_3)
# ====================================================================
ETA = 2 / 9          # Donnelly eta invariant
K = 2 / 3            # Koide ratio
D1 = 6               # first-excited degeneracy
LAMBDA_1 = 5         # first non-trivial eigenvalue
P = 3                # orbifold order
N_CX = 3             # complex dimensions

# Derived
SIGMA = D1 / (D1 + LAMBDA_1)  # 6/11 = spectral impedance
R = (ETA * K) ** 2             # 16/729 = reactor amplitude (two-wall)
EK = ETA * K                   # 4/27  = single fold-wall amplitude
EPSILON = np.sqrt(D1 + LAMBDA_1)  # sqrt(11) = total ell=1 spectral content

# PDG 2024 (NuFIT 5.2, normal ordering)
SIN2_13_PDG = 0.02200
SIN2_12_PDG = 0.307
SIN2_23_PDG = 0.546
DELTA_CP_PDG_DEG = 195


def democratic_matrix() -> np.ndarray:
    """D_ij = 1/3 for all i,j. Z3 symmetric, TBM zeroth order."""
    return np.ones((3, 3)) / 3


def point_side_face_T(
    eta: float = ETA,
    k: float = K,
    sigma: float = SIGMA,
    d1: int = D1,
    p: int = P,
    include_cp: bool = True,
) -> np.ndarray:
    """
    Tunneling overlap matrix from point/side/face geometry.

    Each neutrino is a DIFFERENT geometric object:
      nu_1 = point (cone tip)    -- 0-dimensional
      nu_2 = side (fold wall)    -- codimension-1
      nu_3 = face (sector bulk)  -- full-dimensional

    Returns complex Hermitian 3x3 matrix.
    """
    r = (eta * k) ** 2
    ek = eta * k

    # Diagonal: self-overlap of each geometric object
    T11 = 2 * sigma / p        # point: 2 polarizations, shared among p sectors
    T22 = -sigma / p            # side: wall depletion (1 of p walls)
    T33 = -(d1 - 2) * r        # face: bulk leaks via 2(n-1) transverse channels

    # Off-diagonal: coupling between different objects
    T12 = -eta                  # point<->side: fold bleed (Donnelly)
    T23 = -k * sigma            # side<->face: harmonic lock x impedance
    T13 = ek                    # point<->face: fold-wall amplitude (positive: 2 hops)

    T = np.array([
        [T11, T12, T13],
        [T12, T22, T23],
        [T13, T23, T33],
    ], dtype=complex)

    if include_cp:
        # CP phase: Z3 accumulation, each wall contributes delta/3
        delta_per_wall = np.arctan(2 * np.pi ** 2 / 9) / 3
        cp_amp = eta * np.sin(delta_per_wall)
        T[0, 1] += 1j * cp_amp; T[1, 0] -= 1j * cp_amp
        T[1, 2] += 1j * cp_amp; T[2, 1] -= 1j * cp_amp
        T[0, 2] -= 1j * cp_amp; T[2, 0] += 1j * cp_amp

    return T


def extract_pmns(U: np.ndarray) -> Tuple[float, float, float, float]:
    """Extract sin^2 angles and Jarlskog from PMNS-like unitary matrix."""
    s13 = np.abs(U[0, 2])
    c13 = np.sqrt(max(0, 1 - s13 ** 2))
    if c13 > 1e-10:
        s12 = np.abs(U[0, 1]) / c13
        s23 = np.abs(U[1, 2]) / c13
    else:
        s12, s23 = 0, 0
    J = np.imag(U[0, 0] * np.conj(U[0, 1]) * np.conj(U[1, 0]) * U[1, 1])
    return s13 ** 2, s12 ** 2, s23 ** 2, J


def main():
    print("=" * 72)
    print("  FORWARD PMNS DERIVATION: ZERO FITTED PARAMETERS")
    print("  Point / Side / Face Tunneling Overlap")
    print("=" * 72)

    D = democratic_matrix()
    T = point_side_face_T()

    print(f"\n  Spectral invariants: d1={D1}, lam1={LAMBDA_1}, K={K:.4f}, "
          f"eta={ETA:.4f}, p={P}")
    print(f"  Derived: sigma={SIGMA:.6f}, r={R:.6f}, eK={EK:.6f}")
    print(f"  Epsilon = sqrt(d1+lam1) = sqrt({D1+LAMBDA_1}) = {EPSILON:.6f}")

    # Build mass matrix
    M = D + EPSILON * T

    # Diagonalize
    eigenvalues, eigenvectors = np.linalg.eigh(M)
    idx = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[idx]
    U = eigenvectors[:, idx]

    sin2_13, sin2_12, sin2_23, J = extract_pmns(U)

    print(f"\n  T matrix (real part):")
    for i in range(3):
        row = "    ["
        for j in range(3):
            row += f" {np.real(T[i,j]):8.5f}"
        row += " ]"
        print(row)

    print(f"\n  M = D + sqrt(11)*T:")
    for i in range(3):
        row = "    ["
        for j in range(3):
            row += f" {np.real(M[i,j]):8.5f}"
        row += " ]"
        print(row)

    # Results
    print("\n" + "-" * 72)
    print("  RESULTS (fitted parameters: ZERO)")
    print("-" * 72)

    data = [
        ("sin2 theta_13", sin2_13, SIN2_13_PDG),
        ("sin2 theta_12", sin2_12, SIN2_12_PDG),
        ("sin2 theta_23", sin2_23, SIN2_23_PDG),
    ]

    print(f"\n  {'Parameter':<22} {'Predicted':>10} {'PDG':>10} {'Error':>10}")
    print(f"  {'-'*22} {'-'*10} {'-'*10} {'-'*10}")
    for name, pred, pdg in data:
        print(f"  {name:<22} {pred:10.6f} {pdg:10.6f} {(pred/pdg-1)*100:+9.2f}%")
    print(f"  {'Jarlskog J':<22} {J:10.6f} {'~0.031':>10}")

    print(f"\n  Eigenvalues: [{eigenvalues[0]:.4f}, {eigenvalues[1]:.4f}, "
          f"{eigenvalues[2]:.4f}]")

    # P24: Leptonic CP phase
    gamma_ckm = np.arctan(2 * np.pi**2 / 9)
    delta_pmns_formula = P * gamma_ckm
    delta_pmns_pdg = np.radians(DELTA_CP_PDG_DEG)
    print(f"\n  P24 Leptonic CP phase (P24 formula):")
    print(f"    delta_CP = p * arctan(2*pi^2/9) = {P} * {np.degrees(gamma_ckm):.2f} deg")
    print(f"    Predicted: {np.degrees(delta_pmns_formula):.2f} deg")
    print(f"    PDG central: {DELTA_CP_PDG_DEG} +/- 50 deg")
    print(f"    Error from central: {(np.degrees(delta_pmns_formula)/DELTA_CP_PDG_DEG - 1)*100:+.2f}%")

    # Compare to individual formulas
    print("\n" + "-" * 72)
    print("  COMPARISON: Forward T vs Individual Formulas")
    print("-" * 72)

    # P24 leptonic CP phase: neutrinos traverse all p=3 walls → delta = p * gamma
    gamma_ckm = np.arctan(2 * np.pi**2 / 9)   # CKM CP phase (1 wall)
    delta_pmns_formula = P * gamma_ckm          # PMNS CP phase (p walls)

    formulas = [
        ("sin2_13", (ETA * K) ** 2, sin2_13, SIN2_13_PDG),
        ("sin2_12", P / np.pi**2, sin2_12, SIN2_12_PDG),
        ("sin2_23", SIGMA, sin2_23, SIN2_23_PDG),
    ]

    print(f"\n  {'Angle':<12} {'Formula':>10} {'Forw T':>10} {'PDG':>10} "
          f"{'F err':>8} {'T err':>8}")
    print(f"  {'-'*12} {'-'*10} {'-'*10} {'-'*10} {'-'*8} {'-'*8}")
    for name, fval, tval, pdg in formulas:
        print(f"  {name:<12} {fval:10.6f} {tval:10.6f} {pdg:10.6f} "
              f"{(fval/pdg-1)*100:+7.2f}% {(tval/pdg-1)*100:+7.2f}%")

    # Physical interpretation
    print("\n" + "=" * 72)
    print("  DICTIONARY: POINT / SIDE / FACE")
    print("=" * 72)

    print("""
  T matrix structure:

  T = [ +2s/p    -eta     +eK  ]    sigma = d1/(d1+lam1) = 6/11
      [  -eta    -s/p    -K*s  ]    eta   = 2/9 (Donnelly)
      [  +eK    -K*s  -(d1-2)r ]    K     = 2/3 (Koide)
                                     eK    = eta*K = 4/27
  + i * eta*sin(delta/3) * antisymm  r     = (eK)^2 = 16/729

  Diagonal (SELF-OVERLAP of each geometric object):
    T_11 = +2sigma/p  = +4/11    POINT: 2 polarizations x impedance / p sectors
    T_22 = -sigma/p   = -2/11    SIDE:  wall depletion (1 of p walls)
    T_33 = -(d1-2)*r  = -64/729  FACE:  bulk leaks via 2(n-1) channels

  Off-diagonal (COUPLING between different objects):
    T_12 = -eta   = -2/9     point<->side: fold bleed (Donnelly invariant)
    T_23 = -K*s   = -4/11    side<->face:  harmonic lock x impedance
    T_13 = +eK    = +4/27    point<->face: fold-wall amplitude (2 hops => +)

  Perturbation strength:
    epsilon = sqrt(d1 + lam1) = sqrt(11): total ell=1 spectral content

  WHY NEUTRINOS HAVE NO KOIDE:
    Charged leptons = three copies of ONE thing (cone point mode)
    -> Circulant symmetry -> K = 2/3 exact
    Neutrinos = three DIFFERENT things (point/side/face)
    -> No circulant -> each has its own formula -> Q_nu != 2/3
""")


if __name__ == "__main__":
    main()
