"""
lotus.math.spectral_math — Pure Python Spectral Mathematics
===========================================================

Self-contained implementations of mathematical functions needed for LOTUS.
No external dependencies beyond Python stdlib (math, cmath, decimal).
"""

import math
import cmath
from decimal import Decimal, getcontext
from typing import Union, Tuple, List

# Set high precision for spectral calculations
getcontext().prec = 50

# Type aliases
Real = Union[int, float, Decimal]
Complex = Union[complex, complex]  # cmath complex

class SpectralMath:
    """Pure Python implementations for spectral geometry calculations."""

    @staticmethod
    def zeta(s: Real, terms: int = 1000) -> Decimal:
        """
        Riemann zeta function via series expansion.

        zeta(s) = sum_{n=1}^∞ n^{-s} for Re(s) > 1

        Args:
            s: Complex argument
            terms: Number of terms in series (default 1000)

        Returns:
            Decimal approximation of zeta(s)
        """
        if isinstance(s, complex):
            # Use complex arithmetic
            result = complex(0, 0)
            for n in range(1, terms + 1):
                result += 1 / (n ** s)
            return Decimal(str(result.real)) + Decimal(str(result.imag)) * 1j
        else:
            # Real case
            result = Decimal(0)
            for n in range(1, terms + 1):
                result += Decimal(1) / (Decimal(n) ** Decimal(s))
            return result

    @staticmethod
    def eta_invariant(manifold_dim: int, orbifold_order: int) -> Decimal:
        """
        Eta invariant for lens spaces L(p; 1,...,1) in C^n.

        For S^{2n+1}/Z_p, eta = (1/p) * sum_{k=1}^{p-1} omega^k * (i*cot(π*k/p))^n
        where omega = exp(2πi/p)

        Args:
            manifold_dim: Dimension of the sphere (2n+1)
            orbifold_order: Order of the orbifold group Z_p

        Returns:
            Decimal value of the eta invariant
        """
        if manifold_dim % 2 != 1:
            raise ValueError("Manifold dimension must be odd (sphere dimension)")

        n = (manifold_dim - 1) // 2  # Complex dimension
        p = orbifold_order

        omega = cmath.exp(2j * math.pi / p)
        total = complex(0, 0)

        for k in range(1, p):
            cot_k = math.cos(math.pi * k / p) / math.sin(math.pi * k / p)
            total += omega ** k * (1j * cot_k) ** n

        result = total / p
        return result

    @staticmethod
    def numerical_integrate(
        func: callable,
        a: Real,
        b: Real,
        method: str = 'simpson',
        n: int = 1000
    ) -> Decimal:
        """
        Pure Python numerical integration.

        Args:
            func: Function to integrate
            a, b: Integration limits
            method: 'simpson' or 'trapezoid'
            n: Number of intervals

        Returns:
            Decimal approximation of the integral
        """
        if method == 'simpson':
            return SpectralMath._simpson_rule(func, a, b, n)
        elif method == 'trapezoid':
            return SpectralMath._trapezoid_rule(func, a, b, n)
        else:
            raise ValueError("Method must be 'simpson' or 'trapezoid'")

    @staticmethod
    def _simpson_rule(func: callable, a: Real, b: Real, n: int) -> Decimal:
        """Simpson's rule integration."""
        if n % 2 != 0:
            n += 1  # Ensure even number of intervals

        h = Decimal(str(b - a)) / Decimal(str(n))
        result = Decimal(0)

        # Endpoints
        result += Decimal(str(func(a)))
        result += Decimal(str(func(b)))

        # Interior points
        for i in range(1, n):
            x = a + i * float(h)
            if i % 2 == 0:
                result += 2 * Decimal(str(func(x)))
            else:
                result += 4 * Decimal(str(func(x)))

        result *= h / Decimal(3)
        return result

    @staticmethod
    def _trapezoid_rule(func: callable, a: Real, b: Real, n: int) -> Decimal:
        """Trapezoid rule integration."""
        h = Decimal(str(b - a)) / Decimal(str(n))
        result = Decimal(0)

        # Endpoints
        result += Decimal(str(func(a))) / 2
        result += Decimal(str(func(b))) / 2

        # Interior points
        for i in range(1, n):
            x = a + i * float(h)
            result += Decimal(str(func(x)))

        result *= h
        return result

    @staticmethod
    def gamma(z: Real) -> Decimal:
        """
        Gamma function using Lanczos approximation.

        Args:
            z: Argument

        Returns:
            Decimal approximation of gamma(z)
        """
        # Lanczos coefficients
        g = 7
        p = [0.99999999999980993, 676.5203681218851, -1259.1392167224028,
             771.32342877765313, -176.61502916214059, 12.507343278686905,
             -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7]

        if isinstance(z, complex):
            z = complex(z)
            if z.real < 0.5:
                return math.pi / (math.sin(math.pi * z) * SpectralMath.gamma(1 - z))

            z -= 1
            x = p[0]
            for i in range(1, len(p)):
                x += p[i] / (z + i)
            t = z + g + 0.5
            result = math.sqrt(2 * math.pi) * (t ** (z + 0.5)) * math.exp(-t) * x
            return Decimal(str(result.real)) + Decimal(str(result.imag)) * 1j
        else:
            z = float(z)
            if z < 0.5:
                return float(math.pi / (math.sin(math.pi * z) * SpectralMath.gamma(1 - z)))

            z -= 1
            x = p[0]
            for i in range(1, len(p)):
                x += p[i] / (z + i)
            t = z + g + 0.5
            return Decimal(str(math.sqrt(2 * math.pi) * (t ** (z + 0.5)) * math.exp(-t) * x))

    @staticmethod
    def bessel_j(n: int, x: Real) -> Decimal:
        """
        Bessel function J_n(x) using series expansion.

        Args:
            n: Order
            x: Argument

        Returns:
            Decimal approximation of J_n(x)
        """
        if abs(x) < 1e-10:
            return Decimal(1) if n == 0 else Decimal(0)

        # Series expansion for small x
        if abs(x) < 10:
            result = Decimal(0)
            term = Decimal(1)
            k = 0
            while abs(term) > 1e-20:
                if k >= 100:  # Prevent infinite loop
                    break
                result += term
                term *= -Decimal(str(x**2)) / (4 * (k + 1) * (k + n + 1))
                k += 1
            result *= (Decimal(str(x)) / 2) ** n
            return result
        else:
            # For large x, use asymptotic expansion (simplified)
            # This is a placeholder - full implementation would need more terms
            return Decimal(str(math.sqrt(2 / (math.pi * abs(x))) * math.cos(abs(x) - math.pi*(n/2 + 1/4))))

    @staticmethod
    def complex_sqrt(z: Complex) -> Complex:
        """Principal square root of a complex number."""
        return cmath.sqrt(z)

    @staticmethod
    def complex_exp(z: Complex) -> Complex:
        """Complex exponential."""
        return cmath.exp(z)

    @staticmethod
    def complex_log(z: Complex) -> Complex:
        """Principal logarithm of a complex number."""
        return cmath.log(z)

# Convenience functions
def zeta(s: Real, terms: int = 1000) -> Decimal:
    """Riemann zeta function."""
    return SpectralMath.zeta(s, terms)

def eta_invariant(manifold_dim: int, orbifold_order: int) -> Decimal:
    """Eta invariant for lens spaces."""
    return SpectralMath.eta_invariant(manifold_dim, orbifold_order)

def numerical_integrate(func: callable, a: Real, b: Real, method: str = 'simpson', n: int = 1000) -> Decimal:
    """Numerical integration."""
    return SpectralMath.numerical_integrate(func, a, b, method, n)

def gamma(z: Real) -> Decimal:
    """Gamma function."""
    return SpectralMath.gamma(z)

def bessel_j(n: int, x: Real) -> Decimal:
    """Bessel function J_n(x)."""
    return SpectralMath.bessel_j(n, x)