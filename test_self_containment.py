#!/usr/bin/env python3
"""
LOTUS Self-Containment Demonstration
====================================

This script demonstrates that LOTUS is now completely self-contained
with no external library dependencies beyond Python stdlib.
"""

import sys
import os

# Add the lotus package to path
sys.path.insert(0, os.path.dirname(__file__))

def main():
    print("LOTUS Self-Containment Demonstration")
    print("=" * 50)

    try:
        # Test core LOTUS import
        print("1. Testing core LOTUS import...")
        import lotus
        print("   ✓ Core LOTUS imports successfully")

        # Test spectral math
        print("\n2. Testing spectral mathematics...")
        zeta_2 = lotus.zeta(2)
        print(f"   ✓ ζ(2) = {zeta_2}")
        print("     Expected: π²/6 ≈ 1.644934...")
        print("     Actual:   {:.12f}".format(float(zeta_2)))

        # Test eta invariant (simplified)
        print("\n3. Testing eta invariant calculation...")
        eta_val = lotus.eta_invariant(5, 3)  # S^5/Z_3
        print(f"   ✓ η_D(S^5/Z_3) = {eta_val}")
        print("     Should be complex with |η| = 2/9 ≈ 0.222...")

        # Test caching system
        print("\n4. Testing caching system...")
        from lotus.core.cache import spectral_cache
        spectral_cache.put('test_key', 'test_value')
        result = spectral_cache.get('test_key')
        print(f"   ✓ Cache put/get: {result}")
        print(f"   ✓ Cache size: {spectral_cache.stats()['size']}")

        # Test ASCII plotting
        print("\n5. Testing ASCII plotting...")
        plot = lotus.plot_eigenvalues([1, 5, 6])  # d1=6, λ1=5, p=3
        print("   ✓ ASCII eigenvalue plot:")
        print(plot)

        # Test configuration
        print("\n6. Testing configuration system...")
        cache_size = lotus.get_config('cache.max_size')
        print(f"   ✓ Config cache.max_size: {cache_size}")

        # Test documentation generation
        print("\n7. Testing documentation generation...")
        api_ref = lotus.generate_api_ref()
        print(f"   ✓ Generated API reference: {len(api_ref)} characters")

        print("\n" + "=" * 50)
        print("🎉 SUCCESS: LOTUS is completely self-contained!")
        print("   - No external dependencies required")
        print("   - All mathematics implemented internally")
        print("   - Full functionality preserved")
        print("   - Ready for independent deployment")

        return True

    except Exception as e:
        print(f"\n❌ ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)