#!/usr/bin/env python3
"""
Test script for new LOTUS functionality
"""
import sys
import os

# Add the project path (go up from tools/tests/ to project root)
project_path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, project_path)

def test_lotus_functionality():
    print("Testing new LOTUS functionality...")
    print("=" * 50)

    try:
        # Test basic import
        import lotus
        print("✓ LOTUS import successful")

        # Test Universe creation
        u = lotus.Universe()
        print(f"✓ Universe created with {len(u.predictions)} predictions")

        # Test new derive method
        result = u.derive('F=ma')
        print(f"✓ Universe.derive() works: {result['name']}")

        # Test equations module
        from lotus.equations import derive as eq_derive
        result2 = eq_derive('E=mc2')
        print(f"✓ equations.derive() works: {result2['name']}")

        # Test CLI-like functionality
        print("\nTesting CLI functionality:")
        print("Available equations:")
        for eq in lotus.equations.all_equations()[:5]:
            print(f"  - {eq['name']}")

        print("\n✓ All tests passed!")

    except Exception as e:
        print(f"✗ Error: {e}")
        import traceback
        traceback.print_exc()
        return False

    return True

if __name__ == "__main__":
    test_lotus_functionality()