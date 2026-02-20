#!/usr/bin/env python3
"""
Test script for LOTUS Docker container.
Run this to verify the containerized environment works correctly.
"""

import sys
import subprocess

def run_command(cmd, description):
    """Run a command and report success/failure."""
    print(f"Testing: {description}")
    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=30)
        if result.returncode == 0:
            print(f"✅ PASS: {description}")
            return True
        else:
            print(f"❌ FAIL: {description}")
            print(f"   Error: {result.stderr}")
            return False
    except subprocess.TimeoutExpired:
        print(f"⏰ TIMEOUT: {description}")
        return False
    except Exception as e:
        print(f"💥 ERROR: {description} - {e}")
        return False

def main():
    """Run all container tests."""
    print("🐳 Testing LOTUS Docker Container")
    print("=" * 40)

    tests_passed = 0
    total_tests = 0

    # Test 1: Import LOTUS
    total_tests += 1
    if run_command("python -c 'import lotus; print(\"LOTUS imported successfully\")'",
                   "LOTUS module import"):
        tests_passed += 1

    # Test 2: Create Universe
    total_tests += 1
    if run_command("python -c 'import lotus; u = lotus.Universe(); print(f\"Universe created with {len(u.predictions)} predictions\")'",
                   "Universe instantiation"):
        tests_passed += 1

    # Test 3: Compile universe
    total_tests += 1
    if run_command("python verification/compile_universe.py",
                   "Universe compilation"):
        tests_passed += 1

    # Test 4: Run basic verification
    total_tests += 1
    if run_command("python verification/EtaInvariant.py",
                   "Basic verification script"):
        tests_passed += 1

    # Test 5: Run falsification tests
    total_tests += 1
    if run_command("python -m pytest tools/falsification/test_universe.py -v",
                   "Basic falsification test"):
        tests_passed += 1

    print("\n" + "=" * 40)
    print(f"Results: {tests_passed}/{total_tests} tests passed")

    if tests_passed == total_tests:
        print("🎉 All tests passed! LOTUS container is working correctly.")
        return 0
    else:
        print("⚠️  Some tests failed. Check the output above for details.")
        return 1

if __name__ == "__main__":
    sys.exit(main())