"""
Pytest configuration for LENG tests.
Ensures public-release root and verification/ are on Python path.
"""

import sys
from pathlib import Path

# Add project root (public-release) to path
project_root = Path(__file__).resolve().parent.parent
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))
# Add verification/ for scripts like leng_replication
verification = project_root / "verification"
if str(verification) not in sys.path:
    sys.path.insert(0, str(verification))
