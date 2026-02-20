"""
Pytest configuration for LENG tests.
Ensures falsification/, verification/, and project root are on Python path.
"""

import sys
from pathlib import Path

# Add project root (public-release) to path (tools/falsification -> parent.parent.parent)
project_root = Path(__file__).resolve().parent.parent.parent
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

# Add verification/ for standalone scripts (alpha_s_constraint, leng_replication, etc.)
verification = project_root / "verification"
if str(verification) not in sys.path:
    sys.path.insert(0, str(verification))

# Add falsification/ so standalone scripts can find leng_test_utils and dictionary_spec
falsification_dir = Path(__file__).resolve().parent
if str(falsification_dir) not in sys.path:
    sys.path.insert(0, str(falsification_dir))
