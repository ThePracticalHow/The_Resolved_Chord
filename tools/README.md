# Tools

Support utilities, tests, and build scripts for The Resolved Chord.

## Subfolders

| Folder | Purpose |
|--------|---------|
| **tests/** | Basic tests: `test_docker.py` (Docker container validation), `test_new_functionality.py` (LOTUS API smoke test) |
| **falsification/** | Adversarial pytest suite (~20 test files). Run: `python tools/falsification/run_falsification.py` or `pytest tools/falsification/ -v` |

## Scripts (root of tools/)

| Script | Purpose |
|--------|---------|
| `convert_tex_pdf_md.py` | Generate .md and .pdf from .tex sources |
| `validate_docs.py` | Validate documentation consistency (run by RUN_EVERYTHING.bat) |
| `release_gate.py` | Robust release checks (required files, stale paths, forbidden absolute links, local link integrity) |
| `auto_clean.py` | Cross-platform cleanup of cache/build artifacts (`--dry-run` by default, `--apply` to delete) |
| `handoff_template.md` | Standardized AI handoff format for cross-agent collaboration |
| `check_sync.py` | Compare papers sync status across formats |
| `config.py` | Conversion settings for LaTeX/Markdown |
| `colab_demo.py` | Simplified Colab-compatible demo |
| `pandoc_utils.py` | Pandoc path detection and utilities |
| `find_pandoc.bat` | Locate Pandoc installation |
| `install_pandoc_simple.bat` | Install Pandoc via winget |
| `install_pandoc_download.bat` | Download and install Pandoc |
| `test_pandoc_install.bat` | Verify Pandoc works |

## Running Tests

```bash
# From project root (public-release/):
python tools/falsification/run_falsification.py -v     # Full falsification suite
python -m pytest tools/falsification/ -v               # Same via pytest
python tools/tests/test_docker.py                     # Docker container test (run inside container)
python tools/tests/test_new_functionality.py           # Quick LOTUS API test
python tools/release_gate.py                           # Robustness release checks
python tools/auto_clean.py                             # Dry-run cleanup preview
python tools/auto_clean.py --apply                     # Apply cleanup
```
