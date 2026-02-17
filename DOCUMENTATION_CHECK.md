# Documentation Check Report

**Date:** February 2026  
**Updated:** February 16, 2026

## Summary

| Category | Count |
|----------|-------|
| Critical (wrong/misleading) | 2 |
| Minor (outdated/incomplete) | 5 |
| OK | — |

---

## Issues Found (and Fixed)

### Critical — FIXED

1. ~~**PUBLICATION_INDEX.md** line 9: Claimed "8 figures" correctly — paper uses 8 figures~~ → Fixed figures/README.md to document all 8 figure generation scripts.
2. ~~**README.md** line 172: `python test_pandoc_install.bat`~~ → Fixed: `test_pandoc_install.bat` (run directly).
3. ~~**README.md** config example~~ → Fixed: aligned with actual config.py structure.

### Minor — FIXED

4. ~~**README.md** line 100: "only generates missing"~~ → Fixed: added --force, mtime-based updates.
5. ~~**GETTING_STARTED.md** folder layout~~ → Fixed: added config.py, check_sync.py.
6. ~~**PUBLICATION_INDEX.md** "resolved-chord"~~ → Fixed: "public-release".
7. ~~**README.md** pip install~~ → Fixed: `pip install -r requirements.txt`.
8. ~~**Conversion tools**~~ → Fixed: --force and check_sync.py in PUBLICATION_INDEX.

---

## Verified OK

- All verification scripts listed in PUBLICATION_INDEX exist
- All listed .bat files exist
- run_verification.bat, run_falsification.bat documented
- companion_guide 6 chapters: consistent across docs
- figures/README now includes all 8 figure scripts (updated 2026-02-16)
- TROUBLESHOOTING.md: correct (test_pandoc_install.bat without python)
- Citation consistent (README, GETTING_STARTED)
- Main paper lotus_emblem figure now has proper label (updated 2026-02-16)
