# The Resolved Chord

**The Spectral Geometry of Everything: Standard Model, Gravity, and Cosmology from S⁵/Z₃**

*Jixiang Leng, February 2026*

---

## 🚀 Quick Start

**First time here?** Follow these steps:

1. **Install dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

2. **Install pandoc** (for LaTeX → Markdown):
   ```bash
   # Windows (run as Administrator):
   winget install JohnMacFarlane.Pandoc
   
   # Or download manually from: https://github.com/jgm/pandoc/releases/latest
   ```

3. **Install MiKTeX** (for LaTeX → PDF):
   - Download: https://miktex.org/download
   - Run MiKTeX Console → Updates → Check for updates

4. **Test installation:**
   ```bash
   python convert_tex_pdf_md.py --dry-run
   ```

5. **Generate missing formats:**
   ```bash
   python convert_tex_pdf_md.py
   ```

---

## What this is

A single compact manifold — the lens space S⁵/Z₃ — determines all 26 parameters of the Standard Model, the gravitational coupling (0.10%), the cosmological constant (1.4%), and the full cosmological history, with zero free parameters.

The electron mass is the unit of measurement, not an input ([Duff, Okun, Veneziano 2002](https://arxiv.org/abs/physics/0110060)). Every dimensionless ratio is fixed by the spectral geometry.

## Key results

| Prediction | Precision | Status | Verification script |
|-----------|-----------|--------|---------------------|
| Koide ratio K = 2/3 | Exact | Theorem | `leng_replication.py` |
| Generations N_g = 3 | Exact | Theorem | `lotus_aps_generation.py` |
| Proton mass 6π⁵ | 10⁻¹¹ | Theorem (Parseval) | `ghost_parseval_proof.py` |
| Fine-structure 1/α = 137.038 | 0.001% | Theorem (APS lag) | `alpha_lag_proof.py` |
| Strong coupling α_s = 0.1187 | 0.56% | Theorem (ghost splitting) | `alpha_s_theorem.py` |
| Higgs VEV | 0.004% | Theorem (α cascade) | `vev_overlap.py` |
| Higgs mass | 0.036% | Theorem (α cascade) | `spectral_action_derivation.py` |
| All 6 quark masses | 0.01–5% | Theorem (ordering) | `downtype_spectral_ordering.py` |
| Full CKM matrix (9 elements) | 0.00–2% | Theorem | `ckm_complete.py` |
| Planck mass | 0.10% | Theorem (5-lock) | `gravity_theorem_proof.py` |
| Cosmological constant | 1.4% | Theorem | `cc_aps_proof.py` |
| Strong CP θ̄ = 0 | Exact | Theorem | — |
| Inflation n_s = 0.968 | 0.8σ | Theorem | `geometric_unification.py` |
| Baryogenesis η_B = 6.3×10⁻¹⁰ | 3% | Theorem | `geometric_unification.py` |
| DM abundance Ω_DM/Ω_B = 16/3 | 0.5% | Theorem | `geometric_unification.py` |

**43 predictions. All 43 at Theorem level. Zero Derived. Zero gaps. Zero free parameters.**

## Quick verification

```bash
pip install -r requirements.txt
# Run all core verification scripts at once:
python run_verification.py

# Or run individually:
python verification/EtaInvariant.py          # eta = 2/9
python verification/ghost_parseval_proof.py   # proton mass to 50 digits
python verification/gravity_theorem_proof.py  # Planck mass (16/16 checks)
python verification/alpha_lag_proof.py        # 1/alpha from APS lag correction
python verification/alpha_s_theorem.py        # alpha_s from ghost splitting
python verification/ckm_complete.py           # full CKM matrix from geometry
python verification/cc_aps_proof.py           # cosmological constant
python verification/vev_stiffness_proof.py   # Higgs VEV theorem (0.005%)
python verification/quantum_gravity_lotus.py  # QG dissolution (5 proofs)
python verification/fold_potential_paper.py   # V(phi) = lambda_H v^4 (phi^2-1)^2/4 (P14-P16)
python verification/vev_overlap_paper.py      # v/m_p overlap interpretation (P14)
python -m pytest falsification/ -v           # full test suite (includes verification scripts)
python -m pytest falsification/ -m "not slow"  # skip slow verification scripts (~faster)
```

## Paper

- **paper/The_Resolved_Chord_v10.tex** — Main paper: The Theorem of Everything (43/43 Theorem, 18+ sections)
- **supplements/Supplement_XII_CompanionGuide.tex** — Topic-based reference (Hurricanes, Equations, LOTUS, Predictions, Firewall, Strange Castles)
- **supplements/** — Supplements I–XII (complete proof chains). X = math-to-physics map, XI = derivation status firewall, XII = companion guide.
- **math-papers/** — Four standalone math papers

**Markdown versions:** Run `python convert_tex_pdf_md.py` to generate `.md` (and `.pdf` when missing) from `.tex`. Uses pandoc for tex→md (best quality) or pymupdf for pdf→md when pandoc is unavailable. Generates missing files; when source (.tex/.pdf) is newer than output, regenerates automatically. Use `--force` to regenerate all. If pandoc is missing: run `install_pandoc_simple.bat` or see `TROUBLESHOOTING.md`.

**Sync from other copies:** If you have multiple `public-release*` folders (e.g. `public-release 2-16-26 4-13`), run `python convert_tex_pdf_md.py --sync` to copy the newest version of each file (by modification time) from `paper/`, `supplements/`, and `math-papers/` into this folder, then convert. Use `--sync-only` to only sync without converting.

**Parallel processing:** For faster conversion of many files, use `python convert_tex_pdf_md.py --parallel`.

## ⚙️ Configuration

The `config.py` file contains all settings. Edit it to customize:

- **Pandoc options:** Change LaTeX→Markdown conversion settings
- **Parallel processing:** Adjust worker count and timeouts
- **Backup settings:** Control automatic backups before conversion
- **Logging:** Set log levels and file locations

Example `config.py` snippet:
```python
# Conversion settings
PANDOC_OPTIONS = [
    '--from=latex',
    '--to=markdown',
    '--wrap=none',
    '--extract-media=.',
]

# Performance
MAX_WORKERS = 4  # For parallel processing
TIMEOUT = 300    # Seconds per file

# Logging
LOG_LEVEL = 'INFO'
LOG_FILE = 'conversion.log'
```

## ✨ New Features

- **Parallel processing:** `--parallel` flag for faster batch conversion
- **Configuration system:** Customizable settings in `config.py`
- **Structured logging:** Detailed logs with `--verbose` flag
- **Automatic backups:** Backup files before conversion (configurable)
- **Sync from multiple folders:** `--sync` flag to merge from other `public-release*` folders
- **Dry run mode:** `--dry-run` to preview changes without executing

## 🔧 Troubleshooting

### Common Issues:

**"pandoc not found"**
```bash
# Install pandoc:
winget install JohnMacFarlane.Pandoc

# Or download manually:
# https://github.com/jgm/pandoc/releases/latest
# Extract to C:\pandoc and add to PATH
```

**"pdflatex not found"**
- Install MiKTeX: https://miktex.org/download
- Run MiKTeX Console → Updates → Check for updates
- Alternative: Use online LaTeX compilers (Overleaf, Papeeria)

**Python packages missing**
```bash
pip install numpy scipy mpmath pymupdf tqdm
```

**Permission errors**
- Run terminal as Administrator
- Or use online converters for testing

### Debug Commands:
```bash
# Test pandoc installation (run directly, not with python)
test_pandoc_install.bat

# Test LaTeX compilation (from public-release/)
pdflatex -interaction=nonstopmode paper/The_Resolved_Chord_v10.tex

# Verbose conversion
python convert_tex_pdf_md.py -v

# Dry run (see what would be done)
python convert_tex_pdf_md.py --dry-run
```

See `TROUBLESHOOTING.md` for detailed solutions.

## How to read this

1. Start with the abstract in **paper/The_Resolved_Chord_v10.tex** (the story)
2. For a quick overview of any topic, read **Supplement XII** (6 chapters: Hurricanes, Equations, LOTUS, Predictions, Derivation Status, Strange Castles)
3. For the explicit math-to-physics mapping, read **supplements/Supplement_X_SpectralAction.tex** (every chain from Tr(f(D^2)) to every prediction)
4. If you doubt a claim, check **supplements/Supplement_XI_DerivationStatus.tex** (the firewall)
5. Run the verification script for any prediction you want to check
6. Run `pytest falsification/` for the full adversarial test suite

## Citation

```bibtex
@article{leng2026resolved,
  author = {Leng, Jixiang},
  title = {The Resolved Chord: The Theorem of Everything},
  year = {2026},
  doi = {10.5281/zenodo.18655472},
  url = {https://doi.org/10.5281/zenodo.18655472}
}
```

## Tools in This Folder

| File | Purpose |
|------|---------|
| `run_verification.py` | Run all core verification scripts (run `python run_verification.py`) |
| `run_verification.bat` | **Double-click** to run verification; writes `verification_results.txt` |
| `run_falsification.py` | Run falsification test suite (run `python run_falsification.py`) |
| `run_falsification.bat` | **Double-click** to run falsification; writes `falsification_results.txt` |
| `convert_tex_pdf_md.py` | Generate .md and .pdf from .tex (run `python convert_tex_pdf_md.py`) |
| `test_pandoc_install.bat` | Check if pandoc is installed |
| `install_pandoc_simple.bat` | Install pandoc (winget/Chocolatey) |
| `install_pandoc_download.bat` | Direct download (when above fails) |
| `find_pandoc.bat` | Locate pandoc on system |
| `debug_pandoc.py` | Python pandoc diagnostic |
| `TROUBLESHOOTING.md` | Pandoc, pdflatex, PATH issues |

## License

MIT
