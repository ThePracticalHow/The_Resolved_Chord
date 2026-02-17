# Getting Started — Public Release

**Open this folder.** Everything you need to read, verify, and work with The Resolved Chord is here.

---

## 1. What You Need

| Requirement | Purpose | How to get it |
|-------------|---------|---------------|
| **Python 3.8+** | Run verification and conversion | [python.org](https://python.org) |
| **numpy, scipy, mpmath** | Physics calculations | `pip install -r requirements.txt` |
| **pymupdf** | PDF→Markdown (included in requirements) | `pip install pymupdf` |
| **pandoc** *(optional)* | LaTeX→Markdown (best quality) | See [TROUBLESHOOTING.md](TROUBLESHOOTING.md) |
| **pdflatex** *(optional)* | Compile .tex to .pdf | MiKTeX or TeX Live |

---

## 2. Quick Start (5 minutes)

```bash
# From this folder (public-release/)
pip install -r requirements.txt
python run_verification.py                   # run all core verification scripts
python -m pytest falsification/ -v           # full test suite (includes verification)
python -m pytest falsification/ -m "not slow" # skip slow scripts (~faster)
```

---

## 3. Folder Layout

```
public-release/
├── GETTING_STARTED.md     ← You are here
├── README.md              ← Paper overview, key results
├── PUBLICATION_INDEX.md   ← Full file index
├── requirements.txt       ← pip install -r requirements.txt
├── config.py              ← Conversion settings (optional)
├── check_sync.py          ← Compare tex/pdf/md sync status
│
├── paper/                 ← Main paper (LaTeX + Markdown)
├── supplements/           ← Supplements I–XII (proof chains + companion guide)
├── math-papers/           ← Four standalone math papers
│
├── verification/          ← Scripts that verify each prediction
├── falsification/         ← Adversarial test suite (pytest)
│
├── run_verification.py    ← Run all core verification scripts
├── run_verification.bat   ← Double-click to run verification
├── run_falsification.py   ← Run falsification test suite
├── run_falsification.bat  ← Double-click to run falsification
├── convert_tex_pdf_md.py   ← Generate .md and .pdf from .tex
├── TROUBLESHOOTING.md     ← Pandoc, pdflatex, PATH issues
│
└── Tools (LaTeX conversion):
    ├── test_pandoc_install.bat    ← Check if pandoc works
    ├── install_pandoc_simple.bat   ← Install pandoc (winget/Chocolatey)
    ├── install_pandoc_download.bat ← Direct download (when above fails)
    ├── find_pandoc.bat             ← Locate pandoc on system
    └── debug_pandoc.py             ← Python pandoc diagnostic
```

---

## 4. Generate Markdown from LaTeX

```bash
python convert_tex_pdf_md.py
```

- **With pandoc:** Best quality (preserves structure)
- **Without pandoc:** Uses pymupdf (pdf→md) when .pdf exists — works but lower quality
- **Options:** `--sync` to pull newest from sibling folders, `--dry-run` to preview

If pandoc is missing: run `install_pandoc_simple.bat` or see [TROUBLESHOOTING.md](TROUBLESHOOTING.md).

---

## 5. Read the Paper

1. **Start:** `paper/The_Resolved_Chord_v10.tex` (or `.md` if generated)
2. **Doubt a claim?** Check `supplements/Supplement_XI_DerivationStatus.tex` — every claim, proof location, verification script
3. **Run verification:** Each prediction has a script in `verification/`

---

## 6. Testing

| Command | Purpose |
|---------|---------|
| **`run_verification.bat`** | **Double-click** to run verification; writes `verification_results.txt` |
| **`run_falsification.bat`** | **Double-click** to run falsification; writes `falsification_results.txt` |
| `python run_verification.py` | Run verification scripts directly |
| `python run_verification.py --quick` | Run only fast verification scripts |
| `python -m pytest falsification/ -v` | Full test suite |
| `python -m pytest falsification/ -m "not slow"` | Skip slow verification scripts |
| `python -m pytest falsification/test_paper_claims.py -v` | Paper formula audit only |

---

## 7. Troubleshooting

| Issue | Solution |
|-------|----------|
| `pandoc not found` | [TROUBLESHOOTING.md](TROUBLESHOOTING.md) — install binary, not `pip install pandoc` |
| `pdflatex not found` | Install MiKTeX or TeX Live; or use existing .pdf files |
| `pymupdf` missing | `pip install pymupdf` |
| PATH not updated | Restart terminal after installing pandoc/pdflatex |
| Tests fail | Run from `public-release/`; ensure `pip install -r requirements.txt` |

---

## 8. Citation

```bibtex
@article{leng2026resolved,
  author = {Leng, Jixiang},
  title = {The Resolved Chord: The Theorem of Everything},
  year = {2026},
  doi = {10.5281/zenodo.18655472},
  url = {https://doi.org/10.5281/zenodo.18655472}
}
```
