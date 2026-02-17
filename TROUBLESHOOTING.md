# Troubleshooting — LaTeX Conversion & Dependencies

This guide covers issues with `convert_tex_pdf_md.py` and its dependencies.

---

## Pandoc: "pandoc not found"

`convert_tex_pdf_md.py` uses **pandoc** for high-quality .tex→.md conversion.

**Common mistake:** `pip install pandoc` installs a Python *wrapper*, not the actual binary. You need the real executable.

---

## Quick Fixes

### Option 1: Winget (Windows 10/11)

```powershell
winget install JohnMacFarlane.Pandoc
```

Restart your terminal, then run `test_pandoc_install.bat`.

### Option 2: Chocolatey

```powershell
choco install pandoc
```

Restart your terminal.

### Option 3: Automated Script

Run `install_pandoc_simple.bat` — it tries winget and Chocolatey for you.

### Option 4: Direct Download (when winget/choco fail)

Run `install_pandoc_download.bat` — downloads from GitHub and extracts to `C:\pandoc\`. Then add to PATH (see "Add to PATH Manually" below).

### Option 5: Manual Download

1. Go to https://github.com/jgm/pandoc/releases/latest
2. Download **pandoc-*-windows-x86_64.zip** or **.msi**
3. **MSI:** Run it — it adds pandoc to PATH automatically
4. **ZIP:** Extract to `C:\pandoc\`, then add `C:\pandoc\` to your PATH
5. Restart your terminal

---

## Verify Installation

Run:

```
test_pandoc_install.bat
```

If it says `[OK] Pandoc is in PATH`, you're done.

---

## If Pandoc Still Won't Install

**Use the fallback:** The convert script uses **pdf→md** (pymupdf) when a `.pdf` exists. So:

1. Ensure you have `.pdf` files (from pdflatex or copied from arxiv_submission)
2. Run `convert_tex_pdf_md.py` — it will use pdf→md instead of tex→md
3. Quality is lower (no LaTeX structure) but it works without pandoc

**Online alternative:** https://pandoc.org/try/ — upload `.tex`, convert to Markdown, download.

---

## Add to PATH Manually

1. Win+R → `sysdm.cpl` → Enter
2. Advanced → Environment Variables
3. Under "User variables" or "System variables", select **Path** → Edit
4. New → `C:\pandoc` (or wherever you extracted pandoc)
5. OK → OK → OK
6. **Close and reopen** your terminal

---

## Diagnostic Scripts

| Script | Purpose |
|--------|---------|
| `test_pandoc_install.bat` | Check if pandoc is installed and in PATH |
| `find_pandoc.bat` | Search for pandoc anywhere on the system |
| `install_pandoc_simple.bat` | Try winget/Chocolatey install |
| `install_pandoc_download.bat` | Direct download when winget/choco fail |

---

## pdflatex: "pdflatex not found"

The script uses **pdflatex** to compile .tex→.pdf when .pdf is missing.

**Install:**
- **MiKTeX:** https://miktex.org/download — run installer
- **TeX Live:** https://tug.org/texlive/
- **Winget:** `winget install MiKTeX.MiKTeX`

**If MiKTeX prompts for updates:** Open MiKTeX Console → Check for updates (or run once to dismiss).

**Without pdflatex:** The script will skip tex→pdf. If you have .pdf files (e.g. from arxiv_submission), it can still generate .md via pdf→md (pymupdf).

---

## PATH Not Updating

After installing pandoc or pdflatex:
1. **Close all terminals** (including Cursor/VS Code integrated terminal)
2. Open a new terminal
3. Run `test_pandoc_install.bat` to verify

If Python still can't find pandoc, the script tries `cmd /c where pandoc` and common install locations automatically.
