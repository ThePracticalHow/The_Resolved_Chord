#!/usr/bin/env python3
"""
Convert .tex files in public-release: generate .md and .pdf from .tex sources.

Folder layout:
  papers/             ← Source .tex files (main paper + supplements at root)
  papers/math-papers/ ← Standalone math papers
  papers/pdf/         ← Generated .pdf files
  papers/markdown/    ← Generated .md files

Dependencies:
- pandoc: for tex->md (install: https://pandoc.org/installing.html)
- pdflatex: for tex->pdf (install TeX Live: https://tug.org/texlive/ or MiKTeX: https://miktex.org/)
- pymupdf: for pdf->md (pip install pymupdf)
- Optional: pytesseract for OCR fallback (pip install pytesseract; install tesseract-ocr)
- Optional: tqdm for progress bar (pip install tqdm)

Rules:
- Source of truth: .tex files in tex/
- Generated .md goes to markdown/ (mirrored subdir)
- Generated .pdf goes to pdf/ (mirrored subdir)
- If .md is newer than .tex: skip (do not overwrite manually edited .md)

Run from public-release/ or pass --root. Rerun anytime; only generates missing/stale files.
Automatically cleans up all backup files after each run.
"""

import argparse
import concurrent.futures
import logging
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path


try:
    from tqdm import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False

# pandoc_utils and config live in the same tools/ directory
_tools_dir = str(Path(__file__).resolve().parent)
if _tools_dir not in sys.path:
    sys.path.insert(0, _tools_dir)

try:
    from pandoc_utils import find_pandoc  # type: ignore[import-untyped]  # resolved via runtime sys.path
except ImportError:
    def find_pandoc():  # type: ignore[misc]
        """Fallback when pandoc_utils is unavailable; tex->md will be skipped."""
        return None

# Import configuration
try:
    from config import config, load_config_from_file
    load_config_from_file()  # Load user config if available
except ImportError:
    # Fallback if config.py doesn't exist
    class Config:
        pandoc_options = ["-f", "latex", "-t", "markdown", "--wrap=none", "--markdown-headings=atx"]
        pdflatex_options = ["-interaction=nonstopmode", "-halt-on-error"]
        max_workers = min(4, os.cpu_count() or 1)
        log_level = "INFO"
        log_file = None
        backup_existing = True
        clean_aux_files = True
        auto_clean_backups = True  # NEW: Auto-clean old backups
        backup_retention_days = 7  # NEW: Days to keep backups
        content_subdirs = ("math-papers",)
        tex_dir = "papers"
        markdown_dir = "papers/markdown"
        pdf_dir = "papers/pdf"
        content_exts = (".tex", ".pdf", ".md")
        pandoc_timeout = 300
        pdflatex_timeout = 120
        pdf_extract_timeout = 60
        enable_ocr = True
        ocr_lang = "eng"
        show_progress = True
        progress_desc = "Converting"
    config = Config()

CONTENT_SUBDIRS = config.content_subdirs
CONTENT_EXTS = config.content_exts

# Setup logging
def setup_logging():
    """Setup logging based on configuration."""
    log_level = getattr(logging, config.log_level.upper(), logging.INFO)

    # Create logger
    logger = logging.getLogger('convert_tex_pdf_md')
    logger.setLevel(log_level)

    # Remove existing handlers
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)

    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(log_level)
    console_formatter = logging.Formatter('%(levelname)s: %(message)s')
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)

    # File handler if configured
    if config.log_file:
        file_handler = logging.FileHandler(config.log_file)
        file_handler.setLevel(log_level)
        file_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(file_formatter)
        logger.addHandler(file_handler)

    return logger

logger = setup_logging()


def cleanup_old_backups(directory: Path, retention_days: int = 7) -> int:
    """
    Clean up old backup files based on retention policy.
    Returns number of files removed.
    """
    if retention_days <= 0:
        return 0  # Keep forever

    cutoff_time = time.time() - (retention_days * 24 * 60 * 60)
    removed_count = 0

    for backup_file in directory.glob("**/*.backup"):
        try:
            if backup_file.stat().st_mtime < cutoff_time:
                backup_file.unlink()
                logger.debug(f"Removed old backup: {backup_file}")
                removed_count += 1
        except (OSError, PermissionError) as e:
            logger.warning(f"Could not remove backup {backup_file}: {e}")

    if removed_count > 0:
        logger.info(f"Cleaned up {removed_count} old backup files")

    return removed_count


def cleanup_all_backups(directory: Path) -> int:
    """
    Clean up ALL backup files immediately, regardless of age.
    Returns number of files removed.
    """
    removed_count = 0

    for backup_file in directory.glob("**/*.backup"):
        try:
            backup_file.unlink()
            logger.debug(f"Removed backup: {backup_file}")
            removed_count += 1
        except (OSError, PermissionError) as e:
            logger.warning(f"Could not remove backup {backup_file}: {e}")

    if removed_count > 0:
        logger.info(f"Cleaned up {removed_count} backup files")

    return removed_count


def cleanup_aux_files(directory: Path) -> int:
    """
    Clean up ALL auxiliary LaTeX files (.aux, .log, .out, .fls, .fdb_latexmk, .synctex.gz, .toc).
    Returns number of files removed.
    """
    removed_count = 0
    aux_exts = ['.aux', '.log', '.out', '.fls', '.fdb_latexmk', '.synctex.gz', '.toc']

    for ext in aux_exts:
        for aux_file in directory.glob(f"**/*{ext}"):
            try:
                aux_file.unlink()
                logger.debug(f"Removed auxiliary file: {aux_file}")
                removed_count += 1
            except (OSError, PermissionError) as e:
                logger.warning(f"Could not remove auxiliary file {aux_file}: {e}")

    if removed_count > 0:
        logger.info(f"Cleaned up {removed_count} auxiliary files")

    return removed_count


def perform_cleanup(root: Path, logger) -> int:
    """
    Perform cleanup of backup and auxiliary files across all content directories.
    Returns total number of files removed.
    """
    logger.info("Cleaning up all backup and auxiliary files...")
    total_cleaned = 0

    # Clean root directory too
    cleaned = cleanup_all_backups(root)
    total_cleaned += cleaned
    cleaned = cleanup_aux_files(root)
    total_cleaned += cleaned

    # Clean format-based directories (tex/, markdown/, pdf/)
    for fmt_dir in (config.tex_dir, config.markdown_dir, config.pdf_dir):
        for subdir in config.content_subdirs:
            subpath = root / fmt_dir / subdir
            if subpath.exists():
                cleaned = cleanup_all_backups(subpath)
                total_cleaned += cleaned
                cleaned = cleanup_aux_files(subpath)
                total_cleaned += cleaned

    logger.info(f"Cleanup complete: {total_cleaned} files removed")
    return total_cleaned


def find_public_release_folders(parent: Path) -> list[Path]:
    """Find all public-release* folders (including public-release) in parent."""
    folders: list[Path] = []
    for d in parent.iterdir():
        if d.is_dir() and d.name.startswith("public-release"):
            folders.append(d)
    return sorted(folders, key=lambda p: p.name)


def sync_versions(target_root: Path, parent: Path, dry_run: bool = False) -> int:
    """
    Scan all public-release* folders for tex/{paper,supplements,math-papers}/.
    For each file, copy the newest version (by mtime) into target_root, overwriting older.
    Returns number of files synced.
    """
    folders = find_public_release_folders(parent)
    if len(folders) < 2:
        if len(folders) == 1:
            print("Only one public-release folder found; nothing to sync.")
        return 0

    # Build map: (format_dir, subdir, filename) -> (path, mtime) for newest
    newest: dict[tuple[str, str, str], tuple[Path, float]] = {}
    format_dirs = (config.tex_dir, config.markdown_dir, config.pdf_dir)
    for folder in folders:
        for fmt_dir in format_dirs:
            for subdir in CONTENT_SUBDIRS:
                subpath = folder / fmt_dir / subdir
                if not subpath.is_dir():
                    continue
                for f in subpath.iterdir():
                    if f.is_file() and not f.name.startswith("."):
                        key = (fmt_dir, subdir, f.name)
                        mtime = f.stat().st_mtime
                        if key not in newest or mtime > newest[key][1]:
                            newest[key] = (f, mtime)

    actions: list[tuple[Path, Path]] = []
    for (fmt_dir, subdir, filename), (src, _) in newest.items():
        dst = target_root / fmt_dir / subdir / filename
        if not dst.exists() or src.stat().st_mtime > dst.stat().st_mtime:
            if src.resolve() != dst.resolve():
                actions.append((src, dst))

    if not actions:
        print("All versions already up to date.")
        return 0

    print(f"Syncing {len(actions)} file(s) from sibling public-release* folders:")
    for src, dst in actions:
        rel = src.relative_to(parent)
        print(f"  {rel} -> {target_root.name}/{dst.relative_to(target_root)}")

    if dry_run:
        return len(actions)

    for src, dst in actions:
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
    print(f"Synced {len(actions)} file(s).")
    return len(actions)


def find_extensions(root: Path) -> dict[tuple[Path, str], set[str]]:
    """Find all .tex, .pdf, .md files; return {(dir, base_name): {extensions}}."""
    found: dict[tuple[Path, str], set[str]] = {}
    for ext in (".tex", ".pdf", ".md"):
        for p in root.rglob(f"*{ext}"):
            if p.name.startswith("."):
                continue
            key = (p.parent, p.stem)
            if key not in found:
                found[key] = set()
            found[key].add(ext)
    return found


def find_tex_sources(root: Path) -> list[dict]:
    """
    Find all .tex files in tex_dir (root + subdirs) and compute output paths.
    Returns list of dicts: {tex: Path, md: Path, pdf: Path, section: str}
    """
    tex_root = root / config.tex_dir
    sources = []

    # Scan root of tex_dir for .tex files (main paper + supplements)
    for tex_path in sorted(tex_root.glob("*.tex")):
        if tex_path.name.startswith("."):
            continue
        base = tex_path.stem
        sources.append({
            "tex": tex_path,
            "md": root / config.markdown_dir / f"{base}.md",
            "pdf": root / config.pdf_dir / f"{base}.pdf",
            "section": "",
        })

    # Scan subdirectories (math-papers, etc.)
    for section in config.content_subdirs:
        section_dir = tex_root / section
        if not section_dir.is_dir():
            continue
        for tex_path in sorted(section_dir.glob("*.tex")):
            if tex_path.name.startswith("."):
                continue
            base = tex_path.stem
            sources.append({
                "tex": tex_path,
                "md": root / config.markdown_dir / section / f"{base}.md",
                "pdf": root / config.pdf_dir / section / f"{base}.pdf",
                "section": section,
            })
    return sources


def run_cmd(cmd: list[str], cwd: Path | None = None, timeout: int = 120) -> bool:
    """Run command with improved error handling and logging."""
    try:
        logger.debug(f"Running: {' '.join(cmd)}")
        start_time = time.time()

        r = subprocess.run(
            cmd,
            cwd=cwd,
            capture_output=True,
            text=True,
            timeout=timeout
        )

        elapsed = time.time() - start_time
        logger.debug(f"Command completed in {elapsed:.1f}s with return code {r.returncode}")

        if r.returncode != 0:
            error_msg = r.stderr[:500] if r.stderr else '(no stderr)'
            logger.warning(f"Command failed: {error_msg}")
            return False

        return True

    except subprocess.TimeoutExpired:
        logger.error(f"Command timed out after {timeout}s: {' '.join(cmd)}")
        return False
    except FileNotFoundError:
        logger.error(f"Command not found: {cmd[0]}")
        return False
    except Exception as e:
        logger.error(f"Unexpected error running command: {e}")
        return False


def tex_to_md(tex_path: Path, md_path: Path) -> bool:
    """Convert .tex to .md via pandoc."""
    pandoc = find_pandoc()
    if not pandoc:
        logger.error("pandoc not found; install pandoc for tex->md")
        logger.info("Run: install_pandoc_simple.bat  (or winget install JohnMacFarlane.Pandoc)")
        logger.info("See: TROUBLESHOOTING.md")
        logger.info("NOTE: pip install pandoc is a wrapper only; you need the real binary.")
        return False

    # Backup existing file if configured
    if config.backup_existing and md_path.exists():
        backup_path = md_path.with_suffix(f"{md_path.suffix}.backup")
        shutil.copy2(md_path, backup_path)
        logger.debug(f"Backed up {md_path} to {backup_path}")

    cmd = [pandoc] + config.pandoc_options + [str(tex_path), "-o", str(md_path)]
    success = run_cmd(cmd, cwd=tex_path.parent, timeout=config.pandoc_timeout)

    if success:
        logger.info(f"Generated {md_path}")
        # NEW: Auto-clean backups if configured
        if config.auto_clean_backups:
            cleanup_old_backups(md_path.parent, config.backup_retention_days)
    else:
        logger.error(f"Failed to convert {tex_path} to markdown")

    return success


def tex_to_pdf(tex_path: Path, pdf_path: Path) -> bool:
    """Compile .tex to .pdf via pdflatex."""
    if not shutil.which("pdflatex"):
        logger.error("pdflatex not found; install TeX Live or MiKTeX for tex->pdf")
        logger.info("MiKTeX: https://miktex.org/download")
        logger.info("Or update MiKTeX: Start -> MiKTeX Console -> Updates -> Check for updates")
        logger.info("ALTERNATIVE: Use online LaTeX compiler:")
        logger.info("- https://www.overleaf.com/")
        logger.info("- https://papeeria.com/")
        logger.info("- Upload .tex file and download .pdf")
        return False

    # Backup existing file if configured
    if config.backup_existing and pdf_path.exists():
        backup_path = pdf_path.with_suffix(f"{pdf_path.suffix}.backup")
        shutil.copy2(pdf_path, backup_path)
        logger.debug(f"Backed up {pdf_path} to {backup_path}")

    cwd = tex_path.parent
    tex_name = tex_path.name

    # Run pdflatex twice for TOC/references
    for run_num in range(1, 3):
        logger.debug(f"PDFLaTeX run {run_num}/2 for {tex_name}")
        cmd = ["pdflatex"] + config.pdflatex_options + [tex_name]
        if not run_cmd(cmd, cwd=cwd, timeout=config.pdflatex_timeout):
            logger.error(f"PDFLaTeX compilation failed on run {run_num}")
            return False

    # pdflatex writes to same dir; move if needed
    built = cwd / tex_path.with_suffix(".pdf").name
    if built.exists() and built != pdf_path:
        shutil.move(str(built), str(pdf_path))
        logger.debug(f"Moved {built} to {pdf_path}")

    # Clean auxiliary files if configured
    if config.clean_aux_files:
        aux_exts = ['.aux', '.log', '.out', '.fls', '.fdb_latexmk', '.synctex.gz']
        for ext in aux_exts:
            aux_file = cwd / tex_path.with_suffix(ext).name
            if aux_file.exists():
                aux_file.unlink()
                logger.debug(f"Cleaned auxiliary file: {aux_file}")

    success = pdf_path.exists()
    if success:
        logger.info(f"Generated {pdf_path}")
        # NEW: Auto-clean backups if configured
        if config.auto_clean_backups:
            cleanup_old_backups(pdf_path.parent, config.backup_retention_days)
    else:
        logger.error(f"PDF file was not created: {pdf_path}")

    return success


def pdf_to_md(pdf_path: Path, md_path: Path) -> bool:
    """Extract text from .pdf to .md via pymupdf. Falls back to OCR if text is minimal."""
    try:
        import fitz  # pymupdf
    except ImportError:
        logger.error("pymupdf not installed; run: pip install pymupdf")
        return False

    # Backup existing file if configured
    if config.backup_existing and md_path.exists():
        backup_path = md_path.with_suffix(f"{md_path.suffix}.backup")
        shutil.copy2(md_path, backup_path)
        logger.debug(f"Backed up {md_path} to {backup_path}")

    try:
        logger.debug(f"Extracting text from {pdf_path}")
        try:
            doc = fitz.open(pdf_path)
        except Exception as e:
            logger.error(f"Failed to open PDF {pdf_path}: {e}")
            return False

        text = []

        for page_num, page in enumerate(doc):
            logger.debug(f"Processing page {page_num + 1}/{len(doc)}")
            try:
                page_text = page.get_text()
            except Exception as e:
                logger.warning(f"Failed to extract text from page {page_num + 1}: {e}")
                continue

            if not page_text.strip() and config.enable_ocr:
                # Try OCR if no text and OCR is enabled
                try:
                    import pytesseract
                    from io import BytesIO
                    from PIL import Image

                    pix = page.get_pixmap()
                    img = Image.open(BytesIO(pix.tobytes("png")))
                    page_text = pytesseract.image_to_string(img, lang=config.ocr_lang)
                    logger.debug(f"Used OCR for page {page_num + 1}")
                except ImportError:
                    logger.warning("pytesseract not installed; OCR skipped")
                except Exception as e:
                    logger.warning(f"OCR error on page {page_num + 1}: {e}")

            text.append(page_text)

        doc.close()

        base = pdf_path.stem
        md_content = f"# {base}\n\n" + "\n\n".join(text)
        md_path.write_text(md_content, encoding="utf-8")

        logger.info(f"Generated {md_path}")
        # NEW: Auto-clean backups if configured
        if config.auto_clean_backups:
            cleanup_old_backups(md_path.parent, config.backup_retention_days)
        return True

    except Exception as e:
        logger.error(f"Error extracting text from {pdf_path}: {e}")
        return False


def process_sequential(actions: list[tuple[str, Path, Path]], show_progress: bool = True) -> tuple[int, int]:
    """Process conversions sequentially."""
    ok = 0
    fail = 0

    if actions and show_progress and HAS_TQDM:
        iterator = tqdm(actions, desc=config.progress_desc)
    else:
        iterator = actions
    for action, src, dst in iterator:
        if not src.exists():
            logger.warning(f"Skip (src missing): {src}")
            fail += 1
            continue

        logger.info(f"{action}: {src.name}")
        if action == "tex->md":
            success = tex_to_md(src, dst)
        elif action == "tex->pdf":
            success = tex_to_pdf(src, dst)
        elif action == "pdf->md":
            success = pdf_to_md(src, dst)
        else:
            logger.error(f"Unknown action: {action}")
            success = False

        if success:
            ok += 1
        else:
            fail += 1

    return ok, fail


def process_parallel(actions: list[tuple[str, Path, Path]], max_workers: int) -> tuple[int, int]:
    """Process conversions in parallel."""
    logger.info(f"Using parallel processing with {max_workers} workers")

    def process_single(action_src_dst):
        action, src, dst = action_src_dst
        if not src.exists():
            logger.warning(f"Skip (src missing): {src}")
            return False

        logger.debug(f"{action}: {src.name}")
        if action == "tex->md":
            return tex_to_md(src, dst)
        elif action == "tex->pdf":
            return tex_to_pdf(src, dst)
        elif action == "pdf->md":
            return pdf_to_md(src, dst)
        else:
            logger.error(f"Unknown action: {action}")
            return False

    ok = 0
    fail = 0

    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Use tqdm for progress tracking
        if HAS_TQDM:
            results = list(tqdm(
                executor.map(process_single, actions),
                total=len(actions),
                desc=config.progress_desc
            ))
        else:
            results = list(executor.map(process_single, actions))

    for success in results:
        if success:
            ok += 1
        else:
            fail += 1

    return ok, fail


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--root",
        type=Path,
        default=Path(__file__).resolve().parent.parent,
        help="Root directory (default: public-release)",
    )
    ap.add_argument("--dry-run", action="store_true", help="Show what would be done")
    ap.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    ap.add_argument(
        "--sync",
        action="store_true",
        help="First sync newest versions from sibling public-release* folders into --root",
    )
    ap.add_argument(
        "--sync-only",
        action="store_true",
        help="Only sync; skip tex/pdf/md conversion",
    )
    ap.add_argument(
        "--parallel",
        action="store_true",
        help="Use parallel processing for conversions",
    )
    ap.add_argument(
        "--workers",
        type=int,
        default=config.max_workers,
        help=f"Number of parallel workers (default: {config.max_workers})",
    )
    ap.add_argument(
        "--force",
        action="store_true",
        help="Regenerate from .tex/.pdf even when outputs exist (bypasses mtime check; still skips if only .md)",
    )
    ap.add_argument(
        "--cleanup-backups",
        action="store_true",
        help="Clean up ALL backup files immediately (not just old ones)",
    )
    args = ap.parse_args()

    # Adjust logging level based on verbose flag
    if args.verbose:
        logger.setLevel(logging.DEBUG)
        for handler in logger.handlers:
            handler.setLevel(logging.DEBUG)

    root = args.root.resolve()
    if not root.exists():
        logger.error(f"Root not found: {root}")
        sys.exit(1)

    logger.info(f"Starting conversion in {root}")

    # Sync from sibling public-release* folders if requested
    if args.sync or args.sync_only:
        parent = root.parent
        sync_count = sync_versions(root, parent, dry_run=args.dry_run)
        if args.dry_run:
            logger.info("(Dry run: no files changed. Run without --dry-run to apply.)")
            return
        if args.sync_only:
            logger.info(f"Sync complete: {sync_count} files")
            # Clean up all backup files and auxiliary files after sync-only
            perform_cleanup(root, logger)
            return
        logger.info(f"Synced {sync_count} files from sibling folders")

    # Handle backup cleanup if requested
    if args.cleanup_backups:
        perform_cleanup(root, logger)
        return

    sources = find_tex_sources(root)
    actions: list[tuple[str, Path, Path]] = []  # (action, src, dst)

    def mtime(p: Path) -> float:
        return p.stat().st_mtime if p.exists() else 0.0

    for src in sources:
        tex_path = src["tex"]
        md_path = src["md"]
        pdf_path = src["pdf"]

        # Ensure output directories exist
        md_path.parent.mkdir(parents=True, exist_ok=True)
        pdf_path.parent.mkdir(parents=True, exist_ok=True)

        tex_mtime = mtime(tex_path)

        # If .md is newer than .tex, skip (do not overwrite manually edited .md)
        if mtime(md_path) > tex_mtime and not args.force:
            if args.verbose:
                logger.debug(f"Skip (md newer): {md_path.name}")
        else:
            need_md = args.force or not md_path.exists() or mtime(md_path) < tex_mtime
            if need_md:
                # Prefer pdf->md if pdf exists (faster), unless --force
                if pdf_path.exists() and not args.force:
                    actions.append(("pdf->md", pdf_path, md_path))
                else:
                    actions.append(("tex->md", tex_path, md_path))

        need_pdf = args.force or not pdf_path.exists() or mtime(pdf_path) < tex_mtime
        if need_pdf:
            actions.append(("tex->pdf", tex_path, pdf_path))

    if not actions:
        logger.info("Nothing to generate. All outputs up to date.")
        # Clean up all backup files and auxiliary files even when nothing to do
        perform_cleanup(root, logger)
        return

    logger.info(f"Will generate {len(actions)} file(s):")
    for action, src, dst in actions:
        logger.info(f"  {action}: {src.relative_to(root)} -> {dst.name}")

    if args.dry_run:
        # Clean up all backup files and auxiliary files even in dry-run mode
        perform_cleanup(root, logger)
        return

    # Choose processing method
    if args.parallel and config.max_workers > 1:
        ok, fail = process_parallel(actions, args.workers)
    else:
        ok, fail = process_sequential(actions, config.show_progress and HAS_TQDM)

    logger.info(f"Conversion complete: {ok} ok, {fail} failed")
    
    # Clean up all backup files and auxiliary files immediately after processing
    perform_cleanup(root, logger)
    
    if fail:
        sys.exit(1)


if __name__ == "__main__":
    main()


# Basic tests (run with pytest if desired)
def test_find_extensions(tmp_path):
    # Create test files
    (tmp_path / "test.tex").write_text("test")
    (tmp_path / "test.pdf").write_text("test")
    (tmp_path / "other.md").write_text("test")
    
    found = find_extensions(tmp_path)
    assert (tmp_path, "test") in found
    assert ".tex" in found[(tmp_path, "test")]
    assert ".pdf" in found[(tmp_path, "test")]
    assert (tmp_path, "other") in found
    assert found[(tmp_path, "other")] == {".md"}
