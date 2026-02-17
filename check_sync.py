#!/usr/bin/env python3
"""Compare .tex, .pdf, .md files in public-release to report sync status."""
from pathlib import Path
from datetime import datetime

# Files generated in quick succession can differ only by filesystem timestamp
# precision. Treat sub-second drift as synced to avoid false positives.
MTIME_EPSILON_SECONDS = 1.0

ROOT = Path(__file__).resolve().parent

def find_extensions(root: Path) -> dict:
    """Find all .tex, .pdf, .md files; return {(dir, base): {ext: mtime}}."""
    found: dict = {}
    for ext in (".tex", ".pdf", ".md"):
        for p in root.rglob(f"*{ext}"):
            if p.name.startswith(".") or "backup" in p.name.lower():
                continue
            key = (p.parent, p.stem)
            if key not in found:
                found[key] = {}
            try:
                found[key][ext] = p.stat().st_mtime
            except OSError:
                found[key][ext] = 0
    return found

def main():
    found = find_extensions(ROOT)
    # Filter to content subdirs only
    content_dirs = {"paper", "supplements", "math-papers"}

    out_of_sync = []
    md_only = []
    missing_outputs = []

    for (dirpath, base), data in sorted(found.items()):
        try:
            rel_dir = dirpath.relative_to(ROOT)
        except ValueError:
            continue
        if rel_dir.parts and rel_dir.parts[0] not in content_dirs:
            continue
        if "verification" in str(rel_dir) or "figures" in str(rel_dir):
            continue

        tex_m = data.get(".tex", 0) or 0
        pdf_m = data.get(".pdf", 0) or 0
        md_m = data.get(".md", 0) or 0

        rel = rel_dir / base if rel_dir != Path(".") else base

        if data.keys() == {".md"}:
            md_only.append(str(rel) + ".md")
            continue

        if ".tex" in data:
            src_mtime = tex_m
            if tex_m - md_m > MTIME_EPSILON_SECONDS and md_m > 0:
                out_of_sync.append((str(rel), "tex", "md", tex_m, md_m))
            if tex_m - pdf_m > MTIME_EPSILON_SECONDS and pdf_m > 0:
                out_of_sync.append((str(rel), "tex", "pdf", tex_m, pdf_m))
            if ".md" not in data:
                missing_outputs.append((str(rel), "md (from tex)"))
            if ".pdf" not in data:
                missing_outputs.append((str(rel), "pdf (from tex)"))
        elif ".pdf" in data:
            src_mtime = pdf_m
            if pdf_m - md_m > MTIME_EPSILON_SECONDS and md_m > 0:
                out_of_sync.append((str(rel), "pdf", "md", pdf_m, md_m))
            if ".md" not in data:
                missing_outputs.append((str(rel), "md (from pdf)"))

    print("=" * 60)
    print("SYNC STATUS: public-release")
    print("=" * 60)

    if out_of_sync:
        print("\n[!] OUT OF SYNC (source newer than output):")
        for rel, src, dst, src_t, dst_t in out_of_sync:
            src_dt = datetime.fromtimestamp(src_t).strftime("%Y-%m-%d %H:%M:%S")
            dst_dt = datetime.fromtimestamp(dst_t).strftime("%Y-%m-%d %H:%M:%S")
            print(f"  {rel}: .{src} ({src_dt}) newer than .{dst} ({dst_dt})")
    else:
        print("\n[OK] No out-of-sync files (all outputs are current or missing)")

    if missing_outputs:
        print("\n[!] MISSING OUTPUTS:")
        for rel, what in missing_outputs:
            print(f"  {rel}: missing {what}")
    else:
        print("\n[OK] No missing outputs")

    if md_only:
        print(f"\n[i] MD-only (skipped by convert): {len(md_only)} file(s)")
        for x in md_only[:5]:
            print(f"  {x}")
        if len(md_only) > 5:
            print(f"  ... and {len(md_only) - 5} more")

    # Summary table for files with multiple formats
    print("\n" + "-" * 60)
    print("MTIME SUMMARY (tex/pdf/md):")
    print("-" * 60)
    for (dirpath, base), data in sorted(found.items()):
        try:
            rel_dir = dirpath.relative_to(ROOT)
        except ValueError:
            continue
        if rel_dir.parts and rel_dir.parts[0] not in content_dirs:
            continue
        if "verification" in str(rel_dir) or "figures" in str(rel_dir):
            continue
        rel = str(rel_dir / base) if str(rel_dir) != "." else base
        parts = []
        for ext in (".tex", ".pdf", ".md"):
            if ext in data and data[ext]:
                t = datetime.fromtimestamp(data[ext]).strftime("%Y-%m-%d %H:%M:%S")
                parts.append(f"{ext[1:]}:{t}")
            elif ext in data:
                parts.append(f"{ext[1:]}:—")
        if len(parts) >= 2:
            print(f"  {rel}: {' | '.join(parts)}")

if __name__ == "__main__":
    main()
