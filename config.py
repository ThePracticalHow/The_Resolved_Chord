#!/usr/bin/env python3
"""
Configuration for convert_tex_pdf_md.py
Edit this file to customize conversion behavior.
"""

import os
from pathlib import Path
from typing import Dict, List, Optional

# Conversion settings
class Config:
    """Configuration for the conversion system."""

    def __init__(self):
        # Pandoc settings
        self.pandoc_options = [
            "-f", "latex",
            "-t", "markdown",
            "--wrap=none",  # Don't wrap lines
            "--markdown-headings=atx",  # Use # headers
        ]

        # PDFLaTeX settings
        self.pdflatex_options = [
            "-interaction=nonstopmode",
            "-halt-on-error",
        ]

        # Parallel processing
        self.max_workers = min(4, os.cpu_count() or 1)  # Use up to 4 cores

        # Logging
        self.log_level = "INFO"  # DEBUG, INFO, WARNING, ERROR
        self.log_file = None  # Set to Path for file logging

        # Output preferences
        self.backup_existing = True  # Backup existing files before overwriting
        self.clean_aux_files = True  # Remove .aux, .log, etc. after PDF generation
        self.auto_clean_backups = True  # Auto-clean old backups after conversion
        self.backup_retention_days = 7  # Days to keep backups (0 = forever)

        # Content directories to process
        self.content_subdirs = ("paper", "supplements", "math-papers")

        # File extensions to track
        self.content_exts = (".tex", ".pdf", ".md")

        # Timeout settings (seconds)
        self.pandoc_timeout = 300
        self.pdflatex_timeout = 120
        self.pdf_extract_timeout = 60

        # OCR settings
        self.enable_ocr = True
        self.ocr_lang = "eng"

        # Progress display
        self.show_progress = True
        self.progress_desc = "Converting"

# Global config instance
config = Config()

def load_config_from_file(config_path: Optional[Path] = None) -> Config:
    """Load configuration from a Python file."""
    if config_path is None:
        # Look for config.py in the same directory as this file
        config_path = Path(__file__).parent / "config.py"

    if config_path.exists():
        try:
            # Execute the config file in a controlled environment
            import builtins
            config_globals = {"__builtins__": builtins, "os": os, "Path": Path}
            config_locals = {}
            exec(config_path.read_text(), config_globals, config_locals)

            # Update our config with any defined variables
            for key, value in config_locals.items():
                if hasattr(config, key):
                    setattr(config, key, value)

        except Exception as e:
            print(f"Warning: Could not load config from {config_path}: {e}")

    return config

def save_default_config(config_path: Path) -> None:
    """Save the default configuration to a file."""
    config_template = '''# Configuration for convert_tex_pdf_md.py
# Edit these values to customize conversion behavior

# Pandoc options for tex->md conversion
pandoc_options = [
    "-f", "latex",
    "-t", "markdown",
    "--wrap=none",
    "--markdown-headings=atx",
]

# PDFLaTeX options for tex->pdf compilation
pdflatex_options = [
    "-interaction=nonstopmode",
    "-halt-on-error",
]

# Parallel processing (0 = sequential)
max_workers = 4

# Logging level: DEBUG, INFO, WARNING, ERROR
log_level = "INFO"

# Set to a filename to enable file logging
log_file = None

# Backup existing files before overwriting
backup_existing = True

# Clean auxiliary files (.aux, .log, etc.) after PDF generation
clean_aux_files = True

# Automatically clean up backup files after successful conversion
auto_clean_backups = True

# Backup retention policy (days to keep backups, 0 = keep forever)
backup_retention_days = 7

# Content directories to process
content_subdirs = ("paper", "supplements", "math-papers")

# File extensions to track
content_exts = (".tex", ".pdf", ".md")

# Timeout settings (seconds)
pandoc_timeout = 300
pdflatex_timeout = 120
pdf_extract_timeout = 60

# OCR settings
enable_ocr = True
ocr_lang = "eng"

# Progress display
show_progress = True
progress_desc = "Converting"
'''

    config_path.write_text(config_template)
    print(f"Default config saved to {config_path}")

if __name__ == "__main__":
    # When run directly, save default config
    save_default_config(Path("config.py"))