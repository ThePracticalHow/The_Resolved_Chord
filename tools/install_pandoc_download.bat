@echo off
setlocal
REM Safe-mode installer helper: no scripted binary download/extraction.
REM This avoids common antivirus false positives from batch-driven web downloads.

echo === Pandoc Manual Install (Safe Mode) ===
echo.
echo This script intentionally does NOT download binaries automatically.
echo Please install Pandoc from the official release page:
echo   https://github.com/jgm/pandoc/releases/latest
echo.
echo Recommended:
echo   1. Download the MSI installer and run it.
echo   2. Restart your terminal.
echo   3. Run: test_pandoc_install.bat
echo.
echo If you use the ZIP package instead of MSI:
echo   - Extract to C:\pandoc\
echo   - Add C:\pandoc\ to PATH
echo.
pause
