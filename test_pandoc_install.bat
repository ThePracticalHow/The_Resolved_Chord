@echo off
REM Test if pandoc is installed and accessible
echo === Pandoc Installation Test ===
echo.

where pandoc >nul 2>&1
if %ERRORLEVEL% EQU 0 (
    echo [OK] Pandoc is in PATH
    echo.
    pandoc --version
    echo.
    echo Location:
    where pandoc
) else (
    echo [FAIL] Pandoc is NOT in PATH
    echo.
    echo You need the actual pandoc executable, NOT the Python package.
    echo   pip install pandoc  = WRONG (just a wrapper)
    echo   Download from GitHub = CORRECT
    echo.
    echo Run: install_pandoc_simple.bat
    echo Or see: TROUBLESHOOTING.md
)

echo.
echo Press any key to continue...
pause >nul
