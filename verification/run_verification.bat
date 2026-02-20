@echo off
setlocal
REM Run verification scripts. Double-click to run.
chcp 65001 >nul 2>&1
set PYTHONIOENCODING=utf-8
cd /d "%~dp0"

set "INSTALL_DEPS=0"
if /I "%~1"=="--install" (
    set "INSTALL_DEPS=1"
)

REM --- Check Python ---
python --version >nul 2>&1
if errorlevel 1 (
    echo  ERROR: Python not found. Install from https://python.org
    pause
    exit /b 1
)

REM --- Optional dependency install ---
if "%INSTALL_DEPS%"=="1" (
    echo Installing dependencies...
    python -m pip install -r requirements.txt
    echo.
)

echo === Verification Scripts ===
echo.

python run_verification.py --discover-all --changed-only --workers 4 --hygiene --hygiene-all --results verification_results.txt
set EXITCODE=%ERRORLEVEL%

echo.
if exist verification_results.txt (
    echo Results saved to verification_results.txt:
    type verification_results.txt
)
echo.
pause
