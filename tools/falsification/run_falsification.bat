@echo off
setlocal
REM Run falsification test suite. Double-click to run.
chcp 65001 >nul 2>&1
set PYTHONIOENCODING=utf-8
cd /d "%~dp0"
cd ..\..

set "INSTALL_DEPS=0"
if /I "%~1"=="--install" set "INSTALL_DEPS=1"

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
) else (
    echo Skipping dependency install ^(safe default^).
    echo To install dependencies, run: run_falsification.bat --install
    echo.
)

echo === Falsification Test Suite ===
echo.

python tools\falsification\run_falsification.py --results falsification_results.txt
set EXITCODE=%ERRORLEVEL%

echo.
if exist falsification_results.txt (
    echo Results saved to falsification_results.txt:
    type falsification_results.txt
)
echo.
pause
