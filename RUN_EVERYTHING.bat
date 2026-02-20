@echo off
setlocal enabledelayedexpansion
REM ============================================================
REM  RUN_EVERYTHING.bat — The Resolved Chord: One-Click Runner
REM  Double-click this file. It does everything.
REM ============================================================
REM  Force UTF-8 so Unicode symbols display correctly everywhere
chcp 65001 >nul 2>&1
set PYTHONIOENCODING=utf-8
REM ============================================================
REM
REM  What this does:
REM    1. Checks for Python
REM    2. Validates documentation
REM    3. Runs release gate robustness checks
REM    4. Optional dependency install (use --install)
REM    5. Runs compile_universe.py
REM    6. Runs verification scripts (auto-discovered + hygiene checks)
REM    7. Runs falsification test suite (adversarial pytest)
REM    8. Saves results to output/results.txt
REM
REM  Requirements: Python 3.9+ (https://python.org)
REM  No admin rights needed. No pandoc or LaTeX needed.
REM ============================================================

cd /d "%~dp0"
title The Resolved Chord — Running Everything

REM Create output directory if it doesn't exist
if not exist output mkdir output

REM Generate timestamp for log file (robust method)
set TIMESTAMP=%DATE:~-4%%DATE:~4,2%%DATE:~7,2%_%TIME:~0,2%%TIME:~3,2%%TIME:~6,2%
set TIMESTAMP=%TIMESTAMP: =0%
set TIMESTAMP=%TIMESTAMP:/=%
set TIMESTAMP=%TIMESTAMP:\=%
set LOGFILE=output\session_%TIMESTAMP%.txt

echo.
echo  ================================================================
echo    The Resolved Chord: The Theorem of Everything
echo    Zero free parameters. One geometry: S^5/Z_3
echo  ================================================================
echo.
echo Full session log will be saved to: %LOGFILE%
echo.

REM Function to log and display output
:log
echo %~1
echo %~1 >> "%LOGFILE%"
goto :eof

REM Start logging
call :log "=== THE RESOLVED CHORD SESSION LOG ==="
call :log "Started at: %TIMESTAMP%"
call :log "Working directory: %CD%"
call :log ""

set "INSTALL_DEPS=0"
if /I "%~1"=="--install" set "INSTALL_DEPS=1"

REM --- Step 0: Check Python (adaptive launcher detection) ---
call :log "[1/7] Checking for Python 3.9+..."
set "PYTHON_CMD="

where py >nul 2>&1
if not errorlevel 1 (
    py -3 -c "import sys; raise SystemExit(0 if sys.version_info >= (3,9) else 1)" >nul 2>&1
    if not errorlevel 1 set "PYTHON_CMD=py -3"
)

if not defined PYTHON_CMD (
    where python >nul 2>&1
    if not errorlevel 1 (
        python -c "import sys; raise SystemExit(0 if sys.version_info >= (3,9) else 1)" >nul 2>&1
        if not errorlevel 1 set "PYTHON_CMD=python"
    )
)

if not defined PYTHON_CMD (
    call :log ""
    call :log "ERROR: Python 3.9+ is not installed or not in PATH."
    call :log ""
    call :log "Please install Python 3.9+ from: https://python.org"
    call :log "IMPORTANT: Check 'Add Python to PATH' during installation."
    call :log ""
    pause
    exit /b 1
)
for /f "tokens=*" %%i in ('%PYTHON_CMD% --version 2^>^&1') do call :log "       Found: %%i"
call :log ""

REM --- Step 0.5: Validate Documentation ---
call :log "[2/7] Validating documentation..."
%PYTHON_CMD% tools\validate_docs.py >> "%LOGFILE%" 2>&1
if errorlevel 1 (
    call :log "       Documentation validation failed!"
    call :log "       Please fix the issues above before proceeding."
    call :log ""
    pause
    exit /b 1
)
call :log "       Documentation OK"
call :log ""

REM --- Step 0.6: Release Gate ---
call :log "[3/7] Running release gate checks..."
%PYTHON_CMD% tools\release_gate.py >> "%LOGFILE%" 2>&1
if errorlevel 1 (
    call :log "       Release gate failed!"
    call :log "       Please fix the issues above before proceeding."
    call :log ""
    pause
    exit /b 1
)
call :log "       Release gate OK"
call :log ""

REM --- Step 1: Install dependencies ---
if "%INSTALL_DEPS%"=="1" (
    call :log "[4/7] Installing dependencies..."
    call :log "       Running: %PYTHON_CMD% -m pip install --upgrade pip"
    %PYTHON_CMD% -m pip install --upgrade pip >> "%LOGFILE%" 2>&1
    call :log "       Running: %PYTHON_CMD% -m pip install -r requirements.txt"
    %PYTHON_CMD% -m pip install -r requirements.txt >> "%LOGFILE%" 2>&1
    if errorlevel 1 (
        call :log ""
        call :log "       WARNING: Some packages may have failed to install."
        call :log "       Trying to continue anyway..."
        call :log ""
    ) else (
        call :log "       All dependencies installed."
    )
) else (
    call :log "[4/7] Skipping dependency install (safe default)."
    call :log "       To install dependencies, run: RUN_EVERYTHING.bat --install"
)
call :log ""

REM --- Step 2: Compile the Universe ---
call :log "================================================================="
call :log "[5/7] COMPILING THE UNIVERSE"
call :log "================================================================="
call :log ""
call :log "Running: %PYTHON_CMD% verification\compile_universe.py"
set TEMP_COMPILE=output\temp_compile.txt
%PYTHON_CMD% verification\compile_universe.py > "%TEMP_COMPILE%" 2>&1
set COMPILE_EXIT=%ERRORLEVEL%
type "%TEMP_COMPILE%" >> "%LOGFILE%"
type "%TEMP_COMPILE%"
del "%TEMP_COMPILE%" 2>nul
call :log ""
call :log "compile_universe.py exit code: %COMPILE_EXIT%"
call :log ""

REM --- Step 3: Verification Scripts ---
call :log "================================================================="
call :log "[6/7] RUNNING VERIFICATION SCRIPTS"
call :log "================================================================="
call :log ""
call :log "Running: %PYTHON_CMD% verification\run_verification.py --discover-all --workers 4 --timeout 180 --hygiene --hygiene-all --strict-hygiene --results output\verification_results.txt"
set TEMP_VERIFY=output\temp_verify.txt
%PYTHON_CMD% verification\run_verification.py --discover-all --workers 4 --timeout 180 --hygiene --hygiene-all --strict-hygiene --results output\verification_results.txt > "%TEMP_VERIFY%" 2>&1
set VERIFY_EXIT=%ERRORLEVEL%
type "%TEMP_VERIFY%" >> "%LOGFILE%"
type "%TEMP_VERIFY%"
del "%TEMP_VERIFY%" 2>nul
call :log ""
call :log "run_verification.py exit code: %VERIFY_EXIT%"
call :log ""

REM --- Step 4: Falsification Test Suite ---
call :log "================================================================="
call :log "[7/7] RUNNING FALSIFICATION TEST SUITE"
call :log "================================================================="
call :log ""
call :log "Running: %PYTHON_CMD% tools\falsification\run_falsification.py -v --results output\falsification_results.txt"
set TEMP_FALSIF=output\temp_falsif.txt
%PYTHON_CMD% tools\falsification\run_falsification.py -v --results output\falsification_results.txt > "%TEMP_FALSIF%" 2>&1
set FALSIF_EXIT=%ERRORLEVEL%
type "%TEMP_FALSIF%" >> "%LOGFILE%"
type "%TEMP_FALSIF%"
del "%TEMP_FALSIF%" 2>nul
call :log ""
call :log "run_falsification.py exit code: %FALSIF_EXIT%"
call :log ""

REM --- Write Results ---
call :log "================================================================="
call :log "RESULTS SUMMARY"
call :log "================================================================="
call :log ""

set /a TOTAL_EXIT=%COMPILE_EXIT% + %VERIFY_EXIT% + %FALSIF_EXIT%

if %TOTAL_EXIT% equ 0 (
    call :log "ALL PASSED. The universe compiles correctly."
    call :log "Verification, hygiene, and falsification checks all passed."
) else (
    call :log "Some steps had issues:"
    if %COMPILE_EXIT% neq 0 call :log "  - compile_universe.py: FAILED"
    if %VERIFY_EXIT% neq 0 call :log "  - Verification scripts: SOME FAILED"
    if %FALSIF_EXIT% neq 0 call :log "  - Falsification tests: SOME FAILED"
)

REM Write summary to results.txt
for /f "tokens=*" %%d in ('%PYTHON_CMD% -c "from datetime import datetime; print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))"') do set TIMESTAMP_RESULTS=%%d
(
    echo %TIMESTAMP_RESULTS% ^| RUN_EVERYTHING results:
    echo   compile_universe.py: exit %COMPILE_EXIT%
    echo   run_verification.py: exit %VERIFY_EXIT%
    echo   run_falsification.py: exit %FALSIF_EXIT%
    if %TOTAL_EXIT% equ 0 (echo   OVERALL: ALL PASSED) else (echo   OVERALL: SOME FAILURES)
) > output\results.txt

call :log ""
call :log "Results saved to output\results.txt"
call :log "Full session log saved to: %LOGFILE%"
call :log "================================================================="
call :log "SESSION COMPLETE"
call :log "================================================================="
call :log "Ended at: %DATE% %TIME%"
call :log ""

REM Also display to console
echo.
echo ================================================================
echo  RESULTS SUMMARY
echo ================================================================
echo.
if %TOTAL_EXIT% equ 0 (
    echo  ALL PASSED. The universe compiles correctly.
    echo  Verification, hygiene, and falsification checks all passed.
) else (
    echo  Some steps had issues:
    if %COMPILE_EXIT% neq 0 echo    - compile_universe.py: FAILED
    if %VERIFY_EXIT% neq 0 echo    - Verification scripts: SOME FAILED
    if %FALSIF_EXIT% neq 0 echo    - Falsification tests: SOME FAILED
)
echo.
echo  Results saved to output\results.txt
echo  Full session log saved to: %LOGFILE%
echo.
echo ================================================================
echo  Next steps:
echo    - Read the paper: papers\The_Resolved_Chord_v11.tex
echo    - Explore individual proofs in verification\
echo    - See DISCOVERY_TIMELINE.md for the full story
echo ================================================================
echo.
pause
