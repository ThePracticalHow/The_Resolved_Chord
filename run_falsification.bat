@echo off
REM Run falsification test suite. Double-click to run.
cd /d "%~dp0"

echo === Falsification Test Suite ===
echo.

python run_falsification.py --results falsification_results.txt
set EXITCODE=%ERRORLEVEL%

echo.
if exist falsification_results.txt (
    echo Results saved to falsification_results.txt:
    type falsification_results.txt
)
echo.
pause
