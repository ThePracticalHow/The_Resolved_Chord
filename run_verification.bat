@echo off
REM Run verification scripts. Double-click to run.
cd /d "%~dp0"

echo === Verification Scripts ===
echo.

python run_verification.py --results verification_results.txt
set EXITCODE=%ERRORLEVEL%

echo.
if exist verification_results.txt (
    echo Results saved to verification_results.txt:
    type verification_results.txt
)
echo.
pause
