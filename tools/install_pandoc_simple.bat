@echo off
setlocal
REM Automated pandoc installer for Windows
echo === Pandoc Installer ===
echo.

set /p CONFIRM=Proceed with package-manager install attempts (winget/choco)? [y/N]: 
if /I not "%CONFIRM%"=="y" (
    echo Aborted by user.
    goto :end
)

REM Try winget first (Windows 10/11)
where winget >nul 2>&1
if %ERRORLEVEL% EQU 0 (
    echo Trying winget...
    winget install --source winget --exact --id JohnMacFarlane.Pandoc --accept-package-agreements
    if %ERRORLEVEL% EQU 0 (
        echo.
        echo [OK] Pandoc installed via winget.
        echo RESTART YOUR TERMINAL, then run: test_pandoc_install.bat
        goto :end
    )
)

REM Try chocolatey
where choco >nul 2>&1
if %ERRORLEVEL% EQU 0 (
    echo Trying Chocolatey...
    choco install pandoc -y
    if %ERRORLEVEL% EQU 0 (
        echo.
        echo [OK] Pandoc installed via Chocolatey.
        echo RESTART YOUR TERMINAL, then run: test_pandoc_install.bat
        goto :end
    )
)

REM Manual instructions
echo.
echo Automatic install failed. Manual install:
echo.
echo 1. Download: https://github.com/jgm/pandoc/releases/latest
echo 2. Get: pandoc-*-windows-x86_64.zip (or .msi)
echo 3. MSI: Run installer (adds to PATH automatically)
echo    ZIP: Extract to C:\pandoc\ then add C:\pandoc\ to PATH
echo.
echo 4. Add to PATH: System Properties -^> Environment Variables -^> Path -^> Edit -^> New -^> C:\pandoc
echo 5. RESTART your terminal
echo.
echo See TROUBLESHOOTING.md for full guide.

:end
echo.
echo Press any key to continue...
pause >nul
