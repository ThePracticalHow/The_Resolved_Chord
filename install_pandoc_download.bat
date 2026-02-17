@echo off
REM Download and install pandoc when winget/Chocolatey fail. Run install_pandoc_simple.bat first.
set "PANDOC_DIR=C:\pandoc"
set "PANDOC_ZIP=%TEMP%\pandoc.zip"
REM Use latest release URL - will automatically redirect to newest version
set "PANDOC_URL=https://github.com/jgm/pandoc/releases/latest/download/pandoc-windows-x86_64.zip"

echo Installing pandoc via direct download...
powershell -Command "& {try { Invoke-WebRequest -Uri '%PANDOC_URL%' -OutFile '%PANDOC_ZIP%' -UseBasicParsing } catch { exit 1 }}" 2>nul
if not exist "%PANDOC_ZIP%" (
    echo Download failed. Manual: https://github.com/jgm/pandoc/releases/latest
    pause
    exit /b 1
)
if exist "%PANDOC_DIR%" rmdir /s /q "%PANDOC_DIR%"
powershell -Command "& { Expand-Archive -Path '%PANDOC_ZIP%' -DestinationPath '%PANDOC_DIR%' -Force }" 2>nul
if not exist "%PANDOC_DIR%\pandoc.exe" (
    echo Extraction failed.
    pause
    exit /b 1
)
echo Add to PATH: %PANDOC_DIR%
echo Restart terminal, then run: test_pandoc_install.bat
del "%PANDOC_ZIP%" 2>nul
pause
