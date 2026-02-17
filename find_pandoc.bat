@echo off
REM Search for pandoc executable anywhere on the system
REM Paths match pandoc_utils.py (Python find_pandoc)
echo === Searching for pandoc ===
echo.

echo Checking PATH...
where pandoc 2>nul
if %ERRORLEVEL% EQU 0 goto :found

echo.
echo Searching common locations...
for %%D in (C:\pandoc C:\Program Files\Pandoc "C:\Program Files (x86)\Pandoc" %LOCALAPPDATA%\Pandoc %USERPROFILE%\AppData\Local\Pandoc) do (
    if exist "%%D\pandoc.exe" (
        echo FOUND: %%D\pandoc.exe
        goto :found
    )
)

echo.
echo Searching Program Files...
dir /s /b "C:\Program Files\pandoc.exe" 2>nul
dir /s /b "C:\Program Files (x86)\pandoc.exe" 2>nul

echo.
echo If pandoc was found above but convert still fails, add its folder to PATH:
echo   System Properties -^> Environment Variables -^> Path -^> Edit -^> New
echo.
goto :end

:found
echo.
echo Pandoc found. If convert_tex_pdf_md.py still fails, restart your terminal.

:end
pause
