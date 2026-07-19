#Requires -Version 5.1
<#
.SYNOPSIS
  Install / build the WR Co-Simulation Studio on native Windows (PowerShell).

.DESCRIPTION
  The Windows counterpart of setup.sh. Does what it safely can automatically:
    * creates a Python virtualenv (.venv) and installs numpy + matplotlib into it
    * builds main.exe via make
    * runs `npm install` for the optional schematic layout (if node is present)
  For system tools it cannot install without a package manager (a C++ compiler,
  make, Xyce, pdflatex, poppler, node) it checks presence and prints an install
  hint for the detected package manager (winget / scoop / choco).

  Note: the build needs a GCC/Clang toolchain that works with the Makefile
  (MSYS2/MinGW or scoop's gcc). Plain MSVC (cl.exe) is not used. If you have WSL2,
  running ./setup.sh inside a Linux distro is the smoother path.

.PARAMETER NoVenv
  Install numpy/matplotlib into the current Python instead of a .venv.

.PARAMETER NoOptional
  Skip the optional schematic/PDF tools (node, pdflatex, poppler).

.EXAMPLE
  ./setup.ps1
  ./setup.ps1 -NoVenv -NoOptional
#>
param(
    [switch]$NoVenv,
    [switch]$NoOptional,
    [switch]$Help
)

$ErrorActionPreference = 'Stop'
Set-Location $PSScriptRoot

if ($Help) {
    Get-Help $PSCommandPath -Detailed
    return
}

# --------------------------------------------------------------------------- #
# pretty output
# --------------------------------------------------------------------------- #
function Hdr ($m) { Write-Host "`n==> $m" -ForegroundColor White }
function Ok   ($m) { Write-Host "  [ok]   $m" -ForegroundColor Green }
function Warn ($m) { Write-Host "  [warn] $m" -ForegroundColor Yellow }
function Err  ($m) { Write-Host "  [MISS] $m" -ForegroundColor Red }
function Note ($m) { Write-Host "  $m" -ForegroundColor DarkGray }
function Have ($c) { [bool](Get-Command $c -ErrorAction SilentlyContinue) }

$MissingRequired = @()

# --------------------------------------------------------------------------- #
# package manager detection + install hints
# --------------------------------------------------------------------------- #
$PKG = if (Have winget) { 'winget' } elseif (Have scoop) { 'scoop' } elseif (Have choco) { 'choco' } else { '' }

# Hint <winget-id> <scoop-pkg> <choco-pkg>: print the install command for the detected manager.
function Hint ($winget, $scoop, $choco) {
    switch ($PKG) {
        'winget' { if ($winget) { Note "install: winget install --id $winget" } else { Note 'no winget package; try scoop/choco' } }
        'scoop'  { if ($scoop)  { Note "install: scoop install $scoop" }         else { Note 'no scoop package; try winget/choco' } }
        'choco'  { if ($choco)  { Note "install: choco install $choco" }         else { Note 'no choco package; try winget/scoop' } }
        default  { Note 'install a package manager first: winget (built-in on Win 10/11), scoop, or choco' }
    }
}

Hdr 'Environment'
Ok "OS: Windows$(if ($PKG) { "  (package manager: $PKG)" })"

# --------------------------------------------------------------------------- #
# core toolchain
# --------------------------------------------------------------------------- #
Hdr 'Core toolchain (required to build + run)'

# C++ compiler (MinGW/Clang; MSVC is not used by the Makefile)
$CXX = @('g++', 'clang++') | Where-Object { Have $_ } | Select-Object -First 1
if ($CXX) {
    Ok "C++ compiler: $CXX ($((& $CXX --version 2>$null | Select-Object -First 1)))"
} else {
    Err 'no MinGW/Clang C++ compiler (need g++ or clang++, C++17)'
    Note 'easiest: MSYS2 -> `pacman -S mingw-w64-ucrt-x86_64-gcc make`, then add its /ucrt64/bin to PATH'
    Hint 'MSYS2.MSYS2' 'gcc' 'mingw'
    $MissingRequired += 'C++ compiler'
}

# make
if (Have make) {
    Ok "make: $((Get-Command make).Source)"
} else {
    Err 'make not found'
    Note 'comes with MSYS2 (`pacman -S make`); or: scoop install make'
    Hint '' 'make' 'make'
    $MissingRequired += 'make'
}

# python (python.exe or the py launcher)
$PY = if (Have python) { 'python' } elseif (Have py) { 'py' } else { '' }
if ($PY) {
    Ok "python: $((& $PY --version 2>&1))"
} else {
    Err 'python not found'
    Hint 'Python.Python.3.12' 'python' 'python'
    $MissingRequired += 'python'
}

# Xyce (external circuit solver; cannot auto-install)
if (Have Xyce) {
    Ok "Xyce: $((Get-Command Xyce).Source)"
} else {
    Err 'Xyce not on PATH -- the solver cannot run without it'
    Note 'download the Windows build from https://xyce.sandia.gov (no winget/scoop/choco package)'
    Note "then add Xyce's bin directory to PATH"
    $MissingRequired += 'Xyce'
}

# --------------------------------------------------------------------------- #
# Python dependencies (numpy + matplotlib)
# --------------------------------------------------------------------------- #
Hdr 'Python dependencies (numpy, matplotlib)'
$PyBin = $PY
$useVenv = -not $NoVenv
if ($PY) {
    if ($useVenv) {
        if (-not (Test-Path .venv)) {
            try { & $PY -m venv .venv; Ok 'created virtualenv .venv' }
            catch { Err 'python -m venv failed'; Warn 'falling back to the current Python'; $useVenv = $false }
        } else {
            Ok 'reusing existing .venv'
        }
    }
    if ($useVenv) {
        $PyBin = Join-Path $PSScriptRoot '.venv\Scripts\python.exe'
        & $PyBin -m pip install --quiet --upgrade pip 2>$null | Out-Null
        & $PyBin -m pip install --quiet numpy matplotlib
        if ($LASTEXITCODE -eq 0) { Ok 'installed numpy + matplotlib into .venv' }
        else { Err 'pip install failed inside .venv'; $MissingRequired += 'numpy/matplotlib' }
    } else {
        & $PY -c 'import numpy, matplotlib' 2>$null
        if ($LASTEXITCODE -eq 0) { Ok 'numpy + matplotlib already importable' }
        else {
            Warn 'numpy/matplotlib not importable in the current Python'
            Note "install: $PY -m pip install numpy matplotlib   (or re-run without -NoVenv)"
            $MissingRequired += 'numpy/matplotlib'
        }
    }
} else {
    Warn 'skipping Python deps (no python)'
}

# --------------------------------------------------------------------------- #
# build main.exe
# --------------------------------------------------------------------------- #
Hdr 'Build main.exe'
if ($CXX -and (Have make)) {
    $log = Join-Path $env:TEMP 'wrcosim_build.log'
    # Makefile defaults CXX to clang++; override to the compiler we actually found.
    & make "CXX=$CXX" *> $log
    if ($LASTEXITCODE -eq 0 -and (Test-Path .\main.exe)) {
        Ok 'built main.exe'
    } else {
        Err "build failed -- see $log"
        Get-Content $log -Tail 15 | ForEach-Object { Write-Host "    $_" }
        $MissingRequired += 'main.exe build'
    }
} else {
    Warn 'skipping build (need a MinGW/Clang compiler + make)'
}

# --------------------------------------------------------------------------- #
# optional: schematic (node + elkjs) and PDF/PNG (pdflatex + poppler)
# --------------------------------------------------------------------------- #
if (-not $NoOptional) {
    Hdr 'Optional: circuit schematic + PDF export'

    if (Have npm) {
        & npm install --silent 2>$null | Out-Null
        if ($LASTEXITCODE -eq 0) { Ok 'npm install (elkjs) done -- schematic layout available' }
        else { Warn 'npm install failed -- schematic layout unavailable (non-fatal)' }
    } else {
        Warn 'node/npm not found -- inline schematic unavailable (non-fatal)'
        Hint 'OpenJS.NodeJS' 'nodejs' 'nodejs'
    }

    if (Have pdflatex) {
        Ok "pdflatex: $((Get-Command pdflatex).Source)  (.tex/.pdf export + schematic PDF)"
    } else {
        Warn 'pdflatex not found -- PDF export + schematic render unavailable (non-fatal)'
        Hint 'MiKTeX.MiKTeX' 'latex' 'miktex'
    }

    # poppler (pdftoppm) rasterizes the schematic PDF -> inline PNG. There is no `sips` on Windows,
    # so without pdftoppm the inline schematic PNG won't render (the .tex/.pdf export still works).
    if (Have pdftoppm) {
        Ok "pdftoppm (poppler): $((Get-Command pdftoppm).Source)  (inline schematic PNG)"
    } else {
        Warn 'pdftoppm (poppler) not found -- inline schematic PNG unavailable (non-fatal)'
        Hint '' 'poppler' 'poppler'
    }
} else {
    Note '(skipping optional schematic/PDF tools: -NoOptional)'
}

# --------------------------------------------------------------------------- #
# summary
# --------------------------------------------------------------------------- #
Hdr 'Summary'
if ($MissingRequired.Count -eq 0) {
    Ok 'all required components are in place'
} else {
    Warn "missing required components: $($MissingRequired -join ', ')"
    Note 'resolve the [MISS] items above, then re-run ./setup.ps1'
}

Write-Host "`nRun the studio:" -ForegroundColor White
if ($useVenv) {
    Write-Host '  .\.venv\Scripts\python sim_ui.py'
} else {
    Write-Host "  $PY sim_ui.py"
}
Write-Host '  then open' -ForegroundColor DarkGray -NoNewline; Write-Host ' http://127.0.0.1:8000'
Write-Host "`nHeadless sweep:" -ForegroundColor White
$pyShow = if ($useVenv) { '.\.venv\Scripts\python' } else { $PY }
Write-Host "  $pyShow sim_ui.py sweep --param L_FEM --min 1e-7 --max 1e-5 --steps 8 --scale log"
Write-Host ''
