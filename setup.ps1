#Requires -Version 5.1
<#
.SYNOPSIS
  Install / build the WR Co-Simulation Studio on native Windows (PowerShell).

.DESCRIPTION
  The Windows counterpart of setup.sh. Does what it safely can automatically:
    * creates a Python virtualenv (.venv) and installs numpy + matplotlib into it
    * builds main.exe via make
    * runs `npm install` for the optional schematic layout (if node is present)
  For system tools (a C++ compiler, make, Xyce, pdflatex, poppler, node) it checks
  presence and prints an install hint for the detected package manager
  (winget / scoop / choco). With -Install it runs those install commands itself.

  Note: the build needs a GCC/Clang toolchain that works with the Makefile
  (MSYS2/MinGW or scoop's gcc). Plain MSVC (cl.exe) is not used. If you have WSL2,
  running ./setup.sh inside a Linux distro is the smoother path.

.PARAMETER Install
  Answer yes to every install prompt: install the missing system tools with the
  detected package manager without asking. Notes:
    * choco needs an elevated shell; scoop must NOT be run elevated.
    * on winget the toolchain comes from MSYS2: the package is installed and then
      bootstrapped with `pacman -S mingw-w64-ucrt-x86_64-gcc make` (a few minutes).
    * Xyce has no package on any manager and is never auto-installed.
    * PATH additions are made for this session only; the permanent line to add is
      printed so you can put it in your user PATH yourself.

.PARAMETER InstallXyce
  Answer yes to the Xyce prompt: download Sandia's prebuilt installer and run it.
  Kept out of -Install because it is a ~41 MB download plus an elevated installer
  run. The archive's SHA-256 is pinned in this script and verified before
  anything is executed; a mismatch aborts.

.PARAMETER NoInstall
  Never install and never prompt; only report what is missing. Redirected input
  (CI) behaves the same way on its own.

.PARAMETER NoVenv
  Install numpy/matplotlib into the current Python instead of a .venv.

.PARAMETER NoOptional
  Skip the optional schematic/PDF tools (node, pdflatex, poppler).

.EXAMPLE
  ./setup.ps1
  ./setup.ps1 -Install
  ./setup.ps1 -Install -InstallXyce
  ./setup.ps1 -NoVenv -NoOptional
#>
param(
    [switch]$Install,
    [switch]$InstallXyce,
    [switch]$NoInstall,
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
$DidInstall = $false   # set once anything was actually installed (drives the final hint)

# Confirm-Install <question> [<auto-yes>]
# Returns $true for yes. The auto-yes flag defaults to -Install, so -Install answers
# every prompt yes; -NoInstall answers no. Redirected input (CI) answers no rather
# than blocking forever on Read-Host.
function Confirm-Install {
    param([Parameter(Mandatory)][string]$Question, [bool]$Auto = $Install.IsPresent)
    if ($Auto) { return $true }
    if ($NoInstall) { return $false }
    if ([Console]::IsInputRedirected) { return $false }
    Write-Host "  $Question [y/N] " -ForegroundColor Yellow -NoNewline
    $reply = Read-Host
    return ($reply -match '^\s*y(es)?\s*$')
}

# Run a native command and return its exit code.
# $ErrorActionPreference is reset to Continue for the call: in Windows PowerShell 5.1 a
# native command that writes to stderr *while stderr is redirected* raises a terminating
# NativeCommandError under 'Stop', which would abort setup on a mere compiler warning.
function Invoke-Native {
    param(
        [Parameter(Mandatory)][string]$Exe,
        [string[]]$Arguments = @(),
        [string]$LogFile,
        [switch]$Show
    )
    $prev = $ErrorActionPreference
    $ErrorActionPreference = 'Continue'
    $global:LASTEXITCODE = 0
    try {
        if     ($LogFile) { & $Exe @Arguments *> $LogFile }
        elseif ($Show)    { & $Exe @Arguments 2>&1 | Out-Host }
        else              { & $Exe @Arguments *> $null }
        return $LASTEXITCODE
    } finally { $ErrorActionPreference = $prev }
}

# A freshly installed tool lands in the machine/user PATH, not in this process's copy.
# Append whatever is new, keeping session-only entries the user already had.
function Update-PathFromRegistry {
    $known = @($env:PATH -split ';' | Where-Object { $_ })
    foreach ($scope in 'Machine', 'User') {
        foreach ($p in ([Environment]::GetEnvironmentVariable('PATH', $scope) -split ';')) {
            if ($p -and ($known -notcontains $p)) { $env:PATH += ";$p"; $known += $p }
        }
    }
}

# --------------------------------------------------------------------------- #
# package manager detection + install hints
# --------------------------------------------------------------------------- #
$PKG = if (Have winget) { 'winget' } elseif (Have scoop) { 'scoop' } elseif (Have choco) { 'choco' } else { '' }
$IsAdmin = ([Security.Principal.WindowsPrincipal] [Security.Principal.WindowsIdentity]::GetCurrent()
           ).IsInRole([Security.Principal.WindowsBuiltInRole]::Administrator)

# Hint <winget-id> <scoop-pkg> <choco-pkg>: print the install command for the detected manager.
function Hint ($winget, $scoop, $choco) {
    switch ($PKG) {
        'winget' { if ($winget) { Note "install: winget install --id $winget" } else { Note 'no winget package; try scoop/choco' } }
        'scoop'  { if ($scoop)  { Note "install: scoop install $scoop" }         else { Note 'no scoop package; try winget/choco' } }
        'choco'  { if ($choco)  { Note "install: choco install $choco" }         else { Note 'no choco package; try winget/scoop' } }
        default  { Note 'install a package manager first: winget (built-in on Win 10/11), scoop, or choco' }
    }
}

# Install-Dep <label> <winget-id> <scoop-pkg> <choco-pkg> [-Probe cmd]: run the detected
# manager's install command. Returns $true when the tool is usable afterwards.
function Install-Dep {
    param([string]$Label, [string]$Winget, [string]$Scoop, [string]$Choco, [string]$Probe)

    if (-not $PKG) { Warn "cannot install $Label -- no winget/scoop/choco on this machine"; return $false }
    $pkg = switch ($PKG) { 'winget' { $Winget } 'scoop' { $Scoop } 'choco' { $Choco } }
    if (-not $pkg) { Warn "cannot install $Label with $PKG -- no package for it"; Hint $Winget $Scoop $Choco; return $false }

    Note "installing $Label via $PKG ($pkg) ..."
    $code = switch ($PKG) {
        'winget' { Invoke-Native winget @('install', '--id', $pkg, '-e', '--source', 'winget',
                                          '--accept-package-agreements', '--accept-source-agreements') -Show }
        'scoop'  { Invoke-Native scoop  @('install', $pkg) -Show }
        'choco'  { Invoke-Native choco  @('install', $pkg, '-y') -Show }
    }
    Update-PathFromRegistry

    $script:DidInstall = $true
    $usable = if ($Probe) { Have $Probe } else { $code -eq 0 }
    if ($usable) { Ok "$Label installed"; return $true }
    Warn "$Label not usable after install (exit $code) -- may need a new shell, admin rights, or a manual install"
    return $false
}

# winget's MSYS2 package ships only the base system: gcc/make come from pacman afterwards,
# and the installer puts neither /ucrt64/bin nor /usr/bin on PATH.
function Install-Msys2Toolchain {
    $root = @("$env:SystemDrive\msys64", "$env:ProgramData\msys64", "$env:LOCALAPPDATA\msys64") |
            Where-Object { Test-Path (Join-Path $_ 'usr\bin\bash.exe') } | Select-Object -First 1
    if (-not $root) { Warn 'MSYS2 not found after install -- cannot bootstrap gcc/make'; return $false }

    Note 'bootstrapping the MSYS2 toolchain (pacman -S gcc make) -- this takes a few minutes'
    $bash = Join-Path $root 'usr\bin\bash.exe'
    Invoke-Native $bash @('-lc', 'pacman -Syu --noconfirm --needed mingw-w64-ucrt-x86_64-gcc make') -Show | Out-Null

    $bins = @((Join-Path $root 'ucrt64\bin'), (Join-Path $root 'usr\bin')) | Where-Object { Test-Path $_ }
    foreach ($d in $bins) { if (($env:PATH -split ';') -notcontains $d) { $env:PATH = "$d;$env:PATH" } }
    if ($bins) { Note "added to PATH for this session only -- to keep it, append: $($bins -join ';')" }
    return (Have g++)
}

# --------------------------------------------------------------------------- #
# Xyce (prebuilt binary from Sandia)
# --------------------------------------------------------------------------- #
# Sandia's download links are opaque per-file IDs on a WordPress site, so they name one
# exact build rather than "the latest". The SHA-256 pins that build: when Sandia publishes
# 7.11 the ID either 404s or serves a different file, and the hash check turns that into a
# loud failure instead of a silent stale install. No checksums are published upstream --
# this one was taken from the file itself, so re-verify by hand when bumping the version.
$XyceVersion = '7.10.0'
$XyceWinUrl  = 'https://xyce.sandia.gov/download/2000/'
$XyceWinSha  = 'E3073B6335A43ECC732F56944806005C99B1A378C7D21E293779AE05453B955A'
$XycePage    = 'https://xyce.sandia.gov/downloads/executables/'

function Install-Xyce {
    $tmp = Join-Path $env:TEMP ("xyce_" + [guid]::NewGuid().ToString('N'))
    New-Item -ItemType Directory -Path $tmp -Force | Out-Null
    $zip = Join-Path $tmp 'xyce.zip'
    try {
        Note "downloading Xyce $XyceVersion for Windows (~41 MB) ..."
        # A plain GET only: the download handler 301-redirects HEAD requests to the site
        # homepage, so any preflight check here would wrongly report failure.
        $prevProgress = $ProgressPreference
        $ProgressPreference = 'SilentlyContinue'   # the IWR progress bar makes big downloads crawl
        try { Invoke-WebRequest -Uri $XyceWinUrl -OutFile $zip -UseBasicParsing -TimeoutSec 900 }
        catch { Err "download failed: $($_.Exception.Message)"; Note "URL: $XyceWinUrl"; return $false }
        finally { $ProgressPreference = $prevProgress }

        $got = (Get-FileHash -Path $zip -Algorithm SHA256).Hash
        if ($got -ne $XyceWinSha) {
            Err 'checksum mismatch -- refusing to run the installer'
            Note "expected $XyceWinSha"
            Note "got      $got"
            Note "Sandia probably published a new release; verify the file yourself at $XycePage"
            Note 'then update $XyceWinUrl/$XyceWinSha in this script'
            return $false
        }
        Ok 'checksum verified'

        Expand-Archive -Path $zip -DestinationPath $tmp -Force
        $exe = Get-ChildItem -Path $tmp -Filter '*.exe' -Recurse | Select-Object -First 1
        if (-not $exe) { Err 'no .exe inside the archive'; return $false }

        # NSIS installer: /S is its silent switch. It needs elevation for Program Files.
        if (-not $IsAdmin) { Note 'the installer needs elevation -- expect a UAC prompt' }
        Note "running $($exe.Name) /S ..."
        $p = Start-Process -FilePath $exe.FullName -ArgumentList '/S' -Wait -PassThru -Verb RunAs
        if ($p.ExitCode -ne 0) { Err "installer exited with code $($p.ExitCode)"; return $false }
        $script:DidInstall = $true

        Update-PathFromRegistry
        if (Have Xyce) { Ok "installed Xyce: $((Get-Command Xyce).Source)"; return $true }

        # Silent NSIS runs do not always touch PATH; look where it lands by default.
        $bin = Get-ChildItem -Path @("$env:ProgramFiles", "${env:ProgramFiles(x86)}", $env:SystemDrive) `
                             -Filter 'Xyce*' -Directory -ErrorAction SilentlyContinue |
               ForEach-Object { Join-Path $_.FullName 'bin' } |
               Where-Object { Test-Path (Join-Path $_ 'Xyce.exe') } | Select-Object -First 1
        if ($bin) {
            $env:PATH = "$bin;$env:PATH"
            Ok "installed Xyce to $bin"
            Note "added to PATH for this session only -- to keep it, append: $bin"
            return $true
        }
        Warn 'installer finished but no Xyce.exe was found -- add its bin directory to PATH by hand'
        return $false
    } finally {
        Remove-Item $tmp -Recurse -Force -ErrorAction SilentlyContinue
    }
}

# Detect a tool by any of $Names; with -Install, try the package manager once and re-detect.
# Returns the resolved command name, or '' when still missing.
function Resolve-Tool {
    param([string[]]$Names, [string]$Label, [string]$Winget, [string]$Scoop, [string]$Choco)
    $found = $Names | Where-Object { Have $_ } | Select-Object -First 1
    # Nothing to offer without a package manager -- let the caller print its hint.
    if (-not $found -and $PKG -and (Confirm-Install "$Label is missing. Install it with $PKG`?")) {
        Install-Dep $Label $Winget $Scoop $Choco -Probe $Names[0] | Out-Null
        $found = $Names | Where-Object { Have $_ } | Select-Object -First 1
    }
    return [string]$found
}

Hdr 'Environment'
Ok "OS: Windows$(if ($PKG) { "  (package manager: $PKG)" })"
if ($NoInstall) {
    Note '(-NoInstall: reporting only, nothing will be installed)'
} elseif (-not $PKG) {
    if ($Install) { Warn '-Install given but no winget/scoop/choco found -- falling back to hints only' }
} elseif ([Console]::IsInputRedirected -and -not $Install) {
    Note '(input is redirected: reporting only -- pass -Install to install without prompts)'
} elseif ($PKG -eq 'choco' -and -not $IsAdmin) {
    Warn 'choco installs need an elevated shell -- re-run as Administrator if they fail'
} elseif ($PKG -eq 'scoop' -and $IsAdmin) {
    Warn 'scoop is not meant to run elevated -- re-run in a normal shell if installs fail'
}

# --------------------------------------------------------------------------- #
# core toolchain
# --------------------------------------------------------------------------- #
Hdr 'Core toolchain (required to build + run)'

# C++ compiler (MinGW/Clang; MSVC is not used by the Makefile).
# On winget the compiler and make both come from MSYS2, so that route is taken first.
$CXX = @('g++', 'clang++') | Where-Object { Have $_ } | Select-Object -First 1
if (-not $CXX -and $PKG -eq 'winget' -and
    (Confirm-Install 'No C++ compiler. Install MSYS2 and bootstrap gcc+make with pacman (several minutes)?')) {
    if (Install-Dep 'MSYS2' 'MSYS2.MSYS2' '' '' -Probe '') { Install-Msys2Toolchain | Out-Null }
    $CXX = @('g++', 'clang++') | Where-Object { Have $_ } | Select-Object -First 1
}
if (-not $CXX) { $CXX = Resolve-Tool @('g++', 'clang++') 'C++ compiler' '' 'gcc' 'mingw' }
if ($CXX) {
    Ok "C++ compiler: $CXX ($((& $CXX --version 2>$null | Select-Object -First 1)))"
} else {
    Err 'no MinGW/Clang C++ compiler (need g++ or clang++, C++17)'
    Note 'easiest: MSYS2 -> `pacman -S mingw-w64-ucrt-x86_64-gcc make`, then add its /ucrt64/bin to PATH'
    Hint 'MSYS2.MSYS2' 'gcc' 'mingw'
    $MissingRequired += 'C++ compiler'
}

# make (winget has no usable package -- it rides along with the MSYS2 bootstrap above)
$MAKE = Resolve-Tool @('make') 'make' '' 'make' 'make'
if ($MAKE) {
    Ok "make: $((Get-Command make).Source)"
} else {
    Err 'make not found'
    Note 'comes with MSYS2 (`pacman -S make`); or: scoop install make'
    Hint '' 'make' 'make'
    $MissingRequired += 'make'
}

# python (python.exe or the py launcher)
$PY = Resolve-Tool @('python', 'py') 'python' 'Python.Python.3.12' 'python' 'python'
if ($PY) {
    Ok "python: $((& $PY --version 2>&1))"
} else {
    Err 'python not found'
    Hint 'Python.Python.3.12' 'python' 'python'
    $MissingRequired += 'python'
}

# Xyce (external circuit solver; no package on any manager -- see Install-Xyce)
# The auto-yes flag is -InstallXyce, not -Install: a plain -Install must never kick
# off a multi-megabyte download and an elevated installer run on its own.
if (-not (Have Xyce) -and
    (Confirm-Install "Xyce is missing. Download Sandia's prebuilt $XyceVersion installer (~41 MB, needs elevation)?" $InstallXyce.IsPresent)) {
    Install-Xyce | Out-Null
}
if (Have Xyce) {
    Ok "Xyce: $((Get-Command Xyce).Source)"
} else {
    Err 'Xyce not on PATH -- the solver cannot run without it'
    if ($Install -and -not $InstallXyce) {
        Note 'no winget/scoop/choco package carries Xyce'
        Note "for Sandia's prebuilt installer, re-run with -InstallXyce"
    }
    Note 'download the Windows build from https://xyce.sandia.gov'
    Note "then add Xyce's bin directory to PATH"
    $MissingRequired += 'Xyce'
}

# --------------------------------------------------------------------------- #
# Python dependencies (numpy + matplotlib)
# --------------------------------------------------------------------------- #
Hdr 'Python dependencies (numpy, matplotlib)'
$PyBin = $PY
$useVenv = $PY -and -not $NoVenv
if ($PY) {
    if ($useVenv) {
        $PyBin = Join-Path $PSScriptRoot '.venv\Scripts\python.exe'
        if (-not (Test-Path $PyBin)) {
            Invoke-Native $PY @('-m', 'venv', '.venv') | Out-Null
            # `python -m venv` is a native call: a non-zero exit does not throw, so the
            # interpreter it should have produced is what actually decides success here.
            if (Test-Path $PyBin) {
                Ok 'created virtualenv .venv'
            } else {
                Err 'python -m venv failed (no .venv\Scripts\python.exe)'
                Warn 'falling back to the current Python'
                $useVenv = $false
                $PyBin = $PY
            }
        } else {
            Ok 'reusing existing .venv'
        }
    }
    if ($useVenv) {
        Invoke-Native $PyBin @('-m', 'pip', 'install', '--quiet', '--upgrade', 'pip') | Out-Null
        $code = Invoke-Native $PyBin @('-m', 'pip', 'install', '--quiet', 'numpy', 'matplotlib') -Show
        if ($code -eq 0) { Ok 'installed numpy + matplotlib into .venv' }
        else { Err 'pip install failed inside .venv'; $MissingRequired += 'numpy/matplotlib' }
    } else {
        $code = Invoke-Native $PY @('-c', 'import numpy, matplotlib')
        if ($code -eq 0) { Ok 'numpy + matplotlib already importable' }
        else {
            Warn 'numpy/matplotlib not importable in the current Python'
            $code = 1
            if (Confirm-Install 'Install numpy + matplotlib into the current Python?') {
                $code = Invoke-Native $PY @('-m', 'pip', 'install', '--quiet', 'numpy', 'matplotlib') -Show
            }
            if ($code -eq 0) { Ok 'installed numpy + matplotlib into the current Python' }
            else {
                Note "install: $PY -m pip install numpy matplotlib   (or re-run without -NoVenv)"
                $MissingRequired += 'numpy/matplotlib'
            }
        }
    }
} else {
    Warn 'skipping Python deps (no python)'
}

# --------------------------------------------------------------------------- #
# build main.exe
# --------------------------------------------------------------------------- #
Hdr 'Build main.exe'
if ($CXX -and $MAKE) {
    $log = Join-Path $env:TEMP 'wrcosim_build.log'
    # Makefile defaults CXX to clang++; override to the compiler we actually found.
    $code = Invoke-Native make @("CXX=$CXX") -LogFile $log
    if ($code -eq 0 -and (Test-Path .\main.exe)) {
        Ok 'built main.exe'
    } else {
        Err "build failed -- see $log"
        Get-Content $log -Tail 15 | ForEach-Object { Write-Host "    $_" }
        $MissingRequired += 'main.exe build'
    }
} else {
    Warn 'skipping build (need a MinGW/Clang compiler + make)'
    $MissingRequired += 'main.exe build'
}

# --------------------------------------------------------------------------- #
# optional: schematic (node + elkjs) and PDF/PNG (pdflatex + poppler)
# --------------------------------------------------------------------------- #
if (-not $NoOptional) {
    Hdr 'Optional: circuit schematic + PDF export'

    if (Resolve-Tool @('npm') 'node/npm' 'OpenJS.NodeJS' 'nodejs' 'nodejs') {
        $code = Invoke-Native npm @('install', '--silent')
        if ($code -eq 0) { Ok 'npm install (elkjs) done -- schematic layout available' }
        else { Warn 'npm install failed -- schematic layout unavailable (non-fatal)' }
    } else {
        Warn 'node/npm not found -- inline schematic unavailable (non-fatal)'
        Hint 'OpenJS.NodeJS' 'nodejs' 'nodejs'
    }

    if (Resolve-Tool @('pdflatex') 'pdflatex' 'MiKTeX.MiKTeX' 'latex' 'miktex') {
        Ok "pdflatex: $((Get-Command pdflatex).Source)  (.tex/.pdf export + schematic PDF)"
    } else {
        Warn 'pdflatex not found -- PDF export + schematic render unavailable (non-fatal)'
        Hint 'MiKTeX.MiKTeX' 'latex' 'miktex'
    }

    # poppler (pdftoppm) rasterizes the schematic PDF -> inline PNG. There is no `sips` on Windows,
    # so without pdftoppm the inline schematic PNG won't render (the .tex/.pdf export still works).
    if (Resolve-Tool @('pdftoppm') 'poppler' '' 'poppler' 'poppler') {
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
    if ($DidInstall) { Note 'some installs only take effect in a new shell -- open one, then re-run ./setup.ps1' }
    elseif ($PKG) { Note 'resolve the [MISS] items above (or re-run as ./setup.ps1 -Install), then re-run' }
    else { Note 'resolve the [MISS] items above, then re-run ./setup.ps1' }
}

Write-Host "`nRun the studio:" -ForegroundColor White
$pyShow = if ($useVenv) { '.\.venv\Scripts\python' } elseif ($PY) { $PY } else { 'python' }
Write-Host "  $pyShow sim_ui.py"
Write-Host '  then open' -ForegroundColor DarkGray -NoNewline; Write-Host ' http://127.0.0.1:8000'
Write-Host "`nHeadless sweep:" -ForegroundColor White
Write-Host "  $pyShow sim_ui.py sweep --param L_FEM --min 1e-7 --max 1e-5 --steps 8 --scale log"
Write-Host ''
