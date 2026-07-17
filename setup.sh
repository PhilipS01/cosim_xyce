#!/usr/bin/env bash
#
# setup.sh -- install / build the WR Co-Simulation Studio on macOS or Linux.
#
# Does what it safely can automatically:
#   * creates a Python virtualenv (.venv) and installs numpy + matplotlib into it
#   * builds ./main via make
#   * runs `npm install` for the optional schematic layout (if node is present)
# For system tools it cannot install without sudo (C++ compiler, make, Xyce,
# pdflatex, node) it checks presence and prints an OS-specific install hint.
#
# Usage:  ./setup.sh [--no-venv] [--no-optional] [-h|--help]
#
set -euo pipefail

cd "$(dirname "$0")"

# --------------------------------------------------------------------------- #
# options
# --------------------------------------------------------------------------- #
USE_VENV=1
DO_OPTIONAL=1
for arg in "$@"; do
    case "$arg" in
        --no-venv)     USE_VENV=0 ;;
        --no-optional) DO_OPTIONAL=0 ;;
        -h|--help)
            cat <<'EOF'
setup.sh -- install / build the WR Co-Simulation Studio on macOS or Linux.

Automatically:
  * create a Python virtualenv (.venv) with numpy + matplotlib
  * build ./main via make
  * npm install the optional schematic layout (if node is present)

For system tools it cannot install (C++ compiler, make, Xyce, pdflatex, node)
it checks presence and prints an OS-specific install hint.

Usage: ./setup.sh [--no-venv] [--no-optional] [-h|--help]
  --no-venv       install numpy/matplotlib into the current Python (no .venv)
  --no-optional   skip the optional schematic/PDF tools (node, pdflatex)
EOF
            exit 0 ;;
        *) echo "unknown option: $arg (see --help)" >&2; exit 2 ;;
    esac
done

# --------------------------------------------------------------------------- #
# pretty output
# --------------------------------------------------------------------------- #
if [ -t 1 ] && [ -z "${NO_COLOR:-}" ]; then
    B=$'\033[1m'; G=$'\033[32m'; Y=$'\033[33m'; R=$'\033[31m'; D=$'\033[2m'; Z=$'\033[0m'
else
    B=; G=; Y=; R=; D=; Z=
fi
hdr()  { printf '\n%s==> %s%s\n' "$B" "$*" "$Z"; }
ok()   { printf '  %s[ok]%s  %s\n'   "$G" "$Z" "$*"; }
warn() { printf '  %s[warn]%s %s\n'  "$Y" "$Z" "$*"; }
err()  { printf '  %s[MISS]%s %s\n'  "$R" "$Z" "$*"; }
note() { printf '  %s%s%s\n' "$D" "$*" "$Z"; }

have() { command -v "$1" >/dev/null 2>&1; }

MISSING_REQUIRED=()   # collected for the final summary

# --------------------------------------------------------------------------- #
# OS + package manager detection
# --------------------------------------------------------------------------- #
OS="$(uname -s)"
PKG=""
case "$OS" in
    Darwin) PKG="brew" ;;
    Linux)
        for m in apt-get dnf pacman zypper; do have "$m" && { PKG="$m"; break; }; done
        ;;
    *) warn "unsupported OS '$OS' -- script targets macOS and Linux; continuing best-effort." ;;
esac

# hint <tool> <brew-pkg> <apt-pkg> <dnf-pkg> <pacman-pkg> <zypper-pkg>
# prints the install command for the detected package manager.
hint() {
    local tool="$1" brew="$2" apt="$3" dnf="$4" pac="$5" zyp="$6"
    case "$PKG" in
        brew)    note "install: brew install $brew" ;;
        apt-get) note "install: sudo apt-get install -y $apt" ;;
        dnf)     note "install: sudo dnf install -y $dnf" ;;
        pacman)  note "install: sudo pacman -S --needed $pac" ;;
        zypper)  note "install: sudo zypper install -y $zyp" ;;
        *)       note "install '$tool' with your system package manager" ;;
    esac
}

hdr "Environment"
ok "OS: $OS${PKG:+  (package manager: $PKG)}"

# --------------------------------------------------------------------------- #
# core toolchain
# --------------------------------------------------------------------------- #
hdr "Core toolchain (required to build + run)"

# C++ compiler
CXX=""
for c in c++ g++ clang++; do have "$c" && { CXX="$c"; break; }; done
if [ -n "$CXX" ]; then
    ok "C++ compiler: $CXX ($("$CXX" --version 2>/dev/null | head -1))"
else
    err "no C++ compiler (need g++ or clang++, C++17)"
    if [ "$OS" = Darwin ]; then
        note "install: xcode-select --install   (Command Line Tools)"
    else
        hint "g++" "gcc" "g++ make" "gcc-c++ make" "base-devel" "gcc-c++ make"
    fi
    MISSING_REQUIRED+=("C++ compiler")
fi

# make
if have make; then
    ok "make: $(command -v make)"
else
    err "make not found"
    [ "$OS" = Darwin ] && note "install: xcode-select --install" \
                       || hint "make" "make" "make" "make" "make" "make"
    MISSING_REQUIRED+=("make")
fi

# python3
if have python3; then
    ok "python3: $(python3 --version 2>&1)"
else
    err "python3 not found"
    hint "python3" "python" "python3 python3-venv python3-pip" "python3 python3-pip" "python" "python3 python3-pip"
    MISSING_REQUIRED+=("python3")
fi

# Xyce (external circuit solver; cannot auto-install)
if have Xyce; then
    ok "Xyce: $(command -v Xyce)"
else
    err "Xyce not on PATH -- the solver cannot run without it"
    note "build/download from https://xyce.sandia.gov  (no standard brew/apt package)"
    note "then ensure the 'Xyce' binary is on your PATH"
    MISSING_REQUIRED+=("Xyce")
fi

# --------------------------------------------------------------------------- #
# Python dependencies (numpy + matplotlib)
# --------------------------------------------------------------------------- #
hdr "Python dependencies (numpy, matplotlib)"
PYBIN="python3"
if have python3; then
    if [ "$USE_VENV" = 1 ]; then
        if [ ! -d .venv ]; then
            if python3 -m venv .venv 2>/dev/null; then
                ok "created virtualenv .venv"
            else
                err "python3 -m venv failed (missing venv module?)"
                hint "python3-venv" "python" "python3-venv" "python3" "python" "python3"
                warn "falling back to the system Python"
                USE_VENV=0
            fi
        else
            ok "reusing existing .venv"
        fi
    fi
    if [ "$USE_VENV" = 1 ]; then
        PYBIN=".venv/bin/python"
        "$PYBIN" -m pip install --quiet --upgrade pip >/dev/null 2>&1 || true
        if "$PYBIN" -m pip install --quiet numpy matplotlib; then
            ok "installed numpy + matplotlib into .venv"
        else
            err "pip install failed inside .venv"
            MISSING_REQUIRED+=("numpy/matplotlib")
        fi
    else
        # no venv: check imports, install into current python only if importable pip exists
        if python3 -c 'import numpy, matplotlib' 2>/dev/null; then
            ok "numpy + matplotlib already importable"
        else
            warn "numpy/matplotlib not importable in the system Python"
            note "install: python3 -m pip install numpy matplotlib   (or re-run without --no-venv)"
            MISSING_REQUIRED+=("numpy/matplotlib")
        fi
    fi
else
    warn "skipping Python deps (no python3)"
fi

# --------------------------------------------------------------------------- #
# build ./main
# --------------------------------------------------------------------------- #
hdr "Build ./main"
if [ -n "$CXX" ] && have make; then
    # Makefile defaults CXX to clang++; override to the compiler we actually found.
    if make CXX="$CXX" >/tmp/wrcosim_build.log 2>&1; then
        ok "built ./main"
    else
        err "build failed -- see /tmp/wrcosim_build.log"
        tail -n 15 /tmp/wrcosim_build.log | sed 's/^/    /'
        MISSING_REQUIRED+=("./main build")
    fi
else
    warn "skipping build (need a C++ compiler + make)"
fi

# --------------------------------------------------------------------------- #
# optional: schematic (node + elkjs) and PDF export (pdflatex)
# --------------------------------------------------------------------------- #
if [ "$DO_OPTIONAL" = 1 ]; then
    hdr "Optional: circuit schematic + PDF export"

    if have npm; then
        if npm install --silent >/dev/null 2>&1; then
            ok "npm install (elkjs) done -- schematic layout available"
        else
            warn "npm install failed -- schematic layout unavailable (non-fatal)"
        fi
    else
        warn "node/npm not found -- inline schematic unavailable (non-fatal)"
        hint "node" "node" "nodejs npm" "nodejs npm" "nodejs npm" "nodejs npm"
    fi

    if have pdflatex; then
        ok "pdflatex: $(command -v pdflatex)  (schematic PNG + .tex/.pdf export)"
    else
        warn "pdflatex not found -- schematic render + PDF export unavailable (non-fatal)"
        if [ "$OS" = Darwin ]; then
            note "install: brew install --cask basictex   (then: sudo tlmgr install circuitikz)"
        else
            hint "pdflatex" "" "texlive-latex-extra texlive-pictures" "texlive-scheme-medium" "texlive-most" "texlive-latex"
        fi
    fi
else
    note "(skipping optional schematic/PDF tools: --no-optional)"
fi

# --------------------------------------------------------------------------- #
# summary
# --------------------------------------------------------------------------- #
hdr "Summary"
if [ "${#MISSING_REQUIRED[@]}" -eq 0 ]; then
    ok "all required components are in place"
else
    warn "missing required components: ${MISSING_REQUIRED[*]}"
    note "resolve the [MISS] items above, then re-run ./setup.sh"
fi

printf '\n%sRun the studio:%s\n' "$B" "$Z"
if [ "$USE_VENV" = 1 ]; then
    printf '  source .venv/bin/activate && python sim_ui.py\n'
    printf '  %s# or without activating:%s  .venv/bin/python sim_ui.py\n' "$D" "$Z"
else
    printf '  python3 sim_ui.py\n'
fi
printf '  %sthen open%s http://127.0.0.1:8000\n' "$D" "$Z"
printf '\n%sHeadless sweep:%s\n' "$B" "$Z"
printf '  %s sim_ui.py sweep --param L_FEM --min 1e-7 --max 1e-5 --steps 8 --scale log\n' \
    "$([ "$USE_VENV" = 1 ] && echo '.venv/bin/python' || echo python3)"
echo
