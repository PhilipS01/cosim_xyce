#!/usr/bin/env bash
#
# setup.sh -- install / build the WR Co-Simulation Studio on macOS or Linux.
#
# Does what it safely can automatically:
#   * creates a Python virtualenv (.venv) and installs numpy + matplotlib into it
#   * builds ./main via make
#   * runs `npm install` for the optional schematic layout (if node is present)
# For system tools (C++ compiler, make, Xyce, pdflatex, poppler, node) it checks
# presence and prints an OS-specific install hint; with --install it runs those
# install commands itself (via sudo where the package manager needs it).
#
# Usage:  ./setup.sh [--install] [--no-venv] [--no-optional] [-h|--help]
#
set -euo pipefail

cd "$(dirname "$0")"

# --------------------------------------------------------------------------- #
# options
# --------------------------------------------------------------------------- #
USE_VENV=1
DO_OPTIONAL=1
DO_INSTALL=0
DO_INSTALL_XYCE=0
ASSUME_NO=0
for arg in "$@"; do
    case "$arg" in
        --install)      DO_INSTALL=1 ;;
        --install-xyce) DO_INSTALL_XYCE=1 ;;
        --no-install)   ASSUME_NO=1 ;;
        --no-venv)     USE_VENV=0 ;;
        --no-optional) DO_OPTIONAL=0 ;;
        -h|--help)
            cat <<'EOF'
setup.sh -- install / build the WR Co-Simulation Studio on macOS or Linux.

Automatically:
  * create a Python virtualenv (.venv) with numpy + matplotlib
  * build ./main via make
  * npm install the optional schematic layout (if node is present)

For system tools (C++ compiler, make, Xyce, pdflatex, poppler, node) it checks
presence and prints an OS-specific install hint.

On a terminal it asks before installing anything that is missing. Piped or
redirected input (CI) never prompts -- it just reports, as --no-install does.

Usage: ./setup.sh [--install] [--install-xyce] [--no-install]
                  [--no-venv] [--no-optional]
  --install       answer yes to every install prompt: install the missing system
                  tools with the detected package manager (brew / apt-get / dnf /
                  pacman / zypper) without asking. Uses sudo where needed; on
                  macOS the compiler comes from `xcode-select --install`.
                  Does NOT cover Xyce -- see --install-xyce.
  --install-xyce  answer yes to the Xyce prompt: download Sandia's prebuilt
                  binary and install it (macOS Arm64 only; the checksum is
                  pinned in this script). Kept out of --install because it is a
                  ~18 MB download plus a sudo installer run. Other platforms get
                  pointed at spack/source.
  --no-install    never install and never prompt; only report what is missing
  --no-venv       install numpy/matplotlib into the current Python (no .venv)
  --no-optional   skip the optional schematic/PDF tools (node, pdflatex, poppler)
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

# ask <question> [auto-yes-flag]
# Returns 0 for yes. The auto-yes flag defaults to --install, so `--install` answers
# every prompt yes; --no-install answers no. Anything other than an interactive
# terminal (pipes, CI, `curl | bash`) answers no rather than blocking on input.
ask() {
    local q="$1" auto="${2:-$DO_INSTALL}" reply
    [ "$auto" = 1 ] && return 0
    [ "$ASSUME_NO" = 1 ] && return 1
    [ -t 0 ] || return 1
    printf '  %s%s%s [y/N] ' "$Y" "$q" "$Z"
    read -r reply || { echo; return 1; }
    case "$reply" in [Yy]|[Yy][Ee][Ss]) return 0 ;; *) return 1 ;; esac
}

MISSING_REQUIRED=()   # collected for the final summary

# --------------------------------------------------------------------------- #
# OS + package manager detection
# --------------------------------------------------------------------------- #
OS="$(uname -s)"
PKG=""
case "$OS" in
    # brew is only claimed when it is actually installed -- otherwise --install
    # would try to run a command that isn't there.
    Darwin) have brew && PKG="brew" ;;
    Linux)
        for m in apt-get dnf pacman zypper; do have "$m" && { PKG="$m"; break; }; done
        ;;
    *) warn "unsupported OS '$OS' -- script targets macOS and Linux; continuing best-effort." ;;
esac

SUDO=""
[ "$(id -u)" != 0 ] && have sudo && SUDO="sudo"
APT_UPDATED=0
DID_INSTALL=0   # set once anything was actually installed (drives the final hint)

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

# install_pkgs <label> <brew> <apt> <dnf> <pacman> <zypper>
# runs the detected manager's install command. Callers must use it in a
# condition (if / ||) so that a failing install does not trip `set -e`.
install_pkgs() {
    local label="$1" brew="$2" apt="$3" dnf="$4" pac="$5" zyp="$6" pkgs=""
    case "$PKG" in
        brew)    pkgs="$brew" ;;
        apt-get) pkgs="$apt"  ;;
        dnf)     pkgs="$dnf"  ;;
        pacman)  pkgs="$pac"  ;;
        zypper)  pkgs="$zyp"  ;;
        *) warn "cannot install $label -- no supported package manager detected"; return 1 ;;
    esac
    if [ -z "$pkgs" ]; then
        warn "cannot install $label with $PKG -- no package for it"
        return 1
    fi
    if [ "$PKG" != brew ] && [ -z "$SUDO" ] && [ "$(id -u)" != 0 ]; then
        warn "cannot install $label -- need root and 'sudo' is not available"
        return 1
    fi
    note "installing $label via $PKG ($pkgs) ..."
    if [ "$PKG" = apt-get ] && [ "$APT_UPDATED" = 0 ]; then
        $SUDO apt-get update -qq || true
        APT_UPDATED=1
    fi
    # shellcheck disable=SC2086  # $pkgs is a deliberate multi-package word list
    case "$PKG" in
        brew)    brew install $pkgs ;;
        apt-get) $SUDO apt-get install -y $pkgs ;;
        dnf)     $SUDO dnf install -y $pkgs ;;
        pacman)  $SUDO pacman -S --needed --noconfirm $pkgs ;;
        zypper)  $SUDO zypper install -y $pkgs ;;
    esac || return 1
    DID_INSTALL=1
}

# --------------------------------------------------------------------------- #
# Xyce (prebuilt binary from Sandia)
# --------------------------------------------------------------------------- #
# Sandia's download links are opaque per-file IDs on a WordPress site, so they name
# one exact build rather than "the latest". The SHA-256 below pins that build: when
# Sandia publishes 7.11 the ID either 404s or serves a different file, and the hash
# check turns that into a loud failure instead of a silent stale install.
# No checksums are published upstream -- these were taken from the files themselves,
# so re-verify by hand when bumping the version.
XYCE_VERSION="7.10.0"
XYCE_MAC_URL="https://xyce.sandia.gov/download/1952/"
XYCE_MAC_SHA256="06cd438e0a3d4dc86948c1d0f6684fa125db3dfa1776ff1d9339ebb77865754b"
XYCE_PAGE="https://xyce.sandia.gov/downloads/executables/"

xyce_manual_note() {
    note "options: spack install xyce   (source build, takes hours)"
    note "         or grab a binary yourself from $XYCE_PAGE"
    note "         or build from https://github.com/Xyce/Xyce"
}

install_xyce() {
    if [ "$OS" != Darwin ]; then
        warn "no scripted Xyce install for $OS"
        note "Sandia ships only a RHEL8 rpm; their own page says it does not work on"
        note "Debian/Ubuntu (so not in WSL2 either) and is broken on Fedora"
        xyce_manual_note
        return 1
    fi
    local arch; arch="$(uname -m)"
    if [ "$arch" != arm64 ]; then
        warn "Sandia's macOS build is Arm64-only -- nothing to install on $arch"
        xyce_manual_note
        return 1
    fi
    if [ -z "$SUDO" ] && [ "$(id -u)" != 0 ]; then
        warn "installing Xyce needs root and 'sudo' is not available"
        return 1
    fi

    local tmp zip pkg got
    tmp="$(mktemp -d)" || return 1
    zip="$tmp/xyce.zip"

    note "downloading Xyce $XYCE_VERSION for macOS Arm64 (~18 MB) ..."
    # A plain GET only: the download handler 301-redirects HEAD requests to the
    # site homepage, so any preflight check here would wrongly report failure.
    if ! curl -fL --progress-bar --max-time 900 -o "$zip" "$XYCE_MAC_URL"; then
        err "download failed"
        note "URL: $XYCE_MAC_URL"
        rm -rf "$tmp"; return 1
    fi

    got="$(shasum -a 256 "$zip" | cut -d' ' -f1)"
    if [ "$got" != "$XYCE_MAC_SHA256" ]; then
        err "checksum mismatch -- refusing to install"
        note "expected $XYCE_MAC_SHA256"
        note "got      $got"
        note "Sandia probably published a new release; verify the file yourself at $XYCE_PAGE"
        note "then update XYCE_MAC_URL/XYCE_MAC_SHA256 in this script"
        rm -rf "$tmp"; return 1
    fi
    ok "checksum verified"

    if ! unzip -q -o "$zip" -d "$tmp"; then
        err "unzip failed"; rm -rf "$tmp"; return 1
    fi
    pkg="$(find "$tmp" -name '*.pkg' -maxdepth 2 2>/dev/null | head -1)"
    if [ -z "$pkg" ]; then
        err "no .pkg inside the archive"; rm -rf "$tmp"; return 1
    fi

    note "installing $(basename "$pkg") (sudo installer) ..."
    if ! $SUDO installer -pkg "$pkg" -target / ; then
        err "installer failed"; rm -rf "$tmp"; return 1
    fi
    rm -rf "$tmp"

    # The pkg drops a versioned directory straight into /usr/local.
    local bin; bin="$(find /usr/local -maxdepth 2 -type d -name bin -path '*Xyce*' 2>/dev/null | head -1)"
    if [ -n "$bin" ] && [ -x "$bin/Xyce" ]; then
        export PATH="$bin:$PATH"
        ok "installed Xyce to $bin"
        note "added to PATH for this session only -- to keep it, append: $bin"
    else
        warn "installed, but no Xyce binary found under /usr/local -- add its bin/ to PATH by hand"
    fi
    # Sandia does not notarize these builds; Gatekeeper blocks the first run.
    warn "macOS will block the first run: approve Xyce under"
    warn "System Settings -> Privacy & Security, then re-run ./setup.sh"
    have Xyce
}

# ensure_tool <probe-cmd> <label> <brew> <apt> <dnf> <pacman> <zypper>
# succeeds when the tool is callable, offering to install it first.
ensure_tool() {
    local probe="$1" label="$2"
    have "$probe" && return 0
    # Nothing to offer without a package manager -- let the caller print its hint.
    [ -n "$PKG" ] || return 1
    ask "$label is missing. Install it with $PKG?" || return 1
    shift 2
    install_pkgs "$label" "$@" || return 1
    have "$probe"
}

hdr "Environment"
ok "OS: $OS${PKG:+  (package manager: $PKG)}"
if [ "$ASSUME_NO" = 1 ]; then
    note "(--no-install: reporting only, nothing will be installed)"
elif [ -z "$PKG" ]; then
    if [ "$DO_INSTALL" = 1 ]; then warn "--install given but no supported package manager found -- hints only"; fi
elif [ ! -t 0 ] && [ "$DO_INSTALL" != 1 ]; then
    note "(input is not a terminal: reporting only -- pass --install to install without prompts)"
elif [ "$PKG" != brew ] && [ -n "$SUDO" ]; then
    note "installs run through sudo -- you may be prompted for your password"
fi

# --------------------------------------------------------------------------- #
# core toolchain
# --------------------------------------------------------------------------- #
hdr "Core toolchain (required to build + run)"

# C++ compiler. On macOS this is the Command Line Tools, not a brew package:
# `xcode-select --install` opens a GUI installer and returns immediately, so the
# compiler only shows up on a later run.
CXX=""
for c in c++ g++ clang++; do have "$c" && { CXX="$c"; break; }; done
if [ -z "$CXX" ]; then
    if [ "$OS" = Darwin ]; then
        if ask "No C++ compiler. Request the Xcode Command Line Tools?"; then
            note "a GUI installer opens ..."
            xcode-select --install 2>/dev/null || true
            DID_INSTALL=1
            note "finish that dialog, then re-run ./setup.sh"
        fi
    elif [ -n "$PKG" ] && ask "No C++ compiler. Install one with $PKG?"; then
        install_pkgs "C++ compiler" "gcc" "g++ make" "gcc-c++ make" "base-devel" "gcc-c++ make" || true
    fi
    for c in c++ g++ clang++; do have "$c" && { CXX="$c"; break; }; done
fi
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

# make (on macOS it rides along with the Command Line Tools)
HAVE_MAKE=0
if [ "$OS" = Darwin ]; then
    have make && HAVE_MAKE=1
elif ensure_tool make "make" "make" "make" "make" "make" "make"; then
    HAVE_MAKE=1
fi
if [ "$HAVE_MAKE" = 1 ]; then
    ok "make: $(command -v make)"
else
    err "make not found"
    if [ "$OS" = Darwin ]; then note "install: xcode-select --install"
    else hint "make" "make" "make" "make" "make" "make"; fi
    MISSING_REQUIRED+=("make")
fi

# python3
if ensure_tool python3 "python3" "python" "python3 python3-venv python3-pip" \
                       "python3 python3-pip" "python" "python3 python3-pip"; then
    ok "python3: $(python3 --version 2>&1)"
else
    err "python3 not found"
    hint "python3" "python" "python3 python3-venv python3-pip" "python3 python3-pip" "python" "python3 python3-pip"
    MISSING_REQUIRED+=("python3")
fi

# Xyce (external circuit solver; no package manager carries it -- see install_xyce)
if ! have Xyce; then
    if [ "$DO_INSTALL_XYCE" = 1 ]; then
        install_xyce || true
    elif [ "$OS" = Darwin ] && [ "$(uname -m)" = arm64 ]; then
        # Only offer where a scripted install exists. The auto-yes flag is
        # --install-xyce, not --install: a plain --install must never kick off a
        # multi-megabyte download and an installer run on its own.
        if ask "Xyce is missing. Download Sandia's prebuilt $XYCE_VERSION build (~18 MB, needs sudo)?" "$DO_INSTALL_XYCE"; then
            install_xyce || true
        fi
    fi
fi
if have Xyce; then
    ok "Xyce: $(command -v Xyce)"
else
    err "Xyce not on PATH -- the solver cannot run without it"
    if [ "$DO_INSTALL" = 1 ] && [ "$DO_INSTALL_XYCE" != 1 ]; then
        note "no brew/apt/dnf/pacman/zypper package carries Xyce"
        note "for Sandia's prebuilt macOS binary, re-run with --install-xyce"
    fi
    note "build/download from https://xyce.sandia.gov"
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
        # Judge the venv by its interpreter, not by the directory: a half-created
        # .venv/ would otherwise be "reused" and every pip call would fail.
        if [ ! -x .venv/bin/python ]; then
            if ! python3 -m venv .venv 2>/dev/null; then
                # Debian/Ubuntu ship the venv module in a separate package.
                warn "python3 -m venv failed -- the venv module looks missing"
                if [ -n "$PKG" ] && ask "Install the python venv module with $PKG?"; then
                    install_pkgs "python venv module" "python" "python3-venv" "python3" "python" "python3" || true
                    python3 -m venv .venv 2>/dev/null || true
                fi
            fi
            if [ -x .venv/bin/python ]; then
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
        # no venv: check imports, then install into the current python if asked
        if python3 -c 'import numpy, matplotlib' 2>/dev/null; then
            ok "numpy + matplotlib already importable"
        else
            warn "numpy/matplotlib not importable in the system Python"
            if ask "Install numpy + matplotlib into the system Python?" \
               && python3 -m pip install --quiet numpy matplotlib; then
                ok "installed numpy + matplotlib into the system Python"
            else
                note "install: python3 -m pip install numpy matplotlib   (or re-run without --no-venv)"
                MISSING_REQUIRED+=("numpy/matplotlib")
            fi
        fi
    fi
else
    warn "skipping Python deps (no python3)"
fi

# --------------------------------------------------------------------------- #
# build ./main
# --------------------------------------------------------------------------- #
hdr "Build ./main"
if [ -n "$CXX" ] && [ "$HAVE_MAKE" = 1 ]; then
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
    MISSING_REQUIRED+=("./main build")
fi

# --------------------------------------------------------------------------- #
# optional: schematic (node + elkjs), PDF export (pdflatex), raster (poppler)
# --------------------------------------------------------------------------- #
if [ "$DO_OPTIONAL" = 1 ]; then
    hdr "Optional: circuit schematic + PDF export"

    if ensure_tool npm "node/npm" "node" "nodejs npm" "nodejs npm" "nodejs npm" "nodejs npm"; then
        if npm install --silent >/dev/null 2>&1; then
            ok "npm install (elkjs) done -- schematic layout available"
        else
            warn "npm install failed -- schematic layout unavailable (non-fatal)"
        fi
    else
        warn "node/npm not found -- inline schematic unavailable (non-fatal)"
        hint "node" "node" "nodejs npm" "nodejs npm" "nodejs npm" "nodejs npm"
    fi

    # pdflatex: on macOS it is a cask, and circuitikz still needs a tlmgr step.
    PDFLATEX_OK=0
    if have pdflatex; then
        PDFLATEX_OK=1
    elif [ "$OS" = Darwin ] && [ "$PKG" = brew ]; then
        # brew has no pdflatex formula, so this branch is the cask or nothing --
        # falling through to ensure_tool would ask a second time for a package
        # that does not exist.
        if ask "pdflatex is missing. Install BasicTeX with brew (~90 MB)?"; then
            note "installing BasicTeX via brew (cask) ..."
            if brew install --cask basictex; then DID_INSTALL=1; fi
            # BasicTeX lands in /Library/TeX/texbin, which is not on PATH until re-login.
            [ -x /Library/TeX/texbin/pdflatex ] && export PATH="/Library/TeX/texbin:$PATH"
            have pdflatex && PDFLATEX_OK=1
            [ "$PDFLATEX_OK" = 1 ] && note "still needed once: sudo tlmgr update --self && sudo tlmgr install circuitikz"
        fi
    elif ensure_tool pdflatex "pdflatex" "" "texlive-latex-extra texlive-pictures" \
                     "texlive-scheme-medium" "texlive-most" "texlive-latex"; then
        PDFLATEX_OK=1
    fi
    if [ "$PDFLATEX_OK" = 1 ]; then
        ok "pdflatex: $(command -v pdflatex)  (schematic PNG + .tex/.pdf export)"
    else
        warn "pdflatex not found -- schematic render + PDF export unavailable (non-fatal)"
        if [ "$OS" = Darwin ]; then
            note "install: brew install --cask basictex   (then: sudo tlmgr install circuitikz)"
        else
            hint "pdflatex" "" "texlive-latex-extra texlive-pictures" "texlive-scheme-medium" "texlive-most" "texlive-latex"
        fi
    fi

    # poppler (pdftoppm) rasterizes the schematic PDF at 300 dpi. macOS can fall back
    # to the built-in `sips`, so there it is a quality upgrade rather than a hard need.
    if ensure_tool pdftoppm "poppler" "poppler" "poppler-utils" "poppler-utils" "poppler" "poppler-tools"; then
        ok "pdftoppm (poppler): $(command -v pdftoppm)  (300 dpi inline schematic PNG)"
    elif [ "$OS" = Darwin ]; then
        warn "pdftoppm not found -- falling back to the built-in sips (lower quality, non-fatal)"
        hint "pdftoppm" "poppler" "poppler-utils" "poppler-utils" "poppler" "poppler-tools"
    else
        warn "pdftoppm (poppler) not found -- inline schematic PNG unavailable (non-fatal)"
        hint "pdftoppm" "poppler" "poppler-utils" "poppler-utils" "poppler" "poppler-tools"
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
    if [ "$DID_INSTALL" = 1 ]; then
        note "some installs only take effect in a new shell -- open one, then re-run ./setup.sh"
    elif [ -n "$PKG" ]; then
        note "resolve the [MISS] items above (or re-run as ./setup.sh --install), then re-run"
    else
        note "resolve the [MISS] items above, then re-run ./setup.sh"
    fi
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
