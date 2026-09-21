#!/usr/bin/env python3
r"""WR contraction factor rho(f) from an already-extracted port impedance x_P(f).

    beta(f) = 1 / (R_ROM  + j*2*pi*f*L_ROM)        ROM (preconditioner) admittance
    Y(f)    = 1 / (R_FEM  + j*2*pi*f*L_FEM)        true field admittance
    rho(f)  = (1 + beta*x_P)^-1 * (beta - Y) * x_P

|rho| < 1 over the band the coupled system actually excites is the a-priori statement that the
waveform-relaxation iteration contracts there; |rho| = 1 is the break-even locus. rho vanishes
identically where beta == Y, i.e. when the ROM matches the field exactly -- the nilpotent case.

Deliberately separate from scripts/xp_extract.py: x_P depends only on the CIRCUIT, so the ROM sweep
below can be re-run with new R_ROM/L_ROM/R_FEM/L_FEM against the same CSV without touching Xyce.

CLI:
    python3 scripts/rho_contraction.py --xp results/xp.csv --plot results/rho.png
    python3 scripts/rho_contraction.py --xp results/xp.csv --l-rom 1e-3 --r-rom 1.0
"""
import argparse
import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

from xp_extract import read_xp_csv          # noqa: E402  (same-directory companion script)

RHO_CSV_HEADER = ["f_Hz", "Re_xP_ohm", "Im_xP_ohm", "abs_xP_ohm",
                  "Re_rho", "Im_rho", "abs_rho"]


def read_sim_config(path):
    """Parse sim_config.txt ('key = value', '#' comments) into {key: float}."""
    cfg = {}
    if not os.path.exists(path):
        return cfg
    with open(path) as fh:
        for line in fh:
            s = line.split("#", 1)[0].strip()
            if "=" not in s:
                continue
            k, v = (x.strip() for x in s.split("=", 1))
            try:
                cfg[k] = float(v)
            except ValueError:
                pass
    return cfg


def impedance(f, r, l):
    """Series R + jwL impedance over a frequency array."""
    return r + 1j * 2.0 * np.pi * np.asarray(f, dtype=float) * l


def contraction(f, xp, r_rom, l_rom, r_fem, l_fem):
    """rho(f) = (1 + beta*x_P)^-1 (beta - Y) x_P. Returns (rho, beta, Y)."""
    beta = 1.0 / impedance(f, r_rom, l_rom)
    y = 1.0 / impedance(f, r_fem, l_fem)
    rho = (beta - y) * xp / (1.0 + beta * xp)
    return rho, beta, y


def unity_crossings(f, rho):
    """Frequencies where |rho| crosses 1, interpolated in log f / log|rho|.

    Log-log interpolation because both axes are swept and plotted logarithmically -- a linear
    interpolation between two decades apart would put the crossing in the wrong place."""
    m = np.abs(rho)
    out = []
    good = (m > 0) & np.isfinite(m) & (f > 0)
    lf, lm = np.log10(f[good]), np.log10(m[good])
    for i in range(len(lm) - 1):
        a, b = lm[i], lm[i + 1]
        if a == 0.0:
            out.append(10.0 ** lf[i])
        elif (a < 0.0) != (b < 0.0):
            out.append(10.0 ** (lf[i] + (0.0 - a) * (lf[i + 1] - lf[i]) / (b - a)))
    if len(lm) and lm[-1] == 0.0:
        out.append(10.0 ** lf[-1])
    return np.array(out)


def write_rho_csv(path, f, xp, rho):
    d = os.path.dirname(os.path.abspath(path))
    if d:
        os.makedirs(d, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(RHO_CSV_HEADER)
        for fi, zi, ri in zip(f, xp, rho):
            w.writerow([f"{fi:.10g}", f"{zi.real:.10g}", f"{zi.imag:.10g}", f"{abs(zi):.10g}",
                        f"{ri.real:.10g}", f"{ri.imag:.10g}", f"{abs(ri):.10g}"])
    return os.path.abspath(path)


def plot_rho(f, xp, rho, path=None, crossings=None, title="Port impedance and WR contraction"):
    r"""|x_P| and |rho| against the SAME log-frequency axis.

    They carry different units (Ohm vs dimensionless), so |rho| gets its own right-hand scale; the
    break-even level |rho| = 1 is drawn on that scale and every crossing is marked and labelled."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if crossings is None:
        crossings = unity_crossings(f, rho)
    mag = np.abs(rho)

    fig, ax = plt.subplots(figsize=(8, 3.6))
    ax.loglog(f, np.abs(xp), lw=1.2, color="tab:blue", label=r"$|x_P(f)|$")
    ax.set_xlabel("frequency (Hz)")
    ax.set_ylabel(r"$|x_P|$ ($\Omega$)", color="tab:blue")
    ax.tick_params(axis="y", labelcolor="tab:blue")
    ax.grid(True, which="both", alpha=.3)

    ax2 = ax.twinx()
    ax2.set_yscale("log")
    ax2.plot(f, np.maximum(mag, 1e-16), lw=1.2, color="tab:red", label=r"$|\rho(f)|$")
    ax2.axhline(1.0, ls="--", lw=1.0, color="tab:red", alpha=.6)
    ax2.set_ylabel(r"$|\rho|$  (contraction)", color="tab:red")
    ax2.tick_params(axis="y", labelcolor="tab:red")
    for k, fc in enumerate(crossings):
        ax2.plot([fc], [1.0], "o", ms=7, mfc="none", mew=1.6, color="k", zorder=5)
        # Stagger above/below: two crossings a fraction of a decade apart would otherwise print
        # their labels on top of each other.
        ax2.annotate(rf"$|\rho|=1$ @ {fc:.4g} Hz", (fc, 1.0), textcoords="offset points",
                     xytext=(6, 8 if k % 2 == 0 else -16), fontsize=8, zorder=6,
                     bbox=dict(fc="white", ec="none", alpha=.7, pad=1.2))
    lines = ax.get_lines()[:1] + ax2.get_lines()[:1]
    ax.legend(lines, [l.get_label() for l in lines], fontsize=8, loc="upper right")
    ax.set_title(title)
    if path:
        os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
        fig.savefig(path, dpi=110, bbox_inches="tight")
    return fig


def summarize(f, rho, crossings=None):
    """Human-readable a-priori verdict lines for the CLI and the studio panel."""
    if crossings is None:
        crossings = unity_crossings(f, rho)
    mag = np.abs(rho)
    i = int(np.argmax(mag))
    msgs = [f"max |rho| = {mag[i]:.6g} at f = {f[i]:.6g} Hz",
            f"|rho| at f_start ({f[0]:.4g} Hz) = {mag[0]:.6g}   "
            f"at f_stop ({f[-1]:.4g} Hz) = {mag[-1]:.6g}"]
    if mag[i] < 1.0:
        msgs.append(f"|rho| < 1 across the whole band -- WR contracts a priori "
                    f"(worst-case factor {mag[i]:.4g} per iteration).")
    else:
        band = f"{f[mag >= 1.0].min():.4g} .. {f[mag >= 1.0].max():.4g} Hz"
        msgs.append(f"|rho| >= 1 over {band} -- no a-priori contraction there; the run converges "
                    f"only if the circuit carries no content in that band.")
    if crossings.size:
        msgs.append("|rho| = 1 crossing(s): " + ", ".join(f"{c:.6g} Hz" for c in crossings))
    else:
        msgs.append("|rho| = 1 is never crossed in the swept band.")
    return msgs


def main(argv=None):
    cfg = read_sim_config(os.path.join(ROOT, "sim_config.txt"))
    ap = argparse.ArgumentParser(
        description="WR contraction factor rho(f) from an x_P CSV (no Xyce re-run)")
    ap.add_argument("--xp", default=os.path.join(ROOT, "results", "xp.csv"),
                    help="x_P CSV written by scripts/xp_extract.py")
    ap.add_argument("--r-rom", type=float, default=cfg.get("R_ROM", 4.59e-4), dest="r_rom",
                    help="ROM resistance (default: R_ROM from sim_config.txt)")
    ap.add_argument("--l-rom", type=float, default=cfg.get("L_ROM", 1.44e-7), dest="l_rom",
                    help="ROM inductance (default: L_ROM from sim_config.txt)")
    ap.add_argument("--r-fem", type=float, default=cfg.get("R_FEM", 5.1e-4), dest="r_fem",
                    help="true field resistance (default: R_FEM from sim_config.txt)")
    ap.add_argument("--l-fem", type=float, default=cfg.get("L_FEM", 1.6e-7), dest="l_fem",
                    help="true field inductance (default: L_FEM from sim_config.txt)")
    ap.add_argument("--out", default=os.path.join(ROOT, "results", "rho.csv"), help="CSV to write")
    ap.add_argument("--plot", default=None, help="also save the |x_P| + |rho| PNG here")
    a = ap.parse_args(argv)

    f, xp = read_xp_csv(a.xp)
    rho, _beta, _y = contraction(f, xp, a.r_rom, a.l_rom, a.r_fem, a.l_fem)
    cross = unity_crossings(f, rho)
    out = write_rho_csv(a.out, f, xp, rho)
    print(f"ROM   : R_ROM = {a.r_rom:g} Ohm, L_ROM = {a.l_rom:g} H")
    print(f"field : R_FEM = {a.r_fem:g} Ohm, L_FEM = {a.l_fem:g} H")
    print(f"rho   -> {out}  ({f.size} points)")
    if a.plot:
        plot_rho(f, xp, rho, a.plot, cross)
        print(f"plot  -> {os.path.abspath(a.plot)}")
    for m in summarize(f, rho, cross):
        print(m)
    return 0


if __name__ == "__main__":
    sys.exit(main())
