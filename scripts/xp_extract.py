#!/usr/bin/env python3
r"""Extract the circuit PORT IMPEDANCE x_P(f) with a Xyce .AC sweep.

x_P(f) is the driving-point impedance the field sees looking INTO the circuit at the
field/circuit interface, with the circuit's own independent sources zeroed:

    x_P(f) = V_port(f) / (1 A injected at the port)

The field branch is NOT part of x_P. The WR deck's interface (Vmeas ammeter, the matched-secant
Bfield, R_ROM/L_ROM) is therefore absent here; a 1 A AC current source is injected at the port node
in its place. What survives from circuit_spec.txt:

    R / L / C            kept verbatim (the AC deck drops the IC= the transient deck carries)
    V-sources (any kind) ZEROED -> emitted as a 0 V source, i.e. a short that keeps the node names
    I-sources (any kind) REMOVED -> an open circuit
    SW (time-gated)      frozen to a plain resistor: Ron if closed at --switch-time, else Roff

Nothing downstream refers to an output COLUMN INDEX: the port voltage is fetched from the .prn
header by the netlist node name (VR(<port>) / VI(<port>)), so renaming or reordering print tokens
cannot silently shift the result.

Companion: scripts/rho_contraction.py turns the CSV written here into the WR contraction factor
rho(f). It is a separate script on purpose -- the ROM sweep is re-runnable without touching Xyce.

CLI:
    python3 scripts/xp_extract.py --fstart 1 --fstop 1e7 --points 200 --out results/xp.csv
"""
import argparse
import csv
import os
import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

# Xyce SWITCH model defaults (Xyce Reference Guide); emitCustomTopology leaves them out of the
# .MODEL when the spec omits them, so the frozen resistor has to fall back to the same numbers.
SW_RON_DEFAULT = 1.0
SW_ROFF_DEFAULT = 1.0e12

# SPICE magnitude suffixes, as accepted by Xyce (and validated by emitCustomTopology's VPWL guard).
_SUFFIX = {"t": 1e12, "g": 1e9, "k": 1e3, "m": 1e-3, "u": 1e-6,
           "n": 1e-9, "p": 1e-12, "f": 1e-15}


def spice_float(tok):
    """Parse a spec value token ('0.5', '5e-6', '63n', '1meg') into a float.

    'meg'/'mil' are checked before the single-letter suffixes so '1meg' is not read as 1 milli."""
    s = str(tok).strip().lower()
    m = re.match(r"^([+-]?(?:\d+\.?\d*|\.\d+)(?:e[+-]?\d+)?)\s*([a-z]*)$", s)
    if not m:
        raise ValueError(f"not a number: {tok!r}")
    val, suf = float(m.group(1)), m.group(2)
    if not suf:
        return val
    if suf.startswith("meg"):
        return val * 1e6
    if suf.startswith("mil"):
        return val * 25.4e-6
    return val * _SUFFIX.get(suf[0], 1.0)


class Element(object):
    """One circuit_spec.txt line, classified for the AC deck.

    kind is the ROLE in an AC analysis, not the spec TYPE: 'R'/'L'/'C' passives, 'V'/'I' for any
    independent source (its waveform is irrelevant once zeroed/removed) and 'SW'."""

    def __init__(self, kind, name, a, b, toks):
        self.kind, self.name, self.a, self.b, self.toks = kind, name, a, b, toks

    @property
    def nodes(self):
        return (self.a, self.b)


_SRC_V = ("VSIN", "VDC", "VPULSE", "VPWM", "VPWL")
_SRC_I = ("ISIN", "IDC", "IPULSE", "IPWM", "IPWL")


def parse_spec(text):
    """circuit_spec.txt text -> [Element]. Unknown TYPEs raise, matching emitCustomTopology."""
    elems = []
    for lineno, line in enumerate(text.splitlines(), 1):
        s = line.split("#", 1)[0].split("*", 1)[0].strip()
        if not s:
            continue
        tok = s.split()
        typ = tok[0].upper()
        if len(tok) < 4:
            raise ValueError(f"circuit_spec line {lineno}: '{typ}' needs at least <name> <a> <b>")
        name, a, b = tok[1], tok[2], tok[3]
        rest = tok[4:]
        if typ in ("R", "L", "C"):
            if not rest:
                raise ValueError(f"circuit_spec line {lineno}: {typ} {name} has no value")
            elems.append(Element(typ, name, a, b, rest))
        elif typ in _SRC_V:
            elems.append(Element("V", name, a, b, rest))
        elif typ in _SRC_I:
            elems.append(Element("I", name, a, b, rest))
        elif typ == "SW":
            if len(rest) < 2:
                raise ValueError(f"circuit_spec line {lineno}: SW {name} needs tclose topen")
            elems.append(Element("SW", name, a, b, rest))
        else:
            raise ValueError(f"circuit_spec line {lineno}: unknown TYPE '{typ}'")
    if not elems:
        raise ValueError("circuit_spec has no elements")
    return elems


def switch_resistance(e, t):
    """Frozen resistance of a SW element at time t, plus a flag for 'caught mid-transition'.

    Mirrors emitCustomTopology's gate: closed during [tclose, topen), with a smooth ramp of width
    tau at each edge. A gate strictly between 0 and 1 means t lands ON a ramp, where the real Xyce
    generic switch is interpolating between Roff and Ron -- an AC analysis of a half-thrown switch
    is not a circuit anybody means, so we snap to the nearer state and hand back the warning."""
    tclose = spice_float(e.toks[0])
    topen = spice_float(e.toks[1])
    ron = spice_float(e.toks[2]) if len(e.toks) > 2 else SW_RON_DEFAULT
    roff = spice_float(e.toks[3]) if len(e.toks) > 3 else SW_ROFF_DEFAULT
    tau = spice_float(e.toks[4]) if len(e.toks) > 4 else 1.0e-5
    if tau <= 0.0:
        tau = 1.0e-9
    up = 1.0 if tclose <= 0.0 else min(max((t - tclose) / tau, 0.0), 1.0)
    dn = min(max((t - topen) / tau, 0.0), 1.0)
    gate = up - dn
    mid = 1e-9 < gate < 1.0 - 1e-9
    return (ron if gate >= 0.5 else roff), mid


def build_ac_netlist(elems, port="p", gnd="0", sweep="dec", points=200,
                     fstart=1.0, fstop=1.0e6, switch_time=0.0, prn="xp_ac.prn"):
    """Emit the AC probe deck. Returns (netlist_text, [warnings])."""
    if fstart <= 0.0 or fstop <= 0.0:
        raise ValueError("AC sweep frequencies must be > 0 (a log sweep cannot start at DC)")
    if fstop <= fstart:
        raise ValueError("--fstop must be greater than --fstart")
    if points < 2:
        raise ValueError("--points must be >= 2")
    sweep = sweep.lower()
    if sweep not in ("dec", "oct", "lin"):
        raise ValueError("--sweep must be dec, oct or lin")

    warn, body = [], []
    for e in elems:
        if e.kind in ("R", "L", "C"):
            # IC= is transient-only state; an AC deck must not carry it.
            body.append(f"{e.kind}{e.name} {e.a} {e.b} {e.toks[0]}")
        elif e.kind == "V":
            # Zeroed, not deleted: a 0 V source is the short the superposition argument calls for
            # and it keeps {e.a}/{e.b} as distinct named nodes, so probes/checks still resolve.
            body.append(f"V{e.name} {e.a} {e.b} 0")
        elif e.kind == "I":
            body.append(f"* I{e.name} {e.a} {e.b} removed (independent current source -> open)")
        elif e.kind == "SW":
            r, mid = switch_resistance(e, switch_time)
            state = "closed" if r <= 1.0 else "open"
            if mid:
                warn.append(f"switch '{e.name}' is mid-transition at t={switch_time:g} s; "
                            f"frozen to its nearer state ({state}). Pick a settled time.")
            # Xyce rejects a trailing '$'/';' comment on a device line -> own '*' line above it.
            body.append(f"* SW {e.name} frozen {state} at t={switch_time:g} s")
            body.append(f"R_sw_{e.name} {e.a} {e.b} {r:.10g}")

    if not any(port in e.nodes for e in elems if e.kind != "I"):
        # Nothing but removed current sources touches the port, so x_P is an open circuit and the
        # AC matrix is empty. Xyce's own message for this ("Empty matrix has been found") says
        # nothing about the cause, so refuse here instead.
        raise ValueError(
            f"port node '{port}' is touched only by independent current sources, which the x_P "
            f"probe removes -- the port is then an open circuit (x_P = infinity) and there is "
            f"nothing to sweep. Add the passives the field actually loads.")

    head = [
        "* x_P probe deck -- generated by scripts/xp_extract.py, DO NOT EDIT BY HAND",
        "* Circuit-side elements from circuit_spec.txt with all independent sources zeroed/removed;",
        f"* the field branch is replaced by a 1 A AC injection into the port node '{port}'.",
        "",
        f"I_xp_inj {gnd} {port} AC 1 0",
        "",
    ]
    tail = [
        "",
        f".AC {sweep.upper()} {points} {fstart:.10g} {fstop:.10g}",
        # VR/VI = real/imaginary part of the node voltage (Xyce Reference Guide, AC output
        # operators). Named by NODE, never by column position.
        f".PRINT AC FILE={prn} FORMAT=STD VR({port}) VI({port})",
        ".END",
        "",
    ]
    return "\n".join(head + body + tail), warn


def _find_col(header, name):
    """Column index of `name` in a .prn header, case- and space-insensitive. -1 if absent."""
    want = name.replace(" ", "").upper()
    for i, h in enumerate(header):
        if h.replace(" ", "").upper() == want:
            return i
    return -1


def read_ac_prn(path, port="p"):
    """Read a Xyce AC .prn -> (f, x_P complex). Columns are located BY NAME from the header."""
    with open(path) as fh:
        header = fh.readline().split()
        ifreq = _find_col(header, "FREQ")
        ire = _find_col(header, f"VR({port})")
        iim = _find_col(header, f"VI({port})")
        if min(ifreq, ire, iim) < 0:
            raise ValueError(f"{os.path.basename(path)}: header {header} lacks FREQ/VR({port})/VI({port})")
        need = max(ifreq, ire, iim) + 1
        f, re_, im_ = [], [], []
        for line in fh:
            s = line.strip()
            if not s or s.lower().startswith("end"):
                continue
            parts = s.split()
            if len(parts) < need:
                continue
            try:
                f.append(float(parts[ifreq]))
                re_.append(float(parts[ire]))
                im_.append(float(parts[iim]))
            except ValueError:
                continue
    if not f:
        raise ValueError(f"{os.path.basename(path)}: no AC data rows")
    return np.array(f), np.array(re_) + 1j * np.array(im_)


def dc_port_resistance(elems, port="p", gnd="0", switch_time=0.0):
    """Analytic DC limit of x_P: solve the resistive network with L shorted and C open.

    Independent of Xyce, so it is a genuine cross-check of the AC sweep's low-frequency end rather
    than a restatement of it. Zeroed voltage sources and inductors are ideal shorts, so their nodes
    are merged (union-find) before the nodal solve. Returns None when the port has no DC path to
    ground (the conductance matrix is then singular -- x_P is an open circuit at DC)."""
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    find(port); find(gnd)
    conds = []                                    # (node_a, node_b, G) after shorts are merged
    for e in elems:
        if e.kind == "I":
            continue
        find(e.a); find(e.b)
        if e.kind in ("L", "V"):                  # ideal short at DC
            union(e.a, e.b)
        elif e.kind == "R":
            if spice_float(e.toks[0]) <= 0.0:
                union(e.a, e.b)
        elif e.kind == "SW":
            if switch_resistance(e, switch_time)[0] <= 0.0:
                union(e.a, e.b)
    for e in elems:
        if e.kind == "R":
            r = spice_float(e.toks[0])
            if r > 0.0:
                conds.append((find(e.a), find(e.b), 1.0 / r))
        elif e.kind == "SW":
            r = switch_resistance(e, switch_time)[0]
            if r > 0.0:
                conds.append((find(e.a), find(e.b), 1.0 / r))
        # C is an open circuit at DC -> contributes nothing.

    ref = find(gnd)
    if find(port) == ref:
        return 0.0                                # port shorted to ground through L/V
    # '-' binds tighter than '|', so the reference node is removed from the UNION, not just
    # from the port singleton -- hence the explicit parentheses.
    unk = sorted(({n for tri in conds for n in tri[:2]} | {find(port)}) - {ref})
    idx = {n: i for i, n in enumerate(unk)}
    G = np.zeros((len(unk), len(unk)))
    for a, b, g in conds:
        if a == b:
            continue
        if a in idx:
            G[idx[a], idx[a]] += g
        if b in idx:
            G[idx[b], idx[b]] += g
        if a in idx and b in idx:
            G[idx[a], idx[b]] -= g
            G[idx[b], idx[a]] -= g
    rhs = np.zeros(len(unk))
    rhs[idx[find(port)]] = 1.0                    # 1 A injected at the port
    try:
        v = np.linalg.solve(G, rhs)
    except np.linalg.LinAlgError:
        return None
    return float(v[idx[find(port)]])


def port_shunt_capacitance(elems, port="p", gnd="0"):
    """Total capacitance sitting DIRECTLY across the port (port <-> ground); 0.0 if none.

    Capacitors in parallel add, so a sum is the right reduction here."""
    return sum(spice_float(e.toks[0]) for e in elems
               if e.kind == "C" and {e.a, e.b} == {port, gnd})


def sanity_checks(f, xp, elems, port="p", gnd="0", switch_time=0.0, rtol=0.05):
    """Low- and high-frequency asymptotics of |x_P|. Returns (ok, [message lines]).

    LOW f : |x_P| -> the network's DC resistance (for a series-R port, R_s), computed independently
            by dc_port_resistance().
    HIGH f: a capacitor directly across the port shorts it, so |x_P| must ROLL OFF as 1/f. Tested as
            the log-log slope over the top decade (-1 = capacitive), NOT as |x_P| ~ 0: at a finite
            f_stop a 1/(2*pi*f*C) tail is still a perfectly visible number."""
    msgs, ok = [], True
    mag = np.abs(xp)

    rdc = dc_port_resistance(elems, port, gnd, switch_time)
    if rdc is None:
        msgs.append("low f : no DC path from the port to ground -- skipped")
    else:
        got = mag[int(np.argmin(f))]
        tol = max(rtol * abs(rdc), 1e-12)
        good = abs(got - rdc) <= tol
        ok &= good
        msgs.append(f"low f : |x_P({f.min():.4g} Hz)| = {got:.6g} Ohm vs DC network resistance "
                    f"{rdc:.6g} Ohm -- {'OK' if good else 'MISMATCH'}")
        if not good and rdc != 0.0:
            msgs.append("        (raise --fstart? a series L still dominates if 2*pi*f*L >> R there)")

    cap = port_shunt_capacitance(elems, port, gnd)
    if cap > 0.0:
        # The asymptote is the capacitor's OWN impedance: everything else at the port is in
        # parallel with it, so |x_P| can only sit at or below |Z_C| once the cap dominates.
        fmax = float(f.max())
        got = float(mag[int(np.argmax(f))])
        zc = 1.0 / (2.0 * np.pi * fmax * cap)
        if got <= 1.25 * zc:
            msgs.append(f"high f: |x_P({fmax:.4g} Hz)| = {got:.6g} Ohm vs the port capacitor's "
                        f"|Z_C| = {zc:.6g} Ohm -- OK (the shunt C shorts the port)")
        else:
            # NOT counted as a failure: the user picked the band, and a sweep that stops below the
            # port resonance simply has not reached the roll-off yet.
            msgs.append(f"high f: |x_P({fmax:.4g} Hz)| = {got:.6g} Ohm still exceeds the port "
                        f"capacitor's |Z_C| = {zc:.6g} Ohm -- INCONCLUSIVE, the sweep has not "
                        f"reached the capacitive asymptote (raise --fstop past the port resonance)")
        top = f >= fmax / 10.0
        if top.sum() >= 2 and np.all(mag[top] > 0):
            slope = float(np.polyfit(np.log10(f[top]), np.log10(mag[top]), 1)[0])
            msgs.append(f"        |x_P| log-log slope over the top decade = {slope:+.3f} "
                        f"(-1 once the shunt C dominates)")
    else:
        msgs.append("high f: no capacitor directly across the port -- no roll-off expected, skipped")
    return ok, msgs


def run_xyce(netlist, workdir, prn="xp_ac.prn", xyce="Xyce", timeout=300):
    """Write + solve the AC deck in `workdir`. Returns the .prn path; raises with Xyce's tail."""
    deck = os.path.join(workdir, "xp_probe.cir")
    with open(deck, "w") as fh:
        fh.write(netlist)
    proc = subprocess.run([xyce, os.path.basename(deck)], cwd=workdir,
                          capture_output=True, text=True, timeout=timeout)
    out = os.path.join(workdir, prn)
    if proc.returncode != 0 or not os.path.exists(out):
        tail = "\n".join(((proc.stderr or "") + (proc.stdout or "")).splitlines()[-25:])
        raise RuntimeError(f"Xyce failed on the x_P deck (exit {proc.returncode}):\n{tail}")
    return out


def extract_xp(spec_text, port="p", gnd="0", sweep="dec", points=200, fstart=1.0,
               fstop=1.0e6, switch_time=0.0, xyce="Xyce", keep_dir=None, deck_path=None):
    """Spec text -> {"f", "xp", "elements", "netlist", "warnings", "checks_ok", "checks"}.

    Runs in a temp directory unless keep_dir is given, so the working tree's WR outputs (which use
    the same .prn naming conventions) are never touched. deck_path (if given) receives a copy of
    the generated netlist BEFORE Xyce runs, so a failed solve still leaves the deck to inspect."""
    elems = parse_spec(spec_text)
    netlist, warn = build_ac_netlist(elems, port, gnd, sweep, points, fstart, fstop, switch_time)
    if deck_path:
        os.makedirs(os.path.dirname(os.path.abspath(deck_path)) or ".", exist_ok=True)
        with open(deck_path, "w") as fh:
            fh.write(netlist)
    d = keep_dir or tempfile.mkdtemp(prefix="xp_ac_")
    try:
        os.makedirs(d, exist_ok=True)
        f, xp = read_ac_prn(run_xyce(netlist, d, xyce=xyce), port)
    finally:
        if keep_dir is None:
            shutil.rmtree(d, ignore_errors=True)
    ok, msgs = sanity_checks(f, xp, elems, port, gnd, switch_time)
    return {"f": f, "xp": xp, "elements": elems, "netlist": netlist,
            "warnings": warn, "checks_ok": ok, "checks": msgs}


XP_CSV_HEADER = ["f_Hz", "Re_xP_ohm", "Im_xP_ohm", "abs_xP_ohm"]


def write_xp_csv(path, f, xp):
    d = os.path.dirname(os.path.abspath(path))
    if d:
        os.makedirs(d, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(XP_CSV_HEADER)
        for fi, zi in zip(f, xp):
            w.writerow([f"{fi:.10g}", f"{zi.real:.10g}", f"{zi.imag:.10g}", f"{abs(zi):.10g}"])
    return os.path.abspath(path)


def read_xp_csv(path):
    """Inverse of write_xp_csv -> (f, x_P complex). Columns located by HEADER NAME."""
    with open(path, newline="") as fh:
        rows = list(csv.reader(fh))
    if len(rows) < 2:
        raise ValueError(f"{path}: no data rows")
    hdr = [h.strip() for h in rows[0]]
    try:
        i_f, i_re, i_im = (hdr.index(c) for c in XP_CSV_HEADER[:3])
    except ValueError:
        raise ValueError(f"{path}: header {hdr} is not an x_P CSV ({XP_CSV_HEADER})")
    f, z = [], []
    for r in rows[1:]:
        if len(r) <= max(i_f, i_re, i_im):
            continue
        f.append(float(r[i_f]))
        z.append(complex(float(r[i_re]), float(r[i_im])))
    return np.array(f), np.array(z)


def plot_xp(f, xp, path=None, title="Circuit port impedance"):
    """|x_P| vs f on a log-x axis. Returns the matplotlib figure."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 3.2))
    ax.loglog(f, np.abs(xp), lw=1.2, label=r"$|x_P(f)|$")
    ax.set_xlabel("frequency (Hz)")
    ax.set_ylabel(r"$|x_P|$ ($\Omega$)")
    ax.set_title(title)
    ax.grid(True, which="both", alpha=.3)
    ax.legend(fontsize=8, loc="upper right")
    if path:
        os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
        fig.savefig(path, dpi=110, bbox_inches="tight")
    return fig


def main(argv=None):
    ap = argparse.ArgumentParser(description="Extract the circuit port impedance x_P(f) via Xyce .AC")
    ap.add_argument("--spec", default=os.path.join(ROOT, "circuit_spec.txt"),
                    help="circuit spec to probe (default: the project's circuit_spec.txt)")
    ap.add_argument("--port", default="p", help="port node name (default: p)")
    ap.add_argument("--gnd", default="0", help="ground node name (default: 0)")
    ap.add_argument("--sweep", default="dec", choices=["dec", "oct", "lin"],
                    help="dec/oct = points PER DECADE/OCTAVE, lin = TOTAL points (default: dec)")
    ap.add_argument("--points", type=int, default=200, help="points per decade (or total, if --sweep lin)")
    ap.add_argument("--fstart", type=float, default=1.0, help="start frequency in Hz (> 0)")
    ap.add_argument("--fstop", type=float, default=1.0e6, help="stop frequency in Hz")
    ap.add_argument("--switch-time", type=float, default=0.0, dest="switch_time",
                    help="time at which SW elements are frozen to Ron/Roff (default: 0)")
    ap.add_argument("--out", default=os.path.join(ROOT, "results", "xp.csv"), help="CSV to write")
    ap.add_argument("--plot", default=None, help="also save a |x_P| vs f PNG here")
    ap.add_argument("--deck", default=None, help="also save the generated AC netlist here")
    ap.add_argument("--xyce", default="Xyce", help="Xyce executable (default: Xyce on PATH)")
    ap.add_argument("--no-check", action="store_true", help="skip the asymptotic sanity checks")
    a = ap.parse_args(argv)

    with open(a.spec) as fh:
        spec_text = fh.read()
    res = extract_xp(spec_text, a.port, a.gnd, a.sweep, a.points, a.fstart, a.fstop,
                     a.switch_time, a.xyce, deck_path=a.deck)
    if a.deck:
        print(f"deck  -> {os.path.abspath(a.deck)}")
    out = write_xp_csv(a.out, res["f"], res["xp"])
    print(f"x_P   -> {out}  ({res['f'].size} points, {a.fstart:g}..{a.fstop:g} Hz, {a.sweep})")
    if a.plot:
        plot_xp(res["f"], res["xp"], a.plot)
        print(f"plot  -> {os.path.abspath(a.plot)}")
    for w in res["warnings"]:
        print(f"warn: {w}", file=sys.stderr)
    if not a.no_check:
        for m in res["checks"]:
            print(f"check: {m}")
        if not res["checks_ok"]:
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
