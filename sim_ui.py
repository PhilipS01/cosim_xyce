#!/usr/bin/env python3
"""
Local web UI to drive the waveform-relaxation field/circuit co-simulation.

Lets you tune the solver parameters (impedances, source, time-stepping, coupling
grids, WR settings), run ./main, and view the resulting waveforms and WR
convergence -- all in the browser, fully offline.

Usage:
    python3 sim_ui.py            # then open http://127.0.0.1:8000
    python3 sim_ui.py --port 9000

Requires: Python 3, matplotlib, numpy. No other dependencies.
The backend writes sim_config.txt, runs ./main (building it first if missing),
parses the .prn outputs, and renders plots as PNGs.
"""

import argparse
import base64
import io
import json
import os
import re
import subprocess
import sys
import time
import traceback
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

import matplotlib
matplotlib.use("Agg")  # headless: no GUI backend needed
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

# --- Parameter spec: (key, label, default, kind, slider min/max/step or None) ---
# kind: "float" or "int". slider tuple => render a range slider alongside the box.
PARAMS = [
    # The circuit side (source + passives + topology) is authored in the circuit-netlist text panel,
    # not here -- see the "Circuit netlist" box. Only field/coupling/timing/solver knobs live here.
    ("L_ROM",                           "L_ROM (H)",                    1.44e-7,  "float", None),
    ("R_ROM",                           "R_ROM (Ohm)",                  4.59e-4,  "float", None),
    ("L_FEM",                           "L_FEM (H, 'true' field)",      1.6e-7,   "float", None),
    ("R_FEM",                           "R_FEM (Ohm, 'true' field)",    5.1e-4,   "float", None),
    ("nonlin_model",                    "FEM nonlinearity",             0,        "choice",
        {0: "linear", 1: "magnetic saturation"}),
    ("I_sat",                           "I_sat saturation current (A)", 100.0,    "float", None),
    ("time_mode",                       "Run duration",                 0,        "choice",
        {0: "source periods", 1: "absolute end time"}),
    ("t_end",                           "End time t_end (s)",           2.0e-2,   "float", None),
    ("N_field_windows",                 "Field windows (total)",        50,       "int",   (1, 400, 1)),
    ("N_periods",                       "Number of source periods",     1,        "int",   (1, 10, 1)),
    ("N_field_steps_per_source_period", "Field steps / source period",  50,       "int",   (2, 200, 1)),
    ("N_field_eval_intervals",          "FEM eval intervals / window",  1,        "int",   (1, 64, 1)),
    ("N_xyce_coupling_intervals",       "Xyce coupling intervals",      100,      "int",   (2, 400, 1)),
    ("WRmaxSteps",                      "WR max iterations",            20,       "int",   (1, 100, 1)),
    ("WR_tolerance",                    "WR tolerance",                 1.0e-3,   "float", None),
    ("wr_convergence_method",           "WR convergence metric",        1,        "choice",
        {0: "waveform L1 (this code)", 1: "terminal scalar (reference)"}),
    ("coupling_mode",                   "Coupling direction",           0,        "choice",
        {0: "voltage-driven", 1: "current-driven"}),
    ("reconstruct_mode",                "Field reconstruction",         0,        "choice",
        {0: "pointwise (secant)", 1: "linear ramp", 2: "average (linear+const)",
         3: "pointwise (central diff)"}),
    ("interface_form",                  "Interface stamping",           0,        "choice",
        {0: "Thevenin (V source)", 1: "Norton (I source)"}),
]
DEFAULTS = {k: d for (k, _l, d, _kind, _s) in PARAMS}
KINDS = {k: kind for (k, _l, _d, kind, _s) in PARAMS}
LABELS = {k: l for (k, l, _d, _kind, _s) in PARAMS}
# Numeric params are sweepable (a "choice" metric switch is not a continuum).
SWEEPABLE = [k for (k, _l, _d, kind, _s) in PARAMS if kind in ("float", "int")]

# Circuit-side presets (hybrid model): a preset seeds the editable primitive fields; the user
# may then tweak any field. Increment 1 covers the three source kinds + series R/L; presets 4-6
# (switches) arrive in increment 2. Keys map to the flat config the C++ generator consumes.
# Presets seed the circuit-netlist text (circuit_spec) + the field/coupling/timing knobs. The circuit
# side is a simplified spec: reserved nodes p=port, 0=gnd; one element/line (VSIN/ISIN/VPULSE/R/L/C and
# SW name a b tclose topen [Ron Roff trise]).
PRESETS = {
    "P1: Sine V + RL": {"time_mode": 0, "coupling_mode": 0,
                        "circuit_spec": "VSIN Bemf s 0 1 50\nR Rs s cm0 6e-3\nL Ls cm0 p 1.6e-7\n"},
    # Bare current source directly on the port: series R/L/C are meaningless for a current drive
    # (the current is forced regardless) and an ideal I-source in series with L is degenerate.
    "P2: Sine I (bare)": {"time_mode": 0, "coupling_mode": 1,
                          "circuit_spec": "ISIN Bemf 0 p 1 50\n"},
    # window 1 straddles the whole ramp edge (stiff transient) -> more WR iters (WRmaxSteps=40)
    "P3: Step/ramp V + RL": {"time_mode": 1, "t_end": 2.0e-2, "N_field_windows": 50,
                             "coupling_mode": 0, "WRmaxSteps": 40,
                             "circuit_spec": "VPULSE Vemf s 0 0 1 0 1e-4\nR Rs s cm0 6e-3\nL Ls cm0 p 1.6e-7\n"},
    # Switch circuits (SW = time-gated resistor, closed during [tclose, topen)). NOTE C/Ron/times are
    # numerical-survival defaults (C small enough for WR to contract; the cap<->coil freewheel needs a
    # damped closed switch Ron~10 or its ~undamped LC ring dt-collapses) -- tune to your field.
    # 2-way (SPDT) switch: wiper w throws between the source branch (node a) and a short branch (node b),
    # both rejoining at the port p; cap w->0. Freewheel path is w-b-short-p (electrically w->p; tiny R).
    "P4: 2-way switch (sine U, C)": {"time_mode": 1, "t_end": 2.0e-2, "N_field_windows": 50,
        "coupling_mode": 0, "WRmaxSteps": 40,
        "circuit_spec": "VSIN Bemf a p 1 50\nC Csw w 0 1e-6\n"
                        "SW drv w a 0 6e-3 10 1e9 1e-5\nSW fw w b 6e-3 1e30 10 1e9 1e-5\nR shrt b p 1e-6\n"},
    "P5: 2-way switch (DC U, C)": {"time_mode": 1, "t_end": 2.0e-2, "N_field_windows": 50,
        "coupling_mode": 0, "WRmaxSteps": 40,
        "circuit_spec": "VDC Vemf a p 1\nC Csw w 0 1e-6\n"
                        "SW drv w a 0 6e-3 10 1e9 1e-5\nSW fw w b 6e-3 1e30 10 1e9 1e-5\nR shrt b p 1e-6\n"},
    "P6: 2-way switch (AC vs R)": {"time_mode": 1, "t_end": 2.0e-2, "N_field_windows": 50,
        "coupling_mode": 0, "WRmaxSteps": 40,
        "circuit_spec": "VSIN Bemf p bac 1 50\nR Rload p br 1e4\n"
                        "SW ac bac 0 0 6e-3 1e-3 1e9 1e-5\nSW rd br 0 6e-3 1e30 1e-3 1e9 1e-5\n"},
}


# ---------------------------------------------------------------------------
# Running the solver
# ---------------------------------------------------------------------------
# Conditional visibility: key -> list of AND-condition dicts; a control is shown iff ANY dict fully
# matches the current control values (OR-of-ANDs). Keys absent here are always visible. The custom
# spec box + SVG editor are handled separately in JS (visible only when circuit_kind == 4).
VISIBLE_WHEN = {
    "N_periods":                       [{"time_mode": [0]}],
    "N_field_steps_per_source_period": [{"time_mode": [0]}],
    "t_end":            [{"time_mode": [1]}],
    "N_field_windows":  [{"time_mode": [1]}],
    "I_sat": [{"nonlin_model": [1]}],
}


# Hover help: key -> HTML shown in a tooltip next to the control's label (a "?" icon).
HELP = {
    "circuit_spec_edit": (
        "<div class='hh'>Circuit side (text)</div>"
        "<div class='hn'>The circuit is authored here, not in the side panel. Presets seed it; edit and "
        "<b>Apply &amp; render</b>. Nodes: <code>p</code>=port (field attaches here), <code>0</code>=ground; "
        "reserved. One element per line:<br>"
        "<code>R/L/C name a b value</code><br>"
        "<code>VSIN/ISIN name a b amp freq</code> &middot; <code>VDC/IDC name a b value</code><br>"
        "<code>VPULSE/IPULSE name a b v1 v2 td tr</code> &middot; <code>VPWL/IPWL name a b t1 v1 t2 v2 …</code><br>"
        "<code>SW name a b tclose topen [Ron Roff trise]</code> &mdash; time-gated switch, closed during "
        "[tclose, topen) (use a big topen e.g. <code>1e30</code> to stay closed to the end; Ron=10, "
        "Roff=1e9, trise=1e-5 default).<br>The WR interface (<code>Vmeas</code> ammeter + <code>Bfield</code> "
        "field ROM) and the <code>.INCLUDE</code>/<code>.print</code>/<code>.end</code> directives are "
        "generated and shown locked around the editable box.</div>"
    ),
    "lcapy_export": (
        "<div class='hh'>LaTeX / PDF export</div>"
        "<div class='hn'>Downloads the schematic as <code>circuit.tex</code> (a standalone circuitikz "
        "document) + a compiled <code>circuit.pdf</code> for thesis figures &mdash; the same render shown "
        "on the right. Layout is by <b>ELK</b> (orthogonal placement + wire routing); it's a <b>physical "
        "view</b> &mdash; the interface is the ammeter (<code>Vmeas</code>) + the field ROM source "
        "(<code>Bfield</code>), and the WR PWL signal carriers are dropped. Needs <code>node</code>+"
        "<code>elkjs</code> (<code>npm install</code>) for the layout and <code>pdflatex</code>+"
        "<code>circuitikz</code> for the PDF; without pdflatex you still get the <code>.tex</code>.</div>"
    ),
    "reconstruct_mode": (
        "<div class='hh'>Field reconstruction</div>"
        "<div class='hn'>How the dummy field reconstructs its output waveform within a window &mdash; the "
        "field <b>current</b> (voltage-driven) or the field <b>voltage</b> (current-driven).</div>"
        "<table>"
        "<tr><th>mode</th><th>reconstruction</th><th>solves/win</th><th>note</th></tr>"
        "<tr><td>pointwise (secant)</td><td>value at each point, accumulated window secant</td>"
        "<td>N_field_eval</td><td>curve-following; matches the Xyce Bfield secant exactly</td></tr>"
        "<tr><td>linear ramp</td><td>straight line carried-start &rarr; window-end</td>"
        "<td><b>1</b></td><td>cheapest; pure coupling reconstruction</td></tr>"
        "<tr><td>average</td><td>line, start = &frac12;(carried + end)</td>"
        "<td><b>1</b></td><td>linear + window-start damping (colleague)</td></tr>"
        "<tr><td>pointwise (central)</td><td>local finite diff: central (current-driven V) / "
        "BDF1 (voltage-driven I)</td>"
        "<td>N_field_eval</td><td>lowest raw RMS, but the derivative is a dummy artifact</td></tr>"
        "</table>"
        "<div class='hn'>Both coupling directions. All modes carry the seam (C0-continuous). Accuracy is "
        "within ~1&ndash;2% across modes; the extra solves buy little. Recommend <b>linear</b> / "
        "<b>average</b> (1 field solve per window).</div>"
    ),
    "interface_form": (
        "<div class='hh'>Interface stamping (Thevenin vs Norton)</div>"
        "<div class='hn'>How the field ROM's linearised V&ndash;I law (the matched secant, "
        "Z = Rrom + Lrom/dt) is put into the circuit. <b>Same fixpoint either way</b> (algebraic duals) &mdash; "
        "only the conditioning differs.<br>"
        "<b>Thevenin</b> (default): <code>Bfield</code> is a <b>voltage</b> source "
        "V(nx)=V(vfprev)+Z&middot;(I(Vmeas)&minus;V(iprev)). Pins V(p); the huge window-start Z multiplies "
        "only the small iteration change &rarr; robust.<br>"
        "<b>Norton</b>: the dual, a behavioural <b>current</b> source "
        "I=V(iprev)+(V(nx)&minus;V(vfprev))/Z with shunt G=1/Z. G&rarr;0 at the window start &rarr; V(p) "
        "weakly tied &rarr; stiffer. Same terminal (window-boundary) fixpoint, but the <b>interior V(p) "
        "waveform differs</b> from Thevenin (measured &gt; the signal amplitude on a 20&nbsp;kHz current "
        "source); it still converges (no dt-collapse seen up to ~MHz). Provided to compare the two.</div>"
    ),
    "coupling_mode": (
        "<div class='hh'>Coupling direction</div>"
        "<div class='hn'>voltage-driven (Dirichlet): circuit sets V(p), field returns current; "
        "matched-secant Bfield.<br>current-driven (Neumann): circuit sets I(Vmeas), field returns "
        "V_field; plain voltage source &mdash; removes the high-frequency window-start V(p) spike on "
        "current-source circuits.</div>"
    ),
    "circuit_kind": (
        "<div class='hh'>Circuit topology</div>"
        "<table>"
        "<tr><th>kind</th><th>circuit (coil always sits between port p and gnd)</th></tr>"
        "<tr><td>simple source</td><td>one source (sine V / sine I / step) + series R/L/C to the port</td></tr>"
        "<tr><td>#4 2-way (sine U,C)</td><td>sine U + cap C. drive [0,t1): U charges C &amp; drives the "
        "coil &rarr; freewheel [t1,&infin;): C &#8741; coil, U off (the open pre-t0 state isn't a throw)</td></tr>"
        "<tr><td>#5 2-way (DC U,C)</td><td>DC U + cap C. drive [0,t1) &rarr; freewheel [t1,&infin;)</td></tr>"
        "<tr><td>#6 2-way (AC vs R)</td><td>V_AC and R both at the port. AC-drive [0,t1): V_AC drives "
        "the coil &rarr; R-damp [t1,&infin;): R &#8741; coil, source off</td></tr>"
        "<tr><td>custom</td><td>free node-graph from the spec / drag-drop editor</td></tr>"
        "</table>"
    ),
    "switch_backend": (
        "<div class='hh'>Switch backend</div>"
        "<div class='hn'>behavioral R: each throw is a resistor R=Roff+(Ron&minus;Roff)&middot;g(t), the "
        "gate g a clamped trapezoid &mdash; no .MODEL, robust.<br>native S: Xyce voltage-controlled "
        "switch S + .MODEL VSWITCH, gated by a PWL control. Same schedule; can be stiffer at throws.</div>"
    ),
    "switch_t1": (
        "<div class='hh'>Switch throw time t1</div>"
        "<div class='hn'>First throw instant. #4/#5: drive &rarr; freewheel at t1. #6: AC-drive &rarr; "
        "R-damp at t1. Must fall inside the run (&lt; end time).</div>"
    ),
    "switch_C": (
        "<div class='hh'>Switch capacitor C (#4/#5)</div>"
        "<div class='hn'>Cap in the drive/freewheel branch (C &#8741; coil during freewheel). A large C "
        "stresses WR convergence (reactive coupling) &mdash; keep it small (presets: 1&micro;F).</div>"
    ),
    "switch_R": (
        "<div class='hh'>Switch resistor R (#6)</div>"
        "<div class='hn'>The damping resistor the switch places in parallel with the coil during the "
        "R-damp phase (presets: 10k&Omega;).</div>"
    ),
    "switch_Ron": (
        "<div class='hh'>Closed-switch resistance Ron</div>"
        "<div class='hn'><b>Gotcha:</b> the #4/#5 cap&harr;coil freewheel is a nearly-undamped LC loop; "
        "with a tiny Ron it rings and the timestep collapses at the freewheel throw. Use Ron &ge; a few "
        "&Omega; (presets: 10). #6 has no LC loop &rarr; Ron=1m&Omega; is fine.</div>"
    ),
    "switch_Roff": (
        "<div class='hh'>Open-switch resistance Roff</div>"
        "<div class='hn'>Resistance of an open throw (large, e.g. 1e9). Kept finite (not &infin;) so every "
        "node retains a well-defined admittance.</div>"
    ),
    "switch_trise": (
        "<div class='hh'>Switch transition time</div>"
        "<div class='hn'><b>Essential:</b> a finite ramp of each throw. An instantaneous throw "
        "disconnects an ideal branch carrying inductive coil current &rarr; voltage kick &rarr; "
        "dt-collapse. The gate ramps over this time at each edge.</div>"
    ),
}


def write_config(params):
    path = os.path.join(HERE, "sim_config.txt")
    with open(path, "w") as f:
        f.write("# generated by sim_ui.py\n")
        for k in DEFAULTS:
            v = params.get(k, DEFAULTS[k])
            if KINDS[k] in ("int", "choice"):
                f.write(f"{k} = {int(round(float(v)))}\n")
            else:
                f.write(f"{k} = {float(v):.10g}\n")
        # The circuit side is always authored via the text spec (circuit_spec.txt), so force the
        # custom node-graph path; source_kind / R_series / switch_* keep their C++ defaults (inert).
        f.write("circuit_kind = 4\n")
    return path


def read_config():
    """Load the existing sim_config.txt into a {key: value} dict (only keys the UI knows),
    so the studio opens on the actual saved configuration instead of the factory defaults.
    Unknown keys are ignored; missing file -> empty dict (fall back to DEFAULTS)."""
    path = os.path.join(HERE, "sim_config.txt")
    cfg = {}
    if not os.path.exists(path):
        return cfg
    with open(path) as f:
        for line in f:
            s = line.split("#", 1)[0].strip()
            if "=" not in s:
                continue
            k, v = (x.strip() for x in s.split("=", 1))
            if k not in DEFAULTS:
                continue
            try:
                cfg[k] = (int(round(float(v))) if KINDS[k] in ("int", "choice")
                          else float(v))
            except ValueError:
                pass
    return cfg


def ensure_built():
    main_path = os.path.join(HERE, "main")
    if not os.path.exists(main_path):
        subprocess.run(["make"], cwd=HERE, check=True,
                       capture_output=True, text=True)


def run_solver():
    ensure_built()
    t0 = time.perf_counter()
    proc = subprocess.run([os.path.join(HERE, "main")], cwd=HERE,
                          capture_output=True, text=True, timeout=600)
    elapsed = time.perf_counter() - t0
    return proc.returncode, proc.stdout, proc.stderr, elapsed


def emit_netlist():
    """Regenerate wr_circuit.cir from the current sim_config.txt without solving
    ('main emit'), so the rendered schematic reflects the current config."""
    ensure_built()
    subprocess.run([os.path.join(HERE, "main"), "emit"], cwd=HERE,
                   capture_output=True, text=True, timeout=60)


# Example custom node-graph spec (circuit_kind=4). Reserved nodes: p (port), 0 (ground).
DEFAULT_SPEC = """\
# Custom node-graph circuit. Reserved nodes: p = port (field attaches here), 0 = ground.
# <TYPE> <name> <nodeA> <nodeB> <params...>
#   R/L/C name a b val | {V,I}SIN name a b amp freq | {V,I}DC name a b val
#   {V,I}PULSE name a b v1 v2 td tr | {V,I}PWL name a b t1 v1 t2 v2 ...  (multi-step)
# Example: source+C  ||  R+L   (two parallel branches p->0)
VSIN src p a 1 50
C    c1  a 0 1e-6
R    r1  p b 10
L    l1  b 0 1.6e-7
"""


def read_spec():
    """Seed the spec editor from circuit_spec.txt if present, else the default example."""
    path = os.path.join(HERE, "circuit_spec.txt")
    if os.path.exists(path):
        with open(path) as f:
            return f.read()
    return DEFAULT_SPEC


def write_spec(params):
    """Write the custom node-graph spec to circuit_spec.txt (only when non-empty), so the C++
    generator reads it for circuit_kind=4. Ignored by the built-in kinds."""
    spec = params.get("circuit_spec")
    if spec:
        with open(os.path.join(HERE, "circuit_spec.txt"), "w") as f:
            f.write(spec)


# ---------------------------------------------------------------------------
# lcapy export (LaTeX / PDF)
# ---------------------------------------------------------------------------
def _read_sim_params():
    """Parse sim_params.inc (.PARAM name = value) into a {name: value-string} dict."""
    import re as _re
    d = {}
    path = os.path.join(HERE, "sim_params.inc")
    if os.path.exists(path):
        for ln in open(path):
            mm = _re.match(r"\s*\.PARAM\s+(\w+)\s*=\s*(\S+)", ln.split("*", 1)[0], _re.I)
            if mm:
                d[mm.group(1)] = mm.group(2)
    return d


def _resolve(tok, params):
    """{name} or a bare param name -> its sim_params value; otherwise the literal token."""
    import re as _re
    tok = tok.strip()
    mm = _re.fullmatch(r"\{(\w+)\}", tok)
    if mm:
        return params.get(mm.group(1), tok)
    return params.get(tok, tok)


def _split_netlist(raw):
    """Split the generated wr_circuit.cir into (head, tail) around the editable circuit side, using the
    '* === CIRCUIT SIDE' / '* === WR INTERFACE' markers WriteCircuitNetlist emits. head = title +
    .INCLUDE sim_params.inc (above the circuit); tail = the WR interface + directives (below). Both are
    shown read-only in the UI; the circuit side between them is edited as the text spec."""
    lines = raw.splitlines()
    ci = next((i for i, l in enumerate(lines) if l.startswith("* === CIRCUIT SIDE")), None)
    ii = next((i for i, l in enumerate(lines) if l.startswith("* === WR INTERFACE")), None)
    head = "\n".join(lines[:ci]).rstrip() if ci is not None else ""
    tail = "\n".join(lines[ii:]).strip() if ii is not None else ""
    return head, tail


def _lcapy_records(parsed, params):
    """Map the parsed Xyce netlist to lcapy device records for the schematic. Each record is
    {"line": "<lcapy device w/o hint>", "nodes": [a, b]}. R/L/C direct; sources -> V/I (sin/DC);
    behavioral gate resistors + native S -> switches. The WR interface is shown FAITHFULLY (not
    collapsed): Vmeas -> a 0 V source (ammeter) on p-nx, Bfield -> a labelled voltage source (the
    field-ROM Thevenin) on nx-0. Only the PWL signal carriers (VFprev/VIprev) and switch control
    sources are dropped -- they sit on isolated reference nodes and aren't physical branches."""
    def sine(expr):                                            # 'amp*sin(2*pi*freq*time)' -> (amp,freq)
        mm = re.search(r"([A-Za-z0-9_.+\-]+)\s*\*\s*sin\(\s*2\s*\*\s*pi\s*\*\s*([A-Za-z0-9_.+\-]+)\s*\*\s*time",
                       expr, re.I)
        return (_resolve(mm.group(1), params), _resolve(mm.group(2), params)) if mm else None

    def val(expr):                                            # R/L/C value: strip IC=, resolve {param}
        return _resolve(re.sub(r"\bIC\s*=\s*\S+", "", expr).strip(), params)

    def num(x):                                               # tidy number for a label (1.00e+00 -> 1)
        try:
            return f"{float(x):g}"
        except (TypeError, ValueError):
            return str(x)

    recs = []
    def emit(line, a, b, label=None):
        recs.append({"line": line, "nodes": [a, b], "label": label})
    for e in parsed["elements"]:
        nm, t, expr = e["name"], e["type"], e["expr"]
        if e["signal"] or nm.startswith("Vctrl"):
            continue                                          # drop PWL signal carriers + switch ctrl
        a, b = e["nodes"]
        u = nm[1:] or nm                                      # user name (strip the C++ type prefix)
        if nm == "Vmeas":                                     # 0 V ammeter on the interface branch
            emit(f"V{u} {a} {b} 0", a, b, "ammeter")
        elif nm == "Bfield":                                  # field-ROM Thevenin (matched secant)
            emit(f"Vfield {a} {b}", a, b, "field ROM")        # symbolic (the expr is not lcapy-drawable)
        elif t == "resistor" and expr.strip().startswith("R="):
            emit(f"SW{u} {a} {b}", a, b)                      # behavioral gate -> switch
        elif t == "resistor":
            emit(f"R{u} {a} {b} {val(expr)}", a, b)
        elif t == "inductor":
            emit(f"L{u} {a} {b} {val(expr)}", a, b)
        elif t == "capacitor":
            emit(f"C{u} {a} {b} {val(expr)}", a, b)
        elif nm[:1].upper() == "S":                           # native VC switch device
            emit(f"SW{u} {a} {b}", a, b)
        elif t in ("behavioral V", "voltage src"):
            s = sine(expr)
            # lcapy labels a sin() source with its DC offset (0), so give it an explicit amp/freq label.
            emit(f"V{u} {a} {b} sin(0 {s[0]} {s[1]})" if s
                 else f"V{u} {a} {b} {val(expr) if '{' in expr or expr.strip()[:1].isdigit() else 'dc 1'}",
                 a, b, f"{num(s[0])}V/{num(s[1])}Hz" if s else None)
        elif t in ("behavioral I", "current src"):
            s = sine(expr)
            emit(f"I{u} {a} {b} sin(0 {s[0]} {s[1]})" if s else f"I{u} {a} {b} dc 1",
                 a, b, f"{num(s[0])}A/{num(s[1])}Hz" if s else None)
    return recs


# --- ELK layout + circuitikz renderer -------------------------------------------------------------
# The schematic is laid out by ELK (Eclipse Layout Kernel via elk_layout.js/elkjs) -- a professional
# orthogonal placement + wire-routing engine -- and drawn with circuitikz (LaTeX), so the inline PNG
# and the LaTeX/PDF export are the SAME render for every topology. Nets are ELK nodes; each 2-terminal
# element is an ELK edge whose routed polyline becomes the wire, with the component symbol placed on it.
def _ck_component(rec):
    """(circuitikz to[] type, label) for a device record."""
    dev = rec["line"].split()[0]
    label = rec.get("label") or ""
    if dev == "Vmeas" or label == "ammeter":
        return "rmeter, t=A", ""                          # 0 V ammeter (A shown inside the meter)
    if dev == "Vfield" or label == "field ROM":
        return "sV", "field ROM"
    if dev.startswith("SW"):
        return "nos", dev[2:]                             # normally-open switch (label = switch name)
    if not label:
        toks = rec["line"].split()
        label = toks[3] if len(toks) > 3 else ""
    return {"V": "sV", "I": "sI", "R": "R", "L": "L", "C": "C"}.get(dev[0].upper(), "generic"), label


def _ck_sanitize(s):
    """Escape the few LaTeX-special characters that appear in element labels."""
    for a, b in (("µ", "\\textmu "), ("μ", "\\textmu "), ("Ω", "\\textohm "),
                 ("&", "\\&"), ("%", "\\%"), ("_", "\\_"), ("#", "\\#")):
        s = s.replace(a, b)
    return s


def _elk_layout(recs, gnd="0", port="p"):
    """Run the circuit graph (nets + 2-terminal element edges) through ELK; return {nets, routes, size}
    (net positions + the routed wire polyline per element), or None if node/elkjs is unavailable."""
    import json as _json
    import subprocess
    script = os.path.join(HERE, "elk_layout.js")
    if not os.path.exists(script) or not recs:
        return None
    comps = [{"id": f"c{i}", "a": r["nodes"][0], "b": r["nodes"][1]} for i, r in enumerate(recs)]
    nets = sorted({n for r in recs for n in r["nodes"]})
    payload = _json.dumps({"components": comps, "nets": nets, "port": port, "gnd": gnd})
    try:
        out = subprocess.run(["node", script], input=payload, capture_output=True, text=True,
                             cwd=HERE, timeout=20)
    except Exception:
        return None
    if out.returncode != 0 or not out.stdout.strip():
        return None
    try:
        d = _json.loads(out.stdout)
    except Exception:
        return None
    return d if d.get("nets") else None


def _circuitikz(recs, gnd="0", port="p"):
    """Emit a circuitikz picture from the ELK layout: each element's symbol on its routed wire, plus
    junction dots (nets with >=3 connections) and a ground symbol. Returns the tikz string, or None."""
    import collections
    import math
    lay = _elk_layout(recs, gnd, port)
    if lay is None:
        return None
    P, R, H, scale = lay["nets"], lay["routes"], lay["size"][1], 1.0 / 40.0

    def T(pt):                                            # ELK px (y-down) -> tikz cm (y-up)
        return (pt[0] * scale, (H - pt[1]) * scale)

    def ps(pt):
        return f"({pt[0]:.2f},{pt[1]:.2f})"

    body = []
    for i, r in enumerate(recs):
        route = R.get(f"c{i}") or [P[r["nodes"][0]], P[r["nodes"][1]]]
        pts = [T(q) for q in route]
        base, label = _ck_component(r)
        opt = base + (f", l={{{_ck_sanitize(label)}}}" if label else "")
        # Place the component on its LONGEST straight segment (most room for the symbol + label),
        # centred at a fixed drawn length with plain wires filling the rest of the route.
        seg = [math.dist(pts[j], pts[j + 1]) for j in range(len(pts) - 1)]
        k = max(range(len(seg)), key=lambda j: seg[j])
        a, b = pts[k], pts[k + 1]
        L = seg[k] or 1e-9
        ux, uy = (b[0] - a[0]) / L, (b[1] - a[1]) / L
        cl = min(L, 1.0)                                  # component drawn length (cm)
        mx, my = (a[0] + b[0]) / 2, (a[1] + b[1]) / 2
        c0 = (mx - ux * cl / 2, my - uy * cl / 2)
        c1 = (mx + ux * cl / 2, my + uy * cl / 2)
        pre = pts[:k + 1] + [c0]
        post = [c1] + pts[k + 1:]
        if len(pre) >= 2:
            body.append("\\draw " + " -- ".join(ps(p) for p in pre) + ";")
        body.append(f"\\draw {ps(c0)} to[{opt}] {ps(c1)};")
        if len(post) >= 2:
            body.append("\\draw " + " -- ".join(ps(p) for p in post) + ";")
    deg = collections.Counter(n for r in recs for n in r["nodes"])
    for n, xy in P.items():
        pt = ps(T(xy))
        if n == gnd:
            body.append(f"\\draw {pt} node[ground]{{}};")
        elif deg[n] >= 3:
            body.append(f"\\draw {pt} node[circ]{{}};")
    return "\\begin{circuitikz}\n" + "\n".join(body) + "\n\\end{circuitikz}"


def _compile_circuitikz(tikz):
    """Compile a circuitikz picture to (full_tex, pdf_bytes|None, png_bytes|None) via pdflatex (+ pdf->png).
    Degrades gracefully: returns the .tex even if pdflatex/converters are missing."""
    import tempfile
    import subprocess
    tex = ("\\documentclass[border=4pt]{standalone}\n\\usepackage{circuitikz}\n"
           "\\usepackage{textcomp}\n\\begin{document}\n" + tikz + "\n\\end{document}\n")
    d = tempfile.mkdtemp()
    with open(os.path.join(d, "c.tex"), "w") as f:
        f.write(tex)
    try:
        subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", "c.tex"],
                       cwd=d, capture_output=True, timeout=60)
    except Exception:
        return tex, None, None
    pdfp = os.path.join(d, "c.pdf")
    if not os.path.exists(pdfp):
        return tex, None, None
    pdf = open(pdfp, "rb").read()
    png, pngp = None, os.path.join(d, "c.png")
    for cmd in (["sips", "-s", "format", "png", pdfp, "--out", pngp],
                ["pdftoppm", "-png", "-r", "150", "-singlefile", pdfp, os.path.join(d, "c")]):
        try:
            subprocess.run(cmd, cwd=d, capture_output=True, timeout=30)
            if os.path.exists(pngp):
                png = open(pngp, "rb").read(); break
        except Exception:
            continue
    return tex, pdf, png


def lcapy_schematic_png(parsed, params):
    """Inline schematic as a PNG data-URI: ELK layout -> circuitikz -> pdflatex -> png (the same render
    as the export). Returns None if node/elkjs/pdflatex are unavailable."""
    import base64
    recs = _lcapy_records(parsed, params)
    tikz = _circuitikz(recs) if recs else None
    if tikz is None:
        return None
    _, _, png = _compile_circuitikz(tikz)
    return "data:image/png;base64," + base64.b64encode(png).decode() if png else None


def export_lcapy(params):
    """Regenerate the netlist from the live config and export the schematic (ELK layout + circuitikz):
    the standalone .tex + a compiled PDF (base64). Degrades gracefully when the layout/toolchain is
    unavailable."""
    import base64
    write_config(params); write_spec(params); emit_netlist()
    parsed = parse_netlist(os.path.join(HERE, "wr_circuit.cir"))
    recs = _lcapy_records(parsed, _read_sim_params())
    tikz = _circuitikz(recs) if recs else None
    if tikz is None:
        return {"ok": False, "error": "schematic layout unavailable -- needs node + elkjs (run 'npm install')."}
    tex, pdf, _png = _compile_circuitikz(tikz)
    return {"ok": True, "tex": tex,
            "pdf_b64": base64.b64encode(pdf).decode() if pdf else None,
            "warn": None if pdf else "PDF compile failed (pdflatex missing?) -- .tex still available."}


def graph_to_spec(graph):
    """Serialize an editor graph (from netlist_to_editor) into the custom circuit_spec line format so a
    preset can seed / be forked into a custom node-graph. Switch elements can't be expressed as a spec;
    they are skipped and flagged. Returns (spec_text, has_switch)."""
    lines, has_sw = [], False
    for c in graph:
        et, nm, (a, b), sub, p = c["et"], c["name"], c["nodes"], c.get("sub"), c.get("p", {})
        if et == "SW":
            has_sw = True
            continue
        if et in ("R", "L", "C"):
            lines.append(f"{et} {nm} {a} {b} {p.get('val', 0)}")
        elif et in ("V", "I"):
            if sub == "DC":
                lines.append(f"{et}DC {nm} {a} {b} {p.get('val', 1)}")
            elif sub == "PULSE":
                lines.append(f"{et}PULSE {nm} {a} {b} {p.get('v1',0)} {p.get('v2',1)} {p.get('td',0)} {p.get('tr',1e-4)}")
            elif sub == "PWL":
                lines.append(f"{et}PWL {nm} {a} {b} {p.get('pts','0 0')}")
            else:  # SIN
                lines.append(f"{et}SIN {nm} {a} {b} {p.get('amp',1)} {p.get('freq',50)}")
    return "\n".join(lines), has_sw


# ---------------------------------------------------------------------------
# Parsing .prn outputs
# ---------------------------------------------------------------------------
def _read_columns(path, ncols):
    """Read a whitespace .prn file (1 header line, optional 'End' trailer)."""
    t, cols = [], [[] for _ in range(ncols - 1)]
    if not os.path.exists(path):
        return np.array(t), [np.array(c) for c in cols]
    with open(path) as f:
        next(f, None)  # header
        for line in f:
            s = line.strip()
            if not s or s.startswith("End"):
                continue
            parts = s.split()
            if len(parts) < ncols:
                continue
            try:
                vals = [float(p) for p in parts[:ncols]]
            except ValueError:
                continue
            t.append(vals[1])              # column 1 = TIME (col 0 = index)
            for i in range(ncols - 1):
                cols[i].append(vals[i + 1])
    # cols[0] is TIME duplicate of t; keep cols[1:] as data
    return np.array(t), [np.array(c) for c in cols[1:]]


def read_outputs():
    out = {}
    # Field_solution.prn: Index TIME V(FIELD) I(FIELD)  -> synchronisation endpoints
    t, c = _read_columns(os.path.join(HERE, "Field_solution.prn"), 4)
    out["field"] = {"t": t, "V": c[0] if c else np.array([]),
                    "I": c[1] if len(c) > 1 else np.array([])}
    # Field_waveform_solution.prn: Index TIME V(FIELD) I(FIELD) -> reconstructed waveforms
    t, c = _read_columns(os.path.join(HERE, "Field_waveform_solution.prn"), 4)
    out["field_wave"] = {"t": t, "V": c[0] if c else np.array([]),
                         "I": c[1] if len(c) > 1 else np.array([])}
    # Circuit_solution.prn: Index TIME V(P) V(NX) I(VMEAS)
    t, c = _read_columns(os.path.join(HERE, "Circuit_solution.prn"), 5)
    out["circuit"] = {"t": t,
                      "Vp":  c[0] if len(c) > 0 else np.array([]),
                      "Vnx": c[1] if len(c) > 1 else np.array([]),
                      "I":   c[2] if len(c) > 2 else np.array([])}
    # WR_error.txt: comma separated "Time, RelErr, N_iter, Converged"
    out["wr"] = read_wr_error(os.path.join(HERE, "WR_error.txt"))
    return out


def read_wr_error(path):
    t, err, nit, conv = [], [], [], []
    if not os.path.exists(path):
        return {"t": np.array([]), "err": np.array([]),
                "nit": np.array([]), "conv": np.array([])}
    with open(path) as f:
        next(f, None)
        for line in f:
            parts = [p.strip() for p in line.split(",")]
            if len(parts) < 4:
                continue
            try:
                t.append(float(parts[0])); err.append(float(parts[1]))
                nit.append(float(parts[2])); conv.append(float(parts[3]))
            except ValueError:
                continue
    return {"t": np.array(t), "err": np.array(err),
            "nit": np.array(nit), "conv": np.array(conv)}


# ---------------------------------------------------------------------------
# Plot rendering -> base64 PNG
# ---------------------------------------------------------------------------
def _png(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=110, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    plt.close(fig)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


def make_plots(data):
    plots = {}
    c, fw, fld, wr = data["circuit"], data["field_wave"], data["field"], data["wr"]

    # 1) Voltage vs time
    fig, ax = plt.subplots(figsize=(8, 3.2))
    if c["t"].size:
        ax.plot(c["t"], c["Vp"], lw=0.9, label="V(p) circuit (port)")
    if fw["t"].size:
        ax.plot(fw["t"], fw["V"], "--", lw=1.1, label="V field (reconstructed)")
    if fld["t"].size:
        ax.plot(fld["t"], fld["V"], "o", ms=4, label="V field (window end)")
    ax.set_xlabel("time (s)"); ax.set_ylabel("voltage (V)")
    ax.set_title("Port / field voltage"); ax.grid(True, alpha=.3); ax.legend(fontsize=8)
    plots["voltage"] = _png(fig)

    # 2) Current vs time
    fig, ax = plt.subplots(figsize=(8, 3.2))
    if c["t"].size:
        ax.plot(c["t"], c["I"], lw=0.9, label="I(Vmeas) circuit")
    if fw["t"].size:
        ax.plot(fw["t"], fw["I"], "--", lw=1.1, label="I field (reconstructed)")
    if fld["t"].size:
        ax.plot(fld["t"], fld["I"], "o", ms=4, label="I field (window end)")
    ax.set_xlabel("time (s)"); ax.set_ylabel("current (A)")
    ax.set_title("Interface current (circuit vs field)")
    ax.grid(True, alpha=.3); ax.legend(fontsize=8)
    plots["current"] = _png(fig)

    # 3) WR convergence
    fig, ax = plt.subplots(figsize=(8, 3.2))
    if wr["t"].size:
        ax.semilogy(wr["t"], np.maximum(wr["err"], 1e-16), "o-", ms=4,
                    color="tab:red", label="final WR rel. error")
        ax.set_xlabel("window end time (s)")
        ax.set_ylabel("WR rel. error", color="tab:red")
        ax.tick_params(axis="y", labelcolor="tab:red")
        ax.grid(True, alpha=.3)
        ax2 = ax.twinx()
        ax2.plot(wr["t"], wr["nit"], "s--", ms=4, color="tab:blue",
                 label="WR iterations")
        ax2.set_ylabel("WR iterations", color="tab:blue")
        ax2.tick_params(axis="y", labelcolor="tab:blue")
    ax.set_title("Waveform-relaxation convergence per window")
    plots["wr"] = _png(fig)
    return plots


def scalar_summary(data):
    fld, wr = data["field"], data["wr"]
    s = {}
    if fld["t"].size:
        s["final_time_s"] = float(fld["t"][-1])
        s["final_V_field"] = float(fld["V"][-1]) if fld["V"].size else None
        s["final_I_field"] = float(fld["I"][-1]) if fld["I"].size else None
    if wr["t"].size:
        s["windows"] = int(wr["t"].size)
        s["max_WR_iterations"] = int(np.max(wr["nit"]))
        s["mean_WR_iterations"] = float(np.mean(wr["nit"]))
        s["total_xyce_solves"] = int(np.sum(wr["nit"]))  # cost proxy: one Xyce run per WR iter
        s["all_converged"] = bool(np.all(wr["conv"] >= 1.0))
        s["worst_WR_error"] = float(np.max(wr["err"]))
    return s


# ---------------------------------------------------------------------------
# Parameter sweep / convergence study
# ---------------------------------------------------------------------------
def sweep_values(lo, hi, steps, scale, kind):
    """Return the list of swept parameter values."""
    steps = max(2, int(steps))
    if scale == "log":
        if lo <= 0 or hi <= 0:
            raise ValueError("log sweep needs strictly positive min/max")
        vals = np.exp(np.linspace(np.log(lo), np.log(hi), steps))
    else:
        vals = np.linspace(lo, hi, steps)
    if kind == "int":
        vals = np.unique(np.round(vals).astype(int)).astype(float)
    return [float(v) for v in vals]


def run_sweep(base_params, key, values):
    """Run the solver once per swept value; collect per-point metrics."""
    rows = []
    for v in values:
        params = dict(base_params)
        params[key] = v
        write_config(params)
        rc, stdout, stderr, elapsed = run_solver()
        row = {"value": v, "solver_seconds": round(elapsed, 4)}
        if rc != 0:
            row["ok"] = False
            row["error"] = f"solver exited {rc}"
            row["log"] = "\n".join((stdout or stderr or "").splitlines()[-6:])
        else:
            row.update(scalar_summary(read_outputs()))
            row["ok"] = True
        rows.append(row)
    return rows


def make_sweep_plots(key, rows):
    """2x2 metric-vs-swept-parameter figure (convergence study)."""
    ok = [r for r in rows if r.get("ok")]
    if not ok:
        return {}
    x = np.array([r["value"] for r in ok], dtype=float)
    label = LABELS.get(key, key)

    def col(name):
        return np.array([r.get(name, np.nan) for r in ok], dtype=float)

    fig, axes = plt.subplots(2, 2, figsize=(10, 6.4))
    span = x.max() - x.min()
    logx = bool(x.min() > 0 and span > 0 and (x.max() / max(x.min(), 1e-300)) >= 50)

    def setx(ax):
        ax.set_xlabel(label)
        if logx:
            ax.set_xscale("log")
        ax.grid(True, alpha=.3)

    # (0,0) WR iterations (max + mean) -- convergence speed
    ax = axes[0, 0]
    ax.plot(x, col("max_WR_iterations"), "o-", color="tab:blue", label="max")
    ax.plot(x, col("mean_WR_iterations"), "s--", color="tab:cyan", label="mean")
    ax.set_ylabel("WR iterations / window"); ax.set_title("Convergence speed")
    ax.legend(fontsize=8); setx(ax)

    # (0,1) cost: total Xyce solves + solver time
    ax = axes[0, 1]
    ax.plot(x, col("total_xyce_solves"), "o-", color="tab:purple", label="Xyce solves")
    ax.set_ylabel("total Xyce solves", color="tab:purple")
    ax.tick_params(axis="y", labelcolor="tab:purple")
    a2 = ax.twinx()
    a2.plot(x, col("solver_seconds"), "^--", color="tab:green", label="solver s")
    a2.set_ylabel("solver time (s)", color="tab:green")
    a2.tick_params(axis="y", labelcolor="tab:green")
    ax.set_title("Cost"); setx(ax)

    # (1,0) worst WR error (log y)
    ax = axes[1, 0]
    ax.semilogy(x, np.maximum(col("worst_WR_error"), 1e-16), "o-", color="tab:red")
    ax.set_ylabel("worst WR rel. error"); ax.set_title("WR accuracy"); setx(ax)

    # (1,1) final interface values
    ax = axes[1, 1]
    ax.plot(x, col("final_I_field"), "o-", color="tab:orange", label="I_field")
    ax.set_ylabel("final I_field (A)", color="tab:orange")
    ax.tick_params(axis="y", labelcolor="tab:orange")
    a2 = ax.twinx()
    a2.plot(x, col("final_V_field"), "s--", color="tab:gray", label="V_field")
    a2.set_ylabel("final V_field (V)", color="tab:gray")
    a2.tick_params(axis="y", labelcolor="tab:gray")
    ax.set_title("Final interface values"); setx(ax)

    fig.suptitle(f"Convergence study: sweep of {label}", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    return {"sweep": _png(fig)}


def sweep_table(key, rows):
    """Compact serialisable table for the UI / CLI."""
    cols = ["value", "ok", "max_WR_iterations", "mean_WR_iterations",
            "total_xyce_solves", "worst_WR_error", "all_converged",
            "solver_seconds", "final_I_field", "final_V_field"]
    return {"key": key, "columns": cols,
            "rows": [{c: r.get(c) for c in cols} for r in rows]}


# ---------------------------------------------------------------------------
# Circuit netlist visualizer
# ---------------------------------------------------------------------------
ELEM_TYPE = {
    "R": ("resistor", "tab:gray"),
    "L": ("inductor", "tab:olive"),
    "C": ("capacitor", "tab:brown"),
    "V": ("voltage src", "tab:red"),
    "I": ("current src", "tab:orange"),
    "B": ("behavioral", "tab:blue"),
    "E": ("VCVS", "tab:red"), "G": ("VCCS", "tab:orange"),
    "F": ("CCCS", "tab:orange"), "H": ("CCVS", "tab:red"),
    "D": ("diode", "tab:green"), "Q": ("BJT", "tab:purple"),
    "M": ("MOSFET", "tab:purple"),
}


def parse_netlist(path):
    """Parse a SPICE netlist into elements + nodes. Joins '+' continuations,
    skips comments (*) and directives (.) but keeps them as 'directives'."""
    raw = ""
    if os.path.exists(path):
        with open(path) as f:
            raw = f.read()

    # Join continuation lines (leading '+') onto the previous logical line.
    logical = []
    for line in raw.splitlines():
        s = line.rstrip()
        if not s:
            continue
        if s.lstrip().startswith("+") and logical:
            logical[-1] += " " + s.lstrip()[1:].strip()
        else:
            logical.append(s)

    # SPICE convention: the first line of a deck is always the title, never an
    # element -- skip it so it isn't parsed as a device.
    title = logical[0].strip() if logical else ""
    body_lines = logical[1:] if logical else []

    elements, directives, nodes = [], [], []
    seen = set()
    for ln in body_lines:
        st = ln.strip()
        if st.startswith("*"):
            continue                      # comment
        if st.startswith("."):
            directives.append(st)
            continue
        toks = st.split()
        if len(toks) < 3:
            continue                      # need name + 2 nodes minimum
        name = toks[0]
        etype, color = ELEM_TYPE.get(name[0].upper(), ("element", "tab:cyan"))
        n1, n2 = toks[1], toks[2]
        rest = " ".join(toks[3:])
        # Behavioral V= / I= flavour.
        body = st.upper()
        if name[0].upper() == "B":
            if "I" in body.split("=")[0].split()[-1:] or " I " in body or "I=" in body.replace(" ", ""):
                etype = "behavioral I"
            else:
                etype = "behavioral V"
        # Compact descriptor.
        desc = rest
        fileref = None
        if "PWL" in body and "FILE" in body:
            etype = "PWL source"
            color = "tab:green"
            import re as _re
            m = _re.search(r'"([^"]+)"', st)
            fileref = m.group(1) if m else None
            desc = f"PWL FILE {fileref}" if fileref else "PWL FILE"
        if len(desc) > 60:
            desc = desc[:57] + "..."
        # PWL-FILE sources are signal carriers (they sit on isolated reference
        # nodes that behavioral expressions read via V(node)), not conductive
        # branches in the main loop.
        is_signal = (etype == "PWL source")
        elements.append({"name": name, "type": etype, "color": color,
                         "nodes": [n1, n2], "desc": desc, "file": fileref,
                         "expr": rest, "signal": is_signal})
        for nd in (n1, n2):
            if nd not in seen:
                seen.add(nd); nodes.append(nd)
    return {"elements": elements, "directives": directives,
            "nodes": nodes, "raw": raw, "title": title}


# ---------------------------------------------------------------------------
# Netlist -> editor graph (increment B: the SVG editor renders every circuit)
# ---------------------------------------------------------------------------
_SIN_RE = re.compile(r"([\w.eE+\-]+)\s*\*\s*sin\s*\(\s*2\s*\*\s*pi\s*\*\s*([\w.eE+\-]+)\s*\*\s*time", re.I)


def _val_tok(expr, params):
    """First whitespace token of an R/L/C value expression, resolved against sim_params
    ({name} or a bare param name -> number). Trailing 'IC=...' etc. is dropped."""
    tok = (expr or "").strip().split()
    return _resolve(tok[0], params) if tok else "0"


def netlist_to_editor(parsed, params):
    """Translate a parsed netlist into a component list the SVG editor (CE.load) can render
    for ANY circuit. Drops the WR coupling infrastructure that WriteCircuitNetlist always
    re-appends (interface Vmeas/Bfield, the vf/i PWL signal sources, native-switch control
    sources), so the graph is exactly the user-editable circuit side. Each entry:
      {et:'V'|'I'|'R'|'L'|'C'|'SW', name, sub, nodes:[a,b], p:{...}, editable:bool}
    Switches (native S or a behavioral gate resistor R={...}) map to 'SW' -- display only,
    because the custom-spec format supports only R/L/C + sources (see emitCustomTopology)."""
    out = []
    for e in parsed["elements"]:
        name, t, nodes, expr = e["name"], e["type"], e["nodes"], e.get("expr", "")
        # --- drop coupling infrastructure (re-appended by the generator) ---
        if name in ("Vmeas", "Bfield"):
            continue
        if e.get("signal") or t == "PWL source":          # VFprev / VIprev
            continue
        if name.startswith("Vctrl"):                       # native-switch control PWL
            continue
        # --- switches ---
        if t == "element" and "SWMOD" in expr.upper():     # native Xyce S: 'ctrl 0 SWMOD'
            out.append({"et": "SW", "name": name, "sub": None, "nodes": nodes,
                        "p": {}, "editable": False}); continue
        if t.startswith("resistor") and re.match(r"\s*R\s*=", expr, re.I):  # behavioral gate resistor R={...}
            out.append({"et": "SW", "name": name, "sub": None, "nodes": nodes,
                        "p": {}, "editable": False}); continue
        # --- passives ---
        if t.startswith("resistor"):
            out.append({"et": "R", "name": name, "sub": None, "nodes": nodes,
                        "p": {"val": _val_tok(expr, params)}, "editable": True}); continue
        if t.startswith("inductor"):
            out.append({"et": "L", "name": name, "sub": None, "nodes": nodes,
                        "p": {"val": _val_tok(expr, params)}, "editable": True}); continue
        if t.startswith("capacitor"):
            out.append({"et": "C", "name": name, "sub": None, "nodes": nodes,
                        "p": {"val": _val_tok(expr, params)}, "editable": True}); continue
        # --- sources ---
        if t in ("behavioral V", "behavioral I", "voltage src", "current src"):
            et = "I" if t in ("behavioral I", "current src") else "V"
            body = expr
            if "=" in body and body.split("=")[0].strip().upper() in ("V", "I"):
                body = body.split("=", 1)[1]               # 'V = { ... }' -> ' { ... }'
            body = body.strip().strip("{}").strip()
            m = _SIN_RE.search(body)
            if m:
                sub, p = "SIN", {"amp": _resolve(m.group(1), params),
                                 "freq": _resolve(m.group(2), params)}
            else:
                up = body.upper()
                if up.startswith("PULSE"):
                    seg = body[body.find("(")+1:] if "(" in body else body[5:]
                    nums = re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", seg)
                    nums = (nums + ["0", "1", "0", "1e-4"])[:4]
                    sub, p = "PULSE", {"v1": nums[0], "v2": nums[1], "td": nums[2], "tr": nums[3]}
                elif up.startswith("PWL"):
                    inner = body[body.find("(")+1:body.rfind(")")] if "(" in body else body[3:]
                    sub, p = "PWL", {"pts": " ".join(inner.replace(",", " ").split())}
                else:
                    tok = body.replace("DC", "").split()
                    sub, p = "DC", {"val": _resolve(tok[0], params) if tok else "0"}
            out.append({"et": et, "name": name, "sub": sub, "nodes": nodes,
                        "p": p, "editable": True}); continue
        # --- anything else: show as a generic (non-editable) box ---
        out.append({"et": "SW", "name": name, "sub": None, "nodes": nodes,
                    "p": {}, "editable": False})
    return out


# ---------------------------------------------------------------------------
# HTTP handler
# ---------------------------------------------------------------------------
class Handler(BaseHTTPRequestHandler):
    def log_message(self, *a):  # quieter console
        pass

    def _send(self, code, body, ctype="application/json"):
        data = body.encode() if isinstance(body, str) else body
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def do_GET(self):
        if self.path in ("/", "/index.html"):
            self._send(200, INDEX_HTML, "text/html; charset=utf-8")
        elif self.path.startswith("/netlist"):
            self._handle_netlist()
        else:
            self._send(404, "not found", "text/plain")

    def _handle_netlist(self, params=None):
        try:
            # Regenerate wr_circuit.cir so the drawn schematic matches the current config.
            # POST carries the live form params (write them first); GET uses the saved config.
            if params is not None:
                write_config(params)
                write_spec(params)
            emit_netlist()
            parsed = parse_netlist(os.path.join(HERE, "wr_circuit.cir"))
            sim_params = _read_sim_params()
            head, tail = _split_netlist(parsed["raw"])
            self._send(200, json.dumps({
                "ok": True,
                "schematic": lcapy_schematic_png(parsed, sim_params),   # inline schematic (data-URI or null)
                "locked_head": head,                                    # includes above the circuit side
                "locked_tail": tail,                                    # WR interface + directives (read-only)
                "elements": [{k: e[k] for k in ("name", "type", "nodes", "desc")}
                             for e in parsed["elements"]],
                "directives": parsed["directives"],
                "nodes": parsed["nodes"],
                "raw": parsed["raw"],
            }))
        except Exception as e:
            self._send(200, json.dumps({
                "ok": False, "error": str(e),
                "trace": traceback.format_exc()[-2000:],
            }))

    def do_POST(self):
        if self.path not in ("/run", "/sweep", "/netlist", "/export"):
            self._send(404, json.dumps({"error": "unknown endpoint"}))
            return
        n = int(self.headers.get("Content-Length", 0))
        body = json.loads(self.rfile.read(n) or b"{}")
        if self.path == "/sweep":
            self._handle_sweep(body)
            return
        if self.path == "/netlist":
            self._handle_netlist(body)
            return
        if self.path == "/export":
            try:
                self._send(200, json.dumps(export_lcapy(body)))
            except Exception as e:
                self._send(200, json.dumps({
                    "ok": False, "error": str(e),
                    "trace": traceback.format_exc()[-2000:],
                }))
            return
        try:
            params = body
            write_config(params)
            write_spec(params)
            rc, stdout, stderr, elapsed = run_solver()
            log_tail = "\n".join((stdout or "").splitlines()[-25:])
            if rc != 0:
                self._send(200, json.dumps({
                    "ok": False,
                    "error": f"solver exited with code {rc}",
                    "solver_seconds": elapsed,
                    "log": log_tail, "stderr": (stderr or "")[-2000:],
                }))
                return
            data = read_outputs()
            summary = scalar_summary(data)
            summary["solver_seconds"] = round(elapsed, 4)
            if summary.get("total_xyce_solves"):
                summary["sec_per_xyce_solve"] = round(elapsed / summary["total_xyce_solves"], 5)
            self._send(200, json.dumps({
                "ok": True,
                "plots": make_plots(data),
                "summary": summary,
                "solver_seconds": elapsed,
                "log": log_tail,
            }))
        except Exception as e:
            self._send(200, json.dumps({
                "ok": False, "error": str(e),
                "trace": traceback.format_exc()[-2000:],
            }))

    def _handle_sweep(self, body):
        try:
            base = body.get("params", {})
            key = body.get("sweep_key")
            if key not in SWEEPABLE:
                raise ValueError(f"'{key}' is not a sweepable parameter")
            values = sweep_values(float(body["min"]), float(body["max"]),
                                  int(body["steps"]), body.get("scale", "linear"),
                                  KINDS[key])
            t0 = time.perf_counter()
            rows = run_sweep(base, key, values)
            total = round(time.perf_counter() - t0, 3)
            n_ok = sum(1 for r in rows if r.get("ok"))
            self._send(200, json.dumps({
                "ok": n_ok > 0,
                "sweep_key": key,
                "n_points": len(rows),
                "n_ok": n_ok,
                "sweep_seconds": total,
                "plots": make_sweep_plots(key, rows),
                "table": sweep_table(key, rows),
            }))
        except Exception as e:
            self._send(200, json.dumps({
                "ok": False, "error": str(e),
                "trace": traceback.format_exc()[-2000:],
            }))


# ---------------------------------------------------------------------------
# Frontend
# ---------------------------------------------------------------------------
def _controls_html():
    # Seed the form from the actual saved sim_config.txt (fall back to factory defaults),
    # so the studio opens on the current working setup, not a blank/trivial config.
    initial = {**DEFAULTS, **read_config()}
    # (The circuit preset "library" lives in the circuit-netlist panel, not here.)
    rows = []
    for k, label, default, kind, slider in PARAMS:
        default = initial.get(k, default)
        help_icon = f'<span class="help" data-help="{k}">?</span>' if k in HELP else ''
        if kind == "choice":
            opts = "".join(
                f'<option value="{val}"{" selected" if val == default else ""}>{text}</option>'
                for val, text in slider.items()  # for "choice", 5th field is the options dict
            )
            rows.append(f"""
        <div class="ctl" id="ctl_{k}">
          <label for="f_{k}">{label}{help_icon}</label>
          <div class="inputs">
            <select id="f_{k}" data-key="{k}" class="choice">{opts}</select>
          </div>
        </div>""")
            continue
        step = "any" if kind == "float" else "1"
        slider_html = ""
        if slider:
            mn, mx, st = slider
            slider_html = (
                f'<input type="range" min="{mn}" max="{mx}" step="{st}" '
                f'value="{default}" data-key="{k}" class="slider" '
                f'oninput="syncFromSlider(this)">'
            )
        rows.append(f"""
        <div class="ctl" id="ctl_{k}">
          <label for="f_{k}">{label}{help_icon}</label>
          <div class="inputs">
            <input type="number" id="f_{k}" step="{step}" value="{default}"
                   data-key="{k}" oninput="syncFromBox(this)">
            {slider_html}
          </div>
        </div>""")
    # The circuit spec is edited in the circuit box (right panel), not here.
    return "\n".join(rows)


def _sweep_options_html():
    return "".join(
        f'<option value="{k}"{" selected" if k == "R_ROM" else ""}>{LABELS[k]}</option>'
        for k in SWEEPABLE
    )


INDEX_HTML = """<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>WR Co-Simulation Studio</title>
<style>
  :root { --bg:#0f1115; --panel:#181b22; --fg:#e6e6e6; --muted:#9aa4b2;
          --accent:#4f9dff; --ok:#39d98a; --err:#ff6b6b; }
  * { box-sizing: border-box; }
  body { margin:0; font:14px/1.45 -apple-system,Segoe UI,Roboto,sans-serif;
         background:var(--bg); color:var(--fg); }
  header { padding:14px 20px; background:var(--panel); border-bottom:1px solid #262b35;
           display:flex; align-items:center; gap:16px; }
  header h1 { font-size:16px; margin:0; font-weight:600; }
  header .sub { color:var(--muted); font-size:12px; }
  .wrap { display:flex; gap:16px; padding:16px; align-items:flex-start; }
  .panel { background:var(--panel); border:1px solid #262b35; border-radius:10px; padding:16px; }
  #controls { width:360px; flex:0 0 360px; position:sticky; top:16px; max-height:calc(100vh - 32px); overflow:auto; }
  #results { flex:1; min-width:0; }
  .ctl { margin-bottom:12px; }
  .ctl label { display:block; color:var(--muted); font-size:12px; margin-bottom:4px; }
  .inputs { display:flex; align-items:center; gap:8px; }
  .inputs input[type=number] { width:120px; background:#0d0f14; color:var(--fg);
       border:1px solid #2c333f; border-radius:6px; padding:6px 8px; font-family:ui-monospace,monospace; }
  .inputs select.choice { flex:1; background:#0d0f14; color:var(--fg);
       border:1px solid #2c333f; border-radius:6px; padding:6px 8px; }
  .slider { flex:1; }
  .btns { display:flex; gap:10px; margin-top:8px; }
  .mini { font-size:11px; color:var(--muted); align-self:center; font-family:ui-monospace,monospace; }
  button { background:var(--accent); color:#06122a; border:0; border-radius:8px;
           padding:10px 16px; font-weight:600; cursor:pointer; }
  button.secondary { background:#2c333f; color:var(--fg); }
  button:disabled { opacity:.5; cursor:default; }
  #status { margin:10px 0; font-size:13px; min-height:18px; }
  #status.ok { color:var(--ok); } #status.err { color:var(--err); }
  .plot { width:100%; border-radius:8px; background:#fff; margin-bottom:14px; }
  #p_circuit { background:#0f1115; border:1px solid #262b35; }
  .summary { display:grid; grid-template-columns:repeat(auto-fit,minmax(150px,1fr)); gap:10px; margin-bottom:14px; }
  .card { background:#0d0f14; border:1px solid #2c333f; border-radius:8px; padding:10px 12px; }
  .card .k { color:var(--muted); font-size:11px; } .card .v { font-size:16px; font-family:ui-monospace,monospace; }
  .card.cost { border-color:var(--accent); background:#10243f; } .card.cost .v { color:var(--accent); }
  pre#log { background:#0d0f14; border:1px solid #2c333f; border-radius:8px; padding:10px;
            color:var(--muted); font-size:11.5px; max-height:220px; overflow:auto; white-space:pre-wrap; }
  details summary { cursor:pointer; color:var(--muted); margin-bottom:8px; }
  .sweepbox { margin-top:16px; padding-top:14px; border-top:1px solid #2c333f; }
  .sweephd { font-weight:600; font-size:13px; margin-bottom:10px; color:var(--accent); }
  .sub2 { font-size:12px; min-height:16px; margin-top:6px; }
  .sub2.ok { color:var(--ok); } .sub2.err { color:var(--err); }
  .note { color:var(--muted); font-size:11px; margin-top:6px; }
  table.sweep { width:100%; border-collapse:collapse; font-size:11.5px;
                font-family:ui-monospace,monospace; margin-bottom:14px; }
  table.sweep th, table.sweep td { border:1px solid #2c333f; padding:4px 7px; text-align:right; }
  table.sweep th { color:var(--muted); font-weight:600; background:#0d0f14; }
  table.sweep tr.bad td { color:var(--err); }
  #circuitBox { margin-bottom:16px; padding-bottom:14px; border-bottom:1px solid #2c333f; }
  .secthd { display:flex; align-items:center; justify-content:space-between; margin-bottom:10px; }
  .secthd span { font-weight:600; font-size:13px; }
  button.small { padding:5px 10px; font-size:12px; }
  table.netlist { width:100%; border-collapse:collapse; font-size:12px;
                  font-family:ui-monospace,monospace; margin-bottom:6px; }
  table.netlist th, table.netlist td { border:1px solid #2c333f; padding:4px 8px; text-align:left; }
  table.netlist th { color:var(--muted); font-weight:600; background:#0d0f14; }
  table.netlist td.nm { color:var(--accent); }
  /* --- circuit preset library + spec editor + inline schematic view --- */
  .library { display:flex; align-items:center; gap:8px; flex-wrap:wrap; margin-bottom:12px;
             padding:8px 10px; background:#12151b; border:1px solid #2c333f; border-radius:8px; }
  .library .lib-lbl { color:var(--muted); font-size:12px; font-weight:600; }
  .library button { padding:5px 10px; font-size:12px; background:#22305a; color:#cfe0ff;
                    border:1px solid #33436e; border-radius:6px; cursor:pointer; }
  .library button:hover { filter:brightness(1.2); }
  .cedit { display:flex; gap:14px; margin-bottom:12px; align-items:flex-start; flex-wrap:wrap; }
  .cedit-l, .cedit-r { flex:1; min-width:280px; }
  .cedit-hd { display:flex; align-items:center; gap:8px; color:var(--muted); font-size:12px;
              font-weight:600; margin-bottom:6px; }
  .cedit-hd button { margin-left:auto; }
  #f_circuit_spec { width:100%; min-height:120px; background:#0d0f14; color:#e6e6e6;
       border:1px solid var(--accent); border-radius:0; padding:8px; resize:vertical;
       font-family:ui-monospace,monospace; font-size:12px; line-height:1.5; display:block; }
  pre.locked { margin:0; padding:6px 8px; background:#12151b; color:#7f8895;
       border:1px solid #2c333f; font-family:ui-monospace,monospace; font-size:11px; line-height:1.5;
       white-space:pre-wrap; overflow-x:auto; }
  #lockHead { border-radius:6px 6px 0 0; border-bottom:none; }
  #lockTail { border-radius:0 0 6px 6px; border-top:none; }
  .cedit-r img { width:100%; border-radius:8px; background:#fff; display:block; }
  .note.warn { color:var(--err); }
  /* --- hover help tooltip --- */
  .help { display:inline-block; margin-left:6px; width:14px; height:14px; border-radius:50%;
          background:#2c333f; color:var(--muted); font-size:10px; line-height:14px; text-align:center;
          cursor:help; user-select:none; }
  #helpTip { display:none; position:fixed; z-index:100; width:540px; max-width:72vw;
             background:#0d0f14; border:1px solid #3a4351; border-radius:8px; padding:10px 12px;
             box-shadow:0 8px 24px rgba(0,0,0,.55); color:var(--fg); font-size:11.5px; line-height:1.45; }
  #helpTip .hh { font-weight:600; margin-bottom:6px; color:var(--accent); }
  #helpTip .hn { color:var(--muted); margin-top:6px; }
  #helpTip table { border-collapse:collapse; width:100%; }
  #helpTip th, #helpTip td { border:1px solid #2c333f; padding:3px 7px; text-align:left; vertical-align:top; }
  #helpTip th { color:var(--muted); background:#12151b; font-weight:600; }
</style></head>
<body>
<header>
  <h1>WR Co-Simulation Studio</h1>
  <span class="sub">voltage-driven field model &middot; Xyce + dummy FEM &middot; waveform relaxation</span>
</header>
<div class="wrap">
  <div class="panel" id="controls">
    __CONTROLS__
    <div class="btns">
      <button id="runBtn" onclick="run()">Run simulation</button>
      <button class="secondary" onclick="resetDefaults()">Reset</button>
    </div>
    <div id="status"></div>

    <div class="sweepbox">
      <div class="sweephd">Convergence study (parameter sweep)</div>
      <div class="ctl">
        <label for="sw_key">Sweep parameter</label>
        <div class="inputs"><select id="sw_key" class="choice">__SWEEP_OPTS__</select></div>
      </div>
      <div class="ctl">
        <label>Range (min, max, steps)</label>
        <div class="inputs">
          <input type="number" id="sw_min" step="any" value="0" title="min">
          <input type="number" id="sw_max" step="any" value="0.1" title="max">
          <input type="number" id="sw_steps" step="1" value="8" title="steps" style="width:70px">
        </div>
      </div>
      <div class="ctl">
        <label for="sw_scale">Spacing</label>
        <div class="inputs">
          <select id="sw_scale" class="choice">
            <option value="linear">linear</option>
            <option value="log">log (positive only)</option>
          </select>
        </div>
      </div>
      <div class="btns">
        <button id="sweepBtn" onclick="runSweep()">Run sweep</button>
      </div>
      <div id="sweepStatus" class="sub2"></div>
      <div class="note">Holds all other fields at their current values and runs
        the solver once per swept value, plotting metrics vs the parameter.</div>
    </div>
  </div>
  <div class="panel" id="results">
    <div id="circuitBox">
      <div class="secthd">
        <span>Circuit netlist (<code>wr_circuit.cir</code>)</span>
        <div class="btns">
          <button class="secondary small" id="toggleCircuitBtn" onclick="toggleCircuit()">Hide</button>
          <button class="secondary small" onclick="loadCircuit()">Refresh</button>
          <button class="secondary small" onclick="exportCircuit()">Export LaTeX/PDF<span class="help" data-help="lcapy_export">?</span></button>
          <span class="mini" id="exportStatus"></span>
        </div>
      </div>
      <div id="circuitContent">
        <div class="library" id="library">
          <span class="lib-lbl">Library:</span>
          <span id="libButtons"></span>
        </div>
        <div class="cedit">
          <div class="cedit-l">
            <div class="cedit-hd">Circuit side<span class="help" data-help="circuit_spec_edit">?</span>
              <button class="small" onclick="applySpec()">Apply &amp; render</button></div>
            <pre class="locked" id="lockHead"></pre>
            <textarea id="f_circuit_spec" spellcheck="false">__SPEC_SEED__</textarea>
            <pre class="locked" id="lockTail"></pre>
            <div class="note">Only the <b>circuit side</b> (white box) is editable — the WR interface
              (<code>Vmeas</code>, <code>Bfield</code>) and directives are generated and locked. Reserved
              nodes <code>p</code>=port, <code>0</code>=gnd. One element/line: <code>R/L/C name a b val</code>,
              <code>{V,I}SIN name a b amp f</code>, <code>{V,I}DC name a b val</code>,
              <code>{V,I}PULSE name a b v1 v2 td tr</code>, <code>{V,I}PWL name a b t1 v1 …</code>.</div>
          </div>
          <div class="cedit-r">
            <div class="cedit-hd">Schematic <span class="mini">(lcapy)</span></div>
            <img id="p_circuit" alt="circuit schematic">
            <div class="note" id="schemNote"></div>
          </div>
        </div>
        <div id="netlistTable"></div>
        <details><summary>Full raw netlist + directives</summary><pre id="netlistRaw"></pre></details>
      </div>
    </div>
    <div class="summary" id="summary"></div>
    <img class="plot" id="p_voltage" hidden>
    <img class="plot" id="p_current" hidden>
    <img class="plot" id="p_wr" hidden>
    <img class="plot" id="p_sweep" hidden>
    <div id="sweepTable"></div>
    <details><summary>Solver log</summary><pre id="log"></pre></details>
  </div>
</div>
<div id="helpTip"></div>
<script>
const DEFAULTS = __DEFAULTS__;
const PRESETS = __PRESETS__;
const VISIBLE_WHEN = __VISIBILITY__;
const HELP = __HELP__;

// hover help: position:fixed tooltip so it escapes the controls panel's overflow clipping.
function helpShow(t){ const key=t.getAttribute('data-help'); if(!HELP[key])return;
  const tip=document.getElementById('helpTip'); tip.innerHTML=HELP[key]; tip.style.display='block';
  const r=t.getBoundingClientRect(); const w=tip.offsetWidth, h=tip.offsetHeight;
  let x=r.right+8; if(x+w>window.innerWidth-8) x=Math.max(8, r.left-w-8);
  let y=r.top;    if(y+h>window.innerHeight-8) y=Math.max(8, window.innerHeight-h-8);
  tip.style.left=x+'px'; tip.style.top=y+'px'; }
function helpHide(){ document.getElementById('helpTip').style.display='none'; }

// --- conditional visibility: hide options made irrelevant by another selection ---
function ctlVal(k){ const el=document.getElementById('f_'+k); return el?Math.round(parseFloat(el.value)):NaN; }
function condMatch(cond){ return cond.some(d => Object.keys(d).every(k => d[k].includes(ctlVal(k)))); }
function applyVisibility(){
  for(const k in VISIBLE_WHEN){ const box=document.getElementById('ctl_'+k);
    if(box) box.style.display = condMatch(VISIBLE_WHEN[k]) ? '' : 'none'; }
}

function buildLibrary(){
  const host=document.getElementById('libButtons'); if(!host) return;
  host.innerHTML='';
  for(const name in PRESETS){
    const b=document.createElement('button');
    b.textContent=name.split(':')[0];   // short chip: "P1", "P4", ...
    b.title=name;
    b.onclick=()=>applyPreset(name);
    host.appendChild(b);
  }
}
function applyPreset(name){
  const p = PRESETS[name];
  if (!p) return;
  for (const k in p){
    const b = document.getElementById('f_'+k);
    if (b){ b.value = p[k]; if (b.tagName !== 'SELECT') syncFromBox(b); }
  }
  setStatus('Loaded: '+name, '');
  applyVisibility();
  loadCircuit();
}

function syncFromSlider(el){
  const box = document.getElementById('f_'+el.dataset.key);
  if (box) box.value = el.value;
}
function syncFromBox(el){
  const sl = document.querySelector('.slider[data-key="'+el.dataset.key+'"]');
  if (sl) sl.value = el.value;
}
function collect(){
  const p = {};
  document.querySelectorAll('input[type=number][data-key]').forEach(b => {
    p[b.dataset.key] = parseFloat(b.value);
  });
  document.querySelectorAll('select[data-key]').forEach(s => {
    p[s.dataset.key] = parseFloat(s.value);
  });
  const ta = document.getElementById('f_circuit_spec');
  if (ta) p['circuit_spec'] = ta.value;   // custom node-graph spec (circuit_kind=custom)
  return p;
}
function resetDefaults(){
  for (const k in DEFAULTS){
    const b = document.getElementById('f_'+k);
    if (b){ b.value = DEFAULTS[k]; if (b.tagName !== 'SELECT') syncFromBox(b); }
  }
  setStatus('Reset to defaults.', '');
  applyVisibility();
}
function setStatus(msg, cls){
  const s = document.getElementById('status'); s.textContent = msg; s.className = cls;
}
function showSummary(sum){
  const el = document.getElementById('summary'); el.innerHTML = '';
  const fmt = v => (typeof v === 'number') ? (Math.abs(v)<1e-3||Math.abs(v)>=1e5 ? v.toExponential(4) : v.toPrecision(6)) : String(v);
  const order = ['solver_seconds','total_xyce_solves','sec_per_xyce_solve','windows','max_WR_iterations','worst_WR_error','all_converged','final_time_s','final_V_field','final_I_field'];
  const labels = {solver_seconds:'solver time (s)', total_xyce_solves:'Xyce solves', sec_per_xyce_solve:'s / Xyce solve', final_time_s:'final time (s)', final_V_field:'final V_field', final_I_field:'final I_field', max_WR_iterations:'max WR iters', worst_WR_error:'worst WR error', all_converged:'all converged'};
  for (const k of order){ if (k in sum){
    const c = document.createElement('div'); c.className='card';
    if (k === 'solver_seconds') c.classList.add('cost');
    c.innerHTML = '<div class="k">'+(labels[k]||k)+'</div><div class="v">'+fmt(sum[k])+'</div>';
    el.appendChild(c);
  }}
}
function setPlot(id, src){ const im = document.getElementById(id); if(src){ im.src=src; im.hidden=false; } else { im.hidden=true; } }

async function run(){
  const btn = document.getElementById('runBtn'); btn.disabled = true;
  setStatus('Running solver (Xyce + FEM + WR)...', '');
  try {
    const res = await fetch('/run', {method:'POST', headers:{'Content-Type':'application/json'},
                                     body: JSON.stringify(collect())});
    const j = await res.json();
    document.getElementById('log').textContent = (j.log||'') + (j.stderr? '\\n--- stderr ---\\n'+j.stderr : '') + (j.trace? '\\n'+j.trace : '');
    const tsec = (typeof j.solver_seconds === 'number') ? j.solver_seconds.toFixed(2)+' s' : '';
    if (!j.ok){ setStatus('Error: '+(j.error||'unknown')+(tsec?' (after '+tsec+')':''), 'err'); }
    else {
      setStatus('Done in '+tsec+'.', 'ok');
      showSummary(j.summary||{});
      setPlot('p_voltage', j.plots.voltage);
      setPlot('p_current', j.plots.current);
      setPlot('p_wr', j.plots.wr);
      loadCircuit();  // the solver regenerated the netlist; refresh the schematic
    }
  } catch(e){ setStatus('Request failed: '+e, 'err'); }
  finally { btn.disabled = false; }
}

function setSweepStatus(msg, cls){
  const s = document.getElementById('sweepStatus'); s.textContent = msg; s.className = 'sub2 '+(cls||'');
}
function fmtCell(v){
  if (v === null || v === undefined) return '-';
  if (typeof v === 'boolean') return v ? 'yes' : 'no';
  if (typeof v === 'number') return (Math.abs(v)>0 && (Math.abs(v)<1e-3||Math.abs(v)>=1e5)) ? v.toExponential(3) : (Number.isInteger(v)? v : v.toPrecision(5));
  return String(v);
}
function showSweepTable(tbl){
  const host = document.getElementById('sweepTable'); host.innerHTML = '';
  if (!tbl || !tbl.rows || !tbl.rows.length) return;
  const head = ['value','iters(max)','iters(mean)','Xyce','worstErr','conv','sec','I_field','V_field'];
  const cols = tbl.columns.filter(c => c !== 'ok');
  let html = '<table class="sweep"><thead><tr><th>'+head.join('</th><th>')+'</th></tr></thead><tbody>';
  for (const r of tbl.rows){
    const bad = (r.ok === false) || (r.all_converged === false);
    html += '<tr'+(bad?' class="bad"':'')+'>';
    for (const c of cols) html += '<td>'+fmtCell(r[c])+'</td>';
    html += '</tr>';
  }
  html += '</tbody></table>';
  host.innerHTML = html;
}
async function runSweep(){
  const btn = document.getElementById('sweepBtn'); btn.disabled = true;
  const key = document.getElementById('sw_key').value;
  setSweepStatus('Running sweep of '+key+'...', '');
  try {
    const payload = {
      params: collect(),
      sweep_key: key,
      min: parseFloat(document.getElementById('sw_min').value),
      max: parseFloat(document.getElementById('sw_max').value),
      steps: parseInt(document.getElementById('sw_steps').value, 10),
      scale: document.getElementById('sw_scale').value,
    };
    const res = await fetch('/sweep', {method:'POST', headers:{'Content-Type':'application/json'},
                                       body: JSON.stringify(payload)});
    const j = await res.json();
    if (!j.ok){ setSweepStatus('Sweep error: '+(j.error||'all points failed'), 'err'); }
    else {
      setSweepStatus('Sweep done: '+j.n_ok+'/'+j.n_points+' points in '+j.sweep_seconds.toFixed(1)+' s.', 'ok');
      setPlot('p_sweep', j.plots ? j.plots.sweep : null);
      showSweepTable(j.table);
    }
  } catch(e){ setSweepStatus('Request failed: '+e, 'err'); }
  finally { btn.disabled = false; }
}

function showNetlist(j){
  const host = document.getElementById('netlistTable'); host.innerHTML = '';
  if (j.elements && j.elements.length){
    let html = '<table class="netlist"><thead><tr><th>element</th><th>type</th><th>nodes</th><th>value / expression</th></tr></thead><tbody>';
    for (const e of j.elements){
      html += '<tr><td class="nm">'+e.name+'</td><td>'+e.type+'</td><td>'+e.nodes.join(' &harr; ')+'</td><td>'+(e.desc||'')+'</td></tr>';
    }
    html += '</tbody></table>';
    host.innerHTML = html;
  }
  document.getElementById('netlistRaw').textContent =
    (j.raw||'') + (j.directives && j.directives.length ? '\\n\\n--- directives ---\\n'+j.directives.join('\\n') : '');
}
function toggleCircuit(){
  const c = document.getElementById('circuitContent');
  const b = document.getElementById('toggleCircuitBtn');
  const hidden = c.style.display === 'none';
  c.style.display = hidden ? '' : 'none';
  b.textContent = hidden ? 'Hide' : 'Show';
}
async function loadCircuit(){
  try {
    // POST the live form so the drawn schematic reflects unsaved edits (server writes
    // config + regenerates wr_circuit.cir via 'main emit' before parsing).
    const res = await fetch('/netlist', {method:'POST', headers:{'Content-Type':'application/json'},
                                         body: JSON.stringify(collect())});
    const j = await res.json();
    if (j.ok){
      const img=document.getElementById('p_circuit'), sn=document.getElementById('schemNote');
      if (j.schematic){ img.src=j.schematic; img.style.display=''; sn.textContent=''; }
      else { img.removeAttribute('src'); img.style.display='none';
             sn.textContent='(no schematic — lcapy/pdflatex unavailable, or layout not supported for this circuit)'; }
      // Locked context around the editable circuit side (read-only). The textarea itself is authored
      // by the user / seeded by presets and is NOT overwritten here.
      document.getElementById('lockHead').textContent = j.locked_head || '';
      document.getElementById('lockTail').textContent = j.locked_tail || '';
      showNetlist(j);
    }
  } catch(e){ /* leave circuit box empty on failure */ }
}
function applySpec(){
  const ta=document.getElementById('f_circuit_spec');
  if (!ta || !ta.value.trim()){ setStatus('Circuit side is empty.','err'); return; }
  setStatus('Rendering circuit…','');
  loadCircuit();
}
function _dl(name, href){ const a=document.createElement('a'); a.href=href; a.download=name; document.body.appendChild(a); a.click(); a.remove(); }
async function exportCircuit(){
  const st=document.getElementById('exportStatus');
  st.textContent='exporting…'; st.style.color='var(--muted)';
  try {
    const res=await fetch('/export',{method:'POST',headers:{'Content-Type':'application/json'},
                                     body:JSON.stringify(collect())});
    const j=await res.json();
    if(!j.ok){ st.textContent='✗ '+(j.error||'export failed'); st.style.color='var(--err)'; return; }
    if(j.tex) _dl('circuit.tex','data:application/x-tex;charset=utf-8,'+encodeURIComponent(j.tex));
    if(j.pdf_b64) _dl('circuit.pdf','data:application/pdf;base64,'+j.pdf_b64);
    if(j.warn){ st.textContent='⚠ .tex only ('+j.warn+')'; st.style.color='var(--err)'; }
    else { st.textContent='✓ circuit.tex + circuit.pdf'; st.style.color='var(--ok)'; }
  } catch(e){ st.textContent='✗ '+e; st.style.color='var(--err)'; }
}

window.addEventListener('load', ()=>{
  const ctrls=document.getElementById('controls');
  if(ctrls){ ctrls.addEventListener('input', applyVisibility); ctrls.addEventListener('change', applyVisibility); }
  // Help hover is document-wide now (help icons live in both the controls and the circuit box).
  document.addEventListener('mouseover', e=>{ if(e.target.classList.contains('help')) helpShow(e.target); });
  document.addEventListener('mouseout',  e=>{ if(e.target.classList.contains('help')) helpHide(); });
  buildLibrary(); applyVisibility(); loadCircuit();
});
</script>
</body></html>
"""

import html as _html
INDEX_HTML = (INDEX_HTML
              .replace("__CONTROLS__", _controls_html())
              .replace("__SWEEP_OPTS__", _sweep_options_html())
              .replace("__SPEC_SEED__", _html.escape(read_spec()))
              .replace("__DEFAULTS__", json.dumps(DEFAULTS))
              .replace("__PRESETS__", json.dumps(PRESETS))
              .replace("__VISIBILITY__", json.dumps(VISIBLE_WHEN))
              .replace("__HELP__", json.dumps(HELP)))


def _png_to_file(data_uri, path):
    """Decode a 'data:image/png;base64,...' URI to a file."""
    b64 = data_uri.split(",", 1)[1]
    with open(path, "wb") as f:
        f.write(base64.b64decode(b64))


def cli_sweep(args):
    """Headless convergence study: run the sweep and print a table + save PNG."""
    if args.param not in SWEEPABLE:
        print(f"error: '{args.param}' not sweepable. Choose from: {', '.join(SWEEPABLE)}")
        return 2
    values = sweep_values(args.min, args.max, args.steps, args.scale, KINDS[args.param])
    print(f"Sweeping {args.param} over {len(values)} values "
          f"({args.scale}) from {values[0]:g} to {values[-1]:g}")
    rows = run_sweep(dict(DEFAULTS), args.param, values)

    hdr = ["value", "iters_max", "iters_mean", "Xyce", "worstErr", "conv", "sec", "I_field"]
    print("  ".join(f"{h:>11}" for h in hdr))
    for r in rows:
        if not r.get("ok"):
            print(f"{r['value']:>11g}  {'FAILED: ' + r.get('error',''):>60}")
            continue
        print("  ".join(f"{v:>11}" for v in [
            f"{r['value']:g}",
            r.get("max_WR_iterations", "-"),
            f"{r.get('mean_WR_iterations', float('nan')):.2f}",
            r.get("total_xyce_solves", "-"),
            f"{r.get('worst_WR_error', float('nan')):.2e}",
            "yes" if r.get("all_converged") else "no",
            f"{r.get('solver_seconds', float('nan')):.2f}",
            f"{r.get('final_I_field', float('nan')):.4g}",
        ]))
    plots = make_sweep_plots(args.param, rows)
    if plots.get("sweep"):
        _png_to_file(plots["sweep"], args.out)
        print(f"\nPlot written to {args.out}")
    return 0


def serve(args):
    srv = ThreadingHTTPServer((args.host, args.port), Handler)
    url = f"http://{args.host}:{args.port}"
    print(f"WR Co-Simulation Studio running at {url}")
    print("Press Ctrl+C to stop.")
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        print("\nShutting down.")
        srv.shutdown()


def main():
    ap = argparse.ArgumentParser(description="WR co-simulation UI + convergence study")
    sub = ap.add_subparsers(dest="cmd")

    ps = sub.add_parser("serve", help="run the web UI (default)")
    ps.add_argument("--port", type=int, default=8000)
    ps.add_argument("--host", default="127.0.0.1")

    pw = sub.add_parser("sweep", help="headless parameter sweep / convergence study")
    pw.add_argument("--param", required=True, help=f"one of: {', '.join(SWEEPABLE)}")
    pw.add_argument("--min", type=float, required=True)
    pw.add_argument("--max", type=float, required=True)
    pw.add_argument("--steps", type=int, default=8)
    pw.add_argument("--scale", choices=["linear", "log"], default="linear")
    pw.add_argument("--out", default="sweep.png", help="output plot file")

    # Back-compat: bare `--port` with no subcommand still serves.
    ap.add_argument("--port", type=int, default=8000, help=argparse.SUPPRESS)
    ap.add_argument("--host", default="127.0.0.1", help=argparse.SUPPRESS)
    args = ap.parse_args()

    if args.cmd == "sweep":
        sys.exit(cli_sweep(args))
    serve(args)


if __name__ == "__main__":
    main()
