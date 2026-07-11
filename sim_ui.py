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
import csv
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
    ("t_end",                           "Sim duration (s)",             2.0e-2,   "float", None),
    ("N_field_windows",                 "Field windows (total)",        50,       "int",   (1, 400, 1)),
    ("N_periods",                       "Number of source periods",     1,        "int",   (1, 10, 1)),
    ("N_field_steps_per_source_period", "Field steps / source period",  50,       "int",   (2, 200, 1)),
    ("N_field_eval_intervals",          "FEM eval intervals / window",  1,        "int",   (1, 64, 1)),
    ("N_xyce_coupling_intervals",       "Xyce coupling intervals",      100,      "int",   (2, 400, 1)),
    ("WRmaxSteps",                      "WR max iterations",            20,       "int",   (1, 100, 1)),
    ("WR_tolerance",                    "WR tolerance",                 1.0e-3,   "float", None),
    ("wr_convergence_method",           "WR convergence metric",        1,        "choice",
        {0: "waveform L1", 1: "terminal scalar"}),
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
# choice params: {key: {code: human-readable label}}, used to export labels not raw codes
CHOICES = {k: s for (k, _l, _d, kind, s) in PARAMS if kind == "choice"}
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
# The properties grid forces absolute-end-time mode (time_mode=1, injected in collect()); "Sim duration"
# = t_end. So the old time_mode-conditional rules are gone. What remains: N_field_eval_intervals only
# matters for the pointwise reconstruction modes (0=secant, 3=central diff) -- linear/average force 1
# solve/window and ignore it -- so it's shown only then.
VISIBLE_WHEN = {
    "N_field_eval_intervals": [{"reconstruct_mode": [0, 3]}],
}


# Hover help: key -> HTML shown in a tooltip next to the control's label (a "?" icon).
HELP = {
    "circuit_spec_edit": (
        "<div class='hh'>Circuit side (text)</div>"
        "<div class='hn'>The WR interface (<code>Vmeas</code> ammeter + <code>Bfield</code> field ROM) and "
        "the <code>.INCLUDE</code>/<code>.print</code>/<code>.end</code> directives are generated and shown "
        "locked around the editable box.</div>"
        "<div class='hn'><b>Reserved nodes:</b> <code>p</code> = port (field attaches here) &middot; "
        "<code>0</code> = ground.</div>"
        "<table>"
        "<tr><th>syntax (one element per line)</th><th>element</th><th>parameters</th></tr>"
        "<tr><td><code>R/L/C name a b value</code></td><td>resistor / inductor / capacitor</td>"
        "<td><code>name</code> label &middot; <code>a b</code> nodes &middot; <code>value</code> "
        "&Omega; (R) / H (L) / F (C)</td></tr>"
        "<tr><td><code>VSIN/ISIN name a b amp freq</code></td><td>sine source</td>"
        "<td><code>a b</code> nodes (+&nbsp;&rarr;&nbsp;&minus;) &middot; <code>amp</code> amplitude (V/A) "
        "&middot; <code>freq</code> frequency (Hz); waveform amp&middot;sin(2&pi;&middot;freq&middot;t)</td></tr>"
        "<tr><td><code>VDC/IDC name a b value</code></td><td>DC source</td>"
        "<td><code>a b</code> nodes &middot; <code>value</code> constant level (V/A)</td></tr>"
        "<tr><td><code>VPULSE/IPULSE name a b v1 v2 td tr</code></td><td>pulse source</td>"
        "<td><code>a b</code> nodes &middot; <code>v1</code> initial level &middot; <code>v2</code> pulsed "
        "level &middot; <code>td</code> delay (s) &middot; <code>tr</code> rise time (s)</td></tr>"
        "<tr><td><code>VPWL/IPWL name a b t1 v1 t2 v2 …</code></td><td>piecewise-linear source</td>"
        "<td><code>a b</code> nodes &middot; <code>t1 v1 t2 v2 …</code> (time&nbsp;s, level&nbsp;V/A) "
        "breakpoints, linearly interpolated between</td></tr>"
        "<tr><td><code>SW name a b tclose topen [Ron Roff trise]</code></td>"
        "<td>time-gated switch</td>"
        "<td><code>a b</code> nodes &middot; <code>tclose</code> close time (s) &middot; <code>topen</code> "
        "open time (s) &mdash; big topen (e.g. <code>1e30</code>) stays closed to the end &middot; "
        "<code>Ron</code> closed R (=10) &middot; <code>Roff</code> open R (=1e9) &middot; "
        "<code>trise</code> transition (=1e-5)</td></tr>"
        "</table>"
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
        "<div class='hn'><b>pointwise (secant)</b> with <code>FEM eval intervals / window = 1</code> "
        "collapses to <b>linear ramp</b>: one interval leaves only the carried start and the window end, "
        "so the accumulated secant is a single straight segment. Raise the eval intervals for it to actually "
        "follow the curve.</div>"
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
    "wr_convergence_method": (
        "<div class='hh'>WR convergence metric</div>"
        "<div class='hn'>What must drop below <code>WR_tolerance</code> to accept a window and stop "
        "iterating. <code>i</code> is the field-current waveform, k the WR iteration.</div>"
        "<div class='hn'><b>waveform L1</b> &mdash; relative L1 norm of the iteration-to-iteration change "
        "of the whole field-current waveform over the window (trapezoidal):"
        "<div class='hf'>&epsilon; = "
        "&int;|i<sup>(k)</sup>&minus;i<sup>(k&minus;1)</sup>|&nbsp;dt &nbsp;/&nbsp; "
        "&int;|i<sup>(k)</sup>|&nbsp;dt</div>"
        "Needs &ge;2 iterations (iteration&nbsp;1 returns sentinel 1.0). Whole-waveform, so it also catches "
        "interior mismatch, not only the endpoints.</div>"
        "<div class='hn'><b>terminal scalar</b> (reference port) &mdash; sum of four window-<i>end</i> "
        "terminal residuals: transmission mismatch (field vs circuit) + iteration change of V and I:"
        "<div class='hf'>&epsilon; = "
        "&Delta;<sub>rel</sub>(I<sub>field</sub>,I<sub>circ</sub>) + "
        "&Delta;<sub>rel</sub>(V<sub>field</sub>,V<sub>circ</sub>) + "
        "&Delta;<sub>rel</sub>(I<sub>field</sub><sup>(k)</sup>,I<sub>field</sub><sup>(k&minus;1)</sup>) + "
        "&Delta;<sub>rel</sub>(V<sub>field</sub><sup>(k)</sup>,V<sub>field</sub><sup>(k&minus;1)</sup>)</div>"
        "Each &Delta;<sub>rel</sub>(a,b)=|a&minus;b| made relative (&divide;|a|) when |a|&gt;0.1, else "
        "absolute. Converges from iteration&nbsp;1. Only the terminal scalars &mdash; cheaper, ignores the "
        "waveform interior.</div>"
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
    line = rec["line"]
    toks = line.split()
    is_sin = "sin(" in line                                # AC (sinusoidal) vs DC source
    if not label:
        if "dc" in toks:                                  # DC source -> label its value, not "dc"
            i = toks.index("dc")
            label = toks[i + 1] if i + 1 < len(toks) else ""
        else:
            label = toks[3] if len(toks) > 3 else ""
    # DC sources get the plain source symbol (V/I); sinusoidal sources the AC symbol (sV/sI).
    sym = {"V": "sV" if is_sin else "V", "I": "sI" if is_sin else "I",
           "R": "R", "L": "L", "C": "C"}.get(dev[0].upper(), "generic")
    return sym, label


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


# Hand-placed layouts for the built-in presets. ELK is a fine general fallback, but for the fixed
# preset topologies a tuned by-hand placement reads far cleaner (source left / field ROM right,
# orthogonal columns). Keyed by the preset's exact set of physical net names (which is stable per
# preset and shared by P1/P3 and P4/P5 -- same topology, different source), so an edited/custom
# circuit whose nets differ simply misses and falls back to ELK. Coordinates are in ELK's pixel
# convention (y DOWN, ~40 px = 1 cm); `routes` are per-element orthogonal polylines keyed by the
# element's unordered node pair (the component symbol lands on the polyline's longest segment).
_MANUAL_LAYOUTS = {
    # P1 / P3: source(left) - R - L - port; interface (ammeter over field ROM) on the right.
    frozenset({"s", "cm0", "p", "nx", "0"}): {
        "nets": {"s": (0, 0), "cm0": (120, 0), "p": (240, 0), "nx": (240, 120), "0": (120, 240)},
        "routes": {frozenset({"s", "0"}): [(0, 0), (0, 240), (120, 240)],
                   frozenset({"nx", "0"}): [(240, 120), (240, 240), (120, 240)]},
        "size": (240, 240)},
    # P2: bare current source(left) - port; interface on the right.
    frozenset({"p", "nx", "0"}): {
        "nets": {"p": (0, 0), "nx": (120, 0), "0": (60, 120)},
        "routes": {frozenset({"0", "p"}): [(0, 0), (0, 120), (60, 120)],
                   frozenset({"nx", "0"}): [(120, 0), (120, 120), (60, 120)]},
        "size": (120, 120)},
    # P4 / P5: SPDT wiper w (cap w->gnd) throws to source branch (a) / short branch (b), rejoining at
    # port p; interface (ammeter -> field ROM) on the far right. Switch + branches left, field right.
    frozenset({"a", "b", "w", "p", "nx", "0"}): {
        "nets": {"w": (40, 120), "a": (160, 40), "b": (160, 200), "p": (300, 120),
                 "nx": (420, 120), "0": (330, 240)},
        # The cap (outer branch) taps the wiper on a short stub (x 40->10) then drops at x=10, so its
        # drop doesn't run on top of fw's drop; drv (up) + fw (down) stay colinear on the x=40 wiper line.
        "routes": {frozenset({"w", "0"}): [(40, 120), (10, 120), (10, 240), (330, 240)],
                   frozenset({"w", "a"}): [(40, 120), (40, 40), (160, 40)],
                   frozenset({"w", "b"}): [(40, 120), (40, 200), (160, 200)],
                   frozenset({"a", "p"}): [(160, 40), (300, 40), (300, 120)],
                   frozenset({"b", "p"}): [(160, 200), (300, 200), (300, 120)],
                   frozenset({"nx", "0"}): [(420, 120), (420, 240), (330, 240)]},
        "size": (420, 240)},
    # P6: two parallel branches from port p to gnd -- (V + switch) and (R + switch) -- plus the
    # interface (ammeter + field ROM) as the third column. Top rail = p, bottom rail = gnd.
    # Verticals are 140 (> the 120 top-rail reach) so each column's device (V / R / ammeter) lands on
    # its vertical segment at the same height, not up on the top rail.
    frozenset({"p", "bac", "br", "nx", "0"}): {
        "nets": {"p": (160, 0), "bac": (40, 140), "br": (160, 140), "nx": (280, 140), "0": (160, 260)},
        "routes": {frozenset({"p", "bac"}): [(160, 0), (40, 0), (40, 140)],
                   frozenset({"bac", "0"}): [(40, 140), (40, 260), (160, 260)],
                   frozenset({"p", "nx"}): [(160, 0), (280, 0), (280, 140)],
                   frozenset({"nx", "0"}): [(280, 140), (280, 260), (160, 260)]},
        "size": (280, 260)},
}


def _manual_layout(recs):
    """Hand-placed layout for a known preset topology, matched by its physical net-name set; returns
    the same {nets, routes, size} shape as `_elk_layout`, or None (unknown/edited circuit -> ELK)."""
    present = frozenset(n for r in recs for n in r["nodes"])
    spec = _MANUAL_LAYOUTS.get(present)
    if spec is None:
        return None
    P = spec["nets"]
    if any(n not in P for r in recs for n in r["nodes"]):
        return None
    routes = {}
    for i, r in enumerate(recs):
        a, b = r["nodes"]
        poly = spec["routes"].get(frozenset((a, b)))
        if poly is None:
            poly = [P[a], P[b]]                            # axis-aligned pair -> straight wire
        elif list(poly[0]) != list(P[a]):
            poly = poly[::-1]                              # orient a->b (keeps source polarity sane)
        routes[f"c{i}"] = [list(pt) for pt in poly]
    return {"nets": {k: list(v) for k, v in P.items()}, "routes": routes, "size": list(spec["size"])}


def _circuitikz(recs, gnd="0", port="p"):
    """Emit a circuitikz picture from the layout (hand-placed for a known preset, else ELK): each
    element's symbol on its routed wire, plus junction dots (nets with >=3 connections) and a ground
    symbol. Returns the tikz string, or None."""
    import collections
    import math
    lay = _manual_layout(recs) or _elk_layout(recs, gnd, port)
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
    # High-res raster: pdftoppm at 300 dpi first (crisp); sips fallback upscales to ~2400px wide.
    for cmd in (["pdftoppm", "-png", "-r", "300", "-singlefile", pdfp, os.path.join(d, "c")],
                ["sips", "-s", "format", "png", "--resampleWidth", "2400", pdfp, "--out", pngp]):
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


def _save_data_uri_png(uri, dest):
    """Decode a ``data:image/png;base64,...`` URI to a PNG file. Return True on write."""
    if not uri or "base64," not in uri:
        return False
    b64 = uri.split("base64,", 1)[1]
    with open(dest, "wb") as f:
        f.write(base64.b64decode(b64))
    return True


def export_run_csv(path, params, summary, plots=None):
    """Append one row (all input params + all result metrics) to a CSV file.

    Params get a ``param_`` prefix, result metrics a ``result_`` prefix, so the
    two namespaces never collide. If the file already exists, the row is
    appended; if new columns appear across runs the whole file is rewritten with
    the unioned header so it stays valid.

    Any plots (name -> data-URI PNG) are written as ``<stem>_<stamp>_<name>.png``
    next to the CSV, and their filenames recorded in ``plot_<name>`` columns so
    each row points at its own images. Returns the absolute CSV path written.

    path    -- destination CSV path (created if missing).
    params  -- the flat parameter dict used for the run (the /run body).
    summary -- scalar_summary(...) plus solver_seconds etc.
    plots   -- optional {name: data-URI PNG} to dump alongside the CSV.
    """
    path = os.path.abspath(os.path.expanduser(path))
    stamp = time.strftime("%Y%m%d-%H%M%S")

    def flat(v):
        # keep CSV cells scalar; JSON-encode anything nested
        if isinstance(v, (dict, list, tuple)):
            return json.dumps(v)
        return v

    row = {"run_time": time.strftime("%Y-%m-%d %H:%M:%S")}
    for k, v in (params or {}).items():
        if k in CHOICES:  # export the human-readable label, not the raw 0/1/2 code
            try:
                v = CHOICES[k].get(int(round(float(v))), v)
            except (TypeError, ValueError):
                pass
        row["param_" + str(k)] = flat(v)
    for k, v in (summary or {}).items():
        row["result_" + str(k)] = flat(v)

    # dump plot PNGs into a <stem>_plots/ subdir beside the CSV, to avoid clutter;
    # the CSV records the relative path so links resolve from the CSV's location
    d = os.path.dirname(path)
    if d:
        os.makedirs(d, exist_ok=True)
    stem = os.path.splitext(os.path.basename(path))[0]
    plotdir = stem + "_plots"
    if plots:
        os.makedirs(os.path.join(d, plotdir) if d else plotdir, exist_ok=True)
    for name, uri in (plots or {}).items():
        fn = "{}_{}.png".format(stamp, name)
        rel = os.path.join(plotdir, fn)
        if _save_data_uri_png(uri, os.path.join(d, rel) if d else rel):
            row["plot_" + str(name)] = rel

    existing = []
    fieldnames = []
    if os.path.exists(path) and os.path.getsize(path) > 0:
        with open(path, newline="") as f:
            rd = csv.DictReader(f)
            fieldnames = list(rd.fieldnames or [])
            existing = list(rd)

    # union of old header + this row's keys, preserving old order then appending new
    for k in row:
        if k not in fieldnames:
            fieldnames.append(k)

    with open(path, "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        wr.writeheader()
        for r in existing:
            wr.writerow(r)
        wr.writerow(row)
    return path


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
        if self.path not in ("/run", "/sweep", "/netlist", "/export", "/export_csv"):
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
        if self.path == "/export_csv":
            try:
                out = export_run_csv(body.get("path") or "results.csv",
                                     body.get("params", {}), body.get("summary", {}),
                                     body.get("plots", {}))
                self._send(200, json.dumps({"ok": True, "path": out}))
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
# Properties layout: explicit rows (each an equal-column grid). time_mode is NOT shown -- the grid is
# absolute-end-time only ("Sim duration" = t_end) and collect() injects time_mode=1. Params not placed
# in a row but still consumed by the config (N_field_windows, N_periods, I_sat) are emitted as hidden
# inputs so presets/reset/collect keep working.
_PROP_ROWS = [
    ["coupling_mode", "reconstruct_mode", "N_field_eval_intervals",
     "N_field_steps_per_source_period", "N_xyce_coupling_intervals"],
    ["t_end", "wr_convergence_method", "WRmaxSteps", "WR_tolerance", "interface_form"],
    ["R_ROM", "L_ROM", "R_FEM", "L_FEM", "nonlin_model"],
]
_PROP_HIDDEN = ["N_field_windows", "N_periods", "I_sat"]  # config-only; time_mode injected in JS


def _controls_html():
    # Seed the form from the actual saved sim_config.txt (fall back to factory defaults),
    # so the studio opens on the current working setup, not a blank/trivial config.
    initial = {**DEFAULTS, **read_config()}
    spec = {k: (k, l, d, kind, s) for (k, l, d, kind, s) in PARAMS}

    def ctl(k):
        _k, label, default, kind, slider = spec[k]
        default = initial.get(k, default)
        help_icon = f'<span class="help" data-help="{k}">?</span>' if k in HELP else ''
        if kind == "choice":
            opts = "".join(
                f'<option value="{val}"{" selected" if val == default else ""}>{text}</option>'
                for val, text in slider.items()
            )
            inner = f'<select id="f_{k}" data-key="{k}" class="choice">{opts}</select>'
        else:
            step = "any" if kind == "float" else "1"
            inner = (f'<input type="number" id="f_{k}" step="{step}" value="{default}" '
                     f'data-key="{k}" oninput="syncFromBox(this)">')
        return (f'<div class="ctl" id="ctl_{k}"><label for="f_{k}">{label}{help_icon}</label>'
                f'<div class="inputs">{inner}</div></div>')

    out = []
    for row in _PROP_ROWS:
        cells = "\n".join(ctl(k) for k in row)
        out.append(f'<div class="prop-row" style="grid-template-columns:'
                   f'repeat({len(row)},minmax(0,1fr))">\n{cells}\n</div>')
    # config-only params: hidden number boxes (still collected + settable by presets/reset)
    hidden = "".join(
        f'<input type="number" id="f_{k}" value="{initial.get(k, spec[k][2])}" data-key="{k}" hidden>'
        for k in _PROP_HIDDEN
    )
    out.append(f'<div style="display:none">{hidden}</div>')
    return "\n".join(out)


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
  /* --- flat white theme; single accent (#3a0ca3) reserved for primary actions only --- */
  :root { --bg:#ffffff; --fg:#1a1a1a; --muted:#6f6f6f; --line:#e2e2e2;
          --line-strong:#cfcfcf; --field:#fafafa; --accent:#3a0ca3;
          --ok:#1a7f37; --err:#b3261e; }
  * { box-sizing:border-box; }
  body { margin:0; font:14px/1.5 -apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;
         background:var(--bg); color:var(--fg); -webkit-font-smoothing:antialiased; }
  a { color:var(--accent); }

  /* header: thin, quiet, one full-width rule under it */
  header { padding:22px 0; border-bottom:1px solid var(--line); }
  .head-in { max-width:1180px; margin:0 auto; padding:0 40px;
             display:flex; align-items:baseline; gap:16px; }
  header h1 { font-size:15px; margin:0; font-weight:600; letter-spacing:.02em; }
  header .sub { color:var(--muted); font-size:12px; letter-spacing:.01em; }

  main { max-width:1180px; margin:0 auto; padding:0 40px 80px; }

  /* sections: flat, separated only by a long thin rule */
  section { padding:34px 0; border-bottom:1px solid var(--line); }
  .sec-hd { display:flex; align-items:center; justify-content:space-between; gap:16px;
            margin-bottom:20px; }
  .sec-hd h2 { font-size:11px; font-weight:600; text-transform:uppercase; letter-spacing:.14em;
               color:var(--muted); margin:0; }
  .sec-hd .hd-tools { display:flex; align-items:center; gap:14px; }

  /* controls: label above input, flat field with thin border, no radius */
  .ctl label { display:block; color:var(--muted); font-size:11px; margin-bottom:5px;
               letter-spacing:.01em; }
  .inputs { display:flex; align-items:center; gap:8px; }
  .ctl.disabled { opacity:.4; }
  input:disabled, select:disabled { cursor:not-allowed; background:var(--line); color:var(--muted); }
  input[type=number], select.choice, input[type=text] {
       width:100%; background:var(--field); color:var(--fg);
       border:1px solid var(--line-strong); border-radius:0; padding:7px 9px;
       font-family:ui-monospace,SFMono-Regular,Menlo,monospace; font-size:12.5px; }
  select.choice { font-family:inherit; }
  input:focus, select:focus, textarea:focus { outline:none; border-color:var(--accent); }
  .hd-tools input.pathin { width:220px; padding:5px 8px; font-size:12px; }
  .hd-tools button, .hd-tools input.pathin { flex-shrink:0; }
  #csvStatus { max-width:320px; overflow:hidden; text-overflow:ellipsis; white-space:nowrap; }

  /* properties: explicit equal-column rows (grid cols set inline per row) */
  .prop-rows { display:flex; flex-direction:column; gap:18px; }
  .prop-row { display:grid; gap:16px 22px; align-items:start; }
  /* the convergence-study sweep block still uses a plain auto-fill grid */
  .properties { display:grid; grid-template-columns:repeat(auto-fill,minmax(190px,1fr));
                gap:16px 22px; }

  /* buttons: flat, square. accent only on the primary run button */
  .btns { display:flex; gap:12px; align-items:center; flex-wrap:wrap; }
  button { background:var(--bg); color:var(--fg); border:1px solid var(--line-strong);
           border-radius:0; padding:9px 18px; font-size:13px; font-weight:500; cursor:pointer;
           letter-spacing:.01em; transition:background .12s,border-color .12s; }
  button:hover { border-color:var(--fg); }
  button.primary { background:var(--accent); color:#fff; border-color:var(--accent); font-weight:600; }
  button.primary:hover { background:#2e0982; border-color:#2e0982; }
  button.small { padding:5px 11px; font-size:12px; }
  button:disabled { opacity:.45; cursor:default; }

  .mini { font-size:11px; color:var(--muted); font-family:ui-monospace,monospace; }
  .note { color:var(--muted); font-size:11px; margin-top:8px; line-height:1.5; }
  .note.warn { color:var(--err); }
  #status, .sub2 { font-size:12.5px; min-height:16px; }
  #status { margin-top:14px; }
  #status.ok, .sub2.ok { color:var(--ok); }
  #status.err, .sub2.err { color:var(--err); }

  /* native collapsibles for preset library + convergence study */
  details.fold { border:1px solid var(--line); margin-top:4px; }
  details.fold > summary { cursor:pointer; list-style:none; padding:10px 14px;
       font-size:11px; text-transform:uppercase; letter-spacing:.12em; color:var(--muted);
       background:var(--field); user-select:none; }
  details.fold > summary::-webkit-details-marker { display:none; }
  details.fold > summary::before { content:"+ "; color:var(--accent); font-weight:600; }
  details.fold[open] > summary::before { content:"– "; }
  details.fold > .fold-body { padding:16px 14px; }
  details.plain > summary { cursor:pointer; color:var(--muted); font-size:12px; margin-top:6px; }

  /* circuit section: input (left) + schematic (right) */
  .cedit { display:flex; gap:28px; align-items:flex-start; flex-wrap:wrap; }
  .cedit-l, .cedit-r { flex:1; min-width:300px; }
  /* right column is capped to the textarea height (set in JS); schematic fills the remainder */
  .cedit-r { display:flex; flex-direction:column; min-height:0; }
  #schemWrap { flex:1 1 auto; min-height:0; display:flex; }
  .cedit-hd { display:flex; align-items:center; gap:8px; color:var(--muted); font-size:11px;
              text-transform:uppercase; letter-spacing:.1em; margin-bottom:8px; }
  .cedit-hd button { margin-left:auto; }
  #f_circuit_spec { width:100%; min-height:240px; height:240px; background:var(--field); color:var(--fg);
       border:1px solid var(--line-strong); border-radius:0; padding:9px; resize:vertical;
       font-family:ui-monospace,monospace; font-size:12px; line-height:1.6; display:block; }
  pre.locked { margin:0; padding:7px 9px; background:#f3f3f3; color:#8a8a8a;
       border:1px solid var(--line); font-family:ui-monospace,monospace; font-size:11px;
       line-height:1.6; white-space:pre-wrap; overflow-x:auto; }
  #lockHead { border-bottom:none; }
  #lockTail { border-top:none; }
  .cedit-r img { width:100%; height:100%; object-fit:contain; object-position:center top;
                 border:1px solid var(--line); background:#fff; display:block; }

  /* preset library chips */
  .library { display:flex; align-items:center; gap:8px; flex-wrap:wrap; }
  .library button { padding:5px 12px; font-size:12px; }

  /* result summary cards: flat, thin border; cost card marked with the accent */
  .summary { display:grid; grid-template-columns:repeat(auto-fill,minmax(150px,1fr));
             gap:0; border:1px solid var(--line); border-bottom:none; margin-bottom:24px; }
  .card { border-bottom:1px solid var(--line); border-right:1px solid var(--line); padding:12px 14px; }
  .card .k { color:var(--muted); font-size:10.5px; text-transform:uppercase; letter-spacing:.08em; }
  .card .v { font-size:16px; font-family:ui-monospace,monospace; margin-top:3px; }
  .card.cost .v { color:var(--accent); font-weight:600; }

  /* plots: two per row, click to enlarge */
  .plots { display:grid; grid-template-columns:repeat(2,minmax(0,1fr)); gap:20px; }
  .plot { width:100%; background:#fff; border:1px solid var(--line); display:none;
          cursor:zoom-in; }
  .plot.wide { grid-column:1 / -1; }

  pre#log { background:var(--field); border:1px solid var(--line); border-radius:0; padding:12px;
            color:var(--muted); font-size:11.5px; max-height:240px; overflow:auto; white-space:pre-wrap; }

  table.sweep, table.netlist { width:100%; border-collapse:collapse; font-size:12px;
                font-family:ui-monospace,monospace; margin-bottom:12px; }
  table.sweep th, table.sweep td { border:1px solid var(--line); padding:5px 8px; text-align:right; }
  table.netlist th, table.netlist td { border:1px solid var(--line); padding:5px 9px; text-align:left; }
  table.sweep th, table.netlist th { color:var(--muted); font-weight:600; background:var(--field);
                font-size:11px; text-transform:uppercase; letter-spacing:.06em; }
  table.sweep tr.bad td { color:var(--err); }
  table.netlist td.nm { color:var(--accent); }

  /* lightbox for enlarged plots */
  #lightbox { display:none; position:fixed; inset:0; z-index:200; background:rgba(255,255,255,.94);
              align-items:center; justify-content:center; cursor:zoom-out; padding:40px; }
  #lightbox.on { display:flex; }
  #lightbox img { max-width:96vw; max-height:92vh; border:1px solid var(--line-strong); background:#fff; }

  /* hover help tooltip */
  .help { display:inline-block; margin-left:6px; width:14px; height:14px; border-radius:50%;
          background:#ececec; color:var(--muted); font-size:10px; line-height:14px; text-align:center;
          cursor:help; user-select:none; }
  #helpTip { display:none; position:fixed; z-index:210; width:540px; max-width:72vw;
             background:#fff; border:1px solid var(--line-strong); border-radius:0; padding:12px 14px;
             box-shadow:0 10px 30px rgba(0,0,0,.14); color:var(--fg); font-size:11.5px; line-height:1.5; }
  #helpTip .hh { font-weight:600; margin-bottom:6px; color:var(--accent); }
  #helpTip .hn { color:var(--muted); margin-top:6px; }
  #helpTip .hf { color:var(--fg); font-family:ui-monospace,SFMono-Regular,Menlo,monospace;
       font-size:12px; margin:5px 0; padding:5px 9px; background:var(--field);
       border-left:2px solid var(--accent); overflow-x:auto; }
  #helpTip table { border-collapse:collapse; width:100%; }
  #helpTip th, #helpTip td { border:1px solid var(--line); padding:4px 7px; text-align:left; vertical-align:top; }
  #helpTip th { color:var(--muted); background:var(--field); font-weight:600; }
</style></head>
<body>
<header>
  <div class="head-in">
    <h1>WR Co-Simulation Studio</h1>
    <span class="sub">voltage-driven field model &middot; Xyce + dummy FEM &middot; waveform relaxation</span>
  </div>
</header>
<main>

  <!-- 1. Circuit input + schematic -->
  <section id="circuitBox">
    <div class="sec-hd">
      <h2>Circuit side<span class="help" data-help="circuit_spec_edit">?</span></h2>
      <div class="hd-tools">
        <span class="mini" id="exportStatus"></span>
        <button class="small primary" onclick="applySpec()">Apply &amp; render</button>
        <button class="small" onclick="loadCircuit()">Refresh</button>
        <button class="small" onclick="exportCircuit()">Export LaTeX/PDF<span class="help" data-help="lcapy_export">?</span></button>
      </div>
    </div>
    <div class="cedit">
      <div class="cedit-l">
        <pre class="locked" id="lockHead"></pre>
        <textarea id="f_circuit_spec" spellcheck="false">__SPEC_SEED__</textarea>
        <pre class="locked" id="lockTail"></pre>
      </div>
      <div class="cedit-r">
        <details class="fold" id="libFold" style="margin-top:0">
          <summary>Preset library</summary>
          <div class="fold-body">
            <div class="library" id="library"><span id="libButtons"></span></div>
          </div>
        </details>
        <div class="cedit-hd" style="margin-top:8px">Schematic <span class="mini">(lcapy)</span></div>
        <div id="schemWrap"><img id="p_circuit" alt="circuit schematic"></div>
        <div class="note" id="schemNote"></div>
      </div>
    </div>
    <div id="netlistTable" style="margin-top:8px"></div>
    <details class="plain"><summary>Full raw netlist + directives</summary><pre id="netlistRaw"></pre></details>
  </section>

  <!-- 2. Properties -->
  <section>
    <div class="sec-hd"><h2>Properties</h2></div>
    <div class="prop-rows">
      __CONTROLS__
    </div>
  </section>

  <!-- 3. Actions -->
  <section>
    <div class="sec-hd"><h2>Run</h2></div>
    <div class="btns">
      <button id="runBtn" class="primary" onclick="run()">Run simulation</button>
      <button onclick="resetDefaults()">Reset</button>
    </div>
    <div id="status"></div>
    <details class="fold" style="margin-top:20px">
      <summary>Convergence study (parameter sweep)</summary>
      <div class="fold-body">
        <div class="properties" style="grid-template-columns:repeat(auto-fill,minmax(150px,1fr))">
          <div class="ctl">
            <label for="sw_key">Sweep parameter</label>
            <div class="inputs"><select id="sw_key" class="choice">__SWEEP_OPTS__</select></div>
          </div>
          <div class="ctl">
            <label for="sw_min">Min</label>
            <div class="inputs"><input type="number" id="sw_min" step="any" value="0"></div>
          </div>
          <div class="ctl">
            <label for="sw_max">Max</label>
            <div class="inputs"><input type="number" id="sw_max" step="any" value="0.1"></div>
          </div>
          <div class="ctl">
            <label for="sw_steps">Steps</label>
            <div class="inputs"><input type="number" id="sw_steps" step="1" value="8"></div>
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
        </div>
        <div class="btns" style="margin-top:16px"><button id="sweepBtn" onclick="runSweep()">Run sweep</button></div>
        <div id="sweepStatus" class="sub2" style="margin-top:10px"></div>
        <div class="note">Holds all other fields at their current values and runs
          the solver once per swept value, plotting metrics vs the parameter.</div>
      </div>
    </details>
  </section>

  <!-- 4. Results -->
  <section style="border-bottom:none">
    <div class="sec-hd">
      <h2>Results</h2>
      <div class="hd-tools">
        <span class="mini" id="csvStatus"></span>
        <input type="text" id="csvPath" class="pathin" placeholder="results.csv" value="results.csv" spellcheck="false">
        <button class="small" onclick="exportCsv()">Export CSV</button>
      </div>
    </div>
    <div class="summary" id="summary"></div>
    <div class="plots">
      <img class="plot" id="p_voltage" onclick="enlarge(this)">
      <img class="plot" id="p_current" onclick="enlarge(this)">
      <img class="plot wide" id="p_wr" onclick="enlarge(this)">
      <img class="plot wide" id="p_sweep" onclick="enlarge(this)">
    </div>
    <div id="sweepTable"></div>
    <details class="plain"><summary>Solver log</summary><pre id="log"></pre></details>
  </section>

</main>
<div id="lightbox" onclick="this.classList.remove('on')"><img id="lightboxImg"></div>
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
  for(const k in VISIBLE_WHEN){
    const on = condMatch(VISIBLE_WHEN[k]);
    const box = document.getElementById('ctl_'+k);
    if(box){
      box.classList.toggle('disabled', !on);   // dim, but stay in place (no reflow)
      box.querySelectorAll('input,select').forEach(el => el.disabled = !on);
    }
  }
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
  p['time_mode'] = 1;   // properties grid is absolute-end-time only ("Sim duration" = t_end)
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
  const fmt = v => (typeof v === 'number') ? (Number.isInteger(v) ? String(v) : (Math.abs(v)<1e-3||Math.abs(v)>=1e5 ? v.toExponential(4) : v.toPrecision(6))) : String(v);
  const order = ['solver_seconds','total_xyce_solves','sec_per_xyce_solve','windows','max_WR_iterations','worst_WR_error','all_converged','final_time_s','final_V_field','final_I_field'];
  const labels = {solver_seconds:'solver time (s)', total_xyce_solves:'Xyce solves', sec_per_xyce_solve:'s / Xyce solve', final_time_s:'final time (s)', final_V_field:'final V_field', final_I_field:'final I_field', max_WR_iterations:'max WR iters', worst_WR_error:'worst WR error', all_converged:'all converged'};
  for (const k of order){ if (k in sum){
    const c = document.createElement('div'); c.className='card';
    if (k === 'solver_seconds') c.classList.add('cost');
    c.innerHTML = '<div class="k">'+(labels[k]||k)+'</div><div class="v">'+fmt(sum[k])+'</div>';
    el.appendChild(c);
  }}
}
function setPlot(id, src){ const im = document.getElementById(id);
  if(src){ im.src=src; im.style.display='block'; } else { im.removeAttribute('src'); im.style.display='none'; } }
function enlarge(im){ if(!im.src) return;
  document.getElementById('lightboxImg').src = im.src;
  document.getElementById('lightbox').classList.add('on'); }

// Cap the right column (preset library + schematic) to the full circuit-side block height
// (locked head + editable input + locked tail): the schematic (flex:1) shrinks when the preset
// library expands, and tracks manual textarea resize.
function syncSchemHeight(){
  const l=document.querySelector('.cedit-l');
  const r=document.querySelector('.cedit-r');
  if(!l || !r) return;
  // l.offsetHeight omits the locked-tail <pre>'s bottom padding (baseline quirk); add it back from
  // the computed style so this stays correct if the padding ever changes (was a hard-coded +7).
  const tail=document.getElementById('lockTail');
  const pad=tail ? parseFloat(getComputedStyle(tail).paddingBottom)||0 : 0;
  r.style.height = (l.offsetHeight + pad) + 'px';
}

async function run(){
  const btn = document.getElementById('runBtn'); btn.disabled = true;
  setStatus('Running solver (Xyce + FEM + WR)...', '');
  const sent = collect();
  try {
    const res = await fetch('/run', {method:'POST', headers:{'Content-Type':'application/json'},
                                     body: JSON.stringify(sent)});
    const j = await res.json();
    document.getElementById('log').textContent = (j.log||'') + (j.stderr? '\\n--- stderr ---\\n'+j.stderr : '') + (j.trace? '\\n'+j.trace : '');
    const tsec = (typeof j.solver_seconds === 'number') ? j.solver_seconds.toFixed(2)+' s' : '';
    if (!j.ok){ setStatus('Error: '+(j.error||'unknown')+(tsec?' (after '+tsec+')':''), 'err'); }
    else {
      setStatus('Done in '+tsec+'.', 'ok');
      lastRun = {params: sent, summary: j.summary||{}, plots: j.plots||{}};
      showSummary(j.summary||{});
      setPlot('p_voltage', j.plots.voltage);
      setPlot('p_current', j.plots.current);
      setPlot('p_wr', j.plots.wr);
      loadCircuit();  // the solver regenerated the netlist; refresh the schematic
    }
  } catch(e){ setStatus('Request failed: '+e, 'err'); }
  finally { btn.disabled = false; }
}

let lastRun = null;
async function exportCsv(){
  const st = document.getElementById('csvStatus');
  if (!lastRun){ st.textContent = 'Run a simulation first.'; return; }
  const path = (document.getElementById('csvPath').value||'').trim();
  if (!path){ st.textContent = 'Enter a path.'; return; }
  st.textContent = 'Saving...';
  try {
    const res = await fetch('/export_csv', {method:'POST', headers:{'Content-Type':'application/json'},
      body: JSON.stringify({path: path, params: lastRun.params, summary: lastRun.summary, plots: lastRun.plots})});
    const j = await res.json();
    st.textContent = j.ok ? ('Saved to '+j.path) : ('Error: '+(j.error||'unknown'));
  } catch(e){ st.textContent = 'Failed: '+e; }
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
      syncSchemHeight();  // locked text just changed the left-block height -> recap the schematic
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
  const ctrls=document.querySelector('.prop-rows');   // (was '#controls', removed in the redesign)
  if(ctrls){ ctrls.addEventListener('input', applyVisibility); ctrls.addEventListener('change', applyVisibility); }
  // Help hover is document-wide now (help icons live in both the controls and the circuit box).
  document.addEventListener('mouseover', e=>{ if(e.target.classList.contains('help')) helpShow(e.target); });
  document.addEventListener('mouseout',  e=>{ if(e.target.classList.contains('help')) helpHide(); });
  buildLibrary(); applyVisibility(); loadCircuit();
  // Keep the schematic column capped to the full circuit-side (left) block height. Observe the block
  // itself so ANY height change resyncs: async locked-text fill, web-font reflow, or manual resize --
  // this fixes the wrong height on first paint (previously only fixed by a window resize).
  const l=document.querySelector('.cedit-l');
  if(l && window.ResizeObserver){ new ResizeObserver(syncSchemHeight).observe(l); }
  const lf=document.getElementById('libFold');
  if(lf){ lf.addEventListener('toggle', syncSchemHeight); }
  window.addEventListener('resize', syncSchemHeight);
  if(document.fonts && document.fonts.ready){ document.fonts.ready.then(syncSchemHeight); }
  syncSchemHeight();
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
