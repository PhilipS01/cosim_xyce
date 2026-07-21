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
import itertools
import json
import math
import shutil
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
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

# Compiled solver binary. On Windows the MinGW/Clang + Make build emits "main.exe" (GCC/Clang append
# .exe when the -o name has no extension); everywhere else it's "main". Resolve it here so the build
# check and the run/emit calls all agree on the name.
MAIN_EXE = "main.exe" if os.name == "nt" else "main"


def main_binary():
    return os.path.join(HERE, MAIN_EXE)

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
    ("t_end",                           "Sim duration (s)",             2.0e-2,   "float", None),
    ("N_field_windows",                 "Field windows (total)",        50,       "int",   (1, 400, 1)),
    ("N_field_eval_intervals",          "FEM eval intervals / window",  1,        "int",   (1, 64, 1)),
    ("N_xyce_samples",                  "Xyce solution samples",        100,      "int",   (1, 400, 1)),
    ("WRmaxSteps",                      "WR max iterations",            20,       "int",   (1, 100, 1)),
    ("WR_tolerance",                    "WR tolerance",                 1.0e-3,   "float", None),
    ("wr_convergence_method",           "WR convergence metric",        1,        "choice",
        {0: "waveform L1", 1: "terminal scalar"}),
    ("coupling_mode",                   "Coupling direction",           0,        "choice",
        {0: "voltage-driven", 1: "current-driven"}),
    ("reconstruct_mode",                "Field reconstruction",         0,        "choice",
        {0: "pointwise (secant)", 1: "linear ramp"}),
    ("interface_form",                  "Interface stamping",           0,        "choice",
        {0: "Thevenin (V source)", 1: "Norton (I source)"}),
    ("use_t_floor",                     "Secant t_floor",               1,        "choice",
        {1: "on (guard 1/0)", 0: "off (bare dt)"}),
    ("t_floor_frac",                    "t_floor / window",             0.01,     "float", None),
    ("seam_average",                    "Window-seam handoff",          0,        "choice",
        {0: "one-sided (V<-ckt, I<-fld)", 1: "midpoint (average)"}),
    ("validation_mode",                 "Validation mode",              0,        "choice",
        {0: "WR co-sim", 1: "monolithic (real R_FEM/L_FEM)"}),
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
    "P1: Sine V + RL": {"coupling_mode": 0,
                        "circuit_spec": "VSIN Bemf s 0 1 50\nR Rs s cm0 6e-3\nL Ls cm0 p 1.6e-7\n"},
    # Bare current source directly on the port: series R/L/C are meaningless for a current drive
    # (the current is forced regardless) and an ideal I-source in series with L is degenerate.
    "P2: Sine I (bare)": {"coupling_mode": 1,
                          "circuit_spec": "ISIN Bemf 0 p 1 50\n"},
    # window 1 straddles the whole ramp edge (stiff transient) -> more WR iters (WRmaxSteps=40)
    "P3: Step/ramp V + RL": {"t_end": 2.0e-2, "N_field_windows": 50,
                             "coupling_mode": 0, "WRmaxSteps": 40,
                             "circuit_spec": "VPULSE Vemf s 0 0 1 0 1e-4\nR Rs s cm0 6e-3\nL Ls cm0 p 1.6e-7\n"},
    # Switch circuits (SW = time-gated switch, closed during [tclose, topen); emitted as a Xyce native
    # Generic Switch, S device + .MODEL SWITCH). 2-way (SPDT) switch: wiper w throws between the source
    # branch (node a, drv closed early) and the freewheel branch (fw closed after), with cap w->0.
    # P4 freewheels directly w->p; P5 via a short branch (node b -> R shrt -> p). Ron/C set the damping:
    # an undamped cap<->coil freewheel driven on resonance rings hard (see the validation-vs-cosim study).
    "P4: 2-way switch (sine U, C)": {"t_end": 1.6e-3, "N_field_windows": 100,
        "coupling_mode": 0, "WRmaxSteps": 40, "L_FEM": 1.68e-7, "R_FEM": 5e-4,
        "circuit_spec": "VSIN Bemf a p 1 20000\nC Csw 0 w 0.37e-3\n"
                        "SW drv w a 0 1e-3 1e-3 1e12 0.1e-3\nSW fw w p 1e-3 1e30 1e-3 1e12 0.1e-3\n"},
    "P5: 2-way switch (DC U, C)": {"t_end": 1.6e-3, "N_field_windows": 100,
        "coupling_mode": 0, "WRmaxSteps": 40, "L_FEM": 1.68e-7, "R_FEM": 5e-4,
        "circuit_spec": "VDC Vemf a p 1\nC Csw 0 w 0.37e-3\n"
                        "SW drv w a 0 1e-3 1e-3 1e12 0.1e-3\nSW fw w p 1e-3 1e30 1e-3 1e12 0.1e-3\n"},
    "P6: 2-way switch (AC vs R)": {"t_end": 2.0e-2, "N_field_windows": 50,
        "coupling_mode": 0, "WRmaxSteps": 40,
        "circuit_spec": "VSIN Bemf p bac 1 50\nR Rload p br 1e4\n"
                        "SW ac bac 0 0 6e-3 1e-3 1e9 1e-5\nSW rd br 0 6e-3 1e30 1e-3 1e9 1e-5\n"},
}


# ---------------------------------------------------------------------------
# Running the solver
# ---------------------------------------------------------------------------
# Conditional visibility: key -> list of AND-condition dicts; a control is shown iff ANY dict fully
# matches the current control values (OR-of-ANDs). Keys absent here are always visible. The circuit
# side is always the custom node-graph spec (circuit_spec.txt); the run is absolute-end-time only
# ("Sim duration" = t_end). What remains conditional: N_field_eval_intervals only matters for the
# pointwise-secant reconstruction (0) -- linear (1) forces 1 solve/window and ignores it -- so it's
# shown only then. t_floor_frac only matters when the secant t_floor guard is on (use_t_floor=1);
# with bare dt (0) the floor is irrelevant. validation_mode=1 (monolithic reference: the true field
# stamped as real R_FEM/L_FEM devices, one Xyce transient, no WR) makes the whole WR/coupling/secant
# machinery irrelevant, so those controls are disabled while it is on -- only t_end/windows/samples
# (the tran grid), R_FEM/L_FEM (the devices) and the nonlinearity fields stay live.
VISIBLE_WHEN = {
    "N_field_eval_intervals": [{"reconstruct_mode": [0], "validation_mode": [0]}],
    "t_floor_frac":           [{"use_t_floor": [1], "validation_mode": [0]}],
    "coupling_mode":          [{"validation_mode": [0]}],
    "reconstruct_mode":       [{"validation_mode": [0]}],
    "interface_form":         [{"validation_mode": [0]}],
    "seam_average":           [{"validation_mode": [0]}],
    "wr_convergence_method":  [{"validation_mode": [0]}],
    "WRmaxSteps":             [{"validation_mode": [0]}],
    "WR_tolerance":           [{"validation_mode": [0]}],
    "use_t_floor":            [{"validation_mode": [0]}],
    "R_ROM":                  [{"validation_mode": [0]}],
    "L_ROM":                  [{"validation_mode": [0]}],
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
        "<tr><td><code>VPULSE/IPULSE name a b v1 v2 td tr</code></td><td>pulse source (single edge)</td>"
        "<td><code>a b</code> nodes &middot; <code>v1</code> initial level &middot; <code>v2</code> pulsed "
        "level &middot; <code>td</code> delay (s) &middot; <code>tr</code> rise time (s); one transition, "
        "then holds <code>v2</code></td></tr>"
        "<tr><td><code>VPWM/IPWM name a b v1 v2 freq duty</code></td><td>PWM pulse train</td>"
        "<td><code>a b</code> nodes &middot; <code>v1</code> low level &middot; <code>v2</code> high level "
        "&middot; <code>freq</code> switching frequency (Hz) &middot; <code>duty</code> duty cycle 0&ndash;1 "
        "(fraction of each period held high)</td></tr>"
        "<tr><td><code>VPWL/IPWL name a b t1 v1 t2 v2 …</code></td><td>piecewise-linear source</td>"
        "<td><code>a b</code> nodes &middot; <code>t1 v1 t2 v2 …</code> (time&nbsp;s, level&nbsp;V/A) "
        "breakpoints, linearly interpolated between</td></tr>"
        "<tr><td><code>SW name a b tclose topen [Ron Roff trise]</code></td>"
        "<td>time-gated switch</td>"
        "<td><code>a b</code> nodes &middot; <code>tclose</code> close time (s) &middot; <code>topen</code> "
        "open time (s) &mdash; big topen (e.g. <code>1e30</code>) stays closed to the end &middot; "
        "<code>Ron</code>/<code>Roff</code> closed/open R &mdash; optional; if omitted Xyce's "
        "defaults apply (<code>RON=1</code>, <code>ROFF=1e12</code>) &middot; "
        "<code>trise</code> transition (=1e-5). Emitted as a Xyce native Generic Switch "
        "(<code>S</code> device + <code>.MODEL SWITCH</code>)</td></tr>"
        "</table>"
    ),
    "probes": (
        "<div class='hh'>Probes</div>"
        "<div class='hn'>Add extra output points plotted after a run. <b>V</b> probes a node "
        "voltage &mdash; against ground (<code>V(a)</code>) or, by picking a second node, the "
        "<i>differential</i> <code>V(a,b)</code> = V(a)&minus;V(b). <b>I</b> probes the branch current "
        "of an R/L/C element. Targets are read from your circuit spec above.</div>"
        "<div class='hn'>Emitted into the netlist's <code>.print</code> and captured per window into "
        "<code>Probes_solution.prn</code>; shown as the <b>Probe voltages</b> / <b>Probe currents</b> "
        "plots in Results. The port <code>V(p)</code> and interface current are already plotted.</div>"
    ),
    "sweep_jobs": (
        "<div class='hh'>Parallel jobs</div>"
        "<div class='hn'>How many sweep points solve at once. Each point runs the solver in its own "
        "temporary working directory, so concurrent runs don't clash. <code>0</code> = auto "
        "(all CPU cores). Capped to the number of points.</div>"
        "<div class='hn'>Speeds up multi-point sweeps roughly linearly until you saturate cores or "
        "memory. Set to <code>1</code> for a serial run (lowest memory / easiest to read logs).</div>"
    ),
    "sweep_export": (
        "<div class='hh'>Export sweep</div>"
        "<div class='hn'>Writes the whole sweep to CSV: one row per run (run index, one column per "
        "swept parameter, then the result metrics). Also writes <code>&lt;stem&gt;_iters.csv</code> "
        "(WR iterations per window, one row per run) and the current plots (with reference lines) into "
        "<code>&lt;stem&gt;_plots/</code>.</div>"
        "<div class='hn'>Server-side path, relative to where the studio runs. Overwrites an existing file.</div>"
    ),
    "sweep_vlines": (
        "<div class='hh'>Reference lines</div>"
        "<div class='hn'>Comma-separated values drawn as dashed vertical lines on the sweep plots "
        "(handy for marking a crossover). Each entry is arithmetic over the config parameters "
        "(<code>R_ROM</code>, <code>L_ROM</code>, &hellip;) <b>and your circuit-spec R/L/C element "
        "names</b> &mdash; e.g. <code>Rs+Ls</code> (series R + L above), <code>R_ROM</code>, or a plain "
        "number like <code>0.01</code>.</div>"
        "<div class='hn'>Positioned where the x-axis parameter equals the value (interpolated); entries "
        "outside the swept range are skipped. On the WR-iterations heatmap and line plots the line marks "
        "that x-location; on grid heatmaps it's drawn on the first (x) parameter axis.</div>"
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
    "window_width": (
        "<div class='hh'>WR window width</div>"
        "<div class='hn'>Width of one waveform-relaxation time window (s). <b>Derived</b>, not a config "
        "key: <code>width = t_end / N_field_windows</code>. Two-way linked &mdash; edit the width and the "
        "window count is set to <code>round(t_end / width)</code> (&ge;1), then the width snaps to the true "
        "<code>t_end / N</code> (it may not divide evenly). Narrower windows = more windows = more "
        "checkpoint/restart seams but easier per-window WR convergence.</div>"
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
        "</table>"
        "<div class='hn'><b>pointwise (secant)</b> with <code>FEM eval intervals / window = 1</code> "
        "collapses to <b>linear ramp</b>: one interval leaves only the carried start and the window end, "
        "so the accumulated secant is a single straight segment. Raise the eval intervals for it to actually "
        "follow the curve.</div>"
        "<div class='hn'>Both coupling directions. Both modes carry the seam (C0-continuous). Accuracy is "
        "within ~1&ndash;2% between them; the extra solves buy little. Recommend <b>linear</b> "
        "(1 field solve per window).</div>"
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
    "seam_average": (
        "<div class='hh'>Window-seam handoff</div>"
        "<div class='hn'>Which terminal value seeds the next window (both solvers agree &lt; WR "
        "tolerance at the seam).<br>"
        "<b>one-sided</b> (default): V0 = circuit V(p), I0 = field I.<br>"
        "<b>midpoint</b>: V0 = &frac12;(V_circuit+V_field), I0 = &frac12;(I_circuit+I_field). Only the "
        "carried seeds are averaged; the Xyce restart checkpoint is unchanged. Seam-blend test.</div>"
    ),
    "use_t_floor": (
        "<div class='hh'>Secant denominator floor (t_floor)</div>"
        "<div class='hn'>The Bfield impedance Z = Rrom + Lrom/dt uses dt = time&minus;t_abs_start.<br>"
        "<b>on</b> (default): dt &rarr; MAX(time&minus;t_abs_start, t_floor); guards the 1/0 at the exact "
        "window start (dt=0 &rarr; Lrom/0 singularity). t_floor is set by <b>t_floor / window</b>.<br>"
        "<b>off</b>: bare dt; Z &rarr; &infin; at the window start. Test only.</div>"
    ),
    "t_floor_frac": (
        "<div class='hh'>t_floor as a fraction of the window</div>"
        "<div class='hn'>Sets the secant-denominator floor: t_floor = (this) &times; t_window (only "
        "used when <b>Secant t_floor</b> is on). Scale-free &mdash; the floor tracks window length. The "
        "floor value is accuracy/stability-neutral (floor-sweep), so this is a conditioning knob, not "
        "physics. Default 0.01 (= t_window/100).</div>"
    ),
    "N_xyce_samples": (
        "<div class='hh'>Xyce solution samples</div>"
        "<div class='hn'>How finely the Xyce solution is sampled: sets the .tran print cadence "
        "dt_print = t_window/N (the raw wr_circuit.cir.prn rows) AND the interface PWL resolution "
        "(V(p)&rarr;vf_prev_k.pwl and I(Vmeas)&rarr;i_prev_k.pwl, <b>both coupling directions</b>). "
        "Decoupled from t_floor (that is now <b>t_floor / window</b>).</div>"
    ),
    "coupling_mode": (
        "<div class='hh'>Coupling direction</div>"
        "<div class='hn'>voltage-driven (Dirichlet): circuit sets V(p), field returns current; "
        "matched-secant Bfield.<br>current-driven (Neumann): circuit sets I(Vmeas), field returns "
        "V_field; plain voltage source &mdash; removes the high-frequency window-start V(p) spike on "
        "current-source circuits.</div>"
    ),
    "validation_mode": (
        "<div class='hh'>Validation mode</div>"
        "<div class='hn'><b>WR co-sim</b> (default) &mdash; the normal coupled run: the field is a "
        "reduced-order matched-secant <code>Bfield</code> source, iterated to a fixpoint per window "
        "(waveform relaxation).</div>"
        "<div class='hn'><b>monolithic</b> &mdash; replaces <code>Bfield</code> with the TRUE field as "
        "real Xyce devices (<code>R_FEM</code> + <code>L_FEM</code> in series on the port branch) and "
        "solves the whole circuit as ONE transient over [0, t_end] &mdash; no WR loop, no field solver, "
        "no coupling. Produces a <i>reference</i>; flip back to WR co-sim and check the coupled solution "
        "reproduces it. WR/coupling/secant knobs are disabled here (irrelevant) and the convergence plot "
        "is empty. Linear field only &mdash; magnetic saturation is ignored (a plain inductor can't "
        "reproduce the flux law).</div>"
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
}


def write_config(params, workdir=HERE):
    path = os.path.join(workdir, "sim_config.txt")
    with open(path, "w") as f:
        f.write("# generated by sim_ui.py\n")
        for k in DEFAULTS:
            v = params.get(k, DEFAULTS[k])
            if KINDS[k] in ("int", "choice"):
                f.write(f"{k} = {int(round(float(v)))}\n")
            else:
                f.write(f"{k} = {float(v):.10g}\n")
        # The circuit side is always authored via the text spec (circuit_spec.txt); the C++ generator
        # reads circuit_spec.txt directly, so no circuit-kind/source keys are needed here.
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
    if not os.path.exists(main_binary()):
        subprocess.run(["make"], cwd=HERE, check=True,
                       capture_output=True, text=True)


def run_solver(workdir=HERE):
    ensure_built()
    t0 = time.perf_counter()
    proc = subprocess.run([main_binary()], cwd=workdir,
                          capture_output=True, text=True, timeout=600)
    elapsed = time.perf_counter() - t0
    return proc.returncode, proc.stdout, proc.stderr, elapsed


def emit_netlist():
    """Regenerate wr_circuit.cir from the current sim_config.txt without solving
    ('main emit'), so the rendered schematic reflects the current config."""
    ensure_built()
    subprocess.run([main_binary(), "emit"], cwd=HERE,
                   capture_output=True, text=True, timeout=60)


# Example custom node-graph spec. Reserved nodes: p (port), 0 (ground).
DEFAULT_SPEC = """\
# Custom node-graph circuit. Reserved nodes: p = port (field attaches here), 0 = ground.
# <TYPE> <name> <nodeA> <nodeB> <params...>
#   R/L/C name a b val | {V,I}SIN name a b amp freq | {V,I}DC name a b val
#   {V,I}PULSE name a b v1 v2 td tr | {V,I}PWM name a b v1 v2 freq duty  (pulse train)
#   {V,I}PWL name a b t1 v1 t2 v2 ...  (multi-step)
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


def write_spec(params, workdir=HERE):
    """Write the custom node-graph spec to circuit_spec.txt (only when non-empty), so the C++
    generator reads it as the circuit side."""
    spec = params.get("circuit_spec")
    if spec:
        with open(os.path.join(workdir, "circuit_spec.txt"), "w") as f:
            f.write(spec)


def write_probes(params, workdir=HERE):
    """Write the user's output probes to probes.txt (one Xyce print token per line, e.g. 'V(a)',
    'I(Rr1)'). The C++ LoadProbes() appends them to the netlist .print and captures a column per
    probe into Probes_solution.prn. Always (re)written so removing all probes clears the file."""
    probes = params.get("probes") or []
    if isinstance(probes, str):
        probes = [p for p in probes.replace(",", "\n").split("\n")]
    toks, seen = [], set()
    for p in probes:
        t = str(p).strip()
        if t and t not in seen:
            seen.add(t); toks.append(t)
    with open(os.path.join(workdir, "probes.txt"), "w") as f:
        f.write("\n".join(toks) + ("\n" if toks else ""))


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
    # P4 / P5: SPDT wiper w (cap w->gnd) throws between the source branch (a, up) and the freewheel
    # straight to the port (fw: w->p); same topology, source differs (sine / DC). Interface on the right.
    frozenset({"a", "w", "p", "nx", "0"}): {
        "nets": {"w": (40, 120), "a": (160, 40), "p": (300, 120),
                 "nx": (420, 120), "0": (330, 240)},
        # w is a 3-way junction: drv up (to a), fw straight right (to p), cap down (to gnd).
        "routes": {frozenset({"w", "0"}): [(40, 120), (40, 240), (330, 240)],
                   frozenset({"w", "a"}): [(40, 120), (40, 40), (160, 40)],
                   frozenset({"a", "p"}): [(160, 40), (300, 40), (300, 120)],
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
    write_config(params); write_spec(params); write_probes(params); emit_netlist()
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


def read_outputs(workdir=HERE):
    out = {}
    # Field_solution.prn: Index TIME V(FIELD) I(FIELD)  -> synchronisation endpoints
    t, c = _read_columns(os.path.join(workdir, "Field_solution.prn"), 4)
    out["field"] = {"t": t, "V": c[0] if c else np.array([]),
                    "I": c[1] if len(c) > 1 else np.array([])}
    # Field_waveform_solution.prn: Index TIME V(FIELD) I(FIELD) -> reconstructed waveforms
    t, c = _read_columns(os.path.join(workdir, "Field_waveform_solution.prn"), 4)
    out["field_wave"] = {"t": t, "V": c[0] if c else np.array([]),
                         "I": c[1] if len(c) > 1 else np.array([])}
    # Circuit_solution.prn: Index TIME V(P) V(NX) I(VMEAS)
    t, c = _read_columns(os.path.join(workdir, "Circuit_solution.prn"), 5)
    out["circuit"] = {"t": t,
                      "Vp":  c[0] if len(c) > 0 else np.array([]),
                      "Vnx": c[1] if len(c) > 1 else np.array([]),
                      "I":   c[2] if len(c) > 2 else np.array([])}
    # WR_error.txt: comma separated "Time, RelErr, N_iter, Converged"
    out["wr"] = read_wr_error(os.path.join(workdir, "WR_error.txt"))
    # Probes_solution.prn: Index TIME <probe1> <probe2> ... (column names from the header)
    out["probes"] = _read_named_columns(os.path.join(workdir, "Probes_solution.prn"))
    return out


def _read_named_columns(path):
    """Read a whitespace .prn whose header names the data columns after Index/TIME.
    Returns {"t": array, "names": [...], "series": {name: array}} (empty if absent)."""
    res = {"t": np.array([]), "names": [], "series": {}}
    if not os.path.exists(path):
        return res
    with open(path) as f:
        header = f.readline().split()
        names = header[2:]                      # drop 'Index' + 'TIME'
        t, cols = [], [[] for _ in names]
        for line in f:
            s = line.strip()
            if not s or s.startswith("End"):
                continue
            parts = s.split()
            if len(parts) < 2 + len(names):
                continue
            try:
                vals = [float(p) for p in parts[:2 + len(names)]]
            except ValueError:
                continue
            t.append(vals[1])
            for i in range(len(names)):
                cols[i].append(vals[2 + i])
    res["t"] = np.array(t)
    res["names"] = names
    res["series"] = {n: np.array(c) for n, c in zip(names, cols)}
    return res


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

    # 4/5) User probes -> a voltage figure (V(...) tokens) + a current figure (I(...) tokens)
    pr = data.get("probes") or {"t": np.array([]), "names": [], "series": {}}
    if pr["t"].size and pr["names"]:
        v_names = [n for n in pr["names"] if n.upper().startswith("V(")]
        i_names = [n for n in pr["names"] if n.upper().startswith("I(")]
        if v_names:
            fig, ax = plt.subplots(figsize=(8, 3.2))
            for n in v_names:
                ax.plot(pr["t"], pr["series"][n], lw=0.9, label=n)
            ax.set_xlabel("time (s)"); ax.set_ylabel("voltage (V)")
            ax.set_title("Probe voltages"); ax.grid(True, alpha=.3); ax.legend(fontsize=8)
            plots["probe_v"] = _png(fig)
        if i_names:
            fig, ax = plt.subplots(figsize=(8, 3.2))
            for n in i_names:
                ax.plot(pr["t"], pr["series"][n], lw=0.9, label=n)
            ax.set_xlabel("time (s)"); ax.set_ylabel("current (A)")
            ax.set_title("Probe currents"); ax.grid(True, alpha=.3); ax.legend(fontsize=8)
            plots["probe_i"] = _png(fig)
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


def export_sweep_csv(path, keys, rows, plots=None):
    """Write a whole parameter sweep: one CSV row per swept point (run index,
    one column per swept parameter, then the result metrics). Also writes a
    companion '<stem>_iters.csv' (WR iterations per window: one row per run) and
    dumps any plot PNGs into '<stem>_plots/'. Returns (csv_path, [written files])."""
    metric_cols = ["ok", "max_WR_iterations", "mean_WR_iterations", "total_xyce_solves",
                   "worst_WR_error", "all_converged", "solver_seconds",
                   "final_I_field", "final_V_field", "final_time_s", "windows", "error"]
    d = os.path.dirname(path)
    if d:
        os.makedirs(d, exist_ok=True)

    header = ["index"] + list(keys) + metric_cols
    with open(path, "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(header)
        for i, r in enumerate(rows):
            vals = r.get("vals") or {keys[0]: r.get("value")}
            wr.writerow([i] + [vals.get(k) for k in keys] + [r.get(c) for c in metric_cols])
    written = [os.path.abspath(path)]

    # Per-window WR iterations: rows = run index, columns = window 1..W (ragged -> blank).
    wmax = max((len(r.get("wr_nit") or []) for r in rows), default=0)
    if wmax:
        ipath = os.path.splitext(path)[0] + "_iters.csv"
        with open(ipath, "w", newline="") as f:
            wr = csv.writer(f)
            wr.writerow(["index"] + list(keys) + [f"w{w+1}" for w in range(wmax)])
            for i, r in enumerate(rows):
                vals = r.get("vals") or {keys[0]: r.get("value")}
                nit = r.get("wr_nit") or []
                wr.writerow([i] + [vals.get(k) for k in keys] +
                            [nit[w] if w < len(nit) else "" for w in range(wmax)])
        written.append(os.path.abspath(ipath))

    if plots:
        stem = os.path.splitext(os.path.basename(path))[0]
        plotdir = os.path.join(d, stem + "_plots") if d else (stem + "_plots")
        os.makedirs(plotdir, exist_ok=True)
        for name, uri in plots.items():
            dest = os.path.join(plotdir, f"{name}.png")
            if _save_data_uri_png(uri, dest):
                written.append(os.path.abspath(dest))
    return os.path.abspath(path), written


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


def _all_config_numeric():
    """Every numeric key=value in sim_config.txt (not just UI-known keys), so
    reference-line expressions can name any config parameter."""
    d = {}
    path = os.path.join(HERE, "sim_config.txt")
    if os.path.exists(path):
        with open(path) as f:
            for line in f:
                s = line.split("#", 1)[0].strip()
                if "=" not in s:
                    continue
                k, v = (x.strip() for x in s.split("=", 1))
                try:
                    d[k] = float(v)
                except ValueError:
                    pass
    return d


def _spec_element_values(spec):
    """Map R/L/C element names in a circuit spec to their numeric value, so
    reference-line expressions can name real circuit elements (e.g. Rs, Ls).
    Line form: '<TYPE> <name> <nodeA> <nodeB> <value> ...'; only R/L/C keep a
    plain scalar value in that slot. `spec` = the circuit-spec text (or None)."""
    d = {}
    if not spec:
        spec = read_spec()
    for line in (spec or "").splitlines():
        s = line.split("#", 1)[0].strip()
        if not s:
            continue
        parts = s.split()
        if len(parts) >= 5 and parts[0] in ("R", "L", "C"):
            try:
                d[parts[1]] = float(parts[4])
            except ValueError:
                pass
    return d


# Last completed sweep, kept so Results-section reference lines can be re-drawn
# without re-solving: {"keys","rows","mode","free_keys","params"}.
_LAST_SWEEP = {}


def _vline_env(params):
    """Build the evaluation namespace for reference-line expressions:
    sim_config.txt numerics, then circuit-spec R/L/C element values by name
    (e.g. Rs, Ls), then the live params (highest precedence), plus math helpers."""
    env = {"__builtins__": {}}
    src = dict(_all_config_numeric())
    src.update(_spec_element_values((params or {}).get("circuit_spec")))
    src.update(params or {})
    for k, v in src.items():
        try:
            env[k] = float(v)
        except (TypeError, ValueError):
            pass
    env.update(abs=abs, min=min, max=max, sqrt=math.sqrt, pi=math.pi)
    return env


def _split_exprs(exprs):
    if not exprs:
        return []
    if isinstance(exprs, str):
        exprs = re.split(r"[,\n]", exprs)
    return [str(e).strip() for e in exprs if str(e).strip()]


def eval_vlines(exprs, params):
    """Turn reference-line expressions into [(value, label), ...]. Each expression
    is evaluated as arithmetic over the config params + circuit-spec element names
    (see _vline_env). Non-numeric / failing expressions are skipped; the original
    text becomes the line label."""
    env = _vline_env(params)
    out = []
    for s in _split_exprs(exprs):
        try:
            out.append((float(eval(s, env)), s))   # arithmetic only: builtins stripped
        except Exception:
            continue
    return out


def eval_vlines_report(exprs, params):
    """Like eval_vlines but reports every expression's outcome for live feedback:
    [{expr, ok, value} | {expr, ok:False, error}]."""
    env = _vline_env(params)
    out = []
    for s in _split_exprs(exprs):
        try:
            out.append({"expr": s, "ok": True, "value": float(eval(s, env))})
        except Exception as e:
            out.append({"expr": s, "ok": False, "error": type(e).__name__})
    return out


def build_sweep_points(specs, mode):
    """Turn a list of per-parameter sweep specs into (keys, points).

    A "free" spec is an independently swept range: {key, min, max, steps, scale}.
    A "linked" spec tracks another swept parameter as a fixed percentage and adds
    no new dimension: {key, link:{base, pct}} -> value = pct/100 * base_value.

    Modes (applied to the FREE specs only):
      "single"   -- one free parameter, a plain 1D sweep.
      "parallel" -- N free parameters advance in lock-step (same step count);
                    point i is (p1[i], p2[i], ...). Lists truncate to min length.
      "grid"     -- every combination of the free parameters (Cartesian product).
    Linked values are appended to each point after the free values.
    """
    free = [s for s in specs if not s.get("link")]
    linked = [s for s in specs if s.get("link")]
    if not free:
        raise ValueError("a sweep needs at least one independent (range) parameter")
    free_keys = [s["key"] for s in free]
    per = [sweep_values(float(s["min"]), float(s["max"]), int(s["steps"]),
                        s.get("scale", "linear"), KINDS[s["key"]]) for s in free]
    if len(free) == 1 or mode == "single":
        base_pts = [(v,) for v in per[0]]
        free_keys = free_keys[:1]
    elif mode == "parallel":
        n = min(len(v) for v in per)
        base_pts = [tuple(v[i] for v in per) for i in range(n)]
    else:  # grid
        base_pts = [tuple(pt) for pt in itertools.product(*per)]

    fidx = {k: i for i, k in enumerate(free_keys)}
    for s in linked:
        if s["link"]["base"] not in fidx:
            raise ValueError(f"linked parameter '{s['key']}' references "
                             f"'{s['link']['base']}', which is not an independently swept parameter")
    keys = free_keys + [s["key"] for s in linked]
    points = []
    for pt in base_pts:
        vals = list(pt)
        for s in linked:
            frac = float(s["link"]["pct"]) / 100.0
            vals.append(pt[fidx[s["link"]["base"]]] * frac)
        points.append(tuple(vals))
    return keys, points, free_keys


def _sweep_point(base_params, keys, pt):
    """Solve one sweep point in an isolated temp directory so concurrent points
    never clobber each other's config/output files. Returns the metrics row."""
    d = tempfile.mkdtemp(prefix="wr_sweep_")
    try:
        # Seed the isolated dir with the current input files, so a caller that
        # doesn't pass circuit_spec/probes (e.g. the CLI) still has them.
        for fn in ("circuit_spec.txt", "probes.txt"):
            src = os.path.join(HERE, fn)
            if os.path.exists(src):
                shutil.copy(src, os.path.join(d, fn))
        params = dict(base_params)
        for k, v in zip(keys, pt):
            params[k] = v
        write_config(params, d)
        write_spec(params, d)
        write_probes(params, d)
        rc, stdout, stderr, elapsed = run_solver(d)
        row = {"vals": {k: v for k, v in zip(keys, pt)},
               "value": pt[0], "solver_seconds": round(elapsed, 4)}
        if len(pt) > 1:
            row["value2"] = pt[1]
        if rc != 0:
            row["ok"] = False
            row["error"] = f"solver exited {rc}"
            row["log"] = "\n".join((stdout or stderr or "").splitlines()[-6:])
        else:
            data = read_outputs(d)
            row.update(scalar_summary(data))
            row["wr_nit"] = [int(x) for x in data["wr"]["nit"]]   # per-window iteration counts
            row["ok"] = True
        return row
    finally:
        shutil.rmtree(d, ignore_errors=True)


def run_sweep(base_params, keys, points, workers=None):
    """Run the solver once per swept point; collect per-point metrics.

    Each point runs in its own temp working directory, so points can execute
    concurrently (`workers` threads; the solver is a subprocess, so the GIL is
    released while it runs). Results are returned in input order regardless of
    completion order. `workers=None` -> os.cpu_count(). For back-compat a single
    string key + list of scalars is also accepted (classic 1-parameter sweep).
    """
    if isinstance(keys, str):  # legacy call: run_sweep(base, key, values)
        keys = [keys]
        points = [(v,) for v in points]
    n = len(points)
    if workers is None:
        workers = os.cpu_count() or 1
    workers = max(1, min(int(workers), n or 1))
    if workers == 1:
        return [_sweep_point(base_params, keys, pt) for pt in points]
    rows = [None] * n
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(_sweep_point, base_params, keys, pt): i
                for i, pt in enumerate(points)}
        for fut in as_completed(futs):
            rows[futs[fut]] = fut.result()
    return rows


def make_sweep_plots(keys, rows, mode="single", free_keys=None, vlines=None):
    """Convergence-study figure. Dispatches on the sweep shape:
      single / parallel -> 2x2 metric-vs-parameter line plots.
      grid              -> 2x2 metric heatmaps over the two swept parameters.
    `keys` may be a single string (legacy) or a list of parameter names.
    `vlines` = list of (value, label) reference lines drawn on the x-axis."""
    if isinstance(keys, str):
        keys = [keys]
    # Dimensionality is set by the INDEPENDENT (free) parameters; percentage-linked
    # ones track a base and add no axis.
    dims = list(free_keys) if free_keys else list(keys)
    if mode == "grid" and len(dims) == 2:
        out = _make_grid_plots(dims[0], dims[1], rows, vlines)
    elif mode == "grid" and len(dims) > 2:
        out = _make_index_plots(dims, rows, vlines)   # >2D: heatmap not meaningful
    else:
        out = _make_line_plots(keys, rows, mode, vlines)
    out.update(_make_iters_heatmap(dims, rows, mode, vlines))  # WR-iterations-per-window colormap
    return out


def _make_line_plots(keys, rows, mode, vlines=None):
    """2x2 metric-vs-swept-parameter figure (single or parallel/lock-step)."""
    key = keys[0]
    ok = [r for r in rows if r.get("ok")]
    if not ok:
        return {}
    x = np.array([r["value"] for r in ok], dtype=float)
    label = key   # short axis label (e.g. "R_ROM"); full detail goes in the title
    title_label = LABELS.get(key, key)
    if mode == "parallel" and len(keys) > 1:
        title_label += " (lock-step with " + ", ".join(keys[1:]) + ")"

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

    for ax in axes.flat:   # x-axis is the parameter value -> value is a data coordinate
        _draw_vlines(ax, vlines, data_x=True)
    fig.suptitle(f"Convergence study: sweep of {title_label}", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    return {"sweep": _png(fig)}


def _make_grid_plots(kx, ky, rows, vlines=None):
    """2x2 heatmaps of key metrics over a 2-parameter grid (all combinations)."""
    ok = [r for r in rows if r.get("ok")]
    if not ok:
        return {}
    xs = sorted({r["vals"][kx] for r in ok})
    ys = sorted({r["vals"][ky] for r in ok})
    ix = {v: i for i, v in enumerate(xs)}
    iy = {v: i for i, v in enumerate(ys)}
    lx, ly = kx, ky   # short axis labels; full names go in the suptitle

    def grid(name, log=False):
        g = np.full((len(ys), len(xs)), np.nan)
        for r in ok:
            val = r.get(name, np.nan)
            if log and val is not None:
                val = max(float(val), 1e-16)
            g[iy[r["vals"][ky]], ix[r["vals"][kx]]] = val
        return g

    # cell edges (midpoints) so pcolormesh centres cells on the sampled values
    def edges(v):
        v = np.array(v, dtype=float)
        if len(v) == 1:
            d = abs(v[0]) * 0.5 or 0.5
            return np.array([v[0] - d, v[0] + d])
        mid = (v[:-1] + v[1:]) / 2
        return np.concatenate([[2 * v[0] - mid[0]], mid, [2 * v[-1] - mid[-1]]])

    xe, ye = edges(xs), edges(ys)
    panels = [
        ("max_WR_iterations", "WR iterations / window (max)", "viridis", False),
        ("total_xyce_solves", "total Xyce solves", "magma", False),
        ("worst_WR_error", "worst WR rel. error (log)", "inferno", True),
        ("final_I_field", "final I_field (A)", "cividis", False),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(10, 6.6))
    for ax, (name, title, cmap, log) in zip(axes.flat, panels):
        g = grid(name, log)
        norm = None
        if log:
            from matplotlib.colors import LogNorm
            finite = g[np.isfinite(g)]
            if finite.size:
                norm = LogNorm(vmin=finite.min(), vmax=finite.max())
        pcm = ax.pcolormesh(xe, ye, g, cmap=cmap, norm=norm, shading="flat")
        fig.colorbar(pcm, ax=ax, fraction=0.046, pad=0.04)
        ax.set_title(title, fontsize=10)
        ax.set_xlabel(lx); ax.set_ylabel(ly)
        _draw_vlines(ax, vlines, data_x=True)   # x-axis is kx in real units

    fig.suptitle(f"Grid sweep: {LABELS.get(kx, kx)}  x  {LABELS.get(ky, ky)}", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    return {"sweep": _png(fig)}


def _make_index_plots(keys, rows, vlines=None):
    """Fallback for a >2-parameter grid: metrics vs flat run index (a heatmap
    needs 2 axes). The table carries the full parameter tuple per run."""
    ok = [r for r in rows if r.get("ok")]
    if not ok:
        return {}
    x = np.arange(len(ok), dtype=float)
    # x is run index -> map a reference value via the first swept param's per-run values
    colvals = [r["vals"][keys[0]] for r in ok] if keys else None

    def col(name):
        return np.array([r.get(name, np.nan) for r in ok], dtype=float)

    fig, axes = plt.subplots(2, 2, figsize=(10, 6.4))
    ax = axes[0, 0]
    ax.plot(x, col("max_WR_iterations"), "o-", color="tab:blue", label="max")
    ax.plot(x, col("mean_WR_iterations"), "s--", color="tab:cyan", label="mean")
    ax.set_ylabel("WR iterations / window"); ax.set_title("Convergence speed")
    ax.legend(fontsize=8); ax.grid(True, alpha=.3)

    ax = axes[0, 1]
    ax.plot(x, col("total_xyce_solves"), "o-", color="tab:purple")
    ax.set_ylabel("total Xyce solves"); ax.set_title("Cost"); ax.grid(True, alpha=.3)

    ax = axes[1, 0]
    ax.semilogy(x, np.maximum(col("worst_WR_error"), 1e-16), "o-", color="tab:red")
    ax.set_ylabel("worst WR rel. error"); ax.set_title("WR accuracy"); ax.grid(True, alpha=.3)

    ax = axes[1, 1]
    ax.plot(x, col("final_I_field"), "o-", color="tab:orange")
    ax.set_ylabel("final I_field (A)"); ax.set_title("Final interface value"); ax.grid(True, alpha=.3)

    for ax in axes.flat:
        ax.set_xlabel("run index")
        _draw_vlines(ax, vlines, colvals=colvals)
    combo = " x ".join(LABELS.get(k, k) for k in keys)
    fig.suptitle(f"Grid sweep: {combo}  ({len(ok)} runs) -- see table for parameter tuples", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    return {"sweep": _png(fig)}


def _xpos_for_value(colvals, value):
    """Map a parameter value to a fractional column position (centres at k+0.5),
    by interpolating against the per-column parameter values. None if out of range."""
    xp = np.asarray(colvals, dtype=float)
    centers = np.arange(len(xp)) + 0.5
    order = np.argsort(xp)
    xs, cs = xp[order], centers[order]
    if not (xs.min() <= value <= xs.max()):
        return None
    return float(np.interp(value, xs, cs))


def _make_iters_heatmap(dims, rows, mode, vlines=None):
    """Colormap of WR iterations per window: y = window index, x = swept point
    (one column per run), colour = iterations that window took to converge.
    Reveals which windows are hard and how difficulty shifts with the parameter."""
    ok = [r for r in rows if r.get("ok") and r.get("wr_nit")]
    if not ok:
        return {}
    ncol = len(ok)
    wmax = max(len(r["wr_nit"]) for r in ok)
    M = np.full((wmax, ncol), np.nan)
    for j, r in enumerate(ok):
        for w, v in enumerate(r["wr_nit"]):
            M[w, j] = v

    one = len(dims) == 1
    if one:
        xlabels = [f"{r['vals'][dims[0]]:g}" for r in ok]
        xlabel = dims[0]
    else:
        xlabels = [str(j) for j in range(ncol)]
        xlabel = "run index (" + " x ".join(dims) + ")"

    fig_w = max(5.0, min(14.0, 1.4 + 0.42 * ncol))
    fig_h = max(3.2, min(9.0, 1.6 + 0.30 * wmax))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    vmax = np.nanmax(M)
    pcm = ax.pcolormesh(np.arange(ncol + 1), np.arange(wmax + 1) + 0.5, M,
                        cmap="viridis", vmin=1, vmax=max(2, vmax), shading="flat")
    cbar = fig.colorbar(pcm, ax=ax)
    cbar.set_label("WR iterations")
    # Thin ticks to at most ~15 so dense sweeps stay readable.
    MAXTICKS = 15
    step = max(1, int(np.ceil(ncol / MAXTICKS)))
    tick_idx = list(range(0, ncol, step))
    dense = one and len(tick_idx) > 6
    ax.set_xticks([i + 0.5 for i in tick_idx])
    ax.set_xticklabels([xlabels[i] for i in tick_idx],
                       rotation=45 if dense else 0,
                       ha="right" if dense else "center", fontsize=9)
    # Optional user reference lines (drawn where the x-axis parameter hits the value).
    colvals = [r["vals"][dims[0]] for r in ok] if dims else None
    _draw_vlines(ax, vlines, colvals=colvals)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("window")
    ax.invert_yaxis()   # window 1 at top, time flowing downward
    ax.set_title("WR iterations per window", fontsize=12)
    fig.tight_layout()
    return {"iters": _png(fig)}


def _draw_vlines(ax, vlines, colvals=None, data_x=False):
    """Draw user reference lines. If data_x, `value` is an x-data coordinate
    (axis already in parameter units); otherwise interpolate to a column position
    using per-column `colvals`. `vlines` = list of (value, label)."""
    if not vlines:
        return
    for value, label in vlines:
        if data_x:
            xpos = value
        else:
            if not colvals:
                continue
            xpos = _xpos_for_value(colvals, value)
            if xpos is None:
                continue
        ax.axvline(xpos, color="#e63946", lw=1.6, ls="--", zorder=5)
        ax.text(xpos, 1.005, label, transform=ax.get_xaxis_transform(),
                color="#e63946", fontsize=8, rotation=90, va="bottom", ha="center")


def sweep_table(keys, rows):
    """Compact serialisable table for the UI / CLI. `keys` may be a string
    (single sweep) or a list; one value column is emitted per swept parameter."""
    if isinstance(keys, str):
        keys = [keys]
    valcols = [f"p:{k}" for k in keys]         # one column per swept parameter
    metric = ["ok", "max_WR_iterations", "mean_WR_iterations",
              "total_xyce_solves", "worst_WR_error", "all_converged",
              "solver_seconds", "final_I_field", "final_V_field"]

    def rowout(r):
        vals = r.get("vals") or {keys[0]: r.get("value")}
        d = {f"p:{k}": vals.get(k) for k in keys}
        for c in metric:
            d[c] = r.get(c)
        return d

    return {"key": keys[0], "keys": keys, "value_cols": valcols,
            "value_labels": [LABELS.get(k, k) for k in keys],
            "columns": valcols + metric,
            "rows": [rowout(r) for r in rows]}


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
                write_probes(params)
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
        if self.path not in ("/run", "/sweep", "/netlist", "/export", "/export_csv",
                             "/eval_vlines", "/sweep_replot", "/sweep_export"):
            self._send(404, json.dumps({"error": "unknown endpoint"}))
            return
        n = int(self.headers.get("Content-Length", 0))
        body = json.loads(self.rfile.read(n) or b"{}")
        if self.path == "/sweep":
            self._handle_sweep(body)
            return
        if self.path == "/eval_vlines":
            self._send(200, json.dumps({"results": eval_vlines_report(
                body.get("vlines"), body.get("params", {}))}))
            return
        if self.path == "/sweep_replot":
            # Re-render the last sweep's plots with new reference lines (no re-solve).
            s = _LAST_SWEEP
            if not s:
                self._send(200, json.dumps({"ok": False, "error": "run a sweep first"}))
                return
            vl_eval = eval_vlines(body.get("vlines"), s["params"])
            self._send(200, json.dumps({
                "ok": True,
                "plots": make_sweep_plots(s["keys"], s["rows"], s["mode"],
                                          s["free_keys"], vl_eval),
                "results": eval_vlines_report(body.get("vlines"), s["params"]),
            }))
            return
        if self.path == "/sweep_export":
            s = _LAST_SWEEP
            if not s:
                self._send(200, json.dumps({"ok": False, "error": "run a sweep first"}))
                return
            try:
                vl_eval = eval_vlines(body.get("vlines"), s["params"])
                plots = make_sweep_plots(s["keys"], s["rows"], s["mode"],
                                         s["free_keys"], vl_eval)
                csv_path, files = export_sweep_csv(body.get("path") or "results/sweeps/sweep.csv",
                                                   s["keys"], s["rows"], plots)
                self._send(200, json.dumps({"ok": True, "path": csv_path, "files": files}))
            except Exception as e:
                self._send(200, json.dumps({"ok": False, "error": str(e),
                                            "trace": traceback.format_exc()[-2000:]}))
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
                out = export_run_csv(body.get("path") or "results/results.csv",
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
            write_probes(params)
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
            mode = body.get("mode", "single")
            # Accept either the classic flat form (sweep_key/min/max/steps/scale)
            # or an explicit list of per-parameter specs.
            specs = body.get("specs")
            if not specs:
                specs = [{"key": body.get("sweep_key"), "min": body["min"],
                          "max": body["max"], "steps": body["steps"],
                          "scale": body.get("scale", "linear")}]
                mode = "single"
            for s in specs:
                if s["key"] not in SWEEPABLE:
                    raise ValueError(f"'{s['key']}' is not a sweepable parameter")
            keylist = [s["key"] for s in specs]
            if len(keylist) != len(set(keylist)):
                raise ValueError("each swept parameter must be distinct")
            keys, points, free_keys = build_sweep_points(specs, mode)
            vlines = eval_vlines(body.get("vlines"), base)
            workers = int(body.get("workers") or 0) or None
            workers_used = max(1, min(workers or (os.cpu_count() or 1), len(points) or 1))
            t0 = time.perf_counter()
            rows = run_sweep(base, keys, points, workers)
            total = round(time.perf_counter() - t0, 3)
            n_ok = sum(1 for r in rows if r.get("ok"))
            # Keep the result so reference lines can be re-drawn (Results section)
            # without re-solving. rows/keys are plain JSON-friendly data.
            global _LAST_SWEEP
            _LAST_SWEEP = {"keys": keys, "rows": rows, "mode": mode,
                           "free_keys": free_keys, "params": base}
            self._send(200, json.dumps({
                "ok": n_ok > 0,
                "sweep_key": keys[0],
                "sweep_keys": keys,
                "mode": mode,
                "n_points": len(rows),
                "n_ok": n_ok,
                "workers": workers_used,
                "sweep_seconds": total,
                "plots": make_sweep_plots(keys, rows, mode, free_keys, vlines),
                "table": sweep_table(keys, rows),
            }))
        except Exception as e:
            self._send(200, json.dumps({
                "ok": False, "error": str(e),
                "trace": traceback.format_exc()[-2000:],
            }))


# ---------------------------------------------------------------------------
# Frontend
# ---------------------------------------------------------------------------
# Properties layout: controls are split into labeled groups to cut the "wall of inputs". Each group is a
# collapsible <details class="fold"> whose body holds one or more equal-column rows (grid). The two common
# groups start open; the rest start folded (open=False). The run is absolute-end-time only ("Sim duration"
# = t_end). window_width is a UI-only derived control (not a config key): WR window width = t_end /
# N_field_windows, two-way linked with N_field_windows via linkWindows() (not collected -- no data-key).
#   group := (title, open_by_default, [rows]);  row := [control keys laid out equal-width]
_PROP_GROUPS = [
    ("Run & windows",      True,  [["t_end", "N_field_windows", "__arrow__", "window_width"]]),
    ("Coupling",           True,  [["validation_mode", "coupling_mode", "N_xyce_samples",
                                    "reconstruct_mode", "N_field_eval_intervals"]]),
    ("WR iteration",       False, [["wr_convergence_method", "WRmaxSteps", "WR_tolerance"]]),
    ("Interface & secant", False, [["use_t_floor", "t_floor_frac", "interface_form", "seam_average"]]),
    ("Field / ROM model",  False, [["R_ROM", "L_ROM", "R_FEM", "L_FEM", "nonlin_model"]]),
]
_PROP_HIDDEN = ["I_sat"]  # config-only (emitted as a hidden input so presets/reset/collect keep working)


def _controls_html():
    # Seed the form from the actual saved sim_config.txt (fall back to factory defaults),
    # so the studio opens on the current working setup, not a blank/trivial config.
    initial = {**DEFAULTS, **read_config()}
    spec = {k: (k, l, d, kind, s) for (k, l, d, kind, s) in PARAMS}

    def ctl(k):
        # Decorative "linked" arrow between N_field_windows and window_width (not a control).
        if k == "__arrow__":
            return ('<div class="link-arrow" aria-hidden="true" '
                    'title="linked: WR window width = t_end / N_field_windows">&harr;</div>')
        # Synthetic UI-only control: WR window width (s), derived from t_end / N_field_windows and
        # two-way linked with N_field_windows. No data-key -> not collected into the config.
        if k == "window_width":
            help_icon = '<span class="help" data-help="window_width">?</span>' if "window_width" in HELP else ''
            inner = ('<input type="number" id="f_window_width" step="any" '
                     'onchange="linkWindows(\'width\')">')
            return ('<div class="ctl" id="ctl_window_width">'
                    f'<label for="f_window_width">WR window width (s){help_icon}</label>'
                    f'<div class="inputs">{inner}</div></div>')
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
            # t_end / N_field_windows drive the linked window-width box live.
            extra = ""
            if k == "t_end":
                extra = "; linkWindows('tend')"
            elif k == "N_field_windows":
                extra = "; linkWindows('N')"
            inner = (f'<input type="number" id="f_{k}" step="{step}" value="{default}" '
                     f'data-key="{k}" oninput="syncFromBox(this){extra}">')
        return (f'<div class="ctl" id="ctl_{k}"><label for="f_{k}">{label}{help_icon}</label>'
                f'<div class="inputs">{inner}</div></div>')

    def row_html(row):
        cells = "\n".join(ctl(k) for k in row)
        # arrow cell is a narrow auto column; real controls share the remaining width equally
        cols = " ".join("auto" if k == "__arrow__" else "minmax(0,1fr)" for k in row)
        return f'<div class="prop-row" style="grid-template-columns:{cols}">\n{cells}\n</div>'

    out = []
    for title, is_open, rows in _PROP_GROUPS:
        body = "\n".join(row_html(r) for r in rows)
        out.append(
            f'<details class="fold propgrp"{" open" if is_open else ""}>'
            f'<summary>{title.replace("&", "&amp;")}</summary>'
            f'<div class="fold-body">\n{body}\n</div></details>'
        )
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


def _sweep_options_js():
    """JS array literal [[value,label],...] for building sweep-parameter selects."""
    return json.dumps([[k, LABELS.get(k, k)] for k in SWEEPABLE])


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
  /* one convergence-study parameter = a single no-wrap row; cells shrink instead of wrapping */
  .sw-row { display:flex; gap:8px; align-items:flex-end; flex-wrap:nowrap; margin-top:8px; }
  .sw-row .ctl { flex:1 1 0; min-width:0; }
  .sw-row .ctl.sw-rm { flex:0 0 auto; }
  .sw-row .sw-del { width:34px; height:34px; padding:0; display:inline-flex;
       align-items:center; justify-content:center; font-size:16px; line-height:1; }
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

  /* properties: labeled collapsible groups (details.fold), each body an equal-column row grid */
  .prop-rows { display:flex; flex-direction:column; gap:10px; }
  .prop-row { display:grid; gap:16px 22px; align-items:start; }
  .propgrp { margin-top:0; }                       /* group = details.fold */
  .propgrp > .fold-body { display:flex; flex-direction:column; gap:16px; }
  /* "linked" arrow between N_field_windows and WR window width; sits in a narrow auto column,
     bottom-aligned with the input boxes (row is align-items:start) */
  .link-arrow { align-self:end; display:flex; align-items:center; justify-content:center;
                height:36px; color:var(--muted); font-size:19px; user-select:none; }
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

  /* Probes panel */
  .probes-add { display:flex; gap:8px; align-items:center; margin-bottom:8px; flex-wrap:wrap; row-gap:10px; }
  .probes-add select.choice { flex:0 1 auto; width:auto; min-width:140px; max-width:220px; padding:6px 8px; }
  .probes-add .pb-lbl { flex:0 0 58px; font-size:11px; color:var(--muted);
       text-transform:uppercase; letter-spacing:.06em; }
  .probes-add .pb-div { flex:0 0 auto; width:1px; align-self:stretch; background:var(--line-strong); margin:0 6px; }
  .probes-add .pb-sep { color:var(--muted); font-size:14px; flex:0 0 auto; }
  .pb-btn { flex:0 0 auto; display:inline-flex; align-items:center; justify-content:center;
       width:34px; height:34px; padding:0; background:var(--bg); border:1px solid var(--line-strong);
       border-radius:0; cursor:pointer; }
  .pb-btn:hover { border-color:var(--fg); background:var(--field); }
  .pb-btn.pb-v { color:var(--accent); }
  .pb-btn.pb-i { color:#b5179e; }
  .probes-list { display:flex; flex-wrap:wrap; gap:6px; margin-top:10px; }
  .probes-list:empty { display:none; }
  .probes-list .chip { display:inline-flex; align-items:center; gap:6px; font-family:ui-monospace,monospace;
       font-size:12px; border:1px solid var(--line-strong); padding:3px 6px 3px 9px; }
  .probes-list .chip.vp { border-left:3px solid var(--accent); }
  .probes-list .chip.ip { border-left:3px solid #b5179e; }
  .probes-list .chip button { border:none; background:none; color:var(--muted); cursor:pointer;
       padding:0 2px; font-size:14px; line-height:1; }
  .probes-list .chip button:hover { color:var(--err); }

  /* lightbox for enlarged plots */
  #lightbox { display:none; position:fixed; inset:0; z-index:200; background:rgba(255,255,255,.94);
              align-items:center; justify-content:center; cursor:zoom-out; padding:40px; }
  #lightbox.on { display:flex; }
  #lightbox img { max-width:96vw; max-height:92vh; border:1px solid var(--line-strong); background:#fff; }

  /* PWL editor: per-line ✎ button overlaid at the end of each PWL row in the spec box */
  #specWrap { position:relative; }
  #pwlOverlay { position:absolute; inset:0; overflow:hidden; pointer-events:none; }
  #pwlOverlay button.pwl-edit { position:absolute; pointer-events:auto; width:20px; height:20px;
       padding:0; line-height:18px; font-size:12px; text-align:center; background:var(--bg);
       color:var(--accent); border:1px solid var(--line-strong); border-radius:0; cursor:pointer; }
  #pwlOverlay button.pwl-edit:hover { border-color:var(--fg); }

  /* PWL editor: modal */
  #pwlModal { display:none; position:fixed; inset:0; z-index:220; background:rgba(20,20,30,.35);
              align-items:center; justify-content:center; padding:32px; }
  #pwlModal.on { display:flex; }
  #pwlModal .panel { background:#fff; border:1px solid var(--line-strong); width:min(860px,94vw);
       max-height:92vh; display:flex; flex-direction:column; box-shadow:0 18px 50px rgba(0,0,0,.22); }
  #pwlModal .pwl-hd { display:flex; align-items:baseline; gap:12px; padding:12px 16px;
       border-bottom:1px solid var(--line); }
  #pwlModal .pwl-hd h3 { margin:0; font-size:14px; color:var(--accent); font-weight:600; }
  #pwlModal .pwl-hd .sub { font-size:12px; color:var(--muted); font-family:ui-monospace,monospace; }
  #pwlModal .pwl-body { padding:14px 16px; overflow:auto; }
  #pwlModal .pwl-tools { display:flex; flex-wrap:wrap; gap:8px; align-items:center; margin-bottom:10px; }
  #pwlModal .pwl-tools .hint { font-size:11px; color:var(--muted); margin-left:auto; }
  #pwlModal .pwl-tools label { font-size:12px; color:var(--fg); display:inline-flex;
       align-items:center; gap:5px; cursor:pointer; }
  #pwlSvg { width:100%; height:420px; border:1px solid var(--line); background:var(--field);
       touch-action:none; user-select:none; display:block; }
  #pwlSvg .grid { stroke:var(--line); stroke-width:1; }
  #pwlSvg .axis { stroke:var(--line-strong); stroke-width:1.2; }
  #pwlSvg .tick { fill:var(--muted); font-size:10px; font-family:ui-monospace,monospace; }
  #pwlSvg .albl { fill:var(--muted); font-size:11px; }
  #pwlSvg .seg { stroke:var(--accent); stroke-width:2; fill:none; }
  #pwlSvg .pt { fill:#fff; stroke:var(--accent); stroke-width:2; cursor:grab; }
  #pwlSvg .pt:hover { fill:var(--accent); }
  #pwlSvg .pt.drag { cursor:grabbing; }
  #pwlSvg .hlbl { fill:var(--accent); font-size:10px; font-family:ui-monospace,monospace; pointer-events:none; }
  #pwlModal .src-fields { display:flex; flex-wrap:wrap; gap:10px 16px; margin-top:12px; }
  #pwlModal .src-fields:empty { display:none; }
  #pwlModal .src-fields .fld { display:flex; flex-direction:column; gap:3px; }
  #pwlModal .src-fields label { font-size:11px; color:var(--muted); }
  #pwlModal .src-fields input { width:130px; font-family:ui-monospace,monospace; font-size:12px;
       background:var(--field); color:var(--fg); border:1px solid var(--line-strong); border-radius:0; padding:6px 8px; }
  #pwlModal .pwl-num { margin-top:10px; }
  #pwlModal .pwl-num summary { cursor:pointer; color:var(--muted); font-size:12px; }
  #pwlModal textarea#pwlNum { width:100%; height:64px; margin-top:6px; font-family:ui-monospace,monospace;
       font-size:12px; border:1px solid var(--line-strong); border-radius:0; padding:8px; resize:vertical; }
  #pwlModal .pwl-ft { display:flex; gap:10px; justify-content:flex-end; padding:12px 16px;
       border-top:1px solid var(--line); }
  #pwlModal .pwl-warn { color:var(--err); font-size:12px; margin-right:auto; align-self:center; }

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
        <div id="specWrap">
          <textarea id="f_circuit_spec" spellcheck="false" wrap="off">__SPEC_SEED__</textarea>
          <div id="pwlOverlay"></div>
        </div>
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
        <div id="sw_params"></div>
        <div class="btns" style="margin-top:10px">
          <button class="small" type="button" onclick="addSweepParam()">+ Add parameter</button>
        </div>
        <div class="btns" style="margin-top:16px;align-items:center;gap:10px">
          <button id="sweepBtn" onclick="runSweep()">Run sweep</button>
          <select id="sw_mode" class="choice" style="width:auto;max-width:260px" onchange="updateSweepMode()">
            <option value="parallel">Parallel (lock-step, same steps)</option>
            <option value="grid">Grid (all combinations)</option>
          </select>
          <label for="sw_jobs" style="color:var(--muted);font-size:11px;white-space:nowrap">Jobs
            <span class="help" data-help="sweep_jobs">?</span></label>
          <input type="number" id="sw_jobs" min="1" step="1" value="0" placeholder="auto"
                 style="width:70px" title="Parallel solver processes (0 = auto / all cores)">
        </div>
        <div id="sweepStatus" class="sub2" style="margin-top:10px"></div>
        <div class="note" id="sweepNote"></div>
      </div>
    </details>
    <details class="fold" style="margin-top:10px" id="probesPanel">
      <summary>Probes<span class="help" data-help="probes">?</span></summary>
      <div class="fold-body">
        <div class="probes-add">
          <span class="pb-lbl">Voltage</span>
          <select id="pb_v1" class="choice"></select>
          <span class="pb-sep">&rarr;</span>
          <select id="pb_v2" class="choice"></select>
          <button class="pb-btn pb-v" title="Add voltage probe  V(a) or V(a,b)" onclick="pbAddV()" aria-label="Add voltage probe">
            <svg viewBox="0 0 24 24" width="17" height="17" fill="none" stroke="currentColor" stroke-width="2"
                 stroke-linecap="round" stroke-linejoin="round"><path d="M2 12c3-8 5-8 7 0s4 8 7 0"/><path d="M20 8v8M16 12h8" stroke-width="2"/></svg>
          </button>
        </div>
        <div class="probes-add">
          <span class="pb-lbl">Current</span>
          <select id="pb_i" class="choice"></select>
          <button class="pb-btn pb-i" title="Add current probe  I(element)" onclick="pbAddI()" aria-label="Add current probe">
            <svg viewBox="0 0 24 24" width="17" height="17" fill="none" stroke="currentColor" stroke-width="2"
                 stroke-linecap="round" stroke-linejoin="round"><path d="M2 12h13M11 7l5 5-5 5"/><path d="M20 8v8M16 12h8" stroke-width="2"/></svg>
          </button>
        </div>
        <div id="pbList" class="probes-list"></div>
        <div class="note">Extra output points plotted after a run (Probe voltages / currents).
          Targets come from your circuit spec above.</div>
      </div>
    </details>
  </section>

  <!-- 4. Results -->
  <section style="border-bottom:none">
    <div class="sec-hd">
      <h2>Results</h2>
      <div class="hd-tools">
        <span class="mini" id="csvStatus"></span>
        <input type="text" id="csvPath" class="pathin" placeholder="results/results.csv" value="results/results.csv" spellcheck="false">
        <button class="small" onclick="exportCsv()">Export CSV</button>
      </div>
    </div>
    <div class="summary" id="summary"></div>
    <div class="plots">
      <img class="plot" id="p_voltage" onclick="enlarge(this)">
      <img class="plot" id="p_current" onclick="enlarge(this)">
      <img class="plot" id="p_wr" onclick="enlarge(this)">
      <img class="plot" id="p_probe_v" onclick="enlarge(this)">
      <img class="plot" id="p_probe_i" onclick="enlarge(this)">
      <img class="plot wide" id="p_sweep" onclick="enlarge(this)">
      <img class="plot wide" id="p_sweep_iters" onclick="enlarge(this)">
    </div>
    <div id="vlineBar" style="display:none">
      <div class="btns" style="margin-top:12px;align-items:center;gap:10px">
        <label for="sw_vlines" style="color:var(--muted);font-size:11px;white-space:nowrap">Reference line(s)
          <span class="help" data-help="sweep_vlines">?</span></label>
        <input type="text" id="sw_vlines" spellcheck="false" placeholder="e.g. Rs+Ls, R_ROM, 0.01"
               style="flex:1;min-width:200px;max-width:360px;font-family:ui-monospace,monospace"
               oninput="vlinePreview()">
        <span id="sw_vlines_val" class="sub2" style="font-family:ui-monospace,monospace;white-space:nowrap"></span>
      </div>
      <div class="btns" style="margin-top:10px;align-items:center;gap:10px">
        <input type="text" id="sweepCsvPath" class="pathin" placeholder="results/sweeps/sweep.csv" value="results/sweeps/sweep.csv" spellcheck="false">
        <button class="small" onclick="exportSweep()">Export sweep<span class="help" data-help="sweep_export">?</span></button>
        <span class="mini" id="sweepCsvStatus"></span>
      </div>
    </div>
    <div id="sweepTable" style="margin-top:16px"></div>
    <details class="plain"><summary>Solver log</summary><pre id="log"></pre></details>
  </section>

</main>
<div id="lightbox" onclick="this.classList.remove('on')"><img id="lightboxImg"></div>
<div id="pwlModal">
  <div class="panel">
    <div class="pwl-hd">
      <h3 id="pwlTitle">Source editor</h3>
      <span class="sub" id="pwlName"></span>
    </div>
    <div class="pwl-body">
      <div class="pwl-tools">
        <button class="small" onclick="pwlFit()">Fit</button>
        <button class="small" id="pwlClearBtn" onclick="pwlClear()">Clear</button>
        <label><input type="checkbox" id="pwlSnap"> snap</label>
        <span class="hint" id="pwlHint"></span>
      </div>
      <svg id="pwlSvg" xmlns="http://www.w3.org/2000/svg"></svg>
      <div id="srcFields" class="src-fields"></div>
      <details class="pwl-num" id="pwlNumBox">
        <summary>Numeric points (t v pairs)</summary>
        <textarea id="pwlNum" spellcheck="false" oninput="pwlNumEdit()"></textarea>
      </details>
    </div>
    <div class="pwl-ft">
      <span class="pwl-warn" id="pwlWarn"></span>
      <button onclick="pwlCancel()">Cancel</button>
      <button class="primary" onclick="pwlSave()">Save</button>
    </div>
  </div>
</div>
<div id="helpTip"></div>
<script>
const DEFAULTS = __DEFAULTS__;
const PRESETS = __PRESETS__;
const VISIBLE_WHEN = __VISIBILITY__;
const HELP = __HELP__;
const SWEEP_OPTS = __SWEEP_OPTS_JS__;

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
  linkWindows('N');
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
// WR window discretization: window width = t_end / N_field_windows, two-way linked.
//   src 'width' -> derive N_field_windows = round(t_end/width) (>=1), then snap the width box to the
//                  true t_end/N (width may not divide t_end evenly).
//   src 'tend'/'N' (or init) -> just refresh the width box from the current t_end and N.
function linkWindows(src){
  const teEl=document.getElementById('f_t_end');
  const nEl=document.getElementById('f_N_field_windows');
  const wEl=document.getElementById('f_window_width');
  if(!teEl||!nEl||!wEl) return;
  const te=parseFloat(teEl.value);
  if(src==='width'){
    const w=parseFloat(wEl.value);
    if(isFinite(te)&&te>0&&isFinite(w)&&w>0) nEl.value=Math.max(1, Math.round(te/w));
  }
  let n=Math.round(parseFloat(nEl.value));
  if(!isFinite(n)||n<1) n=1;
  nEl.value=n;                                  // keep the count integer/valid
  if(isFinite(te)&&te>0) wEl.value=(te/n).toExponential(3);
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
  if (ta) p['circuit_spec'] = ta.value;   // custom node-graph spec (authored circuit side)
  p['probes'] = probes.slice();           // user output probes -> probes.txt
  return p;
}
function resetDefaults(){
  for (const k in DEFAULTS){
    const b = document.getElementById('f_'+k);
    if (b){ b.value = DEFAULTS[k]; if (b.tagName !== 'SELECT') syncFromBox(b); }
  }
  setStatus('Reset to defaults.', '');
  applyVisibility();
  linkWindows('N');
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
      setPlot('p_probe_v', j.plots.probe_v);
      setPlot('p_probe_i', j.plots.probe_i);
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
  const valcols = tbl.value_cols || ['value'];
  const vlabels = tbl.value_labels || valcols;
  const head = ['#'].concat(vlabels, ['iters(max)','iters(mean)','Xyce','worstErr','conv','sec','I_field','V_field']);
  const cols = valcols.concat(['max_WR_iterations','mean_WR_iterations','total_xyce_solves',
    'worst_WR_error','all_converged','solver_seconds','final_I_field','final_V_field']);
  let html = '<table class="sweep"><thead><tr><th>'+head.join('</th><th>')+'</th></tr></thead><tbody>';
  tbl.rows.forEach((r, i) => {
    const bad = (r.ok === false) || (r.all_converged === false);
    html += '<tr'+(bad?' class="bad"':'')+'>';
    html += '<td>'+i+'</td>';
    for (const c of cols) html += '<td>'+fmtCell(r[c])+'</td>';
    html += '</tr>';
  });
  html += '</tbody></table>';
  host.innerHTML = html;
}

// ---- dynamic N-parameter sweep rows ----
let sweepRowSeq = 0;
function sweepOptionsHtml(sel){
  return SWEEP_OPTS.map(o => '<option value="'+o[0]+'"'+(o[0]===sel?' selected':'')+'>'+o[1]+'</option>').join('');
}
function firstUnusedParam(){
  const used = Array.from(document.querySelectorAll('#sw_params .sw-key')).map(s => s.value);
  for (const o of SWEEP_OPTS) if (used.indexOf(o[0]) < 0) return o[0];
  return SWEEP_OPTS[0][0];
}
function addSweepParam(key){
  const host = document.getElementById('sw_params');
  const id = ++sweepRowSeq;
  const sel = key || firstUnusedParam();
  const row = document.createElement('div');
  row.className = 'sw-row';
  row.dataset.id = id;
  const baseDefault = SWEEP_OPTS.find(o => o[0] !== sel);
  row.innerHTML =
    '<div class="ctl"><label>Parameter</label><div class="inputs">'+
      '<select class="choice sw-key">'+sweepOptionsHtml(sel)+'</select></div></div>'+
    '<div class="ctl"><label>Type</label><div class="inputs">'+
      '<select class="choice sw-type" onchange="updateSweepMode()">'+
      '<option value="range">Range</option><option value="link">% of&hellip;</option></select></div></div>'+
    // range cells
    '<div class="ctl sw-range"><label>Min</label><div class="inputs">'+
      '<input type="number" class="sw-min" step="any" value="0"></div></div>'+
    '<div class="ctl sw-range"><label>Max</label><div class="inputs">'+
      '<input type="number" class="sw-max" step="any" value="0.1"></div></div>'+
    '<div class="ctl sw-range sw-steps-ctl"><label>Steps</label><div class="inputs">'+
      '<input type="number" class="sw-steps" step="1" value="8" oninput="updateSweepMode()"></div></div>'+
    '<div class="ctl sw-range"><label>Spacing</label><div class="inputs">'+
      '<select class="choice sw-scale"><option value="linear">linear</option>'+
      '<option value="log">log (positive only)</option></select></div></div>'+
    // linked cells (percentage of another swept parameter)
    '<div class="ctl sw-link" style="display:none"><label>Percent</label><div class="inputs">'+
      '<input type="number" class="sw-pct" step="any" value="50"></div></div>'+
    '<div class="ctl sw-link" style="display:none"><label>of</label><div class="inputs">'+
      '<select class="choice sw-base">'+sweepOptionsHtml(baseDefault?baseDefault[0]:sel)+'</select></div></div>'+
    '<div class="ctl sw-rm"><label>&nbsp;</label><div class="inputs">'+
      '<button type="button" class="small sw-del" title="Remove parameter" '+
      'onclick="removeSweepParam(this)">&times;</button></div></div>';
  host.appendChild(row);
  updateSweepMode();
}
function removeSweepParam(btn){
  const rows = document.querySelectorAll('#sw_params .sw-row');
  if (rows.length <= 1) return;   // keep at least one parameter
  btn.closest('.sw-row').remove();
  updateSweepMode();
}
function rowIsLink(row){ return row.querySelector('.sw-type').value === 'link'; }
function updateSweepMode(){
  const mode = document.getElementById('sw_mode').value;
  const rows = Array.from(document.querySelectorAll('#sw_params .sw-row'));
  const rangeRows = rows.filter(r => !rowIsLink(r));
  const firstRange = rangeRows[0] || null;
  const step0 = firstRange ? firstRange.querySelector('.sw-steps').value : '';
  const freeKeys = rangeRows.map(r => r.querySelector('.sw-key').value);
  let nLink = 0;
  rows.forEach((row, i) => {
    const linked = rowIsLink(row);
    if (linked) nLink++;
    // Show only the cells for this row's type.
    row.querySelectorAll('.sw-range').forEach(c => c.style.display = linked ? 'none' : '');
    row.querySelectorAll('.sw-link').forEach(c => c.style.display = linked ? '' : 'none');
    // Column labels only on the first row; the rest align underneath it.
    row.querySelectorAll('label').forEach(l => l.style.display = (i === 0) ? '' : 'none');
    if (!linked){
      // Parallel locks every free parameter to the first range row's step count.
      const inp = row.querySelector('.sw-steps');
      const lock = (mode === 'parallel' && row !== firstRange);
      inp.disabled = lock;
      if (lock) inp.value = step0;
    } else {
      // Base dropdown lists the independently swept (range) parameters, minus self.
      const bsel = row.querySelector('.sw-base'), cur = bsel.value;
      const selfKey = row.querySelector('.sw-key').value;
      const opts = freeKeys.filter(k => k !== selfKey);
      bsel.innerHTML = opts.length
        ? opts.map(k => { const o = SWEEP_OPTS.find(x => x[0] === k);
            return '<option value="'+k+'"'+(k===cur?' selected':'')+'>'+(o?o[1]:k)+'</option>'; }).join('')
        : '<option value="">(add a range parameter)</option>';
    }
    row.querySelector('.sw-del').style.visibility = (rows.length > 1) ? 'visible' : 'hidden';
  });
  const nFree = rangeRows.length;
  const note = document.getElementById('sweepNote');
  const linkNote = nLink ? '  Linked (% of) parameters track their base and add no dimension.' : '';
  if (nFree <= 1)
    note.textContent = 'Single parameter: holds all other fields fixed and runs the solver once per swept value.' + linkNote;
  else if (mode === 'parallel')
    note.textContent = 'Parallel: all range parameters advance together over the same number of steps (Min→Max), one solver run per step.' + linkNote;
  else
    note.textContent = 'Grid: runs every combination of the range parameters (step counts may differ). 2 range params -> heatmaps; more -> metric-vs-run-index (see table for tuples).' + linkNote;
}
function collectSweepSpecs(mode){
  const rows = Array.from(document.querySelectorAll('#sw_params .sw-row'));
  const rangeRows = rows.filter(r => !rowIsLink(r));
  const firstRange = rangeRows[0] || null;
  const step0 = firstRange ? parseInt(firstRange.querySelector('.sw-steps').value, 10) : 8;
  return rows.map(row => {
    const key = row.querySelector('.sw-key').value;
    if (rowIsLink(row))
      return { key, link: { base: row.querySelector('.sw-base').value,
                            pct: parseFloat(row.querySelector('.sw-pct').value) } };
    return {
      key,
      min: parseFloat(row.querySelector('.sw-min').value),
      max: parseFloat(row.querySelector('.sw-max').value),
      steps: (mode === 'parallel' && row !== firstRange) ? step0
             : parseInt(row.querySelector('.sw-steps').value, 10),
      scale: row.querySelector('.sw-scale').value,
    };
  });
}
async function exportSweep(){
  const st = document.getElementById('sweepCsvStatus');
  const path = (document.getElementById('sweepCsvPath').value || '').trim();
  if (!path){ st.textContent = 'Enter a path.'; return; }
  const vlines = (document.getElementById('sw_vlines').value || '').trim();
  st.textContent = 'Saving...';
  try {
    const res = await fetch('/sweep_export', {method:'POST', headers:{'Content-Type':'application/json'},
      body: JSON.stringify({path: path, vlines: vlines})});
    const j = await res.json();
    st.textContent = j.ok ? ('Saved '+(j.files ? j.files.length : 1)+' file(s): '+j.path)
                          : ('Error: '+(j.error||'unknown'));
  } catch(e){ st.textContent = 'Failed: '+e; }
}

let vlineDeb = null;
function showVlineVals(results){
  const out = document.getElementById('sw_vlines_val');
  out.innerHTML = (results && results.length) ? ('= ' + results.map(r =>
    r.ok ? '<span style="color:var(--ok)">'+fmtCell(r.value)+'</span>'
         : '<span style="color:var(--err)">'+r.error+'</span>').join(', ')) : '';
}
function vlinePreview(){
  clearTimeout(vlineDeb);
  vlineDeb = setTimeout(async () => {
    const raw = (document.getElementById('sw_vlines').value || '').trim();
    try {
      // Re-draw the stored sweep's plots with the new reference lines (no re-solve).
      const res = await fetch('/sweep_replot', {method:'POST', headers:{'Content-Type':'application/json'},
        body: JSON.stringify({vlines: raw})});
      const j = await res.json();
      if (j.ok){
        setPlot('p_sweep', j.plots ? j.plots.sweep : null);
        setPlot('p_sweep_iters', j.plots ? j.plots.iters : null);
        showVlineVals(j.results);
      }
    } catch(e){ /* leave plots as-is */ }
  }, 200);
}
async function runSweep(){
  const btn = document.getElementById('sweepBtn'); btn.disabled = true;
  const mode = document.getElementById('sw_mode').value;
  const specs = collectSweepSpecs(mode);
  const keys = specs.map(s => s.key);
  if (keys.length !== new Set(keys).size){
    setSweepStatus('Each swept parameter must be distinct.', 'err');
    btn.disabled = false; return;
  }
  const label = keys.join((mode === 'grid' && keys.length > 1) ? ' x ' : ' + ');
  const jobs = parseInt(document.getElementById('sw_jobs').value, 10) || 0;
  const vlines = (document.getElementById('sw_vlines').value || '').trim();
  setSweepStatus('Running sweep of '+label+'...', '');
  try {
    const payload = { params: collect(), mode: mode, specs: specs, workers: jobs, vlines: vlines };
    const res = await fetch('/sweep', {method:'POST', headers:{'Content-Type':'application/json'},
                                       body: JSON.stringify(payload)});
    const j = await res.json();
    if (!j.ok){ setSweepStatus('Sweep error: '+(j.error||'all points failed'), 'err'); }
    else {
      const jw = j.workers ? (' on '+j.workers+' job'+(j.workers>1?'s':'')) : '';
      setSweepStatus('Sweep done: '+j.n_ok+'/'+j.n_points+' points'+jw+' in '+j.sweep_seconds.toFixed(1)+' s.', 'ok');
      setPlot('p_sweep', j.plots ? j.plots.sweep : null);
      setPlot('p_sweep_iters', j.plots ? j.plots.iters : null);
      showSweepTable(j.table);
      document.getElementById('vlineBar').style.display = '';   // enable post-processing lines
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
      pwlRenderButtons();  // refresh the ✎ PWL edit buttons in the spec box
      pbRefresh();         // refresh probe target dropdowns from the (possibly preset-loaded) spec
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

// ================= Probes ===================================================
// User output probes: V(node) or I(element). Targets parsed client-side from the spec
// textarea. Stored as Xyce print tokens in `probes` and shipped to probes.txt via collect().
let probes = [];
function pbSpecElements(){
  const ta=document.getElementById('f_circuit_spec'); if(!ta) return [];
  const TYPES={R:1,L:1,C:1,VSIN:1,ISIN:1,VDC:1,IDC:1,VPULSE:1,IPULSE:1,VPWM:1,IPWM:1,VPWL:1,IPWL:1,SW:1};
  const out=[];
  for(const raw of ta.value.split('\\n')){
    const s=raw.trim(); if(!s || s[0]==='#') continue;
    const t=s.split(/[\\s,]+/), ty=t[0].toUpperCase();
    if(TYPES[ty]) out.push({type:ty, name:t[1]||'', a:t[2]||'', b:t[3]||''});
  }
  return out;
}
function pbNodes(){
  const list=[], seen={};
  for(const e of pbSpecElements()) for(const nd of [e.a,e.b]){
    if(!nd || nd==='0' || nd.toLowerCase()==='nx') continue;
    if(!seen[nd]){ seen[nd]=1; list.push(nd); }
  }
  return list;
}
function pbSetOptions(sel, opts){   // opts: [[value,label],...]; preserve current selection if still present
  const prev=sel.value;
  sel.innerHTML=opts.map(o=>'<option value="'+o[0]+'">'+pwlEsc(o[1])+'</option>').join('');
  if(prev){ for(const o of sel.options) if(o.value===prev){ sel.value=prev; break; } }
}
function pbFillTargets(){
  const nodes=pbNodes();
  const nodeOpts=nodes.map(n=>[n,n]);
  pbSetOptions(document.getElementById('pb_v1'), nodeOpts);                       // V: node
  pbSetOptions(document.getElementById('pb_v2'), [['0','gnd (0)']].concat(nodeOpts)); // V: reference
  const els=[];
  for(const e of pbSpecElements()) if(e.type==='R'||e.type==='L'||e.type==='C')
    els.push(['I('+e.type+e.name+')', e.name+' ('+e.type+')']);
  pbSetOptions(document.getElementById('pb_i'), els);
}
function pbPush(tok){
  if(probes.indexOf(tok)>=0){ setStatus('Probe '+tok+' already added.',''); return; }
  probes.push(tok); pbRenderList(); loadCircuit();   // re-emit so the .print line updates live
}
function pbAddV(){
  const n1=document.getElementById('pb_v1').value; if(!n1) return;
  const n2=document.getElementById('pb_v2').value;
  pbPush((n2 && n2!=='0' && n2!==n1) ? 'V('+n1+','+n2+')' : 'V('+n1+')');   // differential vs ground
}
function pbAddI(){
  const tok=document.getElementById('pb_i').value; if(tok) pbPush(tok);
}
function pbRemove(i){ probes.splice(i,1); pbRenderList(); loadCircuit(); }
function pbRenderList(){
  const host=document.getElementById('pbList'); if(!host) return;
  let html='';
  probes.forEach((tk,i)=>{ const cls=tk.toUpperCase().indexOf('I(')===0?'ip':'vp';
    html+='<span class="chip '+cls+'">'+pwlEsc(tk)+'<button title="remove" onclick="pbRemove('+i+')">×</button></span>'; });
  host.innerHTML=html;
}
function pbRefresh(){ pbFillTargets(); pbRenderList(); }

// ================= Source graph editor =====================================
// Fully client-side: scan the spec textarea for source lines (SIN/DC/PULSE/PWL),
// overlay an ✎ edit button on each row, and open a modal SVG editor. PWL is a
// draggable point editor; SIN/DC/PULSE show the analytic waveform with draggable
// parameter handles + synced numeric fields. Writes back into the one spec line.
let pwlState = null, pwlDeb = null;
const PWL_PAD = {l:56, r:16, t:16, b:38};

function pwlEsc(s){ return String(s).replace(/[&<>"]/g, c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c])); }
function srcNum(x,d){ return isFinite(x)?x:d; }

function pwlFmt(x){
  if(!isFinite(x)) return '0';
  if(x===0) return '0';
  const a=Math.abs(x);
  let s=(a<1e-3||a>=1e5) ? x.toExponential(4) : x.toPrecision(6);
  if(s.indexOf('e')<0 && s.indexOf('.')>=0) s=s.replace(/0+$/,'').replace(/\\.$/,'');
  return s;
}

// parse a spec source line -> {kind:'V'|'I', type:'SIN'|'DC'|'PULSE'|'PWL', name,a,b, (pts|params)}
function srcParseLine(line){
  const s=(line||'').trim();
  if(!s || s[0]==='#') return null;
  const tok=s.split(/[\\s,]+/);
  const m=tok[0].toUpperCase().match(/^([VI])(SIN|DC|PULSE|PWM|PWL)$/);
  if(!m) return null;
  const base={kind:m[1], type:m[2], name:tok[1]||'', a:tok[2]||'', b:tok[3]||''};
  const nums=tok.slice(4).map(parseFloat);
  if(base.type==='PWL'){
    const g=nums.filter(x=>!isNaN(x)), pts=[];
    for(let i=0;i+1<g.length;i+=2) pts.push([g[i], g[i+1]]);
    base.pts=pts;
  } else if(base.type==='SIN'){ base.params={amp:srcNum(nums[0],1), freq:srcNum(nums[1],50)}; }
  else if(base.type==='DC'){ base.params={value:srcNum(nums[0],1)}; }
  else if(base.type==='PULSE'){ base.params={v1:srcNum(nums[0],0), v2:srcNum(nums[1],1), td:srcNum(nums[2],0), tr:srcNum(nums[3],1e-4)}; }
  else if(base.type==='PWM'){ base.params={v1:srcNum(nums[0],0), v2:srcNum(nums[1],1), freq:srcNum(nums[2],1e3), duty:srcNum(nums[3],0.5)}; }
  return base;
}
function pwlScan(){
  const ta=document.getElementById('f_circuit_spec');
  if(!ta) return [];
  const out=[];
  ta.value.split('\\n').forEach((ln,i)=>{ const p=srcParseLine(ln); if(p){ p.lineIdx=i; out.push(p); } });
  return out;
}
// Overlay one ✎ button at the end (right edge) of each source row in the spec
// textarea. wrap="off" => one spec line is exactly one visual row, so line index
// maps to a row via line-height; re-run on input/scroll/resize.
function pwlRenderButtons(){
  const ta=document.getElementById('f_circuit_spec'), ov=document.getElementById('pwlOverlay');
  if(!ta || !ov) return;
  const cs=getComputedStyle(ta);
  let lh=parseFloat(cs.lineHeight); if(!isFinite(lh)) lh=parseFloat(cs.fontSize)*1.6;
  const padTop=parseFloat(cs.paddingTop)||0, bh=20;
  const right=(ta.scrollHeight>ta.clientHeight)?16:7;   // clear the vertical scrollbar
  let html='';
  for(const s of pwlScan()){
    const top=padTop + s.lineIdx*lh - ta.scrollTop + (lh-bh)/2;
    html += '<button class="pwl-edit" style="top:'+top+'px;right:'+right+'px" '
          + 'title="edit '+s.kind+s.type+' '+pwlEsc(s.name)+'" onclick="pwlOpen('+s.lineIdx+')">✎</button>';
  }
  ov.innerHTML=html;
}

// ---- waveform model (SIN/DC/PULSE) ----
function srcWave(type,p,xmin,xmax,n){
  if(type==='DC') return [[xmin,p.value],[xmax,p.value]];
  if(type==='SIN'){ const out=[]; for(let i=0;i<=n;i++){ const t=xmin+(xmax-xmin)*i/n; out.push([t, p.amp*Math.sin(2*Math.PI*p.freq*t)]); } return out; }
  if(type==='PULSE'){
    const tr=p.tr||1e-30;
    const f=t=> t<p.td ? p.v1 : (t<p.td+tr ? p.v1+(p.v2-p.v1)*(t-p.td)/tr : p.v2);
    const bk=[xmin,xmax,p.td,p.td+tr].filter(t=>t>=xmin&&t<=xmax).sort((a,b)=>a-b);
    const uniq=[]; for(const t of bk){ if(!uniq.length||Math.abs(t-uniq[uniq.length-1])>1e-15) uniq.push(t); }
    return uniq.map(t=>[t,f(t)]);
  }
  if(type==='PWM'){
    const per=p.freq>0?1/p.freq:(xmax-xmin);
    const duty=Math.min(1,Math.max(0,p.duty));
    const val=t=>{ const ph=t/per-Math.floor(t/per); return ph<duty ? p.v2 : p.v1; };  // high at period start
    const edges=[];                              // rising (k*per) and falling (k*per+duty*per) edges
    for(let s=Math.floor(xmin/per)*per; s<=xmax+per && edges.length<6000; s+=per){ edges.push(s); edges.push(s+duty*per); }
    const trans=edges.filter(t=>t>xmin && t<xmax).sort((a,b)=>a-b);
    const out=[[xmin, val(xmin+per*1e-9)]];
    for(const t of trans){ out.push([t, val(t-per*1e-9)]); out.push([t, val(t+per*1e-9)]); }  // vertical edge
    out.push([xmax, val(xmax-per*1e-9)]);
    return out;
  }
  return [];
}
function srcRange(type,p){
  if(type==='DC') return [p.value,p.value];
  if(type==='SIN'){ const a=Math.abs(p.amp); return [-a,a]; }
  if(type==='PULSE'||type==='PWM') return [Math.min(p.v1,p.v2), Math.max(p.v1,p.v2)];
  return [-1,1];
}
function srcHandles(){
  const p=pwlState.params, v=pwlState.view, type=pwlState.type;
  const cl=t=>Math.min(v.xmax, Math.max(v.xmin, t));
  if(type==='DC') return [{id:'value', t:(v.xmin+v.xmax)/2, v:p.value, label:'value'}];
  if(type==='SIN') return [
    {id:'amp',  t:cl(p.freq>0?1/(4*p.freq):(v.xmin+v.xmax)/2), v:p.amp, label:'amp'},
    {id:'freq', t:cl(p.freq>0?1/p.freq:(v.xmax-v.xmin)),       v:0,     label:'freq (1 period)'} ];
  if(type==='PULSE') return [
    {id:'v1', t:cl((v.xmin+p.td)/2),          v:p.v1, label:'v1'},
    {id:'v2', t:cl((p.td+p.tr+v.xmax)/2),     v:p.v2, label:'v2'},
    {id:'td', t:p.td,                         v:p.v1, label:'td'},
    {id:'tr', t:p.td+p.tr,                    v:p.v2, label:'tr'} ];
  if(type==='PWM'){ const per=p.freq>0?1/p.freq:(v.xmax-v.xmin), duty=Math.min(1,Math.max(0,p.duty));
    return [
    {id:'v2',   t:cl(duty*per/2),        v:p.v2, label:'high'},
    {id:'v1',   t:cl(duty*per+(1-duty)*per/2), v:p.v1, label:'low'},
    {id:'freq', t:cl(per),               v:(p.v1+p.v2)/2, label:'freq (1 period)'},
    {id:'duty', t:cl(duty*per),          v:p.v2, label:'duty'} ]; }
  return [];
}
function srcApplyHandle(id,t,y){
  const p=pwlState.params;
  if(id==='value') p.value=y;
  else if(id==='amp') p.amp=y;
  else if(id==='freq'){ if(t>1e-15) p.freq=1/t; }
  else if(id==='v1') p.v1=y;
  else if(id==='v2') p.v2=y;
  else if(id==='td') p.td=Math.max(0,t);
  else if(id==='tr') p.tr=Math.max(1e-15, t-p.td);
  else if(id==='duty'){ const per=p.freq>0?1/p.freq:1; p.duty=Math.min(1,Math.max(0, t/per)); }
}
function srcRenderFields(){
  const host=document.getElementById('srcFields'); if(!host) return;
  if(!pwlState || pwlState.type==='PWL'){ host.innerHTML=''; return; }
  const unit=pwlState.kind==='V'?'V':'A';
  const defs={
    SIN:[['amp','amplitude ('+unit+')'],['freq','frequency (Hz)']],
    DC:[['value','value ('+unit+')']],
    PULSE:[['v1','v1 ('+unit+')'],['v2','v2 ('+unit+')'],['td','delay td (s)'],['tr','rise tr (s)']],
    PWM:[['v1','low ('+unit+')'],['v2','high ('+unit+')'],['freq','frequency (Hz)'],['duty','duty (0-1)']]
  }[pwlState.type];
  const p=pwlState.params;
  let html='';
  for(const kv of defs)
    html+='<div class="fld"><label>'+kv[1]+'</label><input type="text" data-k="'+kv[0]+'" value="'+pwlFmt(p[kv[0]])+'" oninput="srcFieldEdit()"></div>';
  host.innerHTML=html;
}
function srcFieldEdit(){   // typed field -> params (does NOT rebuild fields, to keep caret)
  if(!pwlState || pwlState.type==='PWL') return;
  document.querySelectorAll('#srcFields input').forEach(inp=>{
    const val=parseFloat(inp.value);
    if(isFinite(val)) pwlState.params[inp.getAttribute('data-k')]=val;
  });
  pwlDraw();
}

// ---- geometry ----
function pwlDims(){
  const svg=document.getElementById('pwlSvg');
  const W=svg.clientWidth||800, H=svg.clientHeight||420;
  return {W,H, x0:PWL_PAD.l, y0:H-PWL_PAD.b, x1:W-PWL_PAD.r, y1:PWL_PAD.t};
}
function pwlWX(t){ const d=pwlDims(), v=pwlState.view; return d.x0+(t-v.xmin)/(v.xmax-v.xmin)*(d.x1-d.x0); }
function pwlWY(y){ const d=pwlDims(), v=pwlState.view; return d.y0+(y-v.ymin)/(v.ymax-v.ymin)*(d.y1-d.y0); }
function pwlPX(px){ const d=pwlDims(), v=pwlState.view; return v.xmin+(px-d.x0)/(d.x1-d.x0)*(v.xmax-v.xmin); }
function pwlPY(py){ const d=pwlDims(), v=pwlState.view; return v.ymin+(py-d.y0)/(d.y1-d.y0)*(v.ymax-v.ymin); }
function pwlStep(min,max,n){
  const span=max-min; if(!(span>0)) return 1;
  const raw=span/n, mag=Math.pow(10,Math.floor(Math.log10(raw))), norm=raw/mag;
  return (norm<1.5?1:norm<3?2:norm<7?5:10)*mag;
}
function pwlTicks(min,max,n){
  const step=pwlStep(min,max,n), out=[];
  for(let t=Math.ceil(min/step)*step; t<=max+step*0.5; t+=step) out.push(Math.abs(t)<step*1e-6?0:t);
  return out;
}
function pwlSnap(v,min,max){ const s=pwlStep(min,max,6); return Math.round(v/s)*s; }

// ---- view controls ----
function pwlResetView(){
  if(!pwlState) return;
  const teEl=document.getElementById('f_t_end');
  const te=teEl?parseFloat(teEl.value):NaN;
  let xmin=0, xmax=(isFinite(te)&&te>0)?te:1, ymin, ymax;
  if(pwlState.type==='PWL'){
    const xs=pwlState.pts.map(p=>p[0]), ys=pwlState.pts.map(p=>p[1]);
    if(xs.length){ xmin=Math.min(xmin,...xs); xmax=Math.max(xmax,...xs); }
    ymin=ys.length?Math.min(...ys):-1; ymax=ys.length?Math.max(...ys):1;
  } else {
    const p=pwlState.params;
    if(pwlState.type==='PULSE') xmax=Math.max(xmax, (p.td+p.tr)*1.3);
    if(pwlState.type==='PWM' && p.freq>0){ const per=1/p.freq; if(xmax/per>12) xmax=per*8; }  // show ~8 cycles
    const r=srcRange(pwlState.type,p); ymin=r[0]; ymax=r[1];
  }
  if(xmax<=xmin) xmax=xmin+1;
  if(ymax<=ymin){ ymin-=1; ymax+=1; }
  const py=(ymax-ymin)*0.15;
  pwlState.view={xmin, xmax, ymin:ymin-py, ymax:ymax+py};
  pwlDraw();
}
function pwlFit(){
  if(!pwlState) return;
  if(pwlState.type!=='PWL' || !pwlState.pts.length){ pwlResetView(); return; }
  const xs=pwlState.pts.map(p=>p[0]), ys=pwlState.pts.map(p=>p[1]);
  let xmin=Math.min(...xs), xmax=Math.max(...xs), ymin=Math.min(...ys), ymax=Math.max(...ys);
  if(xmax<=xmin) xmax=xmin+1;
  if(ymax<=ymin){ ymin-=1; ymax+=1; }
  pwlState.view={xmin:xmin-(xmax-xmin)*0.08, xmax:xmax+(xmax-xmin)*0.08,
                 ymin:ymin-(ymax-ymin)*0.12, ymax:ymax+(ymax-ymin)*0.12};
  pwlDraw();
}
function pwlZoomAt(f, c){
  const v=pwlState.view;
  const cx=c?c.tx:(v.xmin+v.xmax)/2, cy=c?c.ty:(v.ymin+v.ymax)/2;
  v.xmin=cx-(cx-v.xmin)/f; v.xmax=cx+(v.xmax-cx)/f;
  v.ymin=cy-(cy-v.ymin)/f; v.ymax=cy+(v.ymax-cy)/f;
  pwlDraw();
}
function pwlClear(){ if(pwlState && pwlState.type==='PWL'){ pwlState.pts=[[0,0]]; pwlDraw(); } }

// ---- draw ----
function pwlDraw(){
  if(!pwlState) return;
  const svg=document.getElementById('pwlSvg'), d=pwlDims(), v=pwlState.view;
  const unit=pwlState.kind==='V'?'V':'A';
  let s='';
  for(const t of pwlTicks(v.xmin,v.xmax,7)){ const x=pwlWX(t);
    s+='<line class="grid" x1="'+x+'" y1="'+d.y1+'" x2="'+x+'" y2="'+d.y0+'"/>';
    s+='<text class="tick" x="'+x+'" y="'+(d.y0+15)+'" text-anchor="middle">'+pwlFmt(t)+'</text>'; }
  for(const y of pwlTicks(v.ymin,v.ymax,6)){ const py=pwlWY(y);
    s+='<line class="grid" x1="'+d.x0+'" y1="'+py+'" x2="'+d.x1+'" y2="'+py+'"/>';
    s+='<text class="tick" x="'+(d.x0-7)+'" y="'+(py+3)+'" text-anchor="end">'+pwlFmt(y)+'</text>'; }
  s+='<line class="axis" x1="'+d.x0+'" y1="'+d.y1+'" x2="'+d.x0+'" y2="'+d.y0+'"/>';
  s+='<line class="axis" x1="'+d.x0+'" y1="'+d.y0+'" x2="'+d.x1+'" y2="'+d.y0+'"/>';
  s+='<text class="albl" x="'+((d.x0+d.x1)/2)+'" y="'+(d.y0+32)+'" text-anchor="middle">time (s)</text>';
  const yc=(d.y0+d.y1)/2, xl=14;
  s+='<text class="albl" x="'+xl+'" y="'+yc+'" text-anchor="middle" transform="rotate(-90 '+xl+' '+yc+')">level ('+unit+')</text>';
  if(pwlState.type==='PWL'){
    if(pwlState.pts.length>1)
      s+='<polyline class="seg" points="'+pwlState.pts.map(p=>pwlWX(p[0])+','+pwlWY(p[1])).join(' ')+'"/>';
    pwlState.pts.forEach((p,i)=>{
      s+='<circle class="pt'+(i===pwlState.drag?' drag':'')+'" data-i="'+i+'" cx="'+pwlWX(p[0])+'" cy="'+pwlWY(p[1])+'" r="5"/>'; });
  } else {
    const wave=srcWave(pwlState.type, pwlState.params, v.xmin, v.xmax, 260);
    if(wave.length>1)
      s+='<polyline class="seg" points="'+wave.map(q=>pwlWX(q[0])+','+pwlWY(q[1])).join(' ')+'"/>';
    for(const h of srcHandles()){ const hx=pwlWX(h.t), hy=pwlWY(h.v);
      s+='<circle class="pt'+(h.id===pwlState.hdrag?' drag':'')+'" data-h="'+h.id+'" cx="'+hx+'" cy="'+hy+'" r="6"/>';
      s+='<text class="hlbl" x="'+(hx+8)+'" y="'+(hy-8)+'">'+h.label+'</text>'; }
  }
  svg.innerHTML=s;
  if(pwlState.type==='PWL') pwlSyncNum();
}
function pwlSyncNum(){
  const ta=document.getElementById('pwlNum'); if(!ta) return;
  if(document.activeElement===ta) return;   // don't clobber while the user types
  ta.value=pwlState.pts.map(p=>pwlFmt(p[0])+' '+pwlFmt(p[1])).join('\\n');
}
function pwlNumEdit(){
  if(!pwlState) return;
  const nums=document.getElementById('pwlNum').value.trim().split(/[\\s,]+/).map(parseFloat).filter(x=>!isNaN(x));
  const pts=[]; for(let i=0;i+1<nums.length;i+=2) pts.push([nums[i], nums[i+1]]);
  if(pts.length){ pwlState.pts=pts; pwlDraw(); }
}
function pwlSort(){ pwlState.pts.sort((a,b)=>a[0]-b[0]); }

// ---- interaction ----
function pwlDown(e){
  if(!pwlState) return;
  const svg=document.getElementById('pwlSvg'), r=svg.getBoundingClientRect();
  const px=e.clientX-r.left, py=e.clientY-r.top, tgt=e.target;
  if(tgt.classList && tgt.classList.contains('pt')){
    if(pwlState.type==='PWL'){
      const i=+tgt.getAttribute('data-i');
      if(e.shiftKey || e.button===2){
        if(pwlState.pts.length>1){ pwlState.pts.splice(i,1); pwlDraw(); }
        e.preventDefault(); return; }
      pwlState.drag=i; svg.setPointerCapture(e.pointerId); e.preventDefault(); pwlDraw(); return;
    } else {
      pwlState.hdrag=tgt.getAttribute('data-h'); svg.setPointerCapture(e.pointerId); e.preventDefault(); pwlDraw(); return;
    }
  }
  pwlState.pan={sx:px, sy:py, moved:false,
                vx:pwlState.view.xmin, vX:pwlState.view.xmax,
                vy:pwlState.view.ymin, vY:pwlState.view.ymax};
  svg.setPointerCapture(e.pointerId);
}
function pwlMove(e){
  if(!pwlState) return;
  const svg=document.getElementById('pwlSvg'), r=svg.getBoundingClientRect();
  const px=e.clientX-r.left, py=e.clientY-r.top, snap=document.getElementById('pwlSnap').checked;
  if(pwlState.type==='PWL' && pwlState.drag>=0){
    let t=pwlPX(px), y=pwlPY(py);
    if(snap){ t=pwlSnap(t,pwlState.view.xmin,pwlState.view.xmax); y=pwlSnap(y,pwlState.view.ymin,pwlState.view.ymax); }
    pwlState.pts[pwlState.drag]=[t,y]; pwlDraw(); return;
  }
  if(pwlState.hdrag){
    let t=pwlPX(px), y=pwlPY(py);
    if(snap){ t=pwlSnap(t,pwlState.view.xmin,pwlState.view.xmax); y=pwlSnap(y,pwlState.view.ymin,pwlState.view.ymax); }
    srcApplyHandle(pwlState.hdrag, t, y); srcRenderFields(); pwlDraw(); return;
  }
  if(pwlState.pan){
    const dx=px-pwlState.pan.sx, dy=py-pwlState.pan.sy;
    if(Math.abs(dx)>3 || Math.abs(dy)>3) pwlState.pan.moved=true;
    if(pwlState.pan.moved){
      const d=pwlDims(), v=pwlState.view;
      const sX=(pwlState.pan.vX-pwlState.pan.vx)/(d.x1-d.x0);   // world/px, >0
      const sY=(pwlState.pan.vY-pwlState.pan.vy)/(d.y1-d.y0);   // world/px, <0
      v.xmin=pwlState.pan.vx-dx*sX; v.xmax=pwlState.pan.vX-dx*sX;
      v.ymin=pwlState.pan.vy-dy*sY; v.ymax=pwlState.pan.vY-dy*sY;
      pwlDraw();
    }
  }
}
function pwlUp(e){
  if(!pwlState) return;
  const svg=document.getElementById('pwlSvg');
  try{ svg.releasePointerCapture(e.pointerId); }catch(_){}
  if(pwlState.type==='PWL' && pwlState.drag>=0){ pwlState.drag=-1; pwlSort(); pwlDraw(); return; }
  if(pwlState.hdrag){ pwlState.hdrag=null; pwlDraw(); return; }
  if(pwlState.pan){
    if(!pwlState.pan.moved && pwlState.type==='PWL'){
      const r=svg.getBoundingClientRect(), px=e.clientX-r.left, py=e.clientY-r.top, d=pwlDims();
      if(px>=d.x0 && px<=d.x1 && py>=d.y1 && py<=d.y0){
        let t=pwlPX(px), y=pwlPY(py);
        if(document.getElementById('pwlSnap').checked){
          t=pwlSnap(t,pwlState.view.xmin,pwlState.view.xmax); y=pwlSnap(y,pwlState.view.ymin,pwlState.view.ymax); }
        pwlState.pts.push([t,y]); pwlSort(); pwlDraw();
      }
    }
    pwlState.pan=null;
  }
}
function pwlWheel(e){
  if(!pwlState) return; e.preventDefault();
  const svg=document.getElementById('pwlSvg'), r=svg.getBoundingClientRect();
  pwlZoomAt(e.deltaY<0?1.15:1/1.15, {tx:pwlPX(e.clientX-r.left), ty:pwlPY(e.clientY-r.top)});
}

// ---- open / save ----
function pwlOpen(idx){
  const ta=document.getElementById('f_circuit_spec'); if(!ta) return;
  const p=srcParseLine((ta.value.split('\\n')[idx])||''); if(!p) return;
  pwlState={lineIdx:idx, kind:p.kind, type:p.type, name:p.name, a:p.a, b:p.b, view:null, drag:-1, hdrag:null, pan:null};
  if(p.type==='PWL') pwlState.pts = p.pts.length?p.pts.slice():[[0,0]];
  else pwlState.params = Object.assign({}, p.params);
  const isPwl=p.type==='PWL';
  document.getElementById('pwlTitle').textContent = p.type+' source editor';
  document.getElementById('pwlName').textContent =
    p.kind+p.type+' '+p.name+'   ('+p.a+' ↔ '+p.b+')   ['+(p.kind==='V'?'V':'A')+']';
  document.getElementById('pwlClearBtn').style.display = isPwl?'':'none';
  document.getElementById('pwlNumBox').style.display   = isPwl?'':'none';
  document.getElementById('pwlHint').textContent = isPwl
    ? 'click empty = add · drag = move · shift/right-click = delete · wheel = zoom · drag bg = pan'
    : 'drag ● handles = edit params · wheel = zoom · drag bg = pan';
  document.getElementById('pwlWarn').textContent='';
  srcRenderFields();
  document.getElementById('pwlModal').classList.add('on');  // show first so clientWidth is real
  pwlResetView();
}
function pwlCancel(){ const m=document.getElementById('pwlModal'); if(m) m.classList.remove('on'); pwlState=null; }
function pwlSave(){
  if(!pwlState) return;
  const w=document.getElementById('pwlWarn');
  let payload;
  if(pwlState.type==='PWL'){
    const pts=pwlState.pts.slice().sort((a,b)=>a[0]-b[0]);
    if(pts.length<2){ w.textContent='Need at least 2 points.'; return; }
    for(let i=1;i<pts.length;i++) if(pts[i][0]<=pts[i-1][0]){
      w.textContent='Times must strictly increase (duplicate/backward t).'; return; }
    payload=pts.map(p=>pwlFmt(p[0])+' '+pwlFmt(p[1])).join(' ');
  } else {
    const p=pwlState.params;
    if(pwlState.type==='SIN'){ if(!(p.freq>0)){ w.textContent='Frequency must be > 0.'; return; }
      payload=pwlFmt(p.amp)+' '+pwlFmt(p.freq); }
    else if(pwlState.type==='DC'){ payload=pwlFmt(p.value); }
    else if(pwlState.type==='PULSE'){ if(p.td<0){ w.textContent='Delay td must be ≥ 0.'; return; }
      if(!(p.tr>0)){ w.textContent='Rise time tr must be > 0.'; return; }
      payload=pwlFmt(p.v1)+' '+pwlFmt(p.v2)+' '+pwlFmt(p.td)+' '+pwlFmt(p.tr); }
    else if(pwlState.type==='PWM'){ if(!(p.freq>0)){ w.textContent='Frequency must be > 0.'; return; }
      const duty=Math.min(1,Math.max(0,p.duty));
      payload=pwlFmt(p.v1)+' '+pwlFmt(p.v2)+' '+pwlFmt(p.freq)+' '+pwlFmt(duty); }
  }
  const label=pwlState.kind+pwlState.type+' '+pwlState.name;
  const line=pwlState.kind+pwlState.type+' '+pwlState.name+' '+pwlState.a+' '+pwlState.b+' '+payload;
  const ta=document.getElementById('f_circuit_spec');
  const lines=ta.value.split('\\n'); lines[pwlState.lineIdx]=line; ta.value=lines.join('\\n');
  pwlCancel();
  pwlRenderButtons();
  if(typeof setStatus==='function') setStatus(label+' updated','');
  loadCircuit();
}

window.addEventListener('load', ()=>{
  const ctrls=document.querySelector('.prop-rows');   // (was '#controls', removed in the redesign)
  if(ctrls){ ctrls.addEventListener('input', applyVisibility); ctrls.addEventListener('change', applyVisibility); }
  // Help hover is document-wide now (help icons live in both the controls and the circuit box).
  document.addEventListener('mouseover', e=>{ if(e.target.classList.contains('help')) helpShow(e.target); });
  document.addEventListener('mouseout',  e=>{ if(e.target.classList.contains('help')) helpHide(); });
  buildLibrary(); applyVisibility(); linkWindows('N'); loadCircuit();
  // PWL editor: chips under the spec + modal SVG interaction wiring.
  pwlRenderButtons();
  addSweepParam('R_ROM'); updateSweepMode();   // seed one sweep-parameter row
  pbRefresh();   // probe target dropdowns
  const spec=document.getElementById('f_circuit_spec');
  if(spec){
    spec.addEventListener('input', ()=>{ clearTimeout(pwlDeb); pwlDeb=setTimeout(()=>{ pwlRenderButtons(); pbFillTargets(); },120); });
    spec.addEventListener('scroll', pwlRenderButtons);
  }
  const psvg=document.getElementById('pwlSvg');
  if(psvg){
    psvg.addEventListener('pointerdown', pwlDown);
    psvg.addEventListener('pointermove', pwlMove);
    psvg.addEventListener('pointerup', pwlUp);
    psvg.addEventListener('wheel', pwlWheel, {passive:false});
    psvg.addEventListener('contextmenu', e=>e.preventDefault());
  }
  const pmod=document.getElementById('pwlModal');
  if(pmod) pmod.addEventListener('pointerdown', e=>{ if(e.target.id==='pwlModal') pwlCancel(); });
  document.addEventListener('keydown', e=>{ if(e.key==='Escape' && pwlState) pwlCancel(); });
  window.addEventListener('resize', ()=>{ pwlRenderButtons(); if(pwlState) pwlDraw(); });
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
              .replace("__SWEEP_OPTS_JS__", _sweep_options_js())
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
    rows = run_sweep(dict(DEFAULTS), args.param, values, getattr(args, "jobs", None) or None)

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
    pw.add_argument("--jobs", type=int, default=0,
                    help="parallel solver processes (0 = auto / all cores)")
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
