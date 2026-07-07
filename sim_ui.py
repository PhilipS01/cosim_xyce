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
    ("source_kind",                     "Circuit source",               0,        "choice",
        {0: "sinusoidal voltage", 1: "sinusoidal current", 2: "step/ramp voltage"}),
    ("frequency",                       "Source frequency f (Hz)",      50.0,     "float", (1, 200, 1)),
    ("amplitude",                       "Source amplitude (V or A)",    1.0,      "float", (0.1, 10, 0.1)),
    ("step_v_initial",                  "Step: initial level",          0.0,      "float", None),
    ("step_v_final",                    "Step: final level",            1.0,      "float", None),
    ("step_delay",                      "Step: onset delay (s)",        0.0,      "float", None),
    ("step_rise",                       "Step: ramp/rise time (s)",     1.0e-4,   "float", None),
    ("R_series",                        "R_series src->port (Ohm)",     6.0e-3,   "float", None),
    ("L_series",                        "L_series src->port (H)",       1.6e-7,   "float", None),
    ("C_series",                        "C_series src->port (F, 0=off)",0.0,      "float", None),
    ("circuit_kind",                    "Circuit topology",             0,        "choice",
        {0: "simple source", 1: "#4 2-way (sine U,C)", 2: "#5 2-way (DC U,C)", 3: "#6 2-way (AC vs R)",
         4: "custom (node-graph spec)"}),
    ("switch_backend",                  "Switch backend",               0,        "choice",
        {0: "behavioral R", 1: "native S"}),
    ("switch_t1",                       "Switch t1 (s)",                6.0e-3,   "float", None),
    ("switch_C",                        "Switch cap C (F)",             1.0e-6,   "float", None),
    ("switch_R",                        "Switch R (Ohm, #6)",           1.0e4,    "float", None),
    ("switch_Ron",                      "Switch Ron closed (Ohm)",      1.0e-3,   "float", None),
    ("switch_Roff",                     "Switch Roff open (Ohm)",       1.0e9,    "float", None),
    ("switch_trise",                    "Switch transition (s)",        1.0e-5,   "float", None),
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
]
DEFAULTS = {k: d for (k, _l, d, _kind, _s) in PARAMS}
KINDS = {k: kind for (k, _l, _d, kind, _s) in PARAMS}
LABELS = {k: l for (k, l, _d, _kind, _s) in PARAMS}
# Numeric params are sweepable (a "choice" metric switch is not a continuum).
SWEEPABLE = [k for (k, _l, _d, kind, _s) in PARAMS if kind in ("float", "int")]

# Circuit-side presets (hybrid model): a preset seeds the editable primitive fields; the user
# may then tweak any field. Increment 1 covers the three source kinds + series R/L; presets 4-6
# (switches) arrive in increment 2. Keys map to the flat config the C++ generator consumes.
PRESETS = {
    "P1: Sine V + RL": {"circuit_kind": 0, "source_kind": 0, "amplitude": 1.0, "frequency": 50.0,
                        "R_series": 6.0e-3, "L_series": 1.6e-7, "C_series": 0.0, "time_mode": 0,
                        "coupling_mode": 0},
    # Bare current source directly on the port: series R/L/C are meaningless for a current drive
    # (the current is forced regardless) and an ideal I-source in series with L is degenerate.
    "P2: Sine I (bare)": {"circuit_kind": 0, "source_kind": 1, "amplitude": 1.0, "frequency": 50.0,
                          "R_series": 0.0, "L_series": 0.0, "C_series": 0.0, "time_mode": 0,
                          "coupling_mode": 1},
    "P3: Step/ramp V + RL": {"circuit_kind": 0, "source_kind": 2, "step_v_initial": 0.0,
                             "step_v_final": 1.0, "step_delay": 0.0, "step_rise": 1.0e-4,
                             "R_series": 6.0e-3, "L_series": 1.6e-7, "C_series": 0.0,
                             "time_mode": 1, "t_end": 2.0e-2, "N_field_windows": 50, "coupling_mode": 0,
                             # window 1 straddles the whole ramp edge (stiff transient) -> more WR iters
                             "WRmaxSteps": 40},
    # --- Increment 2: switch topologies. NOTE the passive values are numerical-survival defaults,
    # not physically tuned to the field scale -- see the notes: C must stay small enough for WR to
    # contract, and the cap<->coil freewheel (P4/P5) needs a damped closed switch (Ron~10) or its
    # ~undamped LC ring dt-collapses. Tune C / Ron / times to your field for a meaningful excitation.
    "P4: 2-way switch (sine U, C)": {"circuit_kind": 1, "switch_backend": 0, "amplitude": 1.0,
                                     "frequency": 50.0, "switch_C": 1.0e-6, "switch_Ron": 10.0,
                                     "switch_t1": 6.0e-3, "WRmaxSteps": 40,
                                     "time_mode": 1, "t_end": 2.0e-2, "N_field_windows": 50, "coupling_mode": 0},
    "P5: 2-way switch (DC U, C)": {"circuit_kind": 2, "switch_backend": 0, "amplitude": 1.0,
                                   "switch_C": 1.0e-6, "switch_Ron": 10.0,
                                   "switch_t1": 6.0e-3, "WRmaxSteps": 40,
                                   "time_mode": 1, "t_end": 2.0e-2, "N_field_windows": 50, "coupling_mode": 0},
    "P6: 2-way switch (AC vs R)": {"circuit_kind": 3, "switch_backend": 0, "amplitude": 1.0,
                                   "frequency": 50.0, "switch_R": 1.0e4, "switch_Ron": 1.0e-3,
                                   "switch_t1": 6.0e-3, "WRmaxSteps": 40,
                                   "time_mode": 1, "t_end": 2.0e-2, "N_field_windows": 50, "coupling_mode": 0},
}


# ---------------------------------------------------------------------------
# Running the solver
# ---------------------------------------------------------------------------
# Conditional visibility: key -> list of AND-condition dicts; a control is shown iff ANY dict fully
# matches the current control values (OR-of-ANDs). Keys absent here are always visible. The custom
# spec box + SVG editor are handled separately in JS (visible only when circuit_kind == 4).
VISIBLE_WHEN = {
    "source_kind": [{"circuit_kind": [0]}],
    "R_series":    [{"circuit_kind": [0]}],
    "L_series":    [{"circuit_kind": [0]}],
    "C_series":    [{"circuit_kind": [0]}],
    "step_v_initial": [{"circuit_kind": [0], "source_kind": [2]}],
    "step_v_final":   [{"circuit_kind": [0], "source_kind": [2]}],
    "step_delay":     [{"circuit_kind": [0], "source_kind": [2]}],
    "step_rise":      [{"circuit_kind": [0], "source_kind": [2]}],
    "amplitude":   [{"circuit_kind": [0, 1, 2, 3]}],
    "frequency":   [{"time_mode": [0]}, {"circuit_kind": [0], "source_kind": [0, 1]},
                    {"circuit_kind": [1, 3]}],
    "switch_backend": [{"circuit_kind": [1, 2, 3]}],
    "switch_t1":      [{"circuit_kind": [1, 2, 3]}],
    "switch_Ron":     [{"circuit_kind": [1, 2, 3]}],
    "switch_Roff":    [{"circuit_kind": [1, 2, 3]}],
    "switch_trise":   [{"circuit_kind": [1, 2, 3]}],
    "switch_C":       [{"circuit_kind": [1, 2]}],
    "switch_R":       [{"circuit_kind": [3]}],
    "N_periods":                       [{"time_mode": [0]}],
    "N_field_steps_per_source_period": [{"time_mode": [0]}],
    "t_end":            [{"time_mode": [1]}],
    "N_field_windows":  [{"time_mode": [1]}],
    "I_sat": [{"nonlin_model": [1]}],
    "reconstruct_mode": [{"coupling_mode": [1]}],
}


# Hover help: key -> HTML shown in a tooltip next to the control's label (a "?" icon).
HELP = {
    "reconstruct_mode": (
        "<div class='hh'>Field reconstruction (current-driven)</div>"
        "<table>"
        "<tr><th>mode</th><th>reconstruction</th><th>solves/win</th><th>note</th></tr>"
        "<tr><td>pointwise (secant)</td><td>V at each field point, accumulated secant</td>"
        "<td>N_field_eval</td><td>curve-following; interior uses the dummy's finite difference</td></tr>"
        "<tr><td>linear ramp</td><td>straight line carried-start &rarr; window-end</td>"
        "<td><b>1</b></td><td>cheapest; pure coupling reconstruction</td></tr>"
        "<tr><td>average</td><td>line, start = &frac12;(carried + end)</td>"
        "<td><b>1</b></td><td>linear + window-start damping (colleague)</td></tr>"
        "<tr><td>pointwise (central)</td><td>V at each point, central difference</td>"
        "<td>N_field_eval</td><td>lowest raw RMS, but the derivative is a dummy artifact</td></tr>"
        "</table>"
        "<div class='hn'>All modes carry the seam (C0-continuous). Accuracy is within ~1&ndash;2% "
        "across modes; the extra solves buy little. Recommend <b>linear</b> / <b>average</b> "
        "(1 field solve per window).</div>"
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


WIRE = "#8b94a3"


def _classify_source(e):
    """Map an element to a schematic symbol kind."""
    t = e["type"]
    if t == "PWL source":
        return "pwl"
    if t == "voltage src" and e.get("expr", "").strip() in ("0", "0.0", "DC 0"):
        return "ammeter"          # a 0 V source is an ammeter
    if t in ("voltage src",):
        return "vsource"
    if t in ("current src",):
        return "isource"
    if t == "behavioral V":
        return "bvsource"         # dependent voltage source (diamond)
    if t == "behavioral I":
        return "bisource"         # dependent current source (diamond)
    if t.startswith("resistor"):
        return "res"
    if t.startswith("inductor"):
        return "ind"
    if t.startswith("capacitor"):
        return "cap"
    return "box"


def _draw_symbol(ax, kind, M, r, d, color):
    """Draw a component symbol centred at M, 'radius' r, oriented along unit vec d
    (the wire direction). Perp p is d rotated 90 deg."""
    import numpy as _np
    cx, cy = M
    dx, dy = d
    px, py = -dy, dx                       # perpendicular
    circ_kinds = ("vsource", "isource", "ammeter", "pwl")
    diamond_kinds = ("bvsource", "bisource")

    def line(a, b, **kw):
        ax.plot([a[0], b[0]], [a[1], b[1]], color=kw.pop("c", color),
                lw=kw.pop("lw", 1.8), zorder=3, **kw)

    if kind in circ_kinds or kind in diamond_kinds:
        if kind in circ_kinds:
            ax.add_patch(plt.Circle(M, r, fill=True, fc="#0f1115",
                                    ec=color, lw=1.8, zorder=3))
        else:
            pts = [(cx + r * dx, cy + r * dy), (cx + r * px, cy + r * py),
                   (cx - r * dx, cy - r * dy), (cx - r * px, cy - r * py)]
            ax.add_patch(plt.Polygon(pts, closed=True, fill=True, fc="#0f1115",
                                     ec=color, lw=1.8, zorder=3))
        # Inner glyph.
        if kind in ("vsource", "bvsource"):       # sine '~'
            t = _np.linspace(-1, 1, 40)
            gx = cx + 0.55 * r * t
            gy = cy + 0.32 * r * _np.sin(_np.pi * t)
            ax.plot(gx, gy, color=color, lw=1.6, zorder=4)
        elif kind in ("isource", "bisource"):     # current arrow along d
            from matplotlib.patches import FancyArrowPatch
            a = (cx - 0.5 * r * dx, cy - 0.5 * r * dy)
            b = (cx + 0.5 * r * dx, cy + 0.5 * r * dy)
            ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>",
                                         mutation_scale=12, lw=1.6,
                                         color=color, zorder=4))
        elif kind == "ammeter":
            ax.text(cx, cy, "A", ha="center", va="center", fontsize=11,
                    color=color, fontweight="bold", zorder=4)
        elif kind == "pwl":                        # small pwl wave
            t = _np.array([-1, -0.4, 0.2, 0.7, 1])
            gx = cx + 0.6 * r * t
            gy = cy + 0.35 * r * _np.array([-1, 0.6, -0.3, 0.8, 0.1])
            ax.plot(gx, gy, color=color, lw=1.4, zorder=4)
        return r
    if kind == "cap":                               # two plates perpendicular
        g = 0.12 * r
        for s in (+1, -1):
            c0 = (cx + s * g * dx, cy + s * g * dy)
            line((c0[0] - r * px, c0[1] - r * py), (c0[0] + r * px, c0[1] + r * py))
        return g
    # Rectangle body (resistor/inductor/generic), long axis along d.
    hl, hw = r, 0.5 * r
    corners = [(cx + hl * dx + hw * px, cy + hl * dy + hw * py),
               (cx + hl * dx - hw * px, cy + hl * dy - hw * py),
               (cx - hl * dx - hw * px, cy - hl * dy - hw * py),
               (cx - hl * dx + hw * px, cy - hl * dy + hw * py)]
    ax.add_patch(plt.Polygon(corners, closed=True, fill=True, fc="#0f1115",
                             ec=color, lw=1.8, zorder=3))
    glyph = {"res": "R", "ind": "L"}.get(kind, "")
    if glyph:
        ax.text(cx, cy, glyph, ha="center", va="center", fontsize=9,
                color=color, fontweight="bold", zorder=4)
    return r


def make_circuit_plot(parsed):
    """Render the netlist as a traditional ladder schematic: non-ground nodes on a
    top line, a ground rail at the bottom, each conductive branch a component on a
    vertical leg (to ground) or a horizontal top segment (node-to-node). PWL signal
    sources sit below as reference inputs with dotted 'reads' arrows."""
    import math
    import re as _re
    import numpy as _np
    from matplotlib.patches import FancyArrowPatch
    nodes, elems = parsed["nodes"], parsed["elements"]
    if not nodes:
        return {}

    cond_elems = [e for e in elems if not e.get("signal")]
    sig_elems = [e for e in elems if e.get("signal")]
    cond_node_set = set()
    for e in cond_elems:
        cond_node_set.update(e["nodes"])
    GND = ("0", "gnd", "GND")
    signal_nodes = [n for n in nodes
                    if n not in cond_node_set
                    and any(n in e["nodes"] for e in sig_elems)]
    top_nodes = [n for n in nodes if n not in signal_nodes and n not in GND]

    # Bus layout: each non-ground node is a VERTICAL bus; ground is the bottom rail. Every 2-terminal
    # element is a rung -- node<->node = horizontal between two buses at its OWN y-level (so nothing
    # overlaps and it never reads as a series rail); node<->ground = vertical from the bus down to the
    # rail. Rungs cross intervening buses without a junction dot (= no connection, standard convention).
    DX, ROW, R = 2.6, 1.15, 0.32
    col = {n: i for i, n in enumerate(top_nodes)}
    X = lambda n: col[n] * DX
    x_lo = -1.2
    x_hi = (len(top_nodes) - 1) * DX + 1.2 if top_nodes else 1.2

    fig, ax = plt.subplots(figsize=(9.6, 6.4))
    fig.patch.set_facecolor("#0f1115")
    ax.set_facecolor("#12151b")
    ax.set_aspect("equal"); ax.axis("off")

    legend_seen = {}
    elem_center = {}

    def wire(a, b, c=WIRE, lw=2.0, ls="-", z=1, alpha=1.0):
        ax.plot([a[0], b[0]], [a[1], b[1]], color=c, lw=lw, ls=ls,
                zorder=z, alpha=alpha, solid_capstyle="round")

    def place_on_segment(e, P0, P1):
        """Draw wires P0->symbol->P1 with the component symbol at the midpoint."""
        x0, y0 = P0; x1, y1 = P1
        mx, my = (x0 + x1) / 2, (y0 + y1) / 2
        L = math.hypot(x1 - x0, y1 - y0) or 1.0
        dx, dy = (x1 - x0) / L, (y1 - y0) / L
        wire(P0, (mx - R * dx, my - R * dy))
        wire((mx + R * dx, my + R * dy), P1)
        _draw_symbol(ax, _classify_source(e), (mx, my), R, (dx, dy), e["color"])
        elem_center[e["name"]] = (mx, my)
        legend_seen[e["type"]] = e["color"]
        pxl, pyl = -dy, dx
        off = R + 0.24
        ax.text(mx + off * pxl, my + off * pyl, e["name"], fontsize=8.5,
                ha="center", va="center", color=e["color"], fontweight="bold", zorder=6)

    # Classify conductive elements: node<->node rungs vs node<->ground legs.
    rungs, glegs = [], {}
    for e in cond_elems:
        a, b = e["nodes"]; ag, bg = a in GND, b in GND
        if ag and bg:
            continue
        if ag ^ bg:
            glegs.setdefault(b if ag else a, []).append(e)
        elif a in col and b in col:
            rungs.append((a, b, e))

    GBAND = 1.0                       # ground legs occupy y in [0, GBAND]
    conn = {n: [] for n in top_nodes}  # node -> y-levels where a rung meets its bus
    top_y = GBAND + 0.4 + max(1, len(rungs)) * ROW
    y = top_y
    rung_rows = []
    for a, b, e in rungs:
        rung_rows.append((a, b, e, y)); conn[a].append(y); conn[b].append(y); y -= ROW

    # Ground rail + symbol.
    wire((x_lo, 0.0), (x_hi, 0.0), c=WIRE, lw=2.2)
    gx = (x_lo + x_hi) / 2
    wire((gx, 0.0), (gx, -0.18), c=WIRE)
    for i, w in enumerate((0.16, 0.10, 0.05)):
        yy = -0.18 - i * 0.07
        wire((gx - w, yy), (gx + w, yy), c=WIRE)
    ax.text(gx + 0.22, -0.30, "0", fontsize=9.5, ha="left", va="center",
            color="#e6e6e6", fontweight="bold")

    # Bus vertical extents (used to draw the buses AND to detect rung crossings).
    bus_ext = {}
    for n in top_nodes:
        cs = conn[n]; grounded = n in glegs
        if not cs and not grounded:
            continue
        y_hi = max(cs) if cs else GBAND
        y_lo = GBAND if grounded else (min(cs) if cs else GBAND)
        bus_ext[n] = (y_lo, y_hi)

    # Vertical buses (one per non-ground node), with a terminal + label at the top.
    for n, (y_lo, y_hi) in bus_ext.items():
        wire((X(n), y_lo), (X(n), y_hi))
        deg = len(conn[n]) + len(glegs.get(n, []))
        if deg >= 3:                                     # junction dots at real T-connections
            for yy in conn[n]:
                ax.plot([X(n)], [yy], marker="o", ms=5, color=WIRE, zorder=5)
        ax.plot([X(n)], [y_hi], marker="o", ms=8, color="#0d0f14",
                mec="tab:cyan", mew=1.7, zorder=5)
        ax.text(X(n), y_hi + 0.24, n, fontsize=10, ha="center", va="bottom",
                color="#e6e6e6", fontweight="bold", zorder=6)

    # node<->node rungs: horizontal at their own y-level, with a semicircular HOP wherever the rung
    # crosses an intervening bus it does NOT connect to (unambiguous 'wires cross, no connection').
    HOP = 0.13

    def draw_rung(e, yy, xa, xb):
        x0, x1 = (xa, xb) if xa <= xb else (xb, xa)
        crosses = sorted(X(n) for n, (lo, hi) in bus_ext.items()
                         if x0 < X(n) < x1 and lo - 1e-6 <= yy <= hi + 1e-6)
        # place the symbol in the widest bus-free gap (so it never sits on a crossing bus)
        posts = [x0] + crosses + [x1]
        gi = max(range(len(posts) - 1), key=lambda i: posts[i + 1] - posts[i])
        mx = (posts[gi] + posts[gi + 1]) / 2.0

        def seg(sa, sb):                                 # straight wire sa->sb at yy, hopping crosses
            cur = sa
            for c in crosses:
                if sa < c < sb:
                    wire((cur, yy), (c - HOP, yy))
                    th = _np.linspace(_np.pi, 0.0, 16)
                    ax.plot(c + HOP * _np.cos(th), yy + HOP * _np.sin(th),
                            color=WIRE, lw=2.0, zorder=1, solid_capstyle="round")
                    cur = c + HOP
            wire((cur, yy), (sb, yy))
        seg(x0, mx - R); seg(mx + R, x1)
        _draw_symbol(ax, _classify_source(e), (mx, yy), R, (1.0, 0.0), e["color"])
        elem_center[e["name"]] = (mx, yy); legend_seen[e["type"]] = e["color"]
        ax.text(mx, yy + R + 0.22, e["name"], fontsize=8.5, ha="center", va="bottom",
                color=e["color"], fontweight="bold", zorder=6)

    for a, b, e, yy in rung_rows:
        draw_rung(e, yy, X(a), X(b))

    # node<->ground legs: vertical from the bus band down to the rail (fan if several on a node).
    for n, es in glegs.items():
        if n not in col:
            continue
        cnt = len(es)
        for i, e in enumerate(es):
            ox = X(n) + (i - (cnt - 1) / 2.0) * 0.75
            if abs(ox - X(n)) > 1e-9:
                wire((X(n), GBAND), (ox, GBAND))
            place_on_segment(e, (ox, GBAND), (ox, 0.0))

    # PWL signal sources: a row beneath the rail, each feeding its reads-arrow.
    sy = -1.25
    sxs = _np.linspace(0.2, max(0.2, x_hi - 1.0), max(1, len(sig_elems)))
    for e, sx in zip(sig_elems, sxs):
        sx = float(sx)
        _draw_symbol(ax, "pwl", (sx, sy), R * 0.85, (1.0, 0.0), e["color"])
        elem_center[e["name"]] = (sx, sy)
        legend_seen["PWL signal source"] = e["color"]
        nd = next((n for n in e["nodes"] if n in signal_nodes), None)
        label = (nd or e["name"])
        ax.text(sx, sy - R - 0.12, f"{e['name']} ({label})", fontsize=8,
                ha="center", va="top", color=e["color"], fontweight="bold")

    # 'reads' arrows: signal node's source -> element whose expr uses V(node).
    for e in sig_elems:
        sn = next((n for n in e["nodes"] if n in signal_nodes), None)
        if not sn or e["name"] not in elem_center:
            continue
        pat = _re.compile(r"V\(\s*" + _re.escape(sn) + r"\s*\)", _re.I)
        for c in cond_elems:
            if pat.search(c.get("expr", "")) and c["name"] in elem_center:
                ax.add_patch(FancyArrowPatch(
                    elem_center[e["name"]], elem_center[c["name"]],
                    connectionstyle="arc3,rad=-0.18", arrowstyle="-|>",
                    mutation_scale=11, lw=1.1, color="#9aa4b2", ls=":",
                    shrinkA=14, shrinkB=16, zorder=0, alpha=0.85))

    # Legend.
    handles = [plt.Line2D([0], [0], color=c, lw=3, label=t)
               for t, c in sorted(legend_seen.items())]
    if sig_elems:
        handles.append(plt.Line2D([0], [0], color="#9aa4b2", lw=1.1, ls=":",
                                  label="reads V(node)"))
    if handles:
        ax.legend(handles=handles, fontsize=8, loc="upper left",
                  framealpha=0.3, facecolor="#181b22", edgecolor="#2c333f",
                  labelcolor="#e6e6e6")

    ax.set_xlim(x_lo - 0.6, x_hi + 0.6)
    ax.set_ylim(sy - 0.9, top_y + 0.7)
    ax.set_title(f"Circuit schematic ({len(cond_elems)} branches, "
                 f"{len(sig_elems)} signal sources, {len(nodes)} nodes)",
                 color="#9aa4b2")
    fig.tight_layout()
    return {"circuit": _png(fig)}


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
            self._send(200, json.dumps({
                "ok": True,
                "plots": make_circuit_plot(parsed),
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
        if self.path not in ("/run", "/sweep", "/netlist"):
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
    # Preset selector (hybrid model): seeds the fields below, then the user may edit them.
    preset_opts = '<option value="">— custom —</option>' + "".join(
        f'<option value="{name}">{name}</option>' for name in PRESETS
    )
    rows = [f"""
        <div class="ctl">
          <label for="presetSel">Circuit preset</label>
          <div class="inputs">
            <select id="presetSel" class="choice" onchange="applyPreset()">{preset_opts}</select>
          </div>
        </div>"""]
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
    # Custom node-graph spec editor (used when circuit_kind = custom). "Refresh circuit" re-renders it.
    import html as _html
    spec_seed = _html.escape(read_spec())
    rows.append(f"""
        <div class="ctl" id="specBox">
          <label for="f_circuit_spec">Custom circuit spec (circuit topology = custom)</label>
          <textarea id="f_circuit_spec" rows="9" spellcheck="false"
                    style="width:100%;background:#0d0f14;color:#e6e6e6;border:1px solid #2c333f;
                           border-radius:6px;padding:8px;font-family:ui-monospace,monospace;font-size:11.5px;"
          >{spec_seed}</textarea>
          <div class="note">Nodes: <code>p</code>=port, <code>0</code>=gnd. One element/line:
            <code>R/L/C name a b val</code>, <code>{{V,I}}SIN amp f</code>, <code>{{V,I}}DC val</code>,
            <code>{{V,I}}PULSE v1 v2 td tr</code>, <code>{{V,I}}PWL t1 v1 t2 v2 …</code>.
            Click "Refresh circuit" to render.</div>
        </div>""")
    return "\n".join(rows)


def _sweep_options_html():
    return "".join(
        f'<option value="{k}"{" selected" if k == "R_series" else ""}>{LABELS[k]}</option>'
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
  /* --- custom circuit editor (SVG drag-drop) --- */
  .cewrap { border:1px solid #2c333f; border-radius:8px; margin-bottom:12px; background:#0d0f14; }
  .cebar { display:flex; flex-wrap:wrap; gap:6px; padding:8px; border-bottom:1px solid #2c333f; align-items:center; }
  .cebar button { padding:5px 9px; font-size:12px; }
  .cebar .pal { background:#22305a; color:#cfe0ff; }
  .cebar .sep { flex:1; }
  #ceSvg { width:100%; height:360px; display:block; background:#12151b; border-radius:0 0 8px 8px; touch-action:none; }
  .ce-comp { cursor:grab; }
  .ce-term { fill:#0d0f14; stroke:var(--accent); stroke-width:1.5; cursor:crosshair; }
  .ce-term.pend { fill:var(--ok); stroke:var(--ok); }
  .ce-sel rect, .ce-sel circle.body, .ce-sel line.plate { stroke:var(--ok) !important; }
  .ce-wire { stroke:#8b94a3; stroke-width:2; }
  .ce-pin { fill:var(--err); }
  .ce-lbl { fill:#e6e6e6; font:10px ui-monospace,monospace; pointer-events:none; }
  .ce-hint { color:var(--muted); font-size:11px; padding:6px 8px; line-height:1.5; }
  .ce-props { padding:8px 10px; border-top:1px solid #2c333f; display:flex; flex-wrap:wrap; gap:8px 14px; align-items:flex-end; }
  .ce-props .cprow { flex-basis:100%; color:var(--muted); font-size:12px; }
  .ce-props label { color:var(--muted); font-size:11px; display:flex; flex-direction:column; gap:3px; }
  .ce-props input, .ce-props select { background:#0d0f14; color:#e6e6e6; border:1px solid #2c333f;
       border-radius:5px; padding:4px 6px; font-family:ui-monospace,monospace; font-size:11.5px; width:112px; }
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
        <span>Circuit (<code>wr_circuit.cir</code>)</span>
        <div class="btns">
          <button class="secondary small" id="toggleCircuitBtn" onclick="toggleCircuit()">Hide</button>
          <button class="secondary small" onclick="loadCircuit()">Refresh circuit</button>
        </div>
      </div>
      <div id="circuitContent">
        <div class="cewrap" id="ceWrap">
          <div class="cebar">
            <button class="pal" onclick="CE.add('V')">+V</button>
            <button class="pal" onclick="CE.add('I')">+I</button>
            <button class="pal" onclick="CE.add('R')">+R</button>
            <button class="pal" onclick="CE.add('L')">+L</button>
            <button class="pal" onclick="CE.add('C')">+C</button>
            <span class="sep"></span>
            <button class="secondary small" onclick="CE.del()">Delete sel</button>
            <button class="secondary small" onclick="CE.clear()">Clear</button>
            <button class="small" onclick="CE.apply()">&rarr; Use as circuit</button>
          </div>
          <svg id="ceSvg"></svg>
          <div class="ce-props" id="ceProps"></div>
          <div class="ce-hint">Palette adds a component &middot; drag bodies to move &middot; click a
            terminal then another to wire &middot; select a component to edit its type/values below
            &middot; select + "Delete sel" removes. Red pins <b>p</b>=port, <b>0</b>=ground.
            "Use as circuit" writes the spec &amp; sets topology = custom.</div>
        </div>
        <img class="plot" id="p_circuit" hidden>
        <div id="netlistTable"></div>
        <details><summary>Raw netlist + directives</summary><pre id="netlistRaw"></pre></details>
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
  const custom = ctlVal('circuit_kind')===4;
  const sb=document.getElementById('specBox'); if(sb) sb.style.display = custom?'':'none';
  const ce=document.getElementById('ceWrap'); if(ce) ce.style.display = custom?'':'none';
}
// smart default: sine simple source -> periods; step/switch/custom -> absolute end time
function suggestTimeMode(){ const ck=ctlVal('circuit_kind'), sk=ctlVal('source_kind'); return (ck===0 && (sk===0||sk===1))?0:1; }
// current-source circuit -> current-driven coupling (avoids the window-start V(p) secant spike)
function suggestCouplingMode(){ const ck=ctlVal('circuit_kind'), sk=ctlVal('source_kind'); return (ck===0 && sk===1)?1:0; }
function onTopoChange(){
  const tm=document.getElementById('f_time_mode'); if(tm) tm.value=String(suggestTimeMode());
  const cm=document.getElementById('f_coupling_mode'); if(cm) cm.value=String(suggestCouplingMode());
  applyVisibility();
}

function applyPreset(){
  const name = document.getElementById('presetSel').value;
  const p = PRESETS[name];
  if (!p) return;
  for (const k in p){
    const b = document.getElementById('f_'+k);
    if (b){ b.value = p[k]; if (b.tagName !== 'SELECT') syncFromBox(b); }
  }
  setStatus('Preset applied: '+name, '');
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
      setPlot('p_circuit', j.plots ? j.plots.circuit : null);
      showNetlist(j);
    }
  } catch(e){ /* leave circuit box empty on failure */ }
}
// ---- Custom circuit editor (increment 3b): hand-rolled SVG drag-drop -> circuit_spec ----
const CE = (function(){
  const NS='http://www.w3.org/2000/svg', GRID=20, HW=30, TR=6;
  const defP={V:{amp:1,freq:50},I:{amp:1,freq:50},R:{val:10000},L:{val:1.6e-7},C:{val:1e-6}};
  let comps=[], wires=[], seq={V:0,I:0,R:0,L:0,C:0}, sel=null, pend=null, drag=null;
  const PIN={'PIN:p':{x:90,y:250,label:'p'},'PIN:0':{x:90,y:315,label:'0'}};
  const S=()=>document.getElementById('ceSvg');
  function E(t,a){const e=document.createElementNS(NS,t);for(const k in a)e.setAttribute(k,a[k]);return e;}
  const snap=v=>Math.round(v/GRID)*GRID;
  const termXY=(c,i)=>({x:c.x+(i?HW:-HW),y:c.y});
  function tpos(id){ if(id in PIN)return PIN[id]; const p=id.split('#'); const c=comps.find(x=>x.id===p[0]); return c?termXY(c,+p[1]):{x:0,y:0}; }
  function add(type){ seq[type]++; const c={id:'k'+Math.random().toString(36).slice(2,7),type,name:type.toLowerCase()+seq[type],x:snap(320),y:snap(70+(comps.length%6)*45),sub:((type==='V'||type==='I')?'SIN':null),p:Object.assign({},defP[type])}; comps.push(c); sel=c.id; render(); renderProps(); }
  function del(){ if(!sel)return; comps=comps.filter(c=>c.id!==sel); wires=wires.filter(w=>w.a.split('#')[0]!==sel&&w.b.split('#')[0]!==sel); sel=null; pend=null; render(); renderProps(); }
  function clr(){ comps=[]; wires=[]; sel=null; pend=null; render(); renderProps(); }
  function label(c){
    if(c.type==='V'||c.type==='I'){
      if(c.sub==='DC')return c.name+' DC '+c.p.val;
      if(c.sub==='PULSE')return c.name+' pulse';
      if(c.sub==='PWL')return c.name+' pwl';
      return c.name+' '+c.p.amp+'/'+c.p.freq+'Hz';
    }
    return c.name+' '+c.p.val;
  }
  function render(){
    const s=S(); if(!s)return; while(s.firstChild)s.removeChild(s.firstChild);
    wires.forEach(w=>{const A=tpos(w.a),B=tpos(w.b); s.appendChild(E('line',{class:'ce-wire',x1:A.x,y1:A.y,x2:B.x,y2:B.y}));});
    for(const id in PIN){const p=PIN[id];
      s.appendChild(E('circle',{class:'ce-term'+(pend===id?' pend':''),cx:p.x,cy:p.y,r:TR,'data-term':id}));
      s.appendChild(E('circle',{class:'ce-pin',cx:p.x,cy:p.y,r:3}));
      const t=E('text',{class:'ce-lbl',x:p.x-20,y:p.y+4}); t.textContent=p.label; s.appendChild(t);}
    comps.forEach(c=>{
      const g=E('g',{class:'ce-comp'+(sel===c.id?' ce-sel':''),'data-comp':c.id});
      g.appendChild(E('line',{x1:c.x-HW,y1:c.y,x2:c.x-16,y2:c.y,stroke:'#8b94a3','stroke-width':2}));
      g.appendChild(E('line',{x1:c.x+16,y1:c.y,x2:c.x+HW,y2:c.y,stroke:'#8b94a3','stroke-width':2}));
      if(c.type==='V'||c.type==='I') g.appendChild(E('circle',{class:'body',cx:c.x,cy:c.y,r:16,fill:'#0f1115',stroke:'var(--accent)','stroke-width':1.8}));
      else if(c.type==='C'){ g.appendChild(E('line',{class:'plate',x1:c.x-4,y1:c.y-13,x2:c.x-4,y2:c.y+13,stroke:'var(--accent)','stroke-width':2})); g.appendChild(E('line',{class:'plate',x1:c.x+4,y1:c.y-13,x2:c.x+4,y2:c.y+13,stroke:'var(--accent)','stroke-width':2})); }
      else g.appendChild(E('rect',{x:c.x-16,y:c.y-10,width:32,height:20,rx:3,fill:'#0f1115',stroke:'var(--accent)','stroke-width':1.8}));
      const gl=E('text',{class:'ce-lbl',x:c.x,y:c.y+4,'text-anchor':'middle'}); gl.textContent=(c.type==='V'?'~':c.type==='I'?'↑':c.type==='C'?'':c.type); g.appendChild(gl);
      [0,1].forEach(i=>{const tp=termXY(c,i); g.appendChild(E('circle',{class:'ce-term'+(pend===c.id+'#'+i?' pend':''),cx:tp.x,cy:tp.y,r:TR,'data-term':c.id+'#'+i}));});
      const nl=E('text',{class:'ce-lbl',x:c.x,y:c.y-15,'text-anchor':'middle'}); nl.textContent=label(c); g.appendChild(nl);
      s.appendChild(g);
    });
  }
  function xy(e){const r=S().getBoundingClientRect(); const t=e.touches&&e.touches[0]; return {x:(t?t.clientX:e.clientX)-r.left,y:(t?t.clientY:e.clientY)-r.top};}
  function onDown(e){
    const term=e.target.getAttribute&&e.target.getAttribute('data-term');
    if(term){ if(pend===null)pend=term; else{ if(pend!==term)wires.push({a:pend,b:term}); pend=null; } render(); e.preventDefault(); return; }
    const g=e.target.closest&&e.target.closest('[data-comp]');
    if(g){ const id=g.getAttribute('data-comp'); sel=id; const c=comps.find(x=>x.id===id); const p=xy(e); drag={id,dx:c.x-p.x,dy:c.y-p.y}; render(); renderProps(); e.preventDefault(); return; }
    sel=null; pend=null; render(); renderProps();
  }
  function onMove(e){ if(!drag)return; const c=comps.find(x=>x.id===drag.id); if(!c)return; const p=xy(e); c.x=snap(p.x+drag.dx); c.y=snap(p.y+drag.dy); render(); e.preventDefault(); }
  function onUp(){ drag=null; }
  // --- properties panel: edit the selected component's type + type-dependent attributes ---
  const subDef={SIN:{amp:1,freq:50},DC:{val:1},PULSE:{v1:0,v2:1,td:0,tr:1e-4},PWL:{pts:'0 0 5e-3 1 15e-3 1 20e-3 0'}};
  const subName={SIN:'sinusoidal',DC:'DC',PULSE:'pulse',PWL:'PWL (multi-step)'};
  function fld(lbl,inner){ return '<label>'+lbl+inner+'</label>'; }
  function inp(key,val){ return '<input data-cp="'+key+'" value="'+val+'">'; }
  function renderProps(){
    const host=document.getElementById('ceProps'); if(!host)return;
    const c=comps.find(x=>x.id===sel);
    if(!c){ host.innerHTML='<span class="ce-hint" style="padding:0">Select a component to edit its type &amp; values.</span>'; return; }
    let h='<div class="cprow"><b>'+c.name+'</b> — '+c.type+'</div>';
    if(c.type==='V'||c.type==='I'){
      const u=c.type==='V'?'V':'A';
      h+='<label>Type<select data-cp="sub">'
        +['SIN','DC','PULSE','PWL'].map(s=>'<option value="'+s+'"'+(c.sub===s?' selected':'')+'>'+subName[s]+'</option>').join('')
        +'</select></label>';
      if(c.sub==='SIN')        h+=fld('Amplitude ('+u+')',inp('amp',c.p.amp))+fld('Frequency (Hz)',inp('freq',c.p.freq));
      else if(c.sub==='DC')    h+=fld('Value ('+u+')',inp('val',c.p.val));
      else if(c.sub==='PULSE') h+=fld('Initial ('+u+')',inp('v1',c.p.v1))+fld('Pulsed ('+u+')',inp('v2',c.p.v2))
                                 +fld('Delay td (s)',inp('td',c.p.td))+fld('Rise tr (s)',inp('tr',c.p.tr));
      else                     h+='<label style="flex-basis:100%">Points &mdash; t v t v … ('+u+', times increasing)'
                                 +'<input data-cp="pts" value="'+c.p.pts+'" style="width:100%"></label>';
    } else {
      h+=fld(c.type+' value',inp('val',c.p.val));
    }
    host.innerHTML=h;
    host.querySelectorAll('[data-cp]').forEach(el2=>{
      el2.addEventListener(el2.tagName==='SELECT'?'change':'input',()=>onProp(el2));
    });
  }
  function onProp(el2){
    const c=comps.find(x=>x.id===sel); if(!c)return;
    const key=el2.getAttribute('data-cp'), v=el2.value;
    if(key==='sub'){ c.sub=v; c.p=Object.assign({},subDef[v]); renderProps(); render(); return; }
    c.p[key]=v.trim(); render();   // live-update the on-canvas label
  }
  function serialize(){
    const par={}, find=x=>{par[x]=par[x]||x; return par[x]===x?x:(par[x]=find(par[x]));}, uni=(a,b)=>{par[find(a)]=find(b);};
    find('PIN:p'); find('PIN:0'); comps.forEach(c=>{find(c.id+'#0');find(c.id+'#1');});
    wires.forEach(w=>uni(w.a,w.b));
    const rn={}; rn[find('PIN:p')]='p'; rn[find('PIN:0')]='0'; let n=1;
    const nf=id=>{const r=find(id); if(!(r in rn))rn[r]='n'+(n++); return rn[r];};
    return comps.map(c=>{const a=nf(c.id+'#0'),b=nf(c.id+'#1');
      if(c.type==='V'||c.type==='I'){
        const P=c.type;   // 'V' or 'I' prefix -> VDC/IDC, VPULSE/IPULSE, VPWL/IPWL, VSIN/ISIN
        if(c.sub==='DC')    return P+'DC '+c.name+' '+a+' '+b+' '+c.p.val;
        if(c.sub==='PULSE') return P+'PULSE '+c.name+' '+a+' '+b+' '+c.p.v1+' '+c.p.v2+' '+c.p.td+' '+c.p.tr;
        if(c.sub==='PWL')   return P+'PWL '+c.name+' '+a+' '+b+' '+c.p.pts;
        return P+'SIN '+c.name+' '+a+' '+b+' '+c.p.amp+' '+c.p.freq;
      }
      return c.type+' '+c.name+' '+a+' '+b+' '+c.p.val;
    }).join('\\n');
  }
  function apply(){
    if(!comps.length){ setStatus('Editor empty — add components first.','err'); return; }
    const ta=document.getElementById('f_circuit_spec');
    if(ta) ta.value='# generated by circuit editor\\n'+serialize()+'\\n';
    const ks=document.getElementById('f_circuit_kind'); if(ks) ks.value='4';
    setStatus('Circuit applied to spec (topology = custom).','');
    loadCircuit();
  }
  function init(){ const s=S(); if(!s)return;
    s.addEventListener('mousedown',onDown); s.addEventListener('mousemove',onMove); window.addEventListener('mouseup',onUp);
    s.addEventListener('touchstart',onDown,{passive:false}); s.addEventListener('touchmove',onMove,{passive:false}); window.addEventListener('touchend',onUp);
    render(); renderProps();
  }
  return {add,del,clear:clr,apply,init};
})();

window.addEventListener('load', ()=>{
  const ctrls=document.getElementById('controls');
  if(ctrls){ ctrls.addEventListener('input', applyVisibility); ctrls.addEventListener('change', applyVisibility);
    ctrls.addEventListener('mouseover', e=>{ if(e.target.classList.contains('help')) helpShow(e.target); });
    ctrls.addEventListener('mouseout',  e=>{ if(e.target.classList.contains('help')) helpHide(); }); }
  const ck=document.getElementById('f_circuit_kind'); if(ck) ck.addEventListener('change', onTopoChange);
  const sk=document.getElementById('f_source_kind'); if(sk) sk.addEventListener('change', onTopoChange);
  applyVisibility(); loadCircuit(); CE.init();
});
</script>
</body></html>
"""

INDEX_HTML = (INDEX_HTML
              .replace("__CONTROLS__", _controls_html())
              .replace("__SWEEP_OPTS__", _sweep_options_html())
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
