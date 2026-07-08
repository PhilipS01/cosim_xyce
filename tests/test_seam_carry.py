#!/usr/bin/env python3
"""Regression: the FEM field solver carries the window seam (C0-continuous).

At N_field_eval_intervals=1 the 'pointwise (secant)' and 'linear ramp' reconstruction modes must
produce a BIT-IDENTICAL field waveform: both use the same two window endpoints (the carried
window-start value + the accumulated-secant window-end value), so with no interior points they
coincide exactly. This guards the seam-carry in FEM_solver_voltage_driven_waveform (field current,
start = carried I0) and FEM_solver_current_driven_waveform (field voltage, start = carried
V_field_last_time). A regression that recomputes the window start instead of carrying it, or that
changes the window-end value between the two modes, breaks this.

Run:  python3 tests/test_seam_carry.py     (needs a Xyce on PATH; builds `main` via sim_ui)
"""
import importlib.util
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
spec = importlib.util.spec_from_file_location("sim_ui", os.path.join(ROOT, "sim_ui.py"))
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

# (name, base params, FEM-output waveform key) -- I for voltage-driven, V for current-driven.
CASES = [
    ("voltage-driven (P1)",
     {**m.DEFAULTS, **m.PRESETS["P1: Sine V + RL"], "coupling_mode": 0, "N_field_eval_intervals": 1},
     "I"),
    ("current-driven (P2)",
     {**m.DEFAULTS, **m.PRESETS["P2: Sine I (bare)"], "coupling_mode": 1, "N_field_eval_intervals": 1},
     "V"),
]


def field_wave(params, key):
    m.write_config(params)
    m.write_spec(params)
    rc, out, err, _ = m.run_solver()
    assert rc == 0, f"solver exited {rc}\n{(err or '')[-500:]}"
    fw = m.read_outputs()["field_wave"]
    return fw[key].copy()


def main():
    m.ensure_built()
    failures = []
    try:
        for name, base, key in CASES:
            w0 = field_wave({**base, "reconstruct_mode": 0}, key)   # pointwise (secant)
            w1 = field_wave({**base, "reconstruct_mode": 1}, key)   # linear ramp
            same = len(w0) == len(w1) and float(np.max(np.abs(w0 - w1))) == 0.0
            delta = "len-mismatch" if len(w0) != len(w1) else f"{float(np.max(np.abs(w0 - w1))):.3e}"
            print(f"[{'PASS' if same else 'FAIL'}] {name}: mode0 vs linear @N=1  max|delta|={delta}")
            if not same:
                failures.append(name)
    finally:
        # Restore the tracked config/netlist to defaults so the working tree isn't left dirtied.
        m.write_config({**m.DEFAULTS})
        m.write_spec({**m.DEFAULTS})
        m.emit_netlist()
    if failures:
        print("SEAM REGRESSION FAILED (seam not carried / window-end differs):", ", ".join(failures))
        sys.exit(1)
    print("OK: seam carried (C0) -- mode 0 and linear coincide at N_field_eval=1")


if __name__ == "__main__":
    main()
