#!/usr/bin/env python3
"""Regression: the FEM field solver carries the window seam (C0-continuous).

At N_field_eval_intervals=1 the field waveform must be exactly ONE straight segment per window:
the window start is the CARRIED previous-window end value (reused, never re-solved) and the window
end is the single backward-Euler step across the window. So Field_waveform_solution.prn holds
N_field_windows+1 rows -- one carried start plus one new node per window -- and every node after the
first equals that window's terminal value in Field_solution.prn.

This is the property that lets a single reconstruction (pointwise local BDF-1) cover what used to be
a separate 'linear ramp' mode: one interval leaves only the carried start and the window end, so the
backward-Euler step IS the ramp. A regression that recomputes the window start instead of carrying
it, or that emits extra interior nodes at N=1, breaks the row count; one that lets the waveform end
somewhere other than the terminal value breaks the value check. Guards both
FEM_solver_voltage_driven_waveform (field current, start = carried I0) and
FEM_solver_current_driven_waveform (field voltage, start = carried V_field_last_time).

Run:  python3 tests/test_seam_carry.py     (needs a Xyce on PATH; builds `main` via sim_ui)
"""
import importlib.util
import os
import shutil
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
spec = importlib.util.spec_from_file_location("sim_ui", os.path.join(ROOT, "sim_ui.py"))
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

# (name, base params, FEM-output waveform key) -- I for voltage-driven, V for current-driven.
CASES = [
    ("voltage-driven (P1)",
     {**m.DEFAULTS, **m.PRESETS["P1: RLC (sine U, series RL, shunt C)"],
      "coupling_mode": 0, "N_field_eval_intervals": 1},
     "I"),
    ("current-driven (P2)",
     {**m.DEFAULTS, **m.PRESETS["P2: Sine I (bare)"], "coupling_mode": 1, "N_field_eval_intervals": 1},
     "V"),
]


def solve(params):
    """Run one case in its own temp directory, so the working tree's config/outputs stay untouched."""
    d = tempfile.mkdtemp(prefix="wr_seam_")
    try:
        m.write_config(params, d)
        m.write_spec(params, d)
        m.write_probes(params, d)
        rc, out, err, _ = m.run_solver(d)
        assert rc == 0, f"solver exited {rc}\n{(err or out or '')[-500:]}"
        data = m.read_outputs(d)
        return data["field_wave"], data["field"]
    finally:
        shutil.rmtree(d, ignore_errors=True)


def main():
    m.ensure_built()
    failures = []
    for name, base, key in CASES:
        fw, fld = solve(base)
        n_win = int(base["N_field_windows"])
        problems = []
        # One carried start + one new node per window.
        if fw["t"].size != n_win + 1:
            problems.append(f"rows={fw['t'].size}, expected {n_win + 1}")
        # Every node after the start is that window's terminal value (Field_solution.prn prints
        # fewer digits than the waveform file, hence the relative tolerance rather than ==).
        elif fld["t"].size != n_win:
            problems.append(f"terminal rows={fld['t'].size}, expected {n_win}")
        else:
            if not np.allclose(fw["t"][1:], fld["t"], rtol=1e-9, atol=0.0):
                problems.append("window-end times differ from Field_solution.prn")
            bad = ~np.isclose(fw[key][1:], fld[key], rtol=1e-5, atol=0.0)
            if bad.any():
                worst = int(np.argmax(np.abs(fw[key][1:] - fld[key])))
                problems.append(f"{int(bad.sum())} window-end {key} values differ "
                                f"(worst at t={fld['t'][worst]:.6g}: "
                                f"{fw[key][1:][worst]:.9g} vs {fld[key][worst]:.9g})")
        print(f"[{'FAIL' if problems else 'PASS'}] {name}: one carried ramp per window @N_field_eval=1"
              + ("  -- " + "; ".join(problems) if problems else f"  ({fw['t'].size} rows)"))
        if problems:
            failures.append(name)
    if failures:
        print("SEAM REGRESSION FAILED (seam not carried / extra nodes / wrong window end):",
              ", ".join(failures))
        sys.exit(1)
    print("OK: seam carried (C0) -- N_field_eval=1 is one straight segment per window")


if __name__ == "__main__":
    main()
