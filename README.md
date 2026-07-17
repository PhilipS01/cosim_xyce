# cosim_xyce — WR field/circuit co-simulation

Waveform-relaxation (WR) co-simulation coupling a **Xyce** SPICE circuit to a
reduced-order / dummy **FEM field** model, with a browser studio to explore it.
Bachelor thesis project.

The circuit and the field are solved **separately** and iterated to a fixpoint on
each time window (waveform relaxation). The circuit is a Xyce netlist; the field
is a reduced-order/dummy solver in C++. They exchange port voltage and interface
current through a matched-secant interface source (`Bfield`) once per iteration,
until the window converges below `WR_tolerance`. The run is split into
`N_field_windows` windows over `[0, t_end]`, each restarted from the previous via
a Xyce checkpoint.

## Requirements

- **Python 3** + `numpy` + `matplotlib` (`pip install numpy matplotlib`).
- **C++17 compiler** + `make` (the `Makefile` uses `clang++`).
- **Xyce** on `PATH` — the circuit solver.
- *Optional:* **Node.js** + `npm install` (elkjs) and **pdflatex** + `circuitikz`
  for the inline schematic and `.tex`/`.pdf` export.

## Quick start

```sh
pip install numpy matplotlib      # core deps
npm install                       # optional: schematic layout
python3 sim_ui.py                 # builds ./main, serves http://127.0.0.1:8000
```

Full usage, parameters, and the circuit-spec format: **[UI_README.md](UI_README.md)**.

## Run the solver directly

```sh
make                 # build ./main (clang++, C++17)
./main [config]      # run WR co-sim (defaults to sim_config.txt)
./main emit [cfg]    # only regenerate wr_circuit.cir (no solve)
```

## Layout

| Path | Role |
|------|------|
| `sim_ui.py` | WR Co-Simulation Studio — browser UI + headless parameter sweeps |
| `src/`, `include/Header.h` | C++ backend (`main`): config, netlist generation, WR driver, dummy/ROM FEM solver, checkpoint restart |
| `Makefile` | builds `main` |
| `circuit_spec.txt` | the authored circuit (custom node-graph; see UI_README) |
| `sim_config.txt` | solver / coupling parameters (`key = value`) |
| `elk_layout.js`, `package.json` | optional schematic layout (elkjs) |

Generated at runtime (not source): `wr_circuit.cir`, `sim_params.inc`,
`restart.inc`, `*.pwl`, `*.prn`, `ckpt_out*`, `restart_state`.
