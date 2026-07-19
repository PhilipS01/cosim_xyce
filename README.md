# cosim_xyce — WR field/circuit co-simulation

Waveform-relaxation (WR) co-simulation coupling a **Xyce** SPICE circuit to a
reduced-order / dummy **FEM field** model, with a browser studio to explore it.
Bachelor thesis project.

The circuit and the field are solved **separately** and iterated to a fixpoint on
each time window (waveform relaxation). The circuit is a Xyce netlist; the field
is a dummy solver in C++. They exchange port voltage and interface
current through a matched-secant interface source (`Bfield`) once per iteration,
until the window converges below `WR_tolerance`. The interface contains a ROM 
(preconditioner) of the arbitrarily complex field model. The run is split into
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
./setup.sh                        # macOS/Linux: venv + deps + build + optional tools
python3 sim_ui.py                 # serves http://127.0.0.1:8000
```

Or by hand:

```sh
pip install numpy matplotlib      # core deps
npm install                       # optional: schematic layout
python3 sim_ui.py                 # builds ./main, serves http://127.0.0.1:8000
```

Full usage, parameters, and the circuit-spec format: **[UI_README.md](UI_README.md)**.

## Windows

`setup.sh` is bash (macOS/Linux). On Windows, two routes:

- **WSL2 (recommended)** — install WSL2 + a Linux distro, then run `./setup.sh`
  inside it; from there it's identical to Linux
  (`sudo apt-get install -y g++ make python3 python3-venv python3-pip poppler-utils`,
  plus the Xyce Linux build). Open `http://127.0.0.1:8000` in the Windows browser.
- **Native** — run `./setup.ps1` in PowerShell. It needs a **MinGW/Clang**
  toolchain + `make` (MSVC is *not* used) — e.g. MSYS2
  (`pacman -S mingw-w64-ucrt-x86_64-gcc make`) or scoop. The build emits
  `main.exe` (the UI resolves the name automatically). The inline schematic PNG
  needs **poppler** (`pdftoppm`) on `PATH` — there is no `sips` on Windows — but
  simulations, plots, and `.tex`/`.pdf` export work without it.

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

## Generated at runtime

None of the following are source — they are (re)written on every run and are
git-ignored (except `wr_circuit.cir` and `sim_params.inc`, which are regenerated
but happen to be tracked). Deleting them is safe; the next run recreates them.

**Netlist + per-window includes** — written by the C++ driver before each Xyce call:

- `wr_circuit.cir` — the SPICE deck generated from `circuit_spec.txt`: the user
  circuit + the fixed WR interface (`Vmeas` ammeter, the matched-secant `Bfield`
  source, the two PWL feedback sources) + `.INCLUDE`s + `.print`/`.tran`. The
  topology is window-invariant, so this is regenerated once per run (and by
  `./main emit`).
- `sim_params.inc` — the per-window `.PARAM` values the deck references:
  `t_start`, `t_stop`, `t_abs_start`, `I0`, `Rrom`, `Lrom`, `dt_print`,
  `t_floor`. Rewritten each window.
- `restart.inc` — the per-window `.OPTIONS RESTART` + `.tran` line. Window 1 is a
  fresh UIC transient from `t=0`; window k>1 restarts from the previous
  checkpoint.

**Interface exchange** — the two coupled waveforms, rewritten every WR iteration
(2-column `time value` PWL tables, read back by the deck as `VFprev`/`VIprev`):

- `vf_prev_k.pwl` — the `V(vfprev)` base of the interface `Bfield`. Voltage-driven
  coupling: the circuit's port voltage `V(p)`; current-driven: the field's own
  computed voltage `V_field`.
- `i_prev_k.pwl` — the `V(iprev)` base. Voltage-driven: the field current
  `I_field`; current-driven: the interface current `I(Vmeas)`.

**Checkpoint / restart** — window-to-window continuity:

- `ckpt_out<time>` — Xyce restart checkpoints (e.g. `ckpt_out0`, `ckpt_out0.02`),
  one per window end (`JOB=ckpt_out`).
- `restart_state` — the newest checkpoint copied out as the restart basis for the
  next window.

**Terminal handoff** — single-line `V I` scalars passed between solver stages:

- `Field.txt` — the field's window-end terminal `(V_field, I_field)`.
- `Circuit.txt` — the circuit's window-end terminal `(V(nx), I(Vmeas))`.

**Solver output (parsed for the plots):**

- `wr_circuit.cir.prn` — Xyce's raw transient output for the current window
  (`Index TIME V(P) V(NX) I(VMEAS)` + any probe columns); parsed by the driver
  each iteration.
- `Circuit_solution.prn` — the stitched circuit solution across all windows
  (`V(P) V(NX) I(VMEAS)` vs absolute time).
- `Field_solution.prn` — the field solution at the window endpoints
  (`V(FIELD) I(FIELD)`).
- `Field_waveform_solution.prn` — the reconstructed/extrapolated field waveforms
  over each window (not only the endpoints).
- `WR_error.txt` — per-window convergence log:
  `Time, WR_TotalRelErr, N_iterations, Converged`.
- `Probes_solution.prn` — the user-requested probe columns (from `probes.txt`),
  captured per converged window. Only written if probes are set.

**Logs:** `xyce_stdout.log`, `xyce_stderr.log` — the last Xyce invocation's output.

The UI also writes three **input** files the solver reads: `sim_config.txt`
(parameters), `circuit_spec.txt` (the circuit), and `probes.txt` (output probes).

## Checkpointing (window-to-window restart)

The run is split into `N_field_windows` windows over `[0, t_end]`. Instead of one
long transient, each window is a **separate Xyce invocation that restarts from
the converged end state of the previous window** — so time is continuous and the
circuit state (capacitor voltages, inductor currents) carries across the seam.
This uses Xyce's native checkpoint/restart, orchestrated per window by the C++
driver. For each window `[t_start, t_stop]`:

1. **Write restart directives** (`restart.inc`, via `WriteRestartDirectives`):
   - Window 1: `.OPTIONS RESTART JOB=ckpt_out INITIAL_INTERVAL=<dt_window>` +
     `.tran ... UIC` — a fresh transient from `t=0` with user initial conditions.
   - Window k>1: `.OPTIONS RESTART FILE=restart_state JOB=ckpt_out
     INITIAL_INTERVAL=<dt_window>` + `.tran ...` (no UIC) — restart from the
     previous window's committed state.

   `INITIAL_INTERVAL = dt_window` makes Xyce drop **exactly one** checkpoint, at
   the (absolute) window end, named `ckpt_out<time>`.
2. **Clear stale candidates** (`ClearCheckpoints("ckpt_out")`) so a leftover
   checkpoint from an earlier window can't be committed by mistake.
3. **Run the WR loop.** Every WR iteration re-solves the *same* window as a full
   transient, restarting from the *same* `restart_state`, changing only the
   interface feedback (`vf_prev_k.pwl` / `i_prev_k.pwl`). Each iteration
   overwrites the same candidate checkpoint, so the last one written holds the
   converged state.
4. **Commit** (`CommitCheckpoint("ckpt_out", "restart_state")`) once the window
   converges: copy the converged checkpoint to `restart_state` — the restart
   basis for the next window.

**Gotcha:** the checkpoint to commit is chosen by the **largest time suffix**,
not file mtime — on restart Xyce also writes a `ckpt_out0` at `t=0`, so the
window-end state is always the one with the greatest sim-time.

The checkpoint carries the **circuit** state; the interface seam values that seed
the *next* window's feedback PWLs (`V0`, `I0`, and the ramp slopes) are handed
over separately via `Field.txt` / `Circuit.txt` (blend controlled by
`seam_average`).
