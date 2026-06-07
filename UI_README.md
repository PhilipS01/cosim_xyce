# WR Co-Simulation Studio

Browser UI to study the waveform-relaxation field/circuit co-simulation by
tuning solver parameters and viewing the resulting waveforms + WR convergence.

## Run

```sh
python3 sim_ui.py            # web UI, then open http://127.0.0.1:8000
python3 sim_ui.py serve --port 9000          # custom port
python3 sim_ui.py sweep --param R_series --min 0 --max 0.05 --steps 6   # headless convergence study
```

Requires Python 3 + `matplotlib` + `numpy` (already installed). `Xyce` must be
on `PATH`. The UI builds `./main` automatically if it is missing.

## What it does

1. You set parameters in the left panel and click **Run simulation**.
2. The backend writes `sim_config.txt`, runs `./main`, parses the `.prn`
   outputs, and renders three plots (rendered server-side as PNGs, so it works
   fully offline — no internet/CDN needed):
   - **Port / field voltage** vs time
   - **Interface current** (circuit `I(Vmeas)` vs field `I_field`) vs time
   - **WR convergence** per window (final rel. error + iteration count)
3. A summary card row shows final values, window count, max WR iterations, and
   whether every window converged.

## Tunable parameters

| Key | Meaning |
|-----|---------|
| `frequency`, `amplitude` | circuit voltage source `Bsrc` |
| `R_series`, `L_series` | series coupling impedance between source EMF and port. `0,0` = one-way (port voltage pinned to source, WR is trivial). Nonzero → **two-way coupling**: `V(p)` reacts to the interface current and WR genuinely iterates. Keep `L_series` ≲ `L_FEM` for numerical stability. |
| `L_ROM`, `R_ROM` | reduced-order model used in the Xyce `Biface` interface condition |
| `L_FEM`, `R_FEM` | "true" field/FEM parameters (may differ from the ROM) |
| `nonlin_model` | distributed-device nonlinearity: `0` = linear, `1` = magnetic saturation (flux `λ(I)=L_FEM·I_sat·atan(I/I_sat)`, so `L(I)=L_FEM/(1+(I/I_sat)²)` drops as the core saturates). Solved by Newton per interval; the Xyce ROM stays linear, so it is a genuine field/ROM mismatch the WR loop must resolve. |
| `I_sat` | saturation current scale (A) for `nonlin_model=1`. `I_sat → ∞` recovers the linear field exactly. Set near/below the operating current to see saturation bite. |
| `N_periods` | number of source periods simulated |
| `N_field_steps_per_source_period` | field/WR windows per source period |
| `N_field_eval_intervals` | FEM evaluation intervals per window (`1` = single ramp) |
| `N_xyce_coupling_intervals` | fixed coupling-grid resolution (print cadence + resample) |
| `WRmaxSteps`, `WR_tolerance` | waveform-relaxation iteration cap and tolerance |
| `wr_convergence_method` | WR convergence metric: `0` = waveform L1 of the field current (this code), `1` = terminal-scalar metric (reference `CoSimulation_WR.cpp`) |

## Circuit visualizer

The top of the results panel shows a live schematic of `wr_circuit.cir`, loaded
automatically on page open (and re-readable with **Refresh circuit** after you
edit the netlist). The backend parses the SPICE deck (joining `+` continuations,
skipping the title/comment lines, keeping `.` directives) and renders a
node-link schematic:

It is drawn as a **traditional ladder schematic**: non-ground nodes sit on a top
line, a **ground rail** runs along the bottom, and every conductive branch is a
component on either a **vertical leg** (node → ground) or a **horizontal top
segment** (node → node), with orthogonal wires.

- **Component symbols** by type: independent voltage source = circle with a sine
  `~`; a 0 V source = ammeter (circle `A`); current source = circle with arrow;
  behavioral/dependent `V=`/`I=` sources = **diamonds** (sine / arrow inside);
  R/L/C = rectangle / plates. Each carries its element name.
- **PWL signal sources** are *not* branches in the loop — they sit **below the
  rail** as small waveform sources on isolated reference nodes, with **dotted
  "reads" arrows** pointing into every element whose expression uses `V(node)`.
- Below the schematic: a table of every element (name, type, nodes, value/
  expression) and a collapsible view of the raw netlist + directives.

For the toy this draws the loop `p`–`nx`–`0`: `Bsrc` (dependent voltage source,
left leg), `Vmeas` (ammeter, top), `Biface` (dependent current source, right
leg), ground rail at the bottom — plus the two PWL signal sources `VFprev`
(`vfprev`) and `VIprev` (`iprev`) below, feeding `Biface` via dotted reads-arrows.
This makes explicit that the field-feedback waveforms are inputs to the interface
condition, not parts of the circuit.

## Convergence study (parameter sweep)

The **Convergence study** panel (bottom of the controls) runs the solver once
per value of one chosen parameter while holding all other fields at their
current values, then plots metrics vs the swept parameter. Use it to find which
parameters slow WR convergence and how cost scales.

1. Pick the **sweep parameter** (any numeric parameter).
2. Set **min, max, steps** and **spacing** (`linear` or `log`; `log` needs
   strictly positive min/max).
3. Click **Run sweep**.

It produces a 2×2 figure:
- **Convergence speed** — max & mean WR iterations per window
- **Cost** — total Xyce solves and wall-clock solver time
- **WR accuracy** — worst WR relative error (log axis)
- **Final interface values** — final `I_field` and `V_field`

plus a per-point table (rows that failed or did not fully converge are red).

### Headless sweep (no browser)

The same study runs from the command line, writing a PNG:

```sh
python3 sim_ui.py sweep --param R_series --min 0 --max 0.05 --steps 6
python3 sim_ui.py sweep --param L_FEM --min 1e-7 --max 1e-5 --steps 8 --scale log --out lfem.png
```

`--param` accepts any numeric config key; `--scale` is `linear` (default) or
`log`; `--out` sets the plot file (default `sweep.png`). The bare
`python3 sim_ui.py` (or `python3 sim_ui.py serve`) still starts the web UI.

> Note: a sweep overwrites `sim_config.txt` with each point's values; the file
> is left holding the last swept value when the sweep finishes.

## Config file

`sim_config.txt` is plain `key = value` (also accepts `key value`; `#`
comments; unknown keys ignored; missing keys keep defaults). The solver reads it
at startup. You can run the solver directly with a custom config:

```sh
./main my_config.txt     # defaults to sim_config.txt if omitted
```

## Things to study

- **Enable two-way coupling** with `R_series`/`L_series` > 0, then watch the WR
  iteration count and `Xyce solves` / solver time grow as the coupling strengthens.
  With `R_series=L_series=0` the problem is one-way and converges in ~2 iterations
  regardless of any other parameter (see `refactor_docs` §"Trivial WR convergence").
- Set `L_FEM`/`R_FEM` ≠ `L_ROM`/`R_ROM` and watch the WR error grow and the
  iteration count rise (the ROM no longer matches the field).
- Raise `N_field_eval_intervals` above `1` to switch the field current from a
  single ramp to a multi-rate piecewise-linear waveform.
- Lower `WR_tolerance` and observe more WR iterations per window.
- Increase `N_periods` to see several source cycles.
- Switch `wr_convergence_method` between the waveform-L1 and terminal-scalar
  metrics and compare the reported WR error / iteration counts for the same run.
- Turn on `nonlin_model=1` (magnetic saturation). Note it only changes the WR
  **rate** under two-way coupling (`R_series`/`L_series` > 0) — in the one-way
  case the field is a leaf with a fixed input, so it still converges in ~2
  iterations regardless of `I_sat`. With `R_series>0`, sweep `I_sat` (log) down
  toward the operating current: once saturation bites (`L(I)` deviates from
  `L_FEM`) the iteration count climbs (e.g. at `R_series=0.02`: `I_sat=100`→10,
  `=10`→22, `=1`→capped). This makes the convergence rate operating-point /
  amplitude dependent, unlike the linear field.
- Use the **Convergence study** sweep (above) to map any of these effects
  automatically — e.g. sweep `R_series` from 0 to 0.05 and watch WR iterations
  jump from ~2 to ~7, or sweep `L_FEM` (log) to see the WR error rise as the ROM
  drifts from the field.
