# WR Co-Simulation Studio

Browser UI to study the waveform-relaxation (WR) field/circuit co-simulation:
tune solver/coupling parameters, author the circuit, and view the resulting
waveforms + WR convergence.

## Requirements

**Core** (run simulations + view plots):

- **Python 3** with **numpy** + **matplotlib** — `pip install numpy matplotlib`.
  Everything else is standard library.
- A **C++17 compiler** + **make** (the `Makefile` uses `clang++`). The UI builds
  `./main` automatically via `make` if it is missing.
- **Xyce** on `PATH` — the circuit solver the backend shells out to.

**Optional** (circuit schematic + `.tex`/`.pdf` export only):

- **Node.js** + `npm install` (pulls `elkjs`) — schematic layout.
- **pdflatex** + `circuitikz` (TeX Live / MiKTeX) — renders the schematic PNG and
  the PDF export.

Without the optional tools the studio still runs and all simulations/plots work;
only the inline schematic and `.tex`/`.pdf` export are unavailable (they degrade
gracefully).

## Run

```sh
python3 sim_ui.py                     # web UI at http://127.0.0.1:8000
python3 sim_ui.py serve --port 9000   # custom port
python3 sim_ui.py sweep --param L_FEM --min 1e-7 --max 1e-5 --steps 8 --scale log   # headless study
```

## What it does

1. You set properties, author the circuit, and click **Run simulation**.
2. The backend writes `sim_config.txt`, `circuit_spec.txt`, `probes.txt`, runs
   `./main`, parses the `.prn` outputs, and renders three plots server-side as
   PNGs (fully offline — no internet/CDN needed):
   - **Port / field voltage** vs time
   - **Interface current** (circuit `I(Vmeas)` vs field `I_field`) vs time
   - **WR convergence** per window (final rel. error + iteration count)
3. A summary card row shows final values, window count, max WR iterations, and
   whether every window converged.

## Circuit authoring

The circuit side is authored as a **custom node-graph** in the "Circuit netlist"
text box (saved to `circuit_spec.txt`), one element per line:

```
<TYPE> <name> <nodeA> <nodeB> <params...>
```

- **Reserved nodes:** `p` = port (the field attaches here), `0` = ground. Any
  other token is a user node. Do **not** use node `nx` — it is reserved by the WR
  interface.
- **Types:**
  `R/L/C name a b val` ·
  `VSIN/ISIN name a b amp freq` ·
  `VDC/IDC name a b val` ·
  `VPULSE/IPULSE name a b v1 v2 td tr` ·
  `VPWM/IPWM name a b v1 v2 freq duty` (pulse train) ·
  `VPWL/IPWL name a b t1 v1 t2 v2 …` (multi-step) ·
  `SW name a b tclose topen [Ron Roff trise]` (time-gated switch).
- `#`/`*` comments and blank lines are ignored.

Presets **P1–P6** (sine V + RL, bare current source, step/ramp + RL, and three
2-way switch circuits) seed the spec; edit freely from there. **Two-way
coupling** comes from series R/L on the source→port path — author them in the
spec (e.g. `R Rs s cm0 6e-3` / `L Ls cm0 p 1.6e-7`). A source sitting directly on
`p` with no series R/L is the trivial one-way case (WR converges in ~2 iters).

The fixed WR interface (`Vmeas`, the matched-secant `Bfield` source, and the PWL
feedback files `vf_prev_k.pwl` / `i_prev_k.pwl`) is appended automatically after
the user circuit.

## Parameters (Properties section → `sim_config.txt`)

| Key | Meaning |
|-----|---------|
| `L_ROM`, `R_ROM` | reduced-order field model in the Xyce `Bfield` interface (matched secant `Z = R_ROM + L_ROM/dt`) |
| `L_FEM`, `R_FEM` | "true" field/FEM parameters (may differ from the ROM → the WR loop must resolve the mismatch) |
| `nonlin_model` | distributed-device nonlinearity: `0` = linear, `1` = magnetic saturation (`λ(I)=L_FEM·I_sat·atan(I/I_sat)`, so `L(I)=L_FEM/(1+(I/I_sat)²)` drops as the core saturates) |
| `I_sat` | saturation current scale (A) for `nonlin_model=1`; `I_sat→∞` recovers the linear field |
| `t_end` | absolute run duration (s) |
| `N_field_windows` | number of WR windows over `[0, t_end]` (`dt_field = t_end / N_field_windows`) |
| `N_field_eval_intervals` | FEM evaluation intervals per window (`1` = single ramp; higher = multi-rate piecewise-linear) |
| `N_xyce_samples` | Xyce solution sampling — sets the `.tran` print cadence and the interface-PWL resolution (both coupling directions) |
| `WRmaxSteps`, `WR_tolerance` | WR iteration cap and tolerance |
| `wr_convergence_method` | `0` = waveform-L1 of the field current, `1` = terminal-scalar metric |
| `coupling_mode` | `0` = voltage-driven (circuit sets `V(p)`, field returns `I`), `1` = current-driven (circuit sets `I(Vmeas)`, field returns `V_field`) |
| `reconstruct_mode` | field reconstruction within a window: `0` pointwise (secant), `1` linear ramp |
| `interface_form` | WR interface stamping: `0` = Thevenin (V source), `1` = Norton (I source) — algebraic duals, same terminal fixpoint |
| `use_t_floor` | `1` = guard the window-start `1/0` in `Z` with `t_floor`, `0` = bare `dt` (test) |
| `t_floor_frac` | `t_floor = t_floor_frac · t_window` (window-scaled secant-denominator floor) |
| `seam_average` | window-seam handoff: `0` = one-sided (V←circuit, I←field), `1` = midpoint |

Output **probes** (extra Xyce `.print` tokens like `V(a)`, `I(Rr1)`) are written
to `probes.txt` and captured per window into `Probes_solution.prn`.

## Circuit visualizer (optional)

The results panel shows a live schematic of `wr_circuit.cir` (regenerated from
the current config via `./main emit`), rendered **ELK layout → circuitikz →
pdflatex → PNG** — the same render as the export. It needs Node + `elkjs` and
`pdflatex` (see **Requirements**); without them the schematic is simply hidden.
**Export .tex/.pdf** produces the standalone circuitikz source plus a compiled
`circuit.pdf` for thesis figures.

## Convergence study (parameter sweep)

The **Convergence study** panel runs the solver once per value of one chosen
numeric parameter (all other fields held at their current values), then plots
metrics vs the swept parameter in a 2×2 figure:

- **Convergence speed** — max & mean WR iterations per window
- **Cost** — total Xyce solves and wall-clock solver time
- **WR accuracy** — worst WR relative error (log axis)
- **Final interface values** — final `I_field` and `V_field`

plus a per-point table (rows that failed / didn't fully converge are red).

Headless (writes a PNG):

```sh
python3 sim_ui.py sweep --param L_FEM --min 1e-7 --max 1e-5 --steps 8 --scale log --out lfem.png
python3 sim_ui.py sweep --param WR_tolerance --min 1e-4 --max 1e-2 --steps 6 --scale log
```

`--param` accepts any numeric config key; `--scale` is `linear` (default) or
`log` (needs strictly positive min/max); `--out` sets the plot file (default
`sweep.png`).

> Note: a sweep overwrites `sim_config.txt` with each point's values and leaves
> it holding the last swept value when it finishes.

## Config file

`sim_config.txt` is plain `key = value` (also accepts `key value`; `#` comments;
unknown keys ignored; missing keys keep defaults). Run the solver directly:

```sh
./main my_config.txt     # run WR co-sim (defaults to sim_config.txt if omitted)
./main emit [config]     # only regenerate wr_circuit.cir (no solve)
```

## Things to study

- **Two-way coupling** — author series R/L on the source→port path in the spec,
  then watch the WR iteration count and `Xyce solves` / solver time grow as the
  coupling strengthens (vs the trivial one-way `R=L=0` ~2-iteration case).
- **ROM/field mismatch** — set `L_FEM`/`R_FEM` ≠ `L_ROM`/`R_ROM` and watch the WR
  error grow and the iteration count rise (the ROM no longer matches the field).
- **Coupling direction** — flip `coupling_mode`; current-driven suits voltage-source-like
  ROMs.
- **Field reconstruction** — raise `N_field_eval_intervals` above `1` and change
  `reconstruct_mode` to alter the field-current waveform and cost.
- **Saturation** — `nonlin_model=1`, then sweep `I_sat` (log) down toward the
  operating current under two-way coupling: once saturation bites the iteration
  count climbs (operating-point / amplitude dependent, unlike the linear field).
- **Interface stamping** — compare `interface_form` Thevenin vs Norton (same
  terminal fixpoint, different interior `V(p)` waveform and timestepping).
- Lower `WR_tolerance` → more WR iterations per window; toggle
  `wr_convergence_method` between waveform-L1 and terminal-scalar and compare.
- **t_floor** – guards the window-start impedance by limiting the secant denominator. Find a sweet-spot. Choice is problem-dependant (source frequency, impedances of field/circuit, etc.)
- and many more dependencies
