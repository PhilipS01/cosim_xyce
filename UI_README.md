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
  `SW name a b tclose topen [Ron Roff trise]` (time-gated switch — emitted as a
  Xyce native Generic Switch: `S` device + `.MODEL SWITCH`).
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
| `reconstruct_mode` | field reconstruction within a window: `0` pointwise (BDF-1: a local backward-Euler step per field node), `1` linear ramp |
| `interface_form` | WR interface stamping: `0` = Thevenin (V source), `1` = Norton (I source) — algebraic duals, same terminal fixpoint |
| `precondition` | interface preconditioner: `1` = on (default, the matched-secant ROM impedance `Z = R_ROM + L_ROM/dt` in `Bfield` — an optimized/Robin transmission that accelerates the WR contraction); `0` = **classical Gauss–Seidel (Dirichlet–Neumann) WR** — drop the `Z` correction, drive the port with a *pure* source of the field's own response (`I_field` current source when voltage-driven, `V_field` voltage source when current-driven). Same fixpoint (the `Z`-term vanishes at convergence), slower contraction, may not converge for stiff coupling. Disables `interface_form`/`use_t_floor`/`R_ROM`/`L_ROM`. **Well-posedness:** an inductive port (series `L` to `p`) needs the voltage source → use `coupling_mode=1` (a pure current source in series with `L` is degenerate and Xyce aborts at `t=0`); a capacitive port is the dual (voltage-driven) |
| `use_t_floor` | secant denominator `dt` in `Z = R_ROM + L_ROM/dt`: `1` = floored `MAX(dt, t_floor)` (default, guards the window-start `1/0`); `0` = bare `dt` (test); `2` = **constant** `t_floor` → fixed interface impedance `Z = R_ROM + L_ROM/t_floor` (a constant Robin/optimized-transmission coefficient instead of the growing-admittance accumulated secant; same fixpoint, different WR rate) |
| `interface_consistency` | **study knob** for the Thevenin interface (`precondition=1`, `interface_form=0`), **both coupling directions**; ignored otherwise; FEM unchanged. `0` = **consistent** (default) — one secant denominator `dt` for both currents, so the two `Lrom/dt` terms cancel at convergence → true fixpoint; `1`/`2`/`3` = **inconsistent** — split it so the live-iterate circuit-current term uses the circuit's own step `dt_C` while the lagged term keeps `time−t_abs_start` (`= dtf`, still per `use_t_floor`). Same emission both directions; the shifted observable differs: **voltage-driven** (`i_prev=I_field`) `V_C = V_F + Rrom·(I_C−I_F) + Lrom·(I_C−I0)/dt_C − Lrom·(I_F−I0)/dtf` → terminal **current** shifts (`I_circuit` vs `I_field`); **current-driven** (`i_prev=` lagged `I(Vmeas)`, base `V_field`) → at the fixpoint `V_C = V_field + Lrom·(I_C−I0)·(1/dt_C − 1/dtf)` → terminal **voltage** shifts. `dt_C` realized as: `1` real `R_ROM`/`L_ROM` Xyce devices on the port (`IC=I0`, exact adaptive step; `Bfield`'s inductive term reads the field-grid BDF-1 derivative `di_F/dt` from `didt_field_k.pwl` directly — exact for >1 FEM eval/window, no anchored secant); `2` `Lrom·DDT(I(Vmeas))` (Xyce supports `DDT`, but this stamping is numerically fragile — tends to step-collapse on the `Bfield` branch and may abort; prefer real RL); `3` over `dt_print = t_window/N_xyce_samples` (a fixed mean step) |
| `t_floor_frac` | `t_floor = t_floor_frac · t_window` (window-scaled); the secant-denominator floor (`use_t_floor=1`) **or** the constant denominator (`use_t_floor=2`) |
| `seam_average` | window-seam handoff: `0` = one-sided (V←circuit, I←field), `1` = midpoint |
| `validation_mode` | `0` = WR co-sim (default). `1` = **monolithic reference**: replace the behavioral `Bfield` with the *true* field as real Xyce devices (`R_FEM` + `L_FEM` in series on the port branch) and solve the whole circuit as one transient over `[0, t_end]` — no WR loop, no field solver, no coupling. Gives a reference to validate the coupled run against; the WR/coupling/secant knobs are disabled and the convergence plot is empty. Linear field only (saturation ignored) |

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

- **Validate the coupling** — set `validation_mode=1` to solve the circuit against
  the *true* field as real `R_FEM`/`L_FEM` devices in one monolithic transient
  (the reference), then flip back to `0` (WR co-sim) and confirm the coupled
  solution reproduces it. Export both runs to CSV / overlay the plots; the
  reference is exact, so any deviation is WR/ROM error. Linear field only.
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
- **Preconditioner vs classical WR** — set `precondition=0` to strip the
  matched-secant ROM impedance and recover plain Gauss–Seidel
  (Dirichlet–Neumann) coupling, then compare WR iteration counts (and whether it
  converges at all) against the default `precondition=1` for the same circuit —
  the fixpoint is identical, only the contraction rate changes. Match the source
  to the port: an inductive port (series `L` to `p`, e.g. presets P1/P3) needs
  `coupling_mode=1` so the circuit sees a voltage source; a pure current source
  in series with `L` is degenerate.
- **Consistent vs inconsistent interface** — with the Thevenin interface
  (`precondition=1`, `interface_form=0`, either coupling direction), flip
  `interface_consistency` from `0` (consistent) to `1`/`2`/`3` (inconsistent: real
  RL / DDT / mean dt). The consistent scheme shares one secant denominator so the
  `Lrom` terms cancel at convergence (true fixpoint); the inconsistent ones give
  the live-iterate and lagged terms different denominators, so they no longer
  cancel and the fixpoint **shifts** — the terminal **current** (`I_circuit` vs
  `I_field`) under voltage-driven, the terminal **voltage** (`V_circuit` vs
  `V_field`) under current-driven. Compare the converged terminal values (and the
  WR iteration count) against `interface_consistency=0` on the same circuit — the
  gap is the cost of the mismatched denominator. The `waveform L1` metric tests
  iteration change (not the cross-solver defect), so it reports "converged" on the
  shifted fixpoint; the `terminal scalar` metric *does* include the field-vs-circuit
  defect and will refuse to converge (then abort the window). Either way the defect
  is reported directly: `WR_error.txt` now carries per-window `relI_FC`/`relV_FC`
  columns (relative field-vs-circuit terminal gap — voltage-driven populates
  `relI_FC`, current-driven `relV_FC`), and the summary card shows **worst
  f-c defect** over all windows. So you can score the inconsistency even
  when the run "converges."
- Lower `WR_tolerance` → more WR iterations per window; toggle
  `wr_convergence_method` between waveform-L1 and terminal-scalar and compare.
- **t_floor** – guards the window-start impedance by limiting the secant denominator. Find a sweet-spot. Choice is problem-dependant (source frequency, impedances of field/circuit, etc.)
- and many more dependencies
