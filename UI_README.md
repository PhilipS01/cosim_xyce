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
   - **Interface voltage** vs time (circuit `V(p)` vs field `V_field`)
   - **Interface current** vs time (circuit `I(Vmeas)` vs field `I_field`)
   - **WR convergence** per window (final rel. error + iteration count)

   With `reference_overlay` on, the monolithic reference is re-solved and drawn on
   the two interface plots as a dotted black curve.

   A fifth plot, **Field-circuit interface defect (within windows)**, appears when
   both waveforms are present. It is the inner-window transmission residual,
   sample by sample rather than only at the window terminals:

   ```
   dV(t) = V_field(t) - V(p)_circuit(t)
   dI(t) = I_field(t) - I(Vmeas)_circuit(t)
   ```

   The field waveform has only `N_field_eval_intervals` nodes per window, so it is
   **linearly interpolated onto the circuit grid** and the defect is reported at
   every circuit sample. That is deliberate: the circuit is driven by the field's
   PWL carriers, which Xyce itself reads as a linear interpolant between field
   nodes, so this is the defect the circuit actually saw — field-grid
   reconstruction error included. `Circuit_solution.prn` carries Xyce's *raw
   adaptive* time points (thousands per run, not the `dt_print` coupling grid), and
   the field nodes are a subset of them — the PWL carriers put a breakpoint at each
   one, which Xyce must step onto — so the interpolation only fills in *between*
   field nodes and never extrapolates at them.

   Values are **signed** on linear twin axes, because the defect typically ramps
   inside a window and resets at the seam and that sawtooth is the content.

   **Only one channel is a real cross-solver defect.** `Field_waveform_solution.prn`
   carries the *circuit's own* waveform, resampled onto the field grid, in the
   other column — so with `coupling_mode = 0` (voltage-driven) `dI` is the
   transmission defect and `dV` is the field-grid **reconstruction error** of
   `V(p)`; with `coupling_mode = 1` the roles swap. The legend names which is
   which for the run you just did. Both are worth reading: the reconstruction
   channel is exactly what changes when you trade WR windows against field
   evaluations.

   Summary scalars `mean_V_defect` / `mean_I_defect` (and `max_*`) are the
   mean/max of `|d·|` over every circuit sample; they also land in the sweep table
   and CSV. Both figures and the defect summary cards carry a hover **?** spelling
   out the interpolation, which channel is the real cross-solver defect for the
   coupling direction in use, and how they differ from **worst f-c defect** (which
   the solver measures only at the window terminals, and normalises only above 0.1).

   Hovering a results plot reveals a **TikZ** button that downloads that one plot
   as `pgfplots` LaTeX (`interface_voltage.tex`, `interface_current.tex`,
   `wr_convergence.tex`, `probe_voltages.tex`, `probe_currents.tex`). The file is
   a bare `tikzpicture` — no `\documentclass` — so it drops into a thesis with
   `\input{...}`; the two preamble lines it needs (`\usepackage{pgfplots}`,
   `\pgfplotsset{compat=1.18}`) and the run's key parameters are written as
   comments at the top. It re-renders from the last run's parsed data (no
   re-solve), and each series is thinned to ~2000 points by min/max decimation so
   narrow spikes survive while `pdflatex` stays fast (`_PGF_MAX_PTS` in
   `sim_ui.py`).
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
| `interface_form` | WR interface stamping: `0` = Thevenin (V source), `1` = Norton (I source) — algebraic duals, same terminal fixpoint |
| `precondition` | interface preconditioner: `1` = on (default, the matched-secant ROM impedance `Z = R_ROM + L_ROM/dt` in `Bfield` — an optimized/Robin transmission that accelerates the WR contraction); `0` = **classical Gauss–Seidel (Dirichlet–Neumann) WR** — drop the `Z` correction, drive the port with a *pure* source of the field's own response (`I_field` current source when voltage-driven, `V_field` voltage source when current-driven). Same fixpoint (the `Z`-term vanishes at convergence), slower contraction, may not converge for stiff coupling. Disables `interface_form`/`use_t_floor`/`R_ROM`/`L_ROM`. **Well-posedness:** an inductive port (series `L` to `p`) needs the voltage source → use `coupling_mode=1` (a pure current source in series with `L` is degenerate and Xyce aborts at `t=0`); a capacitive port is the dual (voltage-driven) |
| `use_t_floor` | secant denominator `dt` in `Z = R_ROM + L_ROM/dt`: `1` = floored `MAX(dt, t_floor)` (default, guards the window-start `1/0`); `0` = bare `dt` (test); `2` = **constant** `t_floor` → fixed interface impedance `Z = R_ROM + L_ROM/t_floor` (a constant Robin/optimized-transmission coefficient instead of the growing-admittance accumulated secant; same fixpoint, different WR rate) |
| `interface_consistency` | **study knob** for the Thevenin interface (`precondition=1`, `interface_form=0`), **both coupling directions**; ignored otherwise; FEM unchanged. `0` = **accumulated** (default) — one secant denominator `dt` for both currents, so the two `Lrom/dt` terms cancel at convergence → true fixpoint; `1` = **BDF-1/BE (naive)** — split it so the live-iterate circuit-current term uses the circuit's own step `dt_C` while the lagged term keeps `time−t_abs_start` (`= dtf`, still per `use_t_floor`). Same emission both directions; the shifted observable differs: **voltage-driven** (`i_prev=I_field`) `V_C = V_F + Rrom·(I_C−I_F) + Lrom·(I_C−I0)/dt_C − Lrom·(I_F−I0)/dtf` → terminal **current** shifts (`I_circuit` vs `I_field`); **current-driven** (`i_prev=` lagged `I(Vmeas)`, base `V_field`) → at the fixpoint `V_C = V_field + Lrom·(I_C−I0)·(1/dt_C − 1/dtf)` → terminal **voltage** shifts. `dt_C` is realized with **real devices**: `R_ROM`/`L_ROM` stamped as real Xyce devices on the port branch (`IC=I0`), so `dt_C` is Xyce's actual adaptive step and the inductive term is the genuine BDF-1/BE backward difference; `Bfield`'s own inductive term reads the field-grid BDF-1 derivative `di_F/dt` from `didt_field_k.pwl` directly — exact for >1 FEM eval/window, no anchored secant |
| `t_floor_frac` | `t_floor = t_floor_frac · t_window` (window-scaled); the secant-denominator floor (`use_t_floor=1`) **or** the constant denominator (`use_t_floor=2`) |
| `seam_average` | window-seam handoff: `0` = one-sided (V←circuit, I←field), `1` = midpoint |
| `validation_mode` | `0` = WR co-sim (default). `1` = **monolithic reference**: replace the behavioral `Bfield` with the *true* field as real Xyce devices (`R_FEM` + `L_FEM` in series on the port branch) and solve the whole circuit as one transient over `[0, t_end]` — no WR loop, no field solver, no coupling. Gives a reference to validate the coupled run against; the WR/coupling/secant knobs are disabled and the convergence plot is empty. Linear field only (saturation ignored) |
| `xyce_max_step` | ceiling (s) on Xyce's **internal adaptive timestep**, emitted as the 4th positional field of the per-window `.tran` (`.TRAN <initial step> <final time> [<start time> [<step ceiling>]] [NOOP] [UIC]`, Xyce RG 2.1.38). Distinct from `N_xyce_samples`, which only sets how often a solved point is *printed*. `0` = off — but off is **not unbounded**: Xyce then applies its own default ceiling of `(t_stop − t_start)/10 = t_window/10`, tightened where breakpoints require ≥10 steps between them. WR windows only; the monolithic solve keeps its own `dt_print` ceiling (a value set here is reported as ignored) |
| `xyce_integration_method` | implicit time-integration scheme via `.OPTIONS TIMEINT` (Xyce UG table 7-3). `0` = trap, Xyce's default variable-order trapezoid (emits nothing); `1` = **Backward-Euler** (`METHOD=trap MAXORD=1`) — 1st order, matches the FEM dummy solver's own BDF-1 stepping; `2` = trap only (`METHOD=trap MINORD=2`); `3` = Gear (`METHOD=gear`); `4` = Gear2 only (`METHOD=gear MINORD=2`). Written into the netlist head, so unlike `xyce_max_step` it **also applies to the monolithic reference** |
| `reference_overlay` | **studio-only** (never reaches the solver — `write_config` round-trips it as a `# ui:` comment). `0` = off, `1` = re-solve the same circuit with `validation_mode=1` in a temp directory and draw its interface voltage/current on the **Interface voltage** / **Interface current** plots as a dotted black curve. Costs one extra full solve per run; a failed reference is dropped with a note in the log. The reference is a linear field, so it is not the truth the coupled run converges to when `nonlin_model=1` |

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
- **WR accuracy** — worst **and mean** WR relative error (log axis)
- **Final interface values** — final `I_field` and `V_field`

`WR_error.txt` holds the *final* relative error each window converged to, one
entry per window. **Worst** is the window that ended furthest from tolerance —
sensitive to one hard window (a switching edge, say); **mean** is the typical
window, which is what moves when accuracy changes across the whole run. Reading
them together separates "one window struggles" from "everything got worse".

A 2-parameter grid sweep draws these as **five heatmaps** (2×3, last cell blank):
max WR iterations, total Xyce solves, final `I_field`, and the accuracy pair —
worst and mean WR relative error, each on its own log colour scale (they differ
by orders of magnitude, so a shared scale would flatten the mean panel).

Every sweep also draws a **Field-circuit interface defect** figure: the
inner-window defect above, collapsed to one number per run. Its shape follows the
sweep —

- **2-parameter grid** → two heatmaps (`mean |dV|`, `mean |dI|`) on a log colour
  scale, one per channel;
- **single / parallel / >2-parameter grid** → two log-y panels of **mean and max**
  against the swept parameter (or the run index for a >2-D grid).

Kept as its own figure rather than folded into the metric grid, because the two
channels are not interchangeable — one is the cross-solver defect and the other
the reconstruction error, and which is which flips with `coupling_mode`, so they
are read as a pair.

Exact zeros can't sit on a log scale. On the heatmaps such a cell is drawn grey
(a `set_under` colour, with the count in the panel title) rather than floored,
which would invent decades of dynamic range that aren't there; the TikZ export
has no equivalent, so it emits those cells as holes and says so in its header. On
the line panels a zero is floored to the smallest positive value in that panel so
the point still plots.

**Heatmap axes follow the sweep's Spacing.** A parameter swept `log` gets a
log-scaled axis, so its cells come out uniform instead of crowding at the small
end (four decades of `WR_tolerance` used to collapse into one unreadable band at
the bottom). Ticks sit on the sampled values either way. This applies to the
five metric heatmaps and the total-WR-iterations heatmap, in the PNG and the TikZ
export alike. The *per-window* colormap is unaffected: its x is a column index —
one column per run — not the parameter value, so its cells are already even.

plus a WR-iterations-per-window colormap and a per-point table (rows that
failed / didn't fully converge are red).

While a sweep runs, a **progress bar** under the buttons shows points completed,
elapsed time, a rough ETA, the failed count, and the number of parallel jobs. The
sweep is a single long request, so the studio polls `/sweep_progress` on a second
connection rather than streaming; the ETA is linear in the mean time per
completed point, which is honest when every point costs about the same and drifts
when the swept parameter itself changes the cost (more windows, more WR
iterations).

A **2-parameter grid sweep** adds one more figure: a heatmap of the **total WR
iterations** (summed over every window) with the first swept parameter on x and
the second on y — the "how expensive is this corner of the parameter space" view.
The per-window colormap above resolves iterations per window but collapses the
two parameters to a flat run index, so the 2-D structure is only visible here.
Cells are labelled with their count while the grid is at most 144 cells; a blank
cell is a combination whose solve failed (a window hit `WRmaxSteps` and the run
aborted).

### Parameter row types

Each row in the panel is one of three kinds, chosen in its **Type** column:

| type | meaning | adds a dimension? |
| --- | --- | --- |
| **Range** | independently swept, Min → Max over Steps (linear or log) | yes — the sweep mode applies to these |
| **% of…** | tracks another swept parameter at a fixed percentage | no |
| **Expression** | computed per point from the other parameters | no |

### Sweeping the source

The circuit side lives in `circuit_spec.txt`, not in `sim_config.txt`, so its
numbers were out of reach of the sweep. When the spec holds **exactly one
source**, its fields now join the parameter dropdown:

| source | fields |
| --- | --- |
| `VSIN` / `ISIN` | `src_amp`, `src_freq` |
| `VDC` / `IDC` | `src_val` |
| `VPULSE` / `IPULSE` | `src_v1`, `src_v2`, `src_td`, `src_tr` |
| `VPWM` / `IPWM` | `src_v1`, `src_v2`, `src_freq`, `src_duty` |

They are labelled with the source's own name (`Bemf frequency (Hz)`), sweep like
any other parameter — range, `% of…`, or **Expression** — and can be named inside
an expression or a reference line. Each point rewrites that one number in the
spec line, leaving comments, spacing and every other element untouched.

Requires exactly one source: with two, "the source's frequency" names nothing and
the sweep would have to guess which line to rewrite, so the options disappear from
the dropdown and an explicit `src_*` key is refused. `{V,I}PWL` counts as a source
but exposes no fields — its `t1 v1 t2 v2 …` list has no stable name for "the third
number". The dropdown refreshes whenever the schematic does, so editing the spec
adds or removes the options live.

An **Expression** row is evaluated once per sweep point. Available names:

- any config parameter — its value *at this point*, so `N_field_windows` inside
  the expression is the value this point is being run at;
- the source fields above (`src_freq`, `src_amp`, …) when the circuit has one source;
- `base_<name>` — that parameter's **pre-sweep** value (whatever the Properties
  panel holds), so you needn't hardcode where the sweep started;
- your circuit-spec `R`/`L`/`C` element names (`Rs`, `Ls`, …);
- `abs`, `min`, `max`, `sqrt`, `pi`.

Integer parameters are rounded. A value below the control's own minimum (or a
name that doesn't resolve, or a division by zero) rejects the whole sweep with a
message naming the row — before anything solves.

#### Holding the Xyce resolution constant

The printed step is

```
dt_print = t_end / (N_field_windows × N_xyce_samples)
```

so sweeping `N_field_windows` silently changes the Xyce sample spacing unless
`N_xyce_samples` falls with it — a `1/x` relationship no percentage can express.
Add `N_xyce_samples` as an **Expression** row:

```
base_N_field_windows*base_N_xyce_samples/N_field_windows
```

That product is the total number of printed samples over the whole run, so the
literal form `20000/N_field_windows` does the same thing. With it in place you
can grid **WR windows (total)** against **Field evals per window** and know the
only things changing are the two you swept:

| `N_field_windows` | `N_field_eval_intervals` | `N_xyce_samples` | `dt_print` |
| --- | --- | --- | --- |
| 25 | 1 · 2 · 4 | 800 | 1e-07 |
| 50 | 1 · 2 · 4 | 400 | 1e-07 |
| 100 | 1 · 2 · 4 | 200 | 1e-07 |

Derived rows appear in the per-point table and the exported CSV like any other
parameter, so the value actually used is always on the record.

Both sweep figures carry the same **TikZ** button as the results plots. Because
the on-screen figures are composites, clicking TikZ on one opens a **panel
picker**: *All N panels (one file)*, or any single panel on its own. The
all-panels file holds one independent `tikzpicture` per panel and `\input`s them
in sequence; a single-panel file carries just that picture, is named after it
(`sweep_<params>_grid_mean_wr_rel_error.tex`), and heads with which panel of
which figure it is — so several panels of the same sweep can sit side by side in
a thesis without colliding. Single-picture figures (the iterations colormaps, the
results plots) download straight away with no menu. The two twin-axis panels
are split so each quantity gets its own axis (six panels: WR iterations, Xyce
solves, solver time, WR accuracy, final `I_field`, final `V_field`); a 2-parameter
grid sweep exports its five heatmaps instead, and the iterations colormap exports
on its own. Error heatmaps carry `log10(value)` as the colour meta with the
colourbar relabelled in powers of ten — pgfplots has no logarithmic colour scale,
and without that a heatmap spanning decades exports as one flat colour. Reference lines are drawn from whatever is in the reference-line box
at export time. Heatmaps use pgfplots' built-in `viridis` for every panel, where
the PNG varies the colormap.

Headless (writes a PNG):

```sh
python3 sim_ui.py sweep --param L_FEM --min 1e-7 --max 1e-5 --steps 8 --scale log --out lfem.png
python3 sim_ui.py sweep --param WR_tolerance --min 1e-4 --max 1e-2 --steps 6 --scale log
```

`--param` accepts any numeric config key; `--scale` is `linear` (default) or
`log` (needs strictly positive min/max); `--out` sets the plot file (default
`sweep.png`).

> Note: every sweep point solves in its own temporary directory (that is what lets
> points run concurrently), so a sweep leaves `sim_config.txt` and the other
> working-tree files exactly as it found them.

## A-priori convergence estimate (port impedance `x_P`)

The **A-priori estimate** fold in the *Run* section answers, before any transient
runs, whether the WR iteration contracts for this circuit and this ROM.

**`x_P(f)`** is the impedance the field sees looking *into* the circuit at the
interface, with the circuit's own independent sources zeroed — the port voltage
response to a 1 A injection at the port. It comes from a Xyce `.AC` sweep of a
probe deck built from your `circuit_spec.txt`:

| spec element | in the probe deck |
| --- | --- |
| `R` / `L` / `C` | kept (without the transient deck's `IC=`) |
| any V source (`VSIN`, `VDC`, `VPULSE`, `VPWM`, `VPWL`) | zeroed → a 0 V source, i.e. a short that keeps the node names |
| any I source (`ISIN`, `IDC`, …) | removed → an open circuit |
| `SW` | frozen to `Ron`/`Roff` at **switch t** — an AC analysis needs a time-invariant circuit |
| the whole field branch (`R_ROM`/`L_ROM`/`Bfield`/`Vmeas`) | absent; replaced by `I_xp_inj 0 p AC 1 0` |

The port voltage is read from the `.prn` **by node name** (`VR(p)`/`VI(p)`), never
by column position.

**`ρ(f)`** is the WR contraction factor built on top of it,

```
β(f) = 1 / (R_ROM + j2πf·L_ROM)      ROM (preconditioner) admittance
Y(f) = 1 / (R_FEM + j2πf·L_FEM)      true field admittance
ρ(f) = (1 + β·x_P)⁻¹ (β − Y) x_P
```

`|ρ| < 1` over the band the circuit actually excites is the a-priori statement
that WR contracts there; `ρ ≡ 0` when the ROM matches the field exactly. The panel
prints `max |ρ|`, the band where `|ρ| ≥ 1`, and every `|ρ| = 1` crossing.

- **Compute x_P (Xyce .AC)** runs the sweep (`f start`/`f stop`, `dec`/`oct`/`lin`,
  points) and then `ρ` on top of it.
- **Update ρ only** recomputes `ρ` from the *stored* `x_P` using the current
  `R_ROM`/`L_ROM`/`R_FEM`/`L_FEM` — no Xyce solve. `x_P` depends only on the
  circuit, so re-tuning the ROM is free.
- Both sweeps are written to the CSV paths in the panel (blank = skip):
  `f_Hz, Re_xP_ohm, Im_xP_ohm, abs_xP_ohm` and the same plus
  `Re_rho, Im_rho, abs_rho`.
- Two plots appear in *Results*: `|x_P(f)|` on a log-log axis, and `|x_P|` with
  `|ρ|` against the same frequency axis (`|ρ|` on the right-hand scale, with the
  `|ρ| = 1` level and its crossings marked). Both carry the same hover **TikZ**
  button as the other results plots, exporting `port_impedance.tex` and
  `wr_contraction.tex` — bare `tikzpicture`s re-rendered from the stored sweep (no
  Xyce re-run), with the sweep settings and the ROM/field values recorded as
  comments. `wr_contraction.tex` reproduces the twin-axis figure: `|x_P|` in ohms
  on the left, `|ρ|` on the right, the dashed `|ρ| = 1` rule and its crossing
  markers. AC sweeps are a few hundred points, so these two are **not** decimated.

**Sanity checks** run automatically and are printed under the buttons:

- *low f* — `|x_P(f_start)|` is compared against the network's DC resistance
  (for a series-R port, `R_s`), solved independently of Xyce by shorting every
  inductor and zeroed source and opening every capacitor. A mismatch usually means
  `f_start` is too high for `2πf·L ≪ R` to still hold.
- *high f* — with a capacitor directly across the port (as in preset P1), the cap
  shorts the port, so `|x_P|` must roll off as `1/f`. Tested against that
  capacitor's *own* `|Z_C| = 1/(2πf·C)` at `f_stop` — everything else at the port
  is in parallel with it, so `|x_P|` can only sit at or below `|Z_C|` once it
  dominates. Not tested as `|x_P| ≈ 0`: at a finite `f_stop` a `1/(2πf·C)` tail is
  still a visible number. A sweep that stops below the port resonance is reported
  **INCONCLUSIVE** (raise `f_stop`), not as a failure — the band is your choice.
  The top-decade log-log slope is printed alongside (`−1` once the cap dominates).

### Outside the studio

The same two steps are standalone CLIs, deliberately split so the `ρ` sweep can be
re-run with new ROM values without touching Xyce:

```sh
python3 scripts/xp_extract.py --fstart 1 --fstop 1e7 --points 200 \
        --out results/xp.csv --plot results/xp.png --deck results/xp_probe.cir
python3 scripts/rho_contraction.py --xp results/xp.csv --plot results/rho.png \
        --r-rom 0.9 --l-rom 0.9e-3        # ROM defaults come from sim_config.txt
```

`xp_extract.py` exits non-zero if a sanity check fails (`--no-check` to skip).
A bare current source on the port (preset P2) has no `x_P`: removing it leaves the
port open-circuit, and the script says so instead of letting Xyce report an empty
matrix.

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
- **Field resolution** — raise `N_field_eval_intervals` above `1` to let the
  dummy field follow the curve within a window instead of taking a single
  backward-Euler step across it (one field solve per interval). At `1` the
  window reduces to the straight carried-start → window-end ramp.
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
  to the port: an inductive port (series `L` to `p`, e.g. preset P5) needs
  `coupling_mode=1` so the circuit sees a voltage source; a pure current source
  in series with `L` is degenerate.
- **Accumulated vs naive BDF-1/BE interface** — with the Thevenin interface
  (`precondition=1`, `interface_form=0`, either coupling direction), flip
  `interface_consistency` from `0` (accumulated) to `1` (BDF-1/BE, naive: `R_ROM`
  and `L_ROM` as real devices on the port). The accumulated scheme shares one
  secant denominator so the `Lrom` terms cancel at convergence (true fixpoint);
  the naive one gives the live-iterate and lagged terms different denominators, so
  they no longer cancel and the fixpoint **shifts** — the terminal **current** (`I_circuit` vs
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
