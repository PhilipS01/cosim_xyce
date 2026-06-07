# Prototype: real BDF devices on the WR interface (companion-network coupling)

Branch `feat/wr-bdf-prototype`. Goal: let circuit-side reactive elements (here `Ls`) be
**real Xyce devices** so Xyce's BDF integrator carries their current + dI/dt natively across
the checkpoint/restart window boundary — instead of the behavioral accumulated-secant source,
which cannot host a real inductor (see the note in the `main`-line `wr_circuit.cir`).

## What changed vs the behavioral version

- **Netlist** (`wr_circuit.cir`): all real devices.
  - `Bemf`/`Rs_d`/`Ls_d`: EMF + series coupling, `Ls_d` a real inductor (`ls_branch.inc`,
    generated: inductor for Ls>0, tiny-resistor short for Ls=0 to avoid a degenerate `L=0`).
  - `Rrom_d`/`Lrom_d`: the field ROM as a **real companion R-L branch** (port → ground). This
    is the WR preconditioner, now BDF-integrated so it is consistent with `Ls_d`. It cancels at
    the fixpoint (`i_Ls = i_rom + icorr → I_FEM`), so it only sets the convergence rate.
  - `Bcorr` + `Vcorr`: **deferred-correction** Norton source `icorr = I_FEM^(k-1) − i_rom^(k-1)`,
    exchanged as a PWL file. `Vmeas` is a 0 V interface-current probe.
- **FEM solver**: anchor-secant → **per-step backward-Euler** (history = previous step, not the
  window anchor), so the field integrates the same differential operator as Xyce's BDF.
  Window-start state seeds the BE history (field-side "restart").
- **Master loop**: builds `icorr` each iteration; seeds/carries it across windows.

## Results

- **Headline (works):** Ls = 1.6e-7 (real inductor), `N_field_eval_intervals = 1` →
  **50/50 windows converge**, no step collapse. The *same* Ls on the behavioral `main` version
  **crashes ~window 27** ("Time step too small") — a BDF inductor in series with the near-ideal
  current-source ROM is stiff/inconsistent. So the prototype does what the behavioral version
  cannot.
- **Residual floor confirmed (the predicted limitation):** with anchor-secant FEM
  (`N_field_eval=1`) the WR residual floors (~5e-3 at Ls=0) above the 1e-3 tol. Switching the
  FEM to BE stepping and refining the grid drops it: Ls=0, `N_field_eval=10` →
  **50/50 converge** (floor ~4e-4).
- **Rough edge:** the floor is **config-dependent and non-monotonic** in `N_field_eval`
  (Ls=0 wants ≥10; Ls>0 wants 1; `=100` floors just above tol near the half-period
  zero-crossing, ~window 26). No single setting nails every config to 1e-3.

## Why the floor is config-dependent (open item)

Three grids interact: Xyce adaptive BDF (~581 pts/window), the coupling grid (100), and the FEM
grid (`N_field_eval`). The deferred correction is resampled between them and lagged one iterate.
Tightening robustly likely needs: consistent/finer correction grids, optional under-relaxation
of `icorr`, and/or a higher-order field integrator. Left as future work — the prototype's purpose
(feasibility of real BDF devices + native restart, plus characterizing the floor) is met.

## Repro

```
make
./main                 # default sim_config.txt: Ls=1.6e-7, N_field_eval=1  -> 50/50
# Ls=0 one-way case needs finer FEM:
sed 's/L_series = .*/L_series = 0/; s/N_field_eval_intervals = .*/N_field_eval_intervals = 10/' \
    sim_config.txt > /tmp/cfg && ./main /tmp/cfg   # -> 50/50
```
