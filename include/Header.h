#ifndef HEADER_H
#define HEADER_H


#include <string>
#include <vector>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cstring>
using namespace std;

#define PRINT( X )   cout<< #X << " =  " << X << flush <<  endl


// Tunable simulation parameters. Loaded from a key=value config file (sim_config.txt)
// by LoadConfig(); any key absent in the file keeps the default set here. The UI writes
// this file before invoking the solver.
struct SimConfig
{
    // Reduced-order model (used inside the Xyce Bfield interface condition)
    double L_ROM = 0.9 * 1.6e-7;
    double R_ROM = 0.9 * 5.1e-4;
    // "True" field/FEM parameters (intentionally different from the ROM)
    double L_FEM = 1.6e-7;
    double R_FEM = 5.1e-4;
    // Nonlinearity of the distributed (FEM) device. 0 = linear (constant L_FEM).
    //   1 = magnetic saturation: flux lambda(I) = L_FEM*I_sat*atan(I/I_sat),
    //       so the small-signal inductance L(I) = L_FEM/(1+(I/I_sat)^2) drops as
    //       the core saturates. I_sat is the saturation current scale (A).
    // The Xyce ROM (Bfield) stays linear, so the nonlinearity is a genuine
    // field/ROM model mismatch that the WR iteration must resolve.
    unsigned nonlin_model = 0;
    double I_sat = 100.0;
    // Run duration (absolute end time): the run spans [0, t_end], split into N_field_windows equal
    // WR windows (dt_field = t_end / N_field_windows). The circuit side (source + passives + switches)
    // is authored as a custom node-graph in circuit_spec.txt, not in this config.
    double t_end = 2.0e-2;
    unsigned N_field_windows = 50;
    // Coupling-grid resolutions
    unsigned N_field_eval_intervals = 1;
    // How finely the Xyce solution is sampled: sets the .tran print cadence dt_print = t_window /
    // N_xyce_samples (-> the raw wr_circuit.cir.prn rows) AND the interface PWL resolution that
    // ReadXyceResults resamples V(p)/I(Vmeas) onto (vf_prev_k.pwl / i_prev_k.pwl, both directions).
    unsigned N_xyce_samples = 100;
    // Secant-denominator floor t_floor as a fraction of the window: t_floor = t_floor_frac * t_window.
    // Only the window-start 1/0 guard in Z = Rrom + Lrom/MAX(dt, t_floor) (see use_t_floor); the floor
    // value is accuracy/stability-neutral (floor-sweep), so a fraction of the window is a scale-free knob.
    double t_floor_frac = 0.01;
    // Waveform relaxation
    unsigned WRmaxSteps = 20;
    double WR_tolerance = 1.0e-3;
    // WR convergence metric: 0 = waveform L1 of the field-current (this codebase),
    //                        1 = terminal-scalar metric of the reference CoSimulation_WR.cpp
    unsigned wr_convergence_method = 0;
    // Coupling direction (Dirichlet vs Neumann):
    //   0 = voltage-driven (default): circuit sets V(p), field returns I; Bfield = matched secant.
    //   1 = current-driven: circuit sets I(Vmeas), field returns V_field; Bfield = plain V source.
    // Current-driven suits current-source circuits (no Lrom/t_floor secant -> no window-start spike).
    unsigned coupling_mode = 0;
    // WR interface stamping (orthogonal to coupling_mode): how the field ROM's linearised V-I law
    // (the matched secant, Z = Rrom + Lrom/dt) enters the circuit. Algebraic DUALS -- same TERMINAL
    // fixpoint (the window-boundary V,I the WR metric checks), but the two stampings give a DIFFERENT
    // interior V(p) WAVEFORM inside each window (and different Xyce timestepping):
    //   0 = Thevenin (default): Bfield is a VOLTAGE source V(nx) = V(vfprev) + Z*(I(Vmeas) - V(iprev)).
    //       Pins V(p); the huge window-start Z multiplies only the small iteration change -> well-conditioned.
    //   1 = Norton: the dual, a behavioral CURRENT source I = V(iprev) + (V(nx) - V(vfprev))/Z with
    //       shunt G = 1/Z. G -> 0 at the window start -> V(p) weakly tied -> stiffer, and the interior
    //       V(p) differs (measured up to > the signal amplitude on a 20 kHz current source). It still
    //       converges (no dt-collapse seen up to ~MHz -- Xyce handles the stiffness), just less cleanly.
    //       Provided to compare the two.
    unsigned interface_form = 0;
    // Interface preconditioner (the matched-secant ROM impedance Z = Rrom + Lrom/dt in Bfield):
    //   1 = on (default): the OPTIMIZED-transmission / Robin interface. Bfield carries the secant
    //       correction Z*(I - i_prev) [Thevenin] or (V - vf_prev)/Z [Norton] -- a Newton linearisation
    //       of the field V-I law that accelerates the WR contraction (this is the preconditioner).
    //   0 = off -> CLASSICAL Gauss-Seidel (Dirichlet-Neumann) WR: drop the Z correction entirely; the
    //       circuit port is driven by a PURE source carrying the field's OWN response variable
    //       (I_field as a current source when voltage-driven; V_field as a voltage source when
    //       current-driven). interface_form is ignored (the base must be the field's output, never the
    //       circuit's prior iterate, else the loop decouples). SAME terminal fixpoint as the
    //       preconditioned scheme (the Z-term vanishes at convergence), but plain Gauss-Seidel
    //       contraction -- slower, and may not converge for stiff/strong coupling.
    //       WELL-POSEDNESS: the pure source must match the port impedance. An INDUCTIVE port (series L
    //       to p) demands the voltage source -> use current-driven (coupling_mode=1); a pure current
    //       source in series with L is degenerate (Xyce aborts "failures at time 0"). A CAPACITIVE port
    //       (shunt C on p) is the dual -> voltage-driven. The finite Z of the preconditioner hides this.
    unsigned precondition = 1;
    // Secant-denominator mode for the Bfield impedance Z = Rrom + Lrom/dt. All variants reuse the
    // t_floor .PARAM (= t_floor_frac * t_window):
    //   1 = floored (default): dt -> MAX(time - t_abs_start, t_floor); the accumulated-secant
    //       denominator, guarding the 1/0 at the exact window start (dt=0 -> Lrom/0 singularity).
    //   0 = bare: dt = (time - t_abs_start); Z -> infinity at the window start. Test only.
    //   2 = constant: dt = t_floor (fixed over the window) -> Z = Rrom + Lrom/t_floor is a CONSTANT
    //       interface impedance (a fixed Robin/optimized-transmission coefficient) rather than the
    //       accumulated secant whose admittance grows across the window. Same terminal fixpoint (the
    //       correction vanishes at convergence); only conditioning / WR contraction rate change. The
    //       FEM solver keeps its own accumulated secant (t_acc = t - t_win_start) -- unchanged.
    //       NB: a constant Z is usually a WEAKER preconditioner (rho -> 1) -> slow; a loose Cauchy
    //       WR_tolerance can then stop early below the true fixpoint (tighten tol / add windows).
    unsigned use_t_floor = 1;
    // Interface CONSISTENCY -- ACCUMULATED (default) vs the NAIVE BDF-1/BE split of the secant
    // denominator, for studying whether the mismatched-denominator scheme is really worse. Scope: the
    // Thevenin form (precondition=1, interface_form=0), BOTH coupling directions; ignored otherwise
    // (falls back to the accumulated emission). The FEM solvers are UNCHANGED (they keep their own
    // accumulated secant and still produce vf_prev/i_prev); only the circuit-side Bfield stamping changes.
    // ACCUMULATED (=0) shares ONE denominator dt for both currents in Z*(I(Vmeas) - i_prev), so the two
    // Lrom/dt terms cancel at convergence -> true terminal fixpoint. NAIVE (=1) splits it so the
    // LIVE-iterate current term uses the circuit's own BDF-1/BE timestep dt_C while the LAGGED term keeps
    // time - t_abs_start (= dtf, still honoring use_t_floor). The device lines are identical for both
    // coupling directions (only the PWL contents differ), so the SAME emission serves both, but the
    // shifted observable differs:
    //   voltage-driven (coupling_mode=0): i_prev = I_field (FEM). Accumulated cancels at I_C = I_field;
    //     the split leaves I_C != I_field -> the terminal CURRENT shifts:
    //     V_C = V_field + Rrom*(I_C - I_field) + Lrom*(I_C - I0)/dt_C - Lrom*(I_field - I0)/dtf.
    //   current-driven (coupling_mode=1): i_prev = I(Vmeas)_{k-1} (lagged circuit current), vf_prev =
    //     V_field. Accumulated cancels at I_C^k = I_C^{k-1}; the split leaves, at the fixpoint,
    //     V_C = V_field + Lrom*(I_C - I0)*(1/dt_C - 1/dtf) -> the terminal VOLTAGE shifts (dual observable).
    // I0 = the carried window-start current (.PARAM I0). Either way the Lrom terms no longer cancel ->
    // the fixpoint SHIFTS (the measured effect). Codes:
    //   0 = ACCUMULATED (default): the single-denominator form, emitted byte-for-byte as before.
    //   1 = NAIVE (BDF-1/BE): Rrom + Lrom stamped as REAL Xyce devices in series on the port branch
    //       (carrying I(Vmeas), Lrom IC={I0}); Bfield carries only the lagged correction. dt_C is then
    //       Xyce's actual (adaptive) circuit timestep -- the genuine BDF-1/BE backward difference, not a
    //       reconstruction. Any nonzero code selects this arm.
    unsigned interface_consistency = 0;
    // Window-seam handoff: how the carried terminal values (V0, I0 seeding the next window) are picked
    // from the converged iterate (the two solvers agree to within WR_tolerance there):
    //   0 = one-sided (default): V0 = circuit V(p), I0 = field I. (Restart checkpoint = Xyce state.)
    //   1 = midpoint: V0 = 0.5*(V_circuit + V_field), I0 = 0.5*(I_circuit + I_field). Seam-blend test
    //       (only the carried seeds are averaged; the Xyce restart checkpoint is unchanged).
    unsigned seam_average = 0;
    // Validation mode: replace the behavioral matched-secant Bfield with the TRUE field as REAL
    // Xyce devices (R_FEM + L_FEM in series on the port branch) and solve the whole circuit as ONE
    // monolithic transient over [0, t_end] -- no WR loop, no dummy field solver, no PWL exchange,
    // no windowing. Produces a reference solution to validate the coupled WR run against.
    //   0 = WR co-sim (default), 1 = monolithic reference.
    // Linear field only: a plain Xyce inductor cannot reproduce the saturation flux law, so
    // nonlin_model=1 is warned + ignored here.
    unsigned validation_mode = 0;
    // Ceiling on Xyce's ADAPTIVE internal timestep, in seconds -- emitted as the 4th POSITIONAL field
    // of the per-window .tran line (Xyce RG 2.1.38: .TRAN <initial step> <final time>
    // [<start time> [<step ceiling>]] [NOOP] [UIC]). Bounds the integrator where the waveform is quiet,
    // independently of the print cadence dt_print (= t_window/N_xyce_samples, which only controls how
    // often a solved point is WRITTEN).
    //   <= 0 (default) = unset: the field is omitted and Xyce keeps its OWN default ceiling, which is
    //   (<final time> - <start time>)/10 = dt_window/10 here, auto-tightened where breakpoints demand
    //   >= 10 steps between them (RG 2.1.38 Comments). So "unset" is NOT "unbounded" -- it is
    //   dt_window/10. A user ceiling overrides any internally generated one.
    //   Unset also means restart.inc is byte-for-byte what it was before this knob existed.
    // NB the initial step is min(<initial step>, <step ceiling>, 1/200 of the time to the next
    // breakpoint), so a ceiling below dt_print shrinks the first step of every window too.
    // The .OPTIONS TIMEINT DELMAX route is deliberately NOT used: per RG table 2-5 it is scoped to
    // ERROPTION=1 and merely combines as min(.TRAN ceiling, DELMAX), so .TRAN is the general knob.
    // Scope: the WR windows only. validation_mode=1 keeps its own hardcoded {dt_print} ceiling (see
    // MonolithicValidationSolve) so the monolithic reference stays as finely resolved as the stitched
    // WR run; a value set here is reported as ignored there.
    double xyce_max_step = 0.0;
    // Xyce time-integration method (User Guide, Table 7-3 "Summary of Xyce-supported time integration
    // methods"). Emitted ONCE as .OPTIONS TIMEINT into the netlist head, so unlike xyce_max_step it
    // applies to BOTH the WR run and the monolithic validation solve (they share wr_circuit.cir) --
    // deliberate, so the reference can be run on the same integrator as the coupled run.
    // Per RG table 2-5, METHOD takes trap (or 7) / gear (or 8); MAXORD caps the order the integrator
    // will attempt, MINORD forces it up to that order -- hence the MAXORD=1 / MINORD=2 spellings.
    //   0 = Xyce default (trap: variable-order trapezoid, dynamically mixing BE and trapezoidal,
    //       MAXORD=2 MINORD=1). Emits NOTHING -> netlist unchanged.
    //   1 = Backward-Euler only   (METHOD=trap MAXORD=1)
    //   2 = Trapezoidal only      (METHOD=trap MINORD=2)
    //   3 = Gear                  (METHOD=gear: backward Euler + 2nd-order Gear)
    //   4 = 2nd-order Gear only   (METHOD=gear MINORD=2)
    // Code 1 matches the circuit side to the FEM dummy solver, which already integrates with BDF-1
    // (see FEM_solver_voltage_driven_waveform) -- i.e. both halves on the same first-order method.
    // Unknown codes fall through to "emit nothing", as elsewhere in the emission chains.
    unsigned xyce_integration_method = 0;
};

extern SimConfig g_cfg;

// Parse key=value lines (also accepts "key value"); '#' starts a comment; unknown keys
// are ignored; missing file leaves all defaults. Returns true if the file was opened.
bool LoadConfig(const string& filename);

// User-requested output probes: fully-resolved Xyce print tokens like "V(a)" or "I(Rr1)".
// Populated by LoadProbes() from probes.txt (one token per line; '#'/blank ignored). Appended
// to the netlist's ".print tran" and captured per converged window into Probes_solution.prn.
extern std::vector<std::string> g_probes;
void LoadProbes(const string& filename);

// Append the probe columns (everything after V(p) V(nx) I(Vmeas)) of a Xyce .prn to the probe
// output file, prefixed by a running index + time. n_probes = g_probes.size().
void appendProbeColumns(const string& xyce_prn, FILE* out, size_t n_probes,
                        unsigned long& global_index, bool skip_first_point);


void MasterProcess();

// Validation-mode driver (validation_mode=1): solves the circuit + true field (real R_FEM/L_FEM
// devices) as ONE monolithic Xyce transient over [0, t_end] and writes the standard output files
// (Circuit_solution.prn / Field_*_solution.prn / WR_error.txt header) so the UI plots/CSV work
// unchanged. No WR loop, no dummy field solver. Reference for validating the coupled run.
void MonolithicValidationSolve();

void FEM_solver_voltage_driven_waveform(double I_win_start, unsigned N_field_eval_intervals);

// Current-driven (Neumann) field solver: reads the interface current waveform I(t) from i_prev_k.pwl,
// computes the field voltage V_field(t) = R_FEM*I + L_FEM*dI/dt (+ saturation) pointwise on the
// field-eval grid (local backward Euler), and writes it to vf_prev_k.pwl for the circuit's
// Bfield voltage source. V_field_last_time = the field voltage carried from the previous window end.
void FEM_solver_current_driven_waveform(double I_win_start, double V_field_last_time,
                                        unsigned N_field_eval_intervals);

// Legacy interface (from the reference CoSimulation_WR.cpp): declared for source compatibility but
// NOT defined in this codebase. Here the circuit side is solved by Xyce (RunXyce) and the interface
// condition is stamped as the netlist Bfield source, so neither of these is called.
void CIRCUIT_solver(const double dt_circuit, const unsigned N_dt_circuit_per_dt_field, const double time_field);

double INTERFACE_condition( const bool reset, const double dt_circuit, const double dt_field,
                            double& V_circuit, double& I_circuit_last_time,
                            double& V_field_last_WR_it, double& V_field_last_time, double& I_field_last_WR_it, double& I_field_last_time);

struct Waveform;
struct CircuitWaveform;

double eval_WR_convergence(const Waveform& i_curr, const Waveform& i_prev_iter, const unsigned WR_iteration);

// Terminal-scalar WR convergence metric (port of the reference CoSimulation_WR.cpp).
// Sums the field-vs-circuit mismatch and the iteration-to-iteration change of the
// terminal V and I (relative when |value| > 0.1, absolute otherwise).
double eval_WR_convergence_terminal(
    double V_field, double I_field,
    double V_circuit, double I_circuit,
    double V_field_last_it, double I_field_last_it,
    const unsigned WR_iteration);

void pushOrReplaceDuplicateTime(Waveform& wf, double time, double value);
void pushOrReplaceDuplicateTime(CircuitWaveform& wf, double time, double vp, double vnx, double i);

void writePWLFile(const string& filename, const Waveform& wf);

Waveform resampleWaveformUniform(const Waveform& raw, double t_start, double t_stop, unsigned N_intervals);

void ReadXyceResults(const string& filename, CircuitWaveform& circuit_raw, double t_start, double t_stop, unsigned N_xyce_samples);

void appendCircuitWaveformXyceStyle(FILE* file, const CircuitWaveform& wf, double t_start, unsigned long& global_index, bool skip_first_point);

void appendFieldWaveformXyceStyle(FILE* file, const Waveform& vf, const Waveform& i, double t_abs_start, unsigned long& global_index, bool skip_first_point);

void RunXyce(const string& filename);

// Generates the full circuit netlist (wr_circuit.cir) from the custom node-graph in circuit_spec.txt
// + the fixed WR interface (Vmeas, Bfield) + includes. Called ONCE before the window loop (and in
// emit mode); the topology is window-invariant. The per-window variable quantities stay in
// sim_params.inc/restart.inc. Written inline (no .INCLUDE of the elements) so the UI netlist parser
// (does not follow .INCLUDEs) can draw the circuit.
void WriteCircuitNetlist(const string& filename);

void WriteSimParams(const string& filename, double t_start, double t_stop, double t_abs_start, double i0, double rrom, double lrom, const unsigned N_xyce_samples);

// Generates restart.inc (included by wr_circuit.cir): the window-specific
// .OPTIONS RESTART and .tran line. Window 1 (first_window=true): fresh UIC transient
// from t=0, writes checkpoints (JOB). Window k>1: restart from committed_file (FILE=),
// writes new checkpoints. Time is absolute: tran stop time = t_stop (absolute).
void WriteRestartDirectives(const string& filename, bool first_window,
                            double dt_window, const string& ckpt_out_prefix,
                            const string& committed_file);

// Deletes stale checkpoint candidates <prefix>* before the WR loop (prevents an
// outdated candidate from being picked as the newest one).
void ClearCheckpoints(const string& prefix);

// Finds the checkpoint <prefix>* with the greatest sim-time (= window end) and copies it to
// committed_file as the restart basis for the next time window.
void CommitCheckpoint(const string& prefix, const string& committed_file);

// Writes the two window-end terminal scalars (V, I) to a one-line text file, used to hand the
// field/circuit end values between solver stages (Field.txt / Circuit.txt).
inline void Write_Terminal_results(const char s[80], const double V, const double I)
{
    FILE* file = fopen(s, "w");
    fprintf(file, " %12.5e   % 12.5e", V, I);
    fflush(file);
    fclose(file);
}

// Reads back the two terminal scalars (V, I) written by Write_Terminal_results.
inline void Read_Terminal_results(const char s[80], double& V, double& I)
{
    FILE* file = fopen(s, "r");
    const int l_char = 256;
    char in[l_char];
    char* ptr_in = &in[0];

    fgets(ptr_in, l_char, file);
    sscanf(ptr_in, "%lg %lg", &V, &I);
    fclose(file);
}


#endif
