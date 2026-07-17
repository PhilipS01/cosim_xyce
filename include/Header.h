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
    // Field-voltage reconstruction within a window (current-driven mode). Every mode CARRIES the
    // window start = the previous window's end V_field (reused -> C0-continuous seam, no solve there):
    //   0 = pointwise (default, secant): interior V from the accumulated secant (N field solves).
    //   1 = linear: straight ramp carried-start -> window-end value (ONE field solve / window).
    //   2 = average: like 1 but start = 0.5*(carried + end) (ONE solve; colleague's blend).
    //   3 = central-diff pointwise: interior V via local central difference (N solves; no window-start lag).
    // The window-end value uses the accumulated window secant (I_end - I0)/dt_win.
    unsigned reconstruct_mode = 0;
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
    // Secant-denominator floor in the Bfield impedance Z = Rrom + Lrom/dt (dt = time - t_abs_start):
    //   1 = on (default): dt -> MAX(time - t_abs_start, t_floor); guards the 1/0 at the exact window
    //       start (time == t_abs_start), where dt would be 0 -> Lrom/0 singularity.
    //   0 = off: dt = (time - t_abs_start) bare; Z -> infinity at the window start. Test only.
    unsigned use_t_floor = 1;
    // Window-seam handoff: how the carried terminal values (V0, I0 seeding the next window) are picked
    // from the converged iterate (the two solvers agree to within WR_tolerance there):
    //   0 = one-sided (default): V0 = circuit V(p), I0 = field I. (Restart checkpoint = Xyce state.)
    //   1 = midpoint: V0 = 0.5*(V_circuit + V_field), I0 = 0.5*(I_circuit + I_field). Seam-blend test
    //       (only the carried seeds are averaged; the Xyce restart checkpoint is unchanged).
    unsigned seam_average = 0;
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

void FEM_solver_voltage_driven_waveform(double I_win_start, unsigned N_field_eval_intervals);

// Current-driven (Neumann) field solver: reads the interface current waveform I(t) from i_prev_k.pwl,
// computes the field voltage V_field(t) = R_FEM*I + L_FEM*dI/dt (+ saturation) on the field-eval grid,
// optionally blends linear+const (reconstruct_mode), and writes it to vf_prev_k.pwl for the circuit's
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
