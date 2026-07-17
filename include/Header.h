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
    // Reduced-order model (used inside the Xyce Biface interface condition)
    double L_ROM = 0.9 * 1.6e-7;
    double R_ROM = 0.9 * 5.1e-4;
    // "True" field/FEM parameters (intentionally different from the ROM)
    double L_FEM = 1.6e-7;
    double R_FEM = 5.1e-4;
    // Nonlinearity of the distributed (FEM) device. 0 = linear (constant L_FEM).
    //   1 = magnetic saturation: flux lambda(I) = L_FEM*I_sat*atan(I/I_sat),
    //       so the small-signal inductance L(I) = L_FEM/(1+(I/I_sat)^2) drops as
    //       the core saturates. I_sat is the saturation current scale (A).
    // The Xyce ROM (Biface) stays linear, so the nonlinearity is a genuine
    // field/ROM model mismatch that the WR iteration must resolve.
    unsigned nonlin_model = 0;
    double I_sat = 100.0;
    // Circuit source on the hot node `s`. Selects the device WriteCircuitNetlist emits:
    //   0 = sinusoidal VOLTAGE source  (Bemf s 0 V={amp*sin(2*pi*f*t)})  -- default / legacy
    //   1 = sinusoidal CURRENT source  (Bemf 0 s I={amp*sin(2*pi*f*t)})
    //   2 = step/ramp VOLTAGE source   (Vemf s 0 PULSE(v_init v_final delay rise ...))
    unsigned source_kind = 0;
    // Sinusoidal source (kinds 0,1): amplitude is Volts (kind 0) or Amps (kind 1).
    double frequency = 50.0;
    double amplitude = 1.0;
    // Step/ramp voltage source (kind 2): initial/final level, onset delay, and ramp (rise) time.
    double step_v_initial = 0.0;
    double step_v_final   = 1.0;
    double step_delay     = 0.0;
    double step_rise      = 1.0e-4;
    // Series R/L/C on the source->port path (only nonzero elements are emitted; all-zero => s==p).
    // Rs = Ls = 0 reproduces the one-way toy (port voltage == source voltage).
    double R_series = 0.0;
    double L_series = 0.0;
    double C_series = 0.0;

    // --- Switch topologies (increment 2). circuit_kind selects the whole circuit side:
    //   0 = simple source (source_kind + series R/L/C above) -- increment 1, default.
    //   1 = #4 two-way switch, sine U, cap C: drive [0,t1) -> freewheel [t1,inf). (open = pre-t0 state)
    //   2 = #5 two-way switch, DC U, cap C:   drive [0,t1) -> freewheel [t1,inf).
    //   3 = #6 two-way switch, AC V_AC vs R:    AC-drive [0,t1) -> R-damp [t1,inf).
    // The field/interface (Vmeas, Bfield) is unchanged; the switch side attaches at port p.
    unsigned circuit_kind = 0;
    // Switch realization: 0 = behavioral resistor R={IF(t..,Ron,Roff)}; 1 = native Xyce S + .MODEL SW.
    unsigned switch_backend = 0;
    // Fixed throw instant (absolute time): drive -> freewheel at t1.
    double switch_t1 = 6.0e-3;
    // Switch-circuit passives: C for #4/#5 (F), R for #6 (Ohm).
    double switch_C = 1.0e-3;
    double switch_R = 1.0e4;
    // Closed / open switch resistances (behavioral gate levels and native SW model RON/ROFF).
    double switch_Ron  = 1.0e-3;
    double switch_Roff = 1.0e9;
    // Switch transition (rise/fall) time. A finite ramp (not an instantaneous jump) is essential:
    // an abrupt throw disconnects an ideal branch carrying inductive field current -> voltage kick
    // -> Xyce dt-collapse at the switch instant. The gate ramps over switch_trise at each edge.
    double switch_trise = 1.0e-5;
    // Time stepping / run duration.
    //   time_mode = 0 (source periods): duration = N_periods / frequency, N_steps_field =
    //               N_field_steps_per_source_period * N_periods  (meaningful for sinusoidal sources).
    //   time_mode = 1 (absolute end time): duration = t_end, N_steps_field = N_field_windows
    //               (for step/switch/custom sources that have no "period").
    unsigned time_mode = 0;
    double   t_end = 2.0e-2;
    unsigned N_field_windows = 50;
    unsigned N_periods = 1;
    unsigned N_field_steps_per_source_period = 50;
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

void CIRCUIT_solver(const double dt_circuit, const unsigned N_dt_circuit_per_dt_field, const double time_field);

double INTERFACE_condition(	const bool reset, const double dt_circuit, const double dt_field,
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

// Generiert die vollstaendige Schaltungs-Netzliste (wr_circuit.cir) aus g_cfg: Quelle (source_kind)
// + serielle R/L/C-Kette Quelle->Port p + feste WR-Schnittstelle (Vmeas, Bfield) + Includes. Wird
// EINMAL vor der Fensterschleife (und im emit-Modus) aufgerufen; Topologie ist fensterinvariant.
// Die pro-Fenster variablen Groessen bleiben in sim_params.inc/restart.inc. Inline (kein .INCLUDE
// der Elemente), damit der UI-Netzlisten-Parser (folgt keinen .INCLUDEs) die Schaltung zeichnen kann.
void WriteCircuitNetlist(const string& filename);

void WriteSimParams(const string& filename, double t_start, double t_stop, double t_abs_start, double i0, double rrom, double lrom, double f_src, double amp_src, double dIdt0, double r_series, double l_series, const unsigned N_xyce_samples);

// Generiert restart.inc (von wr_circuit.cir inkludiert): die fenster-spezifische
// .OPTIONS RESTART und .tran Zeile. Fenster 1 (first_window=true): frischer UIC-Transient
// ab t=0, schreibt Checkpoints (JOB). Fenster k>1: Restart aus committed_file (FILE=),
// schreibt neue Checkpoints. Zeit ist absolut: tran-Stoppzeit = t_stop (absolut).
void WriteRestartDirectives(const string& filename, bool first_window,
                            double dt_window, const string& ckpt_out_prefix,
                            const string& committed_file);

// Loescht alte Checkpoint-Kandidaten <prefix>* vor der WR-Schleife (verhindert, dass ein
// veralteter Kandidat als neuester ausgewaehlt wird).
void ClearCheckpoints(const string& prefix);

// Sucht den neuesten (nach mtime) Checkpoint <prefix>* und kopiert ihn nach committed_file
// als Restart-Basis fuer das naechste Zeitfenster.
void CommitCheckpoint(const string& prefix, const string& committed_file);

inline void Write_Terminal_results(const char s[80], const double V, const double I)
{
	FILE* file = fopen(s, "w");
	fprintf(file, " %12.5e   % 12.5e", V, I);
	fflush(file);
	fclose(file);
}

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
