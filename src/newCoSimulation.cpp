#include "../include/Header.h"
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <sstream>

#ifdef _WIN32
// RunXyce spawns Xyce once per WR iteration (windows x WR steps), so the per-call
// process cost matters. system() would route every one of those through cmd.exe --
// two CreateProcess calls instead of one, plus a console window per spawn.
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#endif


// Global tunable configuration (defaults defined in SimConfig). Overwritten by LoadConfig.
SimConfig g_cfg;

// User-requested output probes (Xyce print tokens, e.g. "V(a)", "I(Rr1)"). See Header.h.
std::vector<std::string> g_probes;

// Load probe tokens from probes.txt (one per line; '#'/blank ignored). Missing file -> no probes.
void LoadProbes(const string& filename)
{
    g_probes.clear();
    ifstream in(filename);
    if (!in) return;
    string line;
    while (getline(in, line)) {
        const size_t hash = line.find('#');
        if (hash != string::npos) line = line.substr(0, hash);
        // trim
        size_t a = line.find_first_not_of(" \t\r\n");
        if (a == string::npos) continue;
        size_t b = line.find_last_not_of(" \t\r\n");
        string tok = line.substr(a, b - a + 1);
        if (!tok.empty()) g_probes.push_back(tok);
    }
}

// Read a Xyce .prn and append its probe columns to the probe output file. The Xyce output
// columns are: Index TIME V(P) V(NX) I(VMEAS) <probe1> ... <probeN>; we copy TIME + the N
// trailing probe columns, re-indexed. skip_first_point drops the first row (shared window edge).
void appendProbeColumns(const string& xyce_prn, FILE* out, size_t n_probes,
                        unsigned long& global_index, bool skip_first_point)
{
    if (!out || n_probes == 0) return;
    FILE* f = fopen(xyce_prn.c_str(), "r");
    if (!f) return;
    char line[4096];
    if (!fgets(line, sizeof(line), f)) { fclose(f); return; }  // header
    const size_t ncol = 5 + n_probes;      // index,time,V(p),V(nx),I(Vmeas),probes...
    std::vector<double> v(ncol);
    bool first = true;
    while (fgets(line, sizeof(line), f)) {
        if (strncmp(line, "End", 3) == 0) break;
        istringstream iss(line);
        size_t got = 0;
        for (; got < ncol; ++got) { if (!(iss >> v[got])) break; }
        if (got < ncol) continue;          // malformed / short row
        if (first && skip_first_point) { first = false; continue; }
        first = false;
        fprintf(out, "%-10lu %-17.8e", global_index++, v[1]);   // index, time
        for (size_t k = 0; k < n_probes; ++k) fprintf(out, " %-17.8e", v[5 + k]);
        fprintf(out, "\n");
    }
    fclose(f);
    fflush(out);
}

// Loads tunable parameters from a key=value config file into the global g_cfg. Accepts "key=value"
// or "key value" ('='/','/tab are treated as separators); '#' starts a comment. Unknown keys are
// logged and ignored; N_xyce_coupling_intervals is a back-compat alias for N_xyce_samples. Returns
// false (and keeps all defaults) if the file cannot be opened.
bool LoadConfig(const string& filename)
{
    ifstream in(filename);
    if (!in) {
        cout << "LoadConfig: '" << filename << "' not found, using defaults." << endl;
        return false;
    }

    string line;
    while (getline(in, line)) {
        // strip comment
        const size_t hash = line.find('#');
        if (hash != string::npos) line = line.substr(0, hash);

        // accept "key=value" or "key value"
        for (char& c : line) if (c == '=' || c == ',' || c == '\t') c = ' ';

        istringstream iss(line);
        string key;
        if (!(iss >> key)) continue; // blank line
        double val;
        if (!(iss >> val)) {
            cout << "LoadConfig: ignoring malformed line for key '" << key << "'" << endl;
            continue;
        }

        if      (key == "L_ROM")                            g_cfg.L_ROM = val;
        else if (key == "R_ROM")                            g_cfg.R_ROM = val;
        else if (key == "L_FEM")                            g_cfg.L_FEM = val;
        else if (key == "R_FEM")                            g_cfg.R_FEM = val;
        else if (key == "nonlin_model")                     g_cfg.nonlin_model = (unsigned)val;
        else if (key == "I_sat")                            g_cfg.I_sat = val;
        else if (key == "t_end")                            g_cfg.t_end = val;
        else if (key == "N_field_windows")                  g_cfg.N_field_windows = (unsigned)val;
        else if (key == "N_field_eval_intervals")           g_cfg.N_field_eval_intervals = (unsigned)val;
        else if (key == "N_xyce_samples")                   g_cfg.N_xyce_samples = (unsigned)val;
        else if (key == "N_xyce_coupling_intervals")        g_cfg.N_xyce_samples = (unsigned)val;  // back-compat alias -> sampling role
        else if (key == "t_floor_frac")                     g_cfg.t_floor_frac = val;
        else if (key == "WRmaxSteps")                       g_cfg.WRmaxSteps = (unsigned)val;
        else if (key == "WR_tolerance")                     g_cfg.WR_tolerance = val;
        else if (key == "wr_convergence_method")            g_cfg.wr_convergence_method = (unsigned)val;
        else if (key == "coupling_mode")                    g_cfg.coupling_mode = (unsigned)val;
        else if (key == "reconstruct_mode")                 g_cfg.reconstruct_mode = (unsigned)val;
        else if (key == "interface_form")                   g_cfg.interface_form = (unsigned)val;
        else if (key == "precondition")                     g_cfg.precondition = (unsigned)val;
        else if (key == "use_t_floor")                      g_cfg.use_t_floor = (unsigned)val;
        else if (key == "interface_consistency")            g_cfg.interface_consistency = (unsigned)val;
        else if (key == "seam_average")                     g_cfg.seam_average = (unsigned)val;
        else if (key == "validation_mode")                  g_cfg.validation_mode = (unsigned)val;
        else cout << "LoadConfig: unknown key '" << key << "' ignored." << endl;
    }
    return true;
}


// A scalar time series (strictly increasing time t[], value y[]); the basic PWL/interface waveform.
struct Waveform {
    vector<double> t;
    vector<double> y;

    // Append (time, value); throws if time is not strictly greater than the last (monotonic guard).
    void push(double time, double value) {
    if (!t.empty() && time <= t.back()) {
        cerr << "Waveform time error: new time = " << time
             << ", last time = " << t.back() << endl;
        throw runtime_error("Waveform times must be strictly increasing.");
    }
        t.push_back(time);
        y.push_back(value);
    }
};

// Append (time, value) to wf, but if `time` equals the last timestamp (within tolerance) overwrite
// the last sample instead of pushing a duplicate. Throws if time goes backwards. Tolerates the
// shared window-edge point Xyce emits twice across a restart boundary.
void pushOrReplaceDuplicateTime(Waveform& wf, double time, double value)
{
    const double absTol = 1e-15;
    const double relTol = 1e-12;

    if (wf.t.empty()) {
        wf.push(time, value);
        return;
    }

    const double lastTime = wf.t.back();
    const double tol = absTol + relTol * std::max(std::abs(time), std::abs(lastTime));

    if (time > lastTime + tol) {
        wf.push(time, value);
        return;
    }

    if (std::abs(time - lastTime) <= tol) {
        // Same time point: keep the newest value.
        wf.y.back() = value;
        wf.t.back() = time;
        return;
    }

    throw std::runtime_error("Xyce output times are decreasing.");
}

// Circuit-side time series: the three printed columns per point -- port voltage vp = V(p),
// interface node vnx = V(nx), interface current i = I(Vmeas).
struct CircuitWaveform {
    vector<double> t;
    vector<double> vp;
    vector<double> vnx;
    vector<double> i;

    // Append one (time, vp, vnx, i) row; throws if time is not strictly increasing.
    void push(double time, double vp_value, double vnx_value, double i_value) {
        if (!t.empty() && time <= t.back()) {
            throw runtime_error("CircuitWaveform times must be strictly increasing.");
        }

        t.push_back(time);
        vp.push_back(vp_value);
        vnx.push_back(vnx_value);
        i.push_back(i_value);
    }
};

// Like the Waveform overload, for CircuitWaveform: append (time, vp, vnx, i) or overwrite the last
// row when `time` coincides with the last timestamp (within tolerance). Throws if time decreases.
void pushOrReplaceDuplicateTime(
    CircuitWaveform& wf,
    double time,
    double vp,
    double vnx,
    double i)
{
    const double absTol = 1e-15;
    const double relTol = 1e-12;

    if (wf.t.empty()) {
        wf.push(time, vp, vnx, i);
        return;
    }

    const double lastTime = wf.t.back();
    const double tol = absTol + relTol * std::max(std::abs(time), std::abs(lastTime));

    if (time > lastTime + tol) {
        wf.push(time, vp, vnx, i);
        return;
    }

    if (std::abs(time - lastTime) <= tol) {
        wf.t.back() = time;
        wf.vp.back() = vp;
        wf.vnx.back() = vnx;
        wf.i.back() = i;
        return;
    }

    throw runtime_error("CircuitWaveform times are decreasing.");
}


// Writes a waveform as a Xyce PWL FILE table (one "time value" pair per line, full double precision).
// Throws if the waveform is empty or its t/y sizes disagree.
void writePWLFile(const string& filename, const Waveform& wf) {
    if (wf.t.size() != wf.y.size() || wf.t.empty()) {
        throw runtime_error("Invalid waveform.");
    }

    ofstream out(filename);
    if (!out) {
        throw runtime_error("Could not open PWL file: " + filename);
    }

    out << scientific << setprecision(16);

    for (size_t i = 0; i < wf.t.size(); ++i) {
        out << wf.t[i] << " " << wf.y[i] << "\n";
    }
}

// Writes a 2-point PWL as a linear ramp from (t_start, value_start) with slope `slope`.
// slope=0 (default) -> constant PWL. slope != 0 sets a consistent initial slope of the
// seed waveform V(iprev)/V(vfprev) at the window start (dIdt0/dVdt0 carried from the previous window).
void WriteInitialPwl(
    const char* filename,
    double t_start,
    double t_stop,
    double value_start,
    double slope = 0.0)
{
    FILE* file = fopen(filename, "w");
    if (!file) {
        throw std::runtime_error("Could not open PWL file.");
    }

    const double value_stop = value_start + slope * (t_stop - t_start);

    fprintf(file, " %.16e %.16e\n", t_start, value_start);
    fprintf(file, " %.16e %.16e\n", t_stop,  value_stop);

    fclose(file);
}

// Backward-difference companion for i_prev_k.pwl, consumed ONLY by the naive/inconsistent "real RL"
// interface arm (interface_consistency=1). On the SAME field grid as `current` (nodes t_r^F) it emits
// the field-grid (h_F) BDF-1 backward difference of the interface current
//     d_r = (i_F(t_r) - i_F(t_{r-1})) / (t_r - t_{r-1})     for r >= 1,
// so the netlist can read L*di_F/dt directly as V(vdidtf) instead of reconstructing an anchored secant
// (correct only for ONE field evaluation per window). The window-start node r=0 has no prior sample;
// we forward-fill the first interval's slope using the carried window-start anchor I0:
//     d_0 = (i_F(t_1) - I0) / (t_1 - t_0).
// In the voltage-driven field solver i_F(t_0) == I0 by construction, so there d_0 == d_1 (a flat first
// interval). Times are absolute, identical to i_prev_k.pwl, so the PWL FILE source aligns with it.
void WriteFieldDidtPwl(const char* filename, const Waveform& current, double I0)
{
    const size_t N = current.t.size();
    Waveform didt;
    if (N < 2) {
        // Degenerate single-sample window: no interval -> zero derivative.
        didt.push(current.t.empty() ? 0.0 : current.t[0], 0.0);
        writePWLFile(filename, didt);
        return;
    }
    // r = 0: forward-filled first-interval slope, anchored to the window-start current I0.
    didt.push(current.t[0], (current.y[1] - I0) / (current.t[1] - current.t[0]));
    // r >= 1: plain backward difference between consecutive field nodes.
    for (size_t r = 1; r < N; ++r) {
        didt.push(current.t[r], (current.y[r] - current.y[r - 1]) / (current.t[r] - current.t[r - 1]));
    }
    writePWLFile(filename, didt);
}

// Reads a "time value" PWL table back into a Waveform (via pushOrReplaceDuplicateTime, so a repeated
// final/edge timestamp collapses). Throws if the file cannot be opened or contains no data.
Waveform readPWLFile(const string& filename)
{
    ifstream in(filename);
    if (!in) {
        throw runtime_error("Could not open PWL file: " + filename);
    }

    Waveform wf;

    double t, y;
    while (in >> t >> y) {
        pushOrReplaceDuplicateTime(wf, t, y);
    }

    if (wf.t.empty()) {
        throw runtime_error("PWL file contains no data: " + filename);
    }

    return wf;
}


//The master process coordinates the simulation
//It can be realized where ever it's easiest, e.g., inside the circuit solver, or inside the FEM solver, or exterior
void MasterProcess()
{
    cout << "Master process started" << endl;

    // Generate the circuit netlist (topology is window-invariant -> once here).
    // Overwrites wr_circuit.cir with the custom node-graph (circuit_spec.txt) + the fixed
    // WR interface. Per-window variable parameters still come from sim_params.inc.
    WriteCircuitNetlist("wr_circuit.cir");

    // File for the field solution at the synchronization points (here equal to the field steps)
    FILE* file_Field = fopen("Field_solution.prn", "w");
    fprintf(file_Field, "Index       TIME              V(FIELD)          I(FIELD)\n");
    fflush(file_Field);

    // File for the reconstructed/extrapolated field solution, i.e. the field data Xyce read in
    FILE* file_Field_waveform = fopen("Field_waveform_solution.prn", "w");
    fprintf(file_Field_waveform, "Index       TIME              V(FIELD)          I(FIELD)\n");
    fflush(file_Field_waveform);
    
    FILE* file_Circuit = fopen("Circuit_solution.prn", "w");
    fprintf(file_Circuit, "Index       TIME              V(P)              V(NX)             I(VMEAS)\n");
    fflush(file_Circuit);

    FILE* file_WR_error = fopen("WR_error.txt", "w");
    fprintf(file_WR_error, "   Time, WR_TotalRelErr, N_iterations, Converged, relI_FC, relV_FC \n");
    fflush(file_WR_error);

    // Optional user-probe output: one column per probe token (see LoadProbes/g_probes).
    FILE* file_Probes = nullptr;
    if (!g_probes.empty()) {
        file_Probes = fopen("Probes_solution.prn", "w");
        fprintf(file_Probes, "Index       TIME");
        for (const string& pr : g_probes) fprintf(file_Probes, "              %s", pr.c_str());
        fprintf(file_Probes, "\n");
        fflush(file_Probes);
    }

    unsigned long global_field_index = 0;
    unsigned long global_field_waveform_index = 0;
    unsigned long global_circuit_index = 0;
    unsigned long global_probe_index = 0;

    // All tunable parameters come from g_cfg (loaded from sim_config.txt; see SimConfig).
    // The impedance values (L_ROM, R_ROM) are a reduced-order model of the field domain used
    // inside the Xyce Bfield interface condition; they need not match the "true" FEM values.
    const double L_ROM = g_cfg.L_ROM;
    const double R_ROM = g_cfg.R_ROM;

    // Time stepping (the circuit side is handled dynamically by Xyce). The run spans [0, t_end],
    // split into N_field_windows equal WR windows: dt_field = t_end / N_field_windows.
    const unsigned N_steps_field = (g_cfg.N_field_windows >= 1) ? g_cfg.N_field_windows : 1u;
    if (g_cfg.t_end <= 0.0) throw runtime_error("MasterProcess: t_end must be > 0.");
    const double dt_field = g_cfg.t_end / N_steps_field;

    // WR parameters
    const unsigned WRmaxSteps = g_cfg.WRmaxSteps;
    const double WR_tolerance = g_cfg.WR_tolerance;
    // Xyce-solution sampling, decoupled from Xyce's adaptive step. Controls TWO things with one value:
    // (1) Xyce print cadence dt_print = dt_field/N (WriteSimParams) and (2) resample resolution of the
    // interface PWL (V(p)->vf_prev_k.pwl / I(Vmeas)->i_prev_k.pwl, both directions). t_floor is now a
    // SEPARATE knob (t_floor_frac). See voltage_driven_refactor.tex.
    const unsigned N_xyce_samples = g_cfg.N_xyce_samples;
    // FEM evaluation intervals per window (FEM-PWL has N+1 points). =1 → 2-point ramp → single
    // linear extrapolation of the field current per window. Higher → multi-rate, piecewise linear.
    const unsigned N_field_eval_intervals = g_cfg.N_field_eval_intervals;

    // Initialisierung
    double V_field = 0;
    double I_field = 0;
    double I0 = 0.0;
    double V0 = 0.0;
    double dIdt_0 = 0.0;
    // dV/dt of the port voltage at the window start (seed slope for vf_prev_k.pwl). At the window
    // end it is carried from the converged V(p) waveform (topology-agnostic). Window 1 uses a neutral
    // seed 0: the circuit is authored from circuit_spec.txt, so no source slope is known here; a wrong
    // slope seed only overshoots and costs WR iterations (the converged result is unaffected).
    double dVdt_0 = 0.0;

    //fprintf(file_Field, "%-10lu %-17.8e %-17.8e %-17.8e\n", global_field_index++, 0.0, V_field, I_field);

    // Outer loop over the time windows
    // here the field interval = WR time window
    for (unsigned step_field = 1; step_field <= N_steps_field; step_field++) {
        // nonlinear behavior of the RL element (dummy)
        //R_ROM *= 1.2;
        //L_ROM *= 1.05;

        const double t_start = (step_field - 1) * dt_field; // window-start absolute time (needed cpp-side); Xyce always simulates from 0 to the stop time when called.
        const double t_stop = step_field * dt_field;

        WriteSimParams("sim_params.inc", 0.0, dt_field, t_start, I0, R_ROM, L_ROM, N_xyce_samples);

        // Generate the window-specific restart/.tran directives (included by wr_circuit.cir).
        // Window 1: fresh UIC transient from 0; window k>1: restart from "restart_state".
        WriteRestartDirectives("restart.inc", step_field == 1, dt_field, "ckpt_out", "restart_state");

        // initialize with the last accepted values. Timestamps are ABSOLUTE ([t_start, t_stop]),
        // because Xyce evaluates the PWL FILE sources at the absolute simulation time (restart starts
        // at t_start, not at 0).
        WriteInitialPwl("vf_prev_k.pwl", t_start, t_stop, V0, dVdt_0);
        // i_prev_k.pwl: linear warm-start ramp (start I0, slope dIdt_0 carried from the previous window end).
        WriteInitialPwl("i_prev_k.pwl",   t_start, t_stop, I0, dIdt_0);
        // Companion for the naive "real RL" arm's V(vdidtf): the initial i_prev ramp has constant slope
        // dIdt_0, so its field-grid backward difference is the constant dIdt_0. Written every window so
        // didt_field_k.pwl already exists at WR iteration 1 (the circuit runs BEFORE the field solver).
        // Only the real-RL arm reads it; harmless (tiny) for the other arms.
        WriteInitialPwl("didt_field_k.pwl", t_start, t_stop, dIdt_0, 0.0);

        // Remove stale checkpoint candidates of this prefix so that CommitCheckpoint after the WR
        // loop is guaranteed to pick the freshly created candidate of this window.
        ClearCheckpoints("ckpt_out");

        CircuitWaveform circuit_sol; // circuit solution array (only needed for output/visualization)

        unsigned WR_iteration;
        double WR_rel_Error = 1.0;
        bool WR_converged = false;
        Waveform i_prev_last_iter; // i_m^(k-1) for the L1 convergence criterion; empty at the window start
        double V_field_last_iter = 0.0; // for the terminal-scalar criterion (reference CoSimulation_WR.cpp)
        double I_field_last_iter = 0.0;
        double V_circuit = 0.0, I_circuit = 0.0; // terminal circuit (V(nx), I(Vmeas)); FC transmission defect
        // WR iteration loop
        for (WR_iteration = 1; WR_iteration <= WRmaxSteps; WR_iteration++) {
            circuit_sol = CircuitWaveform{};

            // Call the circuit solver: Bfield returns I(Vmeas) = INTERFACE_condition(V(p), V(vfprev), V(iprev)).
            // ReadXyceResults writes the V(p) waveform to vf_prev_k.pwl for the FEM input.
            RunXyce("wr_circuit.cir");
            ReadXyceResults("wr_circuit.cir.prn", circuit_sol, t_start, t_stop, N_xyce_samples);

            // Call the FEM solver (dummy). Coupling direction per coupling_mode:
            //   0 voltage-driven: reads vf_prev_k.pwl (V(p)), writes i_prev_k.pwl (I_field).
            //   1 current-driven: reads i_prev_k.pwl (I(Vmeas)), writes vf_prev_k.pwl (V_field).
            // writes the end values to Field.txt
            if (g_cfg.coupling_mode == 1)
                FEM_solver_current_driven_waveform(I0, V0, N_field_eval_intervals);
            else
                FEM_solver_voltage_driven_waveform(I0, N_field_eval_intervals);
            Read_Terminal_results("Field.txt", V_field, I_field);
            // Terminal circuit values (V(nx), I(Vmeas)) written by ReadXyceResults above. Read every
            // iteration so the field-vs-circuit transmission defect can be reported regardless of the
            // gating metric (method 0 never reads them otherwise).
            Read_Terminal_results("Circuit.txt", V_circuit, I_circuit);

            // i_prev = field-current waveform of this iteration (FEM output)
            Waveform i_prev = readPWLFile("i_prev_k.pwl");
            std::cout << "PWL points: " << i_prev.t.size() << std::endl;

            // Check the convergence criterion. Method selectable via g_cfg.wr_convergence_method:
            //   0 = waveform L1 of the field current (this codebase)
            //   1 = terminal-scalar (reference CoSimulation_WR.cpp): field-vs-circuit +
            //       iteration-to-iteration change of the terminal values V,I.
            bool can_converge;
            if (g_cfg.wr_convergence_method == 1) {
                WR_rel_Error = eval_WR_convergence_terminal(
                    V_field, I_field, V_circuit, I_circuit,
                    V_field_last_iter, I_field_last_iter, WR_iteration);
                // Require >= 2 iterations (like method 0). The reference allows convergence from
                // iteration 1 on the FC terms alone, but near a window-start zero-crossing both solvers
                // start ~0 so relI_FC+relV_FC < tol trivially -> a spurious "converged" after ONE solve
                // with no actual coupling iteration. Gating on WR_iteration > 1 forces at least one real
                // exchange before the window can be accepted.
                can_converge = (WR_iteration > 1);
            } else {
                // L1 norm (hybrid relative/absolute) of (i_m^(k) - i_m^(k-1)); needs >=2 iterations
                // (iteration 1 returns the sentinel 1.0, no previous waveform).
                WR_rel_Error = eval_WR_convergence(i_prev, i_prev_last_iter, WR_iteration);
                can_converge = (WR_iteration > 1);
            }

            if (can_converge && WR_rel_Error < WR_tolerance) {
                WR_converged = true;
                break;
            }

            // keep the values of this iteration for the next convergence check
            i_prev_last_iter = i_prev;
            V_field_last_iter = V_field;
            I_field_last_iter = I_field;
        }

        // Screen output
        //
        // Write results to file

        PRINT(t_stop);
        PRINT(WR_rel_Error);

        // Cross-solver TRANSMISSION DEFECT at the window terminal, reported independently of the gating
        // metric (method 0 is blind to it; the inconsistent-interface fixpoint keeps it non-zero). Same
        // hybrid-relative convention as eval_WR_convergence_terminal: normalise by the field value unless
        // it is tiny (< 0.1). voltage-driven -> relI_FC is the shifted quantity; current-driven -> relV_FC.
        const double dI_FC   = fabs(I_field - I_circuit);
        const double relI_FC = (fabs(I_field) > 0.1) ? dI_FC / fabs(I_field) : dI_FC;
        const double dV_FC   = fabs(V_field - V_circuit);
        const double relV_FC = (fabs(V_field) > 0.1) ? dV_FC / fabs(V_field) : dV_FC;

        fprintf(file_WR_error, "%12.3e ,  %12.3e ,  %d ,  %d ,  %12.3e ,  %12.3e \n",
                t_stop, WR_rel_Error, WR_iteration, int(WR_converged), relI_FC, relV_FC);
        fflush(file_WR_error);


        if (!WR_converged) {
            throw runtime_error("WR did not converge on current time window. Stop simulation.");
        } else {
            PRINT("WR converged!");
            PRINT(V_field);
            PRINT(I_field);

            // Store the field solution at the endpoints (here equal to the synchronization points) in Xyce format
            fprintf(file_Field, "%-10lu %-17.8e %-17.8e %-17.8e\n", global_field_index++, t_stop, V_field, I_field);
            fflush(file_Field);

            // Save the converged end state as the restart basis for the next window.
            // (All WR iterations wrote the same candidate ckpt_out<t_stop>; the
            // newest one is the converged one.)
            CommitCheckpoint("ckpt_out", "restart_state");

            // Append the Xyce solution of this time window to the previous ones.
            // Waveforms now carry ABSOLUTE time -> offset 0.0 (no further shifting).
            const bool skip_first_point = (step_field > 1);
            appendCircuitWaveformXyceStyle(file_Circuit, circuit_sol, 0.0, global_circuit_index, skip_first_point);

            // user probes: copy the extra .print columns of the converged Xyce solve for this window
            if (file_Probes)
                appendProbeColumns("wr_circuit.cir.prn", file_Probes, g_probes.size(),
                                   global_probe_index, skip_first_point);

            // We also store the field WAVEFORMS (i.e. not only endpoints) in Xyce format;
            // for that we read in the converged waveforms (last iteration)
            Waveform vf_conv = readPWLFile("vf_prev_k.pwl"); // port voltage V(p) from Xyce (coupling grid)
            Waveform i_conv  = readPWLFile("i_prev_k.pwl");  // field current I_field from FEM (FEM grid)

            // Re-sample both onto a common FEM-eval grid (N_field_eval_intervals+1),
            // so that appendFieldWaveformXyceStyle sees matching timestamps.
            // (vf_conv has N_xyce_samples+1 nodes, i_conv has N_field_eval_intervals+1)
            Waveform vf_endpoints = resampleWaveformUniform(
                vf_conv, t_start, t_stop, N_field_eval_intervals
            );
            Waveform i_endpoints = resampleWaveformUniform(
                i_conv, t_start, t_stop, N_field_eval_intervals
            );

            // and append them to the previous time windows (waveforms already absolute -> offset 0.0)
            appendFieldWaveformXyceStyle(file_Field_waveform, vf_endpoints, i_endpoints, 0.0, global_field_waveform_index, skip_first_point);

            // Record the end values as the initial values of the next time window.
            // Default (one-sided): V0 = circuit V(p), I0 = field I. seam_average: midpoint of both
            // terminals (they agree < WR_tolerance, so the blend stays within tol of either side).
            const double V_circ_end = vf_conv.y.back();  // circuit port voltage at window end
            if (g_cfg.seam_average) {
                double V_circ_dummy, I_circ_end;
                Read_Terminal_results("Circuit.txt", V_circ_dummy, I_circ_end);  // I(Vmeas) at end
                V0 = 0.5 * (V_circ_end + V_field);
                I0 = 0.5 * (I_circ_end + I_field);
            } else {
                V0 = V_circ_end;  // port voltage at window end (= V(p) at t=dt_field)
                I0 = I_field;     // field current at window end (FEM output)
            }

            // dI/dt at window end = dI/dt at the start of the next window (continuity)
            const size_t n_i = i_conv.t.size();
            if (n_i >= 2) {
                const double dt_end = i_conv.t[n_i - 1] - i_conv.t[n_i - 2];
                dIdt_0 = (i_conv.y[n_i - 1] - i_conv.y[n_i - 2]) / dt_end;
            } else {
                dIdt_0 = 0.0;
            }

            // dV/dt of the port voltage at window end -> seed slope of the next window.
            // Measured from the converged V(p) waveform (topology-agnostic): correct also
            // when V(p) != Vsrc (series Rs/Ls), where the pure source slope would be wrong.
            const size_t n_v = vf_conv.t.size();
            if (n_v >= 2) {
                const double dt_end_v = vf_conv.t[n_v - 1] - vf_conv.t[n_v - 2];
                dVdt_0 = (vf_conv.y[n_v - 1] - vf_conv.y[n_v - 2]) / dt_end_v;
            } else {
                dVdt_0 = 0.0;
            }
        }

    }

    fprintf(file_Circuit, "End\n");
    fprintf(file_Field,   "End\n");
    fclose(file_Circuit);
    fclose(file_Field);
    fclose(file_WR_error);
    fclose(file_Field_waveform);
    if (file_Probes) { fprintf(file_Probes, "End\n"); fclose(file_Probes); }

}

// Validation-mode driver (validation_mode=1). The true field is stamped as real Xyce devices
// (R_FEM + L_FEM series on the port branch; see the validation branch of WriteCircuitNetlist) and
// the whole circuit is solved as ONE Xyce transient over [0, t_end] -- no WR loop, no dummy field
// solver, no PWL exchange, no windowing. Writes the SAME output files as MasterProcess so the UI
// plots / CSV export / summary work unchanged: V(FIELD)=V(nx), I(FIELD)=I(Vmeas) (they equal the
// circuit port quantities, so the reference "circuit" and "field" curves coincide). WR_error.txt
// gets only its header -> the convergence plot renders empty and the WR summary cards are omitted
// (both degrade gracefully). Linear field only (nonlin_model warned + ignored in the netlist).
void MonolithicValidationSolve()
{
    cout << "Validation mode: monolithic reference solve (real R_FEM + L_FEM devices)." << endl;
    if (g_cfg.nonlin_model != 0)
        cout << "  WARNING: nonlin_model=" << g_cfg.nonlin_model
             << " (saturation) is IGNORED in validation mode -- L_FEM is a linear inductor." << endl;

    if (g_cfg.t_end <= 0.0) throw runtime_error("MonolithicValidationSolve: t_end must be > 0.");

    // Netlist with the real FEM devices (validation branch inside WriteCircuitNetlist).
    WriteCircuitNetlist("wr_circuit.cir");

    // Reference print grid = the WR run's TOTAL resolution: N_field_windows * N_xyce_samples
    // samples over the whole [0, t_end] (dt_print = t_end / N_eff), so the reference and the coupled
    // run line up point-for-point when compared (plots / CSV export).
    const unsigned N_win  = (g_cfg.N_field_windows >= 1) ? g_cfg.N_field_windows : 1u;
    const unsigned N_samp = (g_cfg.N_xyce_samples  >= 1) ? g_cfg.N_xyce_samples  : 1u;
    const unsigned N_eff  = N_win * N_samp;

    // sim_params.inc provides dt_print (= t_end/N_eff) and t_stop that the .tran line references.
    // Rrom/Lrom are written but unused by the validation netlist.
    WriteSimParams("sim_params.inc", 0.0, g_cfg.t_end, 0.0, 0.0, g_cfg.R_ROM, g_cfg.L_ROM, N_eff);

    // One fresh UIC transient over the whole [0, t_end] -- no restart/checkpoint (single solve).
    // The step CEILING (4th .tran field) = dt_print bounds Xyce's adaptive stepper so this single
    // long transient is resolved on the same grid as the stitched WR run; without it the integrator
    // coasts over quiet stretches and the monolithic reference comes out far coarser than the
    // per-window WR solves.
    {
        ofstream rst("restart.inc");
        if (!rst) throw runtime_error("MonolithicValidationSolve: could not open restart.inc");
        rst << ".tran {dt_print} {t_stop} {t_abs_start} {dt_print} UIC\n";
    }

    RunXyce("wr_circuit.cir");

    // Parse the solve onto the uniform grid (reuses the WR parser; its Circuit.txt / *.pwl
    // side-writes are harmless leftovers in this mode).
    CircuitWaveform circuit_sol;
    ReadXyceResults("wr_circuit.cir.prn", circuit_sol, 0.0, g_cfg.t_end, N_eff);
    if (circuit_sol.t.empty())
        throw runtime_error("MonolithicValidationSolve: no circuit data parsed from wr_circuit.cir.prn.");

    // Field terminal = the real-device branch: V(FIELD) = V(nx), I(FIELD) = I(Vmeas).
    Waveform vfield, ifield;
    for (size_t n = 0; n < circuit_sol.t.size(); ++n) {
        vfield.push(circuit_sol.t[n], circuit_sol.vnx[n]);
        ifield.push(circuit_sol.t[n], circuit_sol.i[n]);
    }

    // Circuit_solution.prn: Index TIME V(P) V(NX) I(VMEAS)  (+ End).
    FILE* file_Circuit = fopen("Circuit_solution.prn", "w");
    fprintf(file_Circuit, "Index       TIME              V(P)              V(NX)             I(VMEAS)\n");
    unsigned long idx_c = 0;
    appendCircuitWaveformXyceStyle(file_Circuit, circuit_sol, 0.0, idx_c, /*skip_first_point=*/false);
    fprintf(file_Circuit, "End\n");
    fclose(file_Circuit);

    // Field_waveform_solution.prn: dense field waveform (no trailing "End", matching MasterProcess).
    FILE* file_Field_waveform = fopen("Field_waveform_solution.prn", "w");
    fprintf(file_Field_waveform, "Index       TIME              V(FIELD)          I(FIELD)\n");
    unsigned long idx_fw = 0;
    appendFieldWaveformXyceStyle(file_Field_waveform, vfield, ifield, 0.0, idx_fw, /*skip_first_point=*/false);
    fclose(file_Field_waveform);

    // Field_solution.prn: the window-end (final) row -> summary card final V/I.
    FILE* file_Field = fopen("Field_solution.prn", "w");
    fprintf(file_Field, "Index       TIME              V(FIELD)          I(FIELD)\n");
    const size_t last = circuit_sol.t.size() - 1;
    fprintf(file_Field, "%-10lu %-17.8e %-17.8e %-17.8e\n",
            0UL, circuit_sol.t[last], circuit_sol.vnx[last], circuit_sol.i[last]);
    fprintf(file_Field, "End\n");
    fclose(file_Field);

    // WR_error.txt: header only (no WR iterations) -> empty convergence plot, no WR summary cards.
    FILE* file_WR_error = fopen("WR_error.txt", "w");
    fprintf(file_WR_error, "   Time, WR_TotalRelErr, N_iterations, Converged, relI_FC, relV_FC \n");
    fclose(file_WR_error);

    // Probes_solution.prn: user probes are in the .print line, so copy their columns from the raw
    // monolithic solve (dense, full resolution), same format as the WR path (MasterProcess). Without
    // this the probe plot would show stale data from an earlier run.
    if (!g_probes.empty()) {
        FILE* file_Probes = fopen("Probes_solution.prn", "w");
        fprintf(file_Probes, "Index       TIME");
        for (const string& pr : g_probes) fprintf(file_Probes, "              %s", pr.c_str());
        fprintf(file_Probes, "\n");
        unsigned long idx_p = 0;
        appendProbeColumns("wr_circuit.cir.prn", file_Probes, g_probes.size(), idx_p,
                           /*skip_first_point=*/false);
        fprintf(file_Probes, "End\n");
        fclose(file_Probes);
    } else {
        // No probes this run: remove any stale Probes_solution.prn so the UI doesn't plot old data.
        remove("Probes_solution.prn");
    }

    cout << "Validation reference written "
         << "(Circuit_solution.prn / Field_solution.prn / Field_waveform_solution.prn"
         << (g_probes.empty() ? "" : " / Probes_solution.prn") << ")." << endl;
}

// FEM solver (voltage-driven): takes the port-voltage waveform V_p(t) from vf_prev_k.pwl and
// computes the field current I_field(t) by integrating the field ODE V_p = R_FEM*I + L_FEM*dI/dt
// with backward Euler (BDF-1) on the field grid: at each node I_ref = the PREVIOUS field node and
// dt = the LOCAL step h_F (reconstruct_mode 0). This is a self-contained dummy field solver -- it
// does NOT mirror the netlist Bfield's window-anchored secant (an earlier version did; that was an
// interface-matching artifact, not physics). At N_field_eval_intervals = 1 the single BE step
// reduces to the window secant. Multi-rate: the FEM works on a coarser grid than Xyce internally.
void FEM_solver_voltage_driven_waveform(double I_win_start, unsigned N_field_eval_intervals)
{
    // Reads the port voltage waveform stored by ReadXyceResults from V(p).
    const Waveform vport = readPWLFile("vf_prev_k.pwl");

    // Time is now absolute (vport.t.front() = absolute window start, not 0). The FEM solver is
    // offset-agnostic: every backward-Euler step below uses a LOCAL difference (t[j] - t[j-1]), which
    // is independent of the absolute start time. Hence no "starts at 0" check.

    const size_t Nv = vport.t.size();
    if (Nv < 2) {
        throw runtime_error("FEM_solver: vf_prev_k.pwl needs at least 2 points.");
    }
    if (N_field_eval_intervals < 1) {
        throw runtime_error("FEM_solver: N_field_eval_intervals must be >= 1.");
    }

    // FEM/reference parameters (from config), intentionally may differ from ROM values.
    const double L_FEM = g_cfg.L_FEM;
    const double R_FEM = g_cfg.R_FEM;

    // Distributed-device nonlinearity. Linear unless nonlin_model selects a model
    // (and I_sat is sane). For saturation the flux lambda(I) is nonlinear so L = dlambda/dI
    // drops with current; the per-interval update is no longer closed-form and needs Newton.
    const bool   saturating = (g_cfg.nonlin_model == 1) && (g_cfg.I_sat > 0.0);
    const double I_sat = g_cfg.I_sat;
    // lambda(I) = L_FEM*I_sat*atan(I/I_sat)  → lambda'(I) = L(I) = L_FEM/(1+(I/I_sat)^2).
    // As I/I_sat → 0 this reduces to lambda ≈ L_FEM*I, L ≈ L_FEM (linear limit).
    auto flux  = [&](double I) { return L_FEM * I_sat * std::atan(I / I_sat); };
    auto L_dyn = [&](double I) { const double r = I / I_sat; return L_FEM / (1.0 + r * r); };

    const double t_start = vport.t.front();
    const double t_end   = vport.t.back();

    // Multi-rate: the FEM re-samples the port voltage onto a coarser FEM grid.
    const Waveform V_eval = resampleWaveformUniform(
        vport, t_start, t_end, N_field_eval_intervals
    );
    const size_t N = V_eval.t.size(); // = N_field_eval_intervals + 1

    // Backward-Euler (BDF-1) time stepping of the field ODE on the field grid. Per step the discrete
    // constitutive relation is
    //     V_p(t_j) = R_FEM*I_j + (lambda(I_j) - lambda(I_{j-1})) / (t_j - t_{j-1})
    // i.e. I_ref = the PREVIOUS field node current and dt = the LOCAL step -> a consistent dI/dt. An
    // earlier version anchored both to the window start (I_ref = I_win_start, dt = t - t_win_start) so
    // the dummy would reproduce the netlist Bfield's accumulated secant; that is an interface-matching
    // artifact, not a field-ODE integration, and is NOT correct for >1 field eval per window.
    // The WR fixpoint stays clean: the default Bfield's correction Z*(I(Vmeas) - V(iprev)) vanishes at
    // convergence regardless of Z's denominator, so switching the FEM to local BDF-1 changes the field
    // response only for N_field_eval_intervals >= 2 (more accurate), not the fixpoint's existence. At
    // N_field_eval_intervals = 1 the single BE step == the window secant (bit-for-bit unchanged).
    const double t_win_start = V_eval.t.front();
    const double t_win_end   = V_eval.t.back();

    // Solve V_p = R_FEM*I + (lambda(I) - lambda(I_ref))/dt for the field current I (linear closed form,
    // then Newton if the core saturates). I_ref/dt are the LOCAL backward-Euler pair (previous node,
    // local step); reconstruct_mode 1 below reuses the same solve as a single BE step over the window.
    auto solve_I = [&](double V_p, double I_ref, double dt) -> double {
        if (dt <= 0.0) return I_ref;                       // singular window start: initial condition
        double I = (V_p + L_FEM * I_ref / dt) / (R_FEM + L_FEM / dt);
        if (saturating) {
            const double lam0 = flux(I_ref);
            for (int it = 0; it < 50; ++it) {
                const double g  = R_FEM * I + (flux(I) - lam0) / dt - V_p;
                const double gp = R_FEM + L_dyn(I) / dt;
                const double dI = g / gp;
                I -= dI;
                if (std::fabs(dI) <= 1e-12 + 1e-10 * std::fabs(I)) break;
            }
        }
        return I;
    };

    // Field-current reconstruction within the window (reconstruct_mode; symmetric to the current-driven
    // field-voltage reconstruction). Both modes CARRY the window start = I_win_start (previous window's
    // end current) -> C0-continuous seam, no solve there.
    //   0 pointwise (BDF-1, DEFAULT): local backward-Euler step per field node (I_ref = previous node,
    //     dt = local h_F). N solves/window; the physically correct dummy-solver field integration.
    //   1 linear : straight ramp from the carried start to the window-end current (ONE solve/window).
    Waveform current;
    if (g_cfg.reconstruct_mode == 1) {
        const double I_end = solve_I(V_eval.y[N - 1], I_win_start, t_win_end - t_win_start);
        for (size_t j = 0; j < N; ++j) {
            const double frac = (t_win_end > t_win_start)
                              ? (V_eval.t[j] - t_win_start) / (t_win_end - t_win_start) : 0.0;
            current.push(V_eval.t[j], I_win_start + frac * (I_end - I_win_start));
        }
    } else {
        current.push(t_win_start, I_win_start);            // window start = carried initial condition
        for (size_t j = 1; j < N; ++j) {
            const double dt = V_eval.t[j] - V_eval.t[j - 1];                        // LOCAL field step h_F
            current.push(V_eval.t[j], solve_I(V_eval.y[j], current.y.back(), dt));  // I_ref = previous node
        }
    }

    // This file is used by Xyce as V(iprev) in the next WR iteration.
    writePWLFile("i_prev_k.pwl", current);

    // Companion for the naive "real RL" arm: the field-grid (h_F) BDF-1 backward difference of I_field,
    // read by Xyce as V(vdidtf). Formed HERE on the field grid (over h_F), not reconstructed via DDT --
    // this is what lets that arm run with >1 field evaluation per window. Anchor = I_win_start (== .PARAM I0).
    WriteFieldDidtPwl("didt_field_k.pwl", current, I_win_start);

    // Terminal values at the window end for convergence propagation into the next window.
    Write_Terminal_results("Field.txt", vport.y.back(), current.y.back());
}

// Current-driven (Neumann) field solver. Reads the interface current I(t) (i_prev_k.pwl, circuit
// output) and returns the field voltage V(t) = R_FEM*I + dlambda/dt in vf_prev_k.pwl -- this is the
// V(vfprev) base of the circuit's matched-secant Bfield. The secant is a Newton linearisation of the
// field V-I characteristic in BOTH coupling directions; current-driven's difference is only that this
// base voltage is the field's OWN directly-computed V (an accurate base), not the circuit's prior port
// voltage -- so the amplified window-start term stays small. reconstruct_mode:
//   0 pointwise (default): V computed at every field-eval point by a LOCAL backward-Euler step of the
//      flux (dlambda/dt over consecutive interface samples) -- symmetric to the voltage-driven solver.
//      Uses V_field_last_time as the carried window-start value (I_win_start only feeds mode 1's V_end).
//   1 linear: replace the interior with a straight ramp from V_field_last_time (previous window end)
//      to this window's end value (the 2-point reconstruction; uses V0).
void FEM_solver_current_driven_waveform(double I_win_start, double V_field_last_time,
                                        unsigned N_field_eval_intervals)
{
    const Waveform iface = readPWLFile("i_prev_k.pwl");   // interface current I(t) on the coupling grid
    if (iface.t.size() < 2)
        throw runtime_error("FEM_solver_current_driven: i_prev_k.pwl needs at least 2 points.");
    if (N_field_eval_intervals < 1)
        throw runtime_error("FEM_solver_current_driven: N_field_eval_intervals must be >= 1.");

    const double L_FEM = g_cfg.L_FEM;
    const double R_FEM = g_cfg.R_FEM;
    const bool   saturating = (g_cfg.nonlin_model == 1) && (g_cfg.I_sat > 0.0);
    const double I_sat = g_cfg.I_sat;
    auto flux = [&](double I) { return L_FEM * I_sat * std::atan(I / I_sat); };  // lambda(I)
    auto dlam = [&](double Ib, double Ia) {                                      // lambda(Ib)-lambda(Ia)
        return saturating ? (flux(Ib) - flux(Ia)) : (L_FEM * (Ib - Ia));
    };

    // Multi-rate: resample the interface current onto the field-eval grid (N_field_eval+1 points).
    const Waveform I_eval = resampleWaveformUniform(
        iface, iface.t.front(), iface.t.back(), N_field_eval_intervals);
    const size_t   N = I_eval.t.size();               // = N_field_eval_intervals + 1
    const double   t_win_start = I_eval.t.front();
    const double   t_end       = I_eval.t.back();

    // The window START field voltage is the carried previous-window END value (V_field_last_time):
    // reused -> the seam is C0-continuous by construction AND no solve is spent there. reconstruct_mode:
    //   1 linear: ONE new evaluation (the window end), straight-line interior => 1 solve/window.
    //   0 BDF-1: pointwise interior via LOCAL backward Euler dlambda/dt (previous node, local step),
    //     N solves, start still carried for seam continuity.
    const double I_end  = I_eval.y.back();
    const double dt_win = t_end - t_win_start;
    // Window-end field voltage via a SINGLE backward-Euler step over the whole window (used by
    // reconstruct_mode 1 only): flux difference lambda(I_end) - lambda(I_win_start) over dt_win.
    const double V_end = (dt_win > 0.0)
                         ? R_FEM * I_end + dlam(I_end, I_win_start) / dt_win
                         : R_FEM * I_end;

    // Window start = carried previous-window end value (all modes; seam-continuous, no solve).
    const double V_start_win = V_field_last_time;

    Waveform vfield;
    if (g_cfg.reconstruct_mode == 1) {
        // 1 solve/window: carried start, straight line to the end.
        for (size_t j = 0; j < N; ++j) {
            const double frac = (dt_win > 0.0) ? (I_eval.t[j] - t_win_start) / dt_win : 0.0;
            vfield.push(I_eval.t[j], V_start_win + frac * (V_end - V_start_win));
        }
    } else {
        // pointwise (0 BDF-1): start carried (reuse the seam); interior/end via LOCAL backward Euler.
        vfield.push(t_win_start, V_start_win);
        for (size_t j = 1; j < N; ++j) {
            const double dl_dt = dlam(I_eval.y[j], I_eval.y[j - 1]) / (I_eval.t[j] - I_eval.t[j - 1]);
            vfield.push(I_eval.t[j], R_FEM * I_eval.y[j] + dl_dt);
        }
    }
    const double V_field_end = vfield.y.back();

    writePWLFile("vf_prev_k.pwl", vfield);                     // field -> circuit: the field voltage
    Write_Terminal_results("Field.txt", V_field_end, I_end);  // (V_field_end, I_end)
}

// WR convergence criterion based on the current waveforms of adjacent iterations:
//
//   ∫|i_m^(k)(t) - i_m^(k-1)(t)| dt   /   ∫|i_m^(k)(t)| dt   ≤   WR_tolerance
//
// Integration via the trapezoidal rule over the coupling grid (i_curr and i_prev_iter share the same
// uniform time grid from resampleWaveformUniform). The norm is HYBRID: relative (as written above) when
// the average |i| over the window exceeds 0.1, absolute (numerator over the window length) when the
// current is ~0 and the relative form is ill-defined. i_curr is the i_prev_k.pwl waveform of this
// iteration -- the FEM field current in voltage-driven coupling, the circuit current I(Vmeas) in
// current-driven; by the transmission condition i_m^(k) = i_c^(k) the two coincide at convergence.
double eval_WR_convergence(const Waveform& i_curr, const Waveform& i_prev_iter, const unsigned WR_iteration)
{
    if (WR_iteration < 2 || i_prev_iter.t.empty()) {
        cout << "WR-step  " << WR_iteration << ": skip (no previous waveform)" << endl;
        return 1.0;
    }

    if (i_curr.t.size() != i_prev_iter.t.size() || i_curr.t.size() < 2) {
        throw runtime_error("eval_WR_convergence: waveform size mismatch or too small.");
    }

    double num = 0.0; // ∫ |i^(k) - i^(k-1)| dt
    double den = 0.0; // ∫ |i^(k)|         dt
    for (size_t n = 0; n + 1 < i_curr.t.size(); ++n) {
        const double dt = i_curr.t[n + 1] - i_curr.t[n];

        const double diff_a = std::fabs(i_curr.y[n]     - i_prev_iter.y[n]);
        const double diff_b = std::fabs(i_curr.y[n + 1] - i_prev_iter.y[n + 1]);
        num += 0.5 * (diff_a + diff_b) * dt;

        const double abs_a = std::fabs(i_curr.y[n]);
        const double abs_b = std::fabs(i_curr.y[n + 1]);
        den += 0.5 * (abs_a + abs_b) * dt;
    }

    // Hybrid abs/rel norm (mirrors eval_WR_convergence_terminal). num = INT|di| dt, den = INT|i| dt;
    // dividing both by the window duration T gives average currents. When the average |i| exceeds 0.1
    // use the relative norm num/den (T cancels -> identical to the plain relative L1); when the current
    // is (near-)zero over the whole window den->0 and the relative norm is ill-defined, so fall back to
    // the ABSOLUTE average change num/T. This replaces the old hard throw on den==0, which aborted the
    // whole run on any zero-current window (open switch, source not yet on, disconnected port).
    const double T = i_curr.t.back() - i_curr.t.front();
    const double avg_mag = (T > 0.0) ? den / T : den;

    double WR_relErr;
    if (avg_mag > 0.1)      WR_relErr = num / den;                 // relative L1 (unchanged behaviour)
    else                    WR_relErr = (T > 0.0) ? num / T : num; // near-zero current: absolute avg change

    cout << "WR-step  " << WR_iteration
         << (avg_mag > 0.1 ? ", L1 rel current error = " : ", L1 abs current error = ")
         << WR_relErr << endl;

    return WR_relErr;
}


// Terminal-scalar WR convergence metric, port of the reference CoSimulation_WR.cpp.
// Checks, at the window-end terminal values, both the transmission mismatch
// (field vs circuit) and the iteration-to-iteration change of V and I. Each term is
// made relative when the reference magnitude exceeds 0.1, absolute otherwise.
double eval_WR_convergence_terminal(
    double V_field, double I_field,
    double V_circuit, double I_circuit,
    double V_field_last_it, double I_field_last_it,
    const unsigned WR_iteration)
{
    // I_field == I_circuit ?  (transmission condition)
    double rel_I_FC = fabs(I_field - I_circuit);
    if (fabs(I_field) > 0.1) rel_I_FC /= fabs(I_field);

    // V_field == V_circuit ?
    double rel_V_FC = fabs(V_field - V_circuit);
    if (fabs(V_field) > 0.1) rel_V_FC /= fabs(V_field);

    // I_field == I_field_last_it ?  /  V_field == V_field_last_it ?  (iteration change)
    double rel_I_it = 0.0;
    double rel_V_it = 0.0;
    if (WR_iteration > 1) {
        rel_I_it = fabs(I_field - I_field_last_it);
        if (fabs(I_field) > 0.1) rel_I_it /= fabs(I_field);
        rel_V_it = fabs(V_field - V_field_last_it);
        if (fabs(V_field) > 0.1) rel_V_it /= fabs(V_field);
    }

    const double WR_relErr = rel_I_FC + rel_V_FC + rel_I_it + rel_V_it;
    cout << "WR-step  " << WR_iteration << " (terminal), errors  "
         << rel_I_FC << " , " << rel_V_FC << " , "
         << rel_I_it << " , " << rel_V_it << endl;
    return WR_relErr;
}


// Reads Xyce .prn output produced by ".print tran V(p) V(nx) I(Vmeas)"
// Columns: Index  time  V(p)  V(nx)  I(Vmeas)
// Resamples the circuit->field quantity per coupling_mode: voltage-driven -> V(p) to
// vf_prev_k.pwl (field reads port voltage); current-driven -> I(Vmeas) to i_prev_k.pwl.
// Writes last (V(nx), I(Vmeas)) to Circuit.txt.
void ReadXyceResults(const string& filename, CircuitWaveform& circuit_raw, double t_start, double t_stop, unsigned N_xyce_eval_points)
{
    Waveform vp_raw;   // V(p): circuit -> field in voltage-driven mode
    Waveform i_raw;    // I(Vmeas): circuit -> field in current-driven mode

    circuit_raw.t.clear();
    circuit_raw.vp.clear();
    circuit_raw.vnx.clear();
    circuit_raw.i.clear();

    FILE* file = fopen(filename.c_str(), "r");
    if (!file)
        throw runtime_error(string("Could not open Xyce results file: ") + filename);

    char line[512];

    // check whether the file is empty, also skips the first line (header)
    if (!fgets(line, sizeof(line), file)) {
        fclose(file);
        throw runtime_error("ReadXyceResults: empty file");
    }

    double last_V = 0.0;
    double last_I = 0.0;
    double idx, time, Vp, Viface, I;
    unsigned step = 0;

    while (fgets(line, sizeof(line), file))
    {
        if (strncmp(line, "End", 3) == 0) break;

        if (sscanf(line, "%lg %lg %lg %lg %lg", &idx, &time, &Vp, &Viface, &I) == 5)
        {
            pushOrReplaceDuplicateTime(vp_raw, time, Vp);
            pushOrReplaceDuplicateTime(i_raw, time, I);
            pushOrReplaceDuplicateTime(circuit_raw, time, Vp, Viface, I);

            last_V = Viface;
            last_I = I;
            step++;
        }
    }

    fclose(file);

    if (step == 0)
        throw runtime_error("ReadXyceResults: no data rows read");

    Write_Terminal_results("Circuit.txt", last_V, last_I);

    // Resample the circuit->field quantity onto the uniform coupling grid.
    //   voltage-driven: V(p) -> vf_prev_k.pwl (field reads the port voltage).
    //   current-driven: I(Vmeas) -> i_prev_k.pwl (field reads the interface current).
    if (g_cfg.coupling_mode == 1) {
        Waveform i_sampled = resampleWaveformUniform(i_raw, t_start, t_stop, N_xyce_eval_points);
        writePWLFile("i_prev_k.pwl", i_sampled);
        cout << "Raw Xyce points: " << i_raw.t.size()
             << ", coupling PWL points: " << i_sampled.t.size() << " (I->field)" << endl;
    } else {
        Waveform vp_sampled = resampleWaveformUniform(vp_raw, t_start, t_stop, N_xyce_eval_points);
        writePWLFile("vf_prev_k.pwl", vp_sampled);
        cout << "Raw Xyce points: " << vp_raw.t.size()
             << ", coupling PWL points: " << vp_sampled.t.size() << " (V->field)" << endl;
    }
}

// Runs the Xyce circuit solver on the given netlist (stdout/stderr redirected to xyce_*.log).
// Throws if Xyce exits non-zero.
#ifdef _WIN32
// Open a log file as an inheritable handle for the child's stdout/stderr, replacing
// the "> xyce_stdout.log" shell redirection that cmd.exe would otherwise do for us.
static HANDLE OpenChildLog(const char* path) {
    SECURITY_ATTRIBUTES sa{};
    sa.nLength = sizeof(sa);
    sa.bInheritHandle = TRUE;
    return CreateFileA(path, GENERIC_WRITE, FILE_SHARE_READ | FILE_SHARE_WRITE,
                       &sa, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, nullptr);
}
#endif

void RunXyce(const string& filename) {
#ifdef _WIN32
    HANDLE hOut = OpenChildLog("xyce_stdout.log");
    HANDLE hErr = OpenChildLog("xyce_stderr.log");
    if (hOut == INVALID_HANDLE_VALUE || hErr == INVALID_HANDLE_VALUE) {
        if (hOut != INVALID_HANDLE_VALUE) CloseHandle(hOut);
        if (hErr != INVALID_HANDLE_VALUE) CloseHandle(hErr);
        throw runtime_error("Could not open xyce_stdout.log / xyce_stderr.log for writing.");
    }

    // CreateProcessA may modify lpCommandLine in place, so it needs a writable buffer.
    // Quote the netlist path: the working directory can sit under "Mobile Documents".
    string cmd = "Xyce \"" + filename + "\"";
    vector<char> cmdline(cmd.begin(), cmd.end());
    cmdline.push_back('\0');

    STARTUPINFOA si{};
    si.cb = sizeof(si);
    si.dwFlags = STARTF_USESTDHANDLES;
    si.hStdInput = GetStdHandle(STD_INPUT_HANDLE);
    si.hStdOutput = hOut;
    si.hStdError = hErr;
    PROCESS_INFORMATION pi{};

    // lpApplicationName = nullptr so that "Xyce" is still resolved via PATH, exactly as
    // the shell did. CREATE_NO_WINDOW suppresses the console that would otherwise flash
    // up once per WR iteration.
    BOOL ok = CreateProcessA(nullptr, cmdline.data(), nullptr, nullptr,
                             TRUE, CREATE_NO_WINDOW, nullptr, nullptr, &si, &pi);
    DWORD spawnErr = ok ? 0 : GetLastError();
    CloseHandle(hOut);
    CloseHandle(hErr);
    if (!ok) {
        throw runtime_error(
            "Could not start Xyce (CreateProcess error " + to_string(spawnErr) +
            "). Is Xyce on PATH?"
        );
    }

    WaitForSingleObject(pi.hProcess, INFINITE);
    DWORD code = 1;
    GetExitCodeProcess(pi.hProcess, &code);
    CloseHandle(pi.hThread);
    CloseHandle(pi.hProcess);

    int ret = static_cast<int>(code);
#else
    string cmd =
        string("Xyce ") + filename +
        " > xyce_stdout.log 2> xyce_stderr.log";

    int ret = system(cmd.c_str());
#endif

    if (ret != 0) {
        throw runtime_error(
            string("Xyce failed with exit code ") + to_string(ret) +
            ". Check xyce_stdout.log and xyce_stderr.log."
        );
    }
}

// Compact full-precision double->string for building netlist expressions (std::to_string loses
// small magnitudes: to_string(1e-7) == "0.000000").
static string fmtg(double x)
{
    std::ostringstream os;
    os << scientific << setprecision(10) << x;
    return os.str();
}

// Absolute end time of the whole run (= t_end). Used by the custom-spec SW element's PWL end time.
static double simEndTime()
{
    return (g_cfg.t_end > 0.0) ? g_cfg.t_end : 1.0;
}

// Emits a user-authored node-graph circuit from circuit_spec.txt. Reserved nodes:
// p (port, field attaches here) and 0 (ground); any other token is a user node. The fixed WR
// interface is appended by WriteCircuitNetlist -- the spec must not touch node nx. Well-posedness is
// the user's responsibility (full generality). Line format, '#'/'*' comment, blanks skipped:
//   <TYPE> <name> <nodeA> <nodeB> <params...>
//   R/L/C name a b val | VSIN/ISIN name a b amp f | VDC name a b val | VPULSE name a b v1 v2 td tr
// Emitted device name = type-letter + user name (valid Xyce device).
static void emitCustomTopology(ofstream& out)
{
    ifstream in("circuit_spec.txt");
    if (!in) {
        throw runtime_error("emitCustomTopology: circuit_spec.txt not found.");
    }
    int emitted = 0, lineno = 0;
    string line;
    while (std::getline(in, line)) {
        ++lineno;
        // strip inline comment (# or *)
        for (char c : {'#', '*'}) { const auto pos = line.find(c); if (pos != string::npos) line.erase(pos); }
        std::istringstream iss(line);
        vector<string> tok; string t;
        while (iss >> t) tok.push_back(t);
        if (tok.empty()) continue;
        string type = tok[0];
        for (char& ch : type) if (ch >= 'a' && ch <= 'z') ch -= 32;   // uppercase
        auto need = [&](size_t n) {
            if (tok.size() < n)
                throw runtime_error("emitCustomTopology: circuit_spec.txt line " + to_string(lineno)
                                    + " (TYPE '" + type + "') needs " + to_string(n - 1) + " fields.");
        };
        if (type == "R" || type == "L" || type == "C") {
            need(5);
            const string& nm = tok[1]; const string& a = tok[2]; const string& b = tok[3];
            const string& val = tok[4];
            if      (type == "R") out << "R" << nm << " " << a << " " << b << " " << val << "\n";
            else if (type == "L") out << "L" << nm << " " << a << " " << b << " " << val << " IC=0\n";
            else                  out << "C" << nm << " " << a << " " << b << " " << val << " IC=0\n";
        } else if (type == "VSIN" || type == "ISIN") {
            need(6);
            const char q = (type == "VSIN") ? 'V' : 'I';
            out << "B" << tok[1] << " " << tok[2] << " " << tok[3] << " " << q
                << " = { " << tok[4] << "*sin(2*pi*" << tok[5] << "*time) }\n";
        } else if (type == "VDC" || type == "IDC") {
            need(5);
            const char dev = (type == "VDC") ? 'V' : 'I';
            out << dev << tok[1] << " " << tok[2] << " " << tok[3] << " " << tok[4] << "\n";
        } else if (type == "VPULSE" || type == "IPULSE") {
            need(8);
            const char dev = (type == "VPULSE") ? 'V' : 'I';
            out << dev << tok[1] << " " << tok[2] << " " << tok[3] << " PULSE("
                << tok[4] << " " << tok[5] << " " << tok[6] << " " << tok[7] << " 0 1e30 1e30)\n";
        } else if (type == "VPWM" || type == "IPWM") {
            // Repeating PWM pulse train: <name> <a> <b> v1 v2 freq duty.
            // Emitted as a periodic Xyce PULSE: high (v2) for duty*period each period.
            need(8);
            const char dev = (type == "VPWM") ? 'V' : 'I';
            double freq = std::strtod(tok[6].c_str(), nullptr);
            double duty = std::strtod(tok[7].c_str(), nullptr);
            if (!(freq > 0.0))
                throw runtime_error("emitCustomTopology: circuit_spec.txt line " + to_string(lineno)
                                    + " (" + type + "): PWM freq must be > 0.");
            if (duty < 0.0) duty = 0.0;
            if (duty > 1.0) duty = 1.0;
            const double per  = 1.0 / freq;
            const double edge = per * 1e-3;          // small but nonzero rise/fall
            double pw = duty * per;                   // high (v2) width
            if (pw < edge)       pw = edge;           // keep 0 < pw < per for Xyce
            if (pw > per - edge) pw = per - edge;
            std::ostringstream ss; ss.setf(std::ios::scientific); ss.precision(9);
            ss << dev << tok[1] << " " << tok[2] << " " << tok[3] << " PULSE("
               << tok[4] << " " << tok[5] << " 0 " << edge << " " << edge << " " << pw << " " << per << ")";
            out << ss.str() << "\n";
        } else if (type == "VPWL" || type == "IPWL") {
            // Piecewise-linear (multi-step) source: <name> <a> <b> t1 v1 t2 v2 ... (>=1 pair).
            need(6);
            if ((tok.size() - 4) % 2 != 0)
                throw runtime_error("emitCustomTopology: circuit_spec.txt line " + to_string(lineno)
                                    + " (" + type + ") needs an even number of (time value) points.");
            // Validate each point is a number -> a clear error instead of Xyce's opaque "Cannot convert
            // 'X' to double" (e.g. a '2-3' typo for '2e-3'). SPICE magnitude suffixes (k/m/u/n/p/f/g/t)
            // are allowed as a trailing letter.
            for (size_t i = 4; i < tok.size(); ++i) {
                char* e = nullptr; std::strtod(tok[i].c_str(), &e);
                bool ok = (e != tok[i].c_str());
                if (ok && *e != '\0') { char c = *e; if (c >= 'a' && c <= 'z') c -= 32;
                    ok = (c=='T'||c=='G'||c=='K'||c=='M'||c=='U'||c=='N'||c=='P'||c=='F'); }
                if (!ok)
                    throw runtime_error("emitCustomTopology: circuit_spec.txt line " + to_string(lineno)
                        + " (" + type + "): PWL point '" + tok[i] + "' is not a number (use e.g. 2e-3, not 2-3).");
            }
            const char dev = (type == "VPWL") ? 'V' : 'I';
            out << dev << tok[1] << " " << tok[2] << " " << tok[3] << " PWL(";
            for (size_t i = 4; i < tok.size(); ++i) out << tok[i] << (i + 1 < tok.size() ? " " : "");
            out << ")\n";
        } else if (type == "SW") {
            // Time-gated switch: CLOSED (Ron) during [tclose, topen), OPEN (Roff) otherwise. Emitted as
            // a Xyce NATIVE Generic Switch (S device + .MODEL SWITCH), which smooths the R transition
            // internally (unlike the plain SPICE switch). CONTROL is a continuous 0->1 gate built from
            // the smooth ramps up/dn (transition time tau) -- Xyce warns that discontinuous CONTROL
            // hurts convergence, so we ramp it rather than hard-switch. With ON=1/OFF=0 the model
            // interpolates Roff->Ron across the full ramp (tau keeps its old trapezoidal-transition
            // meaning). topen>=tEnd (e.g. 1e30) = stays closed to the end; tclose<=0 = closed from start.
            //   SW name a b tclose topen [Ron Roff trise]
            need(6);
            const string& nm = tok[1]; const string& a = tok[2]; const string& b = tok[3];
            const double tclose = std::strtod(tok[4].c_str(), nullptr);
            const double topen  = std::strtod(tok[5].c_str(), nullptr);
            // Ron/Roff are optional: if the spec omits them, leave them OUT of the .MODEL so Xyce uses
            // its own SWITCH defaults (RON=1, ROFF=1/GMIN=1e12) rather than forcing values here.
            const bool  has_Ron  = (tok.size() > 6);
            const bool  has_Roff = (tok.size() > 7);
            const double Ron  = has_Ron  ? std::strtod(tok[6].c_str(), nullptr) : 1.0;    // Xyce default
            const double Roff = has_Roff ? std::strtod(tok[7].c_str(), nullptr) : 1.0e12; // Xyce default
            double tau        = (tok.size() > 8) ? std::strtod(tok[8].c_str(), nullptr) : 1.0e-5;
            if (tau <= 0.0) tau = 1.0e-9;
            const double tEnd = simEndTime();
            const string up = (tclose <= 0.0) ? string("1")
                : ("MIN(MAX((TIME-" + fmtg(tclose) + ")/" + fmtg(tau) + ",0),1)");
            string dn = (topen >= tEnd) ? string("0")
                : ("MIN(MAX((TIME-" + fmtg(topen) + ")/" + fmtg(tau) + ",0),1)");
            // A literal-constant CONTROL (={1 - 0}, closed for the whole run) SEGFAULTS Xyce's generic
            // switch (a control with no TIME/variable dependence). Keep it TIME-referenced: use an
            // "opens at t_end" ramp for dn -- numerically 0 over the run (switch stays closed), but not
            // a compile-time constant.
            if (up == "1" && dn == "0")
                dn = "MIN(MAX((TIME-" + fmtg(tEnd) + ")/" + fmtg(tau) + ",0),1)";
            // Generic Switch: S<nm> a b <model> CONTROL={0..1 gate}; model interpolates Roff..Ron.
            // Emit only the RON/ROFF the spec specified; otherwise let Xyce default them (RON=1, ROFF=1e12).
            out << "S" << nm << " " << a << " " << b << " SWMOD" << nm
                << " CONTROL={" << up << " - " << dn << "}\n";
            out << ".MODEL SWMOD" << nm << " SWITCH (ON=1 OFF=0";
            if (has_Ron)  out << " RON="  << fmtg(Ron);
            if (has_Roff) out << " ROFF=" << fmtg(Roff);
            out << ")\n";
        } else {
            throw runtime_error("emitCustomTopology: circuit_spec.txt line " + to_string(lineno)
                                + " unknown TYPE '" + type + "'.");
        }
        ++emitted;
    }
    if (emitted == 0)
        throw runtime_error("emitCustomTopology: circuit_spec.txt has no elements.");
}

// Generates the full circuit netlist (wr_circuit.cir) from the custom node-graph in
// circuit_spec.txt. The topology is window-invariant -> called ONCE before the WR loop (and in
// emit mode). Elements are written INLINE (no .INCLUDE of the devices) so the UI netlist parser
// (does not follow .INCLUDEs) can draw the circuit. The WR interface (Vmeas/Bfield) is fixed and
// appended after the user circuit; per-window quantities stay as {param} refs in sim_params.inc.
void WriteCircuitNetlist(const string& filename)
{
    LoadProbes("probes.txt");   // populate g_probes for the .print line (both `emit` and solve paths)

    ofstream out(filename);
    if (!out) {
        throw runtime_error("WriteCircuitNetlist: could not open " + filename);
    }
    out << scientific << setprecision(10);

    out << "Generated circuit netlist (WriteCircuitNetlist from sim_config.txt) -- DO NOT EDIT BY HAND\n";
    out << ".INCLUDE sim_params.inc\n\n";

    out << "* === CIRCUIT SIDE (custom node-graph from circuit_spec.txt) ===\n";
    emitCustomTopology(out);
    out << "\n";

    if (g_cfg.validation_mode) {
        // Validation mode: the TRUE field as REAL Xyce devices (R_FEM + L_FEM in series on the
        // port branch), replacing the behavioral matched-secant Bfield. Solved monolithically
        // (MonolithicValidationSolve) -> a reference for the coupled WR run. Vmeas is kept so the
        // print columns V(p)/V(nx)/I(Vmeas) are unchanged; with Vmeas a 0 V ammeter V(p)=V(nx), so
        // the "circuit" and "field" curves coincide by construction. Linear field only.
        out << "* === VALIDATION MODE: true field as real devices (R_FEM + L_FEM series on port) ===\n";
        if (g_cfg.nonlin_model != 0)
            out << "* NOTE: nonlin_model=" << g_cfg.nonlin_model
                << " (saturation) is IGNORED here -- Lfem is a linear inductor.\n";
        out << "Vmeas p nx 0\n";
        out << "Rfem  nx nfem " << g_cfg.R_FEM << "\n";
        out << "Lfem  nfem 0 "  << g_cfg.L_FEM << "\n\n";
    } else {
    // Fixed WR interface. Coupling direction per coupling_mode.
    // Both coupling directions use the SAME matched-secant Bfield: the field's Thevenin equivalent, a
    // voltage source V(vfprev) behind the ROM impedance Z = Rrom + Lrom/dt. That IS a Newton/secant
    // linearisation of the field's V-I characteristic V = Z_field(I) about the previous iterate
    // (V_prev, I_prev), with Z as the Jacobian estimate; the correction Z*(I(Vmeas) - V(iprev)) vanishes
    // at convergence. The device lines are identical -- only which point (V_prev, I_prev) the base sits
    // on differs, i.e. how the FEM evaluated the field (set by the FEM solver + ReadXyceResults target):
    //   voltage-driven (Dirichlet): vf_prev = V(p) [ReadXyce], i_prev = I_field [FEM]. Base = the
    //     circuit's prior port voltage (= the field's voltage only at convergence); the residual is the
    //     cross-solver current defect I_circuit - I_field, so at the window start Lrom/t_floor amplifies
    //     it -> the ~f V(p) spike.
    //   current-driven (Neumann):  vf_prev = V_field [FEM], i_prev = I(Vmeas) [ReadXyce]. Base = the
    //     field's OWN directly-computed voltage (accurate); the residual is the iteration change
    //     I_k - I_{k-1} of a single variable, so the amplified term stays small -> no window-start spike.
    if (g_cfg.coupling_mode == 1)
        out << "* === WR INTERFACE (current-driven): matched-secant Bfield "
               "(vf_prev = V_field from FEM, i_prev = I(Vmeas)) ===\n";
    else
        out << "* === WR INTERFACE (voltage-driven): matched-secant Bfield "
               "(vf_prev = V(p), i_prev = I_field from FEM) ===\n";
    out << "VFprev vfprev 0 PWL FILE \"vf_prev_k.pwl\"\n";
    out << "VIprev iprev  0 PWL FILE \"i_prev_k.pwl\"\n";
    out << "Vmeas p nx 0\n";
    // Secant denominator dt in Z = Rrom + Lrom/dt. use_t_floor selects how it is formed (all reuse the
    // t_floor .PARAM = t_floor_frac * t_window):
    //   1 = floored (default): MAX(time - t_abs_start, t_floor); the accumulated-secant denominator,
    //       guarded against the 1/0 at the exact window start (time==t_abs_start, dt=0 -> Lrom/0).
    //   0 = bare: (time - t_abs_start); Z -> infinity at the window start (test only).
    //   2 = constant: t_floor; the denominator is FIXED over the window, so Z = Rrom + Lrom/t_floor is a
    //       constant interface impedance (a fixed Robin/optimized-transmission coefficient) instead of
    //       the accumulated secant whose admittance grows across the window. Same terminal fixpoint (the
    //       correction vanishes at convergence -- see interface_form/precondition); only conditioning
    //       and the WR contraction rate change. The FEM keeps its own (BDF-1) field integration, unchanged.
    string dt_expr;
    if      (g_cfg.use_t_floor == 2) dt_expr = "t_floor";
    else if (g_cfg.use_t_floor == 1) dt_expr = "MAX(time - t_abs_start, t_floor)";
    else                             dt_expr = "(time - t_abs_start)";
    if (g_cfg.precondition == 0) {
        // CLASSICAL Gauss-Seidel (Dirichlet-Neumann) WR: NO ROM preconditioner. Drop the matched-secant
        // correction (Z*(I - i_prev) / (V - vf_prev)/Z) entirely; the circuit port is driven by a PURE
        // source carrying the field's OWN response variable (never the circuit's prior iterate, which
        // would decouple the loop). Same terminal fixpoint as the preconditioned scheme (the Z-term
        // vanishes at convergence), but plain Gauss-Seidel contraction. interface_form is ignored.
        if (g_cfg.coupling_mode == 1) {
            // current-driven: the field returns V_field (vf_prev) -> pure Dirichlet voltage source.
            out << "* precondition = 0 (classical WR): pure voltage source V(nx) = V_field, no ROM secant\n";
            out << "Bfield nx 0 V = { V(vfprev) }\n\n";
        } else {
            // voltage-driven: the field returns I_field (i_prev) -> pure Neumann current source.
            out << "* precondition = 0 (classical WR): pure current source I(nx->0) = I_field, no ROM secant\n";
            out << "Bfield nx 0 I = { V(iprev) }\n\n";
        }
    } else if (g_cfg.interface_form == 1) {
        // Norton (dual of the Thevenin): a behavioral CURRENT source with shunt G = 1/Z. Same TERMINAL
        // fixpoint, but G -> 0 at the window start -> V(p) weakly tied -> stiffer, and the interior V(p)
        // waveform differs from the Thevenin one (converges anyway; no dt-collapse observed up to ~MHz).
        out << "* interface_form = Norton: current source I = V(iprev) + (V(nx) - V(vfprev))/Z, Z = Rrom + Lrom/dt\n";
        out << "Bfield nx 0 I = {\n";
        out << "+ V(iprev) + (V(nx) - V(vfprev)) / (Rrom + Lrom/" << dt_expr << ")\n";
        out << "+ }\n\n";
    } else if (g_cfg.interface_consistency != 0) {
        // INCONSISTENT interface (Thevenin; BOTH coupling directions). Split the single secant
        // denominator Z*(I(Vmeas) - V(iprev)) so the LIVE-iterate current term uses the circuit's own
        // timestep dt_C while the LAGGED term keeps dtf (= dt_expr = time - t_abs_start, per use_t_floor).
        // The netlist device lines are identical for both directions; only what the PWLs carry differs
        // (FEM/ReadXyce side, untouched), so the SAME emission serves both -- but the shifted quantity
        // differs:
        //   voltage-driven (coupling_mode=0): V(iprev)=I_field (FEM), I(Vmeas)=live circuit current.
        //     Consistent cancels the two Lrom terms at I_C=I_F; the split leaves I_C != I_F at the
        //     fixpoint -> the terminal CURRENT is shifted (I_circuit vs I_field gap).
        //   current-driven (coupling_mode=1): V(iprev)=I(Vmeas)_{k-1} (lagged circuit current, ReadXyce),
        //     V(vfprev)=V_field (FEM base). Consistent cancels at I_C^k=I_C^{k-1}; the split leaves
        //     V_C = V_field + Lrom*(I_C-I0)*(1/dt_C - 1/dtf) at the fixpoint -> the terminal VOLTAGE is
        //     shifted (V_circuit vs V_field gap). Same mechanism, dual observable.
        // R stays a consistent difference. I0 = the carried window-start current (.PARAM I0). The
        // live-current term's dt_C is realized three ways (interface_consistency = 1/2/3):
        if (g_cfg.interface_consistency == 1) {
            // 1 = REAL RL: Rrom + Lrom as REAL Xyce devices in series on the port branch, carrying
            // I(Vmeas). The real Lromc (IC={I0}) contributes Lrom*(I_C - I0)/dt_C at Xyce's actual step;
            // Bfield then carries ONLY the lagged correction (V_field base minus the lagged-current
            // terms). Series branch is p-Vmeas-nx-Rromc-nrc-Lromc-nlc-Bfield-0, so I(Vmeas)=I(Lromc)=I_C.
            // Window-1 UIC + restart carry Lromc's state across windows like the authored circuit inductors.
            out << "* interface_consistency=1 (inconsistent, real RL): Rrom+Lrom as real devices on port; Bfield = lagged correction only\n";
            out << "* vdidtf = the field-grid (h_F) BDF-1 backward difference di_F/dt, precomputed in the field\n";
            out << "* solver into didt_field_k.pwl (NOT reconstructed here via DDT or an anchored secant). The\n";
            out << "* field-side inductive correction reads it directly, so it stays exact for >1 field eval per\n";
            out << "* window. No real inductor sits on this Bfield branch: the real Lromc above carries only the\n";
            out << "* circuit-side current I(Vmeas) (h_C derivative); Bfield supplies just the field-side term.\n";
            out << "Vdidtf vdidtf 0 PWL FILE \"didt_field_k.pwl\"\n";
            out << "Rromc nx nrc {Rrom}\n";
            out << "Lromc nrc nlc {Lrom} IC={I0}\n";
            out << "Bfield nlc 0 V = { V(vfprev) - Rrom*V(iprev) - Lrom*V(vdidtf) }\n\n";
        } else if (g_cfg.interface_consistency == 2) {
            // 2 = DDT: circuit inductive term = Lrom*DDT(I(Vmeas)), Xyce's time derivative at the circuit
            // step (I0 constant -> DDT(I_C) = DDT(I_C - I0)). Requires DDT() support in a B-source.
            out << "* interface_consistency=2 (inconsistent, DDT): circuit term Lrom*DDT(I(Vmeas)), lagged term over dtf\n";
            out << "Bfield nx 0 V = {\n";
            out << "+ V(vfprev) + Rrom*(I(Vmeas) - V(iprev))\n";
            out << "+ + Lrom*DDT(I(Vmeas))\n";
            out << "+ - Lrom*(V(iprev) - I0)/" << dt_expr << "\n";
            out << "+ }\n\n";
        } else {
            // 3 = MEAN dt: circuit inductive term over dt_print = t_window/N_xyce_samples (existing .PARAM),
            // a fixed representative mean step instead of the true adaptive dt_C.
            out << "* interface_consistency=3 (inconsistent, mean dt): circuit term over dt_print, lagged term over dtf\n";
            out << "Bfield nx 0 V = {\n";
            out << "+ V(vfprev) + Rrom*(I(Vmeas) - V(iprev))\n";
            out << "+ + Lrom*(I(Vmeas) - I0)/dt_print\n";
            out << "+ - Lrom*(V(iprev) - I0)/" << dt_expr << "\n";
            out << "+ }\n\n";
        }
    } else {
        // Thevenin (default): a behavioral VOLTAGE source that pins V(nx)=V(p).
        out << "Bfield nx 0 V = {\n";
        out << "+ V(vfprev) + (Rrom + Lrom/" << dt_expr << ") * (I(Vmeas) - V(iprev))\n";
        out << "+ }\n\n";
    }
    } // end WR-interface / validation-mode branch

    out << ".INCLUDE restart.inc\n";
    out << ".print tran V(p) V(nx) I(Vmeas)";
    for (const string& pr : g_probes) out << " " << pr;    // user probes -> extra trailing columns
    out << "\n";
    out << ".end\n";
}

// Writes sim_params.inc: the per-window .PARAM lines the netlist references (t_start/t_stop/t_abs_start,
// I0, Rrom/Lrom) plus the derived dt_print = t_window/N_xyce_samples (print cadence + interface PWL
// resolution) and t_floor = t_floor_frac*t_window (secant-denominator guard).
// Throws if t_window <= 0 or N_xyce_samples < 1.
void WriteSimParams(
    const string& filename,
    double t_start, // print start time (absolute) = t_abs_start; Xyce prints from here on
    double t_window, // window length dt_field
    double t_abs_start, // absolute window start; now used for (time - t_abs_start) in the secant terms
    double i0,
    double rrom,
    double lrom,
    const unsigned N_xyce_samples
)
{
    if (t_window <= 0.0) {
        throw runtime_error("WriteSimParams: t_window (stopping time for xyce) must be positive.");
    }

    if (N_xyce_samples < 1) {
        throw runtime_error("WriteSimParams: N_xyce_samples must be >= 1.");
    }

    // Xyce print cadence dt_print: output at every re-sampling point (see resampleWaveformUniform,
    // which reduces the Xyce points). This is the raw resolution of the interface PWL.
    const double dt_print = t_window / static_cast<double>(N_xyce_samples);

    // Floor on the accumulated secant time t_acc. ONLY a 1/0 safeguard at the exact
    // window-start eval point (time == t_abs_start): Lrom/MAX(t_acc, t_floor) instead of Lrom/t_acc.
    // The actual stability fix for the real Ls_d is NOT this floor, but that
    // Bfield uses a SMOOTH secant (no more IF(...) hard switch on the constant I0):
    // the IF branch forced an ideal current source (admittance 0) at the window start in series with
    // the real BDF Ls_d -> degenerate/stiff -> crash. The smooth secant keeps the admittance
    // finite (tiny, but >0). Floor sweep (factor 0..20) confirms: accuracy/stability-neutral.
    // The floor is now its own knob: t_floor = t_floor_frac * t_window (window-scaled).
    const double t_floor = g_cfg.t_floor_frac * t_window;

    // Time is absolute/continuous (checkpoint/restart): the tran stop time is the
    // ABSOLUTE window-end time, not the window length.
    const double t_stop_abs = t_abs_start + t_window;

    ofstream out(filename);
    out << scientific << setprecision(16);
    out << ".PARAM t_start      = " << t_start      << "\n";
    out << ".PARAM t_stop       = " << t_stop_abs   << "\n";
    out << ".PARAM t_abs_start  = " << t_abs_start  << "\n";
    out << ".PARAM I0           = " << i0           << "\n";
    out << ".PARAM Rrom         = " << rrom         << "\n";
    out << ".PARAM Lrom         = " << lrom         << "\n";
    out << ".PARAM dt_print     = " << dt_print     << "\n";
    out << ".PARAM t_floor      = " << t_floor      << "\n";
}

// Writes restart.inc: the window-specific .OPTIONS RESTART and .tran line.
// Time is absolute -> tran stop time = {t_stop} (absolute), print start = {t_abs_start}.
//   Window 1: fresh UIC transient from 0, writes checkpoints (JOB=...).
//   Window k>1: restart from committed_file (FILE=...), writes new checkpoints.
// INITIAL_INTERVAL = window length -> exactly one checkpoint at the (absolute) window end.
// Written as a literal, since .OPTIONS should not rely on {param} expansion.
void WriteRestartDirectives(
    const string& filename,
    bool first_window,
    double dt_window,
    const string& ckpt_out_prefix,
    const string& committed_file)
{
    ofstream out(filename);
    if (!out) {
        throw runtime_error("WriteRestartDirectives: could not open " + filename);
    }

    out << scientific << setprecision(16);

    if (first_window) {
        out << ".OPTIONS RESTART PACK=0 JOB=" << ckpt_out_prefix
            << " INITIAL_INTERVAL=" << dt_window << "\n";
        // Fresh transient from 0; print from {t_abs_start} (=0 in the first window).
        out << ".tran {dt_print} {t_stop} {t_abs_start} UIC\n";
    } else {
        out << ".OPTIONS RESTART FILE=" << committed_file
            << " JOB=" << ckpt_out_prefix
            << " INITIAL_INTERVAL=" << dt_window << "\n";
        // Restart sets the integrator to the checkpoint time; no UIC.
        out << ".tran {dt_print} {t_stop} {t_abs_start}\n";
    }
}

// Deletes stale checkpoint candidates <prefix>* in the working directory. Prevents an
// outdated candidate from a previous window being committed by mistake.
void ClearCheckpoints(const string& prefix)
{
    namespace fs = std::filesystem;
    for (const auto& entry : fs::directory_iterator(fs::current_path())) {
        if (!entry.is_regular_file()) continue;
        const string name = entry.path().filename().string();
        if (name.rfind(prefix, 0) == 0) { // starts with prefix
            std::error_code ec;
            fs::remove(entry.path(), ec);
        }
    }
}

// Finds the checkpoint <prefix>* with the GREATEST simulation time (= window end) and
// copies it to committed_file as the restart basis for the next time window.
// The sim-time is the suffix in the file name (Xyce: JOB + time, e.g. ckpt_out0.02).
// Max-time (instead of mtime) is robust, because on restart Xyce additionally writes a checkpoint
// at t=0 (ckpt_out0); the desired end state always has the greatest time.
void CommitCheckpoint(const string& prefix, const string& committed_file)
{
    namespace fs = std::filesystem;
    fs::path best;
    double best_time = -1.0;
    bool found = false;

    for (const auto& entry : fs::directory_iterator(fs::current_path())) {
        if (!entry.is_regular_file()) continue;
        const string name = entry.path().filename().string();
        if (name.rfind(prefix, 0) != 0) continue; // not our prefix

        // Parse the time suffix after the prefix (e.g. "0", "0.02", "4e-04"). Require the ENTIRE suffix
        // to be a clean number (*end == '\0'): a genuine Xyce checkpoint name is prefix + pure numeric
        // time, whereas an iCloud conflict copy is "ckpt_out0.0002 2" -> strtod stops at " 2" and would
        // otherwise be accepted with the SAME time as the real file, letting a stale copy win the commit
        // nondeterministically. Rejecting trailing junk keeps only the real checkpoints.
        const string suffix = name.substr(prefix.size());
        char* end = nullptr;
        const double t = std::strtod(suffix.c_str(), &end);
        if (end == suffix.c_str() || *end != '\0') continue; // non-numeric or trailing junk -> skip

        if (!found || t > best_time) {
            best = entry.path();
            best_time = t;
            found = true;
        }
    }

    if (!found) {
        throw runtime_error(
            "CommitCheckpoint: no checkpoint file '" + prefix +
            "*' produced by Xyce. Restart write failed?");
    }
    const fs::path& newest = best;

    std::error_code ec;
    fs::copy_file(newest, committed_file,
                  fs::copy_options::overwrite_existing, ec);
    if (ec) {
        throw runtime_error("CommitCheckpoint: failed to copy " +
                            newest.string() + " -> " + committed_file +
                            ": " + ec.message());
    }
}

// Linear resampler onto a uniform grid with N_intervals+1 nodes over [t_start, t_stop].
// Purpose: reduce the adaptive Xyce output to a fixed grid before it is reused as a PWL source
// (V(vfprev) or V(iprev), depending on coupling_mode) in the next WR iteration (see the call in
// ReadXyceResults).
Waveform resampleWaveformUniform(
    const Waveform& raw,
    double t_start,
    double t_stop,
    unsigned N_intervals)
{
    if (raw.t.size() != raw.y.size() || raw.t.empty()) {
        throw runtime_error("Invalid raw waveform.");
    }

    if (N_intervals < 1) {
        throw runtime_error("N_intervals must be >= 1.");
    }

    Waveform sampled;

    size_t k = 0;

    for (unsigned j = 0; j <= N_intervals; ++j) {
        const double tau =
            t_start
            + (t_stop - t_start) * static_cast<double>(j)
            / static_cast<double>(N_intervals);

        while (k + 1 < raw.t.size() && raw.t[k + 1] < tau) {
            ++k;
        }

        double value = raw.y.back();

        if (tau <= raw.t.front()) {
            value = raw.y.front();
        } else if (tau >= raw.t.back()) {
            value = raw.y.back();
        } else {
            const double t0 = raw.t[k];
            const double t1 = raw.t[k + 1];
            const double y0 = raw.y[k];
            const double y1 = raw.y[k + 1];

            const double alpha = (tau - t0) / (t1 - t0);
            value = y0 + alpha * (y1 - y0);
        }

        sampled.push(tau, value);
    }

    return sampled;
}

// Appends a circuit waveform (Index TIME V(P) V(NX) I(VMEAS)) to an open .prn file with a running
// global_index and absolute time (t_abs_start offset). skip_first_point drops the shared window-edge
// sample when stitching consecutive windows.
void appendCircuitWaveformXyceStyle(
    FILE* file,
    const CircuitWaveform& wf,
    double t_abs_start,
    unsigned long& global_index,
    bool skip_first_point)
{
    if (!file) {
        throw runtime_error("appendCircuitWaveformXyceStyle: file is null.");
    }

    const size_t start = skip_first_point ? 1 : 0;

    for (size_t n = start; n < wf.t.size(); ++n) {
        const double t_abs = t_abs_start + wf.t[n];

        fprintf(
            file,
            "%-10lu %-17.8e %-17.8e %-17.8e %-17.8e\n",
            global_index++,
            t_abs,
            wf.vp[n],
            wf.vnx[n],
            wf.i[n]
        );
    }

    fflush(file);
}

// Appends the field waveforms (Index TIME V(FIELD) I(FIELD)) to an open .prn file with a running
// global_index and absolute time. vf and i must share the same time grid (checked). skip_first_point
// drops the shared window-edge sample between consecutive windows.
void appendFieldWaveformXyceStyle(
    FILE* file,
    const Waveform& vf,
    const Waveform& i,
    double t_abs_start,
    unsigned long& global_index,
    bool skip_first_point)
{
    if (!file) {
        throw runtime_error("appendFieldWaveformXyceStyle: file is null.");
    }

    if (vf.t.size() != i.t.size()) {
        throw runtime_error("appendFieldWaveformXyceStyle: size mismatch.");
    }

    const size_t start = skip_first_point ? 1 : 0;

    for (size_t n = start; n < vf.t.size(); ++n) {
        if (std::abs(vf.t[n] - i.t[n]) > 1e-12) {
            throw runtime_error("appendFieldWaveformXyceStyle: time grid mismatch.");
        }

        const double t_abs = t_abs_start + vf.t[n];

        fprintf(
            file,
            "%-10lu %-17.8e %-17.8e %-17.8e\n",
            global_index++,
            t_abs,
            vf.y[n],
            i.y[n]
        );
    }

    fflush(file);
}
