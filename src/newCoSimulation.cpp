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


// Global tunable configuration (defaults defined in SimConfig). Overwritten by LoadConfig.
SimConfig g_cfg;

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
        else if (key == "source_kind")                      g_cfg.source_kind = (unsigned)val;
        else if (key == "frequency")                        g_cfg.frequency = val;
        else if (key == "amplitude")                        g_cfg.amplitude = val;
        else if (key == "step_v_initial")                   g_cfg.step_v_initial = val;
        else if (key == "step_v_final")                     g_cfg.step_v_final = val;
        else if (key == "step_delay")                       g_cfg.step_delay = val;
        else if (key == "step_rise")                        g_cfg.step_rise = val;
        else if (key == "R_series")                         g_cfg.R_series = val;
        else if (key == "L_series")                         g_cfg.L_series = val;
        else if (key == "C_series")                         g_cfg.C_series = val;
        else if (key == "circuit_kind")                     g_cfg.circuit_kind = (unsigned)val;
        else if (key == "switch_backend")                   g_cfg.switch_backend = (unsigned)val;
        else if (key == "switch_t1")                        g_cfg.switch_t1 = val;
        else if (key == "switch_t2")                        g_cfg.switch_t2 = val;
        else if (key == "switch_C")                         g_cfg.switch_C = val;
        else if (key == "switch_R")                         g_cfg.switch_R = val;
        else if (key == "switch_Ron")                       g_cfg.switch_Ron = val;
        else if (key == "switch_Roff")                      g_cfg.switch_Roff = val;
        else if (key == "switch_trise")                     g_cfg.switch_trise = val;
        else if (key == "time_mode")                        g_cfg.time_mode = (unsigned)val;
        else if (key == "t_end")                            g_cfg.t_end = val;
        else if (key == "N_field_windows")                  g_cfg.N_field_windows = (unsigned)val;
        else if (key == "N_periods")                        g_cfg.N_periods = (unsigned)val;
        else if (key == "N_field_steps_per_source_period")  g_cfg.N_field_steps_per_source_period = (unsigned)val;
        else if (key == "N_field_eval_intervals")           g_cfg.N_field_eval_intervals = (unsigned)val;
        else if (key == "N_xyce_coupling_intervals")        g_cfg.N_xyce_coupling_intervals = (unsigned)val;
        else if (key == "WRmaxSteps")                       g_cfg.WRmaxSteps = (unsigned)val;
        else if (key == "WR_tolerance")                     g_cfg.WR_tolerance = val;
        else if (key == "wr_convergence_method")            g_cfg.wr_convergence_method = (unsigned)val;
        else if (key == "coupling_mode")                    g_cfg.coupling_mode = (unsigned)val;
        else if (key == "reconstruct_mode")                 g_cfg.reconstruct_mode = (unsigned)val;
        else cout << "LoadConfig: unknown key '" << key << "' ignored." << endl;
    }
    return true;
}


struct Waveform {
    vector<double> t;
    vector<double> y;

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

struct CircuitWaveform {
    vector<double> t;
    vector<double> vp;
    vector<double> vnx;
    vector<double> i;

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

// Schreibt 2-Punkt PWL als lineare Rampe von (t_start, value_start) mit Steigung `slope`.
// slope=0 (default) → konstante PWL. slope ≠ 0 sorgt für konsistente Anfangsableitung
// an V(iprev)/V(vfprev) am Fensteranfang (passend zur dIdt0-Verwendung in Biface IF-Guard).
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

    // Schaltungs-Netzliste aus g_cfg generieren (Topologie fensterinvariant -> einmal hier).
    // Ueberschreibt wr_circuit.cir mit Quelle (source_kind) + serieller R/L/C-Kette + fixer
    // WR-Schnittstelle. Pro-Fenster variable Parameter kommen weiter aus sim_params.inc.
    WriteCircuitNetlist("wr_circuit.cir");

    // Datei für Feldlösung an Synchronisationszeitpunkten (entspricht hier Feldschritten)
    FILE* file_Field = fopen("Field_solution.prn", "w");
    fprintf(file_Field, "Index       TIME              V(FIELD)          I(FIELD)\n");
    fflush(file_Field);

    // Datei für die rekonstruierte/extrapolierte Feldlösung, d.h. die Felddaten die Xyce eingelesen hat
    FILE* file_Field_waveform = fopen("Field_waveform_solution.prn", "w");
    fprintf(file_Field_waveform, "Index       TIME              V(FIELD)          I(FIELD)\n");
    fflush(file_Field_waveform);
    
    FILE* file_Circuit = fopen("Circuit_solution.prn", "w");
    fprintf(file_Circuit, "Index       TIME              V(P)              V(NX)             I(VMEAS)\n");
    fflush(file_Circuit);

    FILE* file_WR_error = fopen("WR_error.txt", "w");
    fprintf(file_WR_error, "   Time, WR_TotalRelErr, N_iterations, Converged \n");
    fflush(file_WR_error);

    unsigned long global_field_index = 0;
    unsigned long global_field_waveform_index = 0;
    unsigned long global_circuit_index = 0;

    // All tunable parameters come from g_cfg (loaded from sim_config.txt; see SimConfig).
    // The impedance values (L_ROM, R_ROM) are a reduced-order model of the field domain used
    // inside the Xyce Biface interface condition; they need not match the "true" FEM values.
    const double L_ROM = g_cfg.L_ROM;
    const double R_ROM = g_cfg.R_ROM;
    // Circuit voltage source Bsrc: single source of truth for amplitude/frequency, written via
    // WriteSimParams to sim_params.inc (f_src, amp_src) and referenced by the netlist Bsrc, so
    // dVdt_0 below cannot drift from the netlist.
    const double V_src_amplitude = g_cfg.amplitude;
    const double Frequency = g_cfg.frequency;

    // Time stepping (the circuit side is handled dynamically by Xyce). Duration per time_mode:
    //   0 (source periods):  T = N_periods/f,  N_steps_field = N_fsp*N_periods (sinusoidal sources).
    //   1 (absolute end):    T = t_end,        N_steps_field = N_field_windows (step/switch/custom).
    const unsigned N_periods = g_cfg.N_periods;
    const unsigned N_field_steps_per_source_period = g_cfg.N_field_steps_per_source_period;
    unsigned N_steps_field;
    double dt_field;
    if (g_cfg.time_mode == 1) {
        N_steps_field = (g_cfg.N_field_windows >= 1) ? g_cfg.N_field_windows : 1u;
        if (g_cfg.t_end <= 0.0) throw runtime_error("MasterProcess: t_end must be > 0 in absolute time_mode.");
        dt_field = g_cfg.t_end / N_steps_field;
    } else {
        N_steps_field = N_field_steps_per_source_period * N_periods;
        dt_field = (1.0 / Frequency) / N_field_steps_per_source_period;
    }

    // WR parameters
    const unsigned WRmaxSteps = g_cfg.WRmaxSteps;
    const double WR_tolerance = g_cfg.WR_tolerance;
    // Fixed coupling-grid resolution, decoupled from Xyce's adaptive step. Controls TWO things
    // with one value: (1) Xyce print cadence dt_print = dt_field/N (WriteSimParams) and
    // (2) resample resolution of V(p) into vf_prev_k.pwl. See voltage_driven_refactor.tex.
    const unsigned N_xyce_coupling_intervals = g_cfg.N_xyce_coupling_intervals;
    // FEM evaluation intervals per window (FEM-PWL has N+1 points). =1 → 2-point ramp → single
    // linear extrapolation of the field current per window. Higher → multi-rate, piecewise linear.
    const unsigned N_field_eval_intervals = g_cfg.N_field_eval_intervals;

    // Initialisierung
    double V_field = 0;
    double I_field = 0;
    double I0 = 0.0;
    double V0 = 0.0;
    double dIdt_0 = 0.0;
    // dV/dt der Portspannung am Fensteranfang (Seed-Slope für vf_prev_k.pwl). Wird am Fensterende
    // aus der konvergierten V(p)-Waveform getragen (topologie-agnostisch, korrekt auch mit Rs/Ls).
    // Fenster 1: quellen-abhaengige Anfangssteigung (bei I0=0,dIdt0=0 ist V(p)≈Vsrc). Die sinus-
    // Formel gilt NUR fuer source_kind 0; fuer Strom-/Step-Quellen ist sie falsch und verlangsamt die
    // WR-Konvergenz des ersten Fensters (die restlichen Fenster tragen die Steigung aus der Loesung).
    double dVdt_0;
    switch (g_cfg.source_kind) {
        case 0: // sinusoidale SPANNUNGsquelle (Legacy): Quellensteigung amplitude*2*pi*f (Regressions-Anker).
            dVdt_0 = V_src_amplitude * 2.0 * M_PI * Frequency * cos(2.0 * M_PI * Frequency * 0.0);
            break;
        default: // Strom-/Step-Quellen: neutraler Seed 0. V(p) startet ~0 (Ls blockiert den Anfangsstrom,
                 // beim Step faellt die volle Quellspannung zunaechst ueber Ls ab -> dV(p)/dt(0)≈0). Ein
                 // Steigungs-Seed aus der reinen Rampensteigung ueberschiesst und KOSTET WR-Iterationen.
            dVdt_0 = 0.0;
            break;
    }

    //fprintf(file_Field, "%-10lu %-17.8e %-17.8e %-17.8e\n", global_field_index++, 0.0, V_field, I_field);

    // Äußere Schleife der Zeitfenster
    // hier ist Feldintervall = WR-Zeitfenster
    for (unsigned step_field = 1; step_field <= N_steps_field; step_field++) {
        // nichtlineares Verhalten des RL-Gliedes (dummy)
        //R_ROM *= 1.2;
        //L_ROM *= 1.05;

        const double t_start = (step_field - 1) * dt_field; // brauchen wir vor allem hier (in cpp), Xyce simuliert immer von 0 bis Stopzeitpunkt wenn man es aufruft. In xyce nutzen wir es nur für die Phase der restlichen Schaltung (z.B. Vsrc wenn Schaltung nur Spannungsquelle ist)
        const double t_stop = step_field * dt_field;

        WriteSimParams("sim_params.inc", 0.0, dt_field, t_start, I0, R_ROM, L_ROM, Frequency, V_src_amplitude, dIdt_0, g_cfg.R_series, g_cfg.L_series, N_xyce_coupling_intervals);

        // Fenster-spezifische Restart/.tran-Direktiven generieren (von wr_circuit.cir inkludiert).
        // Fenster 1: frischer UIC-Transient ab 0; Fenster k>1: Restart aus "restart_state".
        WriteRestartDirectives("restart.inc", step_field == 1, dt_field, "ckpt_out", "restart_state");

        // initialisiere mit zuletzt akzeptierten Werten. Zeitstempel sind ABSOLUT ([t_start, t_stop]),
        // da Xyce die PWL-FILE-Quellen an der absoluten Simulationszeit auswertet (Restart startet
        // bei t_start, nicht bei 0).
        WriteInitialPwl("vf_prev_k.pwl", t_start, t_stop, V0, dVdt_0);
        // i_prev_k.pwl: lineare Rampe mit Steigung dIdt_0, passend zur akkumulierten Sekante im FEM.
        WriteInitialPwl("i_prev_k.pwl",   t_start, t_stop, I0, dIdt_0);

        // Veraltete Checkpoint-Kandidaten dieses Prefixes entfernen, damit CommitCheckpoint
        // nach der WR-Schleife garantiert den frisch erzeugten Kandidaten dieses Fensters waehlt.
        ClearCheckpoints("ckpt_out");

        CircuitWaveform circuit_sol; // circuit solution array (nur für Ausgabe/Visualisierung gebraucht)

        unsigned WR_iteration;
        double WR_rel_Error = 1.0;
        bool WR_converged = false;
        Waveform i_prev_last_iter; // i_m^(k-1) für L1-Konvergenzkriterium; leer am Fensteranfang
        double V_field_last_iter = 0.0; // für terminal-skalar Kriterium (Referenz CoSimulation_WR.cpp)
        double I_field_last_iter = 0.0;
        // WR-Iterationsschleife
        for (WR_iteration = 1; WR_iteration <= WRmaxSteps; WR_iteration++) {
            circuit_sol = CircuitWaveform{};

            // Circuit Solver aufrufen: Biface liefert I(Vmeas) = INTERFACE_condition(V(p), V(vfprev), V(iprev)).
            // ReadXyceResults schreibt V(p)-Waveform nach vf_prev_k.pwl für FEM-Eingang.
            RunXyce("wr_circuit.cir");
            ReadXyceResults("wr_circuit.cir.prn", circuit_sol, t_start, t_stop, N_xyce_coupling_intervals);

            // FEM Solver aufrufen (dummy). Kopplungsrichtung per coupling_mode:
            //   0 voltage-driven: liest vf_prev_k.pwl (V(p)), schreibt i_prev_k.pwl (I_field).
            //   1 current-driven: liest i_prev_k.pwl (I(Vmeas)), schreibt vf_prev_k.pwl (V_field).
            // schreibt Endwerte in Field.txt
            if (g_cfg.coupling_mode == 1)
                FEM_solver_current_driven_waveform(I0, V0, N_field_eval_intervals);
            else
                FEM_solver_voltage_driven_waveform(I0, N_field_eval_intervals);
            Read_Terminal_results("Field.txt", V_field, I_field);

            // i_prev = Feldstrom-Waveform dieser Iteration (FEM-Ausgang)
            Waveform i_prev = readPWLFile("i_prev_k.pwl");
            std::cout << "PWL points: " << i_prev.t.size() << std::endl;

            // Prüfe Konvergenzkriterium. Methode wählbar via g_cfg.wr_convergence_method:
            //   0 = Waveform-L1 des Feldstroms (dieses Codebase)
            //   1 = terminal-skalar (Referenz CoSimulation_WR.cpp): Feld-vs-Schaltung +
            //       Iterations-zu-Iterations-Aenderung der Terminalwerte V,I.
            bool can_converge;
            if (g_cfg.wr_convergence_method == 1) {
                double V_circuit, I_circuit;
                Read_Terminal_results("Circuit.txt", V_circuit, I_circuit);
                WR_rel_Error = eval_WR_convergence_terminal(
                    V_field, I_field, V_circuit, I_circuit,
                    V_field_last_iter, I_field_last_iter, WR_iteration);
                // Referenz erlaubt Konvergenz ab Iteration 1 (FC-Terme allein koennen reichen).
                can_converge = true;
            } else {
                // L1-rel-Norm von (i_m^(k) - i_m^(k-1)) / i_m^(k); braucht >=2 Iterationen
                // (Iteration 1 liefert Sentinel 1.0, keine Vorgaenger-Waveform).
                WR_rel_Error = eval_WR_convergence(i_prev, i_prev_last_iter, WR_iteration);
                can_converge = (WR_iteration > 1);
            }

            if (can_converge && WR_rel_Error < WR_tolerance) {
                WR_converged = true;
                break;
            }

            // Werte dieser Iteration für die nächste Konvergenzprüfung aufheben
            i_prev_last_iter = i_prev;
            V_field_last_iter = V_field;
            I_field_last_iter = I_field;
        }

        // Screen Output
        //
        // Resultate in Datei schreiben

        PRINT(t_stop);
        PRINT(WR_rel_Error);

        fprintf(file_WR_error, "%12.3e ,  %12.3e ,  %d ,  %d \n", t_stop, WR_rel_Error, WR_iteration, int(WR_converged));
        fflush(file_WR_error);


        if (!WR_converged) {
            throw runtime_error("WR did not converge on current time window. Stop simulation.");
        } else {
            PRINT("WR converged!");
            PRINT(V_field);
            PRINT(I_field);

            // Feldlösung an den Endpunkten (hier gleich den Synchronisationspunkten) im Xyce Format speichern
            fprintf(file_Field, "%-10lu %-17.8e %-17.8e %-17.8e\n", global_field_index++, t_stop, V_field, I_field);
            fflush(file_Field);

            // Konvergierten Endzustand als Restart-Basis fuers naechste Fenster sichern.
            // (Alle WR-Iterationen schrieben denselben Kandidaten ckpt_out<t_stop>; der
            // neueste ist der konvergierte.)
            CommitCheckpoint("ckpt_out", "restart_state");

            // die Xycelösung des Zeitfensters an die vorigen anhängen.
            // Waveforms tragen jetzt ABSOLUTE Zeit → Offset 0.0 (kein erneutes Verschieben).
            const bool skip_first_point = (step_field > 1);
            appendCircuitWaveformXyceStyle(file_Circuit, circuit_sol, 0.0, global_circuit_index, skip_first_point);

            // auch die Feld-WAVEFORMS (also nicht nur Endpunkte) speichern wir im Xyce Format ab
            // dafür lesen wir die konvergierten Waveforms ein (letzte Iteration)
            Waveform vf_conv = readPWLFile("vf_prev_k.pwl"); // port voltage V(p) from Xyce (coupling grid)
            Waveform i_conv  = readPWLFile("i_prev_k.pwl");  // field current I_field from FEM (FEM grid)

            // Beide auf gemeinsames FEM-Eval-Grid (N_field_eval_intervals+1) re-sampeln,
            // damit appendFieldWaveformXyceStyle übereinstimmende Zeitstempel sieht.
            // (vf_conv hat N_xyce_coupling_intervals+1 Knoten, i_conv hat N_field_eval_intervals+1)
            Waveform vf_endpoints = resampleWaveformUniform(
                vf_conv, t_start, t_stop, N_field_eval_intervals
            );
            Waveform i_endpoints = resampleWaveformUniform(
                i_conv, t_start, t_stop, N_field_eval_intervals
            );

            // und hängen sie an die bisherigen Zeitfenster an (Waveforms bereits absolut → Offset 0.0)
            appendFieldWaveformXyceStyle(file_Field_waveform, vf_endpoints, i_endpoints, 0.0, global_field_waveform_index, skip_first_point);

            // Festhalten der Endwerte für die Anfangswerte des nächsten Zeitfensters
            V0 = vf_conv.y.back(); // Portspannung am Fensterende (= V(p) bei t=dt_field)
            I0 = I_field;          // Feldstrom am Fensterende (FEM-Ausgang)

            // dI/dt am Fensterende = dI/dt am Anfang des nächsten Fensters (Stetigkeit)
            const size_t n_i = i_conv.t.size();
            if (n_i >= 2) {
                const double dt_end = i_conv.t[n_i - 1] - i_conv.t[n_i - 2];
                dIdt_0 = (i_conv.y[n_i - 1] - i_conv.y[n_i - 2]) / dt_end;
            } else {
                dIdt_0 = 0.0;
            }

            // dV/dt der Portspannung am Fensterende → Seed-Slope des nächsten Fensters.
            // Aus der konvergierten V(p)-Waveform gemessen (topologie-agnostisch): korrekt auch
            // wenn V(p) ≠ Vsrc (serielle Rs/Ls), wo die reine Quellensteigung falsch wäre.
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

}

// FEM-Solver (voltage-driven): nimmt Portspannungs-Waveform V_p(t) aus vf_prev_k.pwl
// und berechnet Feldstrom I_field(t) via akkumulierter Sekante. Konsistent mit Biface-Ausdruck
// in wr_circuit.cir. Multi-rate: FEM arbeitet auf gröberem Grid als Xyce-intern.
void FEM_solver_voltage_driven_waveform(double I_win_start, unsigned N_field_eval_intervals)
{
    // Reads the port voltage waveform stored by ReadXyceResults from V(p).
    const Waveform vport = readPWLFile("vf_prev_k.pwl");

    // Zeit ist jetzt absolut (vport.t.front() = absoluter Fensteranfang, nicht 0). Der
    // FEM-Solver ist offset-agnostisch: t_acc wird unten als (t[j] - front) gemessen, also
    // korrekt unabhaengig vom absoluten Startzeitpunkt. Daher keine "startet bei 0"-Pruefung.

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

    // Multi-rate: FEM re-sampelt Portspannung auf groebers FEM-Grid.
    const Waveform V_eval = resampleWaveformUniform(
        vport, t_start, t_end, N_field_eval_intervals
    );
    const size_t N = V_eval.t.size(); // = N_field_eval_intervals + 1

    // KRITISCH: I_field-Berechnung muss exakt dasselbe akkumulierte Sekanten-Modell wie
    // Biface in wr_circuit.cir verwenden, um konsistenten WR-Fixpunkt zu garantieren.
    // Biface: I_c = V(iprev) + (V(p) - V(vfprev)) / (Rrom + Lrom/time)
    // FEM-Herleitung: V_p = R_FEM*I_f + L_FEM*(I_f - I0)/t_acc
    //   → I_f = (V_p + L_FEM*I0/t_acc) / (R_FEM + L_FEM/t_acc)  für t_acc > 0
    //   → I_f = I_win_start                                     für t_acc = 0
    const double t_win_start = V_eval.t.front();

    Waveform current;
    for (size_t j = 0; j < N; ++j) {
        double I_j;
        if (j == 0) {
            // Singulärer Punkt t=0: Anfangsbedingung.
            I_j = I_win_start;
        } else {
            // Akkumulierte Sekante: V_p = R*I_f + L*(I_f - I0)/t_acc
            const double t_acc = V_eval.t[j] - t_win_start;
            // Linear closed-form solution; also the Newton seed for the nonlinear case.
            I_j = (V_eval.y[j] + L_FEM * I_win_start / t_acc)
                / (R_FEM + L_FEM / t_acc);

            if (saturating) {
                // Implicit residual with the accumulated flux secant (same window-start
                // reference as the linear case, so consistent with the Xyce Biface grid):
                //   g(I) = R_FEM*I + (lambda(I) - lambda(I_win_start))/t_acc - V_p = 0
                //   g'(I) = R_FEM + L(I)/t_acc
                const double lam0 = flux(I_win_start);
                for (int it = 0; it < 50; ++it) {
                    const double g  = R_FEM * I_j + (flux(I_j) - lam0) / t_acc - V_eval.y[j];
                    const double gp = R_FEM + L_dyn(I_j) / t_acc;
                    const double dI = g / gp;
                    I_j -= dI;
                    if (std::fabs(dI) <= 1e-12 + 1e-10 * std::fabs(I_j)) break;
                }
            }
        }
        current.push(V_eval.t[j], I_j);
    }

    // This file is used by Xyce as V(iprev) in the next WR iteration.
    writePWLFile("i_prev_k.pwl", current);

    // Terminal-Werte am Fensterende für Konvergenz-Propagation ins nächste Fenster.
    Write_Terminal_results("Field.txt", vport.y.back(), current.y.back());
}

// Current-driven (Neumann) field solver. Reads the interface current I(t) (i_prev_k.pwl, circuit
// output) and returns the field voltage V(t) = R_FEM*I + dlambda/dt in vf_prev_k.pwl for the circuit's
// Bfield voltage source. reconstruct_mode:
//   0 pointwise (default): V computed at every field-eval point from I(t) -- symmetric to the
//      voltage-driven solver, needs only I0 (V_field_last_time unused); follows the current's curve.
//   1 linear: replace the interior with a straight ramp from V_field_last_time (previous window end)
//      to this window's end value (the colleague's 2-point reconstruction; uses V0).
//   2 average: 0.5*(linear + const) -- window-start raised to 0.5*(V_field_last_time + V_end); uses V0.
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

    // Pointwise field voltage V(t) = R_FEM*I + dlambda/dt (needs only I0). dlambda/dt via the
    // accumulated secant from the window start; at the first point (t_acc=0 -> 0/0) use a forward
    // difference over the first interval (the same 0/0 guard the voltage-driven solver / t_floor use).
    Waveform vfield;
    for (size_t j = 0; j < N; ++j) {
        const double t = I_eval.t[j];
        double dl_dt;
        if (j == 0) {
            const double dt01 = I_eval.t[1] - t_win_start;
            dl_dt = (dt01 > 0.0) ? dlam(I_eval.y[1], I_win_start) / dt01 : 0.0;
        } else {
            dl_dt = dlam(I_eval.y[j], I_win_start) / (t - t_win_start);
        }
        vfield.push(t, R_FEM * I_eval.y[j] + dl_dt);
    }
    const double V_field_end = vfield.y.back();

    // Optional colleague reconstruction: replace the interior with a straight line (modes 1/2 use V0).
    if (g_cfg.reconstruct_mode == 1 || g_cfg.reconstruct_mode == 2) {
        const double V_start = (g_cfg.reconstruct_mode == 2)
                               ? 0.5 * (V_field_last_time + V_field_end)
                               : V_field_last_time;
        const double dt_win = t_end - t_win_start;
        for (size_t j = 0; j < N; ++j) {
            const double frac = (dt_win > 0.0) ? (I_eval.t[j] - t_win_start) / dt_win : 0.0;
            vfield.y[j] = V_start + frac * (V_field_end - V_start);
        }
    }
    const double I_end = I_eval.y.back();

    writePWLFile("vf_prev_k.pwl", vfield);                     // field -> circuit: the field voltage
    Write_Terminal_results("Field.txt", V_field_end, I_end);  // (V_field_end, I_end)
}

// WR-Konvergenzkriterium auf Basis der Stromwaveforms benachbarter Iterationen:
//
//   ∫|i_m^(k)(t) - i_m^(k-1)(t)| dt   /   ∫|i_m^(k)(t)| dt   ≤   WR_tolerance
//
// Integration via Trapezregel über das Kopplungsraster (i_curr und i_prev_iter teilen sich
// dasselbe uniforme Zeitgrid aus resampleWaveformUniform). Wegen Transmissionsbedingung
// i_m^(k) = i_c^(k) wird hier der konvergierte/gerade berechnete Schaltungsstrom verwendet.
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

	if (den <= 0.0) {
		throw runtime_error("eval_WR_convergence: zero denominator (current waveform identically zero?).");
	}

	const double WR_relErr = num / den;
	cout << "WR-step  " << WR_iteration
	     << ", L1 rel current error = " << WR_relErr << endl;

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
// Writes V(p) resampled to vf_prev_k.pwl (port voltage for FEM voltage-driven input).
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

    // prüfe ob Datei leer ist, überspringt auch die erste Zeile (Header)
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

void RunXyce(const string& filename) {
    string cmd =
        string("Xyce ") + filename +
        " > xyce_stdout.log 2> xyce_stderr.log";

    int ret = system(cmd.c_str());

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

// Absolute end time of the whole run (matches MasterProcess): absolute mode -> t_end,
// periods mode -> N_periods/frequency. Used by the switch-control PWL end time.
static double simEndTime()
{
    if (g_cfg.time_mode == 1) return (g_cfg.t_end > 0.0) ? g_cfg.t_end : 1.0;
    return (g_cfg.frequency > 0.0) ? (double)g_cfg.N_periods / g_cfg.frequency : 1.0;
}

// --- Switch-throw emitter (increment 2) -------------------------------------------------------
// Emits a two-terminal connection between n1 and n2 that is CLOSED (low R) during [ta, tb) and
// OPEN (high R) otherwise. Backend per g_cfg.switch_backend:
//   0 behavioral: resistor with an IF(time) value (Ron inside the window, Roff outside) -- no
//                 .MODEL, robust, plays well with the WR windows.
//   1 native:     Xyce voltage-controlled switch S + a PWL control source = 1 in [ta,tb)
//                 (shared .MODEL SWMOD emitted once by emitSwitchTopology).
// Open-ended-to-end: pass tb >= tEnd (1e30). From-start: pass ta <= 0.
static void emitThrow(ofstream& out, const string& tag, const string& n1, const string& n2,
                      double ta, double tb, double tEnd)
{
    // Trapezoidal "closed" gate g(t) in [0,1]: rises over [ta, ta+tau] (or =1 from the start if
    // ta<=0), falls over [tb, tb+tau] (or stays 1 to the end if tb>=tEnd). tau = switch_trise.
    const double tau = (g_cfg.switch_trise > 0.0) ? g_cfg.switch_trise : 1.0e-9;
    if (g_cfg.switch_backend == 1) {
        // native VC switch + PWL control (ramps 0<->1 over tau at each edge -> S device slides
        // between ROFF and RON continuously, no instantaneous topology jump).
        out << "S" << tag << " " << n1 << " " << n2 << " ctrl" << tag << " 0 SWMOD\n";
        out << "Vctrl" << tag << " ctrl" << tag << " 0 PWL(";
        if (ta <= 0.0) out << "0 1";
        else           out << "0 0 " << ta << " 0 " << (ta + tau) << " 1";
        if (tb < tEnd) out << " " << tb << " 1 " << (tb + tau) << " 0 " << tEnd << " 0";
        else           out << " " << tEnd << " 1";
        out << ")\n";
    } else {
        // behavioral resistor R = Roff + (Ron-Roff)*g(t), g a clamped trapezoid (continuous, so
        // the integrator steps through the throw instead of colliding with a discontinuity).
        //   up = 1                          if ta<=0   else clamp((t-ta)/tau, 0, 1)
        //   dn = 0                          if tb>=tEnd else clamp((t-tb)/tau, 0, 1)
        //   g  = up - dn
        const string up = (ta <= 0.0) ? string("1")
            : ("MIN(MAX((TIME-" + fmtg(ta) + ")/" + fmtg(tau) + ",0),1)");
        const string dn = (tb >= tEnd) ? string("0")
            : ("MIN(MAX((TIME-" + fmtg(tb) + ")/" + fmtg(tau) + ",0),1)");
        out << "R" << tag << " " << n1 << " " << n2 << " R={" << g_cfg.switch_Roff
            << " + (" << g_cfg.switch_Ron << "-" << g_cfg.switch_Roff << ")*("
            << up << " - " << dn << ")}\n";
    }
}

// Emits the switch-circuit side (circuit_kind 1/2/3) at port p (see SimConfig::circuit_kind).
// The field/interface (Vmeas, Bfield) stays common and is emitted by WriteCircuitNetlist.
static void emitSwitchTopology(ofstream& out)
{
    const double tEnd = simEndTime();
    const double t1 = g_cfg.switch_t1;
    const double t2 = g_cfg.switch_t2;

    if (g_cfg.switch_backend == 1) {
        out << ".MODEL SWMOD VSWITCH(RON=" << g_cfg.switch_Ron << " ROFF=" << g_cfg.switch_Roff
            << " VON=0.5 VOFF=0.4)\n";
    }

    switch (g_cfg.circuit_kind) {
        case 1: // #4 three-way, sine U, cap C: drive[0,t1) -> freewheel[t1,t2) -> open[t2,inf)
            out << "Bemf p a V = { amp_src*sin(2*pi*f_src*time) }\n";   // U: + at p, - at a
            out << "Csw w 0 " << g_cfg.switch_C << " IC=0\n";           // C: wiper -> gnd
            emitThrow(out, "drv", "w", "a", 0.0, t1,   tEnd);          // pos1 drive     (w<->U-)
            emitThrow(out, "fw",  "w", "p", t1,  t2,   tEnd);          // pos2 freewheel (w<->p)
            // [t2,inf): both throws open -> pos0 (C isolated, field open).
            break;
        case 2: // #5 two-way, DC U, cap C: drive[0,t1) -> freewheel[t1,inf)
            out << "Vemf p a {amp_src}\n";                              // U_DC: + at p, - at a
            out << "Csw w 0 " << g_cfg.switch_C << " IC=0\n";
            emitThrow(out, "drv", "w", "a", 0.0, t1,   tEnd);          // drive     (w<->U-)
            emitThrow(out, "fw",  "w", "p", t1,  1e30, tEnd);          // freewheel (w<->p)
            break;
        case 3: // #6 two-way, AC vs R: AC-drive[0,t1) -> R-damp[t1,inf)
            out << "Bemf p bac V = { amp_src*sin(2*pi*f_src*time) }\n"; // V_AC: + at p, - at bac
            out << "Rload p br " << g_cfg.switch_R << "\n";            // R: p -> br
            emitThrow(out, "ac", "bac", "0", 0.0, t1,   tEnd);        // ground V_AC bottom (drive)
            emitThrow(out, "rd", "br",  "0", t1,  1e30, tEnd);        // ground R bottom (damp)
            break;
        default:
            throw runtime_error("emitSwitchTopology: unknown circuit_kind");
    }
}

// Emits a user-authored node-graph circuit (circuit_kind==4) from circuit_spec.txt. Reserved nodes:
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
        throw runtime_error("emitCustomTopology: circuit_kind=4 (custom) but circuit_spec.txt not found.");
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
        } else {
            throw runtime_error("emitCustomTopology: circuit_spec.txt line " + to_string(lineno)
                                + " unknown TYPE '" + type + "'.");
        }
        ++emitted;
    }
    if (emitted == 0)
        throw runtime_error("emitCustomTopology: circuit_spec.txt has no elements.");
}

// Generiert die vollstaendige Schaltungs-Netzliste (wr_circuit.cir) aus g_cfg. Topologie ist
// fensterinvariant -> EINMAL vor der WR-Schleife (und im emit-Modus) aufgerufen. Elemente werden
// INLINE geschrieben (kein .INCLUDE der Devices), damit der UI-Netzlisten-Parser (folgt keinen
// .INCLUDEs) die Schaltung zeichnen kann. Pro-Fenster variable Groessen (amp_src/f_src/Rs/Ls/...)
// bleiben als {param}-Referenzen in sim_params.inc; C_series und Step-Parameter werden als Literale
// geschrieben (konfig-konstant). Die WR-Schnittstelle (Vmeas/Bfield) ist fix und unveraendert.
void WriteCircuitNetlist(const string& filename)
{
    ofstream out(filename);
    if (!out) {
        throw runtime_error("WriteCircuitNetlist: could not open " + filename);
    }
    out << scientific << setprecision(10);

    out << "Generated circuit netlist (WriteCircuitNetlist from sim_config.txt) -- DO NOT EDIT BY HAND\n";
    out << ".INCLUDE sim_params.inc\n\n";

    if (g_cfg.circuit_kind == 4) {
        // --- Custom node-graph (increment 3a): user-authored circuit from circuit_spec.txt. ---
        out << "* === CIRCUIT SIDE (custom node-graph from circuit_spec.txt) ===\n";
        emitCustomTopology(out);
        out << "\n";
    } else if (g_cfg.circuit_kind != 0) {
        // --- Switch topology (increment 2): #4/#5/#6 attach their own network at port p. ---
        out << "* === CIRCUIT SIDE (switch topology circuit_kind=" << g_cfg.circuit_kind
            << ", backend=" << (g_cfg.switch_backend ? "native-S" : "behavioral") << ") ===\n";
        emitSwitchTopology(out);
        out << "\n";
    } else {
    // --- Simple source (increment 1): source_kind + series R/L/C chain source->port. ---
    // Aktive serielle Elemente auf dem Pfad Quelle->Port (nur nonzero). Reihenfolge R,L,C.
    vector<char> series;
    if (g_cfg.R_series != 0.0) series.push_back('R');
    if (g_cfg.L_series != 0.0) series.push_back('L');
    if (g_cfg.C_series != 0.0) series.push_back('C');
    // Ohne serielle Elemente sitzt die Quelle direkt auf dem Port p (Einweg-Toy: V(p)==Vsrc).
    const string hot = series.empty() ? string("p") : string("s");

    out << "* === CIRCUIT SIDE (source_kind=" << g_cfg.source_kind << ") ===\n";
    switch (g_cfg.source_kind) {
        case 1: // sinusoidale STROMquelle: treibt den Schleifenstrom direkt in den Port (Vorzeichen
                // via I(Vmeas) pruefen). HINWEIS: serielle R/L/C sind bei einer Stromquelle physikalisch
                // sinnlos (der Strom ist erzwungen; serielle Elemente floaten nur den Quellknoten) UND
                // eine ideale Stromquelle in Reihe mit L ist entartet (Nulldurchgang -> singulaere
                // Jacobi -> dt-Kollaps). Daher: bare Quelle direkt auf p (Preset P2 setzt R/L/C=0).
            out << "Bemf 0 " << hot << " I = { amp_src*sin(2*pi*f_src*time) }\n";
            break;
        case 2: // Step/Rampen-SPANNUNGsquelle: einzelne steigende Flanke (tf=0, pw/per gross)
            out << "Vemf " << hot << " 0 PULSE("
                << g_cfg.step_v_initial << " " << g_cfg.step_v_final << " "
                << g_cfg.step_delay << " " << g_cfg.step_rise << " 0 1e30 1e30)\n";
            break;
        case 0: // sinusoidale SPANNUNGsquelle (Legacy-Default / Regressions-Anker)
        default:
            out << "Bemf " << hot << " 0 V = { amp_src*sin(2*pi*f_src*time) }\n";
            break;
    }
    out << "\n";

    if (!series.empty()) {
        out << "* series R/L/C  source->port\n";
        string node = "s";
        for (size_t i = 0; i < series.size(); ++i) {
            const string next = (i + 1 == series.size()) ? string("p") : ("cm" + to_string(i));
            switch (series[i]) {
                case 'R': out << "Rs_d " << node << " " << next << " {Rs}\n";       break;
                case 'L': out << "Ls_d " << node << " " << next << " {Ls} IC=0\n";   break;
                case 'C': out << "Cs_d " << node << " " << next << " "
                              << g_cfg.C_series << " IC=0\n";                        break;
            }
            node = next;
        }
        out << "\n";
    }
    }

    // Feste WR-Schnittstelle. Kopplungsrichtung per coupling_mode.
    if (g_cfg.coupling_mode == 1) {
        // Current-driven (Neumann): the field returns its voltage V_field (from the FEM, in
        // vf_prev_k.pwl); the circuit's Bfield is a plain voltage source = V_field (no secant, no
        // Lrom/t_floor amplifier). The circuit's I(Vmeas) is read by the driver and fed to the FEM.
        out << "* === WR INTERFACE (current-driven): field voltage source, ammeter ===\n";
        out << "VFprev vfprev 0 PWL FILE \"vf_prev_k.pwl\"\n";
        out << "Vmeas p nx 0\n";
        out << "Bfield nx 0 V = { V(vfprev) }\n\n";
    } else {
        // Voltage-driven (Dirichlet): THEVENIN matched-secant field ROM (default).
        out << "* === WR INTERFACE (voltage-driven): prev waveforms, ammeter, matched-secant Bfield ===\n";
        out << "VFprev vfprev 0 PWL FILE \"vf_prev_k.pwl\"\n";
        out << "VIprev iprev  0 PWL FILE \"i_prev_k.pwl\"\n";
        out << "Vmeas p nx 0\n";
        out << "Bfield nx 0 V = {\n";
        out << "+ V(vfprev) + (Rrom + Lrom/MAX(time - t_abs_start, t_floor)) * (I(Vmeas) - V(iprev))\n";
        out << "+ }\n\n";
    }

    out << ".INCLUDE restart.inc\n";
    out << ".print tran V(p) V(nx) I(Vmeas)\n";
    out << ".end\n";
}

void WriteSimParams(
    const string& filename,
    double t_start, // Print-Startzeit (absolut) = t_abs_start; Xyce gibt ab hier aus
    double t_window, // Fensterlaenge dt_field
    double t_abs_start, // absoluter Fensteranfang; jetzt für (time - t_abs_start) in den Sekanten-Termen
    double i0,
    double rrom,
    double lrom,
    double f_src,   // Frequenz der Schaltungs-Spannungsquelle (Bsrc); single source of truth
    double amp_src, // Amplitude der Schaltungs-Spannungsquelle (Bsrc)
    double dIdt0,   // Anfangs-dI/dt (für Ls-Sekanten-Guard bei time<=eps_t)
    double r_series, // serielle Kopplungsimpedanz Quelle→Port (two-way coupling)
    double l_series,
    const unsigned N_coupling_intervals
)
{
    if (t_window <= 0.0) {
        throw runtime_error("WriteSimParams: t_window (stopping time for xyce) must be positive.");
    }

    if (N_coupling_intervals < 1) {
        throw runtime_error("WriteSimParams: N_coupling_intervals must be >= 1.");
    }
    
    // Für xyce output/print times (dt_print)
    // hier: output bei jedem re-sampling Punkt (siehe resampleWaveformUniform, welcher die xyce-Punkte reduziert)
    const double h_coupling = t_window / static_cast<double>(N_coupling_intervals); 
    
    // Zum vermeiden von Singularitäten in den Ableitungen
    const double eps_t = 1.0e-15;

    // Boden auf die akkumulierte Sekantenzeit t_acc = h_coupling. NUR eine 1/0-Absicherung
    // am exakten Fensterstart-Eval-Punkt (time == t_abs_start): Lrom/MAX(t_acc, t_floor)
    // statt Lrom/t_acc. Der eigentliche Stabilitaets-Fix fuer das echte Ls_d ist NICHT
    // dieser Boden, sondern dass Biface eine GLATTE Sekante nutzt (kein IF(...)-Hartschalter
    // mehr auf den konstanten I0): der IF-Zweig erzwang am Fensterstart eine ideale
    // Stromquelle (Admittanz 0) in Reihe mit dem realen BDF-Ls_d → entartet/steif → Crash.
    // Die glatte Sekante haelt die Admittanz endlich (winzig, aber >0). Boden-Sweep
    // (Faktor 0..20) bestaetigt: der Boden-Wert ist accuracy/stability-neutral.
    const double t_floor = h_coupling;

    // Zeit ist absolut/kontinuierlich (Checkpoint/Restart): tran-Stoppzeit ist der
    // ABSOLUTE Fensterende-Zeitpunkt, nicht die Fensterlaenge.
    const double t_stop_abs = t_abs_start + t_window;

    ofstream out(filename);
    out << scientific << setprecision(16);
    out << ".PARAM t_start      = " << t_start      << "\n";
    out << ".PARAM t_stop       = " << t_stop_abs   << "\n";
    out << ".PARAM t_abs_start  = " << t_abs_start  << "\n";
    out << ".PARAM I0           = " << i0           << "\n";
    out << ".PARAM Rrom         = " << rrom         << "\n";
    out << ".PARAM Lrom         = " << lrom         << "\n";
    out << ".PARAM eps_t        = " << eps_t        << "\n";
    out << ".PARAM dt_print     = " << h_coupling   << "\n";
    out << ".PARAM f_src        = " << f_src        << "\n";
    out << ".PARAM amp_src      = " << amp_src      << "\n";
    out << ".PARAM dIdt0        = " << dIdt0        << "\n";
    out << ".PARAM Rs           = " << r_series     << "\n";
    out << ".PARAM Ls           = " << l_series     << "\n";
    out << ".PARAM t_floor      = " << t_floor      << "\n";
}

// Schreibt restart.inc: die fenster-spezifische .OPTIONS RESTART und .tran Zeile.
// Zeit ist absolut → tran-Stoppzeit = {t_stop} (absolut), Print-Start = {t_abs_start}.
//   Fenster 1: frischer UIC-Transient ab 0, schreibt Checkpoints (JOB=...).
//   Fenster k>1: Restart aus committed_file (FILE=...), schreibt neue Checkpoints.
// INITIAL_INTERVAL = Fensterlaenge → genau ein Checkpoint am (absoluten) Fensterende.
// Als Literal geschrieben, da .OPTIONS sich nicht auf {param}-Expansion verlassen soll.
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
        // Frischer Transient ab 0; Print ab {t_abs_start} (=0 im ersten Fenster).
        out << ".tran {dt_print} {t_stop} {t_abs_start} UIC\n";
    } else {
        out << ".OPTIONS RESTART FILE=" << committed_file
            << " JOB=" << ckpt_out_prefix
            << " INITIAL_INTERVAL=" << dt_window << "\n";
        // Restart setzt Integrator auf die Checkpoint-Zeit; kein UIC.
        out << ".tran {dt_print} {t_stop} {t_abs_start}\n";
    }
}

// Loescht alte Checkpoint-Kandidaten <prefix>* im Arbeitsverzeichnis. Verhindert, dass
// ein veralteter Kandidat aus einem frueheren Fenster faelschlich committed wird.
void ClearCheckpoints(const string& prefix)
{
    namespace fs = std::filesystem;
    for (const auto& entry : fs::directory_iterator(fs::current_path())) {
        if (!entry.is_regular_file()) continue;
        const string name = entry.path().filename().string();
        if (name.rfind(prefix, 0) == 0) { // beginnt mit prefix
            std::error_code ec;
            fs::remove(entry.path(), ec);
        }
    }
}

// Sucht den Checkpoint <prefix>* mit der GROESSTEN Simulationszeit (= Fensterende) und
// kopiert ihn nach committed_file als Restart-Basis fuers naechste Zeitfenster.
// Die Sim-Zeit steht als Suffix im Dateinamen (Xyce: JOB + Zeit, z.B. ckpt_out0.02).
// Max-Zeit (statt mtime) ist robust, weil Xyce beim Restart zusaetzlich einen Checkpoint
// bei t=0 schreibt (ckpt_out0); der gewollte Endzustand hat immer die groesste Zeit.
void CommitCheckpoint(const string& prefix, const string& committed_file)
{
    namespace fs = std::filesystem;
    fs::path best;
    double best_time = -1.0;
    bool found = false;

    for (const auto& entry : fs::directory_iterator(fs::current_path())) {
        if (!entry.is_regular_file()) continue;
        const string name = entry.path().filename().string();
        if (name.rfind(prefix, 0) != 0) continue; // nicht unser Prefix

        // Zeit-Suffix hinter dem Prefix parsen (z.B. "0", "0.02", "4e-04").
        const string suffix = name.substr(prefix.size());
        char* end = nullptr;
        const double t = std::strtod(suffix.c_str(), &end);
        if (end == suffix.c_str()) continue; // kein numerisches Suffix → ueberspringen

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

// Linearer Resampler auf uniformes Grid mit N_intervals+1 Stuetzstellen ueber [t_start, t_stop].
// Zweck: Xyce-Adaptiv-Output auf festes Raster reduzieren bevor er als V(iprev) PWL-Quelle in der
// naechsten WR-Iteration wiederverwendet wird (siehe Aufruf in ReadXyceResults).
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
