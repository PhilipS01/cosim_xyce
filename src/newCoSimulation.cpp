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
        else if (key == "frequency")                        g_cfg.frequency = val;
        else if (key == "amplitude")                        g_cfg.amplitude = val;
        else if (key == "R_series")                         g_cfg.R_series = val;
        else if (key == "L_series")                         g_cfg.L_series = val;
        else if (key == "N_periods")                        g_cfg.N_periods = (unsigned)val;
        else if (key == "N_field_steps_per_source_period")  g_cfg.N_field_steps_per_source_period = (unsigned)val;
        else if (key == "N_field_eval_intervals")           g_cfg.N_field_eval_intervals = (unsigned)val;
        else if (key == "N_xyce_coupling_intervals")        g_cfg.N_xyce_coupling_intervals = (unsigned)val;
        else if (key == "WRmaxSteps")                       g_cfg.WRmaxSteps = (unsigned)val;
        else if (key == "WR_tolerance")                     g_cfg.WR_tolerance = val;
        else if (key == "wr_convergence_method")            g_cfg.wr_convergence_method = (unsigned)val;
        else if (key == "bfield_deriv")                     g_cfg.bfield_deriv = (unsigned)val;
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

    // Time stepping (the circuit side is handled dynamically by Xyce)
    const unsigned N_periods = g_cfg.N_periods;
    const unsigned N_field_steps_per_source_period = g_cfg.N_field_steps_per_source_period;
    const unsigned N_steps_field = N_field_steps_per_source_period * N_periods;
    const double dt_field = (1 / Frequency) / N_field_steps_per_source_period;

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
    // Fenster 1: Quellen-Steigung (bei I0=0,dIdt0=0 ist V(p)≈Vsrc).
    double dVdt_0 = V_src_amplitude * 2.0 * M_PI * Frequency * cos(2.0 * M_PI * Frequency * 0.0);

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

            // FEM Solver aufrufen (dummy, voltage-driven):
            // liest vf_prev_k.pwl (Portspannung V(p) von Xyce)
            // schreibt i_prev_k.pwl (Feldstrom I_field via akkumulierter Sekante)
            // schreibt Endwerte in Field.txt
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
    Waveform vp_raw;

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
            pushOrReplaceDuplicateTime(circuit_raw, time, Vp, Viface, I);

            last_V = Viface;
            last_I = I;
			step++;
		}
	}

	fclose(file);

	if (step == 0)
		throw runtime_error("ReadXyceResults: no data rows read");

    // V(p) auf uniformes Kopplungsraster re-sampeln → vf_prev_k.pwl als FEM-Eingang (Portspannung).
    Waveform vp_sampled = resampleWaveformUniform(
        vp_raw,
        t_start,
        t_stop,
        N_xyce_eval_points
    );

    Write_Terminal_results("Circuit.txt", last_V, last_I);
    writePWLFile("vf_prev_k.pwl", vp_sampled);
    cout << "Raw Xyce points: " << vp_raw.t.size()
         << ", coupling PWL points: " << vp_sampled.t.size()
         << endl;
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
    // Bfield inductive-term denominator selector (0 = accumulated secant /t_acc, 1 = fixed FD /t_floor).
    out << ".PARAM use_fd       = " << g_cfg.bfield_deriv << "\n";
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
