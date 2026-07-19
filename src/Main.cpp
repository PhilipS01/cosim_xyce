#include <iostream>
#include "../include/Header.h"

#pragma region main
// Entry point. Loads the config (default sim_config.txt, or the path in argv), then either emits
// ONLY the netlist (`emit` mode, for the UI to draw before solving) or runs the full WR
// co-simulation via MasterProcess. Returns 0 on success.
int main( int No_Arguments, char* Arguments[  ] )
{

    std::cout << "Started Waveform Relaxaion Toy Programm " << std::endl;
    printf(" last code-compilation:  %s %s\n\n", __TIME__, __DATE__);

    // Modes:
    //   main                 -> solve (default config sim_config.txt)
    //   main <config>        -> solve with a different config file
    //   main emit [config]   -> ONLY generate the netlist wr_circuit.cir from the config, no solve
    //                           (called by the UI before drawing the circuit).
    const bool emit_only = (No_Arguments > 1) && (string(Arguments[1]) == "emit");
    const string config_path =
        emit_only ? ((No_Arguments > 2) ? Arguments[2] : "sim_config.txt")
                  : ((No_Arguments > 1) ? Arguments[1] : "sim_config.txt");
    LoadConfig(config_path);

    if (emit_only) {
        WriteCircuitNetlist("wr_circuit.cir");
        cout << "Emitted wr_circuit.cir from " << config_path << endl;
        return 0;
    }

    //CALL MASTER PROCESS
    // validation_mode=1: monolithic reference solve (true field as real R_FEM/L_FEM devices, one
    // Xyce transient, no WR) to validate the coupled run against. Else the full WR co-simulation.
    if (g_cfg.validation_mode)
        MonolithicValidationSolve();
    else
        MasterProcess();

    cout << endl <<"Finished WR toy successfully" << endl;
    return 0;
}
#pragma endregion main


