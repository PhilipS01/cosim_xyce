#include <iostream>
#include "../include/Header.h"

#pragma region main
int main( int No_Arguments, char* Arguments[  ] )
{

    std::cout << "Started Waveform Relaxaion Toy Programm " << std::endl;
	printf(" last code-compilation:  %s %s\n\n", __TIME__, __DATE__);

	// Modi:
	//   main                 -> Solve (Standard-Config sim_config.txt)
	//   main <config>        -> Solve mit anderer Config-Datei
	//   main emit [config]   -> NUR die Netzliste wr_circuit.cir aus der Config generieren, kein Solve
	//                           (von der UI vor dem Zeichnen der Schaltung aufgerufen).
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
	MasterProcess();

	cout << endl <<"Finished WR toy successfully" << endl;
	return 0;
}
#pragma endregion main


