#include <iostream>
#include "../include/Header.h"

#pragma region main
int main( int No_Arguments, char* Arguments[  ] )
{

    std::cout << "Started Waveform Relaxaion Toy Programm " << std::endl;
	printf(" last code-compilation:  %s %s\n\n", __TIME__, __DATE__);

	// Optional config-file path as first argument (default: sim_config.txt)
	const string config_path = (No_Arguments > 1) ? Arguments[1] : "sim_config.txt";
	LoadConfig(config_path);

	//CALL MASTER PROCESS
	MasterProcess();

	cout << endl <<"Finished WR toy successfully" << endl;
	return 0; 
}
#pragma endregion main


