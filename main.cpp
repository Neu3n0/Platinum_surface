#include "md.h"

int main() {
	Space* space = new Space;
	int amount_steps;
	std::cin >> amount_steps;
	Print_config(*space, amount_steps);
	space->Init_space();
	Print_atoms(*space);
	Print_mol(*space);

	/******************************************* Macro *******************************************/
	std::ofstream Kout;
	Kout.open("C:\\Drive\\Code\\MolDynamics\\md_v1.4\\Output\\Kin_energy.dat", std::ios::app);
	std::ofstream Pout;
	Pout.open("C:\\Drive\\Code\\MolDynamics\\md_v1.4\\Output\\Pot_energy.dat", std::ios::app);
	std::ofstream Fout;
	Fout.open("C:\\Drive\\Code\\MolDynamics\\md_v1.4\\Output\\Full_energy.dat", std::ios::app);
	std::ofstream Rout;
	Rout.open("C:\\Drive\\Code\\MolDynamics\\md_v1.4\\Output\\Rms_vel.dat", std::ios::app);
	std::ofstream Tout;
	Tout.open("C:\\Drive\\Code\\MolDynamics\\md_v1.4\\Output\\Temp.dat", std::ios::app);
	std::ofstream fouten;
	fouten.open("C:\\Drive\\Code\\MolDynamics\\md_v1.4\\Output\\get_energy.txt", std::ios::app);
	/**********************************************************************************************/

	auto start = steady_clock::now();
	for (int step = 100000; step < amount_steps + 100000; ++step) {
		if (step % 100 == 0) {
			auto tmp = steady_clock::now();
			std::cout << step << '\t' << duration_cast<milliseconds>(tmp - start).count() / 1000 << std::endl;
		}
		space->MDStep();
		Get_energy(fouten, *space, Kout, Pout, Fout, Rout, Tout, step);		
		if (step % 100 == 0) {
			WriteVTK(space);
		}
	}
	auto finish = steady_clock::now();
	std::cout << "Program calculation time: " << duration_cast<milliseconds>(finish - start).count() << " ms\n";

	Print_last(*space);

	fouten.close();
	Kout.close();
	Pout.close();
	Fout.close();
	Rout.close();
	Tout.close();

	delete space;
	_CrtDumpMemoryLeaks();
	return 0;
}