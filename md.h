#include "header.h"

class Atom;
class Molecule_N2;
class Space;
class Cell;

class Atom {
public:
	Atom() : Kinetic_en(0), LJ_potential(0), m(0), type(0) {}
	double coord[3];
	double vel[3];
	double power[3];
	double LJ_potential;
	double Kinetic_en;
	double m;
	int type;
	Atom* at2;
	void Coord_shift(Space* space);		//Coordinates shift
	void Power_shift(Atom* mol_prob, const double* shift);		//Power shift
	void Veloc_shift(Space* space);		//Velocitie shift
	int Change_Cell(Space* space, int i0, int j0, int k0, int l);		//Change cell
};

class Molecule_N2 {
public:
	Molecule_N2() : KX_potential(0), KX_power(0) {}
	Atom* atom[2];
	double KX_power;
	double KX_potential;
};

class Cell {
public:
	Cell() : amount_atoms(0), atoms() {}
	~Cell() {
		while (amount_atoms) {
			delete atoms[amount_atoms - 1];
			--amount_atoms;
		}
	}
	Atom* atoms[MAX_ATOMS_IN_CELL];
	int amount_atoms;
};

class Space {
public:
	Space();
	double spacesize, cellsize, T, p, dt, Kin_En, Pot_En, Energy;
	int total_mol_N2, total_at_Pt, amount_cells;
	double m_N, m_Pt;
	Cell cells[MAX_CELLS][MAX_CELLS][MAX_CELLS];
	Molecule_N2 N_2[MAX_MOLS];
	void Init_space();		//Initialize molecules
	void MDStep();		//MDStep
};

void Print_config(const Space& space, const int& amount_steps);		//Print surface parameters
void Print_atoms(const Space& space);		//Print atoms
void Print_last(const Space& space);		//Print last
void Print_mol(const Space& space);			//Print molec
double LJ_F(const double& r);		//Lennard jones power
double LJ_P(const double& r);		//Lennard jones potential
double KX_P(const double& r);		//KX potential
double KX_F(const double& r);		//KX power
void SetNullMacro(Space* space);		//Set null macro
int WriteVTK(Space* space);		//VTK
void Get_energy(std::ofstream& fouten, const Space& space, std::ofstream& Kout, std::ofstream& Pout, std::ofstream& Fout, std::ofstream& Rout, std::ofstream& Tout, const int& step);		//Get energy
double Kin_En(Space* space, Atom* mol);		//Kinetic energy 
double Temp(Space& space);
double Rms_Vel(const Space& space);