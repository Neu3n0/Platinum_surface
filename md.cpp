#include "md.h"

Space::Space() : N_2() {
	ifstream fin;
	fin.open("config.txt");
	fin >> this->total_mol_N2 >> total_at_Pt >> this->T >> this->p >> dt >> this->spacesize >> this->cellsize;
	amount_cells = static_cast<int>(this->spacesize / this->cellsize);
	m_N = (0.001 / 6.022) / 1000.0;	//кг * 10^(-20)
	m_Pt = 3.24e-5;
	m_N = 1.673e-7;
	Kin_En = 0;
	Pot_En = 0;
	Energy = 0;
	fin.close();
}

void Space::Init_space() {
	ifstream fin;
	fin.open("init_atoms.txt");
	double coord[3]{ 0 };
	double v[3]{ 0 };

	//N2
	for (int l = 0; l < this->total_mol_N2; ++l) {
		for (int at = 0; at < 2; ++at) {
			fin >> coord[0] >> coord[1] >> coord[2] >> v[0] >> v[1] >> v[2];
			int i = static_cast<int>(coord[0] / this->cellsize);
			int j = static_cast<int>(coord[1] / this->cellsize);
			int k = static_cast<int>(coord[2] / this->cellsize);
			if ((i >= 0) && (j >= 0) && (k >= 0) && (i < this->amount_cells) && (j < this->amount_cells) && (k < this->amount_cells)) {
				Atom* atom = new Atom; 
				atom->coord[0] = coord[0];
				atom->coord[1] = coord[1];
				atom->coord[2] = coord[2];
				atom->vel[0] = v[0];
				atom->vel[1] = v[1];
				atom->vel[2] = v[2];
				atom->m = m_N;
				atom->type = 2;
				this->cells[i][j][k].atoms[this->cells[i][j][k].amount_atoms] = atom;
				this->N_2[l].atom[at] = atom;
				this->cells[i][j][k].amount_atoms++;
			}
			else {
				printf("Incorrect data\n");
			}
		}
		this->N_2[l].atom[0]->at2 = this->N_2[l].atom[1];
		this->N_2[l].atom[1]->at2 = this->N_2[l].atom[0];
	}
	//Pt
	for (int l = 0; l < this->total_at_Pt; ++l) {
		fin >> coord[0] >> coord[1] >> coord[2] >> v[0] >> v[1] >> v[2];
		int i = static_cast<int>(coord[0] / this->cellsize);
		int j = static_cast<int>(coord[1] / this->cellsize);
		int k = static_cast<int>(coord[2] / this->cellsize);
		if ((i >= 0) && (j >= 0) && (k >= 0) && (i < this->amount_cells) && (j < this->amount_cells) && (k < this->amount_cells)) {
			Atom* atom = new Atom; 
			atom->coord[0] = coord[0];
			atom->coord[1] = coord[1];
			atom->coord[2] = coord[2];
			atom->vel[0] = v[0];
			atom->vel[1] = v[1];
			atom->vel[2] = v[2];
			atom->m = m_Pt;
			atom->type = 1;
			this->cells[i][j][k].atoms[this->cells[i][j][k].amount_atoms] = atom;
			this->cells[i][j][k].amount_atoms++;
		}
		else {
			printf("Incorrect data\n");
		}		
	}
	fin.close();

	for (int i = 0; i < this->amount_cells; ++i) {
		for (int j = 0; j < this->amount_cells; ++j) {
			for (int k = 0; k < this->amount_cells; ++k) {
				for (int l = 0; l < this->cells[i][j][k].amount_atoms; ++l) {
					Atom* atom = this->cells[i][j][k].atoms[l];
					for (int n = 0; n < 3; ++n) {
						atom->power[n] = 0;
					}
				}
			}
		}
	}
}

void Space::MDStep() {
	//Coordinates and Cell shift

	Atom* atom;
	for (int i = 0; i < this->amount_cells; ++i) {
		for (int j = 0; j < this->amount_cells; ++j) {
			for (int k = 0; k < this->amount_cells; ++k) {
				for (int l = 0; l < this->cells[i][j][k].amount_atoms; ++l) {
					atom = this->cells[i][j][k].atoms[l];
					atom->Coord_shift(this);
					if (atom->Change_Cell(this, i, j, k, l) == 1) {
						--l;
					}
				}
			}
		}
	}

	//Set null macro
	SetNullMacro(this);

	//Power shift
	Atom* at1;
	Atom* at2;
	for (int i = 0; i < this->amount_cells; ++i) {
		for (int j = 0; j < this->amount_cells; ++j) {
			for (int k = 0; k < this->amount_cells; ++k) {
				for (int l = 0; l < this->cells[i][j][k].amount_atoms; ++l) {
					for (int i2 = -1; i2 < 2; ++i2) {
						for (int j2 = -1; j2 < 2; ++j2) {
							for (int k2 = -1; k2 < 2; ++k2) {
								int ii = i + i2, jj = j + j2, kk = k + k2;
								double shift[3] = { 0, 0, 0 };
								if (ii >= this->amount_cells) { ii = 0; shift[0] = this->spacesize; }
								if (jj >= this->amount_cells) { jj = 0; shift[1] = this->spacesize; }
								if (kk >= this->amount_cells) { kk = 0; shift[2] = this->spacesize; }
								if (ii < 0) { ii = this->amount_cells - 1; shift[0] = -(this->spacesize); }
								if (jj < 0) { jj = this->amount_cells - 1; shift[1] = -(this->spacesize); }
								if (kk < 0) { kk = this->amount_cells - 1; shift[2] = -(this->spacesize); }
								shift[2] = 0;
								for (int l2 = 0; l2 < this->cells[ii][jj][kk].amount_atoms; ++l2) {
									if (!((l2 == l) && (ii == i) && (jj == j) && (kk == k))) {
										at1 = this->cells[i][j][k].atoms[l];
										at2 = this->cells[ii][jj][kk].atoms[l2];
										if (at1->at2 != at2) {
											at1->Power_shift(at2, shift);
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	
	//KX power shift
	for (int i = 0; i < this->total_mol_N2; ++i) {		
		Molecule_N2 mol = this->N_2[i];
		at1 = mol.atom[0];
		at2 = mol.atom[1]; 
		double shift[3] = { 0 };
		if (static_cast<int>(at1->coord[0] / this->cellsize) == 0 && static_cast<int>(at2->coord[0] / this->cellsize) == this->amount_cells - 1) { shift[0] = - this->spacesize; }
		if (static_cast<int>(at2->coord[0] / this->cellsize) == 0 && static_cast<int>(at1->coord[0] / this->cellsize) == this->amount_cells - 1) { shift[0] = this->spacesize; }
		if (static_cast<int>(at1->coord[1] / this->cellsize) == 0 && static_cast<int>(at2->coord[1] / this->cellsize) == this->amount_cells - 1) { shift[1] = - this->spacesize; }
		if (static_cast<int>(at2->coord[1] / this->cellsize) == 0 && static_cast<int>(at1->coord[1] / this->cellsize) == this->amount_cells - 1) { shift[1] = this->spacesize; }
		if (static_cast<int>(at1->coord[2] / this->cellsize) == 0 && static_cast<int>(at2->coord[2] / this->cellsize) == this->amount_cells - 1) { shift[2] = - this->spacesize; }
		if (static_cast<int>(at2->coord[2] / this->cellsize) == 0 && static_cast<int>(at1->coord[2] / this->cellsize) == this->amount_cells - 1) { shift[2] = this->spacesize; }
		double r = 0;
		for (int a = 0; a < 3; ++a) {
			r += pow((at1->coord[a] - at2->coord[a] - shift[a]), 2);
		}
		r = sqrt(r);
		mol.KX_potential = KX_P(r);
		//Добавить в потенциальную энергию!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		mol.KX_power = KX_F(r);
		for (int n = 0; n < 3; ++n) {
			at1->power[n] += (at1->coord[n] - at2->coord[n] - shift[n]) / r * mol.KX_power;
			at2->power[n] += (at2->coord[n] - at1->coord[n] + shift[n]) / r * mol.KX_power;
		}
	}

	//Velocities shift
	for (int i = 0; i < this->amount_cells; ++i) {
		for (int j = 0; j < this->amount_cells; ++j) {
			for (int k = 0; k < this->amount_cells; ++k) {
				for (int l = 0; l < this->cells[i][j][k].amount_atoms; ++l) {
					atom = this->cells[i][j][k].atoms[l];
					atom->Veloc_shift(this);
				}
			}
		}
	}

	//Energy
	this->Pot_En /= 2;
	this->Energy = this->Pot_En + this->Kin_En;

	Temp(*this);
}

void Atom::Coord_shift(Space* space) {
	for (int n = 0; n < 3; ++n) {
		this->coord[n] += (this->vel[n] * space->dt);
		if (this->coord[n] < 0 && n != 2) {
			this->coord[n] += space->spacesize;
		}
		if (this->coord[n] > space->spacesize && n != 2) {
			this->coord[n] -= space->spacesize;
		}
	}
}

void Atom::Power_shift(Atom* mol_prob, const double* shift) {
	double r = 0;
	for (int a = 0; a < 3; ++a) {
		r = r + pow(this->coord[a] - mol_prob->coord[a] - shift[a], 2);
	}
	r = sqrt(r);
	double modforce;
	double potential;
	modforce = LJ_F(r);
	potential = LJ_P(r);
	this->LJ_potential += potential;
	for (int n = 0; n < 3; ++n) {
		this->power[n] += (this->coord[n] - mol_prob->coord[n] - shift[n]) / r * modforce;
	}
}

void Atom::Veloc_shift(Space* space) {
	for (int n = 0; n < 3; ++n) {
		this->vel[n] += (this->power[n] * (space->dt / this->m));
	}
	this->Kinetic_en = Kin_En(space, this);
	space->Kin_En += this->Kinetic_en;
	space->Pot_En += this->LJ_potential;
}

int Atom::Change_Cell(Space* space, int i0, int j0, int k0, int l) {
	int i = static_cast<int>(this->coord[0] / space->cellsize);
	int j = static_cast<int>(this->coord[1] / space->cellsize);
	int k = static_cast<int>(this->coord[2] / space->cellsize);
	if ((i != i0) || (j != j0) || (k != k0)) {
		if ((i >= 0) && (j >= 0) && (k >= 0) && (i < space->amount_cells) && (j < space->amount_cells) && (k < space->amount_cells)) {
			////Добавляем в новую ячейку
			//Atom* atom;
			//atom = space->cells[i][j][k].atoms[space->cells[i][j][k].amount_atoms];
			//space->cells[i][j][k].atoms[space->cells[i][j][k].amount_atoms] = this;
			//space->cells[i][j][k].amount_atoms++;
			////Убираем из старой ячейки
			//space->cells[i0][j0][k0].amount_atoms--;
			//space->cells[i0][j0][k0].atoms[l] = space->cells[i0][j0][k0].atoms[space->cells[i0][j0][k0].amount_atoms];
			//space->cells[i0][j0][k0].atoms[space->cells[i0][j0][k0].amount_atoms] = atom;

			swap(space->cells[i0][j0][k0].atoms[l], space->cells[i][j][k].atoms[space->cells[i][j][k].amount_atoms]);
			space->cells[i][j][k].amount_atoms++;
			space->cells[i0][j0][k0].amount_atoms--;
			swap(space->cells[i0][j0][k0].atoms[l], space->cells[i0][j0][k0].atoms[space->cells[i0][j0][k0].amount_atoms]);
		}
		else {
			assert((i >= 0) && (j >= 0) && (k >= 0) && (i < space->amount_cells) && (j < space->amount_cells) && (k < space->amount_cells));
			printf("Going beyond the calculation area\n");
		}
		return 1;
	}
	else {
		return 0;
	}
}

void SetNullMacro(Space* space) {
	space->Kin_En = 0;
	space->Pot_En = 0;
	space->Energy = 0;
	for (int i = 0; i < space->amount_cells; ++i) {
		for (int j = 0; j < space->amount_cells; ++j) {
			for (int k = 0; k < space->amount_cells; ++k) {
				for (int l = 0; l < space->cells[i][j][k].amount_atoms; ++l) {
					for (int n = 0; n < 3; ++n) {
						space->cells[i][j][k].atoms[l]->power[n] = 0;
					}
					space->cells[i][j][k].atoms[l]->LJ_potential = 0;
					space->cells[i][j][k].atoms[l]->Kinetic_en = 0;
				}
			}
		}
	}
}

double KX_P(const double& r) {
	return k_N * pow((r0 - r), 2);
}

double KX_F(const double& r) {
	return - 2 * k_N * (r - r0);
}

double LJ_F(const double& r) {
	return A / pow(r, 13) - B / pow(r, 7);
}

double LJ_P(const double& r) {
	return AA / pow(r, 12) - BB / pow(r, 6);
}

double Kin_En(Space* space, Atom* mol) {
	return mol->m * (pow(mol->vel[0], 2) + pow(mol->vel[1], 2) + pow(mol->vel[2], 2)) / 2;
}

void Get_energy(std::ofstream& fouten, const Space& space, std::ofstream& Kout, std::ofstream& Pout, std::ofstream& Fout, std::ofstream& Rout, std::ofstream& Tout, const int& step) {
	fouten << "Full energy = " << space.Energy << "\t" << "Kinetic energy = " 
		<< space.Kin_En << "\t" << "Potential energy = " << space.Pot_En << "\n";
	Kout << step << "\t" <<  space.Kin_En << endl;
	Pout << step << "\t" << space.Pot_En << endl;
	Fout << step << "\t" << space.Energy << endl;
	Tout << step << "\t" << space.T << endl;
	Rout << step << "\t" << Rms_Vel(space) << endl;
}

void Print_atoms(const Space& space) {
	ofstream fout;
	fout.open("C:\\Drive\\Code\\MolDynamics\\md_v1.4\\Output\\get_atoms.txt");
	Atom* atom;
	for (int i = 0; i < space.amount_cells; ++i) {
		for (int j = 0; j < space.amount_cells; ++j) {
			for (int k = 0; k < space.amount_cells; ++k) {
				for (int l = 0; l < space.cells[i][j][k].amount_atoms; ++l) {
					atom = space.cells[i][j][k].atoms[l];
					fout << atom << "   " << space.cells[i][j][k].amount_atoms << '\t' << i << " " << j << " " << k << '\t'
						<< atom->coord[0] << " " << atom->coord[1] << " " << atom->coord[2] << '\t'
						<< atom->vel[0] << " " << atom->vel[1] << " " << atom->vel[2] << '\n';
				}
			}
		}
	}
	fout.close();
}

void Print_last(const Space& space) {
	ofstream fout;
	fout.open("C:\\Drive\\Code\\MolDynamics\\md_v1.4\\Output\\get_last.txt");
	Atom* atom;
	for (int i = 0; i < space.amount_cells; ++i) {
		for (int j = 0; j < space.amount_cells; ++j) {
			for (int k = 0; k < space.amount_cells; ++k) {
				for (int l = 0; l < space.cells[i][j][k].amount_atoms; ++l) {
					atom = space.cells[i][j][k].atoms[l];
					fout << atom->coord[0] << " " << atom->coord[1] << " " << atom->coord[2] << '\t'
						<< atom->vel[0] << " " << atom->vel[1] << " " << atom->vel[2] << '\n';
				}
			}
		}
	}
	fout.close();
}

void Print_mol(const Space& space) {
	ofstream fout;
	fout.open("C:\\Drive\\Code\\MolDynamics\\md_v1.4\\Output\\get_atoms2.txt");
	Atom* atom;
	for (int i = 0; i < space.total_mol_N2; ++i) {
		for (int at = 0; at < 2; ++at) {
			atom = space.N_2[i].atom[at];
			fout << atom << "   " << atom->coord[0] << " " << atom->coord[1] << " " << atom->coord[2] << '\t'
				<< atom->vel[0] << " " << atom->vel[1] << " " << atom->vel[2] << '\n';
		}
	}
	fout.close();
}

void Print_config(const Space& space, const int& amount_steps) {
	ofstream fout;
	fout.open("C:\\Drive\\Code\\MolDynamics\\md_v1.4\\Output\\get_config.txt");
	fout << space.total_mol_N2 * 2 + space.total_at_Pt << '\n' << space.T << '\n' << space.p << '\n' << space.dt << '\n' << amount_steps << '\n'
		<< space.spacesize << '\n' << space.cellsize << '\n' << space.amount_cells << '\n';
	fout.close();
}

double Temp(Space& space) {
	double mV2 = 0;
	double a = 0;
	double res = 0;
	for (int i = 0; i < space.amount_cells; ++i) {
		for (int j = 0; j < space.amount_cells; ++j) {
			for (int k = 0; k < space.amount_cells; ++k) {
				for (int l = 0; l < space.cells[i][j][k].amount_atoms; ++l) {
					Atom* at = space.cells[i][j][k].atoms[l];
					a = at->m * (pow(at->vel[0], 2) + pow(at->vel[1], 2) + pow(at->vel[2], 2));
					mV2 += a;
				}
			}
		}
	}
	res = mV2 / ((space.total_at_Pt + space.total_mol_N2) * (3 * K_b));
	space.T = res;
	return res;
}

double Rms_Vel(const Space& space) {
	double res = 0;
	res = sqrt(3 * K_b * space.T / space.m_Pt);
	return res;
}

int VTK_num = 1000;
//VTK
int WriteVTK(Space* space)
{
	int tmn = 0;
	for (int i = 0; i < (*space).amount_cells; i++)
		for (int j = 0; j < (*space).amount_cells; j++)
			for (int k = 0; k < (*space).amount_cells; k++)
			{
				tmn += (*space).cells[i][j][k].amount_atoms;
			}
	char fname[100];
	sprintf(fname, "C:\\Drive\\Code\\MolDynamics\\md_v1.4\\Output\\VTK\\./state_%010d.vtk", VTK_num);
	FILE* f;
	f = fopen(fname, "w");
	fprintf(f, "# vtk DataFile Version 2.0\nMolecules states\nASCII\nDATASET POLYDATA\nPOINTS %d float\n", tmn);
	for (int i = 0; i < (*space).amount_cells; i++)
		for (int j = 0; j < (*space).amount_cells; j++)
			for (int k = 0; k < (*space).amount_cells; k++)
			{
				for (int n = 0; n < (*space).cells[i][j][k].amount_atoms; n++)
				{
					fprintf(f, "%lf %lf %lf\n", (*space).cells[i][j][k].atoms[n]->coord[0],
						(*space).cells[i][j][k].atoms[n]->coord[1], (*space).cells[i][j][k].atoms[n]->coord[2]);
				}
			}
	fprintf(f, "POINT_DATA %d\nSCALARS MoleculeType float 1\nLOOKUP_TABLE default\n", tmn);
	for (int i = 0; i < (*space).amount_cells; i++)
		for (int j = 0; j < (*space).amount_cells; j++)
			for (int k = 0; k < (*space).amount_cells; k++)
			{
				for (int n = 0; n < (*space).cells[i][j][k].amount_atoms; n++)
				{
					fprintf(f, "%d ", 1);
				}
			}
	VTK_num++;
	fclose(f);
	return 0;
}

