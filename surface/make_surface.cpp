#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <fstream>
#include <stdio.h>
#include <ctime>

using namespace std;

const double PI = 3.141592653589793;

double RandNorm() {
    double res = 0;
    while (res < 0.0000000001 || res > 0.9999999999) {
        res = ((double)rand()) / (double)(RAND_MAX);
    }
    return res;
}

double RandG() {
    return sqrt(-log(RandNorm())) * sin(2 * PI * RandNorm());
}

int main() {
    double temp = 0;
    srand(time(nullptr));
    cout << "enter temp: \n";
    cin >> temp;
    double a;
    cout << "enter size: \n";
    cin >> a;
    int I, J, K;
    cout << "enter amount atoms x-dir, y-dir, z-dir: \n";
    cin >> I >> J >> K;
    vector<vector<double>> coord;
    vector<int> type;
    double shift = 0.25;
    double shift2 = 3.5; // 13.72
    for (int k = 0; k < K / 2 + 1; ++k) {
        for (int i = 0; i < I; ++i) {
            for (int j = 0; j < J; ++j) {
                cout << i << "\t" << j << "\t" << k << endl;
                cout << i + 1.0 / 2 << "\t" << j + 1.0 / 2 << "\t" << k << endl;
                coord.push_back({(i + shift) * a, (j + shift) * a, (k + shift2) * a});
                coord.push_back({(i + 1.0 / 2 + shift) * a, (j + 1.0 / 2 + shift) * a, (k + shift2) * a});
                type.push_back(1);
                type.push_back(1);
            }
        }
        if (k != K / 2) {
            for (int i = 0; i < I; ++i) {
                for (int j = 0; j < J; ++j) {
                    cout << i << "\t" << j + 1.0 / 2 << "\t" << k + 1.0 / 2 << endl;
                    cout << i + 1.0 / 2 << "\t" << j << "\t" << k + 1.0 / 2 << endl;
                    coord.push_back({(i + shift) * a, (j + 1.0 / 2 + shift) * a, (k + 1.0 / 2 + shift2) * a});
                    coord.push_back({(i + 1.0 / 2 + shift) * a, (j + shift) * a, (k + 1.0 / 2 + shift2) * a});
                    type.push_back(2);
                    type.push_back(2);
                }
            }
        }
    }
    int VTK_num = 0;
    int tmn = coord.size();
	char fname[100];
	sprintf(fname, "./state_%010d.vtk", VTK_num);
	FILE* f;
	f = fopen(fname, "w");
    fprintf(f, "# vtk DataFile Version 2.0\nMolecules states\nASCII\nDATASET POLYDATA\nPOINTS %d float\n", tmn);
    for (auto x : coord) {
        fprintf(f, "%lf %lf %lf\n", x[0], x[1], x[2]);
    }
    fprintf(f, "POINT_DATA %d\nSCALARS MoleculeType float 1\nLOOKUP_TABLE default\n", tmn);
    for (int i = 0; i < coord.size(); ++i) {
        fprintf(f, "%d ", type[i]);
    }
    VTK_num++;
	fclose(f);

    ofstream fout;
    fout.open("input.txt");
    
    double m_Pt = 3.24e-5;
    const double K_b = 0.00138;
    double betta = sqrt(m_Pt / K_b / temp / 2.0);
    for (auto x : coord) {
        for (int i = 0; i < 3; ++i) {
            fout << x[i] << " ";
        }
        for (int j = 0; j < 3; ++j) {
            double new_vel = RandG() / betta;
            fout << new_vel << " ";
        }
        fout << endl;
    }
    fout.close();

    return 0;
}
