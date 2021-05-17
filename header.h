//#pragma once

#include <iostream>
#include <fstream>
#include <algorithm>
#include <ctime>
#include <cmath>
#include <chrono>
#include <assert.h>
#include <vector>
#pragma warning (disable : 4996)
////////////////////Memory Leaks/////////////////////////
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
/////////////////////////////////////////////////////////

using namespace std;
using namespace std::chrono;

const int MAX_MOLS = 100000;
const int MAX_ATOMS_IN_CELL = 2000;
const int MAX_CELLS = 30;
const double K_b = 0.00138;						//Boltzmann's constant

/**********   H   **********/
// const double A = 12 * pow(2.64, 12) * 1.53 * pow(10, -2);					//A = 12 * pow(SIGMA, 12) * EPS;
// const double B = 12 * pow(2.64, 6) * 1.53 * pow(10, -2);						//B = 12 * pow(SIGMA, 6) * EPS;
// const double AA = pow(2.64, 12) * 1.53 * pow(10, -2);					//AA = pow(SIGMA, 12) * EPS;
// const double BB = 2 * pow(2.64, 6) * 1.53 * pow(10, -2);					//BB = 2 * pow(SIGMA, 6) * EPS;
// const double A = 61625.553346;					//A = 12 * pow(SIGMA, 12) * EPS;
// const double B = 89.788794;						//B = 12 * pow(SIGMA, 6) * EPS;
// const double AA = 5135.462778;					//AA = pow(SIGMA, 12) * EPS;
// const double BB = 14.964799;					//BB = 2 * pow(SIGMA, 6) * EPS;

/**********   Pt   **********/
const double eps_div_k = 3771.5;
const double eps = eps_div_k * K_b;
const double sigma = 2.523;
const double A = 12 * 4 * pow(sigma, 12) * eps;					
const double B = 6 * 4 * pow(sigma, 6) * eps;						
const double AA = 4 * pow(sigma, 12) * eps;					
const double BB = 4 * pow(sigma, 6) * eps;					

//Поменять)))))
const double k_N = 326;                           //Stretch constant
const double r0 = 1.4;                            //Equilibrium length


// g++ main.cpp md.cpp md.h header.h -o main
// g++ main.cpp md.cpp md.h header.h -O3 -o main

// 1.5301 * 2.64^6 * 10^-2

//A = 12 * pow(SIGMA, 12) * EPS; 
// 12 * 2.64^12 * 1.53 * 10^-2