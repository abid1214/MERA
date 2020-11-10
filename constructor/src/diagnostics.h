#ifndef DIAGNOSTICS
#define DIAGNOSTICS

#include "diagnostics.h"
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include "utility.h"
#include "mpo.h"
#include <cmath>
#include "simpleWegnerFlow.h"
#include <time.h>
#include <stdlib.h>
#include "wegnerFlow.h"
#include "mps.h"

using namespace std;

MPS load_mps(string data_dir, double W, int L, int l, unsigned seed, double epsilon);

vector<double> setDisorderConfig(int seed, int L, double WJ, double Wh);

double get_energy(MPS psi, double W, int seed);

double get_variance(MPS psi, double W, int seed);

double bipartite_spin_fluctuation(MPS psi);

double participation_entropy(MPS psi, int q);

void run_diagnostics(string data_dir, double W, int L, int l, unsigned seed, double epsilon);

#endif

