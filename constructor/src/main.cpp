#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <algorithm>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <time.h>
#include <omp.h>
#include <cstdint>
#include <string>
#include <ctime>
#include <random>
#include <chrono>

#ifdef Use_MKL
#define EIGEN_USE_MKL_ALL
#endif

#include <Eigen/Dense>
#include <limits>

using namespace Eigen;
using namespace std;

typedef Eigen::MatrixXd Mxd;

#include "observables.h"
#include "utility.h"
#include "EDtoMPS.h"
#include "EDtoMPO.h"
#include "simpleWegnerFlow.h"
#include "sdMERA.h"
#include "qBlockSolver.h"
#include "tauSolver.h"
#include "tests.h"
#include <algorithm>
#include <math.h>
#include <random>
#include "diagnostics.h"


void Abidtest(double W, int L, int l, unsigned seed, double epsilon)
{
    bool uniform = false;
	int pD = 2;
	int bD = 5;
	double WJ = 0;
	double Wh = W;

    int t = (int) round(epsilon*L);
    char* opts = new char [L];
    for(int j = 0; j < L; j++)
        opts[L-1-j] = (j < t) ? 'H' : 'L';
    shuffle(&opts[0], &opts[L], std::default_random_engine(seed));
    shuffle(&opts[0], &opts[L], std::default_random_engine(seed));
    shuffle(&opts[0], &opts[L], std::default_random_engine(seed));
    shuffle(&opts[0], &opts[L], std::default_random_engine(seed));
    shuffle(&opts[0], &opts[L], std::default_random_engine(seed));

    for(int j = 0; j < L; ++j)
        cout<<opts[j];
    cout<<endl;

	sdMERA sdM;
	sdM.setMaxBondDim(100);
	sdM.setWJ(WJ);
	sdM.setWh(Wh);
	sdM.setRandomSeed(seed);
	sdM.setMaxSearchLength(l,l);
	sdM.setInitialMPO(L,pD,bD, uniform);
    sdM.setOpts(opts);

    MPO h;
    h.clearMPO();
    h.setMPO(L, pD, bD, 0);
    h.buildHeisenberg(&(sdM.dJ[0]), &(sdM.dh[0]));


    MPO hS;
    hS.clearMPO();
    hS.setMPO(L, pD, bD, 0);
    hS.buildHeisenberg(&(sdM.dJ[0]), &(sdM.dh[0]));
    hS.square();

    MPO Sz;
    Sz.setMPO(L, pD, bD, 0);
    Sz.buildSz(L/2);

	cout<<"st\top\ted\tphy"<<endl;
    for(int i = 0; i < L; i++)
    {
        sdM.unitaryDecimateMPO(opts);
    }

    MPS psi;
    sdM.buildMPS(psi);
    psi.RC();

    double E = psiHphi(psi,h,psi);
    double SE = psiHphi(psi,hS,psi);

    cout<<"MPS Info: "<<E<<" "<<SE-E*E<<" ";
    psi.EE(false);
    cout<<endl;
    psi.EE(true);
    cout<<endl;

    char fname [50];
    sprintf(fname, "mps_data/mps_W_%2.4f_L_%d_l_%d_d_%d_e_%0.2f.txt", W, L, l, seed, epsilon);
    psi.writeMPS(fname, 10);
    
}

int main (int argc, char const *argv[])
{
	double       W = argc>1 ? atof(argv[1]) : 0;
    int          L = argc>2 ? atof(argv[2]) : 32;
    int          l = argc>3 ? atof(argv[3]) : 4;
    unsigned  seed = argc>4 ? atof(argv[4]) : 0;
    double       n = argc>5 ? atof(argv[5]) : 0;
    double       N = argc>6 ? atof(argv[6]) : 1;

    double epsilon =  n/N;
    //run_tests();
    //Abidtest(W, L, l, seed, epsilon);
    //disorder_average_overlap(W, L, l, epsilon);
    //for(int site=0; site<L; site++)
    //overlap_hist(W, L, l, epsilon, seed);
    //single_overlap( W, L, l, seed, epsilon);
    run_diagnostics("mps_data/", W, L, l, seed, epsilon);
	return 0;
}
