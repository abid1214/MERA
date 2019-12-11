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


void Abidtest(double W, int L, unsigned seed, double epsilon)
{
    bool uniform = false;
	int pD = 2;
	int bD = 5;
	double WJ = 0;
	double Wh = W;

    //unsigned seed = 31;
	std::mt19937 generator((unsigned int)time(NULL));
	std::uniform_real_distribution<double> distribution(0.0,1.0);
    char* opts = new char [L];
    for(int j = 0; j < L; ++j)
        opts[L-j-1] = distribution(generator) > epsilon ? 'L' : 'H';
        //opts[L-j-1] = 'L';


	sdMERA sdM;
	sdM.setMaxBondDim(100);
	sdM.setWJ(WJ);
	sdM.setWh(Wh);
	sdM.setRandomSeed(seed);
	sdM.setMaxSearchLength(4,4);
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

	cout<<"st\top\ted\tphy"<<endl;
    for(int i = 0; i < L; i++)
    {
        sdM.unitaryDecimateMPO(opts[i]);
    }

    MPS psi;
    sdM.buildMPS(psi);
    psi.RC();
    double E = psiHphi(psi,h,psi);
    double SE = psiHphi(psi,hS,psi);
    cout<<"MPS Info: "<<E<<" "<<SE-E*E<<" ";
    psi.EE();

    /*
    if(sdM.good_RG_flow)
        cout<<1;
    else
        cout<<0;
    
    cout<<endl<<endl;
    */
}


int main (int argc, char const *argv[])
{
	double      W  = argc>1 ? atof(argv[1]) : 0;
    int         L  = argc>2 ? atof(argv[2]) : 32;
    unsigned seed  = argc>3 ? atof(argv[3]) : (unsigned int)time(NULL);
    double n          = argc>4 ? atof(argv[4]) : 0;
    double N          = argc>5 ? atof(argv[5]) : 1;
    double epsilon =  n/N;
    cout<<"epsilon = "<<epsilon<<endl;
    Abidtest(W, L, seed, epsilon);
	return 0;
}
