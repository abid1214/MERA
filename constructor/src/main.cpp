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


void Abidtest(double W, int L, int l, unsigned seed, double epsilon)
{
    bool uniform = false;
	int pD = 2;
	int bD = 5;
	double WJ = 0;
	double Wh = W;

	std::mt19937 generator(seed);
	std::uniform_real_distribution<double> distribution(0.0,1.0);
    char* opts = new char [L];
    for(int j = 0; j < L; ++j)
    {
        opts[L-j-1] = distribution(generator) > epsilon ? 'L' : 'H';
    }


	sdMERA sdM;
	sdM.setMaxBondDim(100);
	sdM.setWJ(WJ);
	sdM.setWh(Wh);
	sdM.setRandomSeed(seed);
	sdM.setMaxSearchLength(1,l);
	sdM.setInitialMPO(L,pD,bD, uniform);
    sdM.setOpts(opts);

    MPO h;
    h.clearMPO();
    h.setMPO(L, pD, bD, 0);
    h.buildHeisenberg(&(sdM.dJ[0]), &(sdM.dh[0]));

    Mxd Hmat;
    effH(h, 0, L, Hmat);
    cout<<Hmat<<endl;

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
        sdM.unitaryDecimateMPO(opts[i]);
    }

    MPS psi;
    sdM.buildMPS(psi);
    psi.RC();

    double E = psiHphi(psi,h,psi);
    double SE = psiHphi(psi,hS,psi);
    double sz = psiHphi(psi, Sz, psi);

    psi.EE(true);
    cout<<endl;

    cout<<"MPS Info: "<<E<<" "<<SE-E*E<<" ";
    psi.EE(false);
    cout<<" "<<sz<<endl;

    char fname [50];
    sprintf(fname, "mps_W_%2.4f_L_%d_l_%d_d_%d_e_%0.2f.txt", W, L, l, seed, epsilon);
    psi.writeMPS(fname, 10);
    

    /*
    if(sdM.good_RG_flow)
        cout<<1;
    else
        cout<<0;
    
    cout<<endl<<endl;
    */
}

MPS load_mps(string data_dir, double W, int L, int l, unsigned seed, double epsilon)
{
    char fname [150];
    sprintf(fname, "%smps_W_%2.4f_L_%d_l_%d_d_%d_e_%0.2f.txt", data_dir.c_str(), W, L, l, seed, epsilon);
    MPS psi;
    psi.readMPS(fname);
    return psi;
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
    //cout<<"epsilon = "<<epsilon<<endl;
    //Abidtest(W, L, l, seed, epsilon);

    string data_dir1 = "/home/aakhan3/scratch/MERA/data/flip_data/MPS_orig/";
    string data_dir2 = "/home/aakhan3/scratch/MERA/data/flip_data/MPS_last_bit_flip/";
    string data_dir3 = "/home/aakhan3/scratch/MERA/data/flip_data/MPS_last_lbit_flip/";

    MPS psi = load_mps(data_dir1, W, L, l, seed, epsilon);
    MPS phi = load_mps(data_dir3, W, L, l, seed, epsilon);
    cout<<psiphi(psi, phi)<<endl;

    char fname [50];
    sprintf(fname, "overlapl_W_%2.4f_L_%d_l_%d_d_%d_e_%0.2f.txt", W, L, l, seed, epsilon);
    ofstream fout;
    fout.precision(8);
    fout.open(fname);

    for(int site=0; site<L; site++)
        fout<<MPS_partial_overlap(psi, phi, site, 1)<<endl;
	return 0;
}
