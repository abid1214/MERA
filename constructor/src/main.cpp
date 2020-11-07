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


void Abidtest(double W, int L, int l, unsigned seed, double epsilon)
{
    bool uniform = false;
	int pD = 2;
	int bD = 5;
	double WJ = 0;
	double Wh = W;

    /*
	std::mt19937 generator(seed);
	std::uniform_real_distribution<double> distribution(0.0,1.0);
    char* opts = new char [L];
    for(int j = 0; j < L; ++j)
        opts[L-j-1] = distribution(generator) > epsilon ? 'L' : 'H';

    //opts[L-1] = opts[L-1] == 'H' ? 'L' : 'H';
    //opts[L/2] = opts[L/2] == 'H' ? 'L' : 'H';
    */

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

    char fname [50];
    sprintf(fname, "mps_data/mps_W_%2.4f_L_%d_l_%d_d_%d_e_%0.2f.txt", W, L, l, seed, epsilon);
    psi.writeMPS(fname, 10);
}

MPS load_mps(string data_dir, double W, int L, int l, unsigned seed, double epsilon)
{
    char fname [150];
    sprintf(fname, "%smps_W_%2.4f_L_%d_l_%d_d_%d_e_%0.2f.txt", data_dir.c_str(), W, L, l, seed, epsilon);
    MPS psi;
    psi.readMPS(fname);
    return psi;
}

double mean(double a[], int n)
{
    double avg = 0;
    for(int i=0; i<n; i++)
        avg += a[i];
    return avg/n;
}

void disorder_average_overlap(double W, int L, int l, double epsilon)
{

    int num_dis = 100;
    double overlap_list [L][num_dis];
    string data_dir1 = "/home/abid/programs/MERA/constructor/data/emera/";
    string data_dir2 = "/home/abid/programs/MERA/constructor/data/fmera/";

    for(int seed=0; seed<num_dis; seed++)
    {
        MPS psi = load_mps(data_dir1, W, L, l, seed, epsilon);
        MPS phi = load_mps(data_dir2, W, L, l, seed, epsilon);
        cout<<psiphi(psi, phi)<<endl;

        for(int site=0; site<L; site++)
            overlap_list[site][seed] = MPS_partial_overlap(psi, phi, site, 1);
    }

    char fname [50];
    sprintf(fname, "dis_avg_overlap_W_%2.4f_L_%d_l_%d_e_%0.2f.txt", W, L, l,epsilon);
    ofstream fout;
    fout.precision(8);
    fout.open(fname);
    for(int site=0; site<L; site++)
    {
        sort(overlap_list[site], overlap_list[site]+num_dis);
        double median = overlap_list[site][num_dis/2];
        double avg = mean(overlap_list[site], num_dis);
        fout<<avg<<endl;
    }
}

void overlap_hist(double W, int L, int l, double epsilon, int site)
{

    int num_dis = 100;
    string data_dir1 = "/home/abid/programs/MERA/constructor/data/orig/";
    string data_dir2 = "/home/abid/programs/MERA/constructor/data/flip/";

    char fname [50];
    sprintf(fname, "hist_overlap_W_%2.4f_L_%d_l_%d_e_%0.2f_s_%d.txt", W, L, l,epsilon, site);
    ofstream fout;
    fout.precision(8);
    fout.open(fname);

    for(int seed=0; seed<num_dis; seed++)
    {
        MPS psi = load_mps(data_dir1, W, L, l, seed, epsilon);
        MPS phi = load_mps(data_dir2, W, L, l, seed, epsilon);
        fout<<MPS_partial_overlap(psi, phi, site, 1)<<endl;
    }

    
}

void single_overlap(double W, int L, int l, unsigned seed, double epsilon)
{
    string data_dir1 = "/home/abid/programs/MERA/constructor/data/orig/";
    string data_dir2 = "/home/abid/programs/MERA/constructor/data/flip/";

    MPS psi = load_mps(data_dir1, W, L, l, seed, epsilon);
    MPS phi = load_mps(data_dir2, W, L, l, seed, epsilon);
    cout<<psiphi(psi, phi)<<endl;

    char fname [50];
    sprintf(fname, "overlap_W_%2.4f_L_%d_l_%d_d_%d_e_%0.2f.txt", W, L, l, seed, epsilon);
    ofstream fout;
    fout.precision(8);
    fout.open(fname);

    for(int site=0; site<L; site++)
        fout<<MPS_partial_overlap(psi, phi, site, 1)<<endl;
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
    Abidtest(W, L, l, seed, epsilon);
    //disorder_average_overlap(W, L, l, epsilon);
    //for(int site=0; site<L; site++)
    //overlap_hist(W, L, l, epsilon, seed);
    //single_overlap( W, L, l, seed, epsilon);
	return 0;
}
