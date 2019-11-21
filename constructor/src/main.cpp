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

void trim(Mxd& m, double tol)
{
	for(int i = 0; i < m.size(); ++i)
	{
		if(abs(m(i))<tol) m(i) = 0;
	}
}

void testsdMERA(double W)
{
	typedef std::chrono::high_resolution_clock myclock;
	myclock::time_point beginning = myclock::now();
	cout.precision(10);
	int L  = 32;
	int pD = 2;
	int bD = 5;
	double WJ = 0;
	double Wh = W;
	int rs = 31;
	
	sdMERA sdM;
	sdM.setMaxBondDim(40);
	sdM.setWJ(WJ);
	sdM.setWh(Wh);
	sdM.setRandomSeed(rs);
	sdM.setMaxSearchLength(1,8);
	sdM.setInitialMPO(L,pD,bD);

	myclock::duration d = myclock::now() - beginning;
	unsigned seed = rs+1;
	std::mt19937 generator (seed);
	std::uniform_real_distribution<double> distribution(0.0,1.0);
	std::cout<<"Random seed = "<<seed<<std::endl;
	
	char* opts = new char [L];
	MPS psi;
	for(int i = 100; i < 101; ++i)
	{
		for(int j = 0; j < L; ++j)
		{
		  if(distribution(generator)>0.5)
				opts[L-j-1] = 'L';
			else
				opts[L-j-1] = 'H';
		}
		sdM.setOpts(opts);
		sdM.renormalize();
		sdM.setInitialMPO(L,pD,bD);
		sdM.buildMPS(psi);
		psi.RC();
		cout<<"Norm of MPS = "<<psiphi(psi,psi)<<endl;
		double E = psiHphi(psi,sdM.H,psi);
		double SE = psiHphi(psi,sdM.HS,psi);
		cout<<"MPS Info: "<<E<<" "<<SE-E*E<<" ";
		psi.EE();
		
		if(sdM.good_RG_flow)
			cout<<1;
		else
			cout<<0;
		
		cout<<endl<<endl;
	}
}

void testqBlockSolver()
{
	int L  = 10;
	int pD = 2;
	int bD = 5;
	double WJ = 0;
	double Wh = 10;
	int rs = 1;
	
	qBlockSolver qbS;
	qbS.setMaxBondDim(40);
	qbS.setWJ(WJ);
	qbS.setWh(Wh);
	qbS.setRandomSeed(rs);
	qbS.setInitialMPO(L,pD,bD);
	
	qbS.diag(4,5);
}

void testTauBit(int _RS)
{
	cout.precision(15);
	int L  = 10;
	int pD = 2;
	int bD = 5;
	double WJ = 0;
	double Wh = 10;
	int rs = _RS;
	
	tauSolver qbS;
	qbS.setMaxBondDim(60);
	qbS.setWJ(WJ);
	qbS.setWh(Wh);
	qbS.setRandomSeed(rs);
	qbS.setInitialMPO(L,pD,bD);
	
	qbS.diag(1,10);
}


void Abidtest(double W)
{
    int L  = 32;
	int pD = 2;
	int bD = 5;
	double WJ = 0;
	double Wh = W;

	unsigned seed = 32;
	std::mt19937 generator (seed);
	std::uniform_real_distribution<double> distribution(0.0,1.0);
    char* opts = new char [L];
    for(int j = 0; j < L; ++j)
        opts[L-j-1] = distribution(generator)>0.5 ? 'L' : 'H';


	sdMERA sdM;
	sdM.setMaxBondDim(40);
	sdM.setWJ(WJ);
	sdM.setWh(Wh);
	sdM.setRandomSeed(seed - 1);
	sdM.setMaxSearchLength(1,4);
	sdM.setInitialMPO(L,pD,bD);
    sdM.setOpts(opts);

    sdM.renormalize();

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
	double W = argc>1 ? atof(argv[1]) : 0;
    //cout<<"W = "<<W<<endl;
    Abidtest(W);
	return 0;
}
