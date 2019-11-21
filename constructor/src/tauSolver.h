#ifndef My_tauSolver_CLASS_H
#define My_tauSolver_CLASS_H

#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <string>
#include "mpo.h"
#include "utility.h"
#include "lapack_wrapper.h"
#include "simpleWegnerFlow.h"
#include "wegnerFlow.h"
#include "ezh5.h"

typedef Eigen::MatrixXd Mxd;

class tauSolver
{
public:
	// MPO variables
	MPO S;
	MPO H;
	int L;
	int pD;
	int bD;
	int max_bD;
	
	// Disorder variables
	std::vector<double> dJ;
	std::vector<double> dh;
	double WJ;
	double Wh;
	int randSeed;
	
	tauSolver(){};
	~tauSolver(){};
	
	void setMaxBondDim(int t);
	
	void setWJ(double t);
	void setWh(double t);
	void setRandomSeed(int t);
	void setDisorderConfig();
	
	void setInitialMPO(int _L, int _pD, int _bD);
	
	void getTauBits(MPO& A, int tpL);
	double getTauBits(MPO& A, int tpL, int& pos);
	
	void diag(int maxIter, int subLen);
};


#endif