#ifndef My_qBlockSolver_CLASS_H
#define My_qBlockSolver_CLASS_H

#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include "mpo.h"
#include "utility.h"
#include "lapack_wrapper.h"
#include "simpleWegnerFlow.h"
#include "wegnerFlow.h"

typedef Eigen::MatrixXd Mxd;

class qBlockSolver
{
public:
	// MPO variables
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
	
	qBlockSolver(){};
	~qBlockSolver(){};
	
	void setMaxBondDim(int t);
	
	void setWJ(double t);
	void setWh(double t);
	void setRandomSeed(int t);
	void setDisorderConfig();
	
	void setInitialMPO(int _L, int _pD, int _bD);
	
	void diag(int maxIter, int subLen);
};


#endif