#ifndef My_SIMPLE_WEGNER_FLOW_H
#define My_SIMPLE_WEGNER_FLOW_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <Eigen/Dense>

#include "lapack_wrapper.h"
#include "utility.h"

typedef Eigen::MatrixXd Mxd;

class SimpleWegner
{
public:
	Mxd A;
    Mxd H;
    Mxd U;
    double tau;
	double tol;
	int iter;
	
	SimpleWegner(){};
	~SimpleWegner(){};
	
	void setH(Mxd& tH);
	
	void setTau(double t);
	void setTol(double t);
	void setIter(int t);
	void diag(Mxd& tU);
};


#endif