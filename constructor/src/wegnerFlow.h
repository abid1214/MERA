#ifndef WegnerFlow_SOLVER_H
#define WegnerFlow_SOLVER_H


#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>

#include <Eigen/Dense>

#include "lapack_wrapper.h"
#include "utility.h"

typedef Eigen::MatrixXd Mxd;


class WegnerDiagonalize
{
public:
	WegnerDiagonalize(){};
	~WegnerDiagonalize(){};
	
    // integration variables
	Mxd dH;
    Mxd H;
    Mxd U;
    double tau;           // integration time
    double deltaTau;      // integration step
    int steps;            // total number of steps taken in diagonalization
  
    //scratch space data
    double tauRun;
    Mxd HRun;
    Mxd URun;
    Mxd DeltaH;
    Mxd DeltaU;
    Mxd DeltaH5;  // these are the high precision values used for stepping
    Mxd DeltaU5;
    Mxd DeltaH4;  // these are the low precision values used for tunning step size
    Mxd DeltaU4;

    //parameters
    double tolleranceRK;
    double tolleranceElem;

    void setH(Mxd &myH);
    void setTolRK(double t_tolleranceRK);
    void setTol(double t_tolleranceElem);
    void setTau(double t_deltaTau);
    // computes F' using HRun and URun
    void computeFP();
    // performs a RK45 step (as defined in wikepedia)
    void RKstep();
    // runs the diagonalization using RK45 steps with adjustable step sizes
    bool diag(Mxd& _U);

};


#endif