#ifndef My_OBSERVABLES_H
#define My_OBSERVABLES_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include "mps.h"
#include "mpo.h"

// using Eigen::MatrixXd;

typedef Eigen::MatrixXd Mxd;

double psiHphi(const MPS& Psi, MPO& HH, const MPS& Phi);

double psiphi(const MPS& psi, const MPS& phi);

#endif