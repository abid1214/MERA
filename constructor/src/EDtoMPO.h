#ifndef MY_ED_TO_MPO_FUNCTION_H
#define MY_ED_TO_MPO_FUNCTION_H

#include <cmath>
#include <cassert>
#include <algorithm>
#include <Eigen/Dense>
#include "mpo.h"
#include "lapack_wrapper.h"

typedef Eigen::MatrixXd Mxd;

void EDtoMPO(Mxd& A, int Len, int pD, MPO& H);


#endif