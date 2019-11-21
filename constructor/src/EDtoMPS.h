#ifndef MY_ED_TO_MPS_FUNCTION_H
#define MY_ED_TO_MPS_FUNCTION_H

#include <cmath>
#include <cassert>
#include <algorithm>
#include <Eigen/Dense>
#include "mps.h"
#include "lapack_wrapper.h"

typedef Eigen::MatrixXd Mxd;

void EDtoMPS(double * vec, int Len, int pD, MPS& psi);


#endif