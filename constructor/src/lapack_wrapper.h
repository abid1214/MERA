#ifndef My_LAPACK_WRAPPERS_H
#define My_LAPACK_WRAPPERS_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include <vector>
typedef Eigen::MatrixXd Mxd;

#ifndef Use_MKL
extern "C" void dgeqrf_(int*, int*, double*, int*, double*, double*, int*, int*);
extern "C" void dorgqr_(int*, int*, int*, double*, int*, double*, double*, int*, int*);
extern "C" void dgesvd_(char*, char*, int*, int*, double*, int*, double*, double*, int*, double*, int*, double*, int*, int*);
extern "C" void dsytrs_(char*,int*,int*,double*,int*,int*,double*,int*,int*); // broken?
extern "C" void dsysv_(char*,int*,int*,double*,int*,int*,double*,int*,double*,int*,int*); // the only safe one to use
extern "C" void dposv_(char*,int*,int*,double*,int*,double*,int*,int*);
extern "C" void dgesv_(int*,int*,double*,int*,int*,double*,int*,int*); // the only safe one to use
extern "C" void dsposv_(char*,int*,int*,double*,int*,double*,int*,double*,int*,double*,float*,int*,int*);
extern "C" void dsysvx_(char*,char*,int*,int*,double*,int*,double*,int*,int*,double*,int*,double*,int*,double*,double*,double*,double*,int*,int*,int*);
extern "C" void dposvx_(char*,char*,int*,int*,double*,int*,double*,int*,char*,double*,double*,int*,double*,int*,double*,double*,double*,double*,int*,int*);
extern "C" void dsyev_(char*,char*,int*,double*,int*,double*,double*,int*,int*);
#endif

#ifdef Use_MKL
#include "mkl.h"
#endif

void diag(double* M, double* evals, int nn);

void rQR(const Mxd& MT, Mxd& Q);

void rSVD(const Mxd& M, int ds, double * sv, Mxd& UM, Mxd& VTM, char direction);

void rSVD(const Mxd& M, int ds, double * sv, Mxd& UM, Mxd& VTM, int truncD, double& svd_error);

void rSVD(const Mxd& M, std::vector<double>& sv, Mxd& UM, Mxd& VTM, double& svd_tol, bool dry_run=true);

void linearSolver(int ord, Mxd& A, Mxd& vec);

void indef_linearSolver(int ord, Mxd& A, Mxd& vec);

int def_linearSolver(int ord, Mxd& A, Mxd& vec);

void geindef_linearSolver(int ord, Mxd& A, Mxd& vec);

int sdef_linearSolver(int ord, Mxd& A, Mxd& vec);

void sindef_linearSolver(int ord, Mxd& A, Mxd& vec);

void sxdef_linearSolver(int ord, Mxd& A, Mxd& vec);

#endif
