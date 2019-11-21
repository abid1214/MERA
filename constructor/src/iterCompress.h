#ifndef My_Iterative_Compression_H
#define My_Iterative_Compression_H

#include <cstdlib>
#include <Eigen/Dense>
#include "mps.h"
#include "mpo.h"

typedef Eigen::MatrixXd Mxd;


void buildR(MPS& H, MPS& Hc, Mxd* CR);

void updateSite(MPS& H, MPS& Hc, Mxd* CL, Mxd* CR, int& site, char& direc, double& dt);

void updateEnvr(MPS& H, MPS& Hc, Mxd* CL, Mxd* CR, int& site, char& direc, double& dt);

void iterCompress(bool is_Random, int maxIter, double tol, int nbD, MPS& H);

void iterCompress(bool is_Random, int maxIter, double tol, int nbD, MPS& H, bool print_info);

////////////////////////////////////////////////////////////

void buildR(MPO& H, MPO& Hc, Mxd* CR);

void updateSite(MPO& H, MPO& Hc, Mxd* CL, Mxd* CR, int& site, char& direc, double& dt);

void updateEnvr(MPO& H, MPO& Hc, Mxd* CL, Mxd* CR, int& site, char& direc, double& dt);

void iterCompress(bool is_Random, int maxIter, double tol, int nbD, MPO& H);

void iterCompress(bool is_Random, int maxIter, double tol, int nbD, MPO& H, bool print_info);

#endif