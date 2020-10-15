#ifndef BLOCK_FINDER_H
#define BLOCK_FINDER_H

#include <Eigen/Dense>
#include "mpo.h"
typedef Eigen::MatrixXd Mxd;
using namespace std;

int tree_block(int L, int L0, int l);

int random_block(int L, int l);

int findMaxGap(MPO& H, int l);

int firstBlock();

std::vector<int> int_to_bin(int n, int l, int pD);

int bin_to_int(vector<int> sv, int l, int pD);

Mxd U_seperator(Mxd& U, int pD, int l);

std::vector<double> get_sv(Mxd& Umat);

double get_E(std::vector<double>& sv);

int findMinEntanglement(MPO& H, int l);

#endif
