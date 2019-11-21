#ifndef My_UTILITY_FUNCTIONS_H
#define My_UTILITY_FUNCTIONS_H

#include <cmath>
#include <cassert>
#include <algorithm>
#include <Eigen/Dense>
#include <string>

#include "mps.h"
#include "mpo.h"
#include "iterCompress.h"
#include "EDtoMPO.h"
// #include "ezh5.h"

typedef Eigen::MatrixXd Mxd;


// B = A1 cross A2
void KroneckerProd(Mxd& A1, Mxd& A2, Mxd& B);

void effH(MPO& H, int site, int L, Mxd& Heff);

void applyMPO(MPO& A, MPO& B, int site, char direc, char op);

void applyMPO(MPO& A, MPS& psi, int site, char op);

void applyGates(Mxd& G, MPO& H, int site, int L, int max_BD);

void applyGates(Mxd& G, MPO& H, int site, int L);

double var(MPO& H);

// Find Permutation -- to speed up Wegner Diagonalization processes
bool findPermutation(Mxd& A);

// Returns true if the matrix contains any NaN entries in it
bool findNaN(Mxd &H);

 // Symmetrizes a matrix 
void mkSymmetric(Mxd& h1);

// Find the max element of the matrix -- including or excluding the diagonal
double findMax(Mxd H, bool withDiagonal);

// Check whether a matrix is close to identity
bool is_indentity(Mxd& U);

// void writeToFile(std::string fn, MPS& psi);
//
// void readFromFile(std::string fn, MPS& psi);
//
// void writeToFile(std::string fn, MPO& H);
//
// void readFromFile(std::string fn, MPO& H);


#endif
