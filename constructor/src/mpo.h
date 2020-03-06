#ifndef My_MPO_CLASS_H
#define My_MPO_CLASS_H

#include <string>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
#include <Eigen/Dense>

#include "lapack_wrapper.h"
#include "mps.h"

typedef Eigen::MatrixXd Mxd;

struct SpM
{
	std::vector<int> r;
	std::vector<int> c;
	std::vector<double> v;
};

class MPO{
public:
	int Len;
	int pD;
	int bD;
	Mxd ** M;
	SpM ** H;
	int *  Dim;
	int *  M_IDs;
	
	double norm;
	double J1;
	double J2;
	double ens;
	
	bool if_init;
	MPO ();
	MPO (int l, int pd, int bd, double EnergyShifted);
	~ MPO ();
	void setMPO(int l, int pd, int bd, double EnergyShifted);
	void clearMPO();
	
	void copyMPO(const MPO& other);
	void decimateCopy(const MPO& other, int site, int phy);
	double trace(std::vector<int>& pos);
	double sumAll();
	void keepDiag();
	void diagToMPS(MPS& psi);
	
	void buildHeisenberg(double* dJ, double* dh);
    void buildIsing(double* J, double* h, double* Jp);
	void buildRTIC(double* dJ, double* dh);
	void buildSPTC(double* dJ, double* dh);
	void buildSz(int site);
	void addMPO(double coeff, const MPO& other);
	
	void buildSpMPO();
	
	void setZero();
	void setZero(std::string s);
	void setZero(int nbD);
	void setRand();
	void setRand(int nbD);
	void setIdentity();
	void setIdentity(int nbD);
	void square();
	
	void LC();               // Left  Canonicalization
	void RC();               // right Canonicalization
	void moveRight(int site);
	void moveLeft(int site);
	void compressL(int nbD); // Compress from left
	void examineBond(double tol, bool dry_run=true);
	void EE();
};

#endif
