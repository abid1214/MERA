#ifndef My_sdMERA_CLASS_H
#define My_sdMERA_CLASS_H

#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include "mpo.h"
#include "utility.h"
#include "lapack_wrapper.h"
#include "simpleWegnerFlow.h"
#include "wegnerFlow.h"

typedef Eigen::MatrixXd Mxd;

class uniGate
{
public:
	Mxd M;
	int L;
	int site;
	
	// if it is decimated
	int idx;
	int phy;
};

class sdMERA
{
public:
	// MPO variables
	MPO H;
	MPO HS;
    int L0;
	int L;
	int pD;
	int bD;
	int max_bD;
	
	// Disorder variables
	std::vector<double> dJ;
	std::vector<double> dh;
	double WJ;
	double Wh;
	int randSeed;
	
	// Gap searching variables
	double ENSC;
	int stLen;
	int max_search_L;
	double max_gap;
	int max_gap_site;
	int max_gap_L;
	Mxd gaps;
	
	// quality check of the RG flow;
	std::vector<double> g_factors;
	bool good_RG_flow;
	
	// Which sector to keep
	std::vector<double> opts;
	
	// Build the MPS from the tensor network
	std::vector<uniGate> uG;
	
	void addContractedSite(int site_st, int site_op, int site_ed, int phy, int Len, int idx, Mxd& A);
	
	void buildMPS(MPS& psi);
	
	sdMERA(){};
	~sdMERA(){};
	
	void setMaxBondDim(int t);
	
	void setWJ(double t);
	void setWh(double t);
	void setRandomSeed(int t);
	void setDisorderConfig(bool uniform);
	
	void setMaxSearchLength(int st, int maxLen);
	
	void setOpts(char* _opts);
	
	void setInitialMPO(int _L, int _pD, int _bD, bool uniform);
	
	void unitaryDecimateMPO(char* opts);
    
	
	void renormalize();
};


#endif
