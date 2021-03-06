// no-quantum-number mps class
#ifndef My_MPS_CLASS_H
#define My_MPS_CLASS_H

#include <vector>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <cmath>
#include <Eigen/Dense>

typedef Eigen::MatrixXd Mxd;
using namespace std;

class MPS{
public:
	int Len;
	int pD;
	int bD;
	Mxd ** M;
	int* Dim;
	
	bool if_init;
	
	MPS ();
	MPS (int l, int pd, int bd);
	MPS (const MPS& other);
	~ MPS ();
	
	void setMPS(int l, int pd, int bd);
	void setNorm(double nm);
	void setZero();
	void setRand();
	void setZero(int newD);
	void setRand(int newD);
	void clear();
	void addSite(int site);
	void setDim();
	
	void RC();
	void LC();
	void moveLeft(int site);
	void moveRight(int site);

	void copyMPS(const MPS& other);
	void addMPS(double coeff, const MPS& other);
	
	void writeMPS(std::string filename, int preci);
	void readMPS(std::string filename);
	
    vector<double> EE(bool verbose);
	double norm();
    Mxd evaluateMPS(int N);
    Mxd partial_trace(int pos, int l);
};

#endif
