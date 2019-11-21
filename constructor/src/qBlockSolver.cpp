#ifndef My_qBlockSolver_CLASS
#define My_qBlockSolver_CLASS

#include "qBlockSolver.h"
typedef Eigen::MatrixXd Mxd;

void qBlockSolver::setMaxBondDim(int t)
{
	max_bD = t;
}

void qBlockSolver::setWJ(double t)
{
	WJ = t;
}

void qBlockSolver::setWh(double t)
{
	Wh = t;
}

void qBlockSolver::setRandomSeed(int t)
{
	randSeed = t;
}

void qBlockSolver::setDisorderConfig()
{
	srand48(137*randSeed);
	dJ.clear();
	dh.clear();
	for(int i = 0; i < L; ++i)
	{
		dJ.push_back(2*WJ*(drand48()-0.5));
		dh.push_back(2*Wh*(drand48()-0.5));
	}
	for(int i = 0; i < L; ++i)
	{
		std::cout<<1+dJ[i]<<" ";
	}
	std::cout<<std::endl;
	for(int i = 0; i < L; ++i)
	{
		std::cout<<dh[i]/2<<" ";
	}
	std::cout<<std::endl;
}

void qBlockSolver::setInitialMPO(int _L, int _pD, int _bD)
{
	L  = _L;
	pD = _pD;
	bD = _bD;
	setDisorderConfig();
	H.clearMPO();
	H.setMPO(L, pD, bD, 0);
	H.buildHeisenberg(&dJ[0],&dh[0]);
	// H.buildRTIC(&dJ[0],&dh[0]);
}

void qBlockSolver::diag(int maxIter, int subLen)
{
	assert(subLen>0&&subLen<L);
	double maxOffDiag = 0;
	int step = 0;
	
	std::cout<<"Variance = "<<var(H)<<std::endl;
	
	while(step<maxIter)
	{
		int st = 0;
		int ed = subLen - 1;
		
		// moving to the right
		while(ed<L)
		{
			Mxd A,U;
			effH(H,st,subLen,A);
			WegnerDiagonalize WD;
			WD.setH(A);
			WD.setTol(1e-6);
			WD.setTau(1e-4);
			WD.diag(U);
			applyGates(U, H, st, subLen, max_bD);
			
			std::cout<<"Variance = "<<var(H)<<std::endl;
			
			++st;
			++ed;
		}
		--st;
		--ed;
		while(st>=0)
		{
			Mxd A,U;
			effH(H,st,subLen,A);
			WegnerDiagonalize WD;
			WD.setH(A);
			WD.setTol(1e-6);
			WD.setTau(1e-4);
			WD.diag(U);
			applyGates(U, H, st, subLen, max_bD);
			
			std::cout<<"Variance = "<<var(H)<<std::endl;
			
			--st;
			--ed;
		}
		++step;
	}
	
}



#endif
