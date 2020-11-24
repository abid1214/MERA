#ifndef My_tauSolver_CLASS
#define My_tauSolver_CLASS

#include "tauSolver.h"
typedef Eigen::MatrixXd Mxd;

void tauSolver::setMaxBondDim(int t)
{
	max_bD = t;
}

void tauSolver::setWJ(double t)
{
	WJ = t;
}

void tauSolver::setWh(double t)
{
	Wh = t;
}

void tauSolver::setRandomSeed(int t)
{
	randSeed = t;
}

void tauSolver::setDisorderConfig()
{
	srand48(131*randSeed+17);
	dJ.clear();
	dh.clear();
	for(int i = 0; i < L; ++i)
	{
		dJ.push_back(2*WJ*(drand48()-0.5));
		dh.push_back(2*Wh*(drand48()-0.5));
	}
	for(int i = 0; i < L; ++i) std::cout<<dh[i]<<" ";
	std::cout<<std::endl;
	//exit(0);
	// for(int i = 0; i < L; ++i)
	// {
	// 	std::cout<<1+dJ[i]<<" ";
	// }
	// std::cout<<std::endl;
	// for(int i = 0; i < L; ++i)
	// {
	// 	std::cout<<dh[i]/2<<" ";
	// }
	// std::cout<<std::endl;
}

void tauSolver::setInitialMPO(int _L, int _pD, int _bD)
{
	L  = _L;
	pD = _pD;
	bD = _bD;
	setDisorderConfig();
	H.clearMPO();
	H.setMPO(L, pD, bD, 0);
	H.buildHeisenberg(&dJ[0],&dh[0]);
	// H.buildRTIC(&dJ[0],&dh[0]);
	
	S.setMPO(L, pD, 1, 0);
	S.buildSz(0);
}

double tauSolver::getTauBits(MPO& A, int tpL, int& pos)
{
	// get tau bits coefficients
	MPO tH;
	tH.copyMPO(A);
	std::vector<int> idx_set;
	idx_set.push_back(pos);
	double tp = tH.trace(idx_set);
	std::cout<<"TBC = "<<tp<<std::endl;
	return tp;
}

void tauSolver::getTauBits(MPO& A, int tpL)
{
	MPO tH;
	tH.copyMPO(A);
	// get tau bits coefficients
	int slen = 1;
	while(slen<=tpL)
	{
		int* idx = new int [slen];
		for(int i = 0; i < slen; ++i)
		{
			idx[i] = i;
		}
		bool not_exhausted = true;
		do
		{
			MPO tpH;
			tpH.copyMPO(tH);
			std::vector<int> idx_set;
			for(int i = 0; i < slen; ++i)
			{
				idx_set.push_back(idx[i]);
			}
			std::cout<<"TBC "<<slen<<" ";
			for(int i = 0; i < slen; ++i) std::cout<<idx[i]<<" ";
			std::cout<<tpH.trace(idx_set)<<std::endl;
			// generate the next index set
			for(int i = 0; i < slen; ++i)
			{
				if (i==slen-1 && idx[i]==tH.Len-1)
					not_exhausted = false;
				else if (i==slen-1 && idx[i]<tH.Len-1)
					{++idx[i]; break;}
				else if (idx[i] < idx[i+1]-1)
					{++idx[i]; break;}
				else if (idx[i] == idx[i+1]-1 && i==0)
					idx[i] = 0;
				else if (idx[i] == idx[i+1]-1 && i>0)
					idx[i] = idx[i-1]+1;
			}
			
		}while(not_exhausted);
		delete [] idx;
		++slen;
	}
}

void tauSolver::diag(int maxIter, int subLen)
{
	assert(subLen>0&&subLen<L);
	double maxOffDiag = 0;
	int step = 0;
	std::vector<int> st_pos;
	std::vector<int> ugLen;
	std::vector<Mxd> ug;
	
	std::cout<<"Variance = "<<var(H)<<std::endl;
	
	std::string fn = "TauBit_L"+std::to_string(L)+"_RS"+std::to_string(randSeed)+"_full.h5";
	ezh5::File fh5W(fn, H5F_ACC_TRUNC);
	fh5W["Sys"]["L"]        = L;
	fh5W["Sys"]["pD"]       = pD;
	fh5W["Sys"]["bD"]       = bD;
	fh5W["Sys"]["max_bD"]   = max_bD;
	fh5W["Sys"]["Wh"]       = Wh;
	fh5W["Sys"]["WJ"]       = WJ;
	fh5W["Sys"]["randSeed"] = randSeed;
	fh5W["Sys"]["maxIter"]  = maxIter;
	fh5W["Sys"]["subLen"]   = subLen;
	
	auto gMPO   = fh5W["MPO"];
	gMPO["Len"] = H.Len;
	gMPO["pD"]  = H.pD;
	
	for(int i = 0; i < H.Len; ++i)
	{
		for(int j = 0; j < H.pD*H.pD; ++j)
		{
			std::string site_name = "Site" + std::to_string(i);
			std::string phy_name  = "M" + std::to_string(j);
			gMPO[site_name][phy_name] = H.M[i][j];
		}
	}
	
	auto gGate = fh5W["uGate"];
	int idx_gate = 0;
	
	///////////////
	// Test
	subLen = H.Len;
	Mxd A,U;
	effH(H,0,subLen,A);
	WegnerDiagonalize WD;
	WD.setH(A);
	WD.setTol(1e-6);
	WD.setTau(1e-4);
	WD.diag(U);

	std::string gate_name = "G"+std::to_string(idx_gate);
	gGate[gate_name] = U;
	++idx_gate;

	applyGates(U, H, 0, subLen, max_bD);
	U.transposeInPlace();
	st_pos.push_back(0);
	ugLen.push_back(subLen);
	ug.push_back(U);

	std::cout<<"Variance = "<<var(H)<<std::endl;
	S.setMPO(L, pD, 1, 0);
	S.buildSz(0);
	for(int i = ug.size()-1; i>= 0; --i)
	  {
	    applyGates(ug[i], S, st_pos[i], ugLen[i], max_bD);
	  }
	getTauBits(S, 1);
	///////////////

	// subLen = 2;
	// while(subLen<=5)
	// {	
	//   step = 0;
	//   std::cout<<subLen<<std::endl;
	  
	//   while(step<maxIter)
	//     {
	// 	int st = 0;
	// 	int ed = subLen - 1;
		
	// 	// moving to the right
	// 	while(ed<L)
	// 	  {
	// 	    Mxd A,U;
	// 	    effH(H,st,subLen,A);
	// 	    WegnerDiagonalize WD;
	// 	    WD.setH(A);
	// 	    WD.setTol(1e-6);
	// 	    WD.setTau(1e-4);
	// 	    WD.diag(U);
		    
	// 	    std::string gate_name = "G"+std::to_string(idx_gate);
	// 	    gGate[gate_name] = U;
	// 	    ++idx_gate;
		    
	// 	    applyGates(U, H, st, subLen, max_bD);
	// 	    // applyGates(U, S, st, subLen, max_bD);
	// 	    U.transposeInPlace();
	// 	    st_pos.push_back(st);
	// 	    ugLen.push_back(subLen);
	// 	    ug.push_back(U);
		    
		    
	// 	    ++st;
	// 	    ++ed;
	// 	  }
	// 	--st;
	// 	--ed;
		
	// 	std::cout<<"Variance = "<<var(H)<<std::endl;		
	// 	S.setMPO(L, pD, 1, 0);
	// 	S.buildSz(0);
	// 	for(int i = ug.size()-1; i>= 0; --i)
	// 	  {
	// 	    applyGates(ug[i], S, st_pos[i], ugLen[i], max_bD);
	// 	  }
	// 	getTauBits(S, 1);
		
	// 	while(st>=0)
	// 	  {
	// 	    Mxd A,U;
	// 	    effH(H,st,subLen,A);
	// 	    WegnerDiagonalize WD;
	// 	    WD.setH(A);
	// 	    WD.setTol(1e-6);
	// 	    WD.setTau(1e-4);
	// 	    WD.diag(U);
	// 	    applyGates(U, H, st, subLen, max_bD);
	// 	    // applyGates(U, S, st, subLen, max_bD);
	// 	    U.transposeInPlace();
	// 	    st_pos.push_back(st);
	// 	    ugLen.push_back(subLen);
	// 	    ug.push_back(U);
		    
	// 	    std::string gate_name = "G"+std::to_string(idx_gate);
	// 	    gGate[gate_name] = U;
	// 	    ++idx_gate;
		    
	// 	    --st;
	// 	    --ed;
	// 	  }
		
	// 	std::cout<<"Variance = "<<var(H)<<std::endl;
	// 	S.setMPO(L, pD, 1, 0);
	// 	S.buildSz(0);
	// 	for(int i = ug.size()-1; i>= 0; --i)
	// 	  {
	// 	    applyGates(ug[i], S, st_pos[i], ugLen[i], max_bD);
	// 	  }
	// 	getTauBits(S, 1);
		
	// 	++step;
	//     }
	//   ++subLen;
	// }
}

#endif
