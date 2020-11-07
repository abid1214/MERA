#ifndef My_sdMERA_CLASS
#define My_sdMERA_CLASS

#include "sdMERA.h"
#include "observables.h"
#include "block_finder.h"

typedef Eigen::MatrixXd Mxd;
using namespace std;

void sdMERA::setMaxBondDim(int t)
{
	max_bD = t;
	ENSC = 9999;
}

void sdMERA::setWJ(double t)
{
	WJ = t;
}

void sdMERA::setWh(double t)
{
	Wh = t;
}

void sdMERA::setRandomSeed(int t)
{
	randSeed = t;
}

void sdMERA::setDisorderConfig(bool uniform)
{
    double h;
	srand48(137*randSeed);
	dJ.clear();
	dh.clear();
    if(uniform) {
        for(int i = 0; i < L; ++i)
        {
            dJ.push_back(WJ);
            dh.push_back(Wh);
        }
    }
    else {
        cout<<"W_LIST: ";
        for(int i = 0; i < L; ++i)
        {
            dJ.push_back(2*WJ*(drand48()-0.5));
            h = 2*Wh*(drand48()-0.5);
            cout << h << " "; 
            dh.push_back(h);
        }
        cout<<endl;
    }
	
}

void sdMERA::setMaxSearchLength(int st, int maxLen)
{
	stLen = st;
	max_search_L = maxLen;
}

void sdMERA::setOpts(char* _opts)
{
	opts.clear();
	for(int i = 0; i < L; ++i)
	{
		opts.push_back(_opts[i]);
	}
	uG.clear();
	
	good_RG_flow = true;
	g_factors.clear();
}

void sdMERA::setInitialMPO(int _L, int _pD, int _bD, bool uniform)
{
    L0 = _L;
	L  = _L;
	pD = _pD;
	bD = _bD;
	setDisorderConfig(uniform);
	
	H.clearMPO();
	H.setMPO(L, pD, bD, 0);
	H.buildHeisenberg(&dJ[0], &dh[0]);

    //Mxd A;
    //effH(H, 0, L, A);
    //std::cout<<A<<std::endl;

	
	HS.clearMPO();
	HS.setMPO(L, pD, bD, 0);
	HS.buildHeisenberg(&dJ[0], &dh[0]);
	HS.square();
}

void sdMERA::addContractedSite(int site_st, int site_op, int site_ed, int phy, int Len, int idx, Mxd& A)
{
	uniGate tp;
	tp.L    = Len;
	tp.site = site_op;
	tp.phy  = phy;
	tp.idx  = idx;
	tp.M    = A;
	uG.push_back(tp);
	bool id_g = is_indentity(A);
	if(id_g)
    {
	    std::cout<<"WARN: Unitary gate is identity"<<std::endl;
		std::cout<<"SITE: "<<site_op<<"\t"<<site_op<<"\t"<<site_op<<"\t"<<phy<<std::endl;
    }
	else
		std::cout<<"SITE: "<<site_st<<"\t"<<site_op<<"\t"<<site_ed<<"\t"<<phy<<std::endl;
}

void sdMERA::buildMPS(MPS& psi)
{
	int s = uG.size();
	assert(s>0);
	psi.clear();
	psi.setMPS(1,H.pD,1);
	std::vector<int> positions;
	MPO H;
	EDtoMPO(uG[s-1].M, uG[s-1].L, pD, H);
	for(int i = 0; i < pD; ++i)
		psi.M[0][i] = H.M[0][i*pD+(uG[s-1].phy)];
	psi.setDim();
	positions.push_back(uG[s-1].site);
	//std::cout<<"Unitary Gate at the two site level: "<<std::endl<<uG[s-2].M<<std::endl;
	for(int i = s-2; i >= 0; --i)
	{
		int site  = 0;
		for(int j = 0; j < positions.size(); ++j)
		{
			if(positions[j]<uG[i].site) site++;
		}
		psi.addSite(site);
		positions.push_back(uG[i].site);
		EDtoMPO(uG[i].M, uG[i].L, pD, H);
		for(int j = 0; j < pD; ++j)
		{
			H.M[uG[i].idx][j*pD+1-uG[i].phy] = H.M[uG[i].idx][j*pD+uG[i].phy];
		}
		applyMPO(H, psi, site-uG[i].idx, 'I');
		MPS phi(psi.Len,psi.pD,std::min(psi.bD,max_bD));
		phi.setZero();
		phi.setDim();
		if(psi.bD>phi.bD) iterCompress(true, 4, 1e-14, max_bD, psi, false);
	}
}


void sdMERA::unitaryDecimateMPO(char* opts)
{
	
    int U_start, U_L;
    U_L     = (stLen <= L) ? stLen : L;
    //cout<<"finding max gap"<<std::endl;
	U_start = findMinEntanglement(H, stLen);
	//U_start = random_block(L, U_L);
	//U_start = findMaxGap(H, stLen);
	//U_start = firstBlock();
    //cout<<"found block "<<U_start<<" to "<<U_start + U_L-1<<endl;
	
	        
	Mxd A,U,D;
    //cout<<"calculating effective hamiltonian"<<std::endl;
    effH(H,U_start,U_L,A);
	bool perm_found = false;
	double Init_Tol = 1e-5;
	double Init_Tau = 1e-4;

    //cout<<"performing Wegner flow"<<endl;
	while(!perm_found)
	{
		WegnerDiagonalize WD;
		WD.setH(A);
		WD.setTol(Init_Tol);
		WD.setTau(Init_Tau);
		perm_found = WD.diag(U);
		Init_Tol /= 1.1;
		Init_Tau *= 1.1;
	}
	D = (U.transpose() * A * U).diagonal();
	
	// Which index to fix?
	// Fix to 0 or 1?
    int start = (U_start == 0) ? 0 : 1;
    int end   = (U_L == L) ? U_L-1: U_L-2;

	int idx = start;
	int phy = 0;
	double sum=9999;
	
	for(int i = start; i <= end; i++)
	{
		std::vector<double> tp(pD);
		for(int j = 0; j < std::pow(pD,U_L); ++j)
		{
			tp[int(j/std::pow(pD,i))%pD] += D(j);
		}
		if(sum > *std::min_element(tp.begin(),tp.end()))
		{
			sum = *std::min_element(tp.begin(),tp.end());
			idx = i;
			phy = (std::min_element(tp.begin(),tp.end())-tp.begin());
		}
	}
	
    char opt = opts[H.M_IDs[idx+U_start]];
    int opt_phy = opt=='L' ? phy : 1 - phy;

	addContractedSite(H.M_IDs[U_start],H.M_IDs[idx+U_start],H.M_IDs[U_L+U_start-1],opt_phy,U_L,idx,U);
	
	if(L>1)
	{

        //cout<<"applying gates"<<endl;
		applyGates(U, H, U_start, U_L, max_bD);
        		
        //std::cout<<"decimating Hamiltonian"<<endl;
		MPO HH(L-1,pD,bD,0);
		HH.decimateCopy(H, U_start+idx, opt_phy);
		H.clearMPO();
		H.copyMPO(HH);
    }
	L--;
}


#endif
