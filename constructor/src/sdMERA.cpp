#ifndef My_sdMERA_CLASS
#define My_sdMERA_CLASS

#include "sdMERA.h"
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
        for(int i = 0; i < L; ++i)
        {
            dJ.push_back(2*WJ*(drand48()-0.5));
            dh.push_back(2*Wh*(drand48()-0.5));
        }
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
	L  = _L;
	pD = _pD;
	bD = _bD;
	setDisorderConfig(uniform);
	
	H.clearMPO();
	H.setMPO(L, pD, bD, 0);
	H.buildHeisenberg(&dJ[0], &dh[0]);
    Mxd A;
	effH(H,0,L,A);
    cout<<A<<endl;

	
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
	//std::cout<<"Unitary gate is identity = "<<id_g<<std::endl;
	//if(id_g)
	//	std::cout<<site_op<<"\t"<<site_op<<"\t"<<site_op<<"\t"<<phy<<std::endl;
	//else
	//	std::cout<<site_st<<"\t"<<site_op<<"\t"<<site_ed<<"\t"<<phy<<std::endl;
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
	{
		psi.M[0][i] = H.M[0][i*pD+(uG[s-1].phy)];
	}
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
		if(psi.bD>phi.bD) iterCompress(true, 4, 1e-14, max_bD, psi, true);
	}
}

void sdMERA::findMaxGap()
{
	gaps.setZero(max_search_L,L);
	max_gap = 0;
	int slen = stLen;
	if(slen>L) slen = L;
		
	Mxd A;
	while(slen<=max_search_L&&slen<=H.Len)
	{
		for(int i = 0; i < L+1-slen; ++i)
		{
			effH(H,i,slen,A);
			
			Eigen::SelfAdjointEigenSolver<Mxd> es(A);
			if (es.info() != Eigen::Success) abort();
			gaps(slen-1,i) = es.eigenvalues()(A.rows()-1)-es.eigenvalues()(0);
			
			if(max_gap < gaps(slen-1,i))
			{
				max_gap_site = i;
				max_gap_L    = slen;
				max_gap      = gaps(slen-1,i);
			}
		}
		++slen;
	}
	ENSC = max_gap;
	//std::cout<<"Max Gap = "<<max_gap<<", at Site "<<max_gap_site<<", length "<<max_gap_L<<std::endl;
}

double sdMERA::getTauBits(Mxd& A, int& tpL, int& pos)
{
	MPO tH;
	EDtoMPO(A, tpL, pD, tH);
	std::vector<int> idx_set;
	idx_set.push_back(pos);
	double tp = tH.trace(idx_set);
	//std::cout<<"TBC = "<<tp<<std::endl;
	return tp;
}

void sdMERA::getTauBits(Mxd& A, int& tpL)
{
	MPO tH;
	EDtoMPO(A, tpL, pD, tH);
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
			//std::cout<<"TBC"<<slen<<" = "<<tpH.trace(idx_set)<<std::endl;
			// generate the next index set
			for(int i = 0; i < slen; ++i)
			{
				if (i==slen-1 && idx[i]==tpL-1)
					not_exhausted = false;
				else if (i==slen-1 && idx[i]<tpL-1)
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

void sdMERA::unitaryDecimateMPO(char opt)
{
	assert(opt=='L'||opt=='H');
	
    //cout<<"finding max gap"<<std::endl;
	findMaxGap();
	
	int pos = max_gap_site;
	int tpL = max_gap_L;
	
	Mxd A,U,D;
    //cout<<"calculating effective hamiltonian"<<std::endl;
	effH(H,max_gap_site,max_gap_L,A);
    cout<<A<<endl;
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
	int idx = 0;
	int phy = 0;
	double sum=9999;
	
	// Looking for the bond to decimate
	for(int i = 0; i < max_gap_L; ++i)
	{
		std::vector<double> tp(pD);
		for(int j = 0; j < std::pow(pD,max_gap_L); ++j)
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
	
    int opt_phy = opt=='L' ? phy : 1 - phy;
	addContractedSite(H.M_IDs[max_gap_site],H.M_IDs[idx+max_gap_site],H.M_IDs[max_gap_L+max_gap_site-1],opt_phy,max_gap_L,idx,U);
	//std::cout<<"Fixing site "<<H.M_IDs[idx+max_gap_site]<<" to "<<opt_phy<<std::endl;
	
    //cout<<"getting tau bits"<<endl;
	D = U.transpose() * A * U;
	double leading_energy = getTauBits(D, max_gap_L, idx);
	getTauBits(D, max_gap_L);
	
	if(L>1)
	{

        if(L<6)
        {
            Eigen::SelfAdjointEigenSolver<Mxd> es(A);
            if (es.info() != Eigen::Success) abort();
            std::cout<<L<<"\t";
		    Mxd evls = es.eigenvalues();
            for(int k = 0; k < std::pow(pD,L); k++)
                std::cout<<evls(k)<<"\t";
            std::cout<<std::endl;
        }


        //cout<<"applying gates"<<endl;
		applyGates(U, H, max_gap_site, max_gap_L, max_bD);
		int idx_new = idx;
		if(max_gap_site-1>=0)
		{
			pos -= 1;
			tpL += 1;
			idx_new += 1;
		}
		if(max_gap_site+max_gap_L<L)
		{
			tpL += 1;
		}

        //cout<<"calculating Heff"<<endl;
		effH(H,pos,tpL,A);
		double g_factor = -99999;
		int idx_i=0, idx_j=0;
		for(int i = 0; i < A.rows(); ++i)
		{
			for(int j = i+1; j < A.cols(); ++j)
			{
                double tp = std::log10( std::abs( A(i,j)/leading_energy ) );
                if(A(i,j)!=0 && A(i,i)!=A(j,j) && g_factor<tp) 
                {
                    g_factor = tp;
                    idx_i = i;
                    idx_j = j;
                }
			}
		}
		//std::cout<<"Max g factor = "<<g_factor<<" "<<A(idx_i,idx_i)<<" "<<A(idx_j,idx_j)<<std::endl;
		if(g_factor>0) good_RG_flow = false;
		g_factors.push_back(g_factor);

        //cout<<"getting eigens"<<endl;
		Eigen::SelfAdjointEigenSolver<Mxd> es(A);
		Mxd evls = es.eigenvalues();
		double mean_gap_ratio = 0;
		for(int i = 0; i < evls.size()-2; ++i)
		  {
		    double d1 = std::abs(evls(i)-evls(i+1));
		    double d2 = std::abs(evls(i+1)-evls(i+2));
		    mean_gap_ratio += std::min(d1,d2)/std::max(d1,d2);
		  }
		//if(evls.size()-2>0) std::cout<<"MGapRatio = "<<mean_gap_ratio/(evls.size()-2)<<std::endl;
		
        //std::cout<<"decimating Hamiltonian"<<endl;
		MPO HH(L-1,pD,bD,0);
		if(opt=='L')
			HH.decimateCopy(H, max_gap_site+idx, phy);
		else
			HH.decimateCopy(H, max_gap_site+idx, 1-phy);
		H.clearMPO();
		H.copyMPO(HH);
	}else
	{
        Eigen::SelfAdjointEigenSolver<Mxd> es(A);
        if (es.info() != Eigen::Success) abort();
        std::cout<<L<<"\t";
		Mxd evls = es.eigenvalues();
        for(int k = 0; k < std::pow(pD,L); k++)
            std::cout<<evls(k)<<"\t";
        std::cout<<std::endl;
	}
	
	L--;
}

void sdMERA::renormalize()
{
	int idx = 0;
	std::cout<<"st\top\ted\tphy"<<std::endl;
	while(L>0)
	{
		unitaryDecimateMPO(opts[idx]);
		/*(L>=2)
		{
			MPO A;
			A.copyMPO(H);
			std::cout<<"calculating EE"<<std::endl;
			std::cout<<"MPO EE = ";
			A.EE();
			std::cout<<" "<<L<<std::endl;
		}
        */
		idx++;
	}
    MPO A;
    A.copyMPO(H);
    std::cout<<"calculating EE"<<std::endl;
    std::cout<<"MPO EE = ";
    A.EE();
    std::cout<<" "<<std::endl;
}

#endif
