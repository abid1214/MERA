#ifndef My_OBSERVABLES
#define My_OBSERVABLES

#include "observables.h"

// using Eigen::MatrixXd;

typedef Eigen::MatrixXd Mxd;

double psiHphi (const MPS& Psi, MPO& HH, const MPS& Phi)
{
	int tid;
	int r1, r2, c1, c2;
	double EnergyExpec = 0;
	int Len = Psi.Len;
	int phy = Psi.pD;
	int DimW = HH.bD;
	HH.buildSpMPO();
	////////////////////////////////////////////
	Mxd ** CRM = new Mxd * [Len];
	for(int i = 0; i < Len; i++)
	{
		CRM[i] = new Mxd [DimW];
	}
	Mxd ** TM1 = new Mxd * [phy];
	Mxd ** TM2 = new Mxd * [phy];
	for(int i = 0; i < phy; i++)
	{
		TM1[i] = new Mxd [DimW];
		TM2[i] = new Mxd [DimW];
	}
	////////////////////////////////////////////////
	// Last site //
	// Building the FR //
	r1=Psi.M[Len-1][0].rows();
	c1=Psi.M[Len-1][0].cols();
	r2=Phi.M[Len-1][0].rows();
	c2=Phi.M[Len-1][0].cols();
	for(int i = 0; i < DimW; i++)
	{
		CRM[Len-1][i].setZero(r1,r2);
	}
	for(int i = 0; i < phy; i++)
	{
		for(int k = 0; k < phy; k++)
		{
			for(int j = 0; j < HH.H[Len-1][phy*i+k].r.size(); j++)
			{
				CRM[Len-1][HH.H[Len-1][phy*i+k].r[j]].noalias() += HH.H[Len-1][phy*i+k].v[j] * (Psi.M[Len-1][i]*Phi.M[Len-1][k].transpose());
			}
		}
	}
	// Last-but-1 to the second site //
	for(int i = Len-2; i >= 0; i--) {
		// std::cout<<"Site "<<i<<std::endl;
		r1=Psi.M[i][0].rows();
		c1=Psi.M[i][0].cols();
		r2=Phi.M[i][0].rows();
		c2=Phi.M[i][0].cols();
		// Building the FR //
		for(int j = 0; j < DimW; j++)
		{
			CRM[i][j].setZero(r1,r2);
			for(int k = 0; k < phy; k++)
			{
				TM1[k][j].setZero(c1,r2);
				TM2[k][j].setZero(c1,r2);
			}
		}
		for(tid = 0; tid < phy; tid++)
		{
			for(int j = 0; j < DimW; j++)
			{
				TM1[tid][j].noalias() = CRM[i+1][j] * Phi.M[i][tid].transpose();
			}
		}
		delete [] CRM[i+1];
		for(tid = 0; tid < phy; tid++)
		{
			for(int k = 0; k < phy; k++)
			{
				for(int l = 0; l < HH.H[i][tid*phy+k].r.size(); l++)
				{
					TM2[tid][HH.H[i][tid*phy+k].r[l]] += HH.H[i][phy*tid+k].v[l] * TM1[k][HH.H[i][tid*phy+k].c[l]];
				}
			}
		}
		if(i!=0)
		{
			for(int j = 0; j < DimW; j++)
			{
				for(int k = 0; k < phy; k++)
				{
					CRM[i][j].noalias() += Psi.M[i][k] * TM2[k][j];
				}
			}
		}else
		{
			// 1st site //
			for(int k = 0; k < phy; k++)
			{
				CRM[i][0] += Psi.M[i][k] * TM2[k][0];
			}
			EnergyExpec = CRM[0][0](0,0);
		}
	}
	////////////////////////////////////////////////
	delete [] CRM[0];
	delete [] CRM;
	for(int i = 0; i < phy; i++)
	{
		delete [] TM1[i];
		delete [] TM2[i];
	}
	delete [] TM1;
	delete [] TM2;
	////////////////////////////////////////////////
	return EnergyExpec;
}

double psiphi (const MPS& psi, const MPS& phi)
{
	int tid, pD, L;
	int r1, r2, c1, c2;
	double Expec = 0.0;
	if(psi.pD!=phi.pD)
	{
		std::cout<<"psiphi() error: Unmatched physical dimensions!"<<std::endl;
		abort();
	}
	pD = psi.pD;
	L = psi.Len;
	////////////////////////////////////////////
	Mxd * CRM = new Mxd [L];
	Mxd * TM1 = new Mxd [pD];
	////////////////////////////////////////////
	// Last site //
	r1=psi.M[L-1][0].rows();
	r2=phi.M[L-1][0].rows();
	CRM[L-1].setZero(r1,r2);
	for(int i = 0; i < pD; i++)
	{
		CRM[L-1] += psi.M[L-1][i] * phi.M[L-1][i].transpose();
	}
	// Last-but-1 to the second site //
	for(int i = L-2; i >= 0; i--) {
		r1=psi.M[i][0].rows();
		c1=psi.M[i][0].cols();
		r2=phi.M[i][0].rows();
		c2=phi.M[i][0].cols();
		CRM[i].setZero(r1,r2);
		for(int j = 0; j < pD; j++)
		{
			TM1[j].setZero(c1,r2);
		}
		// #pragma omp parallel num_threads(4) private(tid)
		for(tid = 0; tid < pD; tid++)
		{
			// tid = omp_get_thread_num();
			TM1[tid].noalias() = CRM[i+1] * phi.M[i][tid].transpose();
		}
		for(int j = 0; j < pD; j++)
		{
			CRM[i] += psi.M[i][j] * TM1[j];
		}
	}
	Expec = CRM[0](0,0);
	////////////////////////////////////////////////
	delete [] CRM;
	delete [] TM1;
	////////////////////////////////////////////////
	return Expec;
}

#endif
