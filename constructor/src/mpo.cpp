#ifndef My_MPO_CLASS
#define My_MPO_CLASS

#include "mpo.h"
#include <cmath>

typedef Eigen::MatrixXd Mxd;

MPO::MPO()
{
	norm = 1;
	if_init = false;
}

MPO::MPO(int l, int pd, int bd, double EnergyShifted)
{
	Len = l;
	pD = pd;
	bD = bd;
	ens = EnergyShifted;
	norm = 1;
	
	M_IDs = new int [Len]();
	for(int i = 0; i < Len; ++i)
	{
		M_IDs[i] = i;
	}
	
	Dim = new int [Len+1]();
	
	Dim[0] = 1;
	for(int i = 1; i < Len; ++i)
	{
		int pw = std::min(i,Len-i);
		Dim[i] = std::pow(pD*pD,pw);
		if(Dim[i]<0||log2(Dim[i])<pw)
		{
			Dim[i] = bD;
		}else
		{
			Dim[i] = std::min(Dim[i],bD);
		}
	}
	Dim[Len] = 1;
	
	M = new Mxd * [Len];
	for(int i = 0; i < Len; ++i)
	{
		M[i] = new Mxd [pD*pD];
	}
	
	H = new SpM * [Len];
	for(int i = 0; i < Len; i++)
	{
		H[i] = new SpM [pD*pD];
	}
	
	if_init = true;
}

void MPO::setMPO(int l, int pd, int bd, double EnergyShifted)
{
	Len = l;
	pD = pd;
	bD = bd;
	ens = EnergyShifted;
	norm = 1;
	
	M_IDs = new int [Len]();
	for(int i = 0; i < Len; ++i)
	{
		M_IDs[i] = i;
	}
	
	Dim = new int [Len+1]();
	
	Dim[0] = 1;
	for(int i = 1; i < Len; ++i)
	{
		int pw = std::min(i,Len-i);
		Dim[i] = std::pow(pD*pD,pw);
		if(Dim[i]<0||log2(Dim[i])<pw)
		{
			Dim[i] = bD;
		}else
		{
			Dim[i] = std::min(Dim[i],bD);
		}
	}
	Dim[Len] = 1;
	
	M = new Mxd * [Len];
	for(int i = 0; i < Len; ++i)
	{
		M[i] = new Mxd [pD*pD];
	}
	
	H = new SpM * [Len];
	for(int i = 0; i < Len; i++)
	{
		H[i] = new SpM [pD*pD];
	}
	
	if_init = true;
}

MPO::~MPO()
{
	if(if_init)
	{
		for(int i = 0; i < Len; ++i)
		{
			delete [] M[i];
			delete [] H[i];
		}
		delete [] M;
		delete [] H;
		delete [] Dim;
		delete [] M_IDs;
		if_init = false;
	}
}

void MPO::clearMPO()
{
	if(if_init)
	{
		for(int i = 0; i < Len; ++i)
		{
			delete [] M[i];
			delete [] H[i];
		}
		delete [] M;
		delete [] H;
		delete [] Dim;
		delete [] M_IDs;
		if_init = false;
	}
}

void MPO::copyMPO(const MPO& other)
{
	Len  = other.Len;
	pD   = other.pD;
	bD   = other.bD;
	norm = other.norm;
	
	if(!if_init)
	{
		M_IDs = new int [Len]();
		
		Dim = new int [Len+1]();
		
		M = new Mxd * [Len];
		for(int i = 0; i < Len; ++i)
		{
			M[i] = new Mxd [pD*pD];
		}
		
		H = new SpM * [Len];
		for(int i = 0; i < Len; i++)
		{
			H[i] = new SpM [pD*pD];
		}
		
		if_init = true;
	}
	
	for(int i = 0; i < Len; ++i)
	{
		M_IDs[i] = other.M_IDs[i];
	}
	
	for(int i = 0; i < Len+1; ++i)
	{
		Dim[i] = other.Dim[i];
	}
	
	for(int i = 0; i < Len; ++i)
	{
		for(int j = 0; j < pD*pD; ++j)
		{
			M[i][j] = other.M[i][j];
		}
	}
}

void MPO::decimateCopy(const MPO& other, int site, int phy)
{
	if(other.Len>1)
	{
		Len  = other.Len - 1;
		pD   = other.pD;
		bD   = other.bD;
		norm = other.norm;
	
		if(!if_init)
		{
			M_IDs = new int [Len]();
			
			Dim = new int [Len+1]();
		
			M = new Mxd * [Len];
			for(int i = 0; i < Len; ++i)
			{
				M[i] = new Mxd [pD*pD];
			}
		
			H = new SpM * [Len];
			for(int i = 0; i < Len; i++)
			{
				H[i] = new SpM [pD*pD];
			}
		
			if_init = true;
		}
		
		for(int i = 0; i < site; ++i)
		{
			M_IDs[i] = other.M_IDs[i];
		}
		
		for(int i = site; i < Len; ++i)
		{
			M_IDs[i] = other.M_IDs[i+1];
		}
		
		// for(int i = 0; i < other.Len; ++i)
		// {
		// 	std::cout<<other.M_IDs[i]<<" ";
		// }
		// std::cout<<std::endl;
		//
		// for(int i = 0; i < Len; ++i)
		// {
		// 	std::cout<<M_IDs[i]<<" ";
		// }
		// std::cout<<std::endl;
	
		for(int i = 0; i < site; ++i)
		{
			for(int j = 0; j < pD*pD; ++j)
			{
				M[i][j] = other.M[i][j];
			}
		}
		for(int i = site+1; i < other.Len; ++i)
		{
			for(int j = 0; j < pD*pD; ++j)
			{
				M[i-1][j] = other.M[i][j];
			}
		}
	
		if(other.M[site][0].rows()>other.M[site][0].cols())
		{
			if(phy>=0)
			{
				Mxd tp;
				for(int j = 0; j < pD*pD; ++j)
				{
					tp.noalias() = M[site-1][j] * other.M[site][phy*pD+phy];
					M[site-1][j] = tp;
				}
			}else
			{
				Mxd tp,tpp;
				tpp = other.M[site][0];
				for(int j = 1; j < pD; ++j)
				{
					tpp += other.M[site][j*pD+j];
				}
				tpp /= pD;
				for(int j = 0; j < pD*pD; ++j)
				{
					tp.noalias() = M[site-1][j] * tpp;
					M[site-1][j] = tp;
				}
			}
		}else
		{
			if(phy>=0)
			{
				Mxd tp;
				for(int j = 0; j < pD*pD; ++j)
				{
					tp.noalias() = other.M[site][phy*pD+phy] * M[site][j] ;
					M[site][j] = tp;
				}
			}else
			{
				Mxd tp,tpp;
				tpp = other.M[site][0];
				for(int j = 1; j < pD; ++j)
				{
					tpp += other.M[site][j*pD+j];
				}
				tpp /= pD;
				for(int j = 0; j < pD*pD; ++j)
				{
					tp.noalias() = tpp * M[site][j] ;
					M[site][j] = tp;
				}
			}

		}
	
		for(int i = 0; i < Len; ++i)
		{
			Dim[i] = M[i][0].rows();
		}
		Dim[Len] =1;
	
		bD = *std::max_element(Dim,Dim+Len+1);
	}
}

void MPO::setZero(std::string s)
{
	norm = 1;
	for(int i = 0; i < pD*pD; ++i)
	{
		M[0][i].setZero(1,bD);
	}
	for(int i = 1; i < Len-1; ++i)
	{
		for(int j = 0; j < pD*pD; ++j)
		{
			M[i][j].setZero(bD,bD);
		}
	}
	for(int i = 0; i < pD*pD; ++i)
	{
		M[Len-1][i].setZero(bD,1);
	}
	Dim[0] = 1;
	for(int i = 1; i < Len; ++i)
	{
		Dim[i] = bD;
	}
	Dim[Len] = 1;
}

void MPO::setZero()
{
	norm = 1;
	for(int i = 0; i < pD*pD; ++i)
	{
		M[0][i].setZero();
	}
	for(int i = 1; i < Len-1; ++i)
	{
		for(int j = 0; j < pD*pD; ++j)
		{
			M[i][j].setZero();
		}
	}
	for(int i = 0; i < pD*pD; ++i)
	{
		M[Len-1][i].setZero();
	}
}

void MPO::setZero(int nbD)
{
	norm = 1;
	bD   = nbD;
	
	Dim[0] = 1;
	for(int i = 1; i < Len; ++i)
	{
		int pw = std::min(i,Len-i);
		Dim[i] = std::pow(pD*pD,pw);
		if(Dim[i]<0||log2(Dim[i])<pw)
		{
			Dim[i] = bD;
		}else
		{
			Dim[i] = std::min(Dim[i],bD);
		}
	}
	Dim[Len] = 1;
	
	for(int i = 0; i < pD*pD; ++i)
	{
		M[0][i].setZero(1,Dim[1]);
	}
	for(int i = 1; i < Len-1; ++i)
	{
		for(int j = 0; j < pD*pD; ++j)
		{
			M[i][j].setZero(Dim[i],Dim[i+1]);
		}
	}
	for(int i = 0; i < pD*pD; ++i)
	{
		M[Len-1][i].setZero(Dim[Len-1],1);
	}
}

void MPO::setIdentity()
{
	setZero();
	for(int i = 0; i < pD; ++i)
	{
		M[0][i*pD+i].setIdentity();
	}
	for(int i = 1; i < Len-1; ++i)
	{
		for(int j = 0; j < pD; ++j)
		{
			M[i][j*pD+j].setIdentity();
		}
	}
	for(int i = 0; i < pD; ++i)
	{
		M[Len-1][i*pD+i].setIdentity();
	}
}

void MPO::setIdentity(int nbD)
{
	setZero(nbD);
	setIdentity();
}

void MPO::setRand()
{
	norm = 1;
	for(int i = 0; i < pD*pD; ++i)
	{
		M[0][i].setRandom();
	}
	for(int i = 1; i < Len-1; ++i)
	{
		for(int j = 0; j < pD*pD; ++j)
		{
			M[i][j].setRandom();
		}
	}
	for(int i = 0; i < pD*pD; ++i)
	{
		M[Len-1][i].setRandom(bD,1);
	}
	///////////////////////////////
	// Make symmetric
	for(int j = 0; j < Len; ++j)
	{
		for(int i = 0; i < pD; ++i)
		{
			for(int k = i+1; k < pD; ++k)
			{
				M[j][i*pD+k] = M[j][k*pD+i];
			}
		}
	}
	///////////////////////////////
	for(int i = 0; i < pD*pD; ++i)
	{
		M[Len-1][i].setRandom();
	}
}

void MPO::setRand(int nbD)
{
	norm = 1;
	bD   = nbD;
	
	Dim[0] = 1;
	for(int i = 1; i < Len; ++i)
	{
		int pw = std::min(i,Len-i);
		Dim[i] = std::pow(pD*pD,pw);
		if(Dim[i]<0||log2(Dim[i])<pw)
		{
			Dim[i] = bD;
		}else
		{
			Dim[i] = std::min(Dim[i],bD);
		}
	}
	Dim[Len] = 1;
	
	for(int i = 0; i < pD*pD; ++i)
	{
		M[0][i].setRandom(1,Dim[1]);
	}
	for(int i = 1; i < Len-1; ++i)
	{
		for(int j = 0; j < pD*pD; ++j)
		{
			M[i][j].setRandom(Dim[i],Dim[i+1]);
		}
	}
	for(int i = 0; i < pD*pD; ++i)
	{
		M[Len-1][i].setRandom(Dim[Len-1],1);
	}
	///////////////////////////////
	// Make symmetric
	for(int j = 0; j < Len; ++j)
	{
		for(int i = 0; i < pD; ++i)
		{
			for(int k = i+1; k < pD; ++k)
			{
				M[j][i*pD+k] = M[j][k*pD+i];
			}
		}
	}
	///////////////////////////////
}

void MPO::buildHeisenberg(double* J, double* h)
{
	setZero("Square");
	
	M[0][0]<<.5*h[0],0,0,.5,1;
	M[0][1]<<0,.5,.5,0,0;
	M[0][2]<<0,.5,-.5,0,0;
	M[0][3]<<-.5*h[0],0,0,-.5,1;
	for(int i = 1; i < Len-1; i++)
	{
		M[i][0]<<1,0,0,0,0,  0,0,0,0,0,   0,0,0,0,0,   .5,0,0,0,0,   .5*h[i],0,0,.5,1;
		M[i][1]<<0,0,0,0,0,  .5,0,0,0,0,  -.5,0,0,0,0,   0,0,0,0,0,   0,.5,.5,0,0;
		M[i][2]<<0,0,0,0,0,  .5,0,0,0,0,   .5,0,0,0,0,   0,0,0,0,0,   0,.5,-.5,0,0;
		M[i][3]<<1,0,0,0,0,  0,0,0,0,0,   0,0,0,0,0,  -.5,0,0,0,0,  -.5*h[i],0,0,-.5,1;
	}
	M[Len-1][0]<<1,0,0,.5,.5*h[Len-1];
	M[Len-1][1]<<0,.5,-.5,0,0;
	M[Len-1][2]<<0,.5,.5,0,0;
	M[Len-1][3]<<1,0,0,-.5,-.5*h[Len-1];
}

void MPO::buildIsing(double* J, double* h, double* Jp)
{
	setZero("Square");
	
	M[0][0]<<h[0],0,Jp[0],1;
	M[0][1]<<0,J[0],0,0;
	M[0][2]<<0,J[0],0,0;
	M[0][3]<<-h[0],0,-Jp[0],1;
	for(int i = 1; i < Len-1; i++)
	{
		M[i][0]<<1,0,0,0,  0,0,0,0,  1,0,0,0,   h[i],0,Jp[i],1;
		M[i][1]<<0,0,0,0,  1,0,0,0,  0,0,0,0,   0,J[i],0,0;
		M[i][2]<<0,0,0,0,  1,0,0,0,  0,0,0,0,   0,J[i],0,0;
		M[i][3]<<1,0,0,0,  0,0,0,0, -1,0,0,0,  -h[i],0,-Jp[i],1;
	}
	M[Len-1][0]<<1,0,1,h[Len-1];
	M[Len-1][1]<<0,1,0,0;
	M[Len-1][2]<<0,1,0,0;
	M[Len-1][3]<<1,0,-1,-h[Len-1];
}

void MPO::buildRTIC(double* dJ, double* dh)
{
	norm = 1;
	double enps = (ens)/double(Len);
	J1 = -1;
	J2 = 0.3;
	double h = 0.6;
	
	setZero("Square");
	
	M[0][0]<<-enps+h+dh[0],1,0,0;
	M[0][1]<<0,0,(J1+dJ[0]),J2;
	M[0][2]<<0,0,(J1+dJ[0]),J2;
	M[0][3]<<-enps-h-dh[0],1,0,0;
	for(int i = 1; i < Len-1; i++)
	{
		M[i][0]<<1,0,0,0,  -enps+h+dh[i],1,0,0,   0,0,0,0,   0,0,1,0;
		M[i][1]<<0,0,0,0,  0,0,(J1+dJ[i]),J2,     1,0,0,0,   0,0,0,0;
		M[i][2]<<0,0,0,0,  0,0,(J1+dJ[i]),J2,     1,0,0,0,   0,0,0,0;
		M[i][3]<<1,0,0,0,  -enps-h-dh[i],1,0,0,   0,0,0,0,   0,0,1,0;
	}
	M[Len-1][0]<<1,-enps+h+dh[Len-1],0,0;
	M[Len-1][1]<<0,0,1,0;
	M[Len-1][2]<<0,0,1,0;
	M[Len-1][3]<<1,-enps-h-dh[Len-1],0,0;
}

// void MPO::buildRTIC(double* dJ, double* dh)
// {
// 	double edgeSz = 0;
// 	norm = 1;
// 	double enps = (ens)/double(Len);
// 	J1 = -1;
// 	J2 = 0.3;
// 	double h = 0.6;
//
// 	setZero("Square");
//
// 	M[0][0]<<-enps+edgeSz,1,(J1+dJ[0]),J2;
// 	M[0][1]<<h+dh[0],0,0,0;
// 	M[0][2]<<h+dh[0],0,0,0;
// 	M[0][3]<<-enps-edgeSz,1,-(J1+dJ[0]),-J2;
// 	for(int i = 1; i < Len-1; i++)
// 	{
// 		M[i][0]<<1,0,0,0,  -enps,1,(J1+dJ[i]),J2,   1,0,0,0,   0,0,1,0;
// 		M[i][1]<<0,0,0,0,  h+dh[i],0,0,0,                 0,0,0,0,   0,0,0,0;
// 		M[i][2]<<0,0,0,0,  h+dh[i],0,0,0,                 0,0,0,0,   0,0,0,0;
// 		M[i][3]<<1,0,0,0,  -enps,1,-(J1+dJ[i]),-J2, -1,0,0,0,  0,0,1,0;
// 	}
// 	M[Len-1][0]<<1,-enps+edgeSz,1,0;
// 	M[Len-1][1]<<0,h+dh[Len-1],0,0;
// 	M[Len-1][2]<<0,h+dh[Len-1],0,0;
// 	M[Len-1][3]<<1,-enps-edgeSz,-1,0;
// }

void MPO::buildSPTC(double* dJ, double* dh)
{
	norm = 1;
	double enps = (ens)/double(Len);
	J1 = 1;
	J2 = 0.2;
	double h = 0.2;

	setZero("Square");

	M[0][0]<<-enps+h+dh[0],1,J2,0,0;
	M[0][1]<<0,0,0,0,J1+dJ[0];
	M[0][2]<<0,0,0,0,J1+dJ[0];
	M[0][3]<<-enps-h-dh[0],1,-J2,0,0;
	for(int i = 1; i < Len-1; i++)
	{
		M[i][0]<<1,0,0,0,0,  -enps+h+dh[i],1,J2,0,0,   1,0,0,0,0,   0,0,0,0,0,    0,0,0,1,0;
		M[i][1]<<0,0,0,0,0,  0,0,0,0,(J1+dJ[i]),       0,0,0,0,0,   1,0,0,0,0,    0,0,0,0,0;
		M[i][2]<<0,0,0,0,0,  0,0,0,0,(J1+dJ[i]),       0,0,0,0,0,   1,0,0,0,0,    0,0,0,0,0;
		M[i][3]<<1,0,0,0,0,  -enps-h-dh[i],1,-J2,0,0,  -1,0,0,0,0,  0,0,0,0,0,    0,0,0,-1,0;
	}
	M[Len-1][0]<<1,-enps+h+dh[Len-1],1,0,0;
	M[Len-1][1]<<0,0,0,1,0;
	M[Len-1][2]<<0,0,0,1,0;
	M[Len-1][3]<<1,-enps-h-dh[Len-1],-1,0,0;
}

// void MPO::buildSPTC(double* dJ, double* dh)
// {
// 	double enps = (ens)/double(Len);
// 	J1 = 1;
// 	J2 = 0.0;
// 	double h = 0.1;
//
// 	setZero("Square");
//
// 	M[0][0]<<-enps,1,0,0,J1+dJ[0];
// 	M[0][1]<<h+dh[0],0,J2,0,0;
// 	M[0][2]<<h+dh[0],0,J2,0,0;
// 	M[0][3]<<-enps,1,0,0,-J1-dJ[0];
// 	for(int i = 1; i < Len-1; i++)
// 	{
// 		M[i][0]<<1,0,0,0,0,  -enps,1,0,0,(J1+dJ[i]),     0,0,0,0,0,    1,0,0,0,0,    0,0,0,0,0;
// 		M[i][1]<<0,0,0,0,0,  h+dh[i],0,J2,0,0,           1,0,0,0,0,    0,0,0,0,0,    0,0,0,1,0;
// 		M[i][2]<<0,0,0,0,0,  h+dh[i],0,J2,0,0,           1,0,0,0,0,    0,0,0,0,0,    0,0,0,1,0;
// 		M[i][3]<<1,0,0,0,0,  -enps,1,0,0,-(J1+dJ[i]),    0,0,0,0,0,   -1,0,0,0,0,    0,0,0,0,0;
// 	}
// 	M[Len-1][0]<<1,-enps,0,1,0;
// 	M[Len-1][1]<<0,h+dh[Len-1],1,0,0;
// 	M[Len-1][2]<<0,h+dh[Len-1],1,0,0;
// 	M[Len-1][3]<<1,-enps,0,-1,0;
// }

void MPO::buildSz(int site)
{
	setZero("Square");
	for(int i = 0; i < Len; ++i)
	{
		if(i!=site)
		{
			M[i][0].setIdentity();
			M[i][3].setIdentity();
		}else
		{
			M[i][0].setIdentity();
			M[i][0] *= 1;
			M[i][3].setIdentity();
			M[i][3] *= -1;
		}

	}
}

void MPO::buildSza(int a)
{
    setZero("Square");
    M[0][0]<< 1, 1;
    M[0][1]<< 0, 0;
    M[0][2]<< 0, 0;
    M[0][3]<<-1, 1;
    for(int i = 1; i < a; i++)
    {
        M[i][0]<<1, 0,  1, 1;
        M[i][1]<<0, 0,  0, 0;
        M[i][2]<<0, 0,  0, 0;
        M[i][3]<<1, 0, -1, 1;
    }
    for(int i = a; i < Len-1; i++)
    {
        M[i][0]<<1, 0,  0, 1;
        M[i][1]<<0, 0,  0, 0;
        M[i][2]<<0, 0,  0, 0;
        M[i][3]<<1, 0,  0, 1;
    }
    if(a == Len)
    {
        M[Len-1][0]<< 1,  1;
        M[Len-1][1]<< 0,  0;
        M[Len-1][2]<< 0,  0;
        M[Len-1][3]<< 1, -1;
    }
    else
    {
        M[Len-1][0]<< 1, 0;
        M[Len-1][1]<< 0, 0;
        M[Len-1][2]<< 0, 0;
        M[Len-1][3]<< 1, 0;

    }
}

void MPO::buildSpMPO()
{
	if(!if_init)
	{
		std::cout<<"MPO not initiated!"<<std::endl;
		abort();
	}
	
	for(int site = 0; site < Len; site++)
	{
		for(int i = 0; i < pD*pD; i++)
		{
			if(H[site][i].v.size()!=0) H[site][i].v.clear();
			if(H[site][i].r.size()!=0) H[site][i].r.clear();
			if(H[site][i].c.size()!=0) H[site][i].c.clear();
		}
		for(int l = 0; l < pD*pD; l++)
		{
			for(int j = 0; j < M[site][0].rows(); j++)
			{
				for(int k = 0; k < M[site][0].cols(); k++)
				{
					if(M[site][l](j,k)!=0)
					{
						H[site][l].r.push_back(j);
						H[site][l].c.push_back(k);
						H[site][l].v.push_back(M[site][l](j,k));
					}
				}
			}
		}
	}
}

void MPO::square()
{
	norm *= norm;
	if(!if_init)
	{
		std::cout<<"MPO not initiated!"<<std::endl;
		abort();
	}
	
	Mxd ** tM = new Mxd * [Len];
	for(int i = 0; i < Len; ++i)
	{
		tM[i] = new Mxd [pD*pD];
	}
	for(int i = 0; i < Len; ++i)
	{
		for(int j = 0; j < pD*pD; ++j)
		{
			tM[i][j]=M[i][j];
		}
	}
	bD *= bD;
	setZero("Square");
	
	for(int i = 0; i < Len+1; ++i)
	{
		Dim[i] = Dim[i]*Dim[i];
	}
	
	for(int i = 0; i < Len; i++)
	{
		int r = tM[i][0].rows();
		int c = tM[i][0].cols();
		for(int p1 = 0; p1 < pD; p1++)
		{
			for(int p3 = 0; p3 < pD; p3++)
			{
				for(int j = 0; j < r; j++)
				{
					for(int k = 0; k < c; k++)
					{
						for(int p2 = 0; p2 < pD; p2++)
						{
							M[i][p1*pD+p3].block(j*r,k*c,r,c) += tM[i][p1*pD+p2](j,k) * tM[i][p2*pD+p3];
						} 
					}
				}
			}
		}
	}
	for(int i = 0; i < Len; ++i)
	{
		delete [] tM[i];
	}
	delete [] tM;
}

void MPO::RC()
{
	if(!if_init)
	{
		std::cout<<"MPO not initiated!"<<std::endl;
		abort();
	}
	if(norm==0) norm = 1;
	int tid;
	int row, col;
	Mxd TM, Q, R;
	for(int i = Len-1; i >= 0; i--)
	{
		row=M[i][0].rows();
		col=M[i][0].cols();
		TM.resize(row,pD*pD*col);
		for(tid=0; tid<pD*pD; tid++)
		//#pragma omp parallel num_threads(pD) private(tid)
		{
			//tid = omp_get_thread_num();
			TM.block(0,tid*col,row,col)=M[i][tid];
		}
		TM.transposeInPlace();
		rQR(TM,Q);
		Q.transposeInPlace();
		R.noalias() = TM.transpose()*Q.transpose();
		for(tid=0; tid<pD*pD; tid++)
		// #pragma omp parallel num_threads(pD) private(tid)
		{
			// tid = omp_get_thread_num();
			M[i][tid].setZero();
			M[i][tid].block(0,0,std::min(row,int(Q.cols())),col)=Q.block(0,tid*col,std::min(row,int(Q.cols())),col);
			if(i!=0)
			{
				Mxd tempM;
				tempM.noalias()=M[i-1][tid]*R;
				M[i-1][tid].setZero();
				M[i-1][tid].block(0,0,tempM.rows(),tempM.cols())=tempM;
			}else
			{
				norm = R(0,0);
			}
		}
	}
}

void MPO::LC()
{
	if(!if_init)
	{
		std::cout<<"MPO not initiated!"<<std::endl;
		abort();
	}
	if(norm==0) norm = 1;
	int tid;
	int row, col;
	Mxd TM, Q, R;
	for(int i = 0; i < Len; i++)
	{
		row=M[i][0].rows();
		col=M[i][0].cols();
		TM.resize(pD*pD*row,col);
		for(tid=0; tid<pD*pD; tid++)
		//#pragma omp parallel num_threads(pD) private(tid)
		{
			//tid = omp_get_thread_num();
			TM.block(tid*row,0,row,col)=M[i][tid];
		}
		rQR(TM,Q);
		R.noalias() = Q.transpose()*TM;
		for(tid=0; tid<pD*pD; tid++)
		// #pragma omp parallel num_threads(pD) private(tid)
		{
			// tid = omp_get_thread_num();
			M[i][tid].setZero();
			M[i][tid].block(0,0,row,std::min(col,int(Q.rows())))=Q.block(tid*row,0,row,std::min(col,int(Q.rows())));
			if(i!=Len-1)
			{
				Mxd tempM;
				tempM.noalias()=R*M[i+1][tid];
				M[i+1][tid].setZero();
				M[i+1][tid].block(0,0,tempM.rows(),tempM.cols())=tempM;
			}else
			{
				norm = R(0,0);
			}
		}
	}
}

void MPO::moveLeft(int site)
{
	if(!if_init)
	{
		std::cout<<"MPS not initiated!"<<std::endl;
		abort();
	}
	
	int tid;
	int row, col;
	Mxd TM, Q, R;
	row=M[site][0].rows();
	col=M[site][0].cols();
	TM.resize(row,pD*pD*col);
	for(tid=0; tid<pD*pD; tid++)
	//#pragma omp parallel num_threads(pD) private(tid)
	{
		//tid = omp_get_thread_num();
		TM.block(0,tid*col,row,col)=M[site][tid];
	}
	TM.transposeInPlace();
	rQR(TM,Q);
	Q.transposeInPlace();
	R.noalias() = TM.transpose()*Q.transpose();
	for(tid=0; tid<pD*pD; tid++)
	// #pragma omp parallel num_threads(pD) private(tid)
	{
		// tid = omp_get_thread_num();
		M[site][tid].setZero();
		M[site][tid].block(0,0,std::min(row,int(Q.cols())),col)=Q.block(0,tid*col,std::min(row,int(Q.cols())),col);
		if(site>0)
		{
			Mxd tempM;
			tempM.noalias()=M[site-1][tid]*R;
			M[site-1][tid].setZero();
			M[site-1][tid].block(0,0,tempM.rows(),tempM.cols())=tempM;
		}
	}
}

void MPO::moveRight(int site)
{
	if(!if_init)
	{
		std::cout<<"MPS not initiated!"<<std::endl;
		abort();
	}
	
	int tid;
	int row, col;
	Mxd TM, Q, R;
	row=M[site][0].rows();
	col=M[site][0].cols();
	TM.resize(pD*pD*row,col);
	for(tid=0; tid<pD*pD; tid++)
	//#pragma omp parallel num_threads(pD) private(tid)
	{
		//tid = omp_get_thread_num();
		TM.block(tid*row,0,row,col)=M[site][tid];
	}
	rQR(TM,Q);
	R.noalias() = Q.transpose()*TM;
	for(tid=0; tid<pD*pD; tid++)
	// #pragma omp parallel num_threads(pD) private(tid)
	{
		// tid = omp_get_thread_num();
		M[site][tid].setZero();
		M[site][tid].block(0,0,row,std::min(col,int(Q.rows())))=Q.block(tid*row,0,row,std::min(col,int(Q.rows())));
		if(site<Len-1)
		{
			Mxd tempM;
			tempM.noalias()=R*M[site+1][tid];
			M[site+1][tid].setZero();
			M[site+1][tid].block(0,0,tempM.rows(),tempM.cols())=tempM;
		}
	}
}

void MPO::compressL(int nbD)
{
	bD = *std::max_element(Dim,Dim+Len+1);
	if(bD>nbD)
	{
		RC();
		int tid;
		int row, col;
		double svd_error = 0;
		double * sv = new double [pD*pD*bD];
		/////////////////////////////////////////////
		Mxd TM, U, V;
		///////////////////////////////////////////
		for(int i = 0; i < Len; i++)
		{
			row=M[i][0].rows();
			col=M[i][0].cols();
			TM.resize(pD*pD*row,col);
			for(tid = 0; tid < pD*pD; tid++)
			{
				// tid = omp_get_thread_num();
				TM.block(tid*row,0,row,col)=M[i][tid];
			}
			rSVD(TM,pD*pD*bD,sv,U,V,nbD,svd_error);
			for(tid = 0; tid < pD*pD; tid++)
			{
				// tid = omp_get_thread_num();
				M[i][tid].setZero(row,std::min(nbD,col));
				M[i][tid].block(0,0,row,U.cols())=U.block(tid*row,0,row,U.cols());
				if(i!=Len-1)
				{
					Mxd tempM;
					tempM.noalias()=V*M[i+1][tid].real();
					M[i+1][tid].setZero(nbD,tempM.cols());
					M[i+1][tid].block(0,0,tempM.rows(),tempM.cols())=tempM;
				}
				// else
				// {
				// 	norm = V(0,0);
				// }
			}
			// cout<<svd_error<<endl;
		}
		delete [] sv;
		bD = nbD;
		for(int i = 0; i < Len; ++i)
		{
			Dim[i] = M[i][0].rows();
		}
		
		for(int i = 0; i < pD*pD; ++i)
		{
			M[0][i] *= norm;
		}
	}
}

void MPO::addMPO(double coeff, const MPO& other)
{
	if(!if_init||!other.if_init)
	{
		std::cout<<"MPO(s) to be added is(are) not initiated!"<<std::endl;
		abort();
	}
	if(Len==other.Len&&pD==other.pD)
	{
		
		Mxd tp;
		
		for(int i = 0; i < pD*pD; ++i)
		{
			tp.setZero(1,bD+other.bD);
			tp.block(0,0,1,bD) = M[0][i];
			tp.block(0,bD,1,other.bD) = coeff * other.M[0][i];
			M[0][i] = tp;
		}
	
		for(int i = 1; i < Len-1; ++i)
		{
			for(int j = 0; j < pD*pD; ++j)
			{
				tp.setZero(bD+other.bD,bD+other.bD);
				tp.block(0,0,bD,bD) = M[i][j];
				tp.block(bD,bD,other.bD,other.bD) = other.M[i][j];
				M[i][j] = tp;
			}
		}
		
		for(int i = 0; i < pD*pD; ++i)
		{
			tp.setZero(bD+other.bD,1);
			tp.block(0,0,bD,1) = M[Len-1][i];
			tp.block(bD,0,other.bD,1) = other.M[Len-1][i];
			M[Len-1][i] = tp;
		}
		
		for(int i = 0; i < Len; ++i)
		{
			Dim[i] = M[i][0].rows();
		}
		
		bD = *std::max_element(Dim,Dim+Len+1);
	}
}

void MPO::EE()
{
	if(!if_init)
	{
		std::cout<<"MPO not initiated!"<<std::endl;
		abort();
	}
	// std::cout<<"Entanglement Entropy:"<<std::endl;
	///////////////////////////////////////////
	RC();
	double * sv = new double [pD*pD*bD]();
	int row, col, tid, DT, tDim;
	Mxd TM, U, V;
	double vEE = 0;
	///////////////////////////////////////////
	for(int i = 0; i < Len-1; i++)
	{
		row=M[i][0].rows();
		col=M[i][0].cols();
		
		TM.resize(pD*pD*row,col);
		
		tDim = std::min(TM.rows(),TM.cols());
		
		for(tid = 0; tid < pD*pD; tid++)
		{
			TM.block(tid*row,0,row,col) = M[i][tid].block(0,0,row,col);
		}
		
		rSVD(TM,tDim,sv,U,V,'r');
		
		double EEn = 0;
		for(int j = 0; j < tDim; ++j)
		{
			if(sv[j]>0) EEn -= sv[j]*sv[j]*log(sv[j]*sv[j]); 
		}
		
		if(vEE<EEn) vEE = EEn;
		// if(i==(Len/2-1)) std::cout<<EEn<<" ";
		
		for(tid = 0; tid < pD*pD; tid++)
		{
			M[i][tid].setZero();
			M[i][tid].block(0,0,row,U.cols()) = U.block(tid*row,0,row,U.cols());
			if(i!=Len-1)
			{
				Mxd tempM;
				tempM.noalias() = V * M[i+1][tid];
				M[i+1][tid].setZero();
				M[i+1][tid].block(0,0,tempM.rows(),tempM.cols())=tempM;
			}
		}
	}
	std::cout<<vEE<<" ";
	delete [] sv;
}

void MPO::examineBond(double tol, bool dry_run)
{
  RC();
  int tid;
  int row, col;
  double svd_tol = tol;
  std::vector<double> sv;
  /////////////////////////////////////////////
  Mxd TM, U, V;
  ///////////////////////////////////////////
  for(int i = 0; i < Len-1; i++)
    {
      row=M[i][0].rows();
      col=M[i+1][0].cols();
      TM.resize(pD*pD*row,pD*pD*col);
      for(int j = 0; j < pD*pD; j++)
	{
	  for(int k = 0; k < pD*pD; ++k)
	    {
	      TM.block(j*row,k*col,row,col)=M[i][j] * M[i+1][k];
	    }
	}
      rSVD(TM,sv,U,V,svd_tol,dry_run);
      for(int j = 0; j < pD*pD; j++)
	{
	  M[i][j]   = U.block(j*row,0,row,U.cols());
	  M[i+1][j] = V.block(0,j*col,V.rows(),col);
	}
    }
  for(int i = 0; i < Len; ++i)
    {
      Dim[i] = M[i][0].rows();
    }
  for(int i = 0; i < pD*pD; ++i)
    {
      M[0][i] *= norm;
    }
}

double MPO::trace(std::vector<int>& pos)
{
	if(!if_init)
	{
		std::cout<<"MPO not initiated!"<<std::endl;
		abort();
	}
	
	for(int i = 0; i < pos.size(); ++i)
	{
		M[pos[i]][pD+pD-1] *= -1;
	}
	
	Mxd A = M[0][0]/2;
	for(int i = 1; i < pD; ++i)
	{
		A += M[0][i*pD+i]/2;
	}
	for(int i = 1; i < Len; ++i)
	{
		Mxd tp = M[i][0]/2;
		for(int j = 1; j < pD; ++j)
		{
			tp += M[i][j*pD+j]/2;
		}
		Mxd ttp = A * tp;
		A = ttp;
	}
	
	return A(0,0);
}

double MPO::sumAll()
{
	if(!if_init)
	{
		std::cout<<"MPO not initiated!"<<std::endl;
		abort();
	}
	
	Mxd A = M[0][0];
	for(int i = 1; i < pD; ++i)
	{
		A += M[0][i*pD+i];
	}
	for(int i = 1; i < Len; ++i)
	{
		Mxd tp = M[i][0];
		for(int j = 1; j < pD; ++j)
		{
			tp += M[i][j*pD+j];
		}
		Mxd ttp = A * tp;
		A = ttp;
	}
	
	return A(0,0);
}

void MPO::keepDiag()
{
	if(!if_init)
	{
		std::cout<<"MPO not initiated!"<<std::endl;
		abort();
	}
	
	for(int s = 0; s < Len; ++s)
	{
		for(int i = 0; i < pD*pD; ++i)
		{
			int p1 = i%pD;
			int p2 = i/pD;
			if(p1!=p2) M[s][i].setZero();
		}
	}
}

void MPO::diagToMPS(MPS& psi)
{
	psi.setMPS(Len, pD, bD);
	for(int s = 0; s < Len; ++s)
	{
		for(int i = 0; i < pD; ++i)
		{
			psi.M[s][i] = M[s][i*pD+i];
		}
	}
	psi.setDim();
}


#endif
