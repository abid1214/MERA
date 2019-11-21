#ifndef My_Iterative_Compression
#define My_Iterative_Compression

#include "iterCompress.h"
typedef Eigen::MatrixXd Mxd;

void buildR(MPS& H, MPS& Hc, Mxd* CR)
{
	// cout<<"buildR"<<endl;
	int pD = H.pD;
	int L  = H.Len;
	Mxd * TM1 = new Mxd [pD];
	////////////////////////////////////////////
	// Last site //
	int r1, r2, c1, c2;
	r1=Hc.M[L-1][0].rows();
	r2=H.M[L-1][0].rows();
	CR[L-1].setZero(r1,r2);
	for(int i = 0; i < pD; i++)
	{
		CR[L-1] += Hc.M[L-1][i] * H.M[L-1][i].transpose();
	}
	// Last-but-1 to the second site //
	for(int i = L-2; i > 0; i--)
	{
		r1=Hc.M[i][0].rows();
		c1=Hc.M[i][0].cols();
		r2=H.M[i][0].rows();
		c2=H.M[i][0].cols();
		CR[i].setZero(r1,r2);
		for(int j = 0; j < pD; j++)
		{
			TM1[j].setZero(c1,r2);
		}
		// #pragma omp parallel num_threads(4) private(tid)
		for(int tid = 0; tid < pD; tid++)
		{
			// tid = omp_get_thread_num();
			TM1[tid].noalias() = CR[i+1] * H.M[i][tid].transpose();
		}
		for(int j = 0; j < pD; j++)
		{
			CR[i] += Hc.M[i][j] * TM1[j];
		}
	}
	
	delete [] TM1;
}

void updateSite(MPS& H, MPS& Hc, Mxd* CL, Mxd* CR, int& site, char& direc, double& dt)
{
	// cout<<"updateSite"<<endl;
	double tp = 0;
	for(int i = 0; i < H.pD; i++)
	{
		Mxd T;
		if(site==0)
		{
			T = Hc.M[site][i];
			Hc.M[site][i] = H.M[site][i] * CR[site+1].transpose();
			tp += (T.transpose()*Hc.M[site][i]).trace();
		}else if(site==H.Len-1)
		{
			T = Hc.M[site][i];
			Hc.M[site][i] = CL[site-1] * H.M[site][i];
			tp += (T.transpose()*Hc.M[site][i]).trace();
		}else
		{
			T = Hc.M[site][i];
			Hc.M[site][i] = CL[site-1] * H.M[site][i] * CR[site+1].transpose();
			tp += (T.transpose()*Hc.M[site][i]).trace();
		}
	}
	
	dt = 1-tp;
	
	if(direc=='R')
		Hc.moveRight(site);
	else
		Hc.moveLeft(site);
}

void updateEnv(MPS& H, MPS& Hc, Mxd* CL, Mxd* CR, int& site, char& direc, double& dt)
{
	// cout<<"updateEnv"<<endl;
	int pD = H.pD;
	int L  = H.Len;
	if(site==0)
	{
		int c1, c2;
		c1=Hc.M[0][0].cols();
		c2=H.M[0][0].cols();
		CL[0].setZero(c1,c2);
		for(int i = 0; i < pD; i++)
		{
			CL[0] += Hc.M[0][i].transpose() * H.M[0][i];
		}
	}else if(site==L-1)
	{
		int r1, r2;
		r1=Hc.M[L-1][0].rows();
		r2=H.M[L-1][0].rows();
		CR[L-1].setZero(r1,r2);
		for(int i = 0; i < pD; i++)
		{
			CR[L-1] += Hc.M[L-1][i] * H.M[L-1][i].transpose();
		}
	}else
	{
		Mxd * TM1 = new Mxd [pD];
		int r1, r2, c1, c2;
		r1=Hc.M[site][0].rows();
		c1=Hc.M[site][0].cols();
		r2=H.M[site][0].rows();
		c2=H.M[site][0].cols();
		if(direc=='L')
		{
			CR[site].setZero(r1,r2);
			for(int j = 0; j < pD; j++)
			{
				TM1[j].setZero(c1,r2);
			}
			// #pragma omp parallel num_threads(4) private(tid)
			for(int tid = 0; tid < pD; tid++)
			{
				// tid = omp_get_thread_num();
				TM1[tid].noalias() = CR[site+1] * H.M[site][tid].transpose();
			}
			for(int j = 0; j < pD; j++)
			{
				CR[site] += Hc.M[site][j] * TM1[j];
			}
		}else
		{
			CL[site].setZero(c1,c2);
			for(int j = 0; j < pD; j++)
			{
				TM1[j].setZero(r1,c2);
			}
			// #pragma omp parallel num_threads(4) private(tid)
			for(int tid = 0; tid < pD; tid++)
			{
				// tid = omp_get_thread_num();
				TM1[tid].noalias() = CL[site-1] * H.M[site][tid];
			}
			for(int j = 0; j < pD; j++)
			{
				CL[site] += Hc.M[site][j].transpose() * TM1[j];
			}
		}
		
		delete [] TM1;
	}
}

void iterCompress(bool is_Random, int maxIter, double tol, int nbD, MPS& H)
{
	std::cout<<"Iterative compression... ("<<H.bD<<" --> "<<nbD<<")"<<std::endl;
	int L = H.Len;
	H.LC();
	MPS Hc(H.Len,H.pD,nbD);
	////////////////////////////////////////////
	Hc.setRand(nbD);
	Hc.LC();
	// if(is_Random)
	// {
	// 	Hc.setRand(nbD);
	// 	Hc.RC();
	// }else
	// {
	// 	Hc.copyMPO(H);
	// 	Hc.compressL(nbD);
	// }
	////////////////////////////////////////////
	Mxd * CRM = new Mxd [L];
	Mxd * CLM = new Mxd [L];
	////////////////////////////////////////////
	buildR(H, Hc, CRM);
	double dt = 1.0;
	////////////////////////////////////////////
	for(int ii = 0; ii < maxIter; ii++)
	{
		char direc;
		int site;
		
		if(ii%2==0)
			direc = 'R';
		else
			direc = 'L';
		
		for(int i = 0; i < L-1; i++)
		{
			if(direc=='R')
				site = i;
			else
				site = L-i-1;
			updateSite(H,Hc,CLM,CRM,site,direc,dt);
			updateEnv(H,Hc,CLM,CRM,site,direc,dt);
		}
		if(std::abs(dt)<tol) break;
	}
	std::cout<<"Compression error = "<<dt<<std::endl;
	// std::cout<<std::endl;
	////////////////////////////////////////////
	delete [] CRM;
	delete [] CLM;
	////////////////////////////////////////////
	H.copyMPS(Hc);
}

void iterCompress(bool is_Random, int maxIter, double tol, int nbD, MPS& H, bool print_info)
{
	if(print_info) std::cout<<"Iterative compression... ("<<H.bD<<" --> "<<nbD<<")"<<std::endl;
	int L = H.Len;
	// for(int j = 0; j < H.Len; ++j)
	// {
	// 	std::cout<<"("<<H.M[j][0].rows()<<" "<<H.M[j][0].cols()<<") ";
	// }
	// std::cout<<std::endl;
	// for(int j = 0; j < H.Len; ++j)
	// {
	// 	std::cout<<"("<<H.M[j][1].rows()<<" "<<H.M[j][1].cols()<<") ";
	// }
	// std::cout<<std::endl;
	// for(int j = 0; j < H.Len+1; ++j)
	// {
	// 	std::cout<<H.Dim[j]<<" ";
	// }
	// std::cout<<std::endl;
	H.LC();
	MPS Hc(H.Len,H.pD,nbD);
	////////////////////////////////////////////
	Hc.setRand(nbD);
	Hc.LC();
	// for(int j = 0; j < Hc.Len+1; ++j)
	// {
	// 	std::cout<<Hc.Dim[j]<<" ";
	// }
	// std::cout<<std::endl;
	// if(is_Random)
	// {
	// 	Hc.setRand(nbD);
	// 	Hc.RC();
	// }else
	// {
	// 	Hc.compressL(nbD);
	// }
	////////////////////////////////////////////
	Mxd * CRM = new Mxd [L];
	Mxd * CLM = new Mxd [L];
	////////////////////////////////////////////
	buildR(H, Hc, CRM);
	double dt = 1.0;
	////////////////////////////////////////////
	for(int ii = 0; ii < maxIter; ii++)
	{
		char direc;
		int site;
		
		if(ii%2==0)
			direc = 'R';
		else
			direc = 'L';
		
		for(int i = 0; i < L-1; i++)
		{
			if(direc=='R')
				site = i;
			else
				site = L-i-1;
			updateSite(H,Hc,CLM,CRM,site,direc,dt);
			updateEnv(H,Hc,CLM,CRM,site,direc,dt);
		}
		// std::cout<<"Compression error = "<<dt<<std::endl;
		if(std::abs(dt)<tol) break;
	}
	if(print_info) std::cout<<"Compression error = "<<dt<<std::endl;
	// std::cout<<std::endl;
	////////////////////////////////////////////
	delete [] CRM;
	delete [] CLM;
	////////////////////////////////////////////
	H.copyMPS(Hc);
}


///////////////////////////////////////////////////////////////////////////////


void buildR(MPO& H, MPO& Hc, Mxd* CR)
{
	// cout<<"buildR"<<endl;
	int pD = H.pD;
	int L  = H.Len;
	Mxd * TM1 = new Mxd [pD*pD];
	////////////////////////////////////////////
	// Last site //
	int r1, r2, c1, c2;
	r1=Hc.M[L-1][0].rows();
	r2=H.M[L-1][0].rows();
	CR[L-1].setZero(r1,r2);
	for(int i = 0; i < pD*pD; i++)
	{
		CR[L-1] += Hc.M[L-1][i] * H.M[L-1][i].transpose();
	}
	// Last-but-1 to the second site //
	for(int i = L-2; i > 0; i--)
	{
		r1=Hc.M[i][0].rows();
		c1=Hc.M[i][0].cols();
		r2=H.M[i][0].rows();
		c2=H.M[i][0].cols();
		CR[i].setZero(r1,r2);
		for(int j = 0; j < pD*pD; j++)
		{
			TM1[j].setZero(c1,r2);
		}
		// #pragma omp parallel num_threads(4) private(tid)
		for(int tid = 0; tid < pD*pD; tid++)
		{
			// tid = omp_get_thread_num();
			TM1[tid].noalias() = CR[i+1] * H.M[i][tid].transpose();
		}
		for(int j = 0; j < pD*pD; j++)
		{
			CR[i] += Hc.M[i][j] * TM1[j];
		}
	}
	
	delete [] TM1;
}

void updateSite(MPO& H, MPO& Hc, Mxd* CL, Mxd* CR, int& site, char& direc, double& dt)
{
	// cout<<"updateSite"<<endl;
	double tp = 0;
	for(int i = 0; i < H.pD*H.pD; i++)
	{
		Mxd T;
		if(site==0)
		{
			T = Hc.M[site][i];
			Hc.M[site][i] = H.M[site][i] * CR[site+1].transpose();
			tp += (T.transpose()*Hc.M[site][i]).trace();
		}else if(site==H.Len-1)
		{
			T = Hc.M[site][i];
			Hc.M[site][i] = CL[site-1] * H.M[site][i];
			tp += (T.transpose()*Hc.M[site][i]).trace();
		}else
		{
			T = Hc.M[site][i];
			Hc.M[site][i] = CL[site-1] * H.M[site][i] * CR[site+1].transpose();
			tp += (T.transpose()*Hc.M[site][i]).trace();
		}
	}
	
	dt = 1-tp;
	
	if(direc=='R')
		Hc.moveRight(site);
	else
		Hc.moveLeft(site);
}

void updateEnv(MPO& H, MPO& Hc, Mxd* CL, Mxd* CR, int& site, char& direc, double& dt)
{
	// cout<<"updateEnv"<<endl;
	int pD = H.pD;
	int L  = H.Len;
	if(site==0)
	{
		int c1, c2;
		c1=Hc.M[0][0].cols();
		c2=H.M[0][0].cols();
		CL[0].setZero(c1,c2);
		for(int i = 0; i < pD*pD; i++)
		{
			CL[0] += Hc.M[0][i].transpose() * H.M[0][i];
		}
	}else if(site==L-1)
	{
		int r1, r2;
		r1=Hc.M[L-1][0].rows();
		r2=H.M[L-1][0].rows();
		CR[L-1].setZero(r1,r2);
		for(int i = 0; i < pD*pD; i++)
		{
			CR[L-1] += Hc.M[L-1][i] * H.M[L-1][i].transpose();
		}
	}else
	{
		Mxd * TM1 = new Mxd [pD*pD];
		int r1, r2, c1, c2;
		r1=Hc.M[site][0].rows();
		c1=Hc.M[site][0].cols();
		r2=H.M[site][0].rows();
		c2=H.M[site][0].cols();
		if(direc=='L')
		{
			CR[site].setZero(r1,r2);
			for(int j = 0; j < pD*pD; j++)
			{
				TM1[j].setZero(c1,r2);
			}
			// #pragma omp parallel num_threads(4) private(tid)
			for(int tid = 0; tid < pD*pD; tid++)
			{
				// tid = omp_get_thread_num();
				TM1[tid].noalias() = CR[site+1] * H.M[site][tid].transpose();
			}
			for(int j = 0; j < pD*pD; j++)
			{
				CR[site] += Hc.M[site][j] * TM1[j];
			}
		}else
		{
			CL[site].setZero(c1,c2);
			for(int j = 0; j < pD*pD; j++)
			{
				TM1[j].setZero(r1,c2);
			}
			// #pragma omp parallel num_threads(4) private(tid)
			for(int tid = 0; tid < pD*pD; tid++)
			{
				// tid = omp_get_thread_num();
				TM1[tid].noalias() = CL[site-1] * H.M[site][tid];
			}
			for(int j = 0; j < pD*pD; j++)
			{
				CL[site] += Hc.M[site][j].transpose() * TM1[j];
			}
		}
		
		delete [] TM1;
	}
}

void iterCompress(bool is_Random, int maxIter, double tol, int nbD, MPO& H)
{
	std::cout<<"Iterative compression... ("<<H.bD<<" --> "<<nbD<<")"<<std::endl;
	int L = H.Len;
	H.RC();
	double norm = H.norm;
	MPO Hc(H.Len,H.pD,nbD,0);
	Hc.copyMPO(H);
	////////////////////////////////////////////
	if(is_Random)
	{
		Hc.setRand(nbD);
		Hc.RC();
	}else
	{
		Hc.copyMPO(H);
		Hc.compressL(nbD);
	}
	////////////////////////////////////////////
	Mxd * CRM = new Mxd [L];
	Mxd * CLM = new Mxd [L];
	////////////////////////////////////////////
	buildR(H, Hc, CRM);
	double dt = 1.0;
	////////////////////////////////////////////
	for(int ii = 0; ii < maxIter; ii++)
	{
		char direc;
		int site;
		
		if(ii%2==0)
			direc = 'R';
		else
			direc = 'L';
		
		for(int i = 0; i < L-1; i++)
		{
			if(direc=='R')
				site = i;
			else
				site = L-i-1;
			updateSite(H,Hc,CLM,CRM,site,direc,dt);
			updateEnv(H,Hc,CLM,CRM,site,direc,dt);
		}
		if(std::abs(dt)<tol) break;
	}
	std::cout<<"Compression error = "<<dt<<std::endl;
	// std::cout<<std::endl;
	////////////////////////////////////////////
	delete [] CRM;
	delete [] CLM;
	////////////////////////////////////////////
	H.copyMPO(Hc);
	H.norm = norm;
	for(int i = 0; i < H.pD*H.pD; ++i)
	{
		H.M[0][i] *= norm;
	}
}

void iterCompress(bool is_Random, int maxIter, double tol, int nbD, MPO& H, bool print_info)
{
	if(print_info) std::cout<<"Iterative compression... ("<<H.bD<<" --> "<<nbD<<")"<<std::endl;
	int L = H.Len;
	H.RC();
	double norm = H.norm;
	MPO Hc(H.Len,H.pD,nbD,0);
	Hc.copyMPO(H);
	////////////////////////////////////////////
	if(is_Random)
	{
		Hc.setRand(nbD);
		Hc.RC();
	}else
	{
		Hc.compressL(nbD);
	}
	////////////////////////////////////////////
	Mxd * CRM = new Mxd [L];
	Mxd * CLM = new Mxd [L];
	////////////////////////////////////////////
	buildR(H, Hc, CRM);
	double dt = 1.0;
	////////////////////////////////////////////
	for(int ii = 0; ii < maxIter; ii++)
	{
		char direc;
		int site;
		
		if(ii%2==0)
			direc = 'R';
		else
			direc = 'L';
		
		for(int i = 0; i < L-1; i++)
		{
			if(direc=='R')
				site = i;
			else
				site = L-i-1;
			updateSite(H,Hc,CLM,CRM,site,direc,dt);
			updateEnv(H,Hc,CLM,CRM,site,direc,dt);
		}
		if(std::abs(dt)<tol) break;
	}
	if(print_info) std::cout<<"Compression error = "<<dt<<std::endl;
	// std::cout<<std::endl;
	////////////////////////////////////////////
	delete [] CRM;
	delete [] CLM;
	////////////////////////////////////////////
	H.copyMPO(Hc);
	H.norm = norm;
	for(int i = 0; i < H.pD*H.pD; ++i)
	{
		H.M[0][i] *= norm;
	}
}


#endif