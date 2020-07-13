// no-quantum-number mps class
#ifndef My_MPS_CLASS
#define My_MPS_CLASS

#include "mps.h"
#include "lapack_wrapper.h"
#include <cmath>


typedef Eigen::MatrixXd Mxd;


MPS::MPS()
{
	if_init = false;
}

MPS::MPS(int l, int pd, int bd)
{
	Len = l;
	pD = pd;
	bD = bd;
	Dim = new int [Len+1];
	if(bD>(1<<10))
	{
		std::cout<<"Bond Dimension is too big! "<<bD<<std::endl;
		abort();
	}
	
	Dim[0] = 1;
	for(int i = 1; i < Len; ++i)
	{
		int pw = std::min(i,Len-i);
		Dim[i] = std::pow(pD,pw);
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
		// std::cout<<i<<std::endl;
		M[i] = new Mxd [pD];
	}
	if_init = true;
}

void MPS::setMPS(int l, int pd, int bd)
{
	if(!if_init)
	{
		Len = l;
		pD = pd;
		bD = bd;
		Dim = new int [Len+1];
		if(bD>(1<<16))
		{
			std::cout<<"Bond Dimension is too big! "<<bD<<std::endl;
			abort();
		}
		Dim[0] = 1;
		for(int i = 1; i < Len; ++i)
		{
			int pw = std::min(i,Len-i);
			Dim[i] = std::pow(pD,pw);
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
			// std::cout<<i<<std::endl;
			M[i] = new Mxd [pD];
		}
		if_init = true;
	}
}

MPS::MPS (const MPS& other)
{
	Len = other.Len;
	pD = other.pD;
	bD = other.bD;
	Dim = new int [Len+1];
	M = new Mxd * [Len];
	for(int i = 0; i < Len; ++i)
	{
		M[i] = new Mxd [pD];
	}
	
	for(int i = 0; i < Len+1; ++i)
	{
		Dim[i] = other.Dim[i];
	}
	
	for(int i = 0; i < Len; ++i)
	{
		for(int j = 0; j < pD; ++j)
		{
			M[i][j] = other.M[i][j];
		}
	}
	if_init = true;
}

void MPS::setNorm(double nm)
{
	if(if_init)
	{
		for(int i = 0; i < pD; ++i)
		{
			M[0][i] *= nm;
		}
	}
}

void MPS::copyMPS(const MPS& other)
{
	Len = other.Len;
	pD = other.pD;
	bD = other.bD;
	
	if(!if_init)
	{
		Dim = new int [Len+1];
		M = new Mxd * [Len];
		for(int i = 0; i < Len; ++i)
		{
			M[i] = new Mxd [pD];
		}
		if_init = true;
	}
	
	for(int i = 0; i < Len+1; ++i)
	{
		Dim[i] = other.Dim[i];
	}
	
	for(int i = 0; i < Len; ++i)
	{
		for(int j = 0; j < pD; ++j)
		{
			M[i][j] = other.M[i][j];
		}
	}
}

void MPS::addMPS(double coeff, const MPS& other)
{
	if(!if_init||!other.if_init)
	{
		std::cout<<"MPS(s) to be added is(are) not initiated!"<<std::endl;
		abort();
	}
	if(Len==other.Len&&pD==other.pD)
	{
		bD += other.bD;
		
		Mxd tp;
		
		for(int i = 0; i < pD; ++i)
		{
			tp.setZero(1,Dim[1]+other.Dim[1]);
			tp.block(0,0,1,Dim[1]) = M[0][i];
			tp.block(0,Dim[1],1,other.Dim[1]) = coeff * other.M[0][i];
			M[0][i] = tp;
		}
	
		for(int i = 1; i < Len-1; ++i)
		{
			for(int j = 0; j < pD; ++j)
			{
				tp.setZero(Dim[i]+other.Dim[i],Dim[i+1]+other.Dim[i+1]);
				tp.block(0,0,Dim[i],Dim[i+1]) = M[i][j];
				tp.block(Dim[i],Dim[i+1],other.Dim[i],other.Dim[i+1]) = other.M[i][j];
				M[i][j] = tp;
			}
		}
		
		for(int i = 0; i < pD; ++i)
		{
			tp.setZero(Dim[Len-1]+other.Dim[Len-1],1);
			tp.block(0,0,Dim[Len-1],1) = M[Len-1][i];
			tp.block(Dim[Len-1],0,other.Dim[Len-1],1) = other.M[Len-1][i];
			M[Len-1][i] = tp;
		}
		
		for(int i = 1; i < Len; ++i)
		{
			Dim[i] += other.Dim[i];
		}
	}
}

MPS::~MPS()
{
	if(if_init)
	{
		for(int i = 0; i < Len; ++i)
		{
			delete [] M[i];
		}
		delete [] M;
		delete [] Dim;
		if_init = false;
	}
}

void MPS::clear()
{
	if(if_init)
	{
		for(int i = 0; i < Len; ++i)
		{
			delete [] M[i];
		}
		delete [] M;
		delete [] Dim;
		if_init = false;
	}
}

void MPS::addSite(int site)
{
	if(if_init)
	{
		if(site==0)
		{
			MPS phi(Len+1,pD,1);
			phi.setZero();
			for(int i = 0; i < pD; ++i) phi.M[0][i].setIdentity();
			for(int i = 0; i < Len; ++i)
			{
				for(int j = 0; j < pD; ++j)
				{
					phi.M[i+1][j] = M[i][j];
				}
			}
			phi.setDim();
			clear();
			copyMPS(phi);
		}else if(site==Len)
		{
			MPS phi(Len+1,pD,1);
			phi.setZero();
			for(int i = 0; i < pD; ++i) phi.M[Len][i].setIdentity();
			for(int i = 0; i < Len; ++i)
			{
				for(int j = 0; j < pD; ++j)
				{
					phi.M[i][j] = M[i][j];
				}
			}
			phi.setDim();
			clear();
			copyMPS(phi);
		}else if(site>0&&site<Len)
		{
			setDim();
			MPS phi(Len+1,pD,bD);
			phi.setZero();
			for(int i = 0; i < pD; ++i) phi.M[site][i].setIdentity(M[site-1][0].cols(),M[site][0].rows());
			for(int i = 0; i < site; ++i)
			{
				for(int j = 0; j < pD; ++j)
				{
					phi.M[i][j] = M[i][j];
				}
			}
			for(int i = site; i < Len; ++i)
			{
				for(int j = 0; j < pD; ++j)
				{
					phi.M[i+1][j] = M[i][j];
				}
			}
			phi.setDim();
			clear();
			copyMPS(phi);
		}
	}
}

void MPS::setDim()
{
	for(int i = 0; i < Len; ++i)
	{
		Dim[i] = M[i][0].rows();
	}
	Dim[Len] = 1;
	bD = *std::max_element(Dim,Dim+Len+1);
}

void MPS::setZero()
{
	for(int i = 0; i < Len; ++i)
	{
		for(int j = 0; j < pD; ++j)
		{
			M[i][j].setZero(Dim[i],Dim[i+1]);
		}
	}
}

void MPS::setZero(int newD)
{
	bD = newD;
	Dim[0] = 1;
	for(int i = 1; i < Len; ++i)
	{
		int pw = std::min(i,Len-i);
		Dim[i] = std::pow(pD,pw);
		if(Dim[i]<0||log2(Dim[i])<pw)
		{
			Dim[i] = bD;
		}else
		{
			Dim[i] = std::min(Dim[i],bD);
		}
	}
	Dim[Len] = 1;
	for(int i = 0; i < Len; ++i)
	{
		for(int j = 0; j < pD; ++j)
		{
			M[i][j].setZero(Dim[i],Dim[i+1]);
		}
	}
}

void MPS::setRand()
{
	for(int i = 0; i < Len; ++i)
	{
		for(int j = 0; j < pD; ++j)
		{
			// std::cout<<i<<" "<<j<<" "<<Dim[i]<<" "<<Dim[i+1]<<std::endl;
			M[i][j].setRandom(Dim[i],Dim[i+1]);
		}
	}
}

void MPS::setRand(int newD)
{
	bD = newD;
	Dim[0] = 1;
	for(int i = 1; i < Len; ++i)
	{
		int pw = std::min(i,Len-i);
		Dim[i] = std::pow(pD,pw);
		if(Dim[i]<0||log2(Dim[i])<pw)
		{
			Dim[i] = bD;
		}else
		{
			Dim[i] = std::min(Dim[i],bD);
		}
	}
	Dim[Len] = 1;
	
	for(int i = 0; i < Len; ++i)
	{
		for(int j = 0; j < pD; ++j)
		{
			M[i][j].setRandom(Dim[i],Dim[i+1]);
		}
	}
}

void MPS::RC()
{
	if(!if_init)
	{
		std::cout<<"MPS not initiated!"<<std::endl;
		abort();
	}
	int tid;
	int row, col;
	Mxd TM, Q, R;
	for(int i = Len-1; i >= 0; i--)
	{
		row=M[i][0].rows();
		col=M[i][0].cols();
		TM.resize(row,pD*col);
		for(tid=0; tid<pD; tid++)
		//#pragma omp parallel num_threads(pD) private(tid)
		{
			//tid = omp_get_thread_num();
			TM.block(0,tid*col,row,col)=M[i][tid];
		}
		TM.transposeInPlace();
		rQR(TM,Q);
		Q.transposeInPlace();
		R.noalias() = TM.transpose()*Q.transpose();
		for(tid=0; tid<pD; tid++)
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
			}
		}
	}
}

void MPS::LC()
{
	if(!if_init)
	{
		std::cout<<"MPS not initiated!"<<std::endl;
		abort();
	}
	int tid;
	int row, col;
	Mxd TM, Q, R;
	for(int i = 0; i < Len; i++)
	{
		// std::cout<<i<<std::endl;
		row=M[i][0].rows();
		col=M[i][0].cols();
		// std::cout<<row<<" "<<col<<std::endl;
		TM.resize(pD*row,col);
		for(tid=0; tid<pD; tid++)
		//#pragma omp parallel num_threads(pD) private(tid)
		{
			//tid = omp_get_thread_num();
			TM.block(tid*row,0,row,col)=M[i][tid];
		}
		rQR(TM,Q);
		R.noalias() = Q.transpose()*TM;
		for(tid=0; tid<pD; tid++)
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
			}
		}
	}
}

void MPS::moveLeft(int site)
{
	if(!if_init)
	{
		std::cout<<"MPS not initiated!"<<std::endl;
		abort();
	}
	if(site==Len/2)
	{
		// Calculate EE at mid bond
		double * sv = new double [pD*bD]();
		int row, col, tid, DT, tDim;
		Mxd TM, U, V;	
		row=M[site][0].rows();
		col=M[site][0].cols();		
		TM.resize(row,pD*col);
		tDim = std::min(TM.rows(),TM.cols());
		for(tid = 0; tid < pD; tid++)
		{
			TM.block(0,tid*col,row,col) = M[site][tid].block(0,0,row,col);
		}
		rSVD(TM,pD*bD,sv,U,V,'l');
		double EE = 0;
		for(int j = 0; j < tDim; ++j)
		{
			EE -= sv[j]*sv[j]*log(sv[j]*sv[j]); 
		}
		// std::cout<<EE<<" ";
		for(tid = 0; tid < pD; tid++)
		{
			M[site][tid].setZero();
			M[site][tid].block(0,0,std::min(row,tDim),col) = V.block(0,tid*col,std::min(row,tDim),col);
			Mxd tempM;
			tempM.noalias() = M[site-1][tid] * U;
			M[site-1][tid].setZero();
			M[site-1][tid].block(0,0,tempM.rows(),tempM.cols())=tempM;
		}
		delete [] sv;
	}else
	{
		int tid;
		int row, col;
		Mxd TM, Q, R;
		row=M[site][0].rows();
		col=M[site][0].cols();
		TM.resize(row,pD*col);
		for(tid=0; tid<pD; tid++)
		//#pragma omp parallel num_threads(pD) private(tid)
		{
			//tid = omp_get_thread_num();
			TM.block(0,tid*col,row,col)=M[site][tid];
		}
		TM.transposeInPlace();
		rQR(TM,Q);
		Q.transposeInPlace();
		R.noalias() = TM.transpose()*Q.transpose();
		for(tid=0; tid<pD; tid++)
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
}

void MPS::moveRight(int site)
{
	if(!if_init)
	{
		std::cout<<"MPS not initiated!"<<std::endl;
		abort();
	}
	if(site==Len/2-1)
	{
		// Calculate EE at mid bond
		double * sv = new double [pD*bD]();
		int row, col, tid, DT, tDim;
		Mxd TM, U, V;	
		row=M[site][0].rows();
		col=M[site][0].cols();		
		TM.resize(pD*row,col);
		tDim = std::min(TM.rows(),TM.cols());
		for(tid = 0; tid < pD; tid++)
		{
			TM.block(tid*row,0,row,col) = M[site][tid].block(0,0,row,col);
		}
		rSVD(TM,pD*bD,sv,U,V,'r');
		double EE = 0;
		for(int j = 0; j < tDim; ++j)
		{
			EE -= sv[j]*sv[j]*log(sv[j]*sv[j]); 
		}
		// std::cout<<EE<<" ";
		for(tid = 0; tid < pD; tid++)
		{
			M[site][tid].setZero();
			M[site][tid].block(0,0,row,std::min(col,tDim)) = U.block(tid*row,0,row,std::min(col,tDim));
			Mxd tempM;
			tempM.noalias() = V * M[site+1][tid];
			M[site+1][tid].setZero();
			M[site+1][tid].block(0,0,tempM.rows(),tempM.cols())=tempM;
		}
		delete [] sv;
	}else
	{
		int tid;
		int row, col;
		Mxd TM, Q, R;
		row=M[site][0].rows();
		col=M[site][0].cols();
		TM.resize(pD*row,col);
		for(tid=0; tid<pD; tid++)
		//#pragma omp parallel num_threads(pD) private(tid)
		{
			//tid = omp_get_thread_num();
			TM.block(tid*row,0,row,col)=M[site][tid];
		}
		rQR(TM,Q);
		R.noalias() = Q.transpose()*TM;
		for(tid=0; tid<pD; tid++)
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
}

void MPS::readMPS(std::string filename)
{
	std::ifstream fin;
	fin.open(filename.c_str());
	if(fin.is_open())
	{
		// Load MPS information
		fin>>Len;
		fin>>pD;
		fin>>bD;
		setMPS(Len,pD,bD);
		for(int i = 0; i < Len; ++i)
		{
			fin>>Dim[i];
		}
		fin>>Dim[Len];
		setZero();
		double temp;
		std::vector<double> vec;
		for(int i = 0; i < Len; ++i)
		{
			for(int j = 0; j < pD; ++j)
			{
				vec.clear();
				for(int k = 0; k < Dim[i]*Dim[i+1]; ++k)
				{
					fin>>temp;
					vec.push_back(temp);
				}
				M[i][j] = Eigen::Map<Mxd>(&vec[0],Dim[i],Dim[i+1]);
			}
		}
		fin.close();
		std::cout<<"Reading MPS done!"<<std::endl;
	}else
	{
		std::cout<<"Reading MPS failed! No such file exists!"<<std::endl;
	}
}

void MPS::writeMPS(std::string filename, int preci)
{
	if(!if_init)
	{
		std::cout<<"MPS not written disk! Because it is not initiated!"<<std::endl;
	}
	std::ofstream fout;
	fout.precision(preci);
	fout.open(filename.c_str());
	// MPS information
	fout<<Len<<" "<<pD<<" "<<bD<<std::endl;
	for(int i = 0; i < Len; ++i)
	{
		fout<<Dim[i]<<" ";
	}
	fout<<Dim[Len]<<std::endl;
	// MPS matrices
	for(int i = 0; i < Len; ++i)
	{
		for(int j = 0; j < pD; ++j)
		{
			for(int k = 0; k < M[i][j].size(); ++k)
			{
				fout<<M[i][j].data()[k]<<" ";
			}
			fout<<std::endl;
		}
	}
	fout.close();
	std::cout<<"Writing MPS done!"<<std::endl;
}

void MPS::EE(bool verbose)
{
	if(!if_init)
	{
		std::cout<<"MPS not initiated!"<<std::endl;
		abort();
	}
	// std::cout<<"Entanglement Entropy:"<<std::endl;
	///////////////////////////////////////////
	RC();
	double * sv = new double [pD*bD]();
	int row, col, tid, DT, tDim;
	Mxd TM, U, V;
	///////////////////////////////////////////
	for(int i = 0; i < Len-1; i++)
	{
		row=M[i][0].rows();
		col=M[i][0].cols();
		
		TM.resize(pD*row,col);
		
		tDim = std::min(TM.rows(),TM.cols());
		
		for(tid = 0; tid < pD; tid++)
		{
			TM.block(tid*row,0,row,col) = M[i][tid].block(0,0,row,col);
		}
		
		rSVD(TM,tDim,sv,U,V,'r');
		
		double EEn = 0;
		for(int j = 0; j < tDim; ++j)
		{
			if(sv[j]>0) EEn -= sv[j]*sv[j]*log(sv[j]*sv[j]); 
		}
		
		if(i==(Len/2-1) && !verbose) std::cout<<EEn<<" ";
        if(verbose) std::cout<<EEn<<" ";
		
		for(tid = 0; tid < pD; tid++)
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
	delete [] sv;
}

double MPS::norm()
{
	if(!if_init)
	{
		std::cout<<"MPS not initiated!"<<std::endl;
		abort();
	}
	
	double nm = 0;
	int tid;
	int row, col;
	Mxd TM, Q, R;
	for(int i = Len-1; i >= 0; i--)
	{
		row=M[i][0].rows();
		col=M[i][0].cols();
		TM.resize(row,pD*col);
		for(tid=0; tid<pD; tid++)
		//#pragma omp parallel num_threads(pD) private(tid)
		{
			//tid = omp_get_thread_num();
			TM.block(0,tid*col,row,col)=M[i][tid];
		}
		TM.transposeInPlace();
		rQR(TM,Q);
		Q.transposeInPlace();
		R.noalias() = TM.transpose()*Q.transpose();
		for(tid=0; tid<pD; tid++)
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
				nm = R(0,0);
			}
		}
	}
	
	return std::abs(nm);
}


Mxd MPS::evaluateMPS(int N)
{
    int i = N%2;
    N /= 2;
    Mxd MT = M[0][i];
    for(int site=1; site<Len; site++)
    {
        i  = N%2;
        N /= 2;
        MT *= M[site][i];
    }
    return MT;
}


Mxd MPS::partial_trace(int pos, int l){
    if(!if_init)
	{
		std::cout<<"MPS not initiated!"<<std::endl;
		abort();
	}
    if(pos+l>Len)
    {
		std::cout<<"out of bounds error"<<std::endl;
        abort();
    }

    // matrix is of form M-M-M-M-M-M-M
    //                   | | | | | | |

    int N = std::pow(pD,l);
    Mxd rhoA(N,N);

    // LC to have        Q-Q-M-M-M-M-M
    //                   | | | | | | |
    for(int i=0; i<pos; i++)
        moveRight(i);

    // RC to have        Q-Q-M-M-Q-Q-Q
    //                   | | | | | | |
    for(int i=Len-1; i>=pos+l; i--)
        moveLeft(i);


    // psi_A is             -M-M-     
    //                       | |      
    MPS psi_A = MPS(l, pD, bD);
    for (int i=0; i<l; i++)
        for (int tid=0; tid<pD; tid++)
        {
            psi_A.M[i][tid] = M[pos+i][tid];
        }
        


    // loop over indices to get rho
    //                      -M-M-    
    //                    /  | |  \
    //                    |  i j  |     
    //                    |       |    
    //                    |  k l  |     
    //                    \  | |  /     
    //                      -M-M-     
    /*
    double n=0;
    for(int a=0; a<N; a++)
    {
        Mxd t = evaluateMPS(a);
        n += (t*(t.transpose())).trace();
    }
    std::cout<<"norm = "<<n<<std::endl;
    */

    double t;
    double tol = 1e-8/(N*N);
    for(int a=0; a<N; a++)
        for(int b=0; b<N; b++)
        {
            Mxd tmpA = psi_A.evaluateMPS(a);
            Mxd tmpB = psi_A.evaluateMPS(b);

            t = (tmpA*(tmpB.transpose())).trace();
            rhoA(a,b) = abs(t) < tol ? 0 : t;
        }
	
    return rhoA;
}

#endif
