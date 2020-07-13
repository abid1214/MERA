#ifndef My_UTILITY_FUNCTIONS
#define My_UTILITY_FUNCTIONS

#include "utility.h"

typedef Eigen::MatrixXd Mxd;

double MPS_partial_overlap(MPS psi, MPS phi, unsigned site, unsigned l)
{
    Mxd rhoA   = psi.partial_trace(site, l);
    Mxd sigmaA = phi.partial_trace(site, l);

    return (rhoA*sigmaA).trace();
}

// B = A1 cross A2
void KroneckerProd(Mxd& A1, Mxd& A2, Mxd& B)
{
	assert(A1.size()>0 && A2.size()>0);
	
	B.setZero(A1.rows()*A2.rows(), A1.cols()*A2.cols());
	
	int r1 = A1.rows();
	int c1 = A1.cols();

	int r2 = A2.rows();
	int c2 = A2.cols();
	
	for(int i = 0; i < r1; ++i)
	{
		for(int j = 0; j < c1; ++j)
		{
			B.block(i*r2, j*c2, r2, c2) = A1(i,j)*A2;
		}
	}
}

void effH(MPO& H, int site, int L, Mxd& Heff)
{
	int Len = H.Len;
	int pD  = H.pD;
	int bD  = H.bD;
	//////////////////////////
	Mxd lEnv, rEnv;
	// Left environment
	if(site==0)
	{
		lEnv.setIdentity(1,1);
	}else
	{
		lEnv.setZero(H.Dim[0],H.Dim[1]);
		for(int i = 0; i < pD; ++i)
		{
			lEnv += H.M[0][i*pD+i]/pD;
		}
		for(int i = 1; i < site; ++i)
		{
			Mxd tp(H.Dim[i],H.Dim[i+1]);
			tp.setZero();
			for(int j = 0; j < pD; ++j)
			{
				tp += H.M[i][j*pD+j];
			}
			Mxd ttp = lEnv * tp / pD;
			lEnv = ttp;
		}
	}
	// Right environment
	if(site+L==Len)
	{
		rEnv.setIdentity(1,1);
	}else
	{
		rEnv.setZero(H.Dim[Len-1],H.Dim[Len]);
		for(int i = 0; i < pD; ++i)
		{
			rEnv += H.M[Len-1][i*pD+i]/pD;
		}
		for(int i = Len-2; i > site+L-1; --i)
		{
			Mxd tp(H.Dim[i],H.Dim[i+1]);
			tp.setZero();
			for(int j = 0; j < pD; ++j)
			{
				tp += H.M[i][j*pD+j];
			}
			Mxd ttp = tp * rEnv / pD;
			rEnv = ttp;
		}
	}
	// Build the dense effective Hamiltnonian
	// The default ordering of the basis vectors:
	// 0: +, 1: -
	// 00...0, 10...0, ...
	int size = std::pow(pD,L);
	// std::cout<<"size "<<size<<std::endl;
	int phyr[L];
	int phyc[L];
	Heff.setZero(size,size);
	for(int row = 0; row < size; ++row)
	{
		for(int j = 0; j < L; ++j)
		{
			phyr[j] = int(row/std::pow(pD,j))%pD;
		}
		for(int col = 0; col < size; ++col)
		{
			for(int j = 0; j < L; ++j)
			{
				phyc[j] = int(col/std::pow(pD,j))%pD;
			}
			Mxd tp(H.Dim[site],H.Dim[site]);
			tp.setIdentity();
			for(int k = 0; k < L; ++k)
			{
				Mxd ttp = tp * H.M[k+site][phyr[k]*pD+phyc[k]];
				tp = ttp;
			}
			Heff(row,col) = (lEnv*tp*rEnv)(0,0);
		}
	}
}

// Apply A to psi
// A is thought of as the unitary gates
// A is applied through sites Site to Site+L(A)-1 of psi
// op = 'I' (kept the same), 'T' (tranposed) 
void applyMPO(MPO& A, MPS& psi, int site, char op)
{
	assert(A.if_init && psi.if_init);
	assert(site+A.Len-1<psi.Len && A.pD==psi.pD);
	assert(direc=='T'||direc=='B');
	assert(op=='I'||op=='T');
	
	int pD = A.pD;
	
	Mxd* tp = new Mxd [pD];
	Mxd temp;
	
	for(int i = 0; i < A.Len; ++i)
	{
		for(int j = 0; j < pD; ++j)
		{
			for(int k = 0; k < pD; ++k)
			{
				if(op=='I')
					KroneckerProd(A.M[i][j*pD+k], psi.M[site+i][k], temp);
				else if(op=='T')
					KroneckerProd(A.M[i][k*pD+j], psi.M[site+i][k], temp);
				
				if(k==0)
					tp[j] = temp;
				else
					tp[j] += temp;
			}
		}
		for(int j = 0; j < pD; ++j)
		{
			psi.M[site+i][j] = tp[j];
		}
	}
	
	psi.setDim();

	delete [] tp;
}

// Apply A to B
// A is thought of as the unitary gates
// A is applied through sites Site to Site+L(A)-1 of B
// from both the top and the bottom direction
// direc =T (top), B (bottom), op = 'I' (kept the same), 'T' (tranposed) 
void applyMPO(MPO& A, MPO& B, int site, char direc, char op)
{
	assert(A.if_init && B.if_init);
	assert(site+A.Len-1<B.Len && A.pD==B.pD);
	assert(direc=='T'||direc=='B');
	assert(op=='I'||op=='T');
	
	int pD = A.pD;
	
	Mxd* tp = new Mxd [pD*pD];
	Mxd temp;
	
	for(int i = 0; i < A.Len; ++i)
	{
		for(int j = 0; j < pD*pD; ++j)
		{
			int topIdx = j/pD; int botIdx = j%pD;
			for(int k = 0; k < pD; ++k)
			{
				if(direc=='T' && op=='I')
					KroneckerProd(A.M[i][topIdx*pD+k], B.M[site+i][k*pD+botIdx], temp);
				else if(direc=='B' && op=='I')
					KroneckerProd(B.M[site+i][topIdx*pD+k], A.M[i][k*pD+botIdx], temp);
				else if(direc=='T' && op=='T')
					KroneckerProd(A.M[i][k*pD+topIdx], B.M[site+i][k*pD+botIdx], temp);
				else if(direc=='B' && op=='T')
					KroneckerProd(B.M[site+i][topIdx*pD+k], A.M[i][botIdx*pD+k], temp);
				
				if(k==0)
					tp[j] = temp;
				else
					tp[j] += temp;
			}
		}
		for(int j = 0; j < pD*pD; ++j)
		{
			B.M[site+i][j] = tp[j];
		}
	}
	
	for(int i = 0; i < A.Len+1; ++i)
	{
		B.Dim[site+i] *= A.Dim[i];
	}
	
	B.bD = *std::max_element(B.Dim,B.Dim+B.Len+1);

	delete [] tp;
}

void applyGates(Mxd& G, MPO& H, int site, int L, int max_BD)
{	
	assert(G.rows()==G.cols());
	assert(G.rows()==std::pow(H.pD,L));
	/////////////////////////////
	// Build Gates as MPO
	MPO U;
	EDtoMPO(G, L, H.pD, U);
	U.examineBond(1E-6,false);
	/////////////////////////////
	// Apply Gates -- direct product
	// (U^T*H*U)
	applyMPO(U, H, site, 'T', 'T'); // apply to the top
	if(H.bD>max_BD) iterCompress(true, 4, 1e-14, max_BD, H, false);
	applyMPO(U, H, site, 'B', 'I'); // apply to the bottom
	if(H.bD>max_BD) iterCompress(true, 4, 1e-14, max_BD, H, false);
	/////////////////////////////
	// Compress if necessary
	// if(H.bD>max_BD) iterCompress(true, 4, 1e-14, max_BD, H);
	/////////////////////////////
}

void applyGates(Mxd& G, MPO& H, int site, int L)
{	
	assert(G.rows()==G.cols());
	assert(G.rows()==std::pow(H.pD,L));
	/////////////////////////////
	// Build Gates as MPO
	MPO U;
	EDtoMPO(G, L, H.pD, U);
	/////////////////////////////
	// Apply Gates -- direct product
	// (U^T*H*U)
	applyMPO(U, H, site, 'T', 'T'); // apply to the top
	applyMPO(U, H, site, 'B', 'I'); // apply to the bottom
	/////////////////////////////
}

double var(MPO& H)
{
	// std::cout<<"Calculating variance..."<<std::endl;
  MPO H1;
  H1.copyMPO(H);
  H1.keepDiag();
  H1.norm = 1;
  H1.RC();
  double v1 = H1.norm*H1.norm;
  
  MPO H2;
  H2.copyMPO(H);
  H2.norm = 1;
  H2.RC();
  // H2.square();
  double v2 = H2.norm*H2.norm;
  
  return v2 - v1;
}

// Find Permutation -- to speed up Wegner Diagonalization processes
bool findPermutation(Mxd& A)
{
	// std::cout<<"Look for permutation"<<std::endl;
	// int temp;
	for(int i = 0; i < A.rows(); ++i)
	{
		int idx = -1;
		double ov = 0;
		for(int j = 0; j < A.cols(); ++j)
		{
			if(ov<std::abs(A(i,j)))
			{
				idx = j;
				ov = std::abs(A(i,j));
			}
		}
		for(int j = 0; j < A.cols(); ++j)
		{
			if(j==idx)
			{
				if(A(i,j)>0)
					A(i,j)=1;
				else
					A(i,j)=-1;
			}else
			{
				A(i,j)=0;
			}
		}
	}
	// std::cout<<A<<std::endl;
	// for(int i = 0; i < A.size(); ++i)
	// {
	// 	A(i) = std::round(A(i));
	// }
	Mxd idm(A.rows(),A.cols());
	idm.setIdentity();
	Mxd tp = A*A.transpose()-idm;
	// if( abs((A*A.transpose()-idm).squaredNorm()) > 0.001 ) cout<<"Permutation not found!"<<endl;
	if( tp.norm() > 0.001 )
	{
		// std::cout<<A<<std::endl;
		return false;
	}else
	{
		return true;
	}
}

// Returns true if the matrix contains any NaN entries in it
bool findNaN(Mxd &H)
{
    bool res=false;
    for (int k=0; k<H.size(); ++k){
        if(std::isnan(H(k))) res=true;
    }
    return res;
};

 // Symmetrizes a matrix 
void mkSymmetric(Mxd& h1){
    Mxd h2(h1.rows(),h1.cols());
    
    h2=h1.transpose();
    h1+=h2;
    h1/=2.0;
}

// Find the max element of the matrix -- including or excluding the diagonal
double findMax(Mxd H, bool withDiagonal)
{
	double currMax = 0;
	for(int i = 0; i < H.rows(); ++i)
	{
		if(withDiagonal)
		{
			for(int j = i; j < H.cols(); ++j)
			{
				if(std::abs(H(i,j))>currMax) currMax = std::abs(H(i,j));
			}
		}else
		{
			for(int j = i+1; j < H.cols(); ++j)
			{
				if(std::abs(H(i,j))>currMax) currMax = std::abs(H(i,j));
			}
		}

	}
	return currMax;
}

// Check whether a matrix is close to identity --  assuming it is a square matrix
bool is_indentity(Mxd& U)
{
  if(U.rows()!=U.cols())
    {
      std::cout<<"Not square matrix"<<std::endl;
      return false;
    }

  int r = U.rows();
  int c = U.cols();
  Mxd tp_id(r,c);
  tp_id.setIdentity();
  tp_id = tp_id - U;
  if(tp_id.norm()<1E-9)
    return true;
  else
    return false;
}

// void writeToFile(std::string fn, MPS& psi);
//
// void readFromFile(std::string fn, MPS& psi);
//
// void writeToFile(std::string fn, MPO& H);
//
// void readFromFile(std::string fn, MPO& H);


#endif
