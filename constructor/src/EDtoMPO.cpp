#ifndef MY_ED_TO_MPO_FUNCTION
#define MY_ED_TO_MPO_FUNCTION

#include "EDtoMPO.h"
typedef Eigen::MatrixXd Mxd;

void EDtoMPO(Mxd& A, int Len, int pD, MPO& H)
{
	assert(std::pow(pD*pD,Len/2)<4096);
	int size = std::pow(pD,Len);
	double* v = new double [size*size]();
	///////////////////////////
	// Need to rewrite the Matrix as an ED vector
	int phyr[Len];
	int phyc[Len];
	int idx;
	for(int row = 0; row < size; ++row)
	{
		for(int j = 0; j < Len; ++j)
		{
			phyr[j] = int(row/std::pow(pD,j))%pD;
		}
		for(int col = 0; col < size; ++col)
		{
			for(int j = 0; j < Len; ++j)
			{
				phyc[j] = int(col/std::pow(pD,j))%pD;
			}
			idx = 0;
			for(int i = 0; i < Len; ++i)
			{
				idx += std::pow(pD*pD,Len-i-1)*(phyr[i]*pD+phyc[i]);
			}
			v[idx] = A(row,col);
		}
	}
	///////////////////////////
	H.setMPO(Len, pD, pow(2*pD,int(Len/2)),0);
	H.setZero();
	Mxd T, Q, R;
	for(int i = 0; i < Len; i++)
	{
		int row, col, s1, s2;
		s1  = pow(pD*pD, Len-i);
		s2  = pow(pD*pD, Len-i-1);
		row = H.Dim[i];
		col = H.Dim[i+1];
		T.setZero(pD*pD*row,s2);
		for(int r = 0; r < row; r++)
		{
			for(int p = 0; p < pD*pD; p++)
			{
				for(int c = 0; c < s2; c++)
				{
					T(p*row+r,c) = v[r*s1+p*s2+c];
				}
			}
		}
		rQR(T,Q);
		R.noalias() = Q.transpose()*T;
		for(int j = 0; j < pD*pD; j++)
		{
			H.M[i][j] = Q.block(j*row,0,row,col);
		}
		R.transposeInPlace();
		std::copy(R.data(),R.data()+R.size(),v);
	}
	// std::cout<<R(0,0)<<std::endl;
	for(int i = 0; i < pD*pD; ++i)
	{
		H.M[0][i] *= R(0,0);
	}
	delete [] v;
}


#endif