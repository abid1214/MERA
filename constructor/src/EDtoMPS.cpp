#ifndef MY_ED_TO_MPS_FUNCTION
#define MY_ED_TO_MPS_FUNCTION

#include "EDtoMPS.h"

typedef Eigen::MatrixXd Mxd;

void EDtoMPS(double * vec, int Len, int pD, MPS& psi)
{
	assert(pow(pD,Len/2)<2048);
	int size = pow(pD,Len);
	int phy[Len];
	double* v = new double [size]();
	for(int i = 0; i < size; i++)
	{
		for(int j = 0; j < Len; ++j)
		{
			phy[j] = int(i/std::pow(pD,j))%pD;
		}
		int idx = 0;
		for(int j = 0; j < Len; ++j)
		{
			idx += std::pow(pD,Len-j-1)*(phy[j]);
		}
		v[idx] = vec[i];
		// v[i] = vec[size-i-1];
		// v[i] = vec[i];
	}
	psi.setMPS(Len, pD, pow(pD,Len/2));
	psi.setZero();
	for(int i = 0; i < Len; i++)
	{
		int row, col, s1, s2;
		s1  = pow(pD, Len-i);
		s2  = pow(pD, Len-i-1);
		row = psi.Dim[i];
		col = psi.Dim[i+1];
		Mxd T, Q, R;
		T.setZero(pD*row,s2);
		for(int r = 0; r < row; r++)
		{
			for(int p = 0; p < pD; p++)
			{
				for(int c = 0; c < s2; c++)
				{
					T(p*row+r,c) = v[r*s1+p*s2+c];
				}
			}
		}
		rQR(T,Q);
		R.noalias() = Q.transpose()*T;
		for(int j = 0; j < pD; j++)
		{
			psi.M[i][j] = Q.block(j*row,0,row,col);
		}
		R.transposeInPlace();
		std::copy(R.data(),R.data()+R.size(),v);
	}
	delete [] v;
}


#endif