#ifndef My_SIMPLE_WEGNER_FLOW
#define My_SIMPLE_WEGNER_FLOW

#include "simpleWegnerFlow.h"
#include <unsupported/Eigen/MatrixFunctions>

typedef Eigen::MatrixXd Mxd;

void SimpleWegner::setTau(double t)
{
	tau = t;
}

void SimpleWegner::setTol(double t)
{
	tol = t;
}

void SimpleWegner::setIter(int t)
{
	iter = t;
}

void SimpleWegner::setH(Mxd& tH)
{
	H = tH;
	A = tH;
	U.setIdentity(H.rows(), H.cols());
}

void SimpleWegner::diag(Mxd& tU)
{
	int step = 0;
	
	Mxd diagH, eta, tpU, tH;
	
	while(step<999999)
	{
		diagH  = (H.diagonal()).asDiagonal();
		eta    = -tau*(diagH*H - H*diagH);
		tpU    = eta.exp();
		tU     = U * tpU;
		U      = tU;
		tH     = tpU.transpose() * H * tpU;
		H      = tH;

		mkSymmetric(H);
		double offDiagMax = findMax(H,false);
		// std::cout<<"Step="<<step<<"; off-diag="<<offDiagMax<<std::endl;
		if(offDiagMax < tol) break;
		++step;
	}
	// for(int i = 0; i < H.size(); ++i)
	// {
	// 	if( std::abs(H(i))<tol ) H(i) = 0;
	// 	if( std::abs(U(i))<tol ) U(i) = 0;
	// }
	// std::cout<<H<<std::endl<<std::endl;
	
    Eigen::SelfAdjointEigenSolver<Mxd> es(A);
    if (es.info() != Eigen::Success) abort();
	Mxd D = es.eigenvalues().asDiagonal();
	Mxd V = es.eigenvectors();
	Mxd P = U.transpose() * V;

	if(!findPermutation(P)) std::cout<<"Permutation not found!\n";
	tU.noalias() = V*P.transpose();
}

#endif
