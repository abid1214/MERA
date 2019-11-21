#ifndef WegnerFLow_SOLVER
#define WegnerFLow_SOLVER

#include "wegnerFlow.h"

typedef Eigen::MatrixXd Mxd;


void WegnerDiagonalize::setH(Mxd &myH) {
	dH=myH;
    H=myH;
}

void WegnerDiagonalize::setTolRK(double t_tolleranceRK) {
    tolleranceRK=t_tolleranceRK;
}

void WegnerDiagonalize::setTol(double t_tolleranceElem){
    tolleranceElem=t_tolleranceElem;
	tolleranceRK  =t_tolleranceElem;
}

void WegnerDiagonalize::setTau(double t_deltaTau) {
    deltaTau=t_deltaTau;
}

// computes F' using HRun and URun
void WegnerDiagonalize::computeFP(){
    Mxd diag = HRun.diagonal().asDiagonal();
    Mxd eta = HRun*diag-diag*HRun;
    DeltaH=-eta*HRun+HRun*eta;
    DeltaU=URun*eta;
}

// performs a RK45 step (as defined in wikepedia)
void WegnerDiagonalize::RKstep()
{
    double bt1[7][7] = {{0.,2,3,4,5,6,7},
                        {1./5., 1./5.,3,4,5,6,7},
                        {3./10., 3./40., 9./40., 4,5,6,7},
                        {4./5., 44./45., -56./15., 32./9., 5,6,7},
                        {8./9., 19372./6561., -25360./2187., 64448./6561., -212./729.,6,7},
                        {1., 9017./3168., -355./33., 46732./5247., 49./176., -5103./18656.,7},
                        {1., 35./384., 0., 500./1113., 125./192., -2187./6784., 11./84.}};
    double bt2a[7] = {35./384., 0., 500./1113., 125./192., -2187./6784., 11./84., 0.};
    double bt2b[7] = {5179./57600., 0., 7571./16695., 393./640., -92097./339200., 187./2100., 1./40.};
    
    std::vector< Mxd > ktH;
    std::vector< Mxd > ktU;
    
    for(int i=0; i<7; i++){
        tauRun=tau+deltaTau*bt1[i][0];
        HRun = H;
        URun = U;
        for(int k=1; k<=i; k++){
            HRun += deltaTau * bt1[i][k] * ktH[k-1];
            URun += deltaTau * bt1[i][k] * ktU[k-1];
        }
        computeFP();
        ktH.push_back(DeltaH);
        ktU.push_back(DeltaU);
    }
    
    DeltaH5=bt2a[0]*ktH[0];
    DeltaU5=bt2a[0]*ktU[0];

    DeltaH4=bt2b[0]*ktH[0];
    DeltaU4=bt2b[0]*ktU[0];

    for(int k=1; k<7; k++){
        DeltaH5+=bt2a[k]*ktH[k];
        DeltaU5+=bt2a[k]*ktU[k];
        
        DeltaH4+=bt2b[k]*ktH[k];
        DeltaU4+=bt2b[k]*ktU[k];
    }
	
	ktH.clear();
	ktU.clear();
}

// runs the diagonalization using RK45 steps with adjustable step sizes
bool WegnerDiagonalize::diag(Mxd& _U) {
    Mxd matrixDeltaH;
    Mxd matrixDeltaU;

    int maxRow;
    int maxCol;
    double newtau;
    double maxDelta;
    double offDiagMax=99999;
    
    tau=0.0;
    steps=0;
    U.setIdentity(H.rows(),H.cols());

    // std::cout<<"Starting to diagonalize, matrix size="<<H.rows()<<std::endl;
    
    while ((offDiagMax>0)&&(steps<99999))
	{
        RKstep();

        matrixDeltaH=DeltaH5-DeltaH4;
        matrixDeltaU=DeltaU5-DeltaU4;
		
		// std::cout<<matrixDeltaH<<std::endl;
		// std::cout<<matrixDeltaH<<std::endl;

        maxDelta=std::max(findMax(matrixDeltaH,true),findMax(matrixDeltaU,true));
        
		// std::cout<<"Step="<<steps<<"  maxDelta="<<maxDelta<<" off-diag="<<offDiagMax<<std::endl;
        if(findNaN(DeltaH4)||findNaN(DeltaH5)||findNaN(DeltaU4)||findNaN(DeltaU5))
		{
            deltaTau=0.1;
        }else
		{
            if(maxDelta<tolleranceRK)
			{
                tau+=deltaTau;
            
                H+=deltaTau*DeltaH5;
                U+=deltaTau*DeltaU5;
				
				// trim elements
				for(int i = 0; i < H.size(); ++i)
				{
					if( std::abs(H(i))<tolleranceElem ) H(i) = 0;
					if( std::abs(U(i))<tolleranceElem ) U(i) = 0;
				}

                mkSymmetric(H);
				
                offDiagMax=findMax(H,false);

                if(0.2*tolleranceRK/maxDelta<1.5) {
                    deltaTau*=0.2*tolleranceRK/maxDelta;
                } else {
                    deltaTau*=1.5;
                }
            } else {
                deltaTau/=2.;
            }
        }
        steps++;
        if((std::isnan(deltaTau))||(deltaTau<=0.)) deltaTau=0.1;
		
		if(offDiagMax<tolleranceElem) break;
    }
	// std::cout<<H.diagonal()<<std::endl;
	//std::cout<<H<<std::endl;
	
	// std::cout<<"Diagonalization"<<std::endl;
    Eigen::SelfAdjointEigenSolver<Mxd> es(dH);
    if (es.info() != Eigen::Success) abort();
	Mxd D = es.eigenvalues().asDiagonal();
	Mxd V = es.eigenvectors();
	Mxd P = U.transpose() * V;

	if(findPermutation(P))
	{
		_U.noalias() = V*P.transpose();
		return true;
	}
	return false;
}


#endif
