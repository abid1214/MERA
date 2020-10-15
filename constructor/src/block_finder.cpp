#ifndef BLOCK_FINDER
#define BLOCK_FINDER

#include "block_finder.h"
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include "utility.h"
#include "mpo.h"
#include <cmath>
#include "simpleWegnerFlow.h"
#include <time.h>
#include <stdlib.h>
#include "wegnerFlow.h"

typedef Eigen::MatrixXd Mxd;
using namespace std;

int random_block(int L, int l)
{
    srand ((time(NULL)));
    return rand() % (L-l-1);
}

int findMaxGap(MPO& H, int l)
{
	double max_gap = 0;
    double gap;
    int max_gap_site = 0;
    int L = H.Len;
    if (l >= L)
        return 0;

	Mxd A;
    for(int i = 0; i < L-(l-1); i++)
    {
        effH(H,i,l,A);
        
        Eigen::SelfAdjointEigenSolver<Mxd> es(A);
        if (es.info() != Eigen::Success) abort();
        gap = es.eigenvalues()(A.rows()-1)-es.eigenvalues()(0);
        
        if(max_gap < gap)
        {
            max_gap_site = i;
            max_gap      = gap;
        }
	}
    return max_gap_site;
}

int firstBlock()
{
    return 0;
}

std::vector<int> int_to_bin(int n, int l, int pD)
{
    std::vector<int> sv(l);
    for(int i = 0; i < l; i++)
    {
        sv[l-1-i] = n%pD;
        n /= pD;
    }
    return sv;
}

int bin_to_int(vector<int> sv, int l, int pD)
{
    int n = 0;
    for(int i = 0; i < l; i++)
        n += sv[l-i-1]*std::pow(pD, i);
    return n;
}


Mxd U_seperator(Mxd& U, int pD, int l)
{
    //split up U
    int rows = std::pow(pD,2*(l-2));
    int cols = std::pow(pD,4);
    Mxd Umat(rows, cols);
    int R, C;
    for(int i = 0; i < U.rows(); i++)
        for(int j = 0; j < U.cols(); j++)
        {
            std::vector<int> bin_i = int_to_bin(i, l, pD);
            std::vector<int> bin_j = int_to_bin(j, l, pD);
            int phya1 = bin_i[0];
            int phya2 = bin_i[l-1];
            int phya3 = bin_j[0];
            int phya4 = bin_j[l-1];
            int t1[] = {phya1, phya2, phya3, phya4};
            int t2[2*(l-2)];
            for(int k = 0; k < l-2; k++)
            {
                t2[k]     = bin_i[k+1];
                t2[l-2+k] = bin_j[k+1];
            }
            std::vector<int> rv (t2, t2 + sizeof(t2) / sizeof(int) );
            std::vector<int> cv (t1, t1 + sizeof(t1) / sizeof(int) );
            C = bin_to_int(cv, 4, pD);
            R = bin_to_int(rv, 2*(l-2), pD);

            Umat(R,C) = U(i, j);
        }
    return Umat;
}

std::vector<double> get_sv(Mxd& Umat)
{
    Eigen::JacobiSVD<Mxd> svd(Umat);
    Mxd svmat = svd.singularValues();
    int n = svmat.rows();
    std::vector<double> sv(n);
    for(int i = 0; i < n; i++)
        sv[i] = svmat(i,0);
    return sv;
}

double get_E(std::vector<double>& sv)
{
    double E = 0;
    double Z = 0;
    for(int i = 0; i < sv.size(); i++)
        Z += sv[i]*sv[i];
    for(int i = 0; i < sv.size(); i++)
        E -= sv[i]*sv[i]/Z*std::log(sv[i]*sv[i]/Z);
    return E;
}

int findMinEntanglement(MPO& H, int l)
{
    int pD = H.pD;
    double min_E = 9999;
    double E;
    int min_E_site = 0;
    int L = H.Len;
    if (l >=L || l <= 2)
        return 0;

	Mxd A, U;
    for(int site = 0; site < L-(l-1); site++)
    {
        effH(H,site,l,A);
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

        Mxd Umat = U_seperator(U, pD, l);
        vector<double> sv = get_sv(Umat);
        E = get_E(sv);
        
        if(E < min_E)
        {
            min_E_site = site;
            min_E      = E;
        }
	}
    return min_E_site;
}

#endif
