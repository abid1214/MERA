#ifndef DIAGNOSTICS_H
#define DIAGNOSTICS_H

#include "mps.h"
#include "mpo.h"
#include "observables.h"
#include "sdMERA.h"

using namespace std;

MPS load_mps(string data_dir, double W, int L, int l, unsigned seed, double epsilon)
{
    char fname [150];
    sprintf(fname, "%smps_W_%2.4f_L_%d_l_%d_d_%d_e_%0.2f.txt", data_dir.c_str(), W, L, l, seed, epsilon);
    MPS psi;
    psi.readMPS(fname);
    return psi;
}

vector<double> setDisorderConfig(int seed, int L, double WJ, double Wh)
{
	srand48(137*seed);

    vector<double> dJ;
    vector<double> dh;
    
	dJ.clear();
	dh.clear();
    for(int i = 0; i < L; ++i)
    {
        dJ.push_back(2*WJ*(drand48()-0.5));
        dh.push_back(2*Wh*(drand48()-0.5));
    }
    return dh;
}

double get_energy(MPS psi, double W, int seed)
{
    bool uniform = false;
	int pD = psi.pD;
	int bD = 5;
	double WJ = 0;
	double Wh = W;
    int L = psi.Len;

    vector<double> dh = setDisorderConfig(seed, psi.Len, WJ, Wh);
    vector<double> dJ(psi.Len, 0.0);

    MPO h;
    h.clearMPO();
    h.setMPO(L, pD, bD, 0);
    h.buildHeisenberg(&(dJ[0]), &(dh[0]));

    return psiHphi(psi,h,psi);
}

double get_variance(MPS psi, double W, int seed)
{

    bool uniform = false;
	int pD = psi.pD;
	int bD = 5;
	double WJ = 0;
	double Wh = W;
    int L = psi.Len;

    vector<double> dh = setDisorderConfig(seed, psi.Len, WJ, Wh);
    vector<double> dJ(psi.Len, 0.0);
    
    MPO h;
    h.clearMPO();
    h.setMPO(L, pD, bD, 0);
    h.buildHeisenberg(&(dJ[0]), &(dh[0]));

    MPO hS;
    hS.clearMPO();
    hS.setMPO(L, pD, bD, 0);
    hS.buildHeisenberg(&(dJ[0]), &(dh[0]));
    hS.square();

    double E = psiHphi(psi,h,psi);
    double SE = psiHphi(psi,hS,psi);

    return SE - E*E;
}

double bipartite_spin_fluctuation(MPS psi)
{
    MPO Sza;
    Sza.setMPO(psi.Len, psi.pD, 2, 0);
    Sza.buildSza(psi.Len/2);

    /*
    Mxd s;
    effH(Sza, 0, psi.Len, s);
    cout<<s<<endl;
    */

    MPO Sza_2;
    Sza_2.setMPO(psi.Len, psi.pD, 2, 0);
    Sza_2.buildSza(psi.Len/2);
    Sza_2.square();

    double M  = psiHphi(psi, Sza  , psi);
    double M2 = psiHphi(psi, Sza_2, psi);

    return M2 - M*M;
}

double participation_entropy(MPS psi, int q)
{
    int pD = psi.pD;
    int L  = psi.Len;
    double N  = pow(pD, L);

    double a;
    double PE = 0;
    if(q==1)
        for(int i = 0; i < N; i++)
        {
            a = (psi.evaluateMPS(i)).trace();
            if ((a*a) > 0) PE -= a*a*log(a*a);
        }
    else
    {
        for(int i = 0; i < N; i++)
        {
            a = (psi.evaluateMPS(i)).trace();
            if ((a*a) > 0) PE += pow(a*a, q);
        }
        PE = 1./(1-q)*log(PE);
    }
    return PE;
}

void run_diagnostics(string data_dir, double W, int L, int l, unsigned seed, double epsilon)
{
    MPS psi = load_mps(data_dir, W, L, l, seed, epsilon);
    char fname [100];
    sprintf(fname, "diagnostics_data/diag_W_%2.4f_L_%d_l_%d_d_%d_e_%0.2f.txt", W, L, l, seed, epsilon);

    ofstream fout;
    fout.precision(8);
    fout.open(fname);

    fout<<"Entanglements: ";
    vector<double> EE_list = psi.EE(false);
    for(int i = 0; i < psi.Len-1; i++)
        fout<<EE_list[i]<<" ";
    fout<<endl;

    fout<<"W: ";
    vector<double> Wlist = setDisorderConfig(seed, L, 0, W);
    for(int i = 0; i < L; i++)
        fout<<Wlist[i]<<" ";
    fout<<endl;

    /*
    fout<<"psi: "<<endl;
    for(int i = 0; i < pow(psi.pD, psi.Len); i++)
        fout<<psi.evaluateMPS(i)<<endl;
    fout<<endl;
    */
    
    //fout<<"Participation Entropy q=1:  "<<participation_entropy(psi, 1)<<endl;
    //fout<<"Participation Entropy q=2:  "<<participation_entropy(psi, 2)<<endl;
    fout<<"energy: "<<get_energy(psi, W, seed)<<endl;

    fout<<"variance: "<<get_variance(psi, W, seed)<<endl;

    fout<<"bipartite spin fluctuation: "<<bipartite_spin_fluctuation(psi)<<endl;
}
        
    
#endif

