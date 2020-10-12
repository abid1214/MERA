#ifndef TESTS
#define TESTS

#include <iostream>
#include "block_finder.h"
#include <Eigen/Dense>
#include <cmath>

typedef Eigen::MatrixXd Mxd;
using namespace std;



int test_int_bin(){
    int n  = 21;
    int l  = 6; 
    int pD = 2;

    cout<<n<<endl;
    std::vector<int> b = int_to_bin(n,l,pD);
    for(int i = 0; i < b.size(); i++)
        cout<<b[i]<<" ";
    cout<<endl;
    int np = bin_to_int(b, l, pD);
    cout<<np<<endl;
    return 0;
}

int test_U_seperator(){
    /*expected output: 
     
     0  1  2  3  4  5  6  7
     8  9 10 11 12 13 14 15
    16 17 18 19 20 21 22 23
    24 25 26 27 28 29 30 31
    32 33 34 35 36 37 38 39
    40 41 42 43 44 45 46 47
    48 49 50 51 52 53 54 55
    56 57 58 59 60 61 62 63

     0  1  4  5  8  9 12 13 32 33 36 37 40 41 44 45
     2  3  6  7 10 11 14 15 34 35 38 39 42 43 46 47
    16 17 20 21 24 25 28 29 48 49 52 53 56 57 60 61
    18 19 22 23 26 27 30 31 50 51 54 55 58 59 62 63

    290.643 29.5072 2.34743e-14 3.44082e-15

    0.307522

    */

    int l = 3;
    int pD = 2;
    int n = std::pow(pD,l);
    Mxd U(n,n);
    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
            U(i,j) = n*i + j;
    cout<<U<<endl; 
    cout<<endl;

    Mxd Umat = U_seperator(U, pD, l);
    cout<<Umat<<endl;

    vector<double> sv = get_sv(Umat);
    cout<<endl;
    for(int i = 0; i < sv.size(); i++)
        cout<<sv[i]<<" ";
    cout<<endl;
    cout<<endl;
    double E = get_E(sv);
    cout<<E<<endl;

    return 0;
}

int run_tests(){
    //test_int_bin();
    test_U_seperator();
    return 0;
}


#endif

