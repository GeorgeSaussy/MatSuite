#include "mslib.h"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <vector>
#include <string>
using namespace std;
struct SqMat randomUnitValElmSqMat(int N) {
    struct SqMat toret;
    if(N>0) {
        initZeroSqMat(&toret, N);
        double theta=0.0;
        struct Complex matVal;
        for(int k=0;k<N;k++) {
            for(int k1=0;k1<N;k1++) {
                theta=((double) rand()/(RAND_MAX))*2.0*M_PI;
                matVal.re=cos(theta);
                matVal.im=sin(theta);
                setVal(&toret,k,k1,matVal);
            }
        }
    }
    return toret;
}
void printFullSqMat(struct SqMat mat) {
    int N=mat.N;
    for(int k=0;k<N;k++) {
        for(int k1=0;k1<N;k1++) {
            cout<<getVal(mat,k,k1).re<<"+i"<<getVal(mat,k,k1).im<<"\t";
        }
        cout<<endl;
    }
}
vector<struct SqMat> genRandUnitValElmSqMatList(int m, int N) {
    vector<struct SqMat> toret;
    if(m>0 && N>0) {
        struct SqMat toadd;
        for(int k=0;k<m;k++) {
            toadd=randomUnitValElmSqMat(N);
            toret.push_back(toadd);
        }
    }
    return toret;
}
