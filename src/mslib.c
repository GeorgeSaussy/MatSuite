#include "mslib.h"
#include "stdlib.h"
#include "math.h"

void ms_HamXVec(Hamiltonian H, Cmplx * v, Cmplx * w) {
    Cell * iter=H.ptr;

    while((unsigned int)(iter-H.ptr)<H.len) {
        (w[iter->row]).re+=(iter->val).re*(v[iter->col]).re-(iter->val).im*(v[iter->col]).im;
        (w[iter->row]).im+=(iter->val).re*(v[iter->col]).im+(iter->val).re*(v[iter->col]).im;

        iter++;
    }
}

Hamiltonian *  ms_generateRandomHamiltonian(int sp, int dim) { // this is not fast
    Hamiltonian * ret=(Hamiltonian*)malloc(sizeof(Hamiltonian));
    // first make dim / sp groupings of sp rows, done with Knuths permtation shuffle
    int numberOfGroups=dim/sp;
    int groups[numberOfGroups][sp];
    int basearray[dim];
    int k, k1, k2, j, secondval; // for ANSI compliance
    for(k=0;k<dim;k++) {
        basearray[k]=k;
    }
    for(k=0;k<dim;k++) {
        j=rand()%(dim-k);
        secondval=basearray[k+j];
        basearray[k+j]=basearray[k];
        basearray[k]=secondval;
    }
    j=0;
    for(k=0;k<numberOfGroups;k++) {
        for(k1=0;k1<sp;k1++) {
            groups[k][k1]=basearray[j];
            j++;
        }
    }
    // generate dim / sp Hamiltonians
    Cmplx ham[numberOfGroups][sp][sp];
    Cmplx vec[sp];
    Cmplx one, zero;
    double norm=0;
    one.re=1;
    one.im=0;
    zero.re=0;
    zero.im=0;
    for(k=0;k<numberOfGroups;k++) {
        for(k1=0;k1<sp;k1++) {
            vec[k1].re=((double)rand())/((double)RAND_MAX);
            vec[k1].im=((double)rand())/((double)RAND_MAX);
        }
        norm=0;
        for(k1=0;k1<sp;k1++) {
            norm+=sqrt(vec[k1].re*vec[k1].re+vec[k1].im*vec[k1].im);
        }
        norm=sqrt(norm);
        for(k1=0;k1<sp;k1++) {
            vec[k1].re/=norm;
            vec[k1].im/=norm;
        }
        for(k1=0;k1<sp;k1++) {
            for(k2=0;k2<sp;k2++) {
                if(k1==k2) {
                    ham[k][k2][k1].re=1-2*vec[k2].re*vec[k1].re;
                    ham[k][k2][k1].re=-2*vec[k2].im*vec[k1].im;
                }
                else {
                    ham[k][k2][k1].re=-2*vec[k2].re*vec[k1].re;
                    ham[k][k2][k1].re=-2*vec[k2].im*vec[k1].im;
                }
            }
        }
    }
    // apply them to the row groups
    ret->N=dim;
    ret->len=sp*sp*numberOfGroups+(dim-sp*numberOfGroups);
    ret->ptr=(Cell*)malloc(ret->len*sizeof(Cell));
    j=0;
    for(k=0;k<numberOfGroups;k++) {
        for(k1=0;k1<sp;k1++) {
            for(k2=0;k2<sp;k2++) {
                ((ret->ptr)[j]).row=groups[k][k1];
                ((ret->ptr)[j]).col=groups[k][k2];
                ((ret->ptr)[j]).val=ham[k][k1][k2];
                j++;
            }
        }
    }
}

void ms_GenKrylovBasis(Hamiltonian H, Cmplx * v, Cmplx ** w, unsigned int m) {
    unsigned int count=1;

    ms_HamXVec(H, v, w[0]);
    while(count<m) {
        ms_HamXVec(H, w[count-1], w[count]);
        count++;
    }
}


void ms_parallel_HamXVec(Hamiltonian H, Cmplx ** v, Cmplx ** w, int num) {
    Cell * iter=H.ptr;
    int count=0;

    while((unsigned int)(iter-H.ptr)<H.len) {
        count=0;
        while(count<num) {
            (w[count][iter->row]).re+=(iter->val).re*(v[count][iter->col]).re-(iter->val).im*(v[count][iter->col]).im;
            (w[count][iter->row]).im+=(iter->val).re*(v[count][iter->col]).im+(iter->val).re*(v[count][iter->col]).im;

            (w[count][iter->col]).re+=(iter->val).re*(v[count][iter->row]).re+(iter->val).im*(v[count][iter->row]).im;
            (w[count][iter->col]).im+=(iter->val).re*(v[count][iter->row]).im-(iter->val).im*(v[count][iter->row]).re;

            count++;
        }
        iter++;
    }
}

/*
void ms_parallel_GenKrylovBasis(Hamiltonian H, Cmplx ** v, Cmplx *** w, unsigned int m, int num) {
    unsigned int count=0;

    ms_parallel_HamXVec(H, v, w[0], num);
    while(count<m) {
        ms_parallel_HamXVec(H, w[count-1], w[count], num);
        count++;
    }
}
*/
