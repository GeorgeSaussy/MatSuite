#include "mslib.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>


void ms_HamXVec(Hamiltonian H, Cmplx * v, Cmplx * w) {
    Cell * iter=H.val;

    while((unsigned int)(iter-H.val)<H.len) {
        // Upper right part
        (w[iter->row]).re+=(iter->val).re*(v[iter->col]).re-(iter->val).im*(v[iter->col]).im;
        (w[iter->row]).im+=(iter->val).re*(v[iter->col]).im+(iter->val).re*(v[iter->col]).im;

        // Conjugate
        (w[iter->col]).re+=(iter->val).re*(v[iter->row]).re+(iter->val).im*(v[iter->row]).im;
        (w[iter->col]).im+=(iter->val).re*(v[iter->row]).im-(iter->val).im*(v[iter->row]).re;

        iter++;
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
    Cell * iter=H.val;
    int count=0;

    while((unsigned int)(iter-H.val)<H.len) {
        count=0;
        while(count<num) {
            (w[count][iter->row]).re+=(iter->val).re*(v[iter->col]).re-(iter->val).im*(v[count][iter->col]).im;
            (w[count][iter->row]).im+=(iter->val).re*(v[iter->col]).im+(iter->val).re*(v[count][iter->col]).im;

            (w[count][iter->col]).re+=(iter->val).re*(v[iter->row]).re+(iter->val).im*(v[count][iter->row]).im;
            (w[count][iter->col]).im+=(iter->val).re*(v[iter->row]).im-(iter->val).im*(v[count][iter->row]).re;

            count++;
        }
        iter++;
    }
}

void ms_parallel_GenKrylovBasis(Hamiltonian H, Cmplx ** v, Cmplx *** w, unsigned int m, int num) {
    unsigned int count=0;

    ms_parallel_HamXVec(H, v, w[0], num);
    while(count<m) {
        ms_parallel_HamXVec(H, w[count-1], w[count], num);
        count++;
    }
}
