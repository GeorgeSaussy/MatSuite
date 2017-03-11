#include "mslib.h"
#include <stdlib.h>
#include <stdio.h>


void ms_HamXVec(Hamiltonian H, Cmplx * v, Cmplx * w) {
    int done=0;
    int iter=0;
    int iter1=0;
    for(iter=0;iter<H.N;iter++) {
        for(iter1=0;iter1<H.num_rows_per_col;iter1++) {
            //printf("%d\t%d\n",iter,iter1);
            w[(H.ptr)[H.num_rows_per_col*iter+iter1].row].re +=
                        (H.ptr)[H.num_rows_per_col*iter+iter1].val.re * v[iter].re
                        - (H.ptr)[H.num_rows_per_col*iter+iter1].val.im * v[iter].im;
            w[(H.ptr)[H.num_rows_per_col*iter+iter1].row].im +=
                        (H.ptr)[H.num_rows_per_col*iter+iter1].val.re * v[iter].im
                        + (H.ptr)[H.num_rows_per_col*iter+iter1].val.im * v[iter].re;
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

/*
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
*/

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
