#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "time.h"
#include "assert.h"

// STRUCTS TO AID IN COMPUTATION

typedef struct Cmplx {
    double re;
    double im;
} Cmplx;

typedef struct Cell {
    unsigned int row;
    Cmplx val;
} Cell;

typedef struct Hamiltonian {
    unsigned int N;
    unsigned int num_rows_per_col;
    Cell * ptr;
} Hamiltonian;

// CORE OF IMPLEMENTATION 2

/**
    Multiply a vector by a Hamiltonian and store the result at a given location
    @param H the Hamiltonian (in the spaceial format)
    @param v the vector, product
    @param w where to store the result, should be maqlloc'd and zeroed
    @return void
*/
void ms_HamXVec(Hamiltonian H, Cmplx * v, Cmplx * w) {
    Cell * iter=H.ptr;
    int len = H.N*H.num_rows_per_col;
    int col = 0;

    while((unsigned int)(iter-H.ptr)<len) {
        (w[iter->row]).re+=(iter->val).re*(v[col]).re-(iter->val).im*(v[col]).im;
        (w[iter->row]).im+=(iter->val).re*(v[col]).im+(iter->val).re*(v[col]).im;

        if((unsigned int)(iter-H.ptr)%H.num_rows_per_col==H.num_rows_per_col-1) {
            col++;
        }
        iter++;
    }
}

/**
    Generate the Krylov basis for a Hamiltonian and vector and store the result
    in a a given location
    @param H the Hamiltonian
    @param v the vector
    @param w where to store the result, should be malloc'd and zeroed
    @param m the number of Kylov basis to calculate, must be > 1, should be > 10
    @return void
*/
void ms_GenKrylovBasis(Hamiltonian H, Cmplx * v, Cmplx ** w, unsigned int m) {
    unsigned int count=1;

    ms_HamXVec(H, v, w[0]);
    while(count<m) {
        ms_HamXVec(H, w[count-1], w[count]);
        count++;
    }
}

// IMPLEMENTATION 2 TESTS

/**
    Functin to test implementation with correct sparcity and an invertable
    matrix, matrix is a tensor of Pauli matricies
    Implementation 2
*/
void ms_test3_GenKrylovBasis(int lowdim, int hidim, int skip, int iter) {
    // void ms_GenKrylovBasis(Hamiltonian H, Cmplx * v, Cmplx ** w, unsigned int m);
    Hamiltonian H;
    Cmplx * v;
    Cmplx ** w;
    Cmplx zero;
    Cmplx one;
    Cmplx negone;
    int m=30;
    clock_t begin;
    clock_t end;
    double avg;
    int nonzero;

    zero.re=0.0;
    zero.im=0.0;
    one.re=1.0;
    one.im=0.0;
    negone.re=-1.0;
    negone.im=0.0;
    for(int dim=lowdim;dim<=hidim;dim+=skip) {
        // malloc everything
        H.N=dim;
        nonzero=dim;
        H.num_rows_per_col=1;
        H.ptr=(Cell*)malloc(H.num_rows_per_col*H.N*sizeof(Cell));
        for(int k=0;k<H.N;k++) {
            (H.ptr)[k*H.num_rows_per_col].row=k;
            if(k%2==0) {
                (H.ptr)[k*H.num_rows_per_col].val=one;
            }
            else {
                (H.ptr)[k*H.num_rows_per_col].val=negone;
            }
        }
        v=(Cmplx*)malloc(dim*sizeof(Cmplx));
        w=(Cmplx**)malloc(m*sizeof(Cmplx*));
        for(int k1=0;k1<m;k1++) {
            w[k1]=(Cmplx*)malloc(dim*sizeof(Cmplx));
        }
        avg=0;
        for(int i=0;i<iter;i++) {
            // init data
            for(int k1=0;k1<dim;k1++) {
                v[k1].re=((double)rand()/(double)RAND_MAX);
                v[k1].im=((double)rand()/(double)RAND_MAX);
            }
            for(int k1=0;k1<m;k1++) {
                for(int k2=0;k2<dim;k2++) {
                    w[k1][k2]=zero;
                }
            }
            begin = clock();
                ms_GenKrylovBasis(H, v, w, m);
            end = clock();
            avg += (double)(end - begin);
        }
        avg /= CLOCKS_PER_SEC * iter;
        printf("%d\t%f\n", dim, avg);
        free(v);
        free(w);
        free(H.ptr);
    }
}

// MAIN

int main(int argc, char** argv) {
    ms_test3_GenKrylovBasis(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4])); // 10 100000 1 10 1000
    return 0;
}
