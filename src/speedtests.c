#include "mslib.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <assert.h>

/**
    Function to test performance of implementation with oversparce matricies
*/
void ms_test1_GenKrylovBasis(int lowdim, int hidim, int skip, int nonzero, int iter) {
    // void ms_GenKrylovBasis(Hamiltonian H, Cmplx * v, Cmplx ** w, unsigned int m);
    Hamiltonian H;
    Cmplx * v;
    Cmplx ** w;
    Cmplx zero;
    int m=30;
    clock_t begin;
    clock_t end;
    double avg;

    zero.re=0.0;
    zero.im=0.0;
    H.len=nonzero;
    H.N=0; // not needed for test
    H.ptr=(Cell*)malloc(H.len*sizeof(Cell));
    for(int dim=lowdim;dim<=hidim;dim+=skip) {
        // malloc everything
        H.N=dim;
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
            for(int k1=0;k1<nonzero;k1++) {
                ((H.ptr)[k1]).col=rand()%dim;
                ((H.ptr)[k1]).row=((H.ptr)[k1]).col+rand()%(dim-((H.ptr)[k1]).col);
                ((H.ptr)[k1]).val.re=((double)rand()/(double)RAND_MAX);
                ((H.ptr)[k1]).val.im=((double)rand()/(double)RAND_MAX);
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
    }
}


/**
    Function to test implementation with correct sparcity density although
    matrix still not invertable
*/
void ms_test2_GenKrylovBasis(int lowdim, int hidim, int skip, int iter) {
    // void ms_GenKrylovBasis(Hamiltonian H, Cmplx * v, Cmplx ** w, unsigned int m);
    Hamiltonian H;
    Cmplx * v;
    Cmplx ** w;
    Cmplx zero;
    int m=30;
    clock_t begin;
    clock_t end;
    double avg;
    int nonzero;

    zero.re=0.0;
    zero.im=0.0;
    H.N=0; // not needed for test
    for(int dim=lowdim;dim<=hidim;dim+=skip) {
        // malloc everything
        H.N=dim;
        H.len=dim;
        nonzero=dim;
        H.ptr=(Cell*)malloc(H.len*sizeof(Cell));
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
            for(int k1=0;k1<nonzero;k1++) {
                ((H.ptr)[k1]).col=rand()%dim;
                ((H.ptr)[k1]).row=((H.ptr)[k1]).col+rand()%(dim-((H.ptr)[k1]).col);
                ((H.ptr)[k1]).val.re=((double)rand()/(double)RAND_MAX);
                ((H.ptr)[k1]).val.im=((double)rand()/(double)RAND_MAX);
            }
            begin = clock();
                ms_GenKrylovBasis(H, v, w, m);
            end = clock();
            avg += (double)(end - begin);
        }
        avg /= CLOCKS_PER_SEC * iter;
        printf("%d\t%f\n", dim, avg);
        free(H.ptr);
        free(v);
        free(w);
    }
}

/**
    Functin to test implementation with correct sparcity and an invertable
    matrix, matrix is a tensor of Pauli matricies
*/
void ms_test3_GenKrylovBasis(int lowdim, int hidim, int skip, int iter) {
    // void ms_GenKrylovBasis(Hamiltonian H, Cmplx * v, Cmplx ** w, unsigned int m);
    Hamiltonian H;
    Cmplx * v;
    Cmplx ** w;
    Cmplx zero;
    int m=30;
    clock_t begin;
    clock_t end;
    double avg;
    int nonzero;

    zero.re=0.0;
    zero.im=0.0;
    H.N=0; // not needed for test
    for(int dim=lowdim;dim<=hidim;dim+=skip) {
        // malloc everything
        H.N=dim;
        H.len=dim;
        nonzero=dim;
        H.ptr=(Cell*)malloc(H.len*sizeof(Cell));
        for(int k=0;k<nonzero;k++) {
            ((H.ptr)[k]).col=k;
            ((H.ptr)[k]).row=k;
            if(k%2==1) {
                ((H.ptr)[k]).val.re=1;
                ((H.ptr)[k]).val.im=0;
            }
            else {
                ((H.ptr)[k]).val.re=-1;
                ((H.ptr)[k]).val.im=0;
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
        free(H.ptr);
        free(v);
        free(w);
    }
}

int main(int argc, char** argv) {
    ms_test3_GenKrylovBasis(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4])); // 10 100000 1 10 1000
    return 0;
}
