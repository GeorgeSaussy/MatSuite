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
    unsigned int col;
    Cmplx val;
} Cell;

typedef struct Hamiltonian {
    unsigned int N;
    unsigned int len;
    Cell * ptr;
} Hamiltonian;

// IMPLEMENTATION 1 FUNCTIONS

/**
    Multiply a vector by a Hamiltonian and store the result at a given location
    @param H the Hamiltonian (in the spaceial format)
    @param v the vector, product
    @param w where to store the result, should be maqlloc'd and zeroed
    @return void
*/
void ms_HamXVec(Hamiltonian H, Cmplx * v, Cmplx * w) {
    Cell * iter=H.ptr;

    while((unsigned int)(iter-H.ptr)<H.len) {
        (w[iter->row]).re+=(iter->val).re*(v[iter->col]).re-(iter->val).im*(v[iter->col]).im;
        (w[iter->row]).im+=(iter->val).re*(v[iter->col]).im+(iter->val).re*(v[iter->col]).im;

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

/**
    Generate a random Hamiltonian
    @param sp the max number of nonzero elements per row
    @param dim the dimention of the matrix
    @return a pointer to a random Hamiltonian
*/
Hamiltonian  ms_generateRandomHamiltonian(int sp, int dim) {
    Hamiltonian ret;
    // first make (dim / sp) groupings of sp rows, done with Knuths permtation shuffle
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
    // generate (dim / sp) Hamiltonians
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
    ret.N=dim;
    ret.len=sp*sp*numberOfGroups+(dim-sp*numberOfGroups);
    ret.ptr=(Cell*)malloc(ret.len*sizeof(Cell));
    j=0;
    for(k=0;k<numberOfGroups;k++) {
        for(k1=0;k1<sp;k1++) {
            for(k2=0;k2<sp;k2++) {
                ((ret.ptr)[j]).row=groups[k][k1];
                ((ret.ptr)[j]).col=groups[k][k2];
                ((ret.ptr)[j]).val=ham[k][k1][k2];
                j++;
            }
        }
    }
    for(j=0;j<ret.N-sp*numberOfGroups;j++) {
        ((ret.ptr)[j]).row=basearray[sp*numberOfGroups+j];
        ((ret.ptr)[j]).col=((ret.ptr)[j]).row;
        ((ret.ptr)[j]).val.re=1;
        ((ret.ptr)[j]).val.im=1;
    }
    return ret;
}

// IMPLEMENTATION 1 TESTS

/**
    Functin to test implementation with correct sparcity and an invertable
    matrix, matrix is a tensor of Pauli matricies
    Implementation 1
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

/**
    Functin to test implementation with correct sparcity and an invertable
    matrix, matrix is a random sparce matrix
    Implementation 1
*/
void ms_test4_GenKrylovBasis(int lowdim, int hidim, int skip, int iter, int sparcity) {
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
        nonzero=dim;
        H=ms_generateRandomHamiltonian(sparcity, dim);
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
                //printf("Hello1\n");
                ms_GenKrylovBasis(H, v, w, m);
                //printf("Hello2\n");
            end = clock();
            avg += (double)(end - begin);
        }
        avg /= CLOCKS_PER_SEC * iter;
        printf("%d\t%f\n", dim, avg);
        //free(H.ptr);
        free(v);
        free(w);
    }
}

// MAIN

int main(int argc, char** argv) {
    //ms_test3_GenKrylovBasis(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4])); // 1000 10000 50 1000
    ms_test4_GenKrylovBasis(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5])); //  1000 10000 50 1000 4
    return 0;
}
