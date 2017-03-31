#ifndef MSLIB_H
#define MSLIB_H

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

/**
    Multiply a vector by a Hamiltonian and store the result at a given location
    @param H the Hamiltonian (in the spaceial format)
    @param v the vector, product
    @param w where to store the result, should be maqlloc'd and zeroed
    @return void
*/
void ms_HamXVec(Hamiltonian * H, Cmplx * v, Cmplx * w);
/**
    Generate the Krylov basis for a Hamiltonian and vector and store the result
    in a a given location
    @param H the Hamiltonian
    @param v the vector
    @param w where to store the result, should be malloc'd and zeroed
    @param m the number of Kylov basis to calculate, must be > 1, should be > 10
    @return void
*/
void ms_GenKrylovBasis(Hamiltonian * H, Cmplx * v, Cmplx ** w, unsigned int m);
/**
    Generate a random Hamiltonian
    @param sp the max number of nonzero elements per row
    @param dim the dimention of the matrix
    @return a pointer to a random Hamiltonian
*/
Hamiltonian ms_generateRandomHamiltonian(int sp, int dim);
/**
    Calculate the norm of a vector
    @param p the vector
    @param dim the dimention of the vector
    @return the norm
*/
double vec_norm(Cmplx * p, int dim);
/**
    Compute a matrix exponential
    @param H the matrix to exponentiate as a (Cmplx**)
    @param F where to store the result
    @param S a workspace to store a matrix
    @param p a workspace to store a vector
    @param dim the dimention of the matrix
    @param iter the number of iterations of the matrix exponential to perform
    (must be >0)
    @return void
*/
void ms_expm(Cmplx ** H, Cmplx ** F, Cmplx ** S, Cmplx * p, int dim, int iter);
/**
    Implement Krylov exponentiation

    @param t time
    @param A a pointer to a Hamiltonian
    @param normA the infinity norm of A (should be precomputed to save time)
    @param v the vector, norm must be == 1
    @paran w where to store the answer ( Dim(ans) == Dim(v) == A->N )

    @param V a workspace (to avoid repeated heap allocations)
    @param H a workspace
    @param F a workspace
    @param S a workspace
    @param p a workspace
    @param workspace a workspace

    @return void
*/
void ms_expv(double t, Hamiltonian * A, double normA, Cmplx * v, int m,
    Cmplx * w, Cmplx ** V, Cmplx ** H, Cmplx ** F, Cmplx ** S, Cmplx * p,
    Cmplx * workspace);
#endif
