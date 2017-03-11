#ifndef MSLIB_H
#define MSLIB_H

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

/**
    Multiply a vector by a Hamiltonian and store the result at a given location
    @param H the Hamiltonian (in the spaceial format)
    @param v the vector, product
    @param w where to store the result, should be maqlloc'd and zeroed
    @return void
*/
void ms_HamXVec(Hamiltonian H, Cmplx * v, Cmplx * w);


/**
    Generate the Krylov basis for a Hamiltonian and vector and store the result
    in a a given location
    @param H the Hamiltonian
    @param v the vector
    @param w where to store the result, should be malloc'd and zeroed
    @param m the number of Kylov basis to calculate, must be > 1, should be > 10
    @return void
*/
void ms_GenKrylovBasis(Hamiltonian H, Cmplx * v, Cmplx ** w, unsigned int m);

/**
    Multiply multiple vectors by a Hamiltonian and store the result at a given
    location
    @param H the Hamiltonian (in the spaceial format)
    @param v the vectors, product
    @param w where to store the result, should be maqlloc'd and zeroed
    @param num the number of vectors to multiply
    @return void
*/
void ms_parallel_HamXVec(Hamiltonian H, Cmplx ** v, Cmplx ** w, int num);

/**
    Generate the Krylov basis for a Hamiltonian and set of vector and store the
    result in a a given location
    @param H the Hamiltonian
    @param v the vectors
    @param w where to store the result, should be malloc'd and zeroed
    @param m the number of Kylov basis to calculate, must be > 1, should be > 10
    @param num the number of vectors for which the Krylov basis should be
    calculated
    @return void
*/
//void ms_parallel_GenKrylovBasis(Hamiltonian H, Cmplx ** v, Cmplx *** w, unsigned int m, int num) {

#endif
