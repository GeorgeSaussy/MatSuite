#ifndef MSLIB_H
#define MSLIB_H
/** @file */
/**
    Struct to hold a matrix
*/
struct SqMatRl {
    int N; // the demension of the matrix
    double * pValue; // the values in the matrix, stored left-to-right top-to-bottom
};
/**
    Copy one matrix to another without worring about memory leaks
    @param mat1 the matrix to copy
    @param mat2 a reference to where the matrix is to be copied
    @return void
*/
void safeCopySqMatRl(struct SqMatRl mat1, struct SqMatRl * mat2);
/**
    Initialize a square matrix to be a matrix of all zeroes
    @param mat a reference to the matrix to initialize
    @param N the size of the matrix
*/
void initZeroSqMatRl(struct SqMatRl * mat, int N);
/**
    Initialize a square matrix to be an identity matrix
    @param mat a reference to the matrix to initialize
    @param N the size of the matrix
*/
void initIdSqMatRl(struct SqMatRl * mat, int N);
/**
    Function to get the value from the matrix
    @param mat a reference to the matrix
    @param i the row of the value
    @param j the column of the value
    @return the value found in the element
    @warning returns 0.0 if indexes not in range
*/
double getValRl(struct SqMatRl mat, int i, int j);
/**
    Function to set a value in the matrix
    @param mat a reference to a matrix
    @param i the row of the value
    @param j the column of the matrix
    @param x the value to set
    @return void
    @warning fails if indexes not in range
*/
void setValRl(struct SqMatRl * mat, int i, int j, double x);
/**
    Calculate the value of a matrix scaled by a scalar
    @param lambda the scalar
    @param mat the matrix
    @return the scaled matrix
*/
struct SqMatRl scaleSqMatRl(struct SqMatRl mat, double lambda);
/**
    Add two matricies
    @param mat1 the first matrix
    @param mat2 the second matrix
    @return the sum
    @warning addition fails if mat1.N!=mat2.N
*/
struct SqMatRl addSqMatRl(struct SqMatRl mat1, struct SqMatRl mat2);
/**
    Function to multiply two matricies
    @param mat1 the first factor
    @param mat2 the second factor
    @return the product
    @warning multiplication fails if mat1.N!=mat2.N
*/
struct SqMatRl multSqMatRl(struct SqMatRl mat1, struct SqMatRl mat2);
/**
    Function to calculate the power of a matrix
    @param mat the base matrix
    @param j the power
    @return mat^j
    @warning must have j>=0
*/
struct SqMatRl powSqMatRl(struct SqMatRl mat, int j);
/**
    Function to calculate the transpose of a square matrix
    @param mat the matrix to transpose
    @return the transpose of mat
*/
struct SqMatRl transSqMatRl(struct SqMatRl mat);
/**
    Swap rows within a matrix
    @param mat a refenrece to the matrix on which to operate
    @param i the first rows index
    @param j the second row index
    @return void
*/
void swapRowsRl(struct SqMatRl * mat, int i, int j);
/**
    Calculate the inverse of a matrix using Gaussian elemination
    @param mat the matrix to invert
    @return the inverse of mat
*/
struct SqMatRl invSqMatRl(struct SqMatRl mat);
/**
    Calculate the one norm of a matrix
    @param mat the matrix argument
    @return double the one norm of mat
*/
double oneNormRl(struct SqMatRl mat);
/**
    Function to calculate the N matrix for the Pade approximation
    @param mat the matrix argument
    @param p for p, q Pade approximation
    @param q for p, q Pade approximation
    @return $ N_{p,q} (mat) $
*/
struct SqMatRl PadeNRl(struct SqMatRl mat, int p, int q);
/**
    Function to calculate the D matrix for the Pade approximation
    @param mat the matrix argument
    @param p first argument for p, q Pade approximation
    @param q second argument for p, q Pade approximation
    @return $ D_{p,q} (mat) $
*/
struct SqMatRl PadeDRl(struct SqMatRl mat, int p, int q);
/**
    Calculate exp(mat) with the p, q Pade method
    @param mat the matrix to exponentiate
    @param p first argument for p, q Pade approximation
    @param q second argument for p, q Pade approximation
    @return p, q Pade approximation of exp(mat)
*/
struct SqMatRl expPadeRl(struct SqMatRl mat, int p, int q);
/**
    Function to perform the Krylov aproximation for matrix mat and vector v
    @param mat a matrix to exponentiate
    @param v a vector by which to multiply exp(t*mat)
    @param tau the step-size
    @param minerr the roundoff error used
    @return another vector of the same demension as v as a pointer to doubles
*/
double * expvKrylovRl(struct SqMatRl mat, double t, double * v, double tau, double minerr);
/**
    Struct to store an imaginary number
*/
struct Complex {
    double re; // the real part of the complex number
    double im; // the imaginary part of the complex number
};
/**
    Function to store a complex matrix
*/
struct SqMat {
    int N; // the sixze of the matrix
    struct Complex * pValue;
};

#endif
