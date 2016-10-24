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
    A struct to store a complex number
*/
typedef struct Complex {
    double re;
    double im;
} Complex;
/**
    Function to add complex numbers
    @param num1 the first complex number
    @param num2 the second complex number
    @return the sum nnum1+num2
*/
Complex addCmplx(Complex num1, Complex num2);
/**
    Function to subtract complex numbers
    @param num1 the first complex number
    @param num2 the second complex number
    @return the difference num1-num2
*/
Complex subCmplx(Complex num1, Complex num2);
/**
    Function to multiply complex numbers
    @param num1 the first factor
    @param num2 the second factor
    @return the product num1*num2
*/
Complex multCmplx(Complex num1, Complex num2);
/**
    Function to divide complex numbers
    @param num1 the numerator
    @param num2 the denominator
    @return the fraction num1/num2
*/
Complex divCmplx(Complex num1, Complex num2);
/**
    A struct to store the a complex matrix
*/
struct SqMat {
    int N;
    Complex * pValue;
};
/**
    Copy one matrix to another without worring about memory leaks
    @param mat1 the matrix to copy
    @param mat2 a reference to where the matrix is to be copied
    @return void
*/
void safeCopySqMat(struct SqMat mat1, struct SqMat * mat2);
/**
    Initialize a square matrix to be a matrix of all zeroes
    @param mat a reference to the matrix to initialize
    @param N the size of the matrix
*/
void initZeroSqMat(struct SqMat * mat, int N);
/**
    Function to zero out a matrix
    @param mat the matrix to zero
    @return void
*/
void zeroSqMat(struct SqMat * mat);
/**
    Initialize a square matrix to be an identity matrix
    @param mat a reference to the matrix to initialize
    @param N the size of the matrix
*/
void initIdSqMat(struct SqMat * mat, int N);
/**
    Function to get the value from the matrix
    @param mat a reference to the matrix
    @param i the row of the value
    @param j the column of the value
    @return the value found in the element
    @warning returns 0.0 if indexes not in range
*/
Complex getVal(struct SqMat mat, int i, int j);
/**
    Function to set a value in the matrix
    @param mat a reference to a matrix
    @param i the row of the value
    @param j the column of the matrix
    @param x the value to set
    @return void
    @warning fails if indexes not in range
*/
void setVal(struct SqMat * mat, int i, int j, Complex x);
/**
    Calculate the value of a matrix scaled by a scalar
    @param lambda the scalar
    @param mat the matrix
    @return the scaled matrix
*/
struct SqMat scaleSqMat(struct SqMat mat, Complex lambda);
/**
    Add two matricies
    @param mat1 the first matrix
    @param mat2 the second matrix
    @return the sum
    @warning addition fails if mat1.N!=mat2.N
*/
struct SqMat addSqMat(struct SqMat mat1, struct SqMat mat2);
/**
    Function to multiply two matricies
    @param mat1 the first factor
    @param mat2 the second factor
    @return the product
    @warning multiplication fails if mat1.N!=mat2.N
    @warning do not use, fails often
*/
struct SqMat multSqMat(struct SqMat mat1, struct SqMat mat2);
/**
    Function to safely perform matrix multiplication
    @param mat1 the first factor
    @param mat2 the second factor
    @return the product
*/
void multSafeSqMat(struct SqMat mat1, struct SqMat mat2, struct SqMat * out);
/**
    Function to calculate the power of a matrix
    @param mat the base matrix
    @param j the power
    @return mat^j
    @warning must have j>=0
*/
struct SqMat powSqMat(struct SqMat mat, int j);
/**
    Function to calculate the transpose of a square matrix
    @param mat the matrix to transpose
    @return the transpose of mat
    @warning must j>=0
*/
struct SqMat transSqMat(struct SqMat mat);
/**
    Swap rows within a matrix
    @param mat a refenrece to the matrix on which to operate
    @param i the first rows index
    @param j the second row index
    @return void
*/
void swapRows(struct SqMat * mat, int i, int j);
/**
    Calculate the inverse of a matrix using Gaussian elemination
    @param mat the matrix to invert
    @return the inverse of mat
*/
struct SqMat invSqMat(struct SqMat mat);
/**
    Calculate the one norm of a matrix
    @param mat the matrix argument
    @return double the one norm of mat
*/
double oneNorm(struct SqMat mat);
#endif
