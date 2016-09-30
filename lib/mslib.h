#ifndef MSLIB
#define MSLIB
// TODO add print matrix
/**
    Struct to hold a matrix
*/
struct SqMat {
    int N; // the demension of the matrix
    double * pValue; // the values in the matrix, stored left-to-right top-to-bottom
}
/**
    Copy one matrix to another without worring about memory leaks
    @param mat1 the matrix to copy
    @param mat2 a reference to where the matrix is to be copied
    @return void
*/
void safeCopySqMat(SqMat mat1, SqMat * mat2);
/**
    Initialize a square matrix to be a matrix of all zeroes
    @param mat a reference to the matrix to initialize
    @param N the size of the matrix
*/
void initZeroSqMat(SqMat * mat, int N);
/**
    Initialize a square matrix to be an identity matrix
    @param mat a reference to the matrix to initialize
    @param N the size of the matrix
*/
void initIdSqMat(SqMat * mat, int N);
/**
    Function to get the value from the matrix
    @param mat a reference to the matrix
    @param i the row of the value
    @param j the column of the value
    @return the value found in the element
    @warning returns 0.0 if indexes not in range
*/
double getVal(SqMat * mat, int i, int j);
/**
    Function to set a value in the matrix
    @param mat a reference to a matrix
    @param i the row of the value
    @param j the column of the matrix
    @param x the value to set
    @return void
    @warning fails if indexes not in range
*/
void setVal(SqMat * mat, int i, int j, double x);
/**
    Calculate the value of a matrix scaled by a scalar
    @param lambda the scalar
    @param mat the matrix
    @return the scaled matrix
*/
SqMat scaleSqMat(SqMat mat, double lambda);
/**
    Add two matricies
    @param mat1 the first matrix
    @param mat2 the second matrix
    @return the sum
    @warning addition fails if mat1.N!=mat2.N
*/
SqMat addSqMat(SqMat mat1, SqMat mat2);
/**
    Function to multiply two matricies
    @param mat1 the first factor
    @param mat2 the second factor
    @return the product
    @warning multiplication fails if mat1.N!=mat2.N
*/
SqMat multSqMat(SqMat mat1, SqMat mat2);
/**
    Function to calculate the power of a matrix
    @param mat the base matrix
    @param j the power
    @return mat^j
    @warning must have j>=0
*/
SqMat powSqMat(SqMat mat, int j);
/**
    Function to calculate the transpose of a square matrix
    @param mat the matrix to transpose
    @return the transpose of mat
*/
SqMat transSqMat(SqMat mat);
/**
    Swap rows within a matrix
    @param mat a refenrece to the matrix on which to operate
    @param i the first rows index
    @param j the second row index
    @return void
*/
void swapRows(SqMat * mat, int i, int j);
/** TODO
    Calculate the inverse of a matrix using Gaussian elemination
    @param mat the matrix to invert
    @return the inverse of mat
*/
// SqMat invSqMAt(SqMat mat);
#endif