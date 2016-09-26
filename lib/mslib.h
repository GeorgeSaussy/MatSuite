#ifndef MSLIB
#define MSLIB
#include <vector>
#include <iostream>
#include <math.h>
using namespace std;
/**
    Prints a matrix to std::out
    @param A a matrix
    @return void
*/
void matPrint(vector<vector<double> > A);
/**
    Function to add two matricies to one another
    @param A the first summand
    @param B the second summand
    @return the sum A + B
    @warning demensions of A and B must match
*/
vector<vector<double> > matAdd(vector<vector<double> > A, vector<vector<double> > B);
/**
    Function to perform matix multiplication
    @param A first matrix factor
    @param B seond matrix factor
    @return the product A*B
    @warning the width of A must equal the height of B
*/
vector<vector<double> > matMult(vector<vector<double> > A, vector<vector<double> > B);
/**
    Function to perform scalar multiplication on a matrix
    @param A a sqare matrix
    @param t a scalar
    @return the value A*t
*/
vector<vector<double> > scalarMult(vector<vector<double> > A, double t);
/**
    Function to calculate the one norm of a matrix
    @param A a square matrix
    @return the one norm of A
*/
double oneNorm(vector<vector<double> > A);
/**
    Function to calculate the traspose of a retangular matrix
    @param P a matrix
    @return the transpose of P
*/
vector<vector<double> > matT(vector<vector<double> > P);
/**
    Function to swap the rows of a matrix
    @param A a matrix
    @param i the first row index to swap
    @param j the second row index to swap
    @return a matrix with rows i and j interchanged
    @warning must 0 <= i,j < A.size()
*/
vector<vector<double> > swapRows(vector<vector<double> >, int i, int j);
/**
    Function to calculate the inverse of a matrix
    @param A a square invertable matrix
    @return the inverse of A
    @warning A cannot be singular
*/
vector<vector<double> > matInv(vector<vector<double> > A);
/**
    Calculate the Q function for Pade's approximation
    @param A a square matrix
    @param p degree of polynomial approximation
    @return a square matrix given by (2.1) in Ward 1977
*/
vector<vector<double> > Q(vector<vector<double> > A, int p);
/**
    First draft function to calculate exp(At), where A is a matrix and t is a
    scalar, following Ward 1977
    @param A a square matrix
    @param t a scalar
    @return exp(At) as a matix
*/
vector<vector<double> > matExp(vector<vector<double> > A);
#endif
