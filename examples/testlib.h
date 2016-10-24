#ifndef TESTLIB_H
#define TESTLIB_H
/** @file */
#include "mslib.h"
#include <vector>
using namespace std;
/**
    Function to generate a random matrix
    @param N the size of the matrix
    @return a rendom matrix with complex values on the unit circle
    @warning must N > 0
*/
struct SqMat ramdomUnitValElmSqMat(int N);
/**
    Function to print a full matrix
    @param mat a matrix to print
    @return void
*/
void printFullSqMat(struct SqMat mat);
/**
    Generate a vector of random matrixes
    @param m the number of matrixes to Generate
    @param N the size of the matrix
    @warning must m > 0 && N > 0
*/
vector<struct SqMat> genRandUnitValElmSqMatList(int m, int N);
#endif
