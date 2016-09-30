#include "mslib.h"
void initZeroSqMat(SqMat * mat, int N) {
    *mat.N=N;
    *mat.pValue=malloc(N*N); // TODO check malloc sucedes
    int k=0;
    for(;k<N*N;k++) {
        *((*mat).pValue+k*sizeof(double))=0;
    }
}
void initIdSqMat(SqMat * mat, int N) {
    *mat.N=N;
    *mat.pValue=malloc(N*N); // TODO check malloc sucedes
    int k=0;
    int k1=0;
    for(;k<N;k++) {
        for(k1=0;k1<N;k1++) {
            if(k1==k) {
                *((*mat).pValue+(k*N+k1)*sizeof(double))=1;
            }
            else {
                *((*mat).pValue+(k*N+k1)*sizeof(double))=0;
            }
        }
    }
}
double getVal(SqMat * mat, int i, int j) {
    double toret=0.0;
    if(0<=i && i<N && 0<=j && j<N) {
        toret=*((*mat).pValue+(i*(*mat).N+j)*sizeof(double));
    }
    return toret;
}
void setVal(SqMat * mat, int i, int j, double x) {
    if(0<=i && i<N && 0<=j && j<N) {
        *((*mat).pValue+(i*(*mat).N+j)*sizeof(double))=x;
    }
}
SqMat scaleSqMat(SqMat mat, double lambda) {
    SqMat toret;
    initZeroSqMat(&toret,mat.N);
    int k=0;
    int k1=0;
    for(;k<toret.N;k++) {
        for(k1=0;k1<toret.N;k1++) {
            setVal(&toret,k,k1,lambda*getVal(&mat,k,k1));
        }
    }
    return toret;
}
SqMat addSqMat(SqMat mat1, SqMat mat2) {
    SqMat toret;
    if(mat1.N==mat2.N) {
        initZeroSqMat(&toret,mat1.N);
        int k=0;
        int k1=0;
        for(;k<mat1.N;k++) {
            for(k1=0;k1<mat1.N;k1++) {
                setVal(&toret,k,k1,getVal(mat1,k,k1)+getVal(mat2,k,k1));
            }
        }
    }
    return toret;
}
SqMat multSqMat(SqMat mat1, SqMat mat2) {
    SqMat toret;
    if(mat1.N==mat2.N) {
        initZeroSqMat(&toret,mat1.N);
        int k=0;
        int k1=0;
        int k2=0;
        double x=0;
        for(;k<mat1.N;k++) {
            for(k1=0;k1<mat1.N;k1++) {
                x=0;
                for(k2=0;k2<mat1.N;k2++) {
                    x+=getVal(mat1,k,k2)*getVal(mat2,k2,k1);
                }
                setVal(&toret,k,k1,x);
            }
        }
    }
}
