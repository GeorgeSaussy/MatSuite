#include "mslib.h"
void safeCopySqMat(SqMat mat1, SqMat * mat2) {
    free((*mat2).pValue);
    (*mat2).N=mat1.N;
    *mat2.pValue=malloc(mat1.N*mat1.N*sizeof(double));
    int k=0;
    for(;k<mat1.N*mat1.N;k++) {
        *((*mat).pValue+k*sizeof(double))=*(mat1.pValue+k*sizeof(double));
    }
}
void initZeroSqMat(SqMat * mat, int N) {
    *mat.N=N;
    *mat.pValue=malloc(N*N*sizeof(double)); // TODO check malloc sucedes
    int k=0;
    for(;k<N*N;k++) {
        *((*mat).pValue+k*sizeof(double))=0;
    }
}
void initIdSqMat(SqMat * mat, int N) {
    *mat.N=N;
    *mat.pValue=malloc(N*N*sizeof(double)); // TODO check malloc sucedes
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
SqMat powSqMat(SqMat mat, int j) {
    SqMat toret;
    if(j>=0) {
        initIdSqMat(&toret,mat.N);
        int k=0;
        for(;k<mat.N;k++) {
            safeCopySqMat(multSqMat(toret,mat),toret);
        }
    }
    return toret;
}
SqMat transSqMat(SqMat mat) {
    SqMat toret;
    initZeroSqMat(&toret,mat.N);
    int k=0;
    int k1=0;
    for(;k<mat.N;k++) {
        for(k1=0;k1<mat.N;k1++) {
            setVal(&toret,k,k1,getVal(mat,k1,k));
        }
    }
}
void swapRows(SqMat * mat, int i, int j) {
    double * tempRow=malloc(mat.N*sizeof(double));
    int k=0;
    for(;k<mat.N;k++) {
        *(tempRow+k*sizeof(double))=getVal(*mat,i,k);
    }
    for(k=0;k<mat.N;k++) {
        setVal(mat,i,k,getVal(*mat,j,k));
    }
    for(k=0;k<mat.N;k++) {
        setVal(mat,j,k,*(tempRow+k*sizeof(double)));
    }
}
SqMat PadeN(SqMat mat, int p, int q) {
    SqMat toret;
    initIdSqMat(&toret,mat.N);
    SqMat ongoing;
    safeCopySqMat(toret,&ongoing);
    int k=1;
    double lambda=1.0;
    for(;k<p;k++) {
        lambda*=(p-k+1)/(p+q-j+1)/j;
        safeCopySqMat(multSqMat(ongoing,mat),&ongoing);
        safeCopySqMat(addSqMat(scaleSqMat(ongoing,lambda),toret),&toret);
    }
    return toret;
}
SqMat PadeD(SqMat mat, int p, int q) {
    SqMat matPrime;
    safeCopySqMat(scaleSqMat(mat,-1.0),matPrime);
    SqMat toret;
    initIdSqMat(&toret,matPrime.N);
    SqMat ongoing;
    safeCopySqMat(toret,&ongoing);
    int k=1;
    double lambda=1.0;
    for(;k<p;k++) {
        lambda*=(q-k+1)/(p+q-j+1)/j;
        safeCopySqMat(multSqMat(ongoing,matPrime),&ongoing);
        safeCopySqMat(addSqMat(scaleSqMat(ongoing,lambda),toret),&toret);
    }
    return toret;
}
