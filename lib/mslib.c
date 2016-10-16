#include "mslib.h"
#include <stdlib.h>
#include <math.h>
void safeCopySqMat(struct SqMat mat1, struct SqMat * mat2) {
    free(mat2->pValue);
    mat2->N=mat1.N;
    mat2->pValue=malloc(mat1.N*mat1.N*sizeof(double));
    int k=0;
    for(;k<mat1.N*mat1.N;k++) {
        *(mat2->pValue+k*sizeof(double))=*(mat1.pValue+k*sizeof(double));
    }
}
void initZeroSqMat(struct SqMat * mat, int N) {
    mat->N=N;
    mat->pValue=malloc(N*N*sizeof(double)); // TODO check malloc succeeds
    int k=0;
    for(;k<N*N;k++) {
        *(mat->pValue+k*sizeof(double))=0;
    }
}
void initIdSqMat(struct SqMat * mat, int N) {
    mat->N=N;
    mat->pValue=malloc(N*N*sizeof(double)); // TODO check malloc succeeds
    int k=0;
    int k1=0;
    for(;k<N;k++) {
        for(k1=0;k1<N;k1++) {
            if(k1==k) {
                *(mat->pValue+(k*N+k1)*sizeof(double))=1;
            }
            else {
                *(mat->pValue+(k*N+k1)*sizeof(double))=0;
            }
        }
    }
}
double getVal(struct SqMat mat, int i, int j) {
    double toret=0.0;
    if(0<=i && i<mat.N && 0<=j && j<mat.N) {
        toret=*(mat.pValue+(i*mat.N+j)*sizeof(double));
    }
    return toret;
}
void setVal(struct SqMat * mat, int i, int j, double x) {
    if(0<=i && i<mat->N && 0<=j && j<mat->N) {
        *(mat->pValue+(i*mat->N+j)*sizeof(double))=x;
    }
}
struct SqMat scaleSqMat(struct SqMat mat, double lambda) {
    struct SqMat toret;
    initZeroSqMat(&toret,mat.N);
    int k=0;
    int k1=0;
    for(;k<toret.N;k++) {
        for(k1=0;k1<toret.N;k1++) {
            setVal(&toret,k,k1,lambda*getVal(mat,k,k1));
        }
    }
    return toret;
}
struct SqMat addSqMat(struct SqMat mat1, struct SqMat mat2) {
    struct SqMat toret;
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
struct SqMat multSqMat(struct SqMat mat1, struct SqMat mat2) {
    struct SqMat toret;
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
struct SqMat powSqMat(struct SqMat mat, int j) {
    struct SqMat toret;
    if(j>=0) {
        initIdSqMat(&toret,mat.N);
        int k=0;
        for(;k<mat.N;k++) {
            safeCopySqMat(multSqMat(toret,mat),&toret);
        }
    }
    return toret;
}
struct SqMat transSqMat(struct SqMat mat) {
    struct SqMat toret;
    initZeroSqMat(&toret,mat.N);
    int k=0;
    int k1=0;
    for(;k<mat.N;k++) {
        for(k1=0;k1<mat.N;k1++) {
            setVal(&toret,k,k1,getVal(mat,k1,k));
        }
    }
}
void swapRows(struct SqMat * mat, int i, int j) {
    double * tempRow=malloc(mat->N*sizeof(double));
    int k=0;
    for(;k<mat->N;k++) {
        *(tempRow+k*sizeof(double))=getVal(*mat,i,k);
    }
    for(k=0;k<mat->N;k++) {
        setVal(mat,i,k,getVal(*mat,j,k));
    }
    for(k=0;k<mat->N;k++) {
        setVal(mat,j,k,*(tempRow+k*sizeof(double)));
    }
}
struct SqMat invSqMat(struct SqMat mat) {
    int n=mat.N;
    struct SqMat toret;
    initIdSqMat(&toret,n);
    struct SqMat ongoing;
    safeCopySqMat(mat,&ongoing);
    int k=0;
    int k1=0;
    int i=0;
    int j=0;
    int iMax=0;
    for(;k<n;k++) {
        iMax=k;
        for(k1=k;k1<n;k1++) {
            if(abs(getVal(mat,k1,k))>abs(getVal(mat,iMax,k))) {
                iMax=k1;
            }
        }
        // warning: A singular if A[iMax][k] == 0
        swapRows(&ongoing,k,iMax);
        swapRows(&toret,k,iMax);
        for(i=k+1;i<n;i++) {
            double f=getVal(ongoing,i,k)/getVal(ongoing,k,k);
            for(j=k+1;j<n;j++) {
                setVal(&ongoing,i,j,getVal(ongoing,i,j)-f*getVal(ongoing,k,j));
                setVal(&toret,i,j,getVal(toret,i,j)-f*getVal(toret,k,j));
            }
            setVal(&ongoing,i,k,0.0);
        }
    }
    // reduce
    for(k=n-1;k>=0;k--) {
        for(i=k-1;i>=0;i--) {
            double f=getVal(ongoing,i,k)/getVal(ongoing,k,k);
            for(j=0;j<n;j++) {
                setVal(&ongoing,i,j,getVal(ongoing,i,j)-f*getVal(ongoing,i,j));
                setVal(&toret,i,j,getVal(toret,i,j)-f*getVal(toret,i,j));
            }
        }
    }
    return toret;
}
double oneNorm(struct SqMat mat) {
    double toret=0;
    int k=0;
    int k1=0;
    for(;k<mat.N;k++) {
        for(k1=0;k<mat.N;k1++) {
            if(abs(get(mat,k,k1))>toret) {
                toret=abs(get(mat,k,k1));
            }
        }
    }
    return toret;
}
struct SqMat PadeN(struct SqMat mat, int p, int q) {
    struct SqMat toret;
    initIdSqMat(&toret,mat.N);
    struct SqMat ongoing;
    safeCopySqMat(toret,&ongoing);
    int k=1;
    double lambda=1.0;
    for(;k<p;k++) {
        lambda*=((double) (p-k+1))/((double) (p+q-k+1))/((double) k);
        safeCopySqMat(multSqMat(ongoing,mat),&ongoing);
        safeCopySqMat(addSqMat(scaleSqMat(ongoing,lambda),toret),&toret);
    }
    return toret;
}
struct SqMat PadeD(struct SqMat mat, int p, int q) {
    struct SqMat toret;
    initIdSqMat(&toret,mat.N);
    struct SqMat ongoing;
    safeCopySqMat(toret,&ongoing);
    int k=1;
    double lambda=1.0;
    for(;k<p;k++) {
        lambda*=-1*((double) (p-k+1))/((double) (p+q-k+1))/((double) k);
        safeCopySqMat(multSqMat(ongoing,mat),&ongoing);
        safeCopySqMat(addSqMat(scaleSqMat(ongoing,lambda),toret),&toret);
    }
    return toret;
}
struct SqMat expPade(struct SqMat mat, int p, int q) {
    struct SqMat toret;
    double norm=oneNorm(mat);
    int m=1;
    if(norm>0.5) {
        m=1+((int)(norm*2));
    }
    struct SqMat matPrime;
    safeCopySqMat(scaleSqMat(mat,1.0/m),&matPrime);
    safeCopySqMat(multSqMat(invSqMat(PadeD(matPrime,p,q)),PadeN(matPrime,p,q)),&toret);
    safeCopySqMat(powSqMat(toret,m),&toret);
    return toret;
}
double * expvKrylov(struct SqMat mat, double t, double * v, double tau, double minerr) {
    int k=0;
    int k1=0;
    int j=0;
    int i=0;
    double tk=0;
    double beta=0;
    double tempvar=0;
    double errloc=0;
    double err1=0;
    double err2=0;
    double matvnorm=0;
    //double tau=gamma*(tol/eta)**(1/r)
    struct SqMat H;
    initZeroSqMat(&H, mat.N+2);
    struct SqMat F;
    initZeroSqMat(&F,mat.N+2);
    double * w=malloc(mat.N*sizeof(double));
    double * v1=malloc(mat.N*sizeof(double));
    double * p=malloc(mat.N*sizeof(double));
    for(k=0;k<mat.N;k++) {
        *(w+k*sizeof(double))=0;
        for(k1=0;k1<mat.N;k1++) {
            *(w+k*sizeof(double))+=getVal(mat,k,k1)*(*(v+k1*sizeof(double)));
        }
    }
    for(k=0;k<mat.N;k++) {
        matvnorm+=(*(w+k*sizeof(double)))*(*(w+k*sizeof(double)));
    }
    matvnorm=sqrt(matvnorm);
    for(k=0;k<mat.N;k++) {
        *(w+k*sizeof(double))=*(v+k*sizeof(double));
    }
    while(tk<t) {
        for(k=0;k<mat.N;k++) {
            *(v+k*sizeof(double))=*(w+k*sizeof(double));
            beta+=(*(v+k*sizeof(double)))*(*(v+k*sizeof(double)));
        }
        beta=sqrt(beta);
        for(k=0;k<mat.N;k++) {
            *(v1+k*sizeof(double))=*(v+k*sizeof(double))/beta;
        }
        for(j=0;j<mat.N;j++) { // Arnoldi process
            for(k=0;k<mat.N;k++) {
                *(p+k*sizeof(double))=0;
                for(k1=0;k1<mat.N;k1++) {
                    *(p+k*sizeof(double))+=getVal(mat,k,k1)*(*(v1+k*sizeof(double)));
                }
            }
            for(i=0;i<mat.N;i++) {
                tempvar=0;
                for(k=0;k<mat.N;k++) {
                    tempvar+=(*(v1+k*sizeof(double)))*(*(p+k*sizeof(double)));
                }
                setVal(&H,i,j,tempvar);
                for(k=0;k<mat.N;k++) {
                    *(p+k*sizeof(double))=*(p+k*sizeof(double))-getVal(H,i,j)*(*(v1+k*sizeof(double)));
                }
            }
            tempvar=0;
            for(k=0;k<mat.N;k++) {
                tempvar+=(*(p+k*sizeof(double)))*(*(p+k*sizeof(double)));
            }
            setVal(&H,mat.N,mat.N-1,tempvar);
            //if(getVal(H,mat.N,mat.N)<=tol_abs(mat)) { // TODO WHAT IS THIS ???
            //    happy_breakdown(); // TODO WHAT IS THIS ???
            //}
            for(k=0;k<mat.N;k++) {
                *(v1+k*sizeof(double))=*(p+k*sizeof(double))/getVal(H,mat.N,mat.N-1);
            }
        }
        setVal(&H,mat.N+1,mat.N,1.0);
        double delta_tol=.01; // TODO add delta_tol
        errloc=delta_tol+1;
        while(errloc>delta_tol) {
            F=expPade(scaleSqMat(H,tau),14,14); // TODO pick exp
            for(k=0;k<mat.N;k++) {
                *(w+k*sizeof(double))=beta*(*(v1+k*sizeof(double)))*getVal(F,k,1);
            }
            // LOCAL TRUNCATION ERROR ESTIMATE
            err1=beta*abs(getval(H,mat.N,mat.N-1)*getVal(F,1,0));
            err2=beta*matvnorm*abs(getval(H,mat.N,mat.N-1)*getVal(F,2,0));
            if(err1>10*err2) {
                errloc=err2;
            }
            else if(err1>err2) {
                errloc=err2/(1-(err2/err1));
            }
            else {
                errloc=err1;
            }
        }
        tk+=tau;
    }
    return w;
}
