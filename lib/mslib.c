#include "mslib.h"
#include <stdlib.h>
#include <math.h>
void safeCopySqMatRl(struct SqMatRl mat1, struct SqMatRl * mat2) {
    free(mat2->pValue);
    mat2->N=mat1.N;
    mat2->pValue=malloc(mat1.N*mat1.N*sizeof(double));
    int k=0;
    for(;k<mat1.N*mat1.N;k++) {
        *(mat2->pValue+k*sizeof(double))=*(mat1.pValue+k*sizeof(double));
    }
}
void initZeroSqMatRl(struct SqMatRl * mat, int N) {
    mat->N=N;
    mat->pValue=malloc(N*N*sizeof(double)); // TODO check malloc succeeds
    int k=0;
    for(;k<N*N;k++) {
        *(mat->pValue+k*sizeof(double))=0;
    }
}
void initIdSqMatRl(struct SqMatRl * mat, int N) {
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
double getValRl(struct SqMatRl mat, int i, int j) {
    double toret=0.0;
    if(0<=i && i<mat.N && 0<=j && j<mat.N) {
        toret=*(mat.pValue+(i*mat.N+j)*sizeof(double));
    }
    return toret;
}
void setValRl(struct SqMatRl * mat, int i, int j, double x) {
    if(0<=i && i<mat->N && 0<=j && j<mat->N) {
        *(mat->pValue+(i*mat->N+j)*sizeof(double))=x;
    }
}
struct SqMatRl scaleSqMatRl(struct SqMatRl mat, double lambda) {
    struct SqMatRl toret;
    initZeroSqMatRl(&toret,mat.N);
    int k=0;
    int k1=0;
    for(;k<toret.N;k++) {
        for(k1=0;k1<toret.N;k1++) {
            setValRl(&toret,k,k1,lambda*getValRl(mat,k,k1));
        }
    }
    return toret;
}
struct SqMatRl addSqMatRl(struct SqMatRl mat1, struct SqMatRl mat2) {
    struct SqMatRl toret;
    if(mat1.N==mat2.N) {
        initZeroSqMatRl(&toret,mat1.N);
        int k=0;
        int k1=0;
        for(;k<mat1.N;k++) {
            for(k1=0;k1<mat1.N;k1++) {
                setValRl(&toret,k,k1,getValRl(mat1,k,k1)+getValRl(mat2,k,k1));
            }
        }
    }
    return toret;
}
struct SqMatRl multSqMatRl(struct SqMatRl mat1, struct SqMatRl mat2) {
    struct SqMatRl toret;
    if(mat1.N==mat2.N) {
        initZeroSqMatRl(&toret,mat1.N);
        int k=0;
        int k1=0;
        int k2=0;
        double x=0;
        for(;k<mat1.N;k++) {
            for(k1=0;k1<mat1.N;k1++) {
                x=0;
                for(k2=0;k2<mat1.N;k2++) {
                    x+=getValRl(mat1,k,k2)*getValRl(mat2,k2,k1);
                }
                setValRl(&toret,k,k1,x);
            }
        }
    }
}
struct SqMatRl powSqMatRl(struct SqMatRl mat, int j) {
    struct SqMatRl toret;
    if(j>=0) {
        initIdSqMatRl(&toret,mat.N);
        int k=0;
        for(;k<mat.N;k++) {
            safeCopySqMatRl(multSqMatRl(toret,mat),&toret);
        }
    }
    return toret;
}
struct SqMatRl transSqMatRl(struct SqMatRl mat) {
    struct SqMatRl toret;
    initZeroSqMatRl(&toret,mat.N);
    int k=0;
    int k1=0;
    for(;k<mat.N;k++) {
        for(k1=0;k1<mat.N;k1++) {
            setValRl(&toret,k,k1,getValRl(mat,k1,k));
        }
    }
}
void swapRowsRl(struct SqMatRl * mat, int i, int j) {
    double * tempRow=malloc(mat->N*sizeof(double));
    int k=0;
    for(;k<mat->N;k++) {
        *(tempRow+k*sizeof(double))=getValRl(*mat,i,k);
    }
    for(k=0;k<mat->N;k++) {
        setValRl(mat,i,k,getValRl(*mat,j,k));
    }
    for(k=0;k<mat->N;k++) {
        setValRl(mat,j,k,*(tempRow+k*sizeof(double)));
    }
}
struct SqMatRl invSqMatRl(struct SqMatRl mat) {
    int n=mat.N;
    struct SqMatRl toret;
    initIdSqMatRl(&toret,n);
    struct SqMatRl ongoing;
    safeCopySqMatRl(mat,&ongoing);
    int k=0;
    int k1=0;
    int i=0;
    int j=0;
    int iMax=0;
    for(;k<n;k++) {
        iMax=k;
        for(k1=k;k1<n;k1++) {
            if(abs(getValRl(mat,k1,k))>abs(getValRl(mat,iMax,k))) {
                iMax=k1;
            }
        }
        // warning: A singular if A[iMax][k] == 0
        swapRowsRl(&ongoing,k,iMax);
        swapRowsRl(&toret,k,iMax);
        for(i=k+1;i<n;i++) {
            double f=getValRl(ongoing,i,k)/getValRl(ongoing,k,k);
            for(j=k+1;j<n;j++) {
                setValRl(&ongoing,i,j,getValRl(ongoing,i,j)-f*getValRl(ongoing,k,j));
                setValRl(&toret,i,j,getValRl(toret,i,j)-f*getValRl(toret,k,j));
            }
            setValRl(&ongoing,i,k,0.0);
        }
    }
    // reduce
    for(k=n-1;k>=0;k--) {
        for(i=k-1;i>=0;i--) {
            double f=getValRl(ongoing,i,k)/getValRl(ongoing,k,k);
            for(j=0;j<n;j++) {
                setValRl(&ongoing,i,j,getValRl(ongoing,i,j)-f*getValRl(ongoing,i,j));
                setValRl(&toret,i,j,getValRl(toret,i,j)-f*getValRl(toret,i,j));
            }
        }
    }
    return toret;
}
double oneNormRl(struct SqMatRl mat) {
    double toret=0;
    int k=0;
    int k1=0;
    for(;k<mat.N;k++) {
        for(k1=0;k<mat.N;k1++) {
            if(abs(getValRl(mat,k,k1))>toret) {
                toret=abs(getValRl(mat,k,k1));
            }
        }
    }
    return toret;
}
struct SqMatRl PadeNRl(struct SqMatRl mat, int p, int q) {
    struct SqMatRl toret;
    initIdSqMatRl(&toret,mat.N);
    struct SqMatRl ongoing;
    safeCopySqMatRl(toret,&ongoing);
    int k=1;
    double lambda=1.0;
    for(;k<p;k++) {
        lambda*=((double) (p-k+1))/((double) (p+q-k+1))/((double) k);
        safeCopySqMatRl(multSqMatRl(ongoing,mat),&ongoing);
        safeCopySqMatRl(addSqMatRl(scaleSqMatRl(ongoing,lambda),toret),&toret);
    }
    return toret;
}
struct SqMatRl PadeDRl(struct SqMatRl mat, int p, int q) {
    struct SqMatRl toret;
    initIdSqMatRl(&toret,mat.N);
    struct SqMatRl ongoing;
    safeCopySqMatRl(toret,&ongoing);
    int k=1;
    double lambda=1.0;
    for(;k<p;k++) {
        lambda*=-1*((double) (p-k+1))/((double) (p+q-k+1))/((double) k);
        safeCopySqMatRl(multSqMatRl(ongoing,mat),&ongoing);
        safeCopySqMatRl(addSqMatRl(scaleSqMatRl(ongoing,lambda),toret),&toret);
    }
    return toret;
}
struct SqMatRl expPadeRl(struct SqMatRl mat, int p, int q) {
    struct SqMatRl toret;
<<<<<<< HEAD
    double norm=oneNorm(mat);
=======
    double norm=oneNormRl(mat);
>>>>>>> 8bdac49bb032b17c016452358caefde307581eb3
    int m=1;
    if(norm>0.5) {
        m=1+((int)(norm*2));
    }
    struct SqMatRl matPrime;
    safeCopySqMatRl(scaleSqMatRl(mat,1.0/m),&matPrime);
    safeCopySqMatRl(multSqMatRl(invSqMatRl(PadeDRl(matPrime,p,q)),PadeNRl(matPrime,p,q)),&toret);
    safeCopySqMatRl(powSqMatRl(toret,m),&toret);
    return toret;
}
double * expvKrylovRl(struct SqMatRl mat, double t, double * v, double tau, double minerr) {
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
    struct SqMatRl H;
    initZeroSqMatRl(&H, mat.N+2);
    struct SqMatRl F;
    initZeroSqMatRl(&F,mat.N+2);
    double * w=malloc(mat.N*sizeof(double));
    double * v1=malloc(mat.N*sizeof(double));
    double * p=malloc(mat.N*sizeof(double));
    for(k=0;k<mat.N;k++) {
        *(w+k*sizeof(double))=0;
        for(k1=0;k1<mat.N;k1++) {
            *(w+k*sizeof(double))+=getValRl(mat,k,k1)*(*(v+k1*sizeof(double)));
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
                    *(p+k*sizeof(double))+=getValRl(mat,k,k1)*(*(v1+k*sizeof(double)));
                }
            }
            for(i=0;i<mat.N;i++) {
                tempvar=0;
                for(k=0;k<mat.N;k++) {
                    tempvar+=(*(v1+k*sizeof(double)))*(*(p+k*sizeof(double)));
                }
                setValRl(&H,i,j,tempvar);
                for(k=0;k<mat.N;k++) {
                    *(p+k*sizeof(double))=*(p+k*sizeof(double))-getValRl(H,i,j)*(*(v1+k*sizeof(double)));
                }
            }
            tempvar=0;
            for(k=0;k<mat.N;k++) {
                tempvar+=(*(p+k*sizeof(double)))*(*(p+k*sizeof(double)));
            }
            setValRl(&H,mat.N,mat.N-1,tempvar);
            //if(getValRl(H,mat.N,mat.N)<=tol_abs(mat)) { // TODO WHAT IS THIS ???
            //    happy_breakdown(); // TODO WHAT IS THIS ???
            //}
            for(k=0;k<mat.N;k++) {
                *(v1+k*sizeof(double))=*(p+k*sizeof(double))/getValRl(H,mat.N,mat.N-1);
            }
        }
        setValRl(&H,mat.N+1,mat.N,1.0);
        double delta_tol=.01; // TODO add delta_tol
        errloc=delta_tol+1;
        while(errloc>delta_tol) {
            F=expPadeRl(scaleSqMatRl(H,tau),14,14); // TODO pick exp
            for(k=0;k<mat.N;k++) {
                *(w+k*sizeof(double))=beta*(*(v1+k*sizeof(double)))*getValRl(F,k,1);
            }
            // LOCAL TRUNCATION ERROR ESTIMATE
            err1=beta*abs(getValRl(H,mat.N,mat.N-1)*getValRl(F,1,0));
            err2=beta*matvnorm*abs(getValRl(H,mat.N,mat.N-1)*getValRl(F,2,0));
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
void safeCopySqMat(struct SqMat mat1, struct SqMat * mat2) {
    free(mat2->pValue);
    mat2->N=mat1.N;
    mat2->pValue=malloc(mat2->N*mat2->N*sizeof(struct Complex));
    int k=0;
    for(;k<mat2->N*mat2->N;k++) {
        *(mat2->pValue+k*sizeof(struct Complex))=*(mat1.pValue+k*sizeof(struct Complex));
    }
}
void initZeroSqMat(struct SqMat * mat, int N) {
    free(mat->pValue);
    mat->pValue=malloc(N*N*sizeof(struct Complex));
    mat->N=N;
    int k=0;
    struct Complex zero;
    zero.re=0.0;
    zero.im=0.0;
    for(;k<N*N;k++) {
        *(mat->pValue+k*sizeof(struct Complex))=zero;
    }
}
void initIdSqMat(struct SqMat * mat, int N) {
    free(mat->pValue);
    mat->pValue=malloc(N*N*sizeof(struct Complex));
    mat->N=N;
    int k=0;
    int k1=0;
    struct Complex zero;
    zero.re=0.0;
    zero.im=0.0;
    struct Complex one;
    one.re=1.0;
    one.im=0.0;
    for(;k<N;k++) {
        for(k1=0;k1<N;k1++) {
            if(k1!=k) {
                *(mat->pValue+(k*N+k1)*sizeof(struct Complex))=zero;
            }
            else {
                *(mat->pValue+(k*N+k1)*sizeof(struct Complex))=one;
            }
        }
    }
}
struct Complex getVal(struct SqMat mat, int i, int j) {
    struct Complex toret;
    toret.re=0.0;
    toret.im=0.0;
    if(i>=0 && j>=0 && i<mat.N && j<mat.N) {
        toret=*(mat.pValue+(i*mat.N+j)*sizeof(struct Complex));
    }
    return toret;
}
void setVal(struct SqMat * mat, int i, int j, struct Complex x) {
    if(i>=0 && j>=0 && i<mat->N && j<mat->N) {
        *(mat->pValue+(i*mat->N+j)*sizeof(struct Complex))=x;
    }
}
struct SqMat scaleSqMat(struct SqMatRl mat, struct Complex lambda);
struct SqMat addSqMat(struct SqMatRl mat1, struct SqMatRl mat2);
struct SqMat multSqMatRl(struct SqMat mat1, struct SqMat mat2);
