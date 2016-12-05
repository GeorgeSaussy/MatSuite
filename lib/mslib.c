#include "mslib.h"
#include <stdlib.h>
#include <math.h>
#include <assert.h>
void safeCopySqMatRl(struct SqMatRl mat1, struct SqMatRl * mat2) {
    free(mat2->pValue);
    mat2->N=mat1.N;
    mat2->pValue=(double*) malloc(mat1.N*mat1.N*sizeof(double));
    assert(mat2->pValue!=NULL);
    int k=0;
    for(;k<mat1.N*mat1.N;k++) {
        *(mat2->pValue+k*sizeof(double))=*(mat1.pValue+k*sizeof(double));
    }
}
void initZeroSqMatRl(struct SqMatRl * mat, int N) {
    mat->N=N;
    mat->pValue=(double*) malloc(N*N*sizeof(double));
    assert(mat->pValue!=NULL);
    int k=0;
    for(;k<N*N;k++) {
        *(mat->pValue+k*sizeof(double))=0;
    }
}
void initIdSqMatRl(struct SqMatRl * mat, int N) {
    mat->N=N;
    mat->pValue=(double*) malloc(N*N*sizeof(double));
    assert(mat->pValue!=NULL);
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
    double * tempRow=(double*) malloc(mat->N*sizeof(double));
    assert(tempRow!=NULL);
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
    double f=0;
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
            f=getValRl(ongoing,i,k)/getValRl(ongoing,k,k);
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
    double norm=oneNormRl(mat);
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
    double * w=(double*) malloc(mat.N*sizeof(double));
    assert(w!=NULL);
    double * v1=(double*) malloc(mat.N*sizeof(double));
    assert(v1!=NULL);
    double * p=(double*) malloc(mat.N*sizeof(double));
    assert(p!=NULL);
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
Complex addCmplx(Complex num1, Complex num2) {
    Complex toret;
    toret.re=num1.re+num2.re;
    toret.im=num1.im+num2.im;
    return toret;
}
Complex subCmplx(Complex num1, Complex num2) {
    Complex toret;
    toret.re=num1.re-num2.re;
    toret.im=num1.im-num2.im;
    return toret;
}
Complex multCmplx(Complex num1, Complex num2) {
    Complex toret;
    toret.re=num1.re*num2.re-num1.im*num2.im;
    toret.im=num1.re*num2.im+num1.im*num2.re;
    return toret;
}
Complex divCmplx(Complex num1, Complex num2) {
    Complex toret;
    toret.re=(num1.re*num2.re+num1.im*num2.im)/(num2.re*num2.re+num2.im*num2.im);
    toret.im=(num1.im*num2.re-num1.re*num2.im)/(num2.re*num2.re+num2.im*num2.im);
    return toret;
}
double normCmplx(Complex num) {
    double toret;
    toret=sqrt(num.re*num.re+num.im*num.im);
    return toret;
}
void safeCopySqMat(struct SqMat mat1, struct SqMat * mat2) {
    int k=0;
    if(mat1.N==mat2->N) {
        for(;k<mat2->N*mat2->N;k++) {
            *(mat2->pValue+k*sizeof(Complex))=*(mat1.pValue+k*sizeof(Complex));
        }
    }
    else {
        free(mat2->pValue);
        mat2->N=mat1.N;
        mat2->pValue=(Complex*) malloc(mat2->N*mat2->N*sizeof(Complex));
        assert(mat2->pValue!=NULL);
        for(;k<mat2->N*mat2->N;k++) {
            *(mat2->pValue+k*sizeof(Complex))=*(mat1.pValue+k*sizeof(Complex));
        }
    }
}
void initZeroSqMat(struct SqMat * mat, int N) {
    //free(mat->pValue);
    mat->pValue=(Complex*) malloc(N*N*sizeof(Complex));
    assert(mat->pValue!=NULL);
    mat->N=N;
    int k=0;
    Complex zero;
    zero.re=0.0;
    zero.im=0.0;
    for(k=0;k<N*N;k++) {
        *(mat->pValue+k*sizeof(Complex))=zero;
    }
}
void zeroSqMat(struct SqMat * mat) {
    int N=mat->N;
    int k=0;
    Complex zero;
    zero.re=0.0;
    zero.im=0.0;
    for(k=0;k<N*N;k++) {
        *(mat->pValue+k*sizeof(Complex))=zero;
    }
}
void initIdSqMat(struct SqMat * mat, int N) {
    free(mat->pValue);
    mat->pValue=(Complex*) malloc(N*N*sizeof(Complex));
    assert(mat->pValue!=NULL);
    mat->N=N;
    int k=0;
    int k1=0;
    Complex zero;
    zero.re=0.0;
    zero.im=0.0;
    Complex one;
    one.re=1.0;
    one.im=0.0;
    for(;k<N;k++) {
        for(k1=0;k1<N;k1++) {
            if(k1!=k) {
                *(mat->pValue+(k*N+k1)*sizeof(Complex))=zero;
            }
            else {
                *(mat->pValue+(k*N+k1)*sizeof(Complex))=one;
            }
        }
    }
}
void idSqMat(struct SqMat * out) {
    int k=0;
    int k1=0;
    Complex one;
    one.re=1.0;
    one.im=0.0;
    Complex zero;
    zero.re=0.0;
    zero.im=0.0;
    for(;k<out->N;k++) {
        for(k1=0;k1<out->N;k++) {
            if(k==k1) {
                setVal(out,k,k1,one);
            }
            else {
                setVal(out,k,k1,zero);
            }
        }
    }
}
Complex getVal(struct SqMat mat, int i, int j) {
    Complex toret;
    toret.re=0.0;
    toret.im=0.0;
    if(i>=0 && j>=0 && i<mat.N && j<mat.N) {
        toret=*(mat.pValue+(i*mat.N+j)*sizeof(Complex));
    }
    return toret;
}
void setVal(struct SqMat * mat, int i, int j, Complex x) {
    if(i>=0 && j>=0 && i<mat->N && j<mat->N) {
        *(mat->pValue+(i*mat->N+j)*sizeof(Complex))=x;
    }
}
struct SqMat scaleSqMat(struct SqMat mat, Complex lambda) {
    struct SqMat toret;
    initZeroSqMat(&toret, mat.N);
    int k=0;
    int k1=0;
    for(;k<mat.N;k++) {
        for(k1=0;k1<mat.N;k1++) {
            setVal(&toret,k,k1,multCmplx(getVal(mat,k,k1),lambda));
        }
    }
    return toret;
}
void addSqMat(struct SqMat mat1, struct SqMat mat2, struct SqMat * out) {
    if(mat1.N==mat2.N && mat2.N==out->N) {
        int k=0;
        int k1=0;
        for(;k<mat1.N;k++) {
            for(k1=0;k1<mat1.N;k1++) {
                setVal(out,k,k1,addCmplx(getVal(mat1,k,k1),getVal(mat2,k,k1)));
            }
        }
    }
}
void multSqMat(struct SqMat mat1, struct SqMat mat2, struct SqMat * out) {
    if(mat1.N==mat2.N && mat2.N==out->N) {
        int k=0;
        int k1=0;
        int k2=0;
        Complex value;
        for(;k<mat1.N;k++) {
            for(k1=0;k1<mat2.N;k1++) {
                value.re=0.0;
                value.im=0.0;
                for(k2=0;k2<mat1.N;k2++) {
                    value=addCmplx(value,multCmplx(getVal(mat1,k,k2),getVal(mat2,k2,k1)));
                }
                setVal(out,k,k1,value);
            }
        }
    }
}
void powSqMat(struct SqMat mat, unsigned int j, struct SqMat * out) { // this can be made faster
    struct SqMat ongoing;
    initZeroSqMat(&ongoing, mat.N);
    struct SqMat * workspace;
    initZeroSqMat(workspace, mat.N);
    struct SqMat * tempmat;
    initZeroSqMat(tempmat, mat.N);
    helper_powSqMat(mat,j,out,ongoing,workspace,tempmat);
}
void helper_powSqMat(struct SqMat mat, unsigned int j, struct SqMat * out, struct SqMat ongoing, struct SqMat * workspace, struct SqMat * tempmat) {
    if(j!=0) {
        unsigned int iter=(unsigned int) log2(1.0*j);
        unsigned int diff=(unsigned int) j-pow(2,iter);
        safeCopySqMat(mat,workspace);
        while(iter>0) {
            multSqMat(*workspace,*workspace,tempmat);
            safeCopySqMat(*tempmat,workspace);
            iter--;
        }
        multSqMat(ongoing,*workspace,tempmat);
        safeCopySqMat(*tempmat,&ongoing);
        helper_powSqMat(mat,diff,out,ongoing,workspace,tempmat);
    }
    else {
        safeCopySqMat(ongoing,out);
    }
}
void transSqMat(struct SqMat mat, struct SqMat * out) {
    int k=0;
    int k1=0;
    for(;k<mat.N;k++) {
        for(k1=0;k1<mat.N;k1++) {
            setVal(out,k,k1,getVal(mat,k1,k));
        }
    }
}
void swapRows(struct SqMat * mat, int i, int j) {
    Complex * tempRow=(Complex*) malloc(mat->N*sizeof(Complex));
    assert(tempRow!=NULL);
    int k=0;
    for(;k<mat->N;k++) {
        *(tempRow+(k+i*mat->N)*sizeof(Complex))=getVal(*mat,i,k);
    }
    for(k=0;k<mat->N;k++) {
        setVal(mat,i,k,getVal(*mat,j,k));
    }
    for(k=0;k<mat->N;k++) {
        setVal(mat,i,k,*(tempRow+(k+i*mat->N)*sizeof(Complex)));
    }
}
void invSqMat(struct SqMat mat, struct SqMat * out) {
    int n=mat.N;
    idSqMat(out);
    struct SqMat ongoing;
    safeCopySqMat(mat,&ongoing);
    int k=0;
    int k1=0;
    int i=0;
    int j=0;
    int iMax=0;
    Complex compVal1;
    Complex compVal2;
    Complex f;
    Complex zero;
    zero.re=0.0;
    zero.im=0.0;
    for(;k<n;k++) {
        iMax=k;
        for(k1=k;k1<n;k1++) {
            compVal1=getVal(mat,k1,k);
            compVal2=getVal(mat,iMax,k);
            if(compVal1.re*compVal1.re+compVal1.im*compVal1.im>compVal2.re*compVal2.re+compVal2.im*compVal2.im) {
                iMax=k1;
            }
        }
        // warning: A singular if A[iMax][k] == 0
        swapRows(&ongoing,k,iMax);
        swapRows(out,k,iMax);
        for(i=k+1;i<n;i++) {
            f=divCmplx(getVal(ongoing,i,k),getVal(ongoing,k,k));
            for(j=k+1;j<n;j++) {
                setVal(&ongoing,i,j,subCmplx(getVal(ongoing,i,j),multCmplx(f,getVal(ongoing,k,j))));
                setVal(out,i,j,subCmplx(getVal(*out,i,j),multCmplx(f,getVal(*out,k,j))));
            }
            setVal(&ongoing,i,k,zero);
        }
    }
    // reduce
    for(k=n-1;k>=0;k--) {
        for(i=k-1;i>=0;i--) {
            f=divCmplx(getVal(ongoing,i,k),getVal(ongoing,k,k));
            for(j=0;j<n;j++) {
                setVal(&ongoing,i,j,subCmplx(getVal(ongoing,i,j),multCmplx(f,getVal(ongoing,i,j))));
                setVal(out,i,j,subCmplx(getVal(*out,i,j),multCmplx(f,getVal(*out,i,j))));
            }
        }
    }
}
double oneNorm(struct SqMat mat) {
    double toret=0.0;
    int k=0;
    int k1=0;
    Complex testVal;
    for(;k<mat.N;k++) {
        for(k1=0;k1<mat.N;k1++) {
            testVal=getVal(mat,k,k1);
            if(testVal.re*testVal.re+testVal.im*testVal.im>toret) {
                toret=testVal.re*testVal.re+testVal.im*testVal.im;
            }
        }
    }
    toret=sqrt(toret);
    return toret;
}
void initZeroMat(struct SqMat * mat, int N, int M) {
    mat->N=N;
    mat->M=M;
    mat->pValue=(Complex*)malloc(N*M*sizeof(Complex));
    Complex zero;
    zero->re=0.0;
    zero->im=0.0;
    for(int k=0;k<M*N;k++) {
        *(mat->pValue+k*sizeof(Complex))=zero;
    }
}
Complex getValMat(struct SqMat mat, int i, int j) {
    Complex toret;
    toret=*(mat.pValue+i*mat.N+j*mat.M);
    return toret;
}
void setValMat(struct Matrix * mat, int i, int j, Complex x) {
    *(mat->pValue+i*mat->N+j*mat->M)=x;
}
void multMat(struct Matrix mat1, struct Matrix mat2, struct Matrix * out) {
    ongoing=zero;
    if(mat1.M==mat2.N) {
        int k=0;
        int k1=0;
        int k2=0;
        Complex value;
        for(k=0;k<mat1.N;k++) {
            value.re=0.0;
            value.im=0.0;
            for(k1=0;k1,mat2.M;k2++) {
                value=addCmplx(value,multCmplx(getValMat(mat1,k,k2),getValMat(mst2,k2,k1));
            }
            setValMat(out,k,k1,value);
        }
    }
}
struct Matrix getColMat(struct Matrix mat, int j) {
    struct Matrix toret;
    initZeroMat(&toret,mat.N,j);
    for(int k=0;k<mat.N;k++) {
        setValMat(&toret,k,1,getVal(mat,k,j));
    }
    return toret;
}
struct Matrix transMat(struct Matrix mat) {
    struct Matrix toret;
    initZeroMat(&toret,mat.M,mat.N);
    for(int i=0;i<toret.N;i++) {
        for(int j=0;j<toret.M;j++) {
            setValMat(&toret,i,j,getValMat(mat,j,i));
        }
    }
    return toret;
}
void expvKrylov(double t, struct SqMat A, Matrix v, double tol, int m1, Matrix * w, err, double * hump) {
    // define some constants
    int n=A.N;
    int m=min(m1,30);
    double anorm=oneNorm(A);
    double mxrej=10, btol=1.0e-7;
    double gamma=0.9, delta=1.2;
    double mb=(double)m;
    double  t_out=abs(t), t_new=0, t_now=0;
    int nstep=0;
    double s_error=0, rndoff=anorm*eps; // FIXME hwat is eps?
    int k1=2;
    double xm=1.0/((double)m);
    double normv, beta, fact;
    double s, sgn;
    normv=0.0;
    for(int iter=0;iter<n;iter++) {
        normv+=normCmplx(*(v+iter*sizeof(Complex)));
    }
    normv=sqrt(normv);
    beta=normv;
    fact=pow(((m+1)/exp(1)),(m+1))*sqrt(2*M_PI*(m+1));
    t_new=(1/anorm)*pow(((fact*tol)/(4*beta*anorm)),xm);
    s=pow(10.0,floor(log10(t_new))-1);
    t_new=ceil(t_new/s)*s;
    sgn=1;
    if(t<0) {
        sgn=-1;
    }
    safeCopyMat(v,w);
    *hump=normv;
    while(t_now<t_out) {
        Matrix V;
        Matrix H;
        nstep+=1;
        t_step=min(t_out-t,t_new);
        initZeroMat(&mat,n,m+1);
        initZeroMat(&mat,m+2,m+2);
        for(int iter=0;iter<A.N;iter++) {
            setValMat(&V,iter,0,*(w+iter*sizeof(Complex))/beta);
        }
        struct Matrix p;
        p->pValue=(Complex*)malloc(V.N*sizeof(Complex));
        for(int j=0;j<m;j++) {
            multMat(A,getColMat(V,j));
            for(int i=0;i<j;i++) {
                struct Matrix singleton;
                singleton->pValue=(Complex*)malloc(sizeof(Complex));
                multMat(transMat(getCol(V,i)),p&singleton);
                setValMat(&H,i,j,getValMat(V,))
                // TODO HERE !!!
            }
        }
    }
}
