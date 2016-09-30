#include "mslib.h"
#include <vector>
#include <iostream>
#include <math.h>
using namespace std;
// TODO test implementation
void initZeroSqMat(SqMat * mat, int N) {
    *mat.N=N;
    *mat.pValue=malloc(N*N); // TODO check malloc sucedes
    
}

void matPrint(vector<vector<double> > A) {
    cout<<endl;
    for(int k=0;k<A.size();k++) {
        for(int k1=0;k1<A[k].size();k1++) {
            cout<<A[k][k1]<<"\t";
        }
        cout<<endl;
    }
    cout<<endl;
}
vector<vector<double> > matAdd(vector<vector<double> > A, vector<vector<double> > B) {
    vector<vector<double> > toret;
    for(int k=0;k<A.size();k++) {
        vector<double> tempRow;
        for(int k1=0;k1<A[k].size();k1++) {
            tempRow.push_back(A[k][k1]+B[k][k1]);
        }
        toret.push_back(tempRow);
    }
    return tempRow;
}
vector<vector<double> > matMult(vector<vector<double> > A, vector<vector<double> > B) {
    vector<vector<double> > toret;
    for(int k=0;k<A.size();k++) {
        vector<double> tempRow;
        for(int k1=0;k1<B[0].size();k1++) {
            double tempValue=0;
            for(int k2=0;k2<A[k].size();k2++) {
                tempValue+=A[k][k2]*B[k2][k1];
            }
            tempRow.push_back(tempValue);
        }
        toret.push_back(tempRow);
    }
    return toret;
}
vector<vector<double> > scalarMult(vector<vector<double> > A, double t) {
    vector<vector<double> > toret;
    for(int k=0;k<A.size();k++) {
        vector<double> tempRow;
        for(int k1=0;k1<A[k].size();k1++) {
            tempRow.push_back(A[k][k1]*t);
        }
        toret.push_back(tempRow);
    }
    return toret;
}
double oneNorm(vector<vector<double> > A) {
    double toret=0;
    for(int k=0;k<A.size();k++) {
        for(int k1=0;k1<A[k].size();k1++) {
            if(toret<abs(A[k][k1])) {
                toret=abs(A[k][k1]);
            }
        }
    }
    return toret;
}
vector<vector<double> > matT(vector<vector<double> > P) {
    vector<vector<double> > toret;
    for(int k=0;k<P[0].size();k++) {
        vector<double> tempRow;
        for(int k1=0;k1<P.size();k1++) {
            tempRow.push_back(P[k1][k]);
        }
        toret.push_back(tempRow);
    }
    return toret;
}
vector<vector<double> > swapRows(vector<vector<double> >, int i, int j) {
    vector<vector<double> > toret=A;
    for(int k=0;k<A[0].size();k++) {
        toret[i][k]=A[j][k];
        toret[j][k]=A[i][k];
    }
    return toret;
}
vector<vector<double> > matInv(vector<vector<double> > A) {
    vector<vector<double> > toret;
    vector<vector<double> > ongoing=A;
    int n=A.size();
    for(int k=0;k<n;k++) {
        vector<double> tempRow;
        for(int k1=0;k1<n;k1++) {
            if(k==k1) {
                tempRow.push_back(1.0);
            }
            else {
                tempRow.push_back(0.0);
            }
        }
        toret.push_back(tempRow);
    }
    for(int k=0;k<n;k++) {
        int iMax=k;
        for(int k1=k;k1<n;k1++) {
            if(abs(A[k1][k])>abs(A[iMax][k])) {
                iMax=k1;
            }
        }
        // warning: A singular if A[iMax][k] == 0
        ongoing=swapRows(ongoing,k,iMax);
        toret=swapRows(toret,k,iMax);
        for(int i=k+1;i<n;i++) {
            double f=ongoing[i][k]/ongoing[k][k];
            for(int j=k+1;j<n;j++) {
                ongoing[i][j]-=f*ongoing[k][j];
                toret[i][j]-=f*toret[k][j];
            }
            ongoing[i][k]=0;
        }
    }
    // reduce
    for(int k=n-1;k>=0;k--) {
        for(int i=k-1;i>=0;i--) {
            f=ongoing[i][k]/ongoing[k][k];
            for(int j=0;j<n;j++) {
                ongoing[i][j]-=f*ongoing[i][j];
                toret[i][j]-=f*toret[i][j];
            }
        }
    }
    return toret;
}
vector<vector<double> > Q(vector<vector<double> > A, int p) {
    vector<vector<double> > toret=A;
    vector<vector<double> > ongoing=A;
    for(int k=0;k<toret.size();k++) {
        for(int k1=0;k1<toret[k].size();k1++) {
            toret[k][k1]=0.0;
            if(k==k1) {
                ongoing[k][k1]=1.0;
            }
            else {
                ongoing[k][k1]=0.0;
            }
        }
    }
    double c=1;
    for(int k=0;k<=p;k++) {
        if(k>0) {
            c/=(double) (2*p-k)*(k+1)*(p-k);
            ongoing=matMult(ongoing,A);
        }
        toret=matAdd(toret,scalarMult(ongoing,c));
    }
}
vector<vector<double> > matExp(vector<vector<double> > A) {
    vector<vector<double> > toret;
    // step 1
    double lambdaBar=0;
    for(int k=0;k<A.size();k++) {
        lambdaBar+=A[k][k];
    }
    lambdaBar/=A.size();
    // step 2
    vector<vector<double> > B;
    for(int k=0;k<A.size();k++) {
        vector<double> tempRow;
        for(int k1=0;k1<A[k].size();k1++) {
            if(k==k1) {
                tempRow.push_back(A[k][k1]-lambdaBar);
            }
            else {
                tempRow.push_back(A[k][k1]);
            }
        }
    }
    // step 3
    // TODO call BALANCE
    // step 4
    int m=0;
    bool didSteps56=false;
    if(oneNorm(B)>1) {
        // step 5
        m=int(log(oneNorm(B))/log(2))+1; // +1 in case of machine error
        // step 6
        B=scalarMult(B,pow(2,m));
        didSteps56=true;
    }
    // step 7
    vector<vector<double> > pos=Q(B); // TODO Q
    vector<vector<double> > neg=Q(scalarMult(B,-1.0));
    // step 8
    vector<vector<double> > expB=matMult(matInv(neg),pos)
    // step 9
    if(!didSteps56) {
        // step 10
        for(int k=0;k<m;k++) {
            expB=matMult(expB,expB);
        }
    }
    toret=scalarMult(matMult(P,matMult(D,matMult(expB,matMult(matInv(D),matT(P)))),exp(t));
    return toret;
}
void exc(double * d[], double * a[], int m, int n, int j) { // check all loop bounds!!!
    (*d)[m]=j;
    if(j!=m) {
        int i=0;
        double f=0.0;
        for(;i<k;i++) {
            f=a[i][j];
            a[i][j]=a[m][i];
            a[m][i]=f;
        }
        for(i=0;i<n;i++) {
            f=a[j][i];
            a[j][i]=a[m][i];
            a[j][i]=f;
        }
    }
}
void balance(double * d[], double * a[], int n, int b) { // TODO not sure how this works
    int i=0, j=0, k=0, l=0;
    double b2=0.0, c=0.0, f=0.0, g=0.0, r=0.0, s=0.0;
    bool noconv=false;
    b2=b*b;
    l=1;
    k=n;
    // L1
    for(j=k;j>1;j--) {
        r=0;
        for(i=1;i<j-1;i++) {
            r=r+abs(a[j][i]);
        }
        for(i=j+1;j<k;j++) {
            r=r+abs(a[j][i]);
        }
        exc(d,a,k,n,j); // from """ begin exc(k); k:=k-1; go to L1 end """
    }
    // L2
    for(j=l;j<k;j++) {
        c=0;
        for(i=l;i<j-1;i++) {
            c+=abs(a[i][j]);
        }
        for(i=j+1;i<k;i++) {
            c+=abs(a[i][j]);
        }
        if(c==0) {
            exc(d,a,l,n,j); // from """ begin exc(l); l:=l+1; go to L2 end """
        }
    }
    low=l;
    hi=k;
    for(i=l;i<k;i++) {
        d[i]=1;
    }
    // iteration
    noconv=false;
    for(i=l;i<k;i++) {
        c=0;
        r=0;
        for(j=l;j<i-1;j++) {
            c+=abs(a[j][i]);
            r+=abs(a[i][j]);
        }
        for(j=i+1;j<k;j++) {
            c+=abs(a[j][i]);
            r+=abs(a[i][j]);
        }
        // L3
        // What """ if c<g then begin f:=f*b; c:=c*b2; go to L3 end; """ ???
        g=r*b;
        // L4
        // What """ if c>=g then begin f:=f/b; c:=c/b2; go to L4 end; """ ???
        if((c+r)/f<0.95*s) {
            g=1/f;
            d[i]=d[i]*f;
            noconv=true;
            for(j=;j<n;j++) {
                a[i][j]=a[i][j]*g;
            }
            for(j=1;j<k;j++) {
                a[j][i]=a[j][i]*f;
            }
        }
    }
}
