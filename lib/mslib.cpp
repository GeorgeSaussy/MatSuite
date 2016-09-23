#include <vector>
#include <iostream>
#include <math.h>
using namespace std;
// TODO test implementation
// TODO implement BALANCE

/**
    Prints a matrix to std::out
    @param A a matrix
    @return void
*/
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
/**
    Function to add two matricies to one another
    @param A the first summand
    @param B the second summand
    @return the sum A + B
    @warning demensions of A and B must match
*/
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
/**
    Function to perform matix multiplication
    @param A first matrix factor
    @param B seond matrix factor
    @return the product A*B
    @warning the width of A must equal the height of B
*/
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
/**
    Function to perform scalar multiplication on a matrix
    @param A a sqare matrix
    @param t a scalar
    @return the value A*t
*/
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
/**
    Function to calculate the one norm of a matrix
    @param A a square matrix
    @return the one norm of A
*/
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
/**
    Function to calculate the traspose of a retangular matrix
    @param P a matrix
    @return the transpose of P
*/
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
/**
    Function to swap the rows of a matrix
    @param A a matrix
    @param i the first row index to swap
    @param j the second row index to swap
    @return a matrix with rows i and j interchanged
    @warning must 0 <= i,j < A.size()
*/
vector<vector<double> > swapRows(vector<vector<double> >, int i, int j) {
    vector<vector<double> > toret=A;
    for(int k=0;k<A[0].size();k++) {
        toret[i][k]=A[j][k];
        toret[j][k]=A[i][k];
    }
    return toret;
}
/**
    Function to calculate the inverse of a matrix
    @param A a square invertable matrix
    @return the inverse of A
    @warning A cannot be singular
*/
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
/**
    Calculate the Q function for Pade's approximation
    @param A a square matrix
    @param p degree of polynomial approximation
    @return a square matrix given by (2.1) in Ward 1977
*/
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
/**
    First draft function to calculate exp(At), where A is a matrix and t is a
    scalar, following Ward 1977
    @param A a square matrix
    @param t a scalar
    @return exp(At) as a matix
*/
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