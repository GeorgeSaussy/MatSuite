#include "mslib.h"
#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "time.h"
#include "assert.h"

void ms_HamXVec(Hamiltonian * H, Cmplx * v, Cmplx * w) {
    Cell * iter=H->ptr;

    while((unsigned int)(iter-H->ptr)<H->len) {
        (w[iter->row]).re+=(iter->val).re*(v[iter->col]).re-(iter->val).im*(v[iter->col]).im;
        (w[iter->row]).im+=(iter->val).re*(v[iter->col]).im+(iter->val).re*(v[iter->col]).im;

        iter++;
    }
}
void ms_GenKrylovBasis(Hamiltonian * H, Cmplx * v, Cmplx ** w, unsigned int m) {
    unsigned int count=1;

    ms_HamXVec(H, v, w[0]);
    while(count<m) {
        ms_HamXVec(H, w[count-1], w[count]);
        count++;
    }
}
Hamiltonian  ms_generateRandomHamiltonian(int sp, int dim) {
    Hamiltonian ret;
    // first make (dim / sp) groupings of sp rows, done with Knuths permtation shuffle
    int numberOfGroups=dim/sp;
    int groups[numberOfGroups][sp];
    int basearray[dim];
    int k, k1, k2, j, secondval; // for ANSI compliance
    for(k=0;k<dim;k++) {
        basearray[k]=k;
    }
    for(k=0;k<dim;k++) {
        j=rand()%(dim-k);
        secondval=basearray[k+j];
        basearray[k+j]=basearray[k];
        basearray[k]=secondval;
    }
    j=0;
    for(k=0;k<numberOfGroups;k++) {
        for(k1=0;k1<sp;k1++) {
            groups[k][k1]=basearray[j];
            j++;
        }
    }
    // generate (dim / sp) Hamiltonians
    Cmplx ham[numberOfGroups][sp][sp];
    Cmplx vec[sp];
    Cmplx one, zero;
    double norm=0;
    one.re=1;
    one.im=0;
    zero.re=0;
    zero.im=0;
    for(k=0;k<numberOfGroups;k++) {
        for(k1=0;k1<sp;k1++) {
            vec[k1].re=((double)rand())/((double)RAND_MAX);
            vec[k1].im=((double)rand())/((double)RAND_MAX);
        }
        norm=0;
        for(k1=0;k1<sp;k1++) {
            norm+=sqrt(vec[k1].re*vec[k1].re+vec[k1].im*vec[k1].im);
        }
        norm=sqrt(norm);
        for(k1=0;k1<sp;k1++) {
            vec[k1].re/=norm;
            vec[k1].im/=norm;
        }
        for(k1=0;k1<sp;k1++) {
            for(k2=0;k2<sp;k2++) {
                if(k1==k2) {
                    ham[k][k2][k1].re=1-2*vec[k2].re*vec[k1].re;
                    ham[k][k2][k1].re=-2*vec[k2].im*vec[k1].im;
                }
                else {
                    ham[k][k2][k1].re=-2*vec[k2].re*vec[k1].re;
                    ham[k][k2][k1].re=-2*vec[k2].im*vec[k1].im;
                }
            }
        }
    }
    // apply them to the row groups
    ret.N=dim;
    ret.len=sp*sp*numberOfGroups+(dim-sp*numberOfGroups);
    ret.ptr=(Cell*)malloc(ret.len*sizeof(Cell));
    j=0;
    for(k=0;k<numberOfGroups;k++) {
        for(k1=0;k1<sp;k1++) {
            for(k2=0;k2<sp;k2++) {
                ((ret.ptr)[j]).row=groups[k][k1];
                ((ret.ptr)[j]).col=groups[k][k2];
                ((ret.ptr)[j]).val=ham[k][k1][k2];
                j++;
            }
        }
    }
    for(j=0;j<ret.N-sp*numberOfGroups;j++) {
        ((ret.ptr)[j]).row=basearray[sp*numberOfGroups+j];
        ((ret.ptr)[j]).col=((ret.ptr)[j]).row;
        ((ret.ptr)[j]).val.re=1;
        ((ret.ptr)[j]).val.im=1;
    }
    return ret;
}
double vec_norm(Cmplx * p, int dim) {
    double ret=0;
    for(int k=0;k<dim;k++) {
        ret+=p[k].re*p[k].re+p[k].im*p[k].im;
    }
    ret=sqrt(ret);
    return ret;
}
void ms_expm(Cmplx ** H, Cmplx ** F, Cmplx ** S, Cmplx * p, int dim, int iter) {
    Cmplx zero; zero.re=0.0; zero.im=0.0;
    int k, k1, k2, i;
    double fac;

    for(k=0;k<iter;k++) {
        for(k1=0;k1<dim;k1++) {
            F[k][k1]=zero;
            S[k][k1]=zero;
            if(k==k1) {
                F[k][k1].re=1.0;
                S[k][k1].re=1.0;
            }
        }
    }
    fac=1;
    for(int i=1;i<iter;i++) {
        fac*=i;
        if(i==1) {
            for(k=0;k<dim;k++) {
                for(k1=0;k1<dim;k1++) {
                    S[k][k1]=H[k][k1];
                }
            }
        }
        else {
            for(k=0;k<dim;k++) {
                for(k1=0;k1<dim;k1++) {
                    p[k1]=zero;
                }
                for(k1=0;k1<dim;k1++) {
                    for(k2=0;k2<dim;k2++) {
                        p[k1].re+=H[k1][k2].re*S[k][k2].re-H[k1][k2].im*S[k][k2].im;
                        p[k1].im+=H[k1][k2].re*S[k][k2].im+H[k1][k2].im*S[k][k2].re;
                    }
                }
                for(k1=0;k1<dim;k1++) {
                    S[k][k1]=p[k1];
                }
            }
        }
        for(k=0;k<dim;k++) {
            for(k1=0;k1<dim;k1++) {
                F[k][k1].re+=S[k][k1].re/fac;
                F[k][k1].im+=S[k][k1].im/fac;
            }
        }
    }
}
void ms_expv(double t, Hamiltonian * A, double normA, Cmplx * v, int m, Cmplx * w,
    Cmplx ** V, Cmplx ** H, Cmplx ** F, Cmplx ** S, Cmplx * p, Cmplx * workspace) {
    // constants
    double btol=0.0000001;
    double tol=btol;
    const double eps=pow(10.0,-15); // fudged machine tolerance
    const int exp_iter=5;

    // reserving some space
    const int dim=A->N;
    Cmplx zero; zero.re=0.0; zero.im=0.0;
    Cmplx one; one.re=0.0; one.im=0.0;
    int K1=2, mxrej=10, mb=m, sgn=abs(t)/t, nstep=0;;
    int ireject, mx;
    double s, phi1, phi2, t_step, t_out, t_new, avnorm, normv=1.0, beta=1.0;
    double delta=1.2, xm=1.0/m, gamma=0.9, hump=1.0, t_now=0.0, err;
    double anorm=normA, rndoff=anorm*eps, s_error=0.0;
    double err_loc; // local error


    for(int k=0;k<dim;k++) {
        w[k]=v[k];
    }
    hump=normv;
    while(t_now<t_out) {
        nstep=nstep+1;
        t_step=t_out-t_now;
        if(t_new<t_step) t_step=t_new;
        // V=zeros(n,m+1);
        for(int k=0;k<m+1;k++) {
            for(int k1=0;k1<dim;k1++) {
                V[k][k1]=zero;
            }
        }
        // H=zeros(m+2,m+2);
        for(int k=0;k<m+2;k++) {
            for(int k1=0;k1<m;k1++) {
                H[k][k1]=zero;
            }
        }
        for(int j=1;j<m;j++) {
            ms_HamXVec(A, V[j], p);
            for(int i=0;i<j;i++) {
                // H[i][j]=V[:][i]''*p;
                H[i][j]=zero;
                for(int k=0;k<dim;k++) {
                    H[i][j].re+=V[i][k].re*p[k].re-V[i][k].im*p[k].im;
                    H[i][j].im+=V[i][k].re*p[k].im+V[i][k].im*p[k].re;
                }

                // p=p-H[i][j]*V[:][i];
                for(int k=0;k<dim;k++) {
                    p[k].re-=H[i][j].re*V[k][i].re-H[i][j].im*V[i][k].im;
                    p[k].im-=H[i][j].re*V[k][i].im+H[i][j].im*V[i][k].re;
                }
            }
            s=vec_norm(p, dim);
            if(s<btol) {
                K1=0;
                mb=j;
                t_step=t_out-t_new;
                break;
            }
            H[j+1][j].re=s;
            H[j+1][j].im=0.0;
            // V[:][j+1]=1/s*p;
            for(int k=0;k<dim;k++) {
                V[j+1][k].re=1/s*p[k].re;
                V[j+1][k].im=1/s*p[k].im;
            }
        }
        if(K1!=0) {
            H[m+1][m+1]=one;
            // avnorm=vec_norm(A*V(:,m+1));
    	    ms_HamXVec(A,V[m],workspace);
    	    avnorm=vec_norm(workspace,dim);
        }
        ireject=0;
        while(ireject<=mxrej) {
            mx=mb+K1;
            // F=expm(sgn*t_step*H(1:,mx,1:mx));
            ms_expm(H,F,S,workspace,mx,exp_iter);
            if(K1==0) {
                err_loc=btol;
                break;
            }
            else {
                phi1 = abs(beta)*sqrt(F[0][m].re*F[0][m].re+F[0][m].im*F[0][m].im);
                phi2 = abs(beta*avnorm)*sqrt(F[0][m+1].re*F[0][m+1].re+F[0][m+1].im*F[0][m+1].im);
                if(phi1>10*phi2) {
                    err_loc=phi2;
                    xm=1.0/m;
                }
                else if(phi1>phi2) {
                    err_loc = (phi1*phi2)/(phi1-phi2);
                    xm = 1.0/m;
                }
                else {
                    err_loc = phi1;
                    xm = 1.0/(m-1);
                }
            }
            if(err_loc<=delta*t_step*tol) {
                break;
            }
            else {
                t_step=pow(gamma*t_step*(t_step*tol/err_loc),xm);
                s=pow(10,(floor(log10(t_step))-1));
                t_step=ceil(t_step/s)*s;
                if(ireject==mxrej) {
                    fprintf(stderr,"The requested tolerance is too high.\n");
                }
                ireject+=1;
            }
        }
        mx = mb;
        if(K1>1) mx+=K1;
        // w = V(:,1:mx)*(beta*F(1:mx,1));
        for(int k=0;k<dim;k++) {
            w[k]=zero;
            for(int k1=0;k1,mx;k1++) {
                w[k].re+=beta*(F[1][k1].re*V[k][k1].re-F[1][k1].im*V[k][k1].im);
                w[k].im+=beta*(F[1][k1].re*V[k][k1].im+F[1][k1].im*V[k][k1].re);
            }
        }
        beta=vec_norm(w,dim);
        if(beta>hump) hump=beta;
        t_now=t_now+t_step;

        t_new=pow(gamma*t_step*(t_step*tol/err_loc),xm);
        s=pow(10,(floor(log10(t_new))-1));
        t_new=ceil(t_new/s)*s;

        if(rndoff>err_loc) err_loc=rndoff;
        s_error=s_error+err_loc;
    }
    err=s_error;
    hump=hump/normv;
}
