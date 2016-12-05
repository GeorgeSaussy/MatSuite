% MATLAB functions for tt-SOFT wavepacket propagation
% Chemistry 470a/570a - Introduction to Quantum Chemistry - Fall 2016
% Yale University - Victor S. Batista, Oct 9, 2016.
%
clear all
% Propagation Parameters
nsteps= 30; % # of propagation steps
tol=0.0001; % tolerance for qtt-ultrafast Fourier Transform
m = 1; % system mass
h = 1; % hbar
om= 1;% Frequency
alpha=m*om; % initial width
x0=-1; % initial position
p0= 0; % initial momentum
q=6;   % quantics
nx=2^q;  % # of grid points in coord representation
qshape=[2*ones(1,q)];
np=nx; % # of grid points in momentum representation
L=10; % simulation box
dx= L/nx; % x grid spacing
dp= 2*pi/L; % spacing in momenta space
tau=0.25; % propagation time step
% Build grids for initial state and propagator
for j=1:nx
    x(j)=(j-nx/2)*dx;
    psi(j)=psi0(x(j),p0,x0,alpha);
    V(j)=pot(x(j));
    potp(j)=exp(-1i*V(j)*tau/2);
    pott(j)=V(j);
end
for j=1:np
    if(j<(np/2+1))
        p(j)=(j-1)*dp;
    else
        p(j)=(j-1-np)*dp;
    end
    kinp(j)=exp(-1i*p(j)^2*tau/(2*m));
    kin(j)=p(j)^2/(2*m);
end
% Build quantic tensor trains
psi_tt = tt_tensor(reshape(psi,qshape));
potp_tt = tt_tensor(reshape(potp,qshape));
kinp_tt = tt_tensor(reshape(kinp,qshape));
kin_tt = tt_tensor(reshape(kin,qshape));
pot_tt = tt_tensor(reshape(pott,qshape));

FT_psi_tt =qtt_fft_v(psi_tt,[q],tol,1)/sqrt(2)^(q-4);
rhop_tt = times(conj(FT_psi_tt),FT_psi_tt);
rhox_tt = times(conj(psi_tt),psi_tt);
ENERGY = real(dot(rhop_tt,kin_tt)*dp+dot(rhox_tt,pot_tt)*dx);

% Propagation
for k=1:nsteps
    % soft propagation
    for j=1:nx
        psi(j)=potp(j)*psi(j);
    end
    FTpsi= FT_vic(psi,x,dx,nx,p,np);
    for j=1:np
        FTpsi(j)=kinp(j)*FTpsi(j);
    end
    psi= IFT_vic(FTpsi,p,dp,np,x,nx);
    for j=1:nx
        psi(j)=potp(j)*psi(j);
    end
    % qtt-soft propagation
    psi_tt = times(psi_tt,potp_tt);
    FT_psi_tt =qtt_fft_v(psi_tt,[q],tol,1)/sqrt(2)^(q-4);
    FT_psi_tt = times(FT_psi_tt,kinp_tt);
    psi_tt =qtt_fft_v(FT_psi_tt,[q],tol,-1)*sqrt(2)^(q-4);
    psi_tt = times(psi_tt,potp_tt);
    % Visualization
    % In coordinate representation
    plot(x,ENERGY+abs(full(psi)));
    axis([x(1) x(nx) 0 2*ENERGY]);
    hold on;
    plot(x,ENERGY+abs(full(psi_tt)),'color','r');
        plot(x,ENERGY+0*pott,'color','b');
    plot(x,V,'color','b');
    pause(.1);
    hold off;
    % In momentum representation
    %     plot(p,abs(full(FTpsi)));
    %     hold on;
    %     plot(p,abs(full(FT_psi_tt)),'color','r');
    %     pause(.1);
    %     hold off;
end

