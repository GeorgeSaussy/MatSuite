function fp = FT_vic( fx,x,dx,nx,p,np)
% Brute force FT
for k=1:np
    fp(k)=0;
    for j=1:nx
        fp(k) = fp(k) + fx(j)*exp(-1i*x(j)*p(k));
    end
    fp(k)=fp(k)*dx/sqrt(2*pi);
end

