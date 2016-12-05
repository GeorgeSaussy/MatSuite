function fx = IFT_vic( fp,p,dp,np,x,nx)
% Brute force FT
for k=1:nx
    fx(k)=0;
    for j=1:np
        fx(k) = fx(k) + fp(j)*exp(1i*p(j)*x(k));
    end
    fx(k)=fx(k)*dp/sqrt(2*pi);
end

