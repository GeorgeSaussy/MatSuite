% MATLAB implementation of test 3 to for speed 
lowdim=1000;
hidim=10000;
skip=50;
iter=1000;
for dimiter=1:(hidim-lowdim)/skip
    dim=lowdim+skip*dimiter;
    H=eye(dim);
    for k=1:dim
        if mod(k,2)==0
            H(k,k)=-1;
        end;
    end;
    avg=0;
    for i=1:iter 
        v=zeros(dim,1);
        w=zeros(dim,30);
        for k=1:dim
            v(k)=rand+1j*rand;
        end;
        tic;
        w(:,1)=H*v;
        for k=2:30
            w(:,k)=H*w(:,k-1);
        end;
        avg=avg+toc;
    end;
    avg=avg/iter;
    fprintf('%i\t%f\n',dim,avg);
end;
