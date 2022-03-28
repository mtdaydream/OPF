function [Qs,MakeJ,Mats]=GetQuads(func,nn)
%   func  -   Vector valued function, each component being a quadratic
%   function of input
%   nn  -   Dimension of input to func

%   Qs - (1+nn) x (1+nn) x size(func,1) tensor such that
%   [func(x)]_i=[1;x]'*Qs(:,:,i)*[1;x]

%   MakeJ - Function that returns Jacobian of func given input x

%   Mats - size(func,1) x nn x (1+nn) tensor such that
%   MakeJ(x) =  Mats(:,:,1)+\sum_{i=1}^nn Mats(:,:,i+1)x_i



FF  =   func(zeros(nn,1));
nf  =   length(FF);
Qs  =   zeros(nn+1,nn+1,nf);

for k=1:nf
    Qs(1,1,k)   =   FF(k);
end

for ii=2:(nn+1)
    FF=func(bsol(ii-1));
    FFa=func(2*bsol(ii-1));
    for k=1:nf
        Qs(ii,ii,k)   =   (FFa(k)-2*FF(k)+Qs(1,1,k))/2;
        Qs(1,ii,k)    =   (FF(k)-Qs(1,1,k)-Qs(ii,ii,k))/2;
        Qs(ii,1,k)    =   Qs(1,ii,k);
    end
end

for ii=2:(nn+1)
    for jj=(ii+1):(nn+1)
        FF=func(bsol(ii-1)+bsol(jj-1));
        for k=1:nf
            Qs(ii,jj,k)   =   (FF(k)-2*(Qs(1,ii,k)+Qs(1,jj,k))-Qs(1,1,k)-Qs(ii,ii,k)-Qs(jj,jj,k))/2;
            Qs(jj,ii,k)   =   Qs(ii,jj,k);
        end
    end
end

Mats        =   zeros(nf,nn,nn+1);
Mats(:,:,1) =   Jac(zeros(nn,1));
for ii=1:nn
    Mats(:,:,ii+1)  =   Jac(sparse(ii,1,1,nn,1))-Mats(:,:,1);
end
    
    function J=Jac(x)
    J   =   zeros(nf,nn);
    for it=1:nf
        J(it,:) =   (2*(Qs(2:end,2:end,it)*x+Qs(1,2:end,it)'))';
    end
    end

 function res=bsol(k)
        res=zeros(nn,1);
        if((k>=1)&&(k<=nn))
            res=sparse(k,1,1,nn,1);
        end
  end
    
    MakeJ=@Jac;
    
end
