clc;
clearvars;
%[PFf,nx]=PFFun('case9');
%[Qs,MakeJ,Mats]=GetQuads(PFf,nx);
Qs=rand(11,11,10);
nx=10;
QQ=zeros(nx,nx,nx);
for i=1:nx
    Qt=rand(nx,nx);
    QQ(:,:,i)=Qt'*Qt;
end
%QQ=Qs(2:end,2:end,:);
M=zeros(nx);
%for i=1:nx
%    M(i,:)=Qs(1,2:end,i);
%end
M=-eye(nx);

cvx_begin
variable Lam(nx,nx);
for i=1:nx
    Q=QQ(:,:,1)*Lam(i,1);
    for j=2:nx
        Q=Q+QQ(:,:,j)*Lam(i,j);  
    end
  %  Q>=0;
    Q==semidefinite(nx,nx);
end
J=Lam*M;
2*diag(J)>=sum(abs(J),2)+1e-2;
J(eye(nx)==0)<=0;
cvx_end
