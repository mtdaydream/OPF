function [PFf,nx]=PFFun(casename)
%MATPOWER STUFF%
define_constants;
%compute parameters%
mpc     =   loadcase(casename);
pq      =   find(mpc.bus(:,BUS_TYPE)==1);
pv      =   find(mpc.bus(:,BUS_TYPE)==2);
sb      =   find(mpc.bus(:,BUS_TYPE)==3);
Y       =   makeYbus(mpc);
n       =   size(Y,1);
nsb     =   setdiff(1:n,sb);


%get base point%
try
    res     =   runopf(mpc);
catch
    res     =   runpf(mpc);
end
V0 = res.bus(:,VM).*exp(1i*res.bus(:,VA)*pi/180);
s0 = V0.*conj(Y*V0);

define_constants;
mpc.bus(:,PD) =   -real(s0)*mpc.baseMVA;
mpc.bus(:,QD) =   -imag(s0)*mpc.baseMVA;
mpc.gen(:,PG) =   0;
mpc.gen(:,QG) =   0;
mpopt=mpoption();
mpopt.verbose=0;
res =   runpf(mpc,mpopt);
Y   =   makeYbus(mpc);
V0  =   res.bus(:,VM).*exp(1i*res.bus(:,VA)*pi/180);

PFf =    @(anyx) PFF(anyx,Y,V0,s0,pv,pq);
nx  =   (length(V0)-1)*2;

  
end

function fun=PFF(x,Y,V0,s0,pv,pq)
V   =   V0;
n   =   length(pv)+length(pq);
V([pv;pq])  =   V0([pv;pq])+x(1:n)+1i*x(n+1:end);
s   =   V.*conj(Y*V)-s0;
fun =   [real(s([pv;pq]));imag(s(pq));abs(V(pv)).^2-abs(V0(pv)).^2];
end
    