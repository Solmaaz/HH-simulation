clc;
clear;
dt = 0.01;
t=0:dt:100;
gbarsyn=0.04; %4*10^(-5)S second edit
tau=5.0; %second edit
Esyn=-75;
V(1)=-65;
gsyn(1)=gbarsyn;
Isyn(1)=gsyn*(V(1)-Esyn);
for i=1:length(t)-1
    i;
    gsyn(i)=gbarsyn*exp(-t(i)/tau);
    G(i)=gsyn(i)*(V(i)-Esyn);
end
Gi=G;
plot(t,Gi);