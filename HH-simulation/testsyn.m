clc;
clear;
Cm = 0.01;
dt = 0.04;
t=0:dt:25;
gbarsyn=40; %ps econd edit
tau=5*10^(-3); %second edit
Esyn=-75; %mv second edit inhibitory
rate=2*10^(-3);%to [ms ] spiking?
for i=1:length(t)-1
    gsyn=gbarsyn*exp(-(i-rate)/tau)*heaviside(i-rate); %second edit
Isyn=gsyn*(V(i)-Esyn); %second edit
end
Ie=I;
plot (t,Ie); %second edit