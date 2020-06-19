clc;
clear;
Cm = 0.01;
dt = 0.01;
t=0:dt:20;

I=0.1;

ENa=55.17;  % mv Na reversal potential
EK=-72.14; % mv K reversal potential
El=-49.42; % mv Leakage reversal potential
gbarNa=1.2; % mS/cm^2 Na conductance
gbarK=0.36; % mS/cm^2 K conductance
gbarl=0.003; % mS/cm^2 Leakage conductance
gbarsyn=40*10^(-2); %S second edit
tau=5; %second edit
Esyn=-70; %mv second edit inhibitory
rate=0.2;%to [ms ] spiking?
V(1)=-60; % Initial Membrane voltage
m(1)=am(V(1))/(am(V(1))+bm1(V(1))); 
% Initial m-value
n(1)=an(V(1))/(an(V(1))+bn(V(1))); 
% Initial n-value
h(1)=ah(V(1))/(ah(V(1))+bh(V(1))); 
% Initial h-value
gsyn(1)=gbarsyn;
Isyn(1)=gsyn*(V(1)-Esyn);
for i=1:length(t)-1
    i;
    %Euler method to find the next m/n/h value
    %rate = rate + 0.015;
    m(i+1)=m(i)+dt*((am(V(i))*(1-m(i)))-(bm1(V(i))*m(i))); 
    n(i+1)=n(i)+dt*((an(V(i))*(1-n(i)))-(bn(V(i))*n(i)));
    h(i+1)=h(i)+dt*((ah(V(i))*(1-h(i)))-(bh(V(i))*h(i)));
    gNa=gbarNa*m(i)^3*h(i);
    gK=gbarK*n(i)^4;
    gl=gbarl;
    %kkk=exp(-(i*10^(-3)-rate)/tau)
    %x=exp(-(i*dt-rate)/tau)
    %y=heaviside(i*dt-rate);
    %gbarsyn*exp(-t(i)/tau)
    gsyn(i)=gbarsyn*exp(-t(i)/tau); %second edit
    INa=gNa*(V(i)-ENa);
    IK=gK*(V(i)-EK);
    Il=gl*(V(i)-El);
    Isyn(i+1)=Isyn(i)+dt*gsyn(i)*(V(i)-Esyn);
    %V(i)
    %dt*gsyn(i)*(V(i)-Esyn)
 % %second edit
    %Euler method to find the next voltage value
    %Isyn(i)
    V(i+1)=V(i)+dt*((1/Cm)*(I-(INa+IK+Il)));%+dt*(1/Cm)*Isyn(i);
end

%Store variables for graphing later
FE=V;
FEm=m;
FEn=n;
FEh=h;
Ie=-Isyn;
G=gsyn;

%Plot the functions
plot(t,FE);
%plot(t,Ie);
%plot(t,G);
hold on
%scatter (t(1:end-1),Ie); %second edit
legend('Forward Euler');
xlabel('Time (ms)');
ylabel('Voltage (mV)');
title('Voltage Change for Hodgkin-Huxley Model');