clc;
clear;
Cm = 1;%microFarad
dt = 0.01;
t=0:dt:100;
f=10;
I1=2;

ENa=55.17;  % mv Na reversal potential
EK=-72.14; % mv K reversal potential
El=-49.42; % mv Leakage reversal potential
gbarNa=120; % mS/cm^2 Na conductance
gbarK=36; % mS/cm^2 K conductance
gbarl=0.3; % mS/cm^2 Leakage conductance

V(1)=-60; % Initial Membrane voltage
m(1)=am(V(1))/(am(V(1))+bm1(V(1))); 
% Initial m-value
n(1)=an(V(1))/(an(V(1))+bn(V(1))); 
% Initial n-value
h(1)=ah(V(1))/(ah(V(1))+bh(V(1))); 
% Initial h-value
%I=zeros(length(t));
t1=10;
t2=40;
I2=whitegaussiannoise;
%I(1)=0;
%creating white noise
%L=5000;
%L=5000
%mu=0;
%Sigma=5;

%I2=Sigma+randn(L,1)+mu;
%I=I1+I2;
%y=wgn(1000,1,0)

%I=8*sin(2*pi*f*t/1000);
%I=wgn(50,1,0)
for i=1:length(t)-1
    
    %Euler method to find the next m/n/h value
     % I(i);pulse
    m(i+1)=m(i)+dt*((am(V(i))*(1-m(i)))-(bm1(V(i))*m(i))); 
    n(i+1)=n(i)+dt*((an(V(i))*(1-n(i)))-(bn(V(i))*n(i)));
    h(i+1)=h(i)+dt*((ah(V(i))*(1-h(i)))-(bh(V(i))*h(i)));
    gNa=gbarNa*m(i)^3*h(i);
    gK=gbarK*n(i)^4;
    gl=gbarl;
    INa=gNa*(V(i)-ENa);
    IK=gK*(V(i)-EK);
    Il=gl*(V(i)-El);
    %Euler method to find the next voltage value
    V(i+1)=V(i)+(dt)*((1/Cm)*(I1-(INa+IK+Il)))+sqrt(dt)*(1/Cm)*I2(i);
end

%Store variables for graphing later
FI=I2;
FE=V;
FEm=m;
FEn=n;
FEh=h;

%Plot the functions
plot(t,FE);
legend('Forward Euler');
xlabel('Time (ms)');
ylabel('voltage (mv)');
title('Voltage Change for Hodgkin-Huxley Model');