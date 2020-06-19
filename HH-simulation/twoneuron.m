clc;
clear;
Cm = 1;%microFarad
dt = 0.01;
t=0:dt:100;
f=10;
I1=10;

ENa=55.17;  % mv Na reversal potential
EK=-72.14; % mv K reversal potential
El=-49.42; % mv Leakage reversal potential
Esyn=-70;
gbarNa=120; % mS/cm^2 Na conductance
gbarK=36; % mS/cm^2 K conductance
gbarl=0.3; % mS/cm^2 Leakage conductance
gsyn=20;
I2=5*whitegaussiannoise;
td=10; %ms decay time
tr=5; %ms rise time
tl=0; % ms latency time
mt=10; %ms membrane time cte
V(1)=-60; % Initial Membrane voltage
m(1)=am(V(1))/(am(V(1))+bm1(V(1))); 
% Initial m-value
n(1)=an(V(1))/(an(V(1))+bn(V(1))); 
% Initial n-value
h(1)=ah(V(1))/(ah(V(1))+bh(V(1))); 
% Initial h-value

V1(1)=-60; % Initial Membrane voltage
m1(1)=am1(V1(1))/(am1(V1(1))+bm2(V1(1))); 
% Initial m-value
n1(1)=an1(V1(1))/(an1(V1(1))+bn1(V1(1))); 
% Initial n-value
h1(1)=ah1(V1(1))/(ah1(V1(1))+bh1(V1(1))); 
% Initial h-value
%s(1)=am(V1(1))/(am(V1(1))+bm1(V1(1)));
Isyn(1)=0;
spiketime = 0;
I0 = 0;
spiked = false;

for i=1:length(t)-1
    if(V(i)<0)
    neuron_spiking=false;
    end
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
    
    m1(i+1)=m1(i)+dt*((am1(V1(i))*(1-m1(i)))-(bm2(V1(i))*m1(i))); 
    n1(i+1)=n1(i)+dt*((an1(V1(i))*(1-n1(i)))-(bn1(V1(i))*n1(i)));
    h1(i+1)=h1(i)+dt*((ah1(V1(i))*(1-h1(i)))-(bh1(V1(i))*h1(i)));
    
    gNa1=gbarNa*m1(i)^3*h1(i);
    gK1=gbarK*n1(i)^4;
    gl=gbarl;
    if (V(i)>0 && neuron_spiking==false)
        i
        spiketime=i*dt
        neuron_spiking=true
        spiked = true
        I0 = Isyn(i);
        Isyn(i+1) = Isyn(i);
    else
    	time = i*dt - spiketime;
    	if ((time < 2*(td+tr+tl)) && spiked)
    		s=((mt/(td-tr))*(exp(-(time-tl)/td)-exp(-(time-tl)/tr)));
        	Isyn(i+1)=gsyn*s*(V1(i)-Esyn)*dt+I0;
        else
        	%spike_time=i*dt;
        	Isyn(i+1)=I0*exp((spiketime+tl+tr+td-(i*dt))/(20*td));
        end
    end
    INa1=gNa1*(V1(i)-ENa);
    IK1=gK1*(V1(i)-EK);
    Il1=gl*(V1(i)-El);
    %Euler method to find the next voltage value
    V1(i+1)=V1(i)+dt*(1/Cm)*Isyn(i)+(dt)*((1/Cm)*(0.0-(INa1+IK1+Il1)));
end

%Store variables for graphing later
FI=Isyn;
FE=V;
FE1=V1;
FEm=m;
FEn=n;
FEh=h;

%Plot the functions
figure();
plot(t,FE1);
legend('Forward Euler');
xlabel('Time (ms)');
ylabel('voltage (mV)');
title('Voltage  Change for second Hodgkin-Huxley Model');
figure();
plot(t,FE);
legend('Forward Euler');
xlabel('Time (ms)');
ylabel('Voltage (mV)');
title('Voltage Change for first Hodgkin-Huxley Model');

figure();
plot(t,Isyn);
xlabel('Time (ms)');
ylabel('Current (microA)');
title('Synaptic Current Change for second Hodgkin-Huxley Model');
%hold on
%plot(FEm,FEn);
