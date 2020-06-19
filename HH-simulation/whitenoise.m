clear all;
clc;
L=5000;
mu=0;
Sigma=5;

X=Sigma+randn(L,1)+mu;
figure(1);
plot(X)