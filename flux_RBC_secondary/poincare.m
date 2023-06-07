clear all;
close all;
clc;

z=linspace(0,1000,1000);
K=1;
k=10; 
phi=0;
delta=1.1;
Phi=K*z.^delta.*cos(k*log(z)+phi);
plot(z,Phi);hold on;
plot(z,z);
