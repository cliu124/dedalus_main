clear all;
close all;
clc;

syms mu L(tau) c0(tau) t;

A=[0,1;
    -mu*L(tau)^2, -c0(tau)*L(tau)^2];

sol=simplify(expm(int(A,tau,0,t)),'steps',100);


