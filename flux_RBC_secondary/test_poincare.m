clear all;
close all;
clc;

syms z0 d gamma rho omega tau x0
% tau=-1/gamma*log(z0/d);
A=[-rho,-omega;
    omega,-rho];
% output=simplify(expm(-1/gamma*log(z0/d)*A)*[x0;0],'steps',100);
output=simplify(expm(tau*A)*[x0;0],'Criterion','preferReal');
output_sim1=simplify(output(1),'All',true,'Criterion','preferReal');
output_sim2=simplify(output(2),'All',true,'Criterion','preferReal');


