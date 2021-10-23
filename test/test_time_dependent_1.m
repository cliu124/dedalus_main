clear all;
close all;
clc;

syms t
assume(t,'real')
A=[-1-9*cos(6*t)^2+6*sin(12*t), 12*cos(6*t)^2+9/2*sin(12*t);
    -12*sin(6*t)^2+9/2*sin(12*t), -1-9*sin(6*t)^2-6*sin(12*t)];
x=exp(2*t)*(cos(6*t)+2*sin(6*t));
y=exp(2*t)*(2*cos(6*t)-sin(6*t));

v=[x;y];
res=simplify(diff(v,t)-A*v);
lambda=simplify(eig(A));
lambda_A_T=simplify(eig(A+A'));

e=simplify(x^2+y^2);

