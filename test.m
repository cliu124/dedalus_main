clear all;
close all;
clc;

obj.Ra_ratio=200
obj.k_opt=(1/2*(-2-obj.Ra_ratio+sqrt(obj.Ra_ratio^2+8*obj.Ra_ratio)))^(1/4);
obj.lambda=sqrt(1/2*(-2-obj.Ra_ratio+sqrt(obj.Ra_ratio^2+8*obj.Ra_ratio)))*(3*obj.Ra_ratio-sqrt(obj.Ra_ratio^2+8*obj.Ra_ratio))/(sqrt(obj.Ra_ratio^2+8*obj.Ra_ratio)-obj.Ra_ratio);
u=2*pi*obj.lambda/obj.k_opt;

syms kx ky kz0 Ra S t T_bar_z;
A=T_bar_z*Ra*(kx^2+ky^2+(kz0-S*kx*t)^2)*(kx^2+ky^2)/((kx^2+ky^2+(kz0-S*kx*t)^2)^3+kx^2+ky^2)-(kx^2+ky^2+(kz0-S*kx*t)^2);

A_int=simplify(int(A,t))


