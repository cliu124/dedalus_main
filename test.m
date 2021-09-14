clear all;
close all;
clc;

obj.Ra_ratio=200
obj.k_opt=(1/2*(-2-obj.Ra_ratio+sqrt(obj.Ra_ratio^2+8*obj.Ra_ratio)))^(1/4);
obj.lambda=sqrt(1/2*(-2-obj.Ra_ratio+sqrt(obj.Ra_ratio^2+8*obj.Ra_ratio)))*(3*obj.Ra_ratio-sqrt(obj.Ra_ratio^2+8*obj.Ra_ratio))/(sqrt(obj.Ra_ratio^2+8*obj.Ra_ratio)-obj.Ra_ratio);
u=2*pi*obj.lambda/obj.k_opt;






