clear all;
close all;
clc;

x=linspace(0,100*pi,100000);
A=1;
B=pi;
y=sin(A*x)+sin(B*x);
plot(x,y);

