clear all;
close all;
clc;

global m c k;
m=1; c=1; k=1;
tspan=[0,10];
v0=[1;2];
[t,v]=ode45(@f,tspan,v0);
x_num=v(:,1);
y_num=v(:,2);

figure(1)
plot(t,x_num,'b'); hold on;
%plot(t,double(subs(x,t)),'--r');

figure(2)
plot(t,y_num,'b'); hold on;
%plot(t,double(subs(y,t)),'--r');


function dv=f(t,v)
global m c k;
    %A=[0,1;
    %    -k/m, -(t)*c/m];
    A=[-1-9*cos(6*t)^2+6*sin(12*t), 12*cos(6*t)^2+9/2*sin(12*t);
    -12*sin(6*t)^2+9/2*sin(12*t), -1-9*sin(6*t)^2-6*sin(12*t)];
    dv=A*v;
end




