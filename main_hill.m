clear all;
close all;
clc;

global sigma rho beta r0 zeta;
r0=1;
zeta=1;
sigma=10;
rho=28;
beta=8/3;

options = odeset('RelTol', 1.0e-10, 'AbsTol', 1.0e-10);
dt=0.01;
time=0:dt:(100);
x0=[0;randn(3,1)];
[~,out_sin_lorenz] = ode45(@(t,y)sin_lorenz(t,y),time,x0,options);
[~,out_y_lorenz] = ode45(@(t,y)y_lorenz(t,y),time,x0,options);
% [~,out_sin_noise] = ode45(@(t,y)sin_noise(t,y),time,x0,options);
%[~,out_y_noise] = ode45(@(t,y)y_noise(t,y),time,x0,options);
figure(1)
plot(time,out_sin_lorenz(:,1)); hold on;
plot(time,out_sin_lorenz(:,2)); hold on;

figure(2)
plot(time,out_y_lorenz(:,1)); hold on;
plot(time,out_y_lorenz(:,2)); hold on;
% plot(time,out_sin_noise(:,1)); hold on;
%plot(time,out_y_noise(:,1)); hold on;

function dy=sin_lorenz(t,y)
global sigma rho beta r0 zeta;
theta=y(1,1);
X=y(2,1);
Y=y(3,1);
Z=y(4,1);
dy=[-sin(theta)+r0+zeta*X;
    sigma*(Y-X);
    X*(rho-Z)-Y;
    X*Y-beta*Z];
end

function dy=y_lorenz(t,y)
global sigma rho beta r0 zeta;
theta=y(1,1);
X=y(2,1);
Y=y(3,1);
Z=y(4,1);
dy=[-(theta)+r0+zeta*X;
    sigma*(Y-X);
    X*(rho-Z)-Y;
    X*Y-beta*Z];
end


function dy=sin_noise(t,y)
global sigma rho beta r0 zeta;
theta=y(1,1);
X=y(2,1);
Y=y(3,1);
Z=y(4,1);
dy=[-sin(theta)+r0+zeta*randn(1,1);
    sigma*(Y-X);
    X*(rho-Z)-Y;
    X*Y-beta*Z];
end

function dy=y_noise(t,y)
global sigma rho beta r0 zeta;
theta=y(1,1);
X=y(2,1);
Y=y(3,1);
Z=y(4,1);
dy=[-(theta)+r0+zeta*randn(1,1);
    sigma*(Y-X);
    X*(rho-Z)-Y;
    X*Y-beta*Z];
end

