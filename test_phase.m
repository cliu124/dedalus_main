clear all;
close all;
clc;


x0=randn(2,1);
N=1000;
x(1:2,1)=x0;
x(3,1)=0;
u_old=8+4*1i;
%x_ref=[1,0];
for i=1:N
    u=x(1,i)+1i*x(2,i);
    eta=x(3,i);
    G=[real(u-u*abs(u)^2+eta*(1i)*u);
        imag(u-u*abs(u)^2+eta*(1i)*u);
        real(1i*u_old)*real(u-u_old)+imag(1i*u_old)*imag(u-u_old)];
    dG=[[1-abs(u)^2-2*real(u)^2, -2*real(u)*imag(u);
        -2*real(u)*imag(u), 1-abs(u)^2-2*imag(u)^2] ,[real((1i)*u);imag((1i)*u)];
        [real(1i*u_old),imag(1i*u_old)], 0];
    x(:,i+1)=x(:,i)-inv(dG)*G;
    x_abs(1,i+1)=sqrt(x(1,i+1)^2+x(2,i+1)^2);
end
plot(x(1,:),x(2,:)); hold on;
plot(x(1,end),x(2,end),'r*')
x_abs(end)



