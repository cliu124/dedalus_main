clear all;
close all;
clc;

% global L;
Amp=1; omega=1;
sigma=-1;
mu=3;
Q=1;
k=1;

dt=0.0001;
t_list=0:dt:3;
aQ=zeros(size(t_list));
aQ(1)=10;%sqrt(mu-Q^2);
eigmax_A_AT=zeros(size(t_list));
int_eigmax_A_AT=zeros(size(t_list));
L_list=zeros(size(t_list));
for t_ind=1:length(t_list)
    t=t_list(t_ind);
    L=Amp*exp(sigma*t);
    L_list(t_ind)=L;
    A=[mu-(Q+k)^2/L^2-2*aQ(t_ind)^2, -aQ(t_ind)^2;
        -aQ(t_ind)^2, mu-(Q-k)^2/L^2-2*aQ(t_ind)^2];
    eigmax_A_AT(t_ind)=max(eig(A+A'));
    
    int_eigmax_A_AT(t_ind)=sum(eigmax_A_AT(1:t_ind));
    d_aQ=(mu-Q^2/L^2)*aQ(t_ind)-aQ(t_ind)^3;
    aQ(t_ind+1)=aQ(t_ind)+dt*d_aQ;
    
end

aQ(end)=[];
int_eigmax_A_AT=int_eigmax_A_AT*dt;


