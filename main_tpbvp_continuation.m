clear all;
close all;

clc;

global C R0 R2 theta_10_fun d_theta_10_fun d_theta_bar_fun kx;
flow='darcy_hewitt';
switch flow
    case 'darcy_hewitt'
    C=2;
    darcy_hewitt='2D';
    switch darcy_hewitt
        case '2D'
            Ra=[10000,20000,40000];
            kx_list=0.48*Ra.^0.4;
%                 R0_list=(1/0.48)^(1/0.4);
            %R2_list=[30,30,30];
            R2_list=[20,25,30,40,50,70];
            R0_list=[1,1,1,1,1,1];
            %R0_list=(Ra-R2_list)./kx_list.^2;
            %R1_list=(Ra-kx_list.^2)./kx_list;
        case '3D'
            Ra=[4000,8000,16000];
            kx_list=0.17*Ra.^0.52;
            R2_list=[0,0,0];%Ra-R0_list*kx_list.^2;
            R0_list=(Ra-R2_list)./kx_list.^2;%(1/0.17)^(1/0.52);
            R1_list=(Ra-kx_list.^2)./kx_list;
    end
end

%%main loop to solve...
for R2_ind=1:length(R2_list)
    R2_final=R2_list(R2_ind);
    R0_final=R0_list(R2_ind);
    R2_con_list=linspace(pi^2*C,R2_final,10);%%The list for performing continuation solving...
    R0_con_list=linspace(0,R0_final,10);
    xmesh=linspace(0,1,100);
    sol=bvpinit(xmesh,@guess_linear);
    for R2_con_ind=1:length(R2_con_list)
        for R0_con_ind=1:length(R0_con_list)
            R2=R2_con_list(R2_con_ind);
            R0=R0_con_list(R0_con_ind);
            options = bvpset('FJacobian',@bvpfcnjac,'BCJacobian',@bvpbcjac);%,'Vectorized','on'
            sol = bvp4c(@bvpfcn, @bcfcn, sol,options);
        end
     end
    bvp_sol{R2_ind}.sol=sol;
end




function dz_dtheta = bvpfcn(z,theta)
global C R2 R0;
dz_dtheta = zeros(size(theta));
dz_dtheta = [theta(2)
            theta(3)*theta(4)
            theta(4)
            1/C*(-R2*theta(3)+R0*theta(3)*theta(2) )];
end

function jac=bvpfcnjac(z,theta)
global C R2 R0;
jac=[0,1,0,0;
    0,0,theta(4), theta(3);
    0,0,0,1;
    0,1/C*(R0*theta(3)),1/C*(-R2+R0*theta(2)),0];
end

function res = bcfcn(ya,yb)
res = [ya(1)
       yb(1)
       ya(3)
       yb(3)];
end

function [dBCdya,dBCdyb]=bvpbcjac(ya,yb)
dBCdya=[1,0,0,0;
        0,0,0,0;
        0,0,1,0;
        0,0,0,0];

dBCdyb=[0,0,0,0;
        1,0,0,0;
        0,0,0,0;
        0,0,1,0];

end

function dz_dtheta = bvpfcn_rotation(z,theta)
global R2 kx;
dz_dtheta = zeros(size(theta));
dz_dtheta = [theta(2)
    -4*theta(3)*theta(4)*theta(2)/(1+2*theta(3)^2)+4*theta(3)*theta(4)/(1+2*theta(3)^2)
    theta(4)
    kx^6*theta(3)-R2*kx^2*(1-theta(2))*theta(3)];

end


function g = guess(x_num)
global R2 theta_10_fun d_theta_10_fun d_theta_bar_fun;;
% for z_ind=1:length(d_theta_bar_num)
%    theta_num(z_ind)=sum(d_theta_bar_num(1:z_ind)); 
% end

g = [-R2/2*sin(2*pi*x_num)
     d_theta_bar_fun(x_num)
     theta_10_fun(x_num)
     d_theta_10_fun(x_num)
     ];
end

function g = guess_zero(x_num)
% for z_ind=1:length(d_theta_bar_num)
%    theta_num(z_ind)=sum(d_theta_bar_num(1:z_ind)); 
% end

g = [0
     0
     0
     0
     ];
end


function g = guess_linear(x_num)
% for z_ind=1:length(d_theta_bar_num)
%    theta_num(z_ind)=sum(d_theta_bar_num(1:z_ind)); 
% end

g = 100*[-1/(8*pi)*sin(2*pi*x_num)
     -1/4*cos(2*pi*x_num)
     sin(pi*x_num)
     pi*cos(pi*x_num)
     ];
end
