clear all;
close all;
clc;
global Ra_T Ra_S kx ky tau dy_T_mean dy_S_mean Nz;
Ra_T=6000; %Ra_list=[10000,20000,40000];
kx=0.48*Ra_T.^0.4;
ky=kx;
Ra_S=0; tau=0.01;
dy_T_mean=-1;
dy_S_mean=-1;
% xmesh=linspace(0,1,400);
% solinit=bvpinit(xmesh,@guess_porous);
% sol=bvp4c(@bvpfcn_porous,@bcfcn_porous,solinit);

% 
% tspan=0:0.01:10;
% Nz=100;
% [x,D1]=chebdif(Nz,1); 
% x_mesh=(x+1)/2;
% x_bc=x_mesh(2:end-1);
% solinit=bvpinit(x_bc,@guess_porous);
% y0=[];
% for y_ind=1:size(solinit.y,1)
%     y0=[y0;solinit.y(y_ind,:)'];
% end
% Z=zeros(Nz-2,Nz-2); I=eye(Nz-2,Nz-2);
% M_porous=blkdiag(Z,Z,Z,Z,Z,I/tau,Z,I,Z,I);
% options = odeset('Mass',M_porous);
% [t,y] = ode23t(@marching_porous,tspan,y0,options);

x_bc=linspace(0,1,200);
options=bvpset('NMax',10000);
solinit=bvpinit(x_bc,@guess_porous);
sol=bvp4c(@bvpfcn_porous,@bcfcn_porous,solinit,options);

% 
%solinit=bvpinit(x_bc,@guess_benard);
%sol=bvp4c(@bvpfcn_benard,@bcfcn_benard,solinit,options);

function dz_dtheta = bvpfcn_porous(z,theta)
global Ra_T Ra_S kx ky tau dy_T_mean dy_S_mean;
dz_dtheta = zeros(size(theta));
w_hat=theta(1,:); p_hat=theta(2,:); 
T_hat=theta(3,:); d_T_hat=theta(4,:);
S_hat=theta(5,:); d_S_hat=theta(6,:);
T_0=theta(7,:); d_T_0=theta(8,:);
S_0=theta(9,:); d_S_0=theta(10,:);


dz_dtheta = [-(kx^2+ky^2)*p_hat;
    -w_hat+Ra_T*T_hat-Ra_S*S_hat;
    d_T_hat;
    w_hat.*d_T_0+w_hat*dy_T_mean+(kx^2+ky^2)*T_hat; %
    d_S_hat;
    1/tau*(w_hat.*d_S_0+w_hat*dy_S_mean)+(kx^2+ky^2)*S_hat;
    d_T_0;
    -2*(kx^2+ky^2)*p_hat.*T_hat+2*w_hat.*d_T_hat;
    d_S_0;
    1/tau*(-2*(kx^2+ky^2)*p_hat.*S_hat+2*w_hat.*d_S_hat)
    ];
end

function res = bcfcn_porous(ya,yb)
res = [ya(1)
       yb(1)
       ya(3)
       yb(3)
       ya(5)
       yb(5)
       ya(7)
       yb(7)
       ya(9)
       yb(9)];
end


function g = guess_porous(x_num)
global Ra_T Ra_S kx ky tau dy_T_mean dy_S_mean;
W0=Ra_T;
g=[W0*sin(pi*x_num);
    W0*pi*cos(pi*x_num)/(-(kx^2+ky^2));
    W0*sin(pi*x_num);
    W0*pi*cos(pi*x_num);
    W0*sin(pi*x_num);
    W0*pi*cos(pi*x_num);
    0;
    0;
    0;
    0];
end



function dz_dtheta = bvpfcn_benard(z,theta)
global Ra_T Ra_S kx ky tau dy_T_mean dy_S_mean;
dz_dtheta = zeros(size(theta));
u_tilde=theta(1,:); d_u_tilde=theta(2,:);
v_tilde=theta(3,:); d_v_tilde=theta(4,:);
w_hat=theta(5,:); p_hat=theta(6,:); 
T_hat=theta(7,:); d_T_hat=theta(8,:);
S_hat=theta(9,:); d_S_hat=theta(10,:);
T_0=theta(11,:); d_T_0=theta(12,:);
S_0=theta(13,:); d_S_0=theta(14,:);

dz_dtheta = [d_u_tilde;
    kx*p_hat+(kx^2+ky^2)*u_tilde;
        d_v_tilde;
    ky*p_hat+(kx^2+ky^2)*v_tilde;
    kx*u_tilde+ky*v_tilde;
    kx*d_u_tilde+ky*d_v_tilde-(kx^2+ky^2)*w_hat+Ra_T*T_hat-Ra_S*S_hat;
    d_T_hat;
    w_hat.*d_T_0+w_hat*dy_T_mean+(kx^2+ky^2)*T_hat; %
    d_S_hat;
    1/tau*(w_hat.*d_S_0+w_hat*dy_S_mean)+(kx^2+ky^2)*S_hat;
    d_T_0;
    2*kx*u_tilde.*T_hat+2*ky*v_tilde.*T_hat+2*w_hat.*d_T_hat;
    d_S_0;
    1/tau*(2*kx*u_tilde.*S_hat+2*ky*v_tilde.*S_hat+2*w_hat.*d_S_hat)
    ];
end

function res = bcfcn_benard(ya,yb)
res = [ya(2)
       yb(2)
       ya(4)
       yb(4)
       ya(5)
       yb(5)
       ya(7)
       yb(7)
       ya(9)
       yb(9)
       ya(11)
       yb(11)
       ya(13)
       yb(13)];
end


function g = guess_benard(x_num)
global Ra_T Ra_S kx ky tau dy_T_mean dy_S_mean;
W0=Ra_T;
g=[W0*sin(pi*x_num);
    W0*pi*cos(pi*x_num);
    W0*sin(pi*x_num);
    W0*pi*cos(pi*x_num);
    W0*sin(pi*x_num);
    W0*pi*cos(pi*x_num)/(-(kx^2+ky^2));
    W0*sin(pi*x_num);
    W0*pi*cos(pi*x_num);
    W0*sin(pi*x_num);
    W0*pi*cos(pi*x_num);
    0;
    0;
    0;
    0];
end

function dt_theta=marching_porous(t,theta)
global Ra_T Ra_S kx ky tau dy_T_mean dy_S_mean Nz;
[x,D1]=chebdif(Nz,1);
D1=D1(2:end-1,2:end-1,1);
D1=2*D1; %This is due to a stretching..
z=1;
D1_large=blkdiag(D1,D1,D1,D1,D1,D1,D1,D1,D1,D1);
theta_mat=reshape(theta,[Nz-2,10]);
bvpfcn_theta=reshape(bvpfcn_porous(z,theta_mat')',[(Nz-2)*10,1]);
dt_theta=D1_large*theta-bvpfcn_theta;
end
