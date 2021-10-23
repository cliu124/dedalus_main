clear all;
close all;
clc;


%construct the A matrix for standard shear flow. Note that
%Re is associated with the shear term

mean_elevator_W_ind=1;
zi=sqrt(-1);
zero_bc=0;
I_bc=1;


%These are specific for the unbounde Couette flow 


R_rho_T2S=0.5;
obj.Pe_T=1;
obj.Pe_S=1;
obj.Pr=10; obj.Re=1/obj.Pr;
obj.Ra_T=1;
obj.Ra_S2T=obj.Ra_T/R_rho_T2S;
obj.dy_T_mean=-1;
obj.dy_S_mean=-1;
obj.tau=0.01;
Ri=2;
A=sqrt(2*obj.Pr*(R_rho_T2S^(-1)-1)/Ri);
f=0.5;

syms t;
A_U=A*cos(f*t);
U_bar=0; %This because we have combine that into the time dependence of D2_bc and D4_bc
d_U_bar=A_U;
dd_U_bar=0;

obj.kx=0.2;
obj.kz=0;
ky0=0;

K2=obj.kx^2+obj.kz^2;
ky_sym=ky0-int(A_U,t);
D2_bc_sym=-(ky_sym)^2;
D4_bc_sym=D2_bc_sym^2;


%These are old code...
dt=0.1;
t_list=0:dt:1000;
d_U_bar=0;
x0=rand(3,1)+1i*rand(3,1);
x(:,1)=x0;
lambda_A_AT_int(1,1)=0;
for t_ind=1:length(t_list)-1
    ky=double(subs(ky_sym,t,t_list(t_ind)));
    ky_time(t_ind)=ky;
    D4_bc=double(subs(D4_bc_sym,t,t_list(t_ind)));
    D2_bc=double(subs(D2_bc_sym,t,t_list(t_ind)));
    
    A11=(D4_bc-2*K2*D2_bc+K2^2*I_bc)+obj.Re*(diag(dd_U_bar)*zi*obj.kx*I_bc-zi*obj.kx*diag(U_bar)*(D2_bc-K2*I_bc)); %%Orr-Sommerfeld operator
    A21= -zi*obj.kz*diag(d_U_bar)*I_bc*obj.Re; %Coulping operator
    A22= -zi*obj.kx*diag(U_bar)*I_bc*obj.Re+(D2_bc-K2*I_bc); %Squire operator  
    %inv_lap=inv([D2_bc-K2*I_bc, zero_bc; zero_bc, I_bc]);
    inv_lap=inv(D2_bc-K2*I_bc);
    A_shear= [inv_lap*A11, zero_bc; A21, A22];

    %%standard shear flow model couple with T and S
    obj.A(:,:,t_ind)=[A_shear, obj.Ra_T*[inv_lap*(-K2*I_bc); zero_bc], -obj.Ra_S2T*[inv_lap*(-K2*I_bc); zero_bc];
      -obj.dy_T_mean*I_bc,zero_bc,-zi*obj.kx*diag(U_bar)*I_bc*obj.Pe_T+(D2_bc-K2*I_bc),zero_bc;
      -obj.dy_S_mean*I_bc,zero_bc,zero_bc,-zi*obj.kx*diag(U_bar)*I_bc*obj.Pe_S+obj.tau*(D2_bc-K2*I_bc)];


    obj.M=[obj.Re*I_bc, zero_bc,zero_bc,zero_bc;
                           zero_bc, obj.Re*I_bc, zero_bc, zero_bc;
                           zero_bc, zero_bc, obj.Pe_T*I_bc, zero_bc;
                           zero_bc, zero_bc, zero_bc, obj.Pe_S*I_bc];
    ind=[1,3,4];
    A_2D=(inv(obj.M(ind,ind))*obj.A(ind,ind,t_ind));
    lambda_A(:,t_ind)=eig(A_2D);
    lambda_A_AT(:,t_ind)=eig(A_2D+transpose(A_2D));
    lambda_A_AT_max_int(:,t_ind)=sum(max(lambda_A_AT))*dt;
    x(:,t_ind+1)=expm(A_2D*dt)*x(:,t_ind);
    u(:,t_ind+1)=-x(1,t_ind+1)*ky/obj.kx;
    e(:,t_ind+1)=(abs(x(1,t_ind+1))^2+abs(u(:,t_ind+1))^2)/4+abs(x(3,t_ind+1)-x(2,t_ind+1))^2*obj.Pr/4/(R_rho_T2S^(-1)-1);
    
end
plot(e);
% plot(x(1,:)); hold on;
% plot(x(2,:)); hold on;
% plot(x(3,:));

% lambda_A_AT_num=subs(lambda_A_AT,t_list);
