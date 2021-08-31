clear all;
close all;
clc;

%%12073090: 32*32 size, elevator mode, without shear...
%%12073116: 96*32 size, elevator mode, without shear...
%%12075689: 96*32 size, elevator mode, with shear...
%%12075690: 32*32 size, elevator mode, with shear

%%These are local folder
% folder_name='C:\Data\dedalus\IFSC_2D_without_shear\';
%folder_name='C:\Data\dedalus\IFSC_2D_with_shear\';
%folder_name='C:\Data\dedalus\dedalus_12073090\IFSC_2D_without_shear\';
%folder_name='C:\Data\dedalus\dedalus_12075689\IFSC_2D_with_shear\';
folder_name='C:\Data\dedalus\dedalus_12080412\IFSC_2D_without_shear\';

% folder_name='C:\Data\dedalus\dedalus_12075689\IFSC_2D_with_shear\';

% folder_name='/rc_scratch/chli3324/dedalus_12073090/';
% file_name='IFSC_2D_without_shear_s1_random';
% file_name='IFSC_2D_without_shear_s1_small_domain';

% file_name='IFSC_2D_without_shear_s1_elevator_short';
% file_name='IFSC_2D_without_shear_s1_elevator_long';
% file_name='IFSC_2D_with_shear_s1';
file_name='IFSC_2D_with_shear_s1';
h5_name=[folder_name,file_name,'.h5'];


% u=h5read(h5_name,'/tasks/u');
% w=h5read(h5_name,'/tasks/w');
% T=h5read(h5_name,'/tasks/T');
S=h5read(h5_name,'/tasks/S');
S_coeff=h5read(h5_name,'/tasks/S_coeff');
S_coeff=S_coeff.r+1i*S_coeff.i;


Ra_ratio=1.1;
lambda_opt=sqrt(1/2*(-2-Ra_ratio+sqrt(Ra_ratio^2+8*Ra_ratio)))*(3*Ra_ratio-sqrt(Ra_ratio^2+8*Ra_ratio))/(sqrt(Ra_ratio^2+8*Ra_ratio)-Ra_ratio);
k_opt=(1/2*(-2-Ra_ratio+sqrt(Ra_ratio^2+8*Ra_ratio)))^(1/4);
ks=k_opt/32;
tau=0.01;
Pr=7;
R_rho=90.9;
Ri=1;
uL=1/tau/ks*sqrt(Pr*(1-1/R_rho)/Ri);


data{1}.x=kx_list/k_opt;
data{1}.y=kz_list(1:Nz/2)/k_opt;
data{1}.z=log10(mean(abs(S_coeff(1:Nz/2,:,:)),3));
plot_config.zlim_list=[1,-3,0];
plot_config.xlim_list=[1,0,2];
plot_config.ylim_list=[1,-2,2];
plot_config.ztick_list=[1,-3,-2,-1,0];
plot_config.print_size=[1,1100,900];
plot_config.label_list={1,'$k/k_{opt}$','$m/k_{opt}$'};
plot_config.colormap='white_zero';
plot_contour(data,plot_config);


Nx=length(x);
Nz=length(z);
for t_ind=1:length(t)
%     u_int(t_ind)=sum(sum(u(:,:,t_ind).^2));
    E.S(t_ind)=sum(sum(S(:,:,t_ind).^2))/Nx/Nz/2;
    %E.T(t_ind)=sum(sum(T(:,:,t_ind).^2))/Nx/Nz/2;
    %E.u(t_ind)=sum(sum(u(:,:,t_ind).^2))/Nx/Nz/2;
    %E.w(t_ind)=sum(sum(w(:,:,t_ind).^2))/Nx/Nz/2;
    %E.kinetic(t_ind)=E.u(t_ind)+E.w(t_ind);
    average.S(t_ind)=mean(mean(S(:,:,t_ind)));
    %average.T(t_ind)=mean(mean(T(:,:,t_ind)));
    %average.u(t_ind)=mean(mean(u(:,:,t_ind)));
    %average.w(t_ind)=mean(mean(w(:,:,t_ind)));

end

[val,max_ind]=max(E.S);
t_grow=t(1:max_ind);

data{1}.x=t;
data{1}.y=E.S;
% data{2}.x=t_grow;
% data{2}.y=E.S(max_ind)*exp(2*lambda_opt*(t_grow-max(t_grow)));
plot_config.label_list={1,'$t$','$E_S$'};
plot_config.legend_list={0};
% plot_config.legend_list={1,'Simulation','Linear stability'};
plot_config.name=[folder_name,file_name,'E_S.png'];
plot_config.Markerindex=3;
plot_config.user_color_style_marker_list={'k-','bo--'};
plot_line(data,plot_config);

% error('1')

% plot(t_grow,E.S(1)*exp(2*lambda_opt*t_grow)); hold on;
% plot(t,E.S);
% [~,t_ind_1200]=min(abs(t-1200));
clear data
