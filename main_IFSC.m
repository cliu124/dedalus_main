clear all;
close all;
clc;

%%These are local folder
% folder_name='C:\Data\dedalus\IFSC_2D_without_shear\';
%folder_name='C:\Data\dedalus\IFSC_2D_with_shear\';
folder_name='/rc_scratch/chli3324/dedalus_12072802/';
% file_name='IFSC_2D_without_shear_s1_random';
% file_name='IFSC_2D_without_shear_s1_small_domain';

% file_name='IFSC_2D_without_shear_s1_elevator_short';
% file_name='IFSC_2D_without_shear_s1_elevator_long';
% file_name='IFSC_2D_with_shear_s1';
file_name='IFSC_2D_without_shear_s1';
h5_name=[folder_name,file_name,'.h5'];
h5disp(h5_name);
x=h5read(h5_name,'/scales/x/1.0');
z=h5read(h5_name,'/scales/z/1.0');
t=h5read(h5_name,'/scales/sim_time');

Lx=max(x)-min(x);
Lz=max(z)-min(z);

u=h5read(h5_name,'/tasks/u');
w=h5read(h5_name,'/tasks/w');
T=h5read(h5_name,'/tasks/T');
S=h5read(h5_name,'/tasks/S');

Nx=length(x);
Nz=length(z);
for t_ind=1:length(t)
%     u_int(t_ind)=sum(sum(u(:,:,t_ind).^2));
    E.S(t_ind)=sum(sum(S(:,:,t_ind).^2))/Nx/Nz/2;
    E.T(t_ind)=sum(sum(T(:,:,t_ind).^2))/Nx/Nz/2;
    E.u(t_ind)=sum(sum(u(:,:,t_ind).^2))/Nx/Nz/2;
    E.w(t_ind)=sum(sum(w(:,:,t_ind).^2))/Nx/Nz/2;
    E.kinetic(t_ind)=E.u(t_ind)+E.w(t_ind);
    average.S(t_ind)=mean(mean(S(:,:,t_ind)));
    average.T(t_ind)=mean(mean(T(:,:,t_ind)));
    average.u(t_ind)=mean(mean(u(:,:,t_ind)));
    average.w(t_ind)=mean(mean(w(:,:,t_ind)));

end

Ra_ratio=1.1;
lambda_opt=sqrt(1/2*(-2-Ra_ratio+sqrt(Ra_ratio^2+8*Ra_ratio)))*(3*Ra_ratio-sqrt(Ra_ratio^2+8*Ra_ratio))/(sqrt(Ra_ratio^2+8*Ra_ratio)-Ra_ratio);
k_opt=(1/2*(-2-Ra_ratio+sqrt(Ra_ratio^2+8*Ra_ratio)))^(1/4);
ks=k_opt/32;
tau=0.01;
Pr=7;
R_rho=90.9;
Ri=1;
uL=1/tau/ks*sqrt(Pr*(1-1/R_rho)/Ri);

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



% plot(t_grow,E.S(1)*exp(2*lambda_opt*t_grow)); hold on;
% plot(t,E.S);
% [~,t_ind_1200]=min(abs(t-1200));
clear data
for t_ind=1:length(t)
    data{1}.x=x/(2*pi/k_opt);
    data{1}.y=z/(2*pi/k_opt);
    data{1}.z=S(:,:,t_ind);
    plot_config.fontsize=28;
    plot_config.zlim_list=[1,-3.5,3.5];
    plot_config.label_list={1,'$x/l_{opt}$','$z/l_{opt}$'};
    plot_config.colormap='bluewhitered';%bluewhitered
    plot_config.print_size=[1,1200,1200];
    plot_config.name=[folder_name,file_name,'S_contour_t',num2str(round(t(t_ind),2)),'.png'];
    plot_contour(data,plot_config);
end