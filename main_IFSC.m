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
% folder_name='C:\Data\dedalus\dedalus_12075689\IFSC_2D_with_shear\';
% folder_name='C:\Data\dedalus\dedalus_12080489\IFSC_2D_without_shear\';

% folder_name='C:\Data\dedalus\dedalus_12075689\IFSC_2D_with_shear\';

% folder_name='/rc_scratch/chli3324/dedalus_12073090/';
% file_name='IFSC_2D_without_shear_s1_random';
% file_name='IFSC_2D_without_shear_s1_small_domain';

% file_name='IFSC_2D_without_shear_s1_elevator_short';
% file_name='IFSC_2D_without_shear_s1_elevator_long';
% file_name='IFSC_2D_with_shear_s1';

% folder_name='C:\Data\dedalus\dedalus_12073090\IFSC_2D_without_shear\';
% file_name='IFSC_2D_without_shear_s1';
% folder_name='C:\Data\dedalus\dedalus_12075689\IFSC_2D_with_shear\';
% file_name='IFSC_2D_with_shear_s1';

slurm_num={'12073090',... %%without shear, 32*32, A_elevator=1, A_noise=0
    '12073116',... %%without shear, 32*96, A_elevator=1, A_noise=0
    '12075689',... %%with shear, 32*96, A_elevator=1, A_noise=0.01
    '12083149',...%%without shear, 32*32, A_elevator=1, A_noise=0.01, A_shear=1, show the effect of shear
    '12083150',...%%without shear, 32*32, A_elevator=0, A_noise=0.01, A_shear=1, show the effect of only shear... then decay
    '12083221',...%%without shear, 32*32, A_elevator=1, A_noise=0.01, show effect of random noise...
    '12083491',...with shear, 96*32, A_elevator=1, A_noise=0.01, ks=0.05
    '12083494',...with shear, 96*32, A_elevator=1, A_noise=0.01, ks=0.0264
    '12084941',...%%without shear, 192*32, A_elevator=1,
    '12085402',... %%% with shear, 32*32 box size, Ra=1.1, very large initial condition in elevator and shear
    '12085400',...%%with shear, 32*96, A_elevator=1, A_noise=0.01, ks: 0.008806549460544114
    '12085401',...%%with shear, 32*96, A_elevator=1, A_noise=0.01, ks: 0.017613098921088227
    '12085887',...%%without shear, 32*192, A_elevator=1, A_noise=0, A_shear=0
    '12086951',...%%without shear, 96*32 domain, A_elevator=1, A_noise=0, A_shear=9
    '12088536',...%% with shear, 32*32 domain, A_elevator=1, A_noise=0, A_shear=0
    '12088537',...%%with shear, 32*32 domain, A_elevator=1, A_noise=0, A_shear=1
    '12088673',...%%without shear, 8*192 domain, A_elevator=1, A_noise=0, A_shear=0
    '12089740',...%%with shear, 32*32, A_elevator=1, A_noise=0, A_shear=1, long time up to 20000
    '12089741',...%%with shear, 32*96, A_elevator=1, A_noise=0.01
    '12090288',...%%with shear 8*8, A_elevator=1, A_noise=0.01, A_shear=0
    '12099609',...%%with shear 8*24, ks=2*2pi/Lz, A_elevator=1, A_noise=0.01, u_L=1, F_sin: 0.001240885014416157
    '12099641',...%%with shear 8*24, ks=2*2pi/Lz, A_elevator=1, A_noise=0.01, u_L=0.1, F_sin: 0.0001240885014416157
    '12099647',...%%with shear 8*24, ks=2*2pi/Lz, A_elevator=1, A_noise=0.01, u_L=0.01, F_sin: 0.00001240885014416157
    '12099955',...%%with shear 8*24, ks=2pi/Lz, A_elevator=1, A_noise=0.01, u_L=0.01, F_sin: 3.1022125360403927e-06, 
    '12099956',...%%with shear 8*24, ks=2pi/Lz, A_elevator=1, A_noise=0.01, u_L=0.1, F_sin: 3.1022125360403927e-05, 
    '12099957',...%%with shear 8*24, ks=2pi/Lz, A_elevator=1, A_noise=0.01, u_L=1, F_sin: 3.1022125360403927e-04, 
    '12100768',...%%with shear 8*24, ks=2pi/Lz, A_elevator=1, A_noise=0.01, u_L=1, F_sin: 3.1022125360403927e-04, F_sin_2ks: 0.001240885014416157, 
    '12100769',...%%with shear 8*24, ks=2pi/Lz, A_elevator=1, A_noise=0.01, u_L=1, F_sin: 0.00031022125360403925, F_sin_2ks: 0.001240885014416157, F_sin_3ks: 0.0027919912824363536, 
    '12100793',...%%with shear 8*24, ks=2pi/Lz, A_elevator=1, A_noise=0.01, u_L=1, F_sin: 0.00031022125360403925, F_sin_2ks: 0.001240885014416157, F_sin_3ks: 0.0027919912824363536, F_sin_4ks: 0.004963540057664628, 
    '12100794',...%%with shear 8*24, ks=2pi/Lz, A_elevator=1, A_noise=0.01, u_L=1, F_sin: 3.1022125360403927e-04, F_sin_2ks: 0.001240885014416157, phase_2ks: pi/2
    '12100795',...%%with shear 8*24, ks=2pi/Lz, A_elevator=1, A_noise=0.01, u_L=1, F_sin: 0.00031022125360403925, F_sin_2ks: 0.001240885014416157, phase_2ks: pi/4
    '12131253',...%%with shear, 8*24, ks: 0.4227143741061175, A_elevator=1, A_noise=0, u_L=1,
    '12131254',...%%with shear, 8*24, ks: 0.10567859352652938, A_elevator=1, A_noise=0, u_L=1,
    '12131256',...%%with shear, 8*24, ks: 0.07045239568435291, A_elevator=1, A_noise=0, u_L=1,
    '12131257',...%%with shear, 8*24, ks: 0.21135718705305875, A_elevator=1, A_noise=0, u_L=1,
    '12131258',...%%with shear, 8*24, ks: 0.05283929676326469, A_elevator=1, A_noise=0, u_L=1,
    '12131259',...%%with shear, 8*24, ks: 0.035226197842176454, A_elevator=1, A_noise=0, u_L=1,
    '12131260',...%%with shear, 8*24, ks: 0.017613098921088227, A_elevator=1, A_noise=0, u_L=1,
    '12132188',...%%with shear, 8*24, ks: 0.017613098921088227, A_elevator=1, A_noise=0.01, u_L=1,
    '12132197',...%with shear, 8*24, ks: 0.017613098921088227, A_elevator=1, A_noise=0.01, u_L=10^{-5}, F_sin: 3.1022125360403926e-09, 
    '12132206',...%with shear, 8*24, ks: 0.017613098921088227, A_elevator=1, A_noise=0.01, u_L=10^{-4}, F_sin: 3.1022125360403926e-08, 
    '12132211',...%with shear 8*24, ks: 0.017613098921088227,A_elevator=1, A_noise=0.01, u_L=10^{-3}, F_sin: 3.1022125360403926e-07, 
    '12132414',...%%with shear 8*24, ks: 0.017613098921088227,A_elevator=1, A_noise=0.01, u_L=1, Ra_ratio=5
    '12132615',...%%with shear, 8*24, ks: 0.017613098921088227,A_elevator=1, A_noise=0.01, u_L=10, Ra_ratio=1.1
    '12132761',...%%with shear, 8*24,  ks: 0.017613098921088227,A_elevator=1, A_noise=0.01, u_L=10, Ra_ratio=2, time up to 1000
    '12132764',...%%with shear, 8*24, ks: 0.017613098921088227,A_elevator=1, A_noise=0.01, u_L=1, Ra_ratio=5, time up to 500
    '12135159',...%%with shear, 8*24, ks=0.0176, A_elevator=1, A_noise=0.01, u_
    '12135442',...%%without shear, 8*24, Ra_ratio=5
    '12135952',...
    '12135952',...
    '12136034',...
    '12136695',...
    '12144003',...
    '12148188',...
    '12148389',...
    '12148590',...
    '12149388',...
    '12149389',...
    '12149390',...
    '12149391',...
    '12150307',...
    '12150308',...
    '12150309',...
    '12150310',...
    '12170693',...
    '12170695',...
    '12170697',...
    '12173787',...
    '12173788',...
    '12173789',...
    '12183054',...
    '12200160',...
    '12203442',...
    '12210427',...
    'end'};
%     '12089742',...

flag.print=1;
flag.video=1;
flag.visible=1;
for slurm_ind=length(slurm_num)-1:length(slurm_num)-1
    %find(strcmp(slurm_num,'12136034'))
    %length(slurm_num)-1:length(slurm_num)-1
    
    %%change the path into D... just store data in the external disk...
    h5_name=['C:\Data\dedalus\dedalus_',...
        slurm_num{slurm_ind},...
        '\analysis\analysis_s1.h5'];

     set(0,'DefaultFigureVisible','on')
     IFSC_post_my{slurm_ind}=IFSC_post(h5_name,flag);
     IFSC_post_my{slurm_ind}=IFSC_post_my{slurm_ind}.S_x_ave();
     %IFSC_post_my{slurm_ind}=IFSC_post_my{slurm_ind}.w_x_ave();
     IFSC_post_my{slurm_ind}=IFSC_post_my{slurm_ind}.u_fluctuation_x_ave();
     %IFSC_post_my{slurm_ind}=IFSC_post_my{slurm_ind}.spectrum_TKE_average();
     IFSC_post_my{slurm_ind}=IFSC_post_my{slurm_ind}.E_S_time(1);
     IFSC_post_my{slurm_ind}=IFSC_post_my{slurm_ind}.spectrum_S_average();
     IFSC_post_my{slurm_ind}=IFSC_post_my{slurm_ind}.u_total_xt_ave();
     IFSC_post_my{slurm_ind}=IFSC_post_my{slurm_ind}.T_total_xt_ave();
     IFSC_post_my{slurm_ind}=IFSC_post_my{slurm_ind}.S_total_xt_ave();
    
     IFSC_post_my{slurm_ind}.print=0; IFSC_post_my{slurm_ind}.visible=0;
     IFSC_post_my{slurm_ind}=IFSC_post_my{slurm_ind}.snapshot_S();
     IFSC_post_my{slurm_ind}=IFSC_post_my{slurm_ind}.spectrum_S_snapshot();
     IFSC_post_my{slurm_ind}=IFSC_post_my{slurm_ind}.u_laminar();
    %IFSC_post_my{slurm_ind}=IFSC_post_my{slurm_ind}.E_TKE_time();
    %IFSC_post_my{slurm_ind}=IFSC_post_my{slurm_ind}.spectrum_TKE_average();
    %IFSC_post_my{slurm_ind}=IFSC_post_my{slurm_ind}.uS_x_ave();
    %IFSC_post_my{slurm_ind}=IFSC_post_my{slurm_ind}.wS_x_ave();
    
end



error('1')
% 
% vslow=VideoWriter('TravelingWaveTwo20181111.avi');
% vslow.FrameRate=10;
% open(vslow);
% writeVideo(vslow,Fxnorm);
% close(vslow);

    % file_name='';
    % h5_name=[folder_name,file_name,'.h5'];


    % u=h5read(h5_name,'/tasks/u');
    % w=h5read(h5_name,'/tasks/w');
    % T=h5read(h5_name,'/tasks/T');
%     S=h5read(h5_name,'/tasks/S');
%     S_coeff=h5read(h5_name,'/tasks/S_coeff');
%     S_coeff=S_coeff.r+1i*S_coeff.i;


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
