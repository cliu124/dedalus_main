clear all;
close all;
clc;

%------setup the main parameter for running
Pr=1; Re=0; Pe=1;
Ra_T=1;
Ra_S2T=0;
tau=1;
R_rho_T2S=0;
mean_elevator_amp_list={'W',[2^5]};%logspace(1,5,20)
% mean_elevator_amp_list{2}=mean_elevator_amp_list{2}(1:2);
darcy=1;
Lz=2*pi*2;
kappa_T_elevator=1;
%%-------------

dy_T_mean=-1;
dy_S_mean=-1;
mean='elevator';
% mean_elevator_amp_list={'T',1};%logspace(-4,1,4);
mean_elevator_kx=1;
% elevator_lambda_balance_bisection=[];

% Ny_full_list=[62,62];
Ny_full_list=[92*ones(1,5),62*ones(1,5),92*ones(1,5),...
     122*ones(1,5)]; %This needs to be 92 for high Pe, other case 62 or 32 is enough
%Up to 0.3 and 0.4 for the Pe=100, Ri=10
%Up to 0.5 and 0.8 for the Pe=100, Ri=1 
%Up to 3.5 and 1.5 for the Pe=10^4, R1
% kx_list=linspace(0,0.5,30); %0.3
kx_list=logspace(-2,-0.5,200);%linspace(0.01,0.5,10);
% kz_list=linspace(0,1.5,60); %0.4
kz_list=0;
solve='finished'; %%or finished if we would like to skip but just load the data..
debug='hewitt2013';
% Ri=1/((Pe/100)^2/(Pr/10));
%-----------------


%construct the obj that is using the shear flow formulation as Radko 2016..
primitive_Radko2013=DDC_LST('unified');
primitive_Radko2013.Pr=Pr;
primitive_Radko2013.Pe_T=Pe;
primitive_Radko2013.Pe_S=0;
primitive_Radko2013.Re=Re;
primitive_Radko2013.Ra_T=Ra_T;
primitive_Radko2013.Ra_S2T=Ra_S2T;
primitive_Radko2013.mean=mean;
primitive_Radko2013.dy_T_mean=dy_T_mean;
primitive_Radko2013.kx_list=kx_list;
primitive_Radko2013.kz_list=kz_list;
primitive_Radko2013.solve=solve;
primitive_Radko2013.debug=debug;
primitive_Radko2013.darcy=darcy;
primitive_Radko2013.kappa_T_elevator=kappa_T_elevator;
% primitive_Radko2013.elevator_lambda_balance_C=elevator_lambda_balance_C;
primitive_Radko2013=primitive_Radko2013.convert_shear();
primitive_Radko2013.grid_diff{2}(2)=Lz;
primitive_Radko2013.mean_elevator_kx=mean_elevator_kx;
%copy four of these.. but just modify their reduced formulation.
primitive_Radko2013.operator='uvwpTS'; %uvwpTS

for DDC_LST_ind=1:length(mean_elevator_amp_list{2})
    DDC_LST_list{DDC_LST_ind}=primitive_Radko2013;
    %the parameter setting that are shared by all cases...
    DDC_LST_list{DDC_LST_ind}.mean_elevator_amp_list=...
        {mean_elevator_amp_list{1},mean_elevator_amp_list{2}(DDC_LST_ind)};
    DDC_LST_list{DDC_LST_ind}.Ny_full=Ny_full_list(DDC_LST_ind);
    DDC_LST_list{DDC_LST_ind}=DDC_LST_list{DDC_LST_ind}.mean_profile_elevator_kx();
    %%compute the linear stability analysis over kx kz
    DDC_LST_list{DDC_LST_ind}=DDC_LST_list{DDC_LST_ind}.solve_kxkz();
%     DDC_LST_list{DDC_LST_ind}=DDC_LST_list{DDC_LST_ind}.convert_IFSC_unit_tuS();
%     DDC_LST_list{DDC_LST_ind}=DDC_LST_list{DDC_LST_ind}.post_eig_kxkz_contour();
    DDC_LST_list{DDC_LST_ind}.post_eig_kx();
%     DDC_LST_list{DDC_LST_ind}=DDC_LST_list{DDC_LST_ind}.elevator_flux;
    DDC_LST_list{DDC_LST_ind}=DDC_LST_list{DDC_LST_ind}.post_eigenvector();

    data{DDC_LST_ind}.x=kx_list;
    data{DDC_LST_ind}.y=real(cell2mat(DDC_LST_list{DDC_LST_ind}.eig_val_max_list));
    data_elevator_W{1}.x(DDC_LST_ind)=mean_elevator_amp_list{2}(DDC_LST_ind);
    [data_elevator_W{1}.y(DDC_LST_ind),max_ind]=max(data{DDC_LST_ind}.y);
    data_elevator_kx{1}.y(DDC_LST_ind)=kx_list(max_ind);
    data_elevator_complex{1}.y(DDC_LST_ind)=DDC_LST_list{DDC_LST_ind}.eig_val_max_list{max_ind};
    
    data_elevator_W{3}.x(DDC_LST_ind)=mean_elevator_amp_list{2}(DDC_LST_ind);
    data_elevator_W{3}.y(DDC_LST_ind)=max(real(cell2mat(DDC_LST_list{DDC_LST_ind}.eig_val_max_A_AT_all_list)));
end

% error('1');
secondary_porous_media=get_secondary_porous_media(DDC_LST_list{DDC_LST_ind});

data_elevator_W{1}.x=data_elevator_W{1}.x;
data_elevator_W{1}.y=data_elevator_W{1}.y;
data_elevator_W{2}.x=secondary_porous_media.sigma(:,1);
data_elevator_W{2}.y=secondary_porous_media.sigma(:,2);
data_elevator_W{3}.x=data_elevator_W{1}.x;
data_elevator_W{4}.x=data_elevator_W{1}.x;
data_elevator_W{4}.y=data_elevator_W{1}.x;
%data_elevator_W{1}.x.^(-4/9);
% data_elevator_W{2}.x=mean_elevator_amp_list{2};
% data_elevator_W{2}.y=DDC_LST_list{1}.elevator_lambda_max*ones(size(data_elevator_W{2}.x));
plot_config.Markerindex=3;
plot_config.user_color_style_marker_list={'k-','b*','r-.','m--'};
plot_config.label_list={1,'$A$','$\lambda_{M}$'};
plot_config.name=['C:/Figure/DDC_LST/DDC_LST_Pr=',num2str(Pr),'_tau=',...
    num2str(tau),'_R_rho_T2S=',num2str(R_rho_T2S),...
    '_dy_T_mean=',num2str(dy_T_mean),...
    '_dy_S_mean=',num2str(dy_S_mean),...
    '_mean=',mean,'_mean_elevator_kx=',num2str(mean_elevator_kx),'_Hewitt2013.png'];
plot_config.loglog=[1,1];
plot_config.xtick_list=[1,1,10,100,1000,10000,100000];
plot_config.xlim_list=[1,1,10^5]; plot_config.ylim_list=0;
plot_config.ytick_list=[1,1,10,100,1000,10000,100000];
plot_config.legend_list={1,'$max\{Re[eig(A)]\}$','Hewitt et al. (2013)','$max\{Re[eig(A+A^*)]\}$','$\lambda_M=A$'};
% plot_config.xlim_list=[1,min(data_elevator_W{2}.x),max(data_elevator_W{2}.x)];
plot_line(data_elevator_W,plot_config);



data_elevator_kx{1}.x=mean_elevator_amp_list{2};
data_elevator_kx{1}.y=data_elevator_kx{1}.y;%.*data_elevator_kz{1}.x.^(2/9);
data_elevator_kx{2}.x=secondary_porous_media.kx(:,1);
data_elevator_kx{2}.y=secondary_porous_media.kx(:,2);

%data_elevator_W{2}.y=DDC_LST_list{1}.elevator_lambda_max*ones(size(data_elevator_W{2}.x));
plot_config.user_color_style_marker_list={'k-','b*','m*','r--'};
plot_config.label_list={1,['$\underline{',mean_elevator_amp_list{1},'}$'],'$\lambda_{M}$'};
plot_config.name=['C:/Figure/DDC_LST/DDC_LST_Pr=',num2str(Pr),'_tau=',...
    num2str(tau),'_R_rho_T2S=',num2str(R_rho_T2S),...
    '_dy_T_mean=',num2str(dy_T_mean),...
    '_dy_S_mean=',num2str(dy_S_mean),...
    '_mean=',mean,'_mean_elevator_kx=',num2str(mean_elevator_kx),'_Hewitt2013_kx.png'];
plot_config.loglog=[1,1];
plot_config.xlim_list=0; plot_config.ylim_list=0;
% plot_config.xlim_list=[1,min(data_elevator_W{2}.x),max(data_elevator_W{2}.x)];
plot_line(data_elevator_kx,plot_config);



