clear all;
close all;
clc;

%------setup the main parameter for running
Pr=10; Re=1/Pr; Pe=1;
tau=0.01;gamma=0.5;
R_rho_T2S=gamma/tau;
Ra_T=R_rho_T2S*tau/(1-R_rho_T2S*tau);
Ra_S2T=Ra_T/R_rho_T2S;
mean_elevator_amp_list={'W',4};
mean_kolmogorov=[1,1];
%-------for validation of Holyer (1984)
% 
% %------------for validation of Radko & Smith (2012
% Pr=10; tau=0.01; Pe=1; Re=1/Pr;
% Ra_T=1;  
% % C=4;
% R_rho_T2S=50;
% R_rho_T2S_list=1.1:0.025:2.5;
% %[1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.4,2.5]; %Ra_S2T=Ra_T/R_rho_T2S;
% % Lz=2*pi/1.01;
% C_list=[1.5,2,3,4];
% %%-------------
dy_T_mean=1;
dy_S_mean=1;
mean='elevator';
% mean_elevator_amp_list={'T',1};%logspace(-4,1,4);
mean_elevator_kx=1;
% elevator_lambda_balance_bisection=[];
Ny_full=62; %This needs to be 92 for high Pe, other case 62 or 32 is enough
%Up to 0.3 and 0.4 for the Pe=100, Ri=10
%Up to 0.5 and 0.8 for the Pe=100, Ri=1 
%Up to 3.5 and 1.5 for the Pe=10^4, R1
% kx_list=linspace(0,0.5,30); %0.3
kx_list=linspace(0.01,0.5,20);
% kz_list=linspace(0,1.5,60); %0.4
kz_list=0;
solve='LST'; %%or finished if we would like to skip but just load the data..
debug='';
% Ri=1/((Pe/100)^2/(Pr/10));
%-----------------


%construct the obj that is using the shear flow formulation as Radko 2016..
primitive_Radko2013=DDC_LST('unified');
primitive_Radko2013.tau=tau;
primitive_Radko2013.Pr=Pr;
primitive_Radko2013.Pe_T=Pe;
primitive_Radko2013.Pe_S=Pe;
primitive_Radko2013.Re=Re;
primitive_Radko2013.Ra_T=Ra_T;
primitive_Radko2013.Ra_S2T=Ra_S2T;
primitive_Radko2013.R_rho_T2S=R_rho_T2S;
primitive_Radko2013.mean=mean;
primitive_Radko2013.dy_T_mean=dy_T_mean;
primitive_Radko2013.dy_S_mean=dy_S_mean;
primitive_Radko2013.kx_list=kx_list;
primitive_Radko2013.kz_list=kz_list;
primitive_Radko2013.Ny_full=Ny_full;
primitive_Radko2013.solve=solve;
primitive_Radko2013.debug=debug;
primitive_Radko2013.mean_kolmogorov(1,1:2)=mean_kolmogorov;
% primitive_Radko2013.elevator_lambda_balance_C=elevator_lambda_balance_C;
primitive_Radko2013=primitive_Radko2013.convert_shear();
% primitive_Radko2013.grid_diff{2}(2)=Lz;
primitive_Radko2013.mean_elevator_kx=mean_elevator_kx;
%copy four of these.. but just modify their reduced formulation.
primitive_Radko2013.operator='uvwpTS'; %uvwpTS
primitive_Radko2013=primitive_Radko2013.mean_profile_elevator_kx();
% primitive_Radko2013.Ra_T=0;
% primitive_Radko2013.Ra_S2T=0;
% primitive_Radko2013.Pe_T=0;
% primitive_Radko2013.Pe_S=0;
% primitive_Radko2013.Re=0;
for DDC_LST_ind=1:length(mean_elevator_amp_list{2})
    DDC_LST_list{DDC_LST_ind}=primitive_Radko2013;
    %the parameter setting that are shared by all cases...
    DDC_LST_list{DDC_LST_ind}.mean_elevator_amp_list=...
        {mean_elevator_amp_list{1},mean_elevator_amp_list{2}(DDC_LST_ind)};
    %%compute the linear stability analysis over kx kz
    DDC_LST_list{DDC_LST_ind}=DDC_LST_list{DDC_LST_ind}.solve_kxkz();
%     DDC_LST_list{DDC_LST_ind}=DDC_LST_list{DDC_LST_ind}.convert_IFSC_unit_tuS();
%     DDC_LST_list{DDC_LST_ind}=DDC_LST_list{DDC_LST_ind}.post_eig_kxkz_contour();
    DDC_LST_list{DDC_LST_ind}.post_eig_kx();
%     DDC_LST_list{DDC_LST_ind}=DDC_LST_list{DDC_LST_ind}.elevator_flux;
    %DDC_LST_list{DDC_LST_ind}=DDC_LST_list{DDC_LST_ind}.post_eigenvector();

    data{DDC_LST_ind}.x=kx_list;
    data{DDC_LST_ind}.y=real(cell2mat(DDC_LST_list{DDC_LST_ind}.eig_val_max_list));
    data_elevator_W{1}.x(DDC_LST_ind)=mean_elevator_amp_list{2}(DDC_LST_ind);
    data_elevator_W{1}.y(DDC_LST_ind)=max(data{DDC_LST_ind}.y);
end
% plot_config.ylim_list=[1,2*min(data{5}.y),2*max(data{6}.y)];
plot_config.label_list={1,'$k_z$','$\lambda$'};
% plot_config.legend_list={1,'Primitive','IFSC','MRBC','Stokes'};
%    'Primitive (IFSC unit), $\nabla^2 S=0$',...
%    'IFSC, $\nabla^2 S=0$',...
%     'MRBC, $\nabla^2 S=0$','Stokes, $\nabla^2 S=0$'};
plot_config.Markerindex=3;
plot_config.ylim_list=[0,2*min(data{1}.y)-0.1,2*max(data{1}.y)+0.1];
plot_config.xlim_list=[1,min(kx_list),max(kx_list)];
plot_config.user_color_style_marker_list={'k-','bsquare','m*','r--'...
                                'k-','bo','mx','r--','k--','bd'};
plot_config.title_list={1,['$Pr$=',num2str(Pr),', $\tau$=',...
    num2str(tau),', $R_\rho$=',num2str(R_rho_T2S),...
    ', $\bar{\mathcal{T}}_z=$',num2str(dy_T_mean),...
    ', $\bar{\mathcal{S}}_z=$',num2str(dy_S_mean)]};
plot_config.name=['C:/Figure/DDC_LST/DDC_LST_Pr=',num2str(Pr),'_tau=',...
    num2str(tau),'_R_rho_T2S=',num2str(R_rho_T2S),...
    '_dy_T_mean=',num2str(dy_T_mean),...
    '_dy_S_mean=',num2str(dy_S_mean),...
    '_mean=',mean,'_mean_elevator_kx=',num2str(mean_elevator_kx),'_all.png'];
plot_config.fontsize_legend=24;
plot_config.fontsize=36;
plot_line(data,plot_config);

error('1');
data_elevator_W{2}.x=mean_elevator_amp_list{2};
data_elevator_W{2}.y=DDC_LST_list{1}.elevator_lambda_max*ones(size(data_elevator_W{2}.x));
plot_config.user_color_style_marker_list={'k-','b--','m*','r--'};
plot_config.label_list={1,['$\underline{',mean_elevator_amp_list{1},'}$'],'$\lambda_{M}$'};
plot_config.name=['C:/Figure/DDC_LST/DDC_LST_Pr=',num2str(Pr),'_tau=',...
    num2str(tau),'_R_rho_T2S=',num2str(R_rho_T2S),...
    '_dy_T_mean=',num2str(dy_T_mean),...
    '_dy_S_mean=',num2str(dy_S_mean),...
    '_mean=',mean,'_mean_elevator_kx=',num2str(mean_elevator_kx),'_all_elevator_W.png'];
plot_config.loglog=[0,0];
plot_config.xlim_list=[1,min(data_elevator_W{2}.x),max(data_elevator_W{2}.x)];
plot_line(data_elevator_W,plot_config);




% primitive_Radko2013=DDC_LST('primitive_Radko2013');
% primitive_Radko2013.tau=tau;
% primitive_Radko2013.Pr=Pr;
% primitive_Radko2013.R_rho_T2S=R_rho_T2S;
% primitive_Radko2013=primitive_Radko2013.convert_shear();
% primitive_Radko2013_0tau=primitive_Radko2013;%make one copy but set tau=0
% primitive_Radko2013_0tau.tau=0;
% 
