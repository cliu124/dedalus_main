clear all;
close all;
clc;

Pr=7;
tau=1/24;
R_rho_T2S=1.5;
dy_T_mean=1;
dy_S_mean=1;
% mean_kolmogorov=[1,1];
mean='no';
solve='LST';
operator='uvwpTS';

primitive_IFSC_unit_tuS=DDC_LST('primitive_IFSC_unit_tuS');
primitive_IFSC_unit_tuS.tau=tau;
primitive_IFSC_unit_tuS.Pr=Pr;
primitive_IFSC_unit_tuS.R_rho_T2S=R_rho_T2S;
primitive_IFSC_unit_tuS.post_super_exponential();
% primitive_IFSC_unit_tuS.operator='v_omega_y';
primitive_IFSC_unit_tuS=primitive_IFSC_unit_tuS.convert_shear();
primitive_IFSC_unit_tuS_0tau=primitive_IFSC_unit_tuS; %make one copy but set tau=0
primitive_IFSC_unit_tuS_0tau.tau=0;

IFSC=DDC_LST('IFSC');
IFSC.Ra_ratio=1/tau/R_rho_T2S;
IFSC=IFSC.convert_shear();
IFSC_0tau=IFSC;%make one copy but set tau=0
IFSC_0tau.tau=0;

MRBC=DDC_LST('MRBC');
MRBC.Ra_ratio=1/tau/R_rho_T2S;
MRBC.Sc=Pr/tau;
MRBC=MRBC.convert_shear();
MRBC_0tau=MRBC;%make one copy but set tau=0
MRBC_0tau.tau=0;


Stokes=DDC_LST('Stokes');
Stokes.Ra_ratio=1/tau/R_rho_T2S;
Stokes.tau=tau;
Stokes=Stokes.convert_shear();
Stokes_0tau=Stokes;%make one copy but set tau=0
Stokes_0tau.tau=0;


% DDC_LST_list={IFSC,IFSC_0tau,MRBC,MRBC_0tau,...
%     primitive_IFSC_unit_tuS,primitive_IFSC_unit_tuS_0tau};
DDC_LST_list={primitive_IFSC_unit_tuS,IFSC,MRBC,Stokes,...
        primitive_IFSC_unit_tuS_0tau,IFSC_0tau,MRBC_0tau,Stokes_0tau};

%    primitive_Radko2013,primitive_Radko2013_0tau,...

% kx_list=1;

% kx_list=0.5;
kx_list=linspace(0.01,3,20);%linspace(0,3,20);
for DDC_LST_ind=1:length(DDC_LST_list)
%     DDC_LST_list{DDC_LST_ind}=DDC_LST_list{DDC_LST_ind}.convert_shear();
    %the parameter setting that are shared by all cases...
    DDC_LST_list{DDC_LST_ind}.kx_list=kx_list;
    DDC_LST_list{DDC_LST_ind}.mean=mean;
    DDC_LST_list{DDC_LST_ind}.dy_T_mean=dy_T_mean;
    DDC_LST_list{DDC_LST_ind}.dy_S_mean=dy_S_mean;
%     DDC_LST_list{DDC_LST_ind}.mean_kolmogorov(1,1:2)=mean_kolmogorov;
    DDC_LST_list{DDC_LST_ind}.solve=solve;
    DDC_LST_list{DDC_LST_ind}.operator=operator;

    %%compute the linear stability analysis over kx kz
    DDC_LST_list{DDC_LST_ind}=DDC_LST_list{DDC_LST_ind}.solve_kxkz();
    DDC_LST_list{DDC_LST_ind}=DDC_LST_list{DDC_LST_ind}.convert_IFSC_unit_tuS();
    data{DDC_LST_ind}.x=kx_list;
    data{DDC_LST_ind}.y=real(cell2mat(DDC_LST_list{DDC_LST_ind}.eig_val_max_list));
end
plot_config.ylim_list=[1,2*min(data{1}.y)-3,2*max(data{1}.y)+3];
plot_config.label_list={1,'$k_x$','$\lambda$'};
plot_config.legend_list={1,'Primitive (IFSC unit)','IFSC','MRBC','Stokes',...
    'Primitive (IFSC unit), $\nabla^2 S=0$',...
    'IFSC, $\nabla^2 S=0$',...
     'MRBC, $\nabla^2 S=0$','Stokes, $\nabla^2 S=0$'};
plot_config.Markerindex=3;
plot_config.user_color_style_marker_list={'k-','bsquare','m*','r--'...
                                'k-','bo','mx','r--'};
plot_config.title_list={1,['$Pr$=',num2str(Pr),', $\tau$=',...
    num2str(tau),', $R_\rho$=',num2str(R_rho_T2S),...
    ', $\bar{\mathcal{T}}_z=$',num2str(dy_T_mean),...
    ', $\bar{\mathcal{S}}_z=$',num2str(dy_S_mean)]};
plot_config.name=['C:/Figure/DDC_LST/DDC_LST_Pr=',num2str(Pr),'_tau=',...
    num2str(tau),'_R_rho_T2S=',num2str(R_rho_T2S),...
    '_dy_T_mean=',num2str(dy_T_mean),...
    '_dy_S_mean=',num2str(dy_S_mean),...
    '_mean=',mean,'.png'];
plot_config.fontsize_legend=24;
plot_line(data,plot_config);






% primitive_Radko2013=DDC_LST('primitive_Radko2013');
% primitive_Radko2013.tau=tau;
% primitive_Radko2013.Pr=Pr;
% primitive_Radko2013.R_rho_T2S=R_rho_T2S;
% primitive_Radko2013=primitive_Radko2013.convert_shear();
% primitive_Radko2013_0tau=primitive_Radko2013;%make one copy but set tau=0
% primitive_Radko2013_0tau.tau=0;
% 
