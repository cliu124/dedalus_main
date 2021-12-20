clear all;
close all;
clc;

%------setup the main parameter for running
Pr=10;
Pe=100;
tau=0.01;
R_rho_T2S=0.5;
Ri=10;
dy_T_mean=-1;
dy_S_mean=-1;
mean_kolmogorov=[1,2*pi];
Lz=1;
mean='kolmogorov';
operator='uvwpTS';
Ny_full=62; %This needs to be 92 for high Pe, other case 62 or 32 is enough
%Up to 0.3 and 0.4 for the Pe=100, Ri=10
%Up to 0.5 and 0.8 for the Pe=100, Ri=1 
%Up to 3.5 and 1.5 for the Pe=10^4, R1
% kx_list=linspace(0.005,0.2,60); %0.3
% kx_list=linspace(0.01,4,60); %0.4
kx_list=2*pi/64;
kz_list=0;
solve='LST'; %%or finished if we would like to skip but just load the data..
debug='kz=0';
% Ri=1/((Pe/100)^2/(Pr/10));
%-----------------


%construct the obj that is using the shear flow formulation as Radko 2016..
shear_Radko2016=DDC_LST('shear_Radko2016');
shear_Radko2016.tau=tau;
shear_Radko2016.Pr=Pr;
shear_Radko2016.Pe=Pe;
shear_Radko2016.Ri=Ri;
shear_Radko2016.R_rho_T2S=R_rho_T2S;
shear_Radko2016.mean=mean;
shear_Radko2016.grid_diff{2}(2)=Lz;
shear_Radko2016.dy_T_mean=dy_T_mean;
shear_Radko2016.dy_S_mean=dy_S_mean;
shear_Radko2016.mean_kolmogorov(1,1:2)=mean_kolmogorov;
shear_Radko2016.operator=operator;
% shear_Radko2016.operator='uvwpTS';
%copy four of these.. but just modify their reduced formulation.
DDC_LST_list={shear_Radko2016,shear_Radko2016,shear_Radko2016,shear_Radko2016};
DDC_LST_list{2}.shear_Radko2016_reduced='IFSC';
DDC_LST_list{3}.shear_Radko2016_reduced='MRBC';
DDC_LST_list{4}.shear_Radko2016_reduced='Stokes';

%put them altogether and run 
for DDC_LST_ind=1:length(DDC_LST_list)
%     DDC_LST_list{DDC_LST_ind}=DDC_LST_list{DDC_LST_ind}.convert_shear();
    %the parameter setting that are shared by all cases...
    DDC_LST_list{DDC_LST_ind}.kx_list=kx_list;
    DDC_LST_list{DDC_LST_ind}.kz_list=kz_list;
    DDC_LST_list{DDC_LST_ind}=DDC_LST_list{DDC_LST_ind}.convert_shear();
    DDC_LST_list{DDC_LST_ind}.Ny_full=Ny_full;
    DDC_LST_list{DDC_LST_ind}.post_eig_Radko2016=1;
    DDC_LST_list{DDC_LST_ind}.solve=solve;
    DDC_LST_list{DDC_LST_ind}.debug=debug;

    %%compute the linear stability analysis over kx kz
    DDC_LST_list{DDC_LST_ind}=DDC_LST_list{DDC_LST_ind}.solve_kxkz();
    DDC_LST_list{DDC_LST_ind}=DDC_LST_list{DDC_LST_ind}.convert_IFSC_unit_tuS();
    %DDC_LST_list{DDC_LST_ind}=DDC_LST_list{DDC_LST_ind}.post_eig_kxkz_contour();
%     DDC_LST_list{DDC_LST_ind}=DDC_LST_list{DDC_LST_ind}.post_eigenvector();

    data{DDC_LST_ind}.x=kx_list;
    data{DDC_LST_ind}.y=real(cell2mat(DDC_LST_list{DDC_LST_ind}.eig_val_max_list));
end
% error('1');
% plot_config.ylim_list=[1,2*min(data{5}.y),2*max(data{6}.y)];
plot_config.label_list={1,'$k_x$','$\lambda$'};
plot_config.legend_list={1,'Primitive','IFSC','MRBC','Stokes'};
%    'Primitive (IFSC unit), $\nabla^2 S=0$',...
%    'IFSC, $\nabla^2 S=0$',...
%     'MRBC, $\nabla^2 S=0$','Stokes, $\nabla^2 S=0$'};
plot_config.Markerindex=3;
plot_config.ylim_list=[1,2*min(data{1}.y)-0.1,2*max(data{1}.y)+0.1];
plot_config.xlim_list=[1,min(kx_list),max(kx_list)];
plot_config.user_color_style_marker_list={'k-','bsquare','m*','r--'...
                                'k-','bo','mx','r--'};
plot_config.title_list={1,['$Pr$=',num2str(Pr),', $\tau$=',...
    num2str(tau),', $R_\rho$=',num2str(R_rho_T2S),...
    ', $Pe=$',num2str(Pe),', $Ri=$',num2str(Ri),...
    ', $\bar{\mathcal{T}}_z=$',num2str(dy_T_mean),...
    ', $\bar{\mathcal{S}}_z=$',num2str(dy_S_mean)]};
plot_config.name=['C:/Figure/DDC_LST/DDC_LST_Pr=',num2str(Pr),'_tau=',...
    num2str(tau),'_R_rho_T2S=',num2str(R_rho_T2S),...
    '_Pe=',num2str(Pe),'_Ri=',num2str(Ri),...
    '_dy_T_mean=',num2str(dy_T_mean),...
    '_dy_S_mean=',num2str(dy_S_mean),...
    '_mean=',mean,'.png'];
plot_config.fontsize_legend=24;
plot_config.fontsize=36;
plot_line(data,plot_config);






% primitive_Radko2013=DDC_LST('primitive_Radko2013');
% primitive_Radko2013.tau=tau;
% primitive_Radko2013.Pr=Pr;
% primitive_Radko2013.R_rho_T2S=R_rho_T2S;
% primitive_Radko2013=primitive_Radko2013.convert_shear();
% primitive_Radko2013_0tau=primitive_Radko2013;%make one copy but set tau=0
% primitive_Radko2013_0tau.tau=0;
% 
