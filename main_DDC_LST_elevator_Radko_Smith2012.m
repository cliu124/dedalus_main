clear all;
close all;
clc;

Pr=7; tau=0.01; Pe=1; Re=1/Pr;
Ra_T=1;  
R_rho_T2S_list=1.1:0.025:2.5;
C_list=[1.5,2,3,4];
dy_T_mean=1;
dy_S_mean=1;
mean='elevator';
mean_elevator_kx='max';
% elevator_lambda_balance_bisection=[];
Ny_full=32; %This needs to be 92 for high Pe, other case 62 or 32 is enough

kx_list=linspace(0,1,20); %0.3
kz_list=0;
solve='LST'; %%or finished if we would like to skip but just load the data..
debug='RadkoSmith2012';
%-----------------


%construct the obj that is using the shear flow formulation as Radko 2016..
primitive_Radko2013=DDC_LST('unified');
primitive_Radko2013.tau=tau;
primitive_Radko2013.Pr=Pr;
primitive_Radko2013.Pe_T=Pe;
primitive_Radko2013.Pe_S=Pe;
primitive_Radko2013.Re=Re;
primitive_Radko2013.Ra_T=Ra_T;
primitive_Radko2013.mean=mean;
primitive_Radko2013.dy_T_mean=dy_T_mean;
primitive_Radko2013.dy_S_mean=dy_S_mean;
primitive_Radko2013.kx_list=kx_list;
primitive_Radko2013.kz_list=kz_list;
primitive_Radko2013.Ny_full=Ny_full;
primitive_Radko2013.solve=solve;
primitive_Radko2013.debug=debug;
% primitive_Radko2013.elevator_lambda_balance_C=elevator_lambda_balance_C;
primitive_Radko2013=primitive_Radko2013.convert_shear();
primitive_Radko2013.mean_elevator_kx=mean_elevator_kx;
%copy four of these.. but just modify their reduced formulation.

% %put them altogether and run 
% DDC_LST_ind=1;
primitive_Radko2013=primitive_Radko2013.elevator_flux();
for C_ind=1:length(C_list)
    for R_rho_T2S_ind=1:length(R_rho_T2S_list)
        R_rho_T2S=R_rho_T2S_list(R_rho_T2S_ind);
        DDC_LST_list{C_ind,R_rho_T2S_ind}=primitive_Radko2013;
        DDC_LST_list{C_ind,R_rho_T2S_ind}.Ra_S2T=DDC_LST_list{C_ind,R_rho_T2S_ind}.Ra_T/R_rho_T2S;
        DDC_LST_list{C_ind,R_rho_T2S_ind}.R_rho_T2S=R_rho_T2S;
        DDC_LST_list{C_ind,R_rho_T2S_ind}=DDC_LST_list{C_ind,R_rho_T2S_ind}.mean_profile_elevator_kx();

        DDC_LST_list{C_ind,R_rho_T2S_ind}.elevator_lambda_balance_C=C_list(C_ind);
        DDC_LST_list{C_ind,R_rho_T2S_ind}=DDC_LST_list{C_ind,R_rho_T2S_ind}.lambda_balance_bisection();
        DDC_LST_list{C_ind,R_rho_T2S_ind}=DDC_LST_list{C_ind,R_rho_T2S_ind}.elevator_flux();
        data{C_ind}.y(R_rho_T2S_ind)=DDC_LST_list{C_ind,R_rho_T2S_ind}.F_T;
    end
    data{C_ind}.x=R_rho_T2S_list;
end
lambda_balance=DDC_LST_list{1,1}.get_lambda_balance_F_T_Radko_Smith2012();
C_var_list={'C1p5','C2','C3','C4'};
for C_var_ind=1:length(C_var_list)
    data{4+C_var_ind}.x=lambda_balance.(C_var_list{C_var_ind})(:,1);
    data{4+C_var_ind}.y=lambda_balance.(C_var_list{C_var_ind})(:,2);
end
plot_config.label_list={1,'$R_\rho$','$F_T$'};
plot_config.Markerindex=3;
plot_config.user_color_style_marker_list={'k-','b--','r-.','m-'...
                                'k^','bo','rx','msquare'};
plot_line(data,plot_config);
