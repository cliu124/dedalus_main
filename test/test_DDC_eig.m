clear all;
close all;
clc;

%%%This is to test the eigenvalue computation and the subcritical regime in
%%%this flow configuration... Triple periodic...
% R_rho_list=8.63;%logspace(-1,2,30); %%This is subcritical regime.
R_rho_list=[1.6736,3.47,8.36,33.47, 83.6, 167.36, 333.47, 1673.6]; %set up based on experiments
kx_list=linspace(0.01,100,300);%no need to change
%1. kx_list, set the second number to be large enough to include the peak
%when plot(kx_list,result_DDC{1}.growth_rate)
%2. 

ky_list=0;%no need to change
kz_list=[0];%logspace(-6,-1,10);%logspace(-6,0,30);
Pr_list=7;%no need to change
tau_list=0.01;%no need to change
Ra_T=4.26*10^7;%no need to change
dy_T_mean=1;
dy_S_mean=1;
for R_rho_ind=1:length(R_rho_list)
    for Pr_ind=1:length(Pr_list)
        for tau_ind=1:length(tau_list)
            for kx_ind=1:length(kx_list)
                 for ky_ind=1:length(ky_list)
                     for kz_ind=1:length(kz_list)
                         Pr=Pr_list(Pr_ind);
                         R_rho=R_rho_list(R_rho_ind);
                         tau=tau_list(tau_ind);
                         kx=kx_list(kx_ind);
                         %ky=ky_list(ky_ind);
                         ky=kx;
                         kz=kz_list(kz_ind);
                        A_D=[-(kx^2+ky^2+kz^2)*Pr, 0,inv(-(kx^2+ky^2+kz^2))*(-(kx^2+ky^2))*[Pr, -Pr/R_rho]*Ra_T;
                            0, -(kx^2+ky^2+kz^2)*Pr, 0,0;
                            -dy_T_mean, 0, -(kx^2+ky^2+kz^2), 0;
                            -dy_S_mean, 0, 0, -tau*(kx^2+ky^2+kz^2)];
                        result_DDC{R_rho_ind,Pr_ind,tau_ind}.growth_rate(kx_ind,ky_ind,kz_ind)=max(real(eig(A_D)));
                     end
                 end
            end
        end
    end
end

plot_config_kx.legend_list={1};
for R_rho_ind=1:length(R_rho_list)
    lambda_max(R_rho_ind)=max(result_DDC{R_rho_ind}.growth_rate)
    data_kx{R_rho_ind}.x=kx_list;
    data_kx{R_rho_ind}.y=result_DDC{R_rho_ind}.growth_rate;
    plot_config_kx.legend_list{1+R_rho_ind}=...
        ['$R_\rho=$',num2str(round(R_rho_list(R_rho_ind),2))];
end
data{1}.x=R_rho_list;
data{1}.y=lambda_max;
plot_config.label_list={1,'$R_\rho$','$\lambda$'};
plot_config.loglog=[0,0];
plot_config.name='Radko_2013_figure_2_2.png';
plot_line(data,plot_config);

plot_config_kx.label_list={1,'$k_x$','$\lambda$'};
plot_config_kx.Markerindex=3;
plot_config_kx.user_color_style_marker_list={'k-','b-','r-','m-','k--','b--','r--','m--'};
plot_config_kx.label_list={1,'$k_x$','$\lambda$'};
plot_config_kx.fontsize_legend=16;
plot_config_kx.fontsize=36;
plot_config_kx.name='growth_rate_kx_R_rho.png';
plot_line(data_kx,plot_config_kx);