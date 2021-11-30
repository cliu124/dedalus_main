clear all;
close all;
clc;

%------setup the main parameter for running
Pr=10;

% tau_list=[0.005,0.1,1,1.5];
tau_list=0.01;
% R_rho_T2S_list=1./[1.5,2,10,50];
R_rho_T2S_list=1/2;
% tau_list=0.01;
% R_rho_T2S_list=0.5;

% Old range 
% Pe_list=logspace(0,4,30);
% Ri_list=logspace(-0.61,2.31,30);

Pe_list=logspace(0,4,100);
Ri_list=logspace(-0.61,2.4,100);
dy_T_mean=-1;
dy_S_mean=-1;
mean_kolmogorov=[1,2*pi];
Lz=1;
mean='kolmogorov';
operator='v_omega_y';
Ny_full_list=62*ones(size(Pe_list));
Ny_full_list(end-length(Pe_list)/6:end)=92;
Ny_full_list(1:length(Pe_list)*2/3)=32;
% Ny_full_list=[32,62];
%This needs to be 92 for high Pe, other case 62 or 32 is enough
%Up to 0.3 and 0.4 for the Pe=100, Ri=10
%Up to 0.5 and 0.8 for the Pe=100, Ri=1 
%Up to 3.5 and 1.5 for the Pe=10^4, R1
kx_list=logspace(-1.5,0,30); %0.3
% kz_list=linspace(0,1.5,60); %0.4
kz_list=0;
solve='LST'; %%or finished if we would like to skip but just load the data..
% debug='new_Ny'; %old data: kz=0
debug='fine';
% Ri=1/((Pe/100)^2/(Pr/10));
%-----------------


%construct the obj that is using the shear flow formulation as Radko 2016..
shear_Radko2016=DDC_LST('shear_Radko2016');
shear_Radko2016.Pr=Pr;
shear_Radko2016.mean=mean;
shear_Radko2016.grid_diff{2}(2)=Lz;
shear_Radko2016.dy_T_mean=dy_T_mean;
shear_Radko2016.dy_S_mean=dy_S_mean;
shear_Radko2016.mean_kolmogorov(1,1:2)=mean_kolmogorov;
shear_Radko2016.operator=operator;
shear_Radko2016.kx_list=kx_list;
shear_Radko2016.kz_list=kz_list;
shear_Radko2016.solve=solve;
shear_Radko2016.debug=debug;
shear_Radko2016.post_eig_Radko2016=1;



for Pe_ind=1:length(Pe_list)
    Pe=Pe_list(Pe_ind)
    shear_Radko2016.Pe=Pe;
    shear_Radko2016.Ny_full=Ny_full_list(Pe_ind);
    for Ri_ind=1:length(Ri_list)
        Ri=Ri_list(Ri_ind);
        shear_Radko2016.Ri=Ri;
%         tic;

        %{
        shear_Radko2016.R_rho_T2S=0.5;
        for tau_ind=1:length(tau_list)
            shear_Radko2016.tau=tau_list(tau_ind);
            shear_Radko2016=shear_Radko2016.convert_shear();
            shear_Radko2016_tau_list{Pe_ind,Ri_ind,tau_ind}=shear_Radko2016;
            shear_Radko2016_tau_list{Pe_ind,Ri_ind,tau_ind}=...
                shear_Radko2016_tau_list{Pe_ind,Ri_ind,tau_ind}.solve_kxkz();
            data_shear_Radko2016_tau{tau_ind}.z(Pe_ind,Ri_ind)=...
                max(real(cell2mat(shear_Radko2016_tau_list{Pe_ind,Ri_ind,tau_ind}.eig_val_max_list)));

            shear_Radko2016_Stokes=shear_Radko2016;
            shear_Radko2016_Stokes.shear_Radko2016_reduced='Stokes';
            shear_Radko2016_Stokes=shear_Radko2016_Stokes.convert_shear();
            shear_Radko2016_Stokes_tau_list{Pe_ind,Ri_ind,tau_ind}=shear_Radko2016_Stokes;
            shear_Radko2016_Stokes_tau_list{Pe_ind,Ri_ind,tau_ind}=...
                shear_Radko2016_Stokes_tau_list{Pe_ind,Ri_ind,tau_ind}.solve_kxkz();
            data_shear_Radko2016_Stokes_tau{tau_ind}.z(Pe_ind,Ri_ind)=...
                max(real(cell2mat(shear_Radko2016_Stokes_tau_list{Pe_ind,Ri_ind,tau_ind}.eig_val_max_list)));
        end
        %}

        shear_Radko2016.tau=0.01;

        for R_rho_T2S_ind=1:length(R_rho_T2S_list)
            shear_Radko2016.R_rho_T2S=R_rho_T2S_list(R_rho_T2S_ind);
            shear_Radko2016=shear_Radko2016.convert_shear();
            shear_Radko2016_R_rho_T2S_list{Pe_ind,Ri_ind,R_rho_T2S_ind}=shear_Radko2016;
            shear_Radko2016_R_rho_T2S_list{Pe_ind,Ri_ind,R_rho_T2S_ind}=...
                shear_Radko2016_R_rho_T2S_list{Pe_ind,Ri_ind,R_rho_T2S_ind}.solve_kxkz();
            data_shear_Radko2016_R_rho_T2S{R_rho_T2S_ind}.z(Pe_ind,Ri_ind)=...
                max(real(cell2mat(shear_Radko2016_R_rho_T2S_list{Pe_ind,Ri_ind,R_rho_T2S_ind}.eig_val_max_list)));

            shear_Radko2016_Stokes=shear_Radko2016;
            shear_Radko2016_Stokes.shear_Radko2016_reduced='Stokes';
            shear_Radko2016_Stokes=shear_Radko2016_Stokes.convert_shear();
            shear_Radko2016_Stokes_R_rho_T2S_list{Pe_ind,Ri_ind,R_rho_T2S_ind}=shear_Radko2016_Stokes;
            shear_Radko2016_Stokes_R_rho_T2S_list{Pe_ind,Ri_ind,R_rho_T2S_ind}=...
                shear_Radko2016_Stokes_R_rho_T2S_list{Pe_ind,Ri_ind,R_rho_T2S_ind}.solve_kxkz();
            data_shear_Radko2016_Stokes_R_rho_T2S{R_rho_T2S_ind}.z(Pe_ind,Ri_ind)=...
                max(real(cell2mat(shear_Radko2016_Stokes_R_rho_T2S_list{Pe_ind,Ri_ind,R_rho_T2S_ind}.eig_val_max_list)));
        end
        
%         toc
    end
end

plot_config.loglog=[1,1];
plot_config.print_size=[1,1100,1000];
plot_config.label_list={1,'$Ri$','$Pe$'};
plot_config.xtick_list=[1,0.25,2,10,200];
plot_config.ytick_list=[1,1,10,100,1000,10000];
plot_config.zlim_list=[1,-6,0];
plot_config.ztick_list=[1,0,-1,-2,-3,-4,-5,-6];


for R_rho_T2S_ind=1:length(R_rho_T2S_list)
    R_rho_T2S=R_rho_T2S_list(R_rho_T2S_ind);
    %%set up local data so I do not need to modify the raw data...
    data{R_rho_T2S_ind}.x=Ri_list;
    data{R_rho_T2S_ind}.y=Pe_list;
    data{R_rho_T2S_ind}.z=data_shear_Radko2016_R_rho_T2S{R_rho_T2S_ind}.z;
    data{R_rho_T2S_ind}.z(find(data{R_rho_T2S_ind}.z<10^(-5)))=NaN;
    data{R_rho_T2S_ind}.z=log10(data{R_rho_T2S_ind}.z);
    plot_config.name=[shear_Radko2016.path_fig,'_Radko_2016_primitive_figure10_tau=0.01_R_rho_T2S=',num2str(R_rho_T2S),'.png'];
    plot_contour(data(R_rho_T2S_ind),plot_config);

    data{R_rho_T2S_ind}.z=data_shear_Radko2016_Stokes_R_rho_T2S{R_rho_T2S_ind}.z;
    data{R_rho_T2S_ind}.z(find(data{R_rho_T2S_ind}.z<10^(-5)))=NaN;
    data{R_rho_T2S_ind}.z=log10(data{R_rho_T2S_ind}.z);
    plot_config.name=[shear_Radko2016.path_fig,'_Radko_2016_Stokes_figure10_tau=0.01_R_rho_T2S=',num2str(R_rho_T2S),'.png'];
    plot_contour(data(R_rho_T2S_ind),plot_config);
end

for tau_ind=1:length(tau_list)
    tau=tau_list(tau_ind);
    data{tau_ind}.x=Ri_list;
    data{tau_ind}.y=Pe_list;
    data{tau_ind}.z=data_shear_Radko2016_tau{tau_ind}.z;
    data{tau_ind}.z(find(data{tau_ind}.z<10^(-5)))=NaN;
    data{tau_ind}.z=log10(data{tau_ind}.z);
    plot_config.name=[shear_Radko2016.path_fig,'_Radko_2016_primitive_figure11_R_rho_T2S=0.5_tau=',num2str(tau),'.png'];
    plot_contour(data(tau_ind),plot_config);

    data{tau_ind}.z=data_shear_Radko2016_Stokes_tau{tau_ind}.z;
    data{tau_ind}.z(find(data{tau_ind}.z<10^(-5)))=NaN;
    data{tau_ind}.z=log10(data{tau_ind}.z);
    plot_config.name=[shear_Radko2016.path_fig,'_Radko_2016_Stokes_figure11_R_rho_T2S=0.5_tau=',num2str(tau),'.png'];
    plot_contour(data(tau_ind),plot_config);
end

error('1')
%put them altogether and run 
for DDC_LST_ind=[1:4]%:length(DDC_LST_list)
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
