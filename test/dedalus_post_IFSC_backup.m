classdef dedalus_post
    %dedalus_POST Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %%Here, also set arbitrary default value to identify data type..
        %%whether number or a string. 
        
        %%These are parameters included in the flag of the IFSC 2D
        %%simulation with or without shear...
        x_list=1;
        z_list=1;
        kx_list=1;
        kz_list=1;
        Nx=1;
        Ny=1;
        Nz=1;
        Lx=1;
        Ly=1;
        Lz=1;
        t_list=1;
        h5_name='test';
        name='test';
        spectral_z='Fourier';
        flow='test';
        
        Ra_ratio=1;
        post_store_dt=20;
        stop_sim_time=2000;
        
        ks=1;
        F_sin=0;
        F_sin_2ks=0;
        F_sin_3ks=0;
        F_sin_4ks=0;
        
        phase_2ks=0;
        phase_3ks=0;
        phase_4ks=0;
        
        current_path='./';

        A_S=0;
        A_elevator=0;
        A_noise=0;
        A_shear=0;
        
        k_elevator=0; %%This is the mode for the setting elevator mode in primitive equations
        lambda_elevator=0; %%Store the growth rate of the elevator mode...
        
        
        dy_T_mean=1;
        dy_S_mean=1;
        
        %%parameter added 2021/09/13 for double diffusive in primitive
        %%equation
        R_rho_T2S=1;
        tau=1;
        Pr=1;
        
        %%This is the optimal wavenumber that I need to compute from
        %%Rayleigh ratio
        k_opt=1;

        %%flag for the post processing
        print=0;
        video=0;
        visible=0;
        
        u_laminar_normalize=1;
        
        %%Data I want to store after post-processing
        S; %%snapshot of salnity
        T; %%temperature...
        w; %%snapshot of u
        u; %%snapshot of w
        v;
        p;
        u_fluctuation; %%this minus the laminar base flow
        
        S_coeff; %%fourier coefficient of S
        w_coeff; %%fourier coefficient of u
        u_coeff; %%fourier coefficient of w
        
        E_S; %%salnity potential energy
        E_T;
        TKE_time; %%turbulence kinematic energy... I need to remove the laminar background flow of Kolmogorov type
        spectrum_TKE;
        
        flow_sub_double_diffusive_shear_2D;
        shear_Radko2016_reduced;
        
        Re=1; %#The Reynolds number appearing in front of the inertial term in momentum
        Pe_T=1; %#The Peclet number appearing in front of the inertial term in temperature
        Pe_S=1; %#The Peclet number appearing in front of the inertial term in salinity
        Ra_T=1; %#The Rayleigh number appearing in front of the temperature term, defined as Ra_T=g\alpha T_z L^4/\nu \kappa_T
        Ra_S2T=1; %#The Rayleigh number appearing in front of the salinity term, this is defined based salintiy over temperature, thus Ra_T=g\beta S_z L^4/\nu \kappa_T
        %tau=1; %This tau has been defined before %#This is the diffusivity ratio, \kappa_S/\kappa_T 
        
        
        wS;
        uS;
        uT;
        wT;
        uw;
    end
    
    methods
        function obj = dedalus_post(h5_name,flag)
            %dedalus_post Construct an instance of this class
            %   Detailed explanation goes here

            %construction function...
            obj.h5_name=h5_name;%name for the h5file
            
            %modify these flag.
            obj.print=flag.print;
            obj.video=flag.video;
            obj.visible=flag.visible;
            
            %display h5 file
            h5disp(h5_name);
            
            %read data for x, z, t, kx, kz, and compute Lx, Lz
            obj.x_list=h5read(h5_name,'/scales/x/1.0');
            obj.Nx=length(obj.x_list);
            obj.z_list=h5read(h5_name,'/scales/z/1.0');
            obj.Nz=length(obj.z_list);
            obj.t_list=h5read(h5_name,'/scales/sim_time');
            obj.kx_list=h5read(h5_name,'/scales/kx');
            obj.kz_list=h5read(h5_name,'/scales/kz');
            
            %be careful about this in computing Lx and Lz
            obj.Lx=max(obj.x_list)-min(obj.x_list)+obj.x_list(2);
            obj.Lz=max(obj.z_list)-min(obj.z_list)+obj.z_list(2);
            
            %read the flag_table.
            flag_table=readtable([h5_name(1:end-14),'flag.txt']);
            for table_ind=1:length(flag_table.x_Test)
               if isnumeric(obj.(flag_table.x_Test{table_ind}(3:end)))
                   obj.(flag_table.x_Test{table_ind}(3:end))=str2num(flag_table.x123_{table_ind}(1:end-1));
               else
                   obj.(flag_table.x_Test{table_ind}(3:end))=flag_table.x123_{table_ind}(1:end-1);
               end
            end
            
            %compute the optimal wavenumber IFSC.. This is only suitable
            %for IFSC.. see 
            obj.k_opt=(1/2*(-2-obj.Ra_ratio+sqrt(obj.Ra_ratio^2+8*obj.Ra_ratio)))^(1/4);

        end
        
        function obj=snapshot(obj,variable_name)
            %%plot the snapshot of salinity and generate video if any
            obj.S=h5read(obj.h5_name,['/tasks/',variable_name]);

            variable_max=max(max(max(obj.(variable_name))));
            variable_min=min(min(min(obj.(variable_name))));
            for t_ind=1:length(obj.t_list)
                data{1}.z=obj.(variable_name)(:,:,t_ind);
                if strcmp(obj.flow(1:7),'IFSC_2D')
                    data{1}.x=obj.x_list/(2*pi/obj.k_opt);
                    data{1}.y=obj.z_list/(2*pi/obj.k_opt);
                    plot_config.label_list={1,'$x/l_{opt}$','$z/l_{opt}$'};
                else
                    data{1}.x=obj.x_list;
                    data{1}.y=obj.z_list;
                    plot_config.label_list={1,'$x$','$z$'};
                end
                plot_config.fontsize=28;
                plot_config.zlim_list=[1,variable_min,variable_max];
                plot_config.colormap='bluewhitered';%bluewhitered
                plot_config.print_size=[1,1200,1200];
                plot_config.name=[obj.h5_name(1:end-3),'_snapshot_',variable_name,'_t_',num2str(round(obj.t_list(t_ind),2)),'.png'];
                plot_config.print=obj.print;
                plot_config.visible=obj.visible;
                snapshot(t_ind)=plot_contour(data,plot_config);
            end
            if obj.video
               plot_config.name=[obj.h5_name(1:end-3),'_snapshot_',variable_name,'_t_video.avi'];
               plot_video(snapshot,plot_config);
            end
        end
        
        function obj=spectrum_snapshot(obj,variable_name)
            %%plot the spectrum of salnity as time varies, also generate
            %%video if any
            coeff=h5read(obj.h5_name,['/tasks/',variable_name,'_coeff']);
            obj.([variable_name,'_coeff'])=coeff.r+1i*coeff.i;
            if obj.video
                for t_ind=1:length(obj.t_list)
                    clear data plot_config;
                    if strcmp(obj.flow(1:7),'IFSC_2D')
                        data{1}.x=obj.kx_list/obj.k_opt;
                        data{1}.y=obj.kz_list(1:obj.Nz/2)/obj.k_opt;
                        plot_config.zlim_list=[0,-3,0];
                        plot_config.xlim_list=[1,0,2];
                        plot_config.ylim_list=[1,0,2];
                        plot_config.xtick_list=[1,0,1,2];
                        plot_config.ytick_list=[1,0,1,2];
                        plot_config.ztick_list=[0,-3,-2,-1,0];
                        plot_config.label_list={1,'$k_x/k_{opt}$','$k_z/k_{opt}$'};
                    else
                        data{1}.x=obj.kx_list;
                        data{1}.y=obj.kz_list(1:obj.Nz/2);
                        plot_config.label_list={1,'$k_x$','$k_z$'};
                    end
                    data{1}.z=log10(abs(obj.([variable_name,'_coeff'])(1:obj.Nz/2,:,t_ind)).^2);
                    %data{2}.x=obj.kx_list/obj.k_opt;
                    %data{2}.y=obj.ks/obj.k_opt*ones(size(obj.kx_list));
                    plot_config.print_size=[1,1100,900];
                    plot_config.loglog=[0,0];
                    plot_config.print=obj.print;
                    plot_config.visible=obj.visible;
                    %plot_config.xtick_list=[0.01,0.1,1,10];

                    plot_config.colormap='white_zero';
                    plot_config.name=[obj.h5_name(1:end-3),'_spectrum_',variable_name,'_2D_t_',num2str(round(obj.t_list(t_ind),2)),'.png'];
                    frame_spectrum_2D(t_ind)=plot_contour(data,plot_config);

                    dx=diff(obj.kx_list); dx=dx(1);
                    dz=diff(obj.kz_list); dz=dz(1);
                    
                    if strcmp(obj.flow(1:7),'IFSC_2D')
                        data{1}.x=obj.kx_list/obj.k_opt;
                        data{2}.x=obj.kz_list(1:obj.Nz/2)/obj.k_opt;

                        if obj.F_sin ~=0
                            data{3}.x=2*obj.ks/obj.k_opt*ones(10,1);
                            data{3}.y=logspace(log10(min([data{1}.y,data{2}.y'])),...
                                       log10(max([data{1}.y,data{2}.y'])),10);
                        end
                        plot_config.label_list={1,'$k_x/k_{opt}$ or $k_z/k_{opt}$',''};
                    else
                        data{1}.x=obj.kx_list;
                        data{2}.x=obj.kz_list(1:obj.Nz/2);
                        plot_config.label_list={1,'$k_x$ or $k_z$',''};
                    end
                    data{1}.y=2*dz*sum(abs(obj.([variable_name,'_coeff'])(1:obj.Nz/2,:,t_ind)).^2,1);
                    data{2}.y=2*dx*sum(abs(obj.([variable_name,'_coeff'])(1:obj.Nz/2,:,t_ind)).^2,2);
                    
    %                 data{3}.y=linspace(min(data{1}.y),10);
                    plot_config.loglog=[1,1];
                    plot_config.ytick_list=[0,0.001,0.01,0.1,1,10,100,1000];
                    plot_config.ylim_list=[0];%,0.1,10];
                    plot_config.xtick_list=[1,0.001,0.01,0.1,1,10,100];
                    plot_config.legend_list={1,['$\int E_',variable_name,'(k_x,k_z)dk_z$'],['$\int E_',variable_name,'(k_x,k_z)d k_x$']};
                    plot_config.name=[obj.h5_name(1:end-3),'_spectrum_',variable_name,'_1D_t_',num2str(round(obj.t_list(t_ind),2)),'.png'];
                    plot_config.print=obj.print;
                    plot_config.visible=obj.visible;
                    frame_spectrum_1D(t_ind)=plot_line(data,plot_config);
                end
               plot_config.name=[obj.h5_name(1:end-3),'_spectrum_',variable_name,'_2D_t_video.avi'];
               plot_video(frame_spectrum_2D,plot_config);
               plot_config.name=[obj.h5_name(1:end-3),'_spectrum_',variable_name,'_1D_t_video.avi'];
               plot_video(frame_spectrum_1D,plot_config);
           end
        end
        
        function obj=spectrum_average(obj,variable_name)
            %%This function plot the 
            %%plot the overall spectrum averaged over time
            coeff=h5read(obj.h5_name,['/tasks/',variable_name,'_coeff']);
            obj.([variable_name,'_coeff'])=coeff.r+1i*coeff.i;
            if strcmp(obj.flow(1:7),'IFSC_2D')
                data{1}.x=obj.kx_list/obj.k_opt;
                data{1}.y=obj.kz_list(1:obj.Nz/2)/obj.k_opt;
                plot_config.zlim_list=[0,-3,0];
                plot_config.xlim_list=[1,0,2];
                plot_config.ylim_list=[1,0,2];
                plot_config.ztick_list=[0,-3,-2,-1,0];
                plot_config.label_list={1,'$k_x/k_{opt}$','$k_z/k_{opt}$'};
                if obj.F_sin ~=0
                    data{3}.x=2*obj.ks/obj.k_opt*ones(10,1);
                    data{3}.x=2*obj.ks*ones(10,1);
                    data{3}.y=logspace(log10(min([data{1}.y,data{2}.y'])),...
                                  log10(max([data{1}.y,data{2}.y'])),10);
                end
            else
                data{1}.x=obj.kx_list;
                data{1}.y=obj.kz_list(1:obj.Nz/2);
                plot_config.label_list={1,'$k_x$','$k_z$'};
            end
            
            
            obj.(variable_name)=h5read(obj.h5_name,['/tasks/',variable_name]);
            for t_ind=1:length(obj.t_list)
                obj.(['E_',variable_name])(t_ind)=sum(sum(obj.(variable_name)(:,:,t_ind).^2))/obj.Nx/obj.Nz/2;
            end
            [val,max_ind]=max(obj.(['E_',variable_name]));
%             t_grow=obj.t_list(1:max_ind);
            
            spectrum_average=mean(abs(obj.([variable_name,'_coeff'])(1:obj.Nz/2,:,max(3*max_ind,30):end)).^2,3);
%             obj.spectrum_S_average=spectrum_S_average;
            data{1}.z=log10(spectrum_average);
            
%             data{2}.x=obj.kx_list/obj.k_opt;
%             data{2}.y=obj.ks/obj.k_opt*ones(size(obj.kx_list));
%             plot_config.user_color_style_marker_list={'k--'};
%             


            plot_config.loglog=[0,0];
            plot_config.print_size=[1,1100,900];
            plot_config.colormap='white_zero';
            plot_config.name=[obj.h5_name(1:end-3),'_spectrum_',variable_name,'_2D_time_average.png'];
            plot_config.print=obj.print;
            plot_config.visible=obj.visible;
            plot_contour(data,plot_config);
            
            dx=diff(obj.kx_list); dx=dx(1);
            dz=diff(obj.kz_list); dz=dz(1);

%             data{1}.x=obj.kx_list/obj.k_opt;
            data{1}.x=obj.kx_list;
            data{1}.y=2*dz*sum(spectrum_average,1);
%             data{2}.x=obj.kz_list(1:obj.Nz/2)/obj.k_opt;
            data{2}.x=obj.kz_list(1:obj.Nz/2);
            data{2}.y=2*dx*sum(spectrum_average,2);
            
            plot_config.loglog=[1,1];
            plot_config.ytick_list=[0,0.001,0.01,0.1,1,10,100,1000];
            plot_config.ylim_list=[0];%,0.1,10];
%             plot_config.label_list={1,'$k_x/k_{opt}$ or $k_z/k_{opt}$',''};
            plot_config.label_list={1,'$k_x$ or $k_z$',''};
            plot_config.legend_list={1,['$\int E_',variable_name,'(k_x,k_z)dk_z$'],['$\int E_',variable_name,'(k_x,k_z)d k_x$']};
            plot_config.name=[obj.h5_name(1:end-3),'_spectrum_',variable_name,'_1D_time_average.png'];
            plot_config.print=obj.print;
            plot_config.visible=obj.visible;
            plot_line(data,plot_config);
        end
        
        
        function obj=spectrum_TKE_average(obj)
            %%This function plot the kinetic energy, (deduct the background laminar flow)
            %%This is average over time... plot the turbulence kinetic
            %%energy spectrum as kx and kz
            
            if strcmp(obj.flow,'IFSC_2D_without_shear') || strcmp(obj.flow,'IFSC_2D_with_shear') || strcmp(obj.flow,'IFSC_2D') || strcmp(obj.flow,'double_diffusive_2D')
                %%This is the post-processing for the TKE in 2D... 
                w_coeff=h5read(obj.h5_name,'/tasks/w_coeff');
                obj.w_coeff=w_coeff.r+1i*w_coeff.i;
                u_coeff=h5read(obj.h5_name,'/tasks/u_coeff');
                obj.u_coeff=u_coeff.r+1i*u_coeff.i;
                for t_ind=1:length(obj.t_list)
                    obj.spectrum_TKE(:,:,t_ind)=abs(obj.u_coeff(:,:,t_ind)).^2+abs(obj.w_coeff(:,:,t_ind)).^2;
                end
            
                obj.u=h5read(obj.h5_name,'/tasks/u');
                obj.w=h5read(obj.h5_name,'/tasks/w');

                for t_ind=1:length(obj.t_list)
                    obj.TKE_time(t_ind)=sum(sum(obj.u(:,:,t_ind).^2+obj.w(:,:,t_ind).^2))/obj.Nx/obj.Nz/2;
                end
                
                [val,max_ind]=max(obj.TKE_time);
                spectrum_TKE_average=mean(abs(obj.spectrum_TKE(1:obj.Nz/2,:,max(3*max_ind,30):end)).^2,3);
            end
            
            
            data{1}.z=log10(spectrum_TKE_average);
            data{1}.x=obj.kx_list/obj.k_opt;
            data{1}.y=obj.kz_list(1:obj.Nz/2)/obj.k_opt;
            
            plot_config.zlim_list=[0,-3,0];
            plot_config.xlim_list=[1,0,2];
            plot_config.ylim_list=[1,0,2];
            plot_config.loglog=[0,0];
            plot_config.ztick_list=[0,-3,-2,-1,0];
            plot_config.print_size=[1,1100,900];
            plot_config.label_list={1,'$k_x/k_{opt}$','$k_z/k_{opt}$'};
            plot_config.colormap='white_zero';
            plot_config.name=[obj.h5_name(1:end-3),'_spectrum_TKE_2D_time_average.png'];
            plot_config.print=obj.print;
            plot_config.visible=obj.visible;
            plot_contour(data,plot_config);
            
            dx=diff(obj.kx_list); dx=dx(1);
            dz=diff(obj.kz_list); dz=dz(1);

            data{1}.x=obj.kx_list/obj.k_opt;
            data{1}.y=2*dz*sum(spectrum_TKE_average,1);
            data{2}.x=obj.kz_list(1:obj.Nz/2)/obj.k_opt;
            data{2}.y=2*dx*sum(spectrum_TKE_average,2);
            if obj.F_sin ~=0
                data{3}.x=2*obj.ks/obj.k_opt*ones(10,1);
                data{3}.y=logspace(log10(min([data{1}.y,data{2}.y'])),...
                               log10(max([data{1}.y,data{2}.y'])),10);
            end
            plot_config.loglog=[1,1];
            plot_config.ytick_list=[0,0.001,0.01,0.1,1,10,100,1000];
            plot_config.ylim_list=[0];%,0.1,10];
            plot_config.label_list={1,'$k_x/k_{opt}$ or $k_z/k_{opt}$',''};
            plot_config.legend_list={1,'$\int E_u(k_x,k_z)dk_z$','$\int E_u(k_x,k_z)d k_x$','$2k_s/k_{opt}$'};
            plot_config.name=[obj.h5_name(1:end-3),'_spectrum_TKE_1D_time_average.png'];
            plot_config.print=obj.print;
            plot_config.visible=obj.visible;
            plot_line(data,plot_config);
        end
        
        function obj=E_time(obj,variable_name,elevator_growth_rate)
            %%Plot the salinity potential energy as a function over time
            if nargin<2 || isempty(variable_name)
                variable_name='S';
            end
            
            if nargin<3 || isempty(elevator_growth_rate)
                elevator_growth_rate=0;    
                %flag.mean='laminar_cou';   %%default value of flag_mean if not given, just set the laminar  flow.
                    %error('The flag_mean is missing.')
            end
            obj.(variable_name)=h5read(obj.h5_name,['/tasks/',variable_name]);
            
            for t_ind=1:length(obj.t_list)
                obj.(['E_',variable_name])(t_ind)=sum(sum(obj.(variable_name)(:,:,t_ind).^2))/obj.Nx/obj.Nz/2;
            end
            data{1}.x=obj.t_list;
            data{1}.y=obj.(['E_',variable_name]);
            plot_config.label_list={1,'$t$',['$E_',variable_name,'$']};
            plot_config.legend_list={0};
            if elevator_growth_rate
                [val,max_ind]=max(obj.(['E_',variable_name]));
                [~,ind_100]=min(abs(obj.t_list-100));
                t_grow=obj.t_list(1:max_ind);
                if max_ind==1
                    data{2}.x=obj.t_list;
                elseif max_ind>ind_100
                    data{2}.x=obj.t_list(1:ind_100);
                else
                    data{2}.x=t_grow;
                end
                if strcmp(obj.flow,'IFSC_2D')
                    lambda_opt=2*pi/obj.k_opt;
                    data{2}.y=obj.(['E_',variable_name])(max_ind)*exp(2*lambda_opt*(t_grow-max(t_grow)));
                    plot_config.legend_list={1,'Simulation','Linear stability'};
                elseif strcmp(obj.flow,'double_diffusive_2D')
                    k2=obj.k_elevator^2;
                    A=[-k2*obj.Pr, obj.Pr, -obj.Pr/obj.R_rho_T2S;
                        -obj.dy_T_mean, -k2, 0;
                        -obj.dy_S_mean, 0, -obj.tau*k2];
                    
                    [vec,lambda]=eig(A);
                    [val,lambda_max_ind]=max(real(diag(lambda)));
                    lambda_max=lambda(lambda_max_ind,lambda_max_ind);
                    vec_max=vec(:,lambda_max_ind);
                    switch variable_name
                        case 'S'
                             S_vec_max=vec_max(3);
                        case 'T'
                             S_vec_max=vec_max(2);
                    end
                    for t_ind=1:length(data{2}.x)
                        S2_LST=(obj.(variable_name)(:,:,1)*real(exp(lambda_max*(data{2}.x(t_ind)-data{2}.x(1)))) ...
                            -obj.(variable_name)(:,:,1)/real(S_vec_max)*imag(S_vec_max)*imag(exp(lambda_max*(data{2}.x(t_ind)-data{2}.x(1))))).^2;
                        data{2}.y(t_ind)=sum(sum(S2_LST))/obj.Nx/obj.Nz/2;
                    end
                    %data{2}.y=obj.E_S(1)*exp(2*lambda_max*(obj.t_list));
                    %data{2}.x=obj.t_list;
                elseif strcmp(obj.flow,'double_diffusive_shear_2D')
                    k2=obj.k_elevator^2;
                    A=[-k2, obj.Ra_T, -obj.Ra_S2T;
                        -obj.dy_T_mean, -k2, 0;
                        -obj.dy_S_mean, 0, -obj.tau*k2];
                    B=diag([obj.Re,obj.Pe_T,obj.Pe_S]);
                    [vec,lambda]=eig(A,B);
                    lambda(isinf(lambda)|isnan(lambda)) = -Inf;

                    [val,lambda_max_ind]=max(real(diag(lambda)));
                    lambda_max=lambda(lambda_max_ind,lambda_max_ind);
                    vec_max=vec(:,lambda_max_ind);
                    switch variable_name
                        case 'S'
                             S_vec_max=vec_max(3);
                        case 'T'
                             S_vec_max=vec_max(2);
                    end
                    for t_ind=1:length(data{2}.x)
                        S2_LST=(obj.(variable_name)(:,:,1)*real(exp(lambda_max*(data{2}.x(t_ind)-data{2}.x(1)))) ...
                            -obj.(variable_name)(:,:,1)/real(S_vec_max)*imag(S_vec_max)*imag(exp(lambda_max*(data{2}.x(t_ind)-data{2}.x(1))))).^2;
                        data{2}.y(t_ind)=sum(sum(S2_LST))/obj.Nx/obj.Nz/2;
                    end
                    plot_config.legend_list={1,'Simulation','Linear stability'};
                end
            end
            plot_config.name=[obj.h5_name(1:end-3),'_E_',variable_name,'.png'];
            plot_config.Markerindex=3;
            plot_config.user_color_style_marker_list={'k-','bo--'};
            plot_config.print=obj.print;
            plot_config.visible=obj.visible;
            plot_line(data,plot_config);
            
            plot_config.name=[obj.h5_name(1:end-3),'_E_',variable_name,'_loglog.png'];
%             plot_config.label_list={1,'$t$','$\textrm{log}_{10}(E_S)$'};
            plot_config.loglog=[0,1];
            plot_line(data,plot_config);

            if obj.video
                data{1}.x=obj.t_list;
                data{1}.y=obj.(['E_',variable_name]);
                for t_ind=1:length(obj.t_list)
                    data{2}.x=data{1}.x(t_ind);
                    data{2}.y=data{1}.y(t_ind);
                    plot_config.Markerindex=3;
                    plot_config.user_color_style_marker_list={'k-','rsquare'};

                    plot_config.fontsize=28;
                    plot_config.print_size=[1,1200,1200];
                    plot_config.print=0;
                    plot_config.visible=0;
                    plot_config.legend_list={0};
                    plot_config.loglog=[0,0];
                    plot_config.label_list={1,'$t$',['$E_',variable_name,'$']};

                    E_time(t_ind)=plot_line(data,plot_config);
                end
               plot_config.name=[obj.h5_name(1:end-3),'_E_',variable_name,'_t_video.avi'];
               plot_video(E_time,plot_config);
            end
            
        end
                
        function obj=E_TKE_time(obj)
            %%Plot the turbulence kinetic energy as a function over time
            obj=obj.u_fluctuation_read;
            obj.w=h5read(obj.h5_name,'/tasks/w');
            for t_ind=1:length(obj.t_list)
                obj.TKE_time(t_ind)=sum(sum(obj.w(:,:,t_ind).^2+obj.u_fluctuation(:,:,t_ind)))/obj.Nx/obj.Nz/2;
            end
            data{1}.x=obj.t_list;
            data{1}.y=obj.TKE_time;
            plot_config.label_list={1,'$t$','$E_u$'};
            plot_config.legend_list={0};
            plot_config.name=[obj.h5_name(1:end-3),'_E_TKE.png'];
            plot_config.Markerindex=3;
            plot_config.user_color_style_marker_list={'k-','bo--'};
            plot_config.print=obj.print;
            plot_config.visible=obj.visible;
            plot_line(data,plot_config);
        end
        
        
        function obj=u_laminar(obj)
            %%plot the laminar base flow of kolmogorov shear flow
            %%Here, also plot u', u'', and |u|'. that is the derivative of
            %%the absolute value of u... Here, normalize these quantity
            %%with the maximum magnitude as 1, if any 
            syms z;
            u_laminar=obj.F_sin/obj.ks^2*sin(obj.ks*z)...
                      +obj.F_sin_2ks/(2*obj.ks)^2*sin(2*obj.ks*z+obj.phase_2ks)...
                      +obj.F_sin_3ks/(3*obj.ks)^2*sin(3*obj.ks*z+obj.phase_3ks)...
                      +obj.F_sin_4ks/(4*obj.ks)^2*sin(4*obj.ks*z+obj.phase_4ks);
            data{1}.x=double(subs(u_laminar,z,obj.z_list));
            data{2}.x=double(subs(diff(u_laminar,z),z,obj.z_list));
            data{3}.x=double(subs(diff(diff(u_laminar,z),z),z,obj.z_list));
            data{4}.x=double(subs(diff(abs(u_laminar),z),z,obj.z_list));
            
            if strcmp(obj.flow(1:7),'IFSC_2D')
                data{1}.y=obj.z_list/(2*pi/obj.k_opt);
                data{2}.y=obj.z_list/(2*pi/obj.k_opt);
                data{3}.y=obj.z_list/(2*pi/obj.k_opt);
                data{4}.y=obj.z_list/(2*pi/obj.k_opt);
                plot_config.label_list={1,'','$z/l_{opt}$'};
            else
                data{1}.y=obj.z_list;
                data{2}.y=obj.z_list;
                data{3}.y=obj.z_list;
                data{4}.y=obj.z_list;
                plot_config.label_list={1,'','$z$'};
            end
            
            if obj.u_laminar_normalize
                data{1}.x=data{1}.x/max(abs(data{1}.x));
                data{2}.x=data{2}.x/max(abs(data{2}.x));
                data{3}.x=data{3}.x/max(abs(data{3}.x));
                data{4}.x=data{4}.x/max(abs(data{4}.x));
            end
            
            plot_config.legend_list={1,'$u_L$','$u_L''$','$u_L''''$','$|u_L|''$'};
            plot_config.print_size=[1,900,1200];
            plot_config.Markerindex=3;
            plot_config.xlim_list=[1,-1,2];
            plot_config.ylim_list=[1,min(data{1}.y),max(data{1}.y)];
            plot_config.user_color_style_marker_list={'k-','b--','r-.','bo'};
            plot_config.name=[obj.h5_name(1:end-3),'_u_laminar_normalized_',num2str(obj.u_laminar_normalize),'.png'];
            plot_line(data,plot_config);
        end
        
        function obj=u_fluctuation_read(obj)
            %%read the u fluctuation... this is done by firstly reading the
            %%total velocity and then construct the laminar base flow.
            %%Then, just remove the laminar base flow. 
            
            syms z;
            u_laminar=obj.F_sin/obj.ks^2*sin(obj.ks*z)...
                      +obj.F_sin_2ks/(2*obj.ks)^2*sin(2*obj.ks*z+obj.phase_2ks)...
                      +obj.F_sin_3ks/(3*obj.ks)^2*sin(3*obj.ks*z+obj.phase_3ks)...
                      +obj.F_sin_4ks/(4*obj.ks)^2*sin(4*obj.ks*z+obj.phase_4ks);
            u_laminar_num=double(subs(u_laminar,z,obj.z_list));
            obj.u=h5read(obj.h5_name,'/tasks/u');
            obj.u_fluctuation=obj.u;
            for z_ind=1:length(u_laminar_num)
                obj.u_fluctuation(z_ind,:,:)=obj.u(z_ind,:,:)-u_laminar_num(z_ind);
            end
        end
        
        function obj=u_fluctuation_x_ave(obj)
            %%plot the streamwise averaged u fluctuations...
            %%as a function of z (vertical axis) and time
            obj=obj.u_fluctuation_read();
%             syms z;
%             u_laminar=obj.F_sin/obj.ks^2*sin(obj.ks*z)...
%                       +obj.F_sin_2ks/(2*obj.ks)^2*sin(2*obj.ks*z+obj.phase_2ks)...
%                       +obj.F_sin_3ks/(3*obj.ks)^2*sin(3*obj.ks*z+obj.phase_3ks)...
%                       +obj.F_sin_4ks/(4*obj.ks)^2*sin(4*obj.ks*z+obj.phase_4ks);
%             u_laminar_num=double(subs(u_laminar,z,obj.z_list));
%             obj.u=h5read(obj.h5_name,'/tasks/u');
%             obj.u_fluctuation=obj.u;
%             for z_ind=1:length(u_laminar_num)
%                 obj.u_fluctuation(z_ind,:,:)=obj.u(z_ind,:,:)-u_laminar_num(z_ind);
%             end
            data{1}.x=obj.t_list;
            if strcmp(obj.flow(1:7),'IFSC_2D')
                data{1}.y=obj.z_list/(2*pi/obj.k_opt);
                plot_config.label_list={1,'$t$','$z/l_{opt}$'};
            else
                data{1}.y=obj.z_list;
                plot_config.label_list={1,'$t$','$z$'};
            end
            data{1}.z=squeeze(mean(obj.u_fluctuation,2));
            plot_config.colormap='bluewhitered';
            plot_config.print_size=[1,1200,1200];
            plot_config.print=obj.print;
            plot_config.name=[obj.h5_name(1:end-3),'_u_fluctuation_x_ave.png'];
            plot_contour(data,plot_config);
            
            data{1}.z=squeeze(mean(abs(obj.u_fluctuation),2));
            plot_config.name=[obj.h5_name(1:end-3),'_u_fluctuation_mag_x_ave.png'];
            plot_contour(data,plot_config);
        end

        
        
%         function obj=w_x_ave(obj)
%             %%plot the streamwise averaged w velocity (vertical velocity)
%             %%as a function of z (vertical axis) and time
% 
%             data{1}.x=obj.t_list;
%             if strcmp(obj.flow(1:7),'IFSC_2D')
%                 data{1}.y=obj.z_list/(2*pi/obj.k_opt);
%                 plot_config.label_list={1,'$t$','$z/l_{opt}$'};
%             else
%                 data{1}.y=obj.z_list;
%                 plot_config.label_list={1,'$t$','$z$'};
%             end
%             obj.w=h5read(obj.h5_name,'/tasks/w');
%             data{1}.z=squeeze(mean(obj.w,2));
%             plot_config.colormap='bluewhitered';
%             plot_config.print_size=[1,1200,1200];
%             plot_config.print=obj.print;
%             plot_config.name=[obj.h5_name(1:end-3),'_w_x_ave.png'];
%             plot_contour(data,plot_config);
%             
%             data{1}.z=squeeze(mean(abs(obj.w),2));
%             plot_config.name=[obj.h5_name(1:end-3),'_w_mag_x_ave.png'];
%             plot_contour(data,plot_config);
%         end
        
        
        
        function obj=x_ave(obj,variable_name)
            %%plot the streamwise averaged salnity
            %%as a function of z (vertical axis) and time

            data{1}.x=obj.t_list;
            if strcmp(obj.flow(1:7),'IFSC_2D')
                data{1}.y=obj.z_list/(2*pi/obj.k_opt);
                plot_config.label_list={1,'$t$','$z/l_{opt}$'};
            else
                data{1}.y=obj.z_list;
                plot_config.label_list={1,'$t$','$z$'};
            end
            switch variable_name
                case {'u','v','w','S','T','p'}
                    obj.(variable_name)=h5read(obj.h5_name,['/tasks/',variable_name]);
                case {'uS','wS','uT','wT','uw'}%%
                    var_1=variable_name(1);
                    var_2=variable_name(2);
                    if strcmp(var_1,'u') %%this require in default, the u is always in the first variable....
                       obj=obj.u_fluctuation_read();
                       var_1_data=obj.u_fluctuation;
                    else
                       var_1_data=h5read(obj.h5_name,['/tasks/',var_1]);
                    end
                    var_2_data=h5read(obj.h5_name,['/tasks/',var_2]);
                    obj.(variable_name)=var_1_data.*var_2_data;
            end
            data{1}.z=squeeze(mean(obj.(variable_name),2));
%             plot_config.label_list={1,'$t$','$z/l_{opt}$'};
            plot_config.colormap='bluewhitered';
            plot_config.print_size=[1,1200,1200];
            plot_config.print=obj.print;
            plot_config.name=[obj.h5_name(1:end-3),'_',variable_name,'_x_ave.png'];
            plot_contour(data,plot_config);
            
        end
%         
%         function obj=uS_x_ave(obj)
%             %%plot the streamwise averaged uS
%             %%as a function of z (vertical axis) and time
%             %%Here, u is fluctuations that need to call u_fluctuation_read.
%             
%             data{1}.x=obj.t_list;
% %             data{1}.y=obj.z_list/(2*pi/obj.k_opt);
%             if strcmp(obj.flow(1:7),'IFSC_2D')
%                 data{1}.y=obj.z_list/(2*pi/obj.k_opt);
%                 plot_config.label_list={1,'$t$','$z/l_{opt}$'};
%             else
%                 data{1}.y=obj.z_list;
%                 plot_config.label_list={1,'$t$','$z$'};
%             end
%             S=h5read(obj.h5_name,'/tasks/S');
%             obj=obj.u_fluctuation_read();
%             obj.uS=S.*obj.u_fluctuation;
%             data{1}.z=squeeze(mean(obj.uS,2));
% %             plot_config.label_list={1,'$t$','$z/l_{opt}$'};
%             plot_config.colormap='bluewhitered';
%             plot_config.print_size=[1,1200,1200];
%             plot_config.print=obj.print;
%             plot_config.name=[obj.h5_name(1:end-3),'_uS_x_ave.png'];
%             plot_contour(data,plot_config);
%         end
%         
%         function obj=wS_x_ave(obj)
%             %%plot the streamwise averaged wS
%             %%as a function of z (vertical axis) and time
%             
%             data{1}.x=obj.t_list;
% %             data{1}.y=obj.z_list/(2*pi/obj.k_opt);
%             if strcmp(obj.flow(1:7),'IFSC_2D')
%                 data{1}.y=obj.z_list/(2*pi/obj.k_opt);
%                 plot_config.label_list={1,'$t$','$z/l_{opt}$'};
%             else
%                 data{1}.y=obj.z_list;
%                 plot_config.label_list={1,'$t$','$z$'};
%             end
%             S=h5read(obj.h5_name,'/tasks/S');
%             w=h5read(obj.h5_name,'/tasks/w');
%             obj.wS=w.*S;
%             data{1}.z=squeeze(mean(obj.wS,2));
% %             plot_config.label_list={1,'$t$','$z/l_{opt}$'};
%             plot_config.colormap='bluewhitered';
%             plot_config.print_size=[1,1200,1200];
%             plot_config.print=obj.print;
%             plot_config.name=[obj.h5_name(1:end-3),'_wS_x_ave.png'];
%             plot_contour(data,plot_config);
%         end
        
        function obj=total_xt_ave(obj,variable_name)
%             data{1}.y=obj.z_list/(2*pi/obj.k_opt);
%             dz=diff(obj.z_list); dz=dz(1);
%             z_list_full=[obj.z_list;obj.z_list(end)+dz];
            switch variable_name
                case {'T','S'}
                    variable_data=h5read(obj.h5_name,['/tasks/',variable_name]);
                    data{1}.x=obj.(['dy_',variable_name,'_mean'])*obj.z_list;
                    data{2}.x=obj.(['Pe_',variable_name])*squeeze(mean(mean(variable_data,2),3))+obj.(['dy_',variable_name,'_mean'])*obj.z_list;
                    plot_config.legend_list={1,['$\bar{\mathcal{',variable_name,'}}_z z$'],['$\bar{\mathcal{',variable_name,'}}_z z+Pe_',variable_name,'\langle ',variable_name,'''\rangle_h$']};
                case {'rho'}
                    %error('not ready');
                    variable_data_T=h5read(obj.h5_name,['/tasks/T']);
                    variable_data_S=h5read(obj.h5_name,['/tasks/S']);

                    R_rho_T2S=obj.Ra_T/obj.Ra_S2T;
                    data{1}.x=-obj.dy_T_mean*obj.z_list+1/R_rho_T2S*obj.dy_S_mean*obj.z_list;
                    data{2}.x=-(obj.Pe_T*squeeze(mean(mean(variable_data_T,2),3))+obj.dy_T_mean*obj.z_list)...
                        +1/R_rho_T2S*(obj.Pe_S*squeeze(mean(mean(variable_data_S,2),3))+obj.dy_S_mean*obj.z_list);
                    plot_config.legend_list={1,['$-\bar{\mathcal{T}}_z z+R_\rho^{-1}\bar{\mathcal{S}}_z z$'],['$-(\bar{\mathcal{T}}_z z+Pe_T \langle ','T','''\rangle_h)+R_\rho^{-1}(\bar{\mathcal{S}}_z z+Pe_S \langle ','S','''\rangle_h)$']};
                    plot_config.fontsize_legend=24;
                case 'u'
                    u=h5read(obj.h5_name,'/tasks/u');
                    syms z;
                    u_laminar=obj.F_sin/obj.ks^2*sin(obj.ks*z)...
                              +obj.F_sin_2ks/(2*obj.ks)^2*sin(2*obj.ks*z+obj.phase_2ks)...
                              +obj.F_sin_3ks/(3*obj.ks)^2*sin(3*obj.ks*z+obj.phase_3ks)...
                              +obj.F_sin_4ks/(4*obj.ks)^2*sin(4*obj.ks*z+obj.phase_4ks);
                    u_laminar_num=double(subs(u_laminar,z,obj.z_list));
                    data{1}.x=u_laminar_num;
                    data{2}.x=squeeze(mean(mean(u,2),3));
                    plot_config.legend_list={1,['$\bar{',variable_name,'}$'],['$\bar{',variable_name,'}_+''\langle ',variable_name,'''\rangle_h$']};

            end
            if strcmp(obj.flow(1:7),'IFSC_2D')
                data{1}.y=obj.z_list/(2*pi/obj.k_opt);
                data{2}.y=obj.z_list/(2*pi/obj.k_opt);
                plot_config.label_list={1,'','$z/l_{opt}$'};
            else
                data{1}.y=obj.z_list;
                data{2}.y=obj.z_list;
                plot_config.label_list={1,'','$z$'};
            end
            plot_config.ylim_list=[1,min(data{1}.y),max(data{1}.y)];
%             plot_config.label_list={1,'','$z/l_{opt}$'};
            plot_config.print_size=[1,1200,1200];
            plot_config.print=obj.print;
            plot_config.name=[obj.h5_name(1:end-3),'_',variable_name,'_total_xt_ave.png'];
            plot_line(data,plot_config);
            
            plot_config.print=0;
            plot_config.visible=0;
            if obj.video
                for t_ind=1:length(obj.t_list)
                    data{2}.x=squeeze(mean(variable_data(:,:,t_ind),2))+obj.z_list;
                    snapshot(t_ind)=plot_line(data,plot_config);
                end
               plot_config.name=[obj.h5_name(1:end-3),'_',variable_name,'_total_xt_ave_video.avi'];
               plot_video(snapshot,plot_config);
            end
            
        end
        
%         function obj=T_total_xt_ave(obj)
% %             data{1}.y=obj.z_list/(2*pi/obj.k_opt);
%             T=h5read(obj.h5_name,'/tasks/T');
% %             data{1}.y=obj.z_list/(2*pi/obj.k_opt);
%             data{1}.x=obj.z_list;
% %             data{2}.y=obj.z_list/(2*pi/obj.k_opt);
%             data{2}.x=squeeze(mean(mean(T,2),3))+obj.z_list;
%             if strcmp(obj.flow(1:7),'IFSC_2D')
%                 data{1}.y=obj.z_list/(2*pi/obj.k_opt);
%                 data{2}.y=obj.z_list/(2*pi/obj.k_opt);
%                 plot_config.label_list={1,'$t$','$z/l_{opt}$'};
%             else
%                 data{1}.y=obj.z_list;
%                 data{2}.y=obj.z_list
%                 plot_config.label_list={1,'$t$','$z$'};
%             end
%             plot_config.ylim_list=[1,min(data{1}.y),max(data{1}.y)];
% %             plot_config.label_list={1,'','$z/l_{opt}$'};
%             plot_config.legend_list={1,'$\bar{T}$','$\bar{T}+T$'}
%             plot_config.print_size=[1,1200,1200];
%             plot_config.print=obj.print;
%             plot_config.name=[obj.h5_name(1:end-3),'_T_total_xt_ave.png'];
%             plot_line(data,plot_config);
%         end
%         
%         function obj=u_total_xt_ave(obj)
%             u=h5read(obj.h5_name,'/tasks/u');
%             syms z;
%             u_laminar=obj.F_sin/obj.ks^2*sin(obj.ks*z)...
%                       +obj.F_sin_2ks/(2*obj.ks)^2*sin(2*obj.ks*z+obj.phase_2ks)...
%                       +obj.F_sin_3ks/(3*obj.ks)^2*sin(3*obj.ks*z+obj.phase_3ks)...
%                       +obj.F_sin_4ks/(4*obj.ks)^2*sin(4*obj.ks*z+obj.phase_4ks);
%             u_laminar_num=double(subs(u_laminar,z,obj.z_list));
%             
%             
% %             data{1}.y=obj.z_list/(2*pi/obj.k_opt);
%             data{1}.x=u_laminar_num;
%             
% %             data{2}.y=obj.z_list/(2*pi/obj.k_opt);
%             data{2}.x=squeeze(mean(mean(u,2),3));
%             if strcmp(obj.flow(1:7),'IFSC_2D')
%                 data{1}.y=obj.z_list/(2*pi/obj.k_opt);
%                 data{2}.y=obj.z_list/(2*pi/obj.k_opt);
%                 plot_config.label_list={1,'$t$','$z/l_{opt}$'};
%             else
%                 data{1}.y=obj.z_list;
%                 data{2}.y=obj.z_list
%                 plot_config.label_list={1,'$t$','$z$'};
%             end
%             plot_config.ylim_list=[1,min(data{1}.y),max(data{1}.y)];
% %             plot_config.label_list={1,'','$z/l_{opt}$'};
%             plot_config.legend_list={1,'$\bar{U}$','$\bar{U}+u$'};
%             plot_config.print_size=[1,1200,1200];
%             plot_config.print=obj.print;
%             plot_config.name=[obj.h5_name(1:end-3),'_u_total_xt_ave.png'];
%             plot_line(data,plot_config);
%         end
        
    end
end





%         function obj=E_T_time(obj,elevator_growth_rate)
%             %%Plot the salinity potential energy as a function over time
%             
%             
%             if nargin<2 || isempty(elevator_growth_rate)
%                 elevator_growth_rate=0;    
%                 %flag.mean='laminar_cou';   %%default value of flag_mean if not given, just set the laminar  flow.
%                     %error('The flag_mean is missing.')
%             end
%             obj.T=h5read(obj.h5_name,'/tasks/T');
%             
%             for t_ind=1:length(obj.t_list)
%                 obj.E_T(t_ind)=sum(sum(obj.T(:,:,t_ind).^2))/obj.Nx/obj.Nz/2;
%             end
%             data{1}.x=obj.t_list;
%             data{1}.y=obj.E_T;
%             plot_config.label_list={1,'$t$','$E_T$'};
%             plot_config.legend_list={0};
%             if elevator_growth_rate
%                 [val,max_ind]=max(obj.E_T);
%                 [~,ind_100]=min(abs(obj.t_list-100));
%                 t_grow=obj.t_list(1:max_ind);
%                 if max_ind==1
%                     data{2}.x=obj.t_list;
%                 elseif max_ind>ind_100
%                     data{2}.x=obj.t_list(1:ind_100);
%                 else
%                     data{2}.x=t_grow;
%                 end
%                 if strcmp(obj.flow,'IFSC_2D')
%                     lambda_opt=2*pi/obj.k_opt;
%                     data{2}.y=obj.E_T(max_ind)*exp(2*lambda_opt*(t_grow-max(t_grow)));
%                     plot_config.legend_list={1,'Simulation','Linear stability'};
%                 elseif strcmp(obj.flow,'double_diffusive_2D')
%                     k2=obj.k_elevator^2;
%                     A=[-k2*obj.Pr, obj.Pr, -obj.Pr/obj.R_rho_T2S;
%                         -obj.dy_T_mean, -k2, 0;
%                         -obj.dy_S_mean, 0, -obj.tau*k2];
%                     
%                     [vec,lambda]=eig(A);
%                     [val,lambda_max_ind]=max(real(diag(lambda)));
%                     lambda_max=lambda(lambda_max_ind,lambda_max_ind);
%                     vec_max=vec(:,lambda_max_ind);
%                     T_vec_max=vec_max(2);
%                     for t_ind=1:length(data{2}.x)
%                         T2_LST=(obj.T(:,:,1)*real(exp(lambda_max*(data{2}.x(t_ind)-data{2}.x(1)))) ...
%                           -obj.T(:,:,1)/real(T_vec_max)*imag(T_vec_max)*imag(exp(lambda_max*(data{2}.x(t_ind)-data{2}.x(t_ind))))).^2;
%                         data{2}.y(t_ind)=sum(sum(T2_LST))/obj.Nx/obj.Nz/2;
%                     end
%                     %data{2}.y=obj.E_S(1)*exp(2*lambda_max*(obj.t_list));
%                     %data{2}.x=obj.t_list;
% 
%                     plot_config.legend_list={1,'Simulation','Linear stability'};
%                 end
%             end
%             plot_config.name=[obj.h5_name(1:end-3),'_E_T.png'];
%             plot_config.Markerindex=3;
%             plot_config.user_color_style_marker_list={'k-','bo--'};
%             plot_config.print=obj.print;
%             plot_config.visible=obj.visible;
%             plot_line(data,plot_config);
%             
%             plot_config.name=[obj.h5_name(1:end-3),'_E_T_loglog.png'];
% %             plot_config.label_list={1,'$t$','$\textrm{log}_{10}(E_T)$'};
%             plot_config.loglog=[0,1];
%             plot_line(data,plot_config);
% 
%             
%             data{1}.x=obj.t_list;
%             data{1}.y=obj.E_T;
%             for t_ind=1:length(obj.t_list)
%                 data{2}.x=data{1}.x(t_ind);
%                 data{2}.y=data{1}.y(t_ind);
%                 plot_config.Markerindex=3;
%                 plot_config.user_color_style_marker_list={'k-','rsquare'};
%             
%                 plot_config.fontsize=28;
%                 plot_config.print_size=[1,1200,1200];
%                 plot_config.print=0;
%                 plot_config.visible=0;
%                 plot_config.legend_list={0};
%                 E_T_time(t_ind)=plot_line(data,plot_config);
%             end
%             if obj.video
%                plot_config.name=[obj.h5_name(1:end-3),'_E_T_t_video.avi'];
%                plot_video(E_T_time,plot_config);
%             end
%             
%         end
% 


