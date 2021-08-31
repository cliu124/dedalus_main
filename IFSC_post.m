classdef IFSC_post
    %IFSC_POST Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        x_list
        z_list
        kx_list
        kz_list
        Nx
        Nz
        Lx
        Lz
        t_list
        h5_name
        
        Ra_ratio
        k_opt
        F_sin
        A_S
        A_elevator
        A_noise
        A_shear
        
        S
        E_S
        S_coeff
    end
    
    methods
        function obj = IFSC_post(h5_name,flag)
            %IFSC_POST Construct an instance of this class
            %   Detailed explanation goes here
%             obj.h5_name=[folder_name,file_name,'.h5'];
%             flag_table=readtable([folder_name,'flag.txt']);
%             obj.folder_name=folder_name;

            %%construction function...
            obj.h5_name=h5_name;
            obj.Ra_ratio=flag.Ra_ratio;
            obj.k_opt=(1/2*(-2-obj.Ra_ratio+sqrt(obj.Ra_ratio^2+8*obj.Ra_ratio)))^(1/4);

            h5disp(h5_name);
            obj.x_list=h5read(h5_name,'/scales/x/1.0');
            obj.Nx=length(obj.x_list);
            obj.z_list=h5read(h5_name,'/scales/z/1.0');
            obj.Nz=length(obj.z_list);
            obj.t_list=h5read(h5_name,'/scales/sim_time');
            obj.kx_list=h5read(h5_name,'/scales/kx');
            obj.kz_list=h5read(h5_name,'/scales/kz');
            obj.Lx=max(obj.x_list)-min(obj.x_list)+obj.x_list(2);
            obj.Lz=max(obj.z_list)-min(obj.z_list)+obj.z_list(2);
        end
        
        function snapshot_S(obj)
            obj.S=h5read(obj.h5_name,'/tasks/S');

            for t_ind=1:length(obj.t_list)
                data{1}.x=obj.x_list/(2*pi/obj.k_opt);
                data{1}.y=obj.z_list/(2*pi/obj.k_opt);
                data{1}.z=obj.S(:,:,t_ind);
                plot_config.fontsize=28;
                plot_config.zlim_list=[1,-3.5,3.5];
                plot_config.label_list={1,'$x/l_{opt}$','$z/l_{opt}$'};
                plot_config.colormap='bluewhitered';%bluewhitered
                plot_config.print_size=[1,1200,1200];
                plot_config.name=[obj.h5_name(1:end-3),'_snapshot_S_t_',num2str(round(obj.t_list(t_ind),2)),'.png'];
                plot_contour(data,plot_config);
            end
        end
        
        function spectrum_S_snapshot(obj)
            S_coeff=h5read(obj.h5_name,'/tasks/S_coeff');
            obj.S_coeff=S_coeff.r+1i*S_coeff.i;
            for t_ind=1:length(obj.t_list)
                data{1}.x=obj.kx_list/obj.k_opt;
                data{1}.y=obj.kz_list(1:obj.Nz/2)/obj.k_opt;
                data{1}.z=log10(abs(obj.S_coeff(1:obj.Nz/2,:,t_ind)));
                plot_config.zlim_list=[1,-3,0];
                plot_config.xlim_list=[1,0,2];
                plot_config.ylim_list=[1,-2,2];
                plot_config.ztick_list=[1,-3,-2,-1,0];
                plot_config.print_size=[1,1100,900];
                plot_config.label_list={1,'$k/k_{opt}$','$m/k_{opt}$'};
                plot_config.colormap='white_zero';
                plot_config.name=[obj.h5_name(1:end-3),'_spectrum_S_t_',num2str(round(obj.t_list(t_ind),2)),'.png'];
                plot_contour(data,plot_config);
            end
           
        end
        
        function spectrum_S_average(obj)
            %%This function plot the 
            %%plot the overall spectrum
            S_coeff=h5read(obj.h5_name,'/tasks/S_coeff');
            obj.S_coeff=S_coeff.r+1i*S_coeff.i;
            data{1}.x=obj.kx_list/obj.k_opt;
            data{1}.y=obj.kz_list(1:obj.Nz/2)/obj.k_opt;
            
            obj.S=h5read(obj.h5_name,'/tasks/S');
            for t_ind=1:length(obj.t_list)
                obj.E_S(t_ind)=sum(sum(obj.S(:,:,t_ind).^2))/obj.Nx/obj.Nz/2;
            end
            [val,max_ind]=max(obj.E_S);
%             t_grow=obj.t_list(1:max_ind);
            
            spectrum_S_average=mean(abs(obj.S_coeff(1:obj.Nz/2,:,3*max_ind:end)),3);
            data{1}.z=log10(spectrum_S_average);
            plot_config.zlim_list=[1,-3,0];
            plot_config.xlim_list=[1,0,2];
            plot_config.ylim_list=[1,0,2];
            plot_config.ztick_list=[1,-3,-2,-1,0];
            plot_config.print_size=[1,1100,900];
            plot_config.label_list={1,'$k/k_{opt}$','$m/k_{opt}$'};
            plot_config.colormap='white_zero';
            plot_config.name=[obj.h5_name(1:end-3),'_spectrum_S_time_average.png'];
            plot_contour(data,plot_config);
            
            dx=diff(obj.kx_list); dx=dx(1);
            dz=diff(obj.kz_list); dz=dz(1);

            data{1}.x=obj.kx_list/obj.k_opt;
            data{1}.y=2*dz*sum(spectrum_S_average,1);
            data{2}.x=obj.kz_list(1:obj.Nz/2)/obj.k_opt;
            data{2}.y=2*dx*sum(spectrum_S_average,2);
            plot_config.loglog=[1,1];
            plot_config.ytick_list=[1,0.001,0.01,0.1,1,10,100,1000];
             plot_config.ylim_list=[0];%,0.1,10];
            plot_config.label_list={1,'$k_x/k_{opt}$ or $k_z/k_{opt}$',''};
            plot_config.legend_list={1,'$\int E_S(k_x,k_z)dk_z$','$\int E_S(k_x,k_z)d k_x$'};
            plot_line(data,plot_config);
        end
        function E_S_time(obj,elevator_growth_rate)
            if nargin<2 || isempty(elevator_growth_rate)
                elevator_growth_rate=0;    
                %flag.mean='laminar_cou';   %%default value of flag_mean if not given, just set the laminar  flow.
                    %error('The flag_mean is missing.')
            end
            obj.S=h5read(obj.h5_name,'/tasks/S');

            for t_ind=1:length(obj.t_list)
                obj.E_S(t_ind)=sum(sum(obj.S(:,:,t_ind).^2))/obj.Nx/obj.Nz/2;
            end
            data{1}.x=obj.t_list;
            data{1}.y=obj.E_S;
            plot_config.label_list={1,'$t$','$E_S$'};
            plot_config.legend_list={0};
            if elevator_growth_rate
                [val,max_ind]=max(obj.E_S);
                t_grow=obj.t_list(1:max_ind);
                data{2}.x=t_grow;
                lambda_opt=2*pi/obj.k_opt;
                data{2}.y=obj.E_S(max_ind)*exp(2*lambda_opt*(t_grow-max(t_grow)));
                plot_config.legend_list={1,'Simulation','Linear stability'};
            end
            plot_config.name=[obj.h5_name(1:end-3),'_E_S.png'];
            plot_config.Markerindex=3;
            plot_config.user_color_style_marker_list={'k-','bo--'};
            plot_line(data,plot_config);

            
        end
    end
end

