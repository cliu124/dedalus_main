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
    end
    
    methods
        function obj = IFSC_post(folder_name,file_name)
            %IFSC_POST Construct an instance of this class
            %   Detailed explanation goes here
            obj.h5_name=[folder_name,file_name,'.h5'];
            flag_table=readtable([folder_name,'flag.txt']);
            
            h5disp(h5_name);
            obj.x_list=h5read(h5_name,'/scales/x/1.0');
            obj.Nx=length(x);
            obj.z_list=h5read(h5_name,'/scales/z/1.0');
            obj.Nz=length(z);
            obj.t_list=h5read(h5_name,'/scales/sim_time');
            obj.kx_list=h5read(h5_name,'/scales/kx');
            obj.kz_list=h5read(h5_name,'/scales/kz');
            obj.Lx=max(x)-min(x)+x(2);
            obj.Lz=max(z)-min(z)+z(2);
        end
        
        function outputArg = S_snapshot(obj)
            for t_ind=1:length(obj.t_list)
                data{1}.x=x/(2*pi/k_opt);
                data{1}.y=z/(2*pi/k_opt);
                data{1}.z=S(:,:,t_ind);
                plot_config.fontsize=28;
                plot_config.zlim_list=[1,-3.5,3.5];
                plot_config.label_list={1,'$x/l_{opt}$','$z/l_{opt}$'};
                plot_config.colormap='bluewhitered';%bluewhitered
                plot_config.print_size=[1,1200,1200];
                plot_config.name=[folder_name,file_name,'S_contour_t',num2str(round(t(t_ind),2)),'.png'];
                plot_contour(data,plot_config);
            end

        end
    end
end

