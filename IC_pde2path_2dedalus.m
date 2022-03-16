clear all;
close all;
clc

group_name='HB_benard_salt_finger_Ra_S2T_IC';
%This is just plot the profile for specific Ra_S2T....
branch_name_list={'tr/bpt1','tr/bpt2','tr/bpt3'};%,'tr/bpt4'
switch group_name
    case {'HB_benard_salt_finger_Ra_S2T','HB_benard_salt_finger_Ra_S2T_IC'}
%         point_list=[1/40,0.1,0.2,0.5,1,2,5,10]*10^5;
        Ra_T=10^5;
        point_list=1300:100:4000;
%         folder_name='pde2path_13266198/salt_finger_Ra_S2T';
        folder_name='salt_finger_Ra_S2T_IC';
        point_ind=4;
%         point_plot=[1,2];
    case 'HB_benard_salt_finger_tau'
%         point_list=1./[0.5,1/3,0.33:-0.02:0.01];
        point_list=[33.967];
%         folder_name='pde2path_13266199/salt_finger_tau';
%         folder_name='salt_finger_tau';
        folder_name='salt_finger_tau_low_Ra_S2T_Nu';
        point_ind=5;
        point_plot=1;
end
my.folder_name=['C:/Data/pde2path/HB_benard_cheb/',folder_name,'/'];
for branch_ind=[1,2,3]
    branch_name=branch_name_list{branch_ind};
    node_name=get_node_name([my.folder_name,branch_name,'/']);
    %load([my.folder_name,branch_name,'/',node_name.pt_last]);
    %data_branch{branch_ind}=p.branch;
    for pt_ind=1:length(node_name.pt_list)
        load([my.folder_name,branch_name,'/',node_name.pt_list{pt_ind}]);
        point=p.u(p.nu+point_ind);
        if any(abs(point_list-point)<0.01)
            ind=find(abs(point_list-point)<0.01);
            mkdir(['./IC/tau_0p01_Ra_S2T',num2str(point_list(ind))]);
            branch_name=branch_name_list{branch_ind};
            p.my.folder_name=my.folder_name;
            p.my.plot_config.visible=0;
            p.my.plot_config.no_ylabel=1;
            p.my.plot_config.print=0;
            p.my.plot_config.post=1;
            p.my.plot_config.branch_name=branch_name;
            p.my.plot_config.point_name=node_name.pt_list{pt_ind}(1:end-4);
            p=userplot(p);
            obj_pde2path{ind}=p.my.obj;
            mat_pde2path=p.mat;
            
            h5_name='analysis_s1_Nx64_Nz128.h5';
            obj_dedalus{ind}.x_list=h5read_complex(h5_name,'/scales/x/1.0');
            obj_dedalus{ind}.z_list=h5read_complex(h5_name,'/scales/z/1.0');
            obj_dedalus{ind}.z_list_cheb=obj_dedalus{ind}.z_list*2-1;

            field_list={'u_tilde','w_hat',...
                'S_0','T_0','S_hat','T_hat'};
            for field_ind=1:length(field_list)
                field=field_list{field_ind};
                obj_pde2path{ind}.(['d_',field])=mat_pde2path.D1*obj_pde2path{ind}.(field);
            end
            obj_pde2path{ind}.p_hat=inv(mat_pde2path.D2-obj_pde2path{ind}.kx^2-obj_pde2path{ind}.ky^2)...
                *(obj_pde2path{ind}.Ra_T*obj_pde2path{ind}.T_hat-obj_pde2path{ind}.Ra_S2T*obj_pde2path{ind}.S_hat);
            field_list={'u_tilde','w_hat',...
                'S_0','T_0','S_hat','T_hat',...
                'd_u_tilde','d_w_hat',...
                'd_S_0','d_T_0','d_S_hat','d_T_hat','p_hat'};
            for field_ind=1:length(field_list)
                field=field_list{field_ind};
                obj_dedalus{ind}.(field)=chebint([0;obj_pde2path{ind}.(field);0],obj_dedalus{ind}.z_list_cheb);
            end
            kx_final=obj_pde2path{ind}.kx;
            kx=kx_final*sqrt(2);
            Lx_old=obj_dedalus{ind}.x_list(end)+obj_dedalus{ind}.x_list(2);
            
            Lx2d=1;
            Lx_new=Lx2d*2*pi/(kx_final*sqrt(2));
            x=obj_dedalus{ind}.x_list'/Lx_old*Lx_new;

            field_list={'u','d_u','w','d_w','p','T','S','d_T','d_S'};
            for field_ind=1:length(field_list)
                field=field_list{field_ind};
                obj_dedalus{ind}.(field)=h5read(h5_name,['/tasks/',field]);
            end

            obj_dedalus{ind}.u(:,:,end)=2*real(1i*obj_dedalus{ind}.u_tilde*exp(1i*kx*x));
            obj_dedalus{ind}.d_u(:,:,end)=2*real(1i*obj_dedalus{ind}.d_u_tilde*exp(1i*kx*x));
            obj_dedalus{ind}.w(:,:,end)=2*real(obj_dedalus{ind}.w_hat*exp(1i*kx*x));
            obj_dedalus{ind}.d_w(:,:,end)=2*real(obj_dedalus{ind}.d_w_hat*exp(1i*kx*x));
            obj_dedalus{ind}.p(:,:,end)=2*real(obj_dedalus{ind}.p_hat*exp(1i*kx*x));
            obj_dedalus{ind}.T(:,:,end)=obj_dedalus{ind}.T_0*ones(size(x))+2*real(obj_dedalus{ind}.T_hat*exp(1i*kx*x));
            obj_dedalus{ind}.S(:,:,end)=obj_dedalus{ind}.S_0*ones(size(x))+2*real(obj_dedalus{ind}.S_hat*exp(1i*kx*x));
            obj_dedalus{ind}.d_T(:,:,end)=obj_dedalus{ind}.d_T_0*ones(size(x))+2*real(obj_dedalus{ind}.d_T_hat*exp(1i*kx*x));
            obj_dedalus{ind}.d_S(:,:,end)=obj_dedalus{ind}.d_S_0*ones(size(x))+2*real(obj_dedalus{ind}.d_S_hat*exp(1i*kx*x));

            h5_name=['analysis_s1_Nx64_Nz128.h5'];
            % field_list={'u','d_u','w','d_w','p','T','S','d_T','d_S'};
            for field_ind=1:length(field_list)
                field=field_list{field_ind};
                h5write(h5_name,['/tasks/',field],obj_dedalus{ind}.(field));
            end
            h5_name_destination=['./IC/tau_0p01_Ra_S2T',num2str(point_list(ind)),'/'...
                h5_name(1:end-3),strrep(branch_name,'/','_'),'_Lx2d',num2str(Lx2d),'.h5'];
            copyfile(h5_name,h5_name_destination);
        end
    end
    
end

error('1')


clear obj_dedalus;
for field_ind=1:length(field_list)
    field=field_list{field_ind};
    obj_dedalus.(field)=h5read(h5_name,['/tasks/',field]);
end

% plot(obj_dedalus.z_list,obj_dedalus.S_0,'r-'); hold on;
% plot(obj_pde2path.z_list,obj_pde2path.S_0,'b--');



