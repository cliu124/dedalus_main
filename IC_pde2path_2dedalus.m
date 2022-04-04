clear all;
close all;
clc

% group_name='HB_benard_salt_finger_Ra_S2T_IC';
group_name='HB_benard_salt_finger_Ra_S2T_hopf';
%This is just plot the profile for specific Ra_S2T....
branch_name_list={'tr/bpt1','tr/bpt2','tr/bpt3'};%,'tr/bpt4'
root_folder_name='C:/Data/pde2path/HB_benard_cheb4c_zonal/';
switch group_name
    case {'HB_benard_salt_finger_Ra_S2T','HB_benard_salt_finger_Ra_S2T_IC'}
        point_list=[1/90,1/80,1/70,1/60,1/50,1/40,1/30,1/20]*10^5;
        %point_list=[1/40,0.1,0.2,0.5,1,2,5,10]*10^5;
        Ra_T=10^5;
%         point_list=1300:100:4000;
%         point_list=1500:500:3000;
%         point_list=3000;
%         folder_name='pde2path_13266198/salt_finger_Ra_S2T';
        folder_name='salt_finger_Ra_S2T';
        ilam=4;
        IC_write_folder_name='./IC/tau_0p01_Ra_S2T_';
        branch_name_list={'tr/bpt1','tr/bpt2','tr/bpt3'};%,'tr/bpt4'
        Lx2d=1;
        
    case 'HB_benard_salt_finger_Ra_S2T_hopf'
        Ra_T=10^5;
        point_list=(0.4:0.1:1)*Ra_T;
        folder_name='salt_finger_Ra_S2T_hopf_2D';
        ilam=4;%:length(point_list);
        IC_write_folder_name='./IC/tau_0p33_Ra_S2T_';
%         branch_name_list={'tr/bpt1','tr/bpt2','tr/bpt3','tr/bpt2/hpt1'};%,'tr/bpt4'
        branch_name_list={'tr/bpt1/bpt3'};
        Lx2d=1;

    case 'HB_benard_salt_finger_tau'
%         point_list=1./[0.5,1/3,0.33:-0.02:0.01];
        point_list=[33.967];
%         folder_name='pde2path_13266199/salt_finger_tau';
%         folder_name='salt_finger_tau';
        folder_name='salt_finger_tau_low_Ra_S2T_Nu';
        ilam=5;
        point_plot=1;
        Lx2d=1;
        IC_write_folder_name='./IC/tau_0p01_Ra_S2T_';
        branch_name_list={'tr/bpt1','tr/bpt2','tr/bpt3'};%,'tr/bpt4'
end
my.folder_name=[root_folder_name,folder_name,'/'];
for branch_ind=1:length(branch_name_list)
    branch_name=branch_name_list{branch_ind};
    node_name=get_node_name([my.folder_name,branch_name,'/']);
    %load([my.folder_name,branch_name,'/',node_name.pt_last]);
    %data_branch{branch_ind}=p.branch;
    for pt_ind=1:length(node_name.pt_list)
        load([my.folder_name,branch_name,'/',node_name.pt_list{pt_ind}]);
        point=p.u(p.nu+ilam);
        p.sol.ineg;
        if any(abs(point_list-point)<0.01)
            ind=find(abs(point_list-point)<0.01);
            mkdir([IC_write_folder_name,num2str(round(point_list(ind)))]);
            branch_name=branch_name_list{branch_ind};
            p.my.folder_name=my.folder_name;
            p.my.plot_config.visible=0;
            p.my.plot_config.no_ylabel=1;
            p.my.plot_config.print=0;
            p.my.plot_config.post=1;
            p.my.plot_config.branch_name=branch_name;
            p.my.plot_config.point_name=node_name.pt_list{pt_ind}(1:end-4);
            p=userplot(p);
            obj_pde2path{branch_ind,ind}=p.my.obj;
            mat_pde2path=p.mat;
            
            h5_name='analysis_s1_Nx64_Nz128.h5';
            obj_dedalus{branch_ind,ind}.x_list=h5read_complex(h5_name,'/scales/x/1.0');
            obj_dedalus{branch_ind,ind}.z_list=h5read_complex(h5_name,'/scales/z/1.0');
            obj_dedalus{branch_ind,ind}.z_list_cheb=obj_dedalus{branch_ind,ind}.z_list*2-1;

            field_list={'u_tilde','w_hat',...
                'S_0','T_0','S_hat','T_hat','U_0'};
            for field_ind=1:length(field_list)
                field=field_list{field_ind};
                obj_pde2path{branch_ind,ind}.(['d_',field])=mat_pde2path.D1*obj_pde2path{branch_ind,ind}.(field);
            end
            obj_pde2path{branch_ind,ind}.p_hat=inv(mat_pde2path.D2-obj_pde2path{branch_ind,ind}.kx^2-obj_pde2path{branch_ind,ind}.ky^2)...
                *(obj_pde2path{branch_ind,ind}.Ra_T*obj_pde2path{branch_ind,ind}.T_hat-obj_pde2path{branch_ind,ind}.Ra_S2T*obj_pde2path{branch_ind,ind}.S_hat);
            field_list={'u_tilde','w_hat',...
                'S_0','T_0','S_hat','T_hat',...
                'd_u_tilde','d_w_hat',...
                'd_S_0','d_T_0','d_S_hat','d_T_hat','p_hat',...
                'U_0','d_U_0'};
            for field_ind=1:length(field_list)
                field=field_list{field_ind};
                obj_dedalus{branch_ind,ind}.(field)=chebint([0;obj_pde2path{branch_ind,ind}.(field)(:,1);0],obj_dedalus{branch_ind,ind}.z_list_cheb);
            end
            kx=obj_pde2path{branch_ind,ind}.kx;
            ky=obj_pde2path{branch_ind,ind}.ky;
            kx_2D=sqrt(kx^2+ky^2);
%             if ky_final
%                 kx=kx_final;
%             else
%                 kx=kx_final^2+;
%             end
            Lx_old=obj_dedalus{branch_ind,ind}.x_list(end)+obj_dedalus{branch_ind,ind}.x_list(2);
            
%             Lx2d=1;
            Lx_new=Lx2d*2*pi/(kx_2D);
            x=obj_dedalus{branch_ind,ind}.x_list'/Lx_old*Lx_new;

            field_list={'u','d_u','w','d_w','p','T','S','d_T','d_S'};
            for field_ind=1:length(field_list)
                field=field_list{field_ind};
                obj_dedalus{branch_ind,ind}.(field)=h5read(h5_name,['/tasks/',field]);
            end

            obj_dedalus{branch_ind,ind}.u(:,:,end)=obj_dedalus{branch_ind,ind}.U_0*ones(size(x))+2*real(1i*obj_dedalus{branch_ind,ind}.u_tilde*exp(1i*kx_2D*x));
            obj_dedalus{branch_ind,ind}.d_u(:,:,end)=obj_dedalus{branch_ind,ind}.d_U_0*ones(size(x))+2*real(1i*obj_dedalus{branch_ind,ind}.d_u_tilde*exp(1i*kx_2D*x));
            obj_dedalus{branch_ind,ind}.w(:,:,end)=2*real(obj_dedalus{branch_ind,ind}.w_hat*exp(1i*kx_2D*x));
            obj_dedalus{branch_ind,ind}.d_w(:,:,end)=2*real(obj_dedalus{branch_ind,ind}.d_w_hat*exp(1i*kx_2D*x));
            obj_dedalus{branch_ind,ind}.p(:,:,end)=2*real(obj_dedalus{branch_ind,ind}.p_hat*exp(1i*kx_2D*x));
            obj_dedalus{branch_ind,ind}.T(:,:,end)=obj_dedalus{branch_ind,ind}.T_0*ones(size(x))+2*real(obj_dedalus{branch_ind,ind}.T_hat*exp(1i*kx_2D*x));
            obj_dedalus{branch_ind,ind}.S(:,:,end)=obj_dedalus{branch_ind,ind}.S_0*ones(size(x))+2*real(obj_dedalus{branch_ind,ind}.S_hat*exp(1i*kx_2D*x));
            obj_dedalus{branch_ind,ind}.d_T(:,:,end)=obj_dedalus{branch_ind,ind}.d_T_0*ones(size(x))+2*real(obj_dedalus{branch_ind,ind}.d_T_hat*exp(1i*kx_2D*x));
            obj_dedalus{branch_ind,ind}.d_S(:,:,end)=obj_dedalus{branch_ind,ind}.d_S_0*ones(size(x))+2*real(obj_dedalus{branch_ind,ind}.d_S_hat*exp(1i*kx_2D*x));

            h5_name=['analysis_s1_Nx64_Nz128.h5'];
            % field_list={'u','d_u','w','d_w','p','T','S','d_T','d_S'};
            for field_ind=1:length(field_list)
                field=field_list{field_ind};
                h5write(h5_name,['/tasks/',field],obj_dedalus{branch_ind,ind}.(field));
            end
            h5_name_destination=[IC_write_folder_name,num2str(round(point_list(ind))),'/'...
                h5_name(1:end-3),'_',strrep(branch_name,'/','_'),'_Lx2d_',num2str(Lx2d),'.h5'];
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



