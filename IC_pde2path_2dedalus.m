clear all;
close all;
clc

% group_name='HB_benard_salt_finger_Ra_S2T_IC';
group_name='HB_benard_salt_finger_Ra_S2T';
group_name='HB_benard_salt_finger_kx';
%This is just plot the profile for specific Ra_S2T....
% branch_name_list={'tr/bpt1','tr/bpt2','tr/bpt3'};%,'tr/bpt4'
root_folder_name='C:/Data/pde2path/HB_benard_full_zonal/';
switch group_name
    case {'HB_benard_salt_finger_Ra_S2T','HB_benard_salt_finger_Ra_S2T_IC'}
        point_list=[1/90,1/80,1/70,1/60,1/50,1/40,1/30,1/20]*10^5;
        %point_list=[3258.1, 5542,8766.4,16289.89,33486.50,59281.41, 102272.9];
        %point_list=[1/40,0.1,0.2,0.5,1,2,5,10]*10^5;
        Ra_T=10^5;
%         point_list=1300:100:4000;
%         point_list=1500:500:3000;
%         point_list=3000;
%         folder_name='pde2path_13266198/salt_finger_Ra_S2T';
        folder_name='salt_finger_Ra_S2T_3D';
        ilam=4;
        IC_write_folder_name='./IC/3D_tau_0p01_Pr_7_kx_yang_Ra_S2T_';
        branch_name_list={'tr/bpt1'};
        %branch_name_list={'tr/bpt1','tr/bpt2','tr/bpt3'};%,'tr/bpt4'
        Lx2d=1;
        
    case 'HB_benard_salt_finger_Ra_S2T_hopf'
        Ra_T=10^5;
        point_list=(0.4:0.1:1)*Ra_T;
        folder_name='salt_finger_Ra_S2T_hopf';
        ilam=4;%:length(point_list);
        IC_write_folder_name='./IC/3D_tau_0p33_Ra_S2T_';
        branch_name_list={'tr/bpt1','tr/bpt2','tr/bpt3','tr/bpt2/hpt1'};%,'tr/bpt4'
%         branch_name_list={'tr/bpt1','tr/bpt1/bpt3'};
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
    case 'HB_benard_salt_finger_kx'
%         folder_name='RBC_Ra_low_Ra_S2T_MM_periodic';
        folder_name='salt_finger_kx_R_rho_90_MM_periodic';
        switch folder_name
            case 'salt_finger_kx_low_Ra_S2T_2D_no_slip'
                branch_name_list={'tr/bpt1','tr/bpt2','tr/bpt3','tr/bpt1/bpt1'};%,'tr/bpt4'
                %branch_name_list={'tr/bpt1/bpt1'};
                IC_write_folder_name='./IC/2D_tau_0p01_Ra_S2T_2500_Pr_7_kx_';
                h5_name='./IC/analysis_s1_Nx128_Nz128.h5';
            case 'salt_finger_kx_low_Ra_S2T_2D_stress_free'
                branch_name_list={'tr/bpt1'};%,'tr/bpt4'
                %branch_name_list={'tr/bpt1/bpt1'};
                IC_write_folder_name='./IC/stress_free_2D_tau_0p01_Ra_S2T_2500_Pr_7_kx_';
                h5_name='./IC/analysis_s1_Nx128_Nz128.h5';
            case 'salt_finger_kx_low_Ra_S2T_low_Pr_2D_no_slip'
                %branch_name_list={'tr/bpt1','tr/bpt1/bpt1','tr/bpt1/bpt2'};
                branch_name_list={'tr/bpt1/bpt2'};
                IC_write_folder_name='./IC/2D_tau_0p01_Ra_S2T_2500_Pr_0p05_kx_';
                h5_name='./IC/analysis_s1_Nx128_Nz128.h5';
            case 'salt_finger_kx_low_Ra_S2T_low_Pr_3D_no_slip'
                branch_name_list={'tr/bpt1','tr/bpt1/bpt1','tr/bpt1/bpt2'};
                IC_write_folder_name='./IC/3D_tau_0p01_Ra_S2T_2500_Pr_0p05_kx_';
                h5_name='./IC/analysis_s1_Nx128_Nz128.h5';
                
            case 'salt_finger_kx_low_Ra_S2T_2D_periodic_test'
                branch_name_list={'tr/bpt1/last_pt/last_pt_p_ds'};
                IC_write_folder_name='./IC/periodic_2D_tau_0p01_Ra_S2T_2500_Pr_7_kx_';
                h5_name='./IC/periodic_analysis_s1_Nx128_Nz128.h5';
            case 'salt_finger_kx_R_rho_10_MM_periodic'
                %error('This is from 2D continuation, need to be finished!!!');
                branch_name_list={'tr/bpt1/last_pt/last_pt_p_ds'};
                IC_write_folder_name='./IC/periodic_MM_tau_R_rho_10_Ra_T_1e5_tau_0p01_Pr_7_kx_';
                h5_name='./IC/periodic_analysis_s1_Nx128_Nz128.h5';
                point_list=[-29,-28,-27,-26,-25,-24,-23,-22,-21,-20,-19,-18];
            case 'salt_finger_kx_R_rho_40_MM_periodic'
                %error('This is from 2D continuation, need to be finished!!!');
                branch_name_list={'tr/bpt1/last_pt/last_pt_p_ds'};
                IC_write_folder_name='./IC/periodic_MM_tau_R_rho_40_Ra_T_1e5_tau_0p01_Pr_7_kx_';
                h5_name='./IC/periodic_analysis_s1_Nx128_Nz128.h5';
                point_list=[-19,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1];

            case 'salt_finger_kx_R_rho_90_MM_periodic'
                %error('This is from 2D continuation, need to be finished!!!');
                branch_name_list={'tr/bpt1/last_pt/last_pt_p_ds'};
                IC_write_folder_name='./IC/periodic_MM_tau_R_rho_90_Ra_T_1e5_tau_0p01_Pr_7_kx_';
                h5_name='./IC/periodic_analysis_s1_Nx128_Nz128.h5';
                point_list=[-10,-9,-8,-7,-6];
            case 'RBC_Ra_low_Ra_S2T_2D_periodic'
                branch_name_list={'tr/bpt1/bpt1'};
                IC_write_folder_name='./IC/flux_T_periodic_RBC_kx_2pi_Ra_T_';
                h5_name='./IC/flux_T_periodic_analysis_s1_Nx128_Nz128.h5';
            case 'RBC_kx_low_Ra_S2T_2D_periodic'
                branch_name_list={'tr/bpt1/bpt1'};
                IC_write_folder_name='./IC/flux_T_periodic_RBC_Ra_T_20000_kx_';
                h5_name='./IC/flux_T_periodic_analysis_s1_Nx128_Nz128.h5';
            case 'RBC_Ra_low_Ra_S2T_MM_periodic'
                branch_name_list={'tr/bpt1','tr/bpt1/bpt1'};
                IC_write_folder_name='./IC/flux_T_periodic_RBC_kx_10_Ra_T_';
                h5_name='./IC/flux_T_periodic_analysis_s1_Nx128_Nz128.h5';
        end
%         point_list=-0.5;
        
        %point_list=[(2:1:9)*10^4,(1:1:9)*10^5,(1:1:9)*10^6];
%         point_list=[2500,5000,7500,10000,12500,15000,17500,20000,22500,25000,27500];
%         point_list=[-18,-16,-14,-12,-10,-8,-6,-4,-2,-1];
%         point_list=-16;
%         point_list=[-18,-16,-14,-12,-10,-8,-6,-4,-2,-1];
        %point_list=[-18,-16,-14,-12,-10,-8,-6,-4,-2,-1];
%         point_list=[-0.5];
        %point_list=[-19,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6.873,-6,-5,-4,-3,-2,-1,-0.01];
        %point_list=[-12];
        ilam=1;
        %point_plot=1;
        Lx2d=1;
%         branch_name_list={'tr/bpt1'};
        %point_list=[-13,-12,-11,-10,-9,-8,-7,-6.873,-6,-5,-4,-3,-2,-1,-0.01];

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
            ind=find(abs(point_list-point)<0.1);
            ind=ind(1);
            if point_list(ind)==-0.5
               mkdir([IC_write_folder_name,'0p5']);
            else
                mkdir([IC_write_folder_name,num2str(round(abs(point_list(ind))))]);
            end
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
            
            %h5_name='analysis_s1_Nx64_Nz128.h5';
            obj_dedalus{branch_ind,ind}.x_list=h5read_complex(h5_name,'/scales/x/1.0');
            obj_dedalus{branch_ind,ind}.z_list=h5read_complex(h5_name,'/scales/z/1.0');
            obj_dedalus{branch_ind,ind}.z_list_cheb=obj_dedalus{branch_ind,ind}.z_list*2-1;
            obj_dedalus{branch_ind,ind}.z_list_fourier=obj_dedalus{branch_ind,ind}.z_list*2*pi;

            switch p.my.variable_version
                 case {'real_minimal',... only have five  variable, w_hat, T_hat, S_hat, T_0, S_0
                'complex_minimal',... only have 8 variable, add imag of the hat, w_hat, T_hat, S_hat, T_0, S_0, w_hat_imag, T_hat_imag, S_hat_imag
                'complex_zonal',... add one variable for the zonal flow, that should be the 9th variable
                'complex_zonal_with_omega_z'... add two more variable, that should be 11 variable, (most comprehensive one), w_hat, T_hat, S_hat, T_0, S_0, w_hat_imag, T_hat_imag, S_hat_imag, U_0, omega_z_hat, omega_z_hat_imag
                'real_with_omega_z',...A short version to test the hexagon planform, %This should be size variable, w_hat, T_hat, S_hat, T_0, S_0, omega_z_hat
                'complex_with_omega_z'...A short version to test the hexagon planform, but  with the complex one. This should have 10 variable, w_hat, T_hat, S_hat, T_0, S_0, w_hat_imag, T_hat_
                }

                field_list={'u_tilde','v_tilde','w_hat',...
                'S_0','T_0','S_hat','T_hat','U_0'};
                for field_ind=1:length(field_list)
                    field=field_list{field_ind};
                    obj_pde2path{branch_ind,ind}.(['d_',field])=mat_pde2path.D1*obj_pde2path{branch_ind,ind}.(field);
                end
                I=eye(size(mat_pde2path.D2));
                Laplacian=mat_pde2path.D2-obj_pde2path{branch_ind,ind}.kx^2*I-obj_pde2path{branch_ind,ind}.ky^2*I;
                Laplacian_bc=Laplacian; 
                Laplacian_bc([1,p.np],:)=0;
                obj_pde2path{branch_ind,ind}.p_hat=pinv(Laplacian_bc)...
                    *(obj_pde2path{branch_ind,ind}.Ra_T*obj_pde2path{branch_ind,ind}.d_T_hat-obj_pde2path{branch_ind,ind}.Ra_S2T*obj_pde2path{branch_ind,ind}.d_S_hat...
                    -2*1i*obj_pde2path{branch_ind,ind}.kx/obj_pde2path{branch_ind,ind}.Pr*obj_pde2path{branch_ind,ind}.w_hat.*obj_pde2path{branch_ind,ind}.d_U_0);

                %Update 2022/06/15, add the hydrostatic pressure
    %             dz=diff(obj_pde2path{branch_ind,ind}.z_list);
    %             dz=dz(1);
                for z_ind=1:length(obj_pde2path{branch_ind,ind}.z_list)
                    obj_pde2path{branch_ind,ind}.p_0(z_ind,1)=sum(mat_pde2path.Iw(1:z_ind).*(obj_pde2path{branch_ind,ind}.Ra_T*obj_pde2path{branch_ind,ind}.T_0(1:z_ind,1)-obj_pde2path{branch_ind,ind}.Ra_S2T*obj_pde2path{branch_ind,ind}.S_0(1:z_ind,1)));
                end
                p_0_int=sum(mat_pde2path.Iw.*obj_pde2path{branch_ind,ind}.p_0);
                obj_pde2path{branch_ind,ind}.p_0=obj_pde2path{branch_ind,ind}.p_0-p_0_int;
                %Update 2022/06/13, rescale the pressure to be the same unit in
                %dedalus

                field_list={'u_tilde','v_tilde','w_hat',...
                    'S_0','T_0','S_hat','T_hat',...
                    'd_u_tilde','d_v_tilde','d_w_hat',...
                    'd_S_0','d_T_0','d_S_hat','d_T_hat','p_hat','p_0'...
                    'U_0','d_U_0'};
                switch p.my.z_basis_mode
                    case 'Chebyshev'
                        for field_ind=1:length(field_list)
                            field=field_list{field_ind};
                            obj_dedalus{branch_ind,ind}.(field)=chebint([obj_pde2path{branch_ind,ind}.(field)(:,1)],obj_dedalus{branch_ind,ind}.z_list_cheb);
                        end
                    case 'Fourier'
                        for field_ind=1:length(field_list)
                            field=field_list{field_ind};
                            obj_dedalus{branch_ind,ind}.(field)=fourint([obj_pde2path{branch_ind,ind}.(field)(:,1)],obj_dedalus{branch_ind,ind}.z_list_fourier);
                        end
                    otherwise
                        error('Wrong p.my.z_basis_mode');
                end
                kx=obj_pde2path{branch_ind,ind}.kx;
                ky=obj_pde2path{branch_ind,ind}.ky;
                kx_2D=sign(kx)*sqrt(kx^2+ky^2);
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
                obj_dedalus{branch_ind,ind}.p(:,:,end)=obj_dedalus{branch_ind,ind}.p_0*ones(size(x))+2*real(obj_dedalus{branch_ind,ind}.p_hat*exp(1i*kx_2D*x));
                obj_dedalus{branch_ind,ind}.T(:,:,end)=obj_dedalus{branch_ind,ind}.T_0*ones(size(x))+2*real(obj_dedalus{branch_ind,ind}.T_hat*exp(1i*kx_2D*x));
                obj_dedalus{branch_ind,ind}.S(:,:,end)=obj_dedalus{branch_ind,ind}.S_0*ones(size(x))+2*real(obj_dedalus{branch_ind,ind}.S_hat*exp(1i*kx_2D*x));
                obj_dedalus{branch_ind,ind}.d_T(:,:,end)=obj_dedalus{branch_ind,ind}.d_T_0*ones(size(x))+2*real(obj_dedalus{branch_ind,ind}.d_T_hat*exp(1i*kx_2D*x));
                obj_dedalus{branch_ind,ind}.d_S(:,:,end)=obj_dedalus{branch_ind,ind}.d_S_0*ones(size(x))+2*real(obj_dedalus{branch_ind,ind}.d_S_hat*exp(1i*kx_2D*x));

             case {'2D_psi_T_S'}
                field_list={'u_2D','w_2D','T_2D','S_2D'};
                for field_ind=1:length(field_list)
                    field=field_list{field_ind};
                    obj_pde2path{branch_ind,ind}.(['d_',field])=mat_pde2path.D1*obj_pde2path{branch_ind,ind}.(field);
                end
                I=eye(size(mat_pde2path.D2));
                D1x=mat_pde2path.D1_2pi_2D*abs(p.u(1+p.nu));
                D1z=mat_pde2path.D1_z_2D;
                u_vel=obj_pde2path{branch_ind,ind}.u_vel;
                w_vel=obj_pde2path{branch_ind,ind}.w_vel;
                
                %source term of the pressure Poisson equation
                %This should be \partial_x u\partial_x w+\partial_x
                %w\partial_z u+\partial_z u \partial_x w+\partial_z
                %w\partial_z w
                PPE_source=diag(D1x*u_vel)*D1x*u_vel...
                    +diag(D1x*w_vel)*D1z*u_vel...
                    +diag(D1z*u_vel)*D1x*w_vel+...
                    +diag(D1z*w_vel)*D1z*w_vel;
                Pr=p.u(6+p.nu);
                obj_pde2path{branch_ind,ind}.p=pinv(mat_pde2path.Laplacian)...
                    *(obj_pde2path{branch_ind,ind}.Ra_T*D1z*obj_pde2path{branch_ind,ind}.T-obj_pde2path{branch_ind,ind}.Ra_S2T*D1z*obj_pde2path{branch_ind,ind}.S...
                    -1/Pr*PPE_source);

                %Update 2022/06/15, add the hydrostatic pressure
    %             dz=diff(obj_pde2path{branch_ind,ind}.z_list);
    %             dz=dz(1);
                for z_ind=1:length(obj_pde2path{branch_ind,ind}.z_list)
                    obj_pde2path{branch_ind,ind}.p_0(z_ind,1)=sum(diag(mat_pde2path.Iw_z(1:z_ind,1:z_ind)).*(obj_pde2path{branch_ind,ind}.Ra_T*obj_pde2path{branch_ind,ind}.T_0(1:z_ind,1)-obj_pde2path{branch_ind,ind}.Ra_S2T*obj_pde2path{branch_ind,ind}.S_0(1:z_ind,1)));
                end
                p_0_int=sum(diag(mat_pde2path.Iw_z).*obj_pde2path{branch_ind,ind}.p_0);
                obj_pde2path{branch_ind,ind}.p_0=obj_pde2path{branch_ind,ind}.p_0-p_0_int;
                %Update 2022/06/13, rescale the pressure to be the same unit in
                %dedalus
                obj_pde2path{branch_ind,ind}.p_2D=...
                    reshape(obj_pde2path{branch_ind,ind}.p,p.nz,p.nx)+obj_pde2path{branch_ind,ind}.p_0*ones(1,p.nx);

                field_list={'u_2D','w_2D','T_2D','S_2D',...
                    'd_u_2D','d_w_2D','d_T_2D','d_S_2D',...
                    'p_2D'};
                switch p.my.z_basis_mode
                    case 'Chebyshev'
                        for field_ind=1:length(field_list)
                            field=field_list{field_ind};
                            obj_dedalus{branch_ind,ind}.(field)(:,x_ind)=chebint([obj_pde2path{branch_ind,ind}.(field)(:,x_ind)],obj_dedalus{branch_ind,ind}.z_list_cheb);
                        end
                    case 'Fourier'
                        for field_ind=1:length(field_list)
                            field=field_list{field_ind};
                            for x_ind=1:p.nx
                                obj_dedalus{branch_ind,ind}.(field)(:,x_ind)=fourint([obj_pde2path{branch_ind,ind}.(field)(:,x_ind)],obj_dedalus{branch_ind,ind}.z_list_fourier);
                            end
                        end
                    otherwise
                        error('Wrong p.my.z_basis_mode');
                end
                kx=obj_pde2path{branch_ind,ind}.kx;
                %ky=obj_pde2path{branch_ind,ind}.ky;
                ky=0;
                kx_2D=sign(kx)*sqrt(kx^2+ky^2);
    %             if ky_final
    %                 kx=kx_final;
    %             else
    %                 kx=kx_final^2+;
    %             end
                Lx_old=obj_dedalus{branch_ind,ind}.x_list(end)+obj_dedalus{branch_ind,ind}.x_list(2);

    %             Lx2d=1;
                Lx_new=Lx2d*2*pi/(kx_2D);
                x_2pi=obj_dedalus{branch_ind,ind}.x_list'/Lx_old*2*pi;
                
                field_list={'u','d_u','w','d_w','p','T','S','d_T','d_S'};
                for field_ind=1:length(field_list)
                    field=field_list{field_ind};
                    obj_dedalus{branch_ind,ind}.(field)=h5read(h5_name,['/tasks/',field]);
                    switch p.my.z_basis_mode
                        case 'Fourier'
                            for z_ind=1:length(obj_dedalus{branch_ind,ind}.z_list)
                                obj_dedalus{branch_ind,ind}.(field)(z_ind,:)=fourint([obj_dedalus{branch_ind,ind}.([field,'_2D'])(z_ind,:)],x_2pi);
                            end
                        case 'Chebyshev'
                            error('Not ready');
                        otherwise
                            error('Wrong z_basis_mode');
                    end
                end

%                 obj_dedalus{branch_ind,ind}.u(:,:,end)=obj_dedalus{branch_ind,ind}.U_0*ones(size(x))+2*real(1i*obj_dedalus{branch_ind,ind}.u_tilde*exp(1i*kx_2D*x));
%                 obj_dedalus{branch_ind,ind}.d_u(:,:,end)=obj_dedalus{branch_ind,ind}.d_U_0*ones(size(x))+2*real(1i*obj_dedalus{branch_ind,ind}.d_u_tilde*exp(1i*kx_2D*x));
%                 obj_dedalus{branch_ind,ind}.w(:,:,end)=2*real(obj_dedalus{branch_ind,ind}.w_hat*exp(1i*kx_2D*x));
%                 obj_dedalus{branch_ind,ind}.d_w(:,:,end)=2*real(obj_dedalus{branch_ind,ind}.d_w_hat*exp(1i*kx_2D*x));
%                 obj_dedalus{branch_ind,ind}.p(:,:,end)=obj_dedalus{branch_ind,ind}.p_0*ones(size(x))+2*real(obj_dedalus{branch_ind,ind}.p_hat*exp(1i*kx_2D*x));
%                 obj_dedalus{branch_ind,ind}.T(:,:,end)=obj_dedalus{branch_ind,ind}.T_0*ones(size(x))+2*real(obj_dedalus{branch_ind,ind}.T_hat*exp(1i*kx_2D*x));
%                 obj_dedalus{branch_ind,ind}.S(:,:,end)=obj_dedalus{branch_ind,ind}.S_0*ones(size(x))+2*real(obj_dedalus{branch_ind,ind}.S_hat*exp(1i*kx_2D*x));
%                 obj_dedalus{branch_ind,ind}.d_T(:,:,end)=obj_dedalus{branch_ind,ind}.d_T_0*ones(size(x))+2*real(obj_dedalus{branch_ind,ind}.d_T_hat*exp(1i*kx_2D*x));
%                 obj_dedalus{branch_ind,ind}.d_S(:,:,end)=obj_dedalus{branch_ind,ind}.d_S_0*ones(size(x))+2*real(obj_dedalus{branch_ind,ind}.d_S_hat*exp(1i*kx_2D*x));

                 
                otherwise
                    error('Wrong p.my.variable_name');
         end
            
%             h5_name=['analysis_s1_Nx64_Nz128.h5'];
            % field_list={'u','d_u','w','d_w','p','T','S','d_T','d_S'};
            for field_ind=1:length(field_list)
                field=field_list{field_ind};
                h5write(h5_name,['/tasks/',field],obj_dedalus{branch_ind,ind}.(field));
            end
            
            %if fixed flux on T, also writeto the dy_T_mean
            if p.my.flux_T==1
                field='dy_T_mean_q';
                obj_dedalus{branch_ind,ind}.(field)=h5read(h5_name,['/tasks/',field]);
                obj_dedalus{branch_ind,ind}.(field)=p.my.obj.dy_T_mean*ones(size(obj_dedalus{branch_ind,ind}.(field)));
                h5write(h5_name,['/tasks/',field],obj_dedalus{branch_ind,ind}.(field));
            end
            
            %if fixed flux on S, also write the dy_S_mean
            if p.my.flux_S==1
                field='dy_S_mean_q';
                obj_dedalus{branch_ind,ind}.(field)=h5read(h5_name,['/tasks/',field]);
                obj_dedalus{branch_ind,ind}.(field)=p.my.obj.dy_S_mean*ones(size(obj_dedalus{branch_ind,ind}.(field)));
                h5write(h5_name,['/tasks/',field],obj_dedalus{branch_ind,ind}.(field));
            end
            
            if point_list(ind)==-0.5
                 h5_name_destination=[IC_write_folder_name,'0p5','/'...
                    h5_name(6:end-3),'_',strrep(branch_name,'/','_'),'_Lx2d_',num2str(Lx2d),'.h5'];            
            else
                h5_name_destination=[IC_write_folder_name,num2str(round(abs(point_list(ind)))),'/'...
                    h5_name(6:end-3),'_',strrep(branch_name,'/','_'),'_Lx2d_',num2str(Lx2d),'.h5'];
            end
            copyfile(h5_name,h5_name_destination);
            
            %{
            %write initial condition to the dedalus simulation of single
            %mode DNS 2022/05/04
            h5_name=['HB_benard_shear_analysis_s1_Nz128.h5'];
            field_list={'u_tilde','v_tilde','w_hat',...
                'S_hat','T_hat',...
                'd_u_tilde','d_v_tilde',...
                'd_S_hat','d_T_hat','p_hat',...
                };
            for field_ind=1:length(field_list)
                field=field_list{field_ind};
                %obj_dedalus{branch_ind,ind}.(field)=chebint([0;obj_pde2path{branch_ind,ind}.(field)(:,1);0],obj_dedalus{branch_ind,ind}.z_list_cheb);
                h5write(h5_name,['/tasks/',field,'_real'],real(obj_dedalus{branch_ind,ind}.(field)));
                h5write(h5_name,['/tasks/',field,'_imag'],imag(obj_dedalus{branch_ind,ind}.(field)));
            end

            field_list={'S_0','T_0',...
                'd_S_0','d_T_0',...
                'U_0','d_U_0'};
            for field_ind=1:length(field_list)
                field=field_list{field_ind};
                %obj_dedalus{branch_ind,ind}.(field)=chebint([0;obj_pde2path{branch_ind,ind}.(field)(:,1);0],obj_dedalus{branch_ind,ind}.z_list_cheb);
                h5write(h5_name,['/tasks/',field],obj_dedalus{branch_ind,ind}.(field));
            end
            
            h5_name_destination=[IC_write_folder_name,num2str(round(abs(point_list(ind)))),'/'...
                h5_name(1:end-3),'_',strrep(branch_name,'/','_'),'_Lx2d_',num2str(Lx2d),'.h5'];
            copyfile(h5_name,h5_name_destination);
            %}
            
            
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



