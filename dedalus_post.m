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
        A_elevator_imag=0;
        A_noise=0;
        A_shear=0;
        
        k_elevator=0; %%This is the mode for the setting elevator mode in primitive equations
        lambda_elevator=0; %%Store the growth rate of the elevator mode...
        
        A_secondary_T=0;
        A_secondary_S=0;
        A_secondary_w=0;
        A_secondary_U0=0;
        k_secondary=0;
        
        A_secondary_phase=0;
        
        
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
        d_S;
        d_T;
        
        S_rms_xt;
        T_rms_xt;
        w_rms_xt;
        u_rms_xt;
        v_rms_xt;
        p_rms_xt;
        
        freq;
        freq_sort; %the frequency sorted based on the spectrum in descending way...
        spec_t;
        
        S_tot;
        T_tot;
        
        
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
        ww;
        
        %new flag added 2021/11/15
        kx=1
        ky=1
        kx_2=0
        ky_2=0
        
        problem='IVP';% #This can be IVP, BVP, EVP depends on the problem you want to solve
        bvp_tolerance=1e-11;% #This is the tolerance for BVP.
        EVP_homogeneous_tolerance=1e-10;
        
        z_bc_w_left='dirichlet';% #This can be also dirichlet
        z_bc_u_v_left='dirichlet';% #This can be periodic, dirichlet, or neumann
        z_bc_T_left='dirichlet';
        z_bc_S_left='dirichlet';
        
        z_bc_w_right='dirichlet';% #This can be also dirichlet
        z_bc_u_v_right='dirichlet';% #This can be periodic, dirichlet, or neumann
        z_bc_T_right='dirichlet';
        z_bc_S_right='dirichlet';
        
        z_basis_mode='Fourier';
        timesteppers='RK443';%
        analysis=0;%
        solver=0;
        initial_dt=0.01;
        continuation=0;
        
        continuation_asymmetric=0
        %variable for harmonic balance
        w_hat;
        p_hat;
        S_hat;
        d_S_hat;
        T_hat;
        d_T_hat;
        w_hat_2;
        p_hat_2;
        S_hat_2;
        d_S_hat_2;
        T_hat_2;
        d_T_hat_2;
        T_0;
        d_T_0;
        S_0;
        d_S_0;
        
        U_0;
        d_U_0;
        
        %These two are useful for the Benard convection. 
        u_hat; %%THese are just converted from the p_hat
        u_hat_2; 

        %These are additiona u, v, and their derivative for the Benard
        %convection of harmonic balance. 
        u_tilde;
        d_u_tilde;
        v_tilde;
        d_v_tilde;
        
        %
        Nu;%Nusselt number 
        Nu_S;
        d_T_0_mid; %d_T_0 at the mid plane
        d_S_0_mid;
        T_rms_mid;
        S_rms_mid;
        u_rms_mid;
        w_rms_mid;        
        
        uvw_hewitt=1; %whether we want to convert to the unit of Hewitt..
    
        T_0_handle=0;
        S_0_handle=0;
        
        IBM_A=1000;
        IBM_z0=1/2;
        IBM_sigma=0.0001
       
        z_T_BL=0;
        z_S_BL=0;
        d_T_0_overshoot=0;
        d_S_0_overshoot=0;
        
        z_T_rms_max=0;
        z_S_rms_max=0;
        
        T_rms_max=0;
        S_rms_max=0;
        
        HB_porous_shear_phi=0;
        
        HB_porous_2_layer_Omega=0;
        
        HB_porous_3_layer_Pi=1;
        HB_porous_3_layer_h=0.01;
        
        initial_kz=2*pi
        
        %%Update, add the residue of the z_list, sG, the residue for BVP,
        %%and sGjac, and sGjac_eig
        z_list_res=0;
        G_sym=0;
        G_num=0;
        Gjac_sym=0;
        Gjac_num=0;
        Gjac_eig=0;
        Gjac_B=0;
        
        Gjac_num_bc=0;
        Gjac_B_bc=0;
        
        var_list=0;
        
        
        EVP_secondary=0;
        
        eigenvalues=0;
        eigenvectors=0;
        eigenvector_lead=struct;
        no_ylabel=0;
        
        title_time=1;
        
        store_variable='all';
        
        spec_kx_z_S=0;
        spec_kx_z_u=0;
        spec_kx_z_w=0;
        
        nx_trunc_num=0;
        nz_trunc_num=0;

        flux_T=0;
        flux_S=0;
        
        damping_1_beta=0;
        
        dy_T_mean_q=0;
        dy_S_mean_q=0;
        
        E_T_upper_bound=0;
        
        Ra_T_no_q=0;
        Ra_S2T_no_q=0;
        
        Nu_T_t=0; %Nusselt number as a function of t
    
        Nu_T_t_full=0;
        
        flow_out;
        checkpoint;
        
        S_active=1;
        
        phase_c_z=0;
        
        z_phase_diagram=0;
        
        A_w_mean=0;
        A_u_mean=0;
        
        A_w_hat=0;
        restart_t0=1;
        
        period_t=0;
        
        %u_coeff;
        %w_coeff;
        T_coeff;
        
        growth_rate;
        
        Nu_average_after;
        
        freq_local_max;
        
        Q0;
        
        Ta_sqrt_z;
    end
    
    methods
        function obj = dedalus_post(h5_name,flag)
            %dedalus_post Construct an instance of this class
            %   Detailed explanation goes here
            if nargin<2 || isempty(flag)
                flag.print=1;
                flag.video=1;
                flag.visible=1;
            end
            %construction function...
            obj.h5_name=h5_name;%name for the h5file
            
            %modify these flag.
            obj.print=flag.print;
            obj.video=flag.video;
            obj.visible=flag.visible;
            obj.no_ylabel=flag.no_ylabel;
            %display h5 file
            %h5disp(h5_name);
            
            %read the flag_table.
            flag_table=readtable([h5_name(1:end-14),'flag.txt']);
            for table_ind=1:length(flag_table.x_Test)
               if isnumeric(obj.(flag_table.x_Test{table_ind}(3:end)))
                   obj.(flag_table.x_Test{table_ind}(3:end))=str2num(flag_table.x123_{table_ind}(1:end-1));
               else
                   obj.(flag_table.x_Test{table_ind}(3:end))=flag_table.x123_{table_ind}(1:end-1);
               end
            end
            
            %add these grid points here. 
            obj.z_list=h5read_complex(h5_name,'/scales/z/1.0');
            if ~strcmp(obj.flow,'HB_benard_shear_periodic')
                obj.x_list=h5read_complex(h5_name,'/scales/x/1.0');
            end
            obj.Nz=length(obj.z_list);
            
            obj.t_list=h5read_complex(h5_name,'/scales/sim_time');
            
        end
        
        function obj = dedalus_post_ivp(obj)
            h5_name=obj.h5_name;
        %read data for x, z, t, kx, kz, and compute Lx, Lz
            obj.x_list=h5read_complex(h5_name,'/scales/x/1.0');
            obj.Nx=length(obj.x_list);
            
            obj.z_list=h5read_complex(h5_name,'/scales/z/1.0');
            obj.Nz=length(obj.z_list);
            
            obj.t_list=h5read_complex(h5_name,'/scales/sim_time');
            obj.kx_list=h5read_complex(h5_name,'/scales/kx');
            if strcmp(obj.z_basis_mode,'Fourier')
                obj.kz_list=h5read_complex(h5_name,'/scales/kz');
            end
            %be careful about this in computing Lx and Lz
            obj.Lx=max(obj.x_list)-min(obj.x_list)+obj.x_list(2);
            obj.Lz=max(obj.z_list)-min(obj.z_list)+obj.z_list(2);
            
            
            %compute the optimal wavenumber IFSC.. This is only suitable
            %for IFSC.. see 
            obj.k_opt=(1/2*(-2-obj.Ra_ratio+sqrt(obj.Ra_ratio^2+8*obj.Ra_ratio)))^(1/4);
        end
        
        function obj = dedalus_post_bvp(obj)
            h5_name=obj.h5_name;
            
            %%Read the eigenvalue problem if performed the secondary
            %%stability analysis
            if obj.EVP_secondary==1 & (strcmp(obj.problem,'BVP') | strcmp(obj.problem,'EVP'))
                obj.eigenvalues=h5read_complex(h5_name,'/eigenvalues');
                obj.eigenvectors=h5read_complex(h5_name,'/eigenvectors');
                obj.z_list=h5read_complex(h5_name,'/scales/z/1.0');
%                 eigenvector_lead=obj.eigenvectors(1,:);
                if size(obj.eigenvectors,2)==1
                   obj.eigenvectors=obj.eigenvectors'; 
                end
                if strcmp(obj.flow,'HB_benard')
                    obj.eigenvector_lead.u_tilde=obj.eigenvectors(1,1:obj.Nz);
                    obj.eigenvector_lead.d_u_tilde=obj.eigenvectors(1,obj.Nz+1:2*obj.Nz);
                    obj.eigenvector_lead.v_tilde=obj.eigenvectors(1,2*obj.Nz+1:3*obj.Nz);
                    obj.eigenvector_lead.d_v_tilde=obj.eigenvectors(1,3*obj.Nz+1:4*obj.Nz);
                    obj.eigenvector_lead.w_hat=obj.eigenvectors(1,4*obj.Nz+1:5*obj.Nz);
                    obj.eigenvector_lead.p_hat=obj.eigenvectors(1,5*obj.Nz+1:6*obj.Nz);
                    obj.eigenvector_lead.T_hat=obj.eigenvectors(1,6*obj.Nz+1:7*obj.Nz);
                    obj.eigenvector_lead.d_T_hat=obj.eigenvectors(1,7*obj.Nz+1:8*obj.Nz);
                    obj.eigenvector_lead.S_hat=obj.eigenvectors(1,8*obj.Nz+1:9*obj.Nz);
                    obj.eigenvector_lead.d_S_hat=obj.eigenvectors(1,9*obj.Nz+1:10*obj.Nz);
                    obj.eigenvector_lead.T_0=obj.eigenvectors(1,10*obj.Nz+1:11*obj.Nz);
                    obj.eigenvector_lead.d_T_0=obj.eigenvectors(1,11*obj.Nz+1:12*obj.Nz);
                    obj.eigenvector_lead.S_0=obj.eigenvectors(1,12*obj.Nz+1:13*obj.Nz);
                    obj.eigenvector_lead.d_S_0=obj.eigenvectors(1,13*obj.Nz+1:14*obj.Nz);

                elseif strcmp(obj.flow,'HB_porous')

                end
            end
            
            
            %read the data for BVP... These are results for harmonic
            %balance...
            obj.z_list=h5read_complex(h5_name,'/scales/z/1.0');
            obj.w_hat=h5read_complex(h5_name,'/tasks/w_hat');
            obj.p_hat=h5read_complex(h5_name,'/tasks/p_hat');
            obj.T_hat=h5read_complex(h5_name,'/tasks/T_hat');
            obj.d_T_hat=h5read_complex(h5_name,'/tasks/d_T_hat');
            obj.S_hat=h5read_complex(h5_name,'/tasks/S_hat');
            obj.d_S_hat=h5read_complex(h5_name,'/tasks/d_S_hat');
            obj.T_0=h5read_complex(h5_name,'/tasks/T_0');
            obj.d_T_0=h5read_complex(h5_name,'/tasks/d_T_0');
            obj.S_0=h5read_complex(h5_name,'/tasks/S_0');
            obj.d_S_0=h5read_complex(h5_name,'/tasks/d_S_0');
            try 
                obj.U_0=h5read_complex(h5_name,'/tasks/U_0');
            end
            
            switch obj.flow
                case 'HB_porous'
                    obj.u_tilde=-obj.kx*obj.p_hat;
                case {'HB_benard','HB_benard_shear'}
                    obj.u_tilde=h5read_complex(h5_name,'/tasks/u_tilde');
                    obj.v_tilde=h5read_complex(h5_name,'/tasks/v_tilde');
                    obj.d_u_tilde=h5read_complex(h5_name,'/tasks/d_u_tilde');
                    obj.d_v_tilde=h5read_complex(h5_name,'/tasks/d_v_tilde');
            end
            
            %read the second harmonic data if any...
            if obj.kx_2~= 0 || obj.ky_2~=0
                obj.w_hat_2=h5read_complex(h5_name,'/tasks/w_hat_2');
                if obj.flow=='HB_porous'
                    obj.p_hat=h5read_complex(h5_name,'/tasks/p_hat_2');
                elseif obj.flow=='HB_benard'
                    error('This is not supported yet');
                    obj.u_tilde_2=h5read_complex(h5_name,'/tasks/u_tilde_2');
                    obj.v_tilde_2=h5read_complex(h5_name,'/tasks/v_tilde_2');
                end
                obj.T_hat_2=h5read_complex(h5_name,'/tasks/T_hat_2');
                obj.d_T_hat_2=h5read_complex(h5_name,'/tasks/d_T_hat_2');
                obj.S_hat_2=h5read_complex(h5_name,'/tasks/S_hat_2');
                obj.d_S_hat_2=h5read_complex(h5_name,'/tasks/d_S_hat_2');
            end
            
            
            if strcmp(obj.flow,'HB_porous_2_layer')
                obj.z_list=h5read_complex(h5_name,'/scales/z/1.0');
                obj.z_list=[obj.z_list;obj.z_list+0.5];
                obj.Nz=2*obj.Nz;
                obj.w_hat=[obj.w_hat;h5read_complex(h5_name,'/tasks/w_hat_top')];
                obj.p_hat=[obj.p_hat;h5read_complex(h5_name,'/tasks/p_hat_top')];
                obj.T_hat=[obj.T_hat;h5read_complex(h5_name,'/tasks/T_hat_top')];
                obj.d_T_hat=[obj.d_T_hat;h5read_complex(h5_name,'/tasks/d_T_hat_top')];
                obj.S_hat=[obj.S_hat;h5read_complex(h5_name,'/tasks/S_hat_top')];
                obj.d_S_hat=[obj.d_S_hat;h5read_complex(h5_name,'/tasks/d_S_hat_top')];
                obj.T_0=[obj.T_0;h5read_complex(h5_name,'/tasks/T_0_top')];
                obj.d_T_0=[obj.d_T_0;h5read_complex(h5_name,'/tasks/d_T_0_top')];
                obj.S_0=[obj.S_0;h5read_complex(h5_name,'/tasks/S_0_top')];
                obj.d_S_0=[obj.d_S_0;h5read_complex(h5_name,'/tasks/d_S_0_top')];
                obj.u_tilde=-obj.kx*obj.p_hat;
            elseif strcmp(obj.flow,'HB_porous_3_layer')
                h=obj.HB_porous_3_layer_h; Lz=obj.Lz;
                z_list=h5read_complex(h5_name,'/scales/z/1.0');
                obj.z_list=[(1-h)/2/Lz*z_list;(1-h)/2+h/Lz*z_list;(1+h)/2+(1-h)/2/Lz*z_list];
                obj.Nz=3*obj.Nz;
                
                obj.w_hat=[obj.w_hat;h5read_complex(h5_name,'/tasks/w_hat_mid');h5read_complex(h5_name,'/tasks/w_hat_top')];
                obj.p_hat=[obj.p_hat;h5read_complex(h5_name,'/tasks/p_hat_mid');h5read_complex(h5_name,'/tasks/p_hat_top')];
                obj.T_hat=[obj.T_hat;h5read_complex(h5_name,'/tasks/T_hat_mid');h5read_complex(h5_name,'/tasks/T_hat_top')];
                obj.d_T_hat=[obj.d_T_hat;h5read_complex(h5_name,'/tasks/d_T_hat_mid');h5read_complex(h5_name,'/tasks/d_T_hat_top')];
                obj.S_hat=[obj.S_hat;h5read_complex(h5_name,'/tasks/S_hat_mid');h5read_complex(h5_name,'/tasks/S_hat_top')];
                obj.d_S_hat=[obj.d_S_hat;h5read_complex(h5_name,'/tasks/d_S_hat_mid');h5read_complex(h5_name,'/tasks/d_S_hat_top')];
                obj.T_0=[obj.T_0;h5read_complex(h5_name,'/tasks/T_0_mid');h5read_complex(h5_name,'/tasks/T_0_top')];
                obj.d_T_0=[obj.d_T_0;h5read_complex(h5_name,'/tasks/d_T_0_mid');h5read_complex(h5_name,'/tasks/d_T_0_top')];
                obj.S_0=[obj.S_0;h5read_complex(h5_name,'/tasks/S_0_mid');h5read_complex(h5_name,'/tasks/S_0_top')];
                obj.d_S_0=[obj.d_S_0;h5read_complex(h5_name,'/tasks/d_S_0_mid');h5read_complex(h5_name,'/tasks/d_S_0_top')];
                
                Pi=obj.HB_porous_3_layer_Pi
                obj.u_tilde=[-obj.kx*h5read_complex(h5_name,'/tasks/p_hat');-obj.kx*Pi*h5read_complex(h5_name,'/tasks/p_hat_mid');-obj.kx*h5read_complex(h5_name,'/tasks/p_hat_top')];
            end
            
            
            if obj.uvw_hewitt & (strcmp(obj.flow,'HB_porous') | strcmp(obj.flow,'HB_porous_2_layer'))
                obj.w_hat=obj.w_hat/obj.Ra_T;
                obj.p_hat=obj.p_hat/obj.Ra_T;
                obj.u_tilde=obj.u_tilde/obj.Ra_T;
                obj.v_tilde=obj.v_tilde/obj.Ra_T;
                obj.w_hat_2=obj.w_hat_2/obj.Ra_T;
                obj.p_hat_2=obj.p_hat_2/obj.Ra_T;
            end
            
            mid_ind=obj.Nz/2;
            obj.Nu=-(obj.d_T_0+obj.dy_T_mean);%Nusselt number 
            obj.Nu_S=-(obj.d_S_0+obj.dy_S_mean); %nussel number for salinity... also add the background one...
%             [T_BL_ind]=find(diff(sign(obj.d_T_0)));%min(abs(obj.d_T_0+obj.dy_T_mean));
%             if ~isempty(T_BL_ind)
%                 obj.z_T_BL=obj.z_list(T_BL_ind(1));
%             end
            [obj.d_T_0_overshoot,d_T_0_max_ind]=max(obj.d_T_0+obj.dy_T_mean);
            obj.z_T_BL=obj.z_list(d_T_0_max_ind);
            if obj.z_T_BL>0.5
                obj.z_T_BL=1-obj.z_T_BL;
            end
            
%             [S_BL_ind]=find(diff(sign(obj.d_S_0)));
%             if ~isempty(S_BL_ind)
%                 obj.z_S_BL=obj.z_list(S_BL_ind(1));
%             end
            [obj.d_S_0_overshoot,d_S_0_max_ind]=max(obj.d_T_0+obj.dy_S_mean);
            obj.z_S_BL=obj.z_list(d_S_0_max_ind);
            if obj.z_S_BL>0.5
                obj.z_S_BL=1-obj.z_S_BL;
            end
            
            obj.T_rms_max=max(obj.T_hat*sqrt(2));
            [~,max_ind]=max(obj.T_hat);
            obj.z_T_rms_max=obj.z_list(max_ind(1));
            if obj.z_T_rms_max>0.5
                obj.z_T_rms_max=1-obj.z_T_rms_max;
            end
            
            obj.S_rms_max=max(obj.S_hat*sqrt(2));
            [~,max_ind]=max(obj.S_hat);
            obj.z_S_rms_max=obj.z_list(max_ind(1));
            if obj.z_S_rms_max>0.5
                obj.z_S_rms_max=1-obj.z_S_rms_max;
            end
            
            obj.d_T_0_mid=obj.d_T_0(mid_ind)+obj.dy_T_mean; %d_T_0 at the mid plane
            obj.d_S_0_mid=obj.d_S_0(mid_ind)+obj.dy_S_mean;
            obj.T_rms_mid=obj.T_hat(mid_ind)*sqrt(2);
            obj.S_rms_mid=obj.S_hat(mid_ind)*sqrt(2);
            obj.w_rms_mid=obj.w_hat(mid_ind)*sqrt(2);
            obj.u_rms_mid=abs(obj.u_tilde(mid_ind))*sqrt(2);
            switch obj.flow
                case 'HB_porous'
                    obj.u_rms_mid=obj.kx*abs(obj.p_hat(mid_ind))*sqrt(2);
                case {'HB_benard','HB_benard_shear'}
                    obj.u_rms_mid=abs(obj.u_tilde(mid_ind))*sqrt(2);
            end
            
            try
                disp('Also get the time series');
                obj.t_list=h5read_complex(h5_name,'/scales/sim_time');
            catch
                disp('No time series here');
            end
            
        end
        
        
        
        function obj=bvp_plot(obj)
            plot_config.print=obj.print;
            plot_config.visible=obj.visible;
            
            %Set up the sign for the background temperature and salintiy
            if obj.dy_T_mean==-1
                T_mean_var=1+obj.dy_T_mean*obj.z_list;
                T_mean_sign='+1-z';
                dy_T_mean_sign='-1';
            elseif obj.dy_T_mean==1
                T_mean_var=obj.dy_T_mean*obj.z_list;
                T_mean_sign='+z';
                dy_T_mean_sign='+1';
            end
            
            if obj.dy_S_mean==-1
                S_mean_var=1+obj.dy_S_mean*obj.z_list;
                S_mean_sign='+1-z';
                dy_S_mean_sign='-1';
            elseif obj.dy_S_mean==1
                S_mean_var=obj.dy_S_mean*obj.z_list;
                S_mean_sign='+z';
                dy_S_mean_sign='+1';
            end
            
            data{1}.x=obj.T_0(:,1)+T_mean_var;
            data{1}.y=obj.z_list;
            if strcmp(obj.problem, 'IVP')
               data{2}.x=obj.T_0(:,end)+T_mean_var;
               data{2}.y=obj.z_list;
            end
            plot_config.fontsize=20;

            plot_config.label_list={1,['$\bar{T}_0',T_mean_sign,'$'], '$z$'};
            plot_config.print_size=[1,500,900];
            plot_config.name=[obj.h5_name(1:end-3),'_HB_','T_0.png'];
            plot_line(data,plot_config);
            clear data
            
            data{1}.x=obj.d_T_0(:,1)+obj.dy_T_mean;
            data{1}.y=obj.z_list;
            plot_config.fontsize=20;
            plot_config.label_list={1,['$\partial_z \bar{T}_0',dy_T_mean_sign,'$'], '$z$'};
            plot_config.print_size=[1,500,900];
            plot_config.name=[obj.h5_name(1:end-3),'_HB_','d_T_0.png'];
            plot_line(data,plot_config);

            data{1}.x=obj.T_0(:,1)+T_mean_var;
            data{1}.y=obj.z_list;
            plot_config.label_list={1,['$\bar{T}_0',T_mean_sign,'$'], '$z$'};
            plot_config.print_size=[1,500,900];
            plot_config.ylim_list=[1,0,0.01];
            plot_config.name=[obj.h5_name(1:end-3),'_HB_','T_0_local.png'];
            plot_line(data,plot_config);
            plot_config.ylim_list=0;

            data{1}.x=obj.T_0(:,1)+T_mean_var;
            data{1}.y=obj.z_list;
            plot_config.label_list={1,['$\bar{T}_0',T_mean_sign,'$'], '$z$'};
            plot_config.print_size=[1,500,900];
            plot_config.ylim_list=[1,0.49,0.51];
%             plot_config.xlim_list=[1,0.49,0.51];
            plot_config.name=[obj.h5_name(1:end-3),'_HB_','T_0_local_core.png'];
            plot_line(data,plot_config);
            plot_config.ylim_list=0; %plot_config.xlim_list=0;

            
            
            data{1}.x=obj.S_0(:,1)+S_mean_var;
            data{1}.y=obj.z_list;
            if strcmp(obj.problem, 'IVP')
               data{2}.x=obj.S_0(:,end)+S_mean_var;
               data{2}.y=obj.z_list;
            end
            plot_config.label_list={1,['$\bar{S}_0',S_mean_sign,'$'], '$z$'};
            plot_config.print_size=[1,500,900];
            plot_config.name=[obj.h5_name(1:end-3),'_HB_','S_0.png'];
            plot_line(data,plot_config);
            clear data

            
            data{1}.x=obj.d_S_0(:,1)+obj.dy_S_mean;
            data{1}.y=obj.z_list;
            plot_config.label_list={1,['$\partial_z \bar{S}_0',dy_S_mean_sign,'$'], '$z$'};
            plot_config.print_size=[1,500,900];
            plot_config.name=[obj.h5_name(1:end-3),'_HB_','d_S_0.png'];
            plot_line(data,plot_config);

            data{1}.x=obj.S_0(:,1)+S_mean_var;
            data{1}.y=obj.z_list;
            plot_config.label_list={1,['$\bar{S}_0',S_mean_sign,'$'], '$z$'};
            plot_config.print_size=[1,500,900];
            plot_config.ylim_list=[1,0,0.01];
            plot_config.name=[obj.h5_name(1:end-3),'_HB_','S_0_local.png'];
            plot_line(data,plot_config);
            plot_config.ylim_list=0;

            data{1}.x=obj.S_0(:,1)+S_mean_var;
            data{1}.y=obj.z_list;
            plot_config.label_list={1,['$\bar{S}_0',S_mean_sign,'$'], '$z$'};
            plot_config.print_size=[1,500,900];
            plot_config.ylim_list=[1,0.49,0.51];
%             plot_config.xlim_list=[1,0.49,0.51];
            plot_config.name=[obj.h5_name(1:end-3),'_HB_','S_0_local_core.png'];
            plot_line(data,plot_config);
            plot_config.ylim_list=0;% plot_config.xlim_list=0;


            data{1}.x=obj.T_hat(:,1);
            data{1}.y=obj.z_list;
            if strcmp(obj.problem, 'IVP')
               data{2}.x=obj.T_hat(:,end);
               data{2}.y=obj.z_list;
            end
            plot_config.label_list={1,'$\widehat{T}$', '$z$'};
            plot_config.print_size=[1,500,900];
            plot_config.name=[obj.h5_name(1:end-3),'_HB_','T_hat.png'];
            plot_line(data,plot_config);
            clear data

            
            data{1}.x=obj.S_hat(:,1);
            data{1}.y=obj.z_list;
            if strcmp(obj.problem, 'IVP')
               data{2}.x=obj.S_hat(:,end);
               data{2}.y=obj.z_list;
            end
            plot_config.label_list={1,'$\widehat{S}$', '$z$'};
            plot_config.print_size=[1,500,900];
            plot_config.name=[obj.h5_name(1:end-3),'_HB_','S_hat.png'];
            plot_line(data,plot_config);
            clear data


            data{1}.x=obj.S_hat(:,1);
            data{1}.y=obj.z_list;
            plot_config.label_list={1,'$\widehat{S}$', '$z$'};
            plot_config.print_size=[1,500,900];
            plot_config.ylim_list=[1,0,0.01];
            plot_config.name=[obj.h5_name(1:end-3),'_HB_','S_hat_local.png'];
            plot_line(data,plot_config);
            plot_config.ylim_list=0;

            data{1}.x=obj.w_hat(:,1);
            data{1}.y=obj.z_list;
            if strcmp(obj.problem, 'IVP')
               data{2}.x=obj.w_hat(:,end);
               data{2}.y=obj.z_list;
            end
            plot_config.label_list={1,'$\widehat{w}$', '$z$'};
            plot_config.print_size=[1,500,900];
            plot_config.name=[obj.h5_name(1:end-3),'_HB_','w_hat.png'];
            plot_line(data,plot_config);
            clear data

            data{1}.x=obj.u_tilde(:,1);
            data{1}.y=obj.z_list;
            if strcmp(obj.problem, 'IVP')
               data{2}.x=obj.u_tilde(:,end);
               data{2}.y=obj.z_list;
            end
            plot_config.label_list={1,'$\widetilde{u}$', '$z$'};
            plot_config.print_size=[1,500,900];
            plot_config.name=[obj.h5_name(1:end-3),'_HB_','u_tilde.png'];
            plot_line(data,plot_config);
            clear data

            %Update 2021/12/06, add streamline...
            data{1}.x=linspace(0,1,10);
            data{1}.y=linspace(0,1,10);
            data{1}.z=NaN*ones(10,10);
            x=linspace(0,2*pi,1000);
            %z_ind_N=length(obj.z_list)/10;
            %z_ind=round(linspace(1,length(obj.z_list),z_ind_N));
            z_ind=1:length(obj.z_list);
            y=obj.z_list(z_ind);
            [data{2}.x,data{2}.y]=meshgrid(x,y);
            %data{2}.y=obj.z_list;
            data{2}.u=obj.u_tilde(z_ind,1)*real(1i*exp(1i*x));
            data{2}.v=obj.w_hat(z_ind,1)*real(exp(1i*x));
            plot_config.xlim_list=[1,0,2*pi];
            plot_config.xtick_list=[1,0,pi/2,pi,3*pi/2,2*pi];
            plot_config.xticklabels_list={1,'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'};
            plot_config.ylim_list=[1,0,1];
            plot_config.label_list={1,'$x k_x$','$z$'};
            plot_config.streamline=1;
            plot_config.user_color_style_marker_list={'k-','b--'};
            plot_config.panel_num=2;
            plot_config.colorbar=0;
            plot_config.name=[obj.h5_name(1:end-3),'_HB_','streamline.png'];
            plot_contour(data,plot_config);
            
            if obj.no_ylabel
                plot_config.label_list{3}='';
                plot_config.name=[obj.h5_name(1:end-3),'_HB_','streamline_no_ylabel.png'];
                plot_contour(data,plot_config);
            end
            
            %Plot the isocontour of T
            clear data plot_config
            plot_config.print=obj.print;
            plot_config.visible=obj.visible;
            
            x=linspace(0,2*pi,1000);
            z_ind=1:length(obj.z_list);
            y=obj.z_list(z_ind);
            [data{1}.x,data{1}.y]=meshgrid(x,y);
            data{1}.z=T_mean_var+obj.T_0(z_ind,1)+obj.T_hat(z_ind,1)*2*cos(x);
            %data{2}.y=obj.z_list;
%             data{1}.u=obj.u_tilde(z_ind)*real(1i*exp(1i*x));
%             data{1}.v=obj.w_hat(z_ind)*real(exp(1i*x));
            plot_config.xlim_list=[1,0,2*pi];
            plot_config.xtick_list=[1,0,pi/2,pi,3*pi/2,2*pi];
            plot_config.xticklabels_list={1,'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'};
            plot_config.ylim_list=[1,0,1];
            plot_config.ytick_list=[1,0,0.2,0.4,0.6,0.8,1];
            plot_config.label_list={1,'$x k_x$','$z$'};
            plot_config.contour_line=0;
            plot_config.colorbar=1;
            %plot_config.user_color_style_marker_list={'k-','b--'};
            plot_config.name=[obj.h5_name(1:end-3),'_HB_','isocontour_T.png'];
            plot_config.print_size=[1,1000,900];
%             plot_config.colormap='bluewhitered';
            plot_config.zlim_list=[1,0,1];
            plot_config.ztick_list=[1,0,0.2,0.4,0.6,0.8,1];
            plot_contour(data,plot_config);
            if obj.no_ylabel
                plot_config.label_list{3}='';
                plot_config.name=[obj.h5_name(1:end-3),'_HB_','isocontour_T_no_ylabel.png'];
                plot_contour(data,plot_config);
            end
            
            clear data plot_config
            plot_config.print=obj.print;
            plot_config.visible=obj.visible;
            
            x=linspace(0,2*pi,1000);
            z_ind=1:length(obj.z_list);
            y=obj.z_list(z_ind);
            [data{1}.x,data{1}.y]=meshgrid(x,y);
            data{1}.z=S_mean_var+obj.S_0(z_ind,1)+obj.S_hat(z_ind,1)*2*cos(x);
            plot_config.xlim_list=[1,0,2*pi];
            plot_config.xtick_list=[1,0,pi/2,pi,3*pi/2,2*pi];
            plot_config.xticklabels_list={1,'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'};
            plot_config.ylim_list=[1,0,1];
            plot_config.ytick_list=[1,0,0.2,0.4,0.6,0.8,1];
            plot_config.label_list={1,'$x k_x$','$z$'};
            plot_config.contour_line=0;
            plot_config.colorbar=1;
%             plot_config.colormap='bluewhitered';
            plot_config.zlim_list=[1,0,1];
            plot_config.ztick_list=[1,0,0.2,0.4,0.6,0.8,1];
            %plot_config.user_color_style_marker_list={'k-','b--'};
            plot_config.name=[obj.h5_name(1:end-3),'_HB_','isocontour_S.png'];
            plot_config.print_size=[1,1000,900];
            plot_contour(data,plot_config);
            if obj.no_ylabel
                plot_config.label_list{3}='';
                plot_config.name=[obj.h5_name(1:end-3),'_HB_','isocontour_S_no_ylabel.png'];
                plot_contour(data,plot_config);
            end
            
%             plot_config.panel_num=4;
            
            %Plot the isocontour of S
        end
        
        
        function obj=dedalus_post_evp(obj)
            
            
        
        end
        
        function obj=evp_plot(obj)
            
        end
        
        
        function obj=snapshot(obj,variable_name,video_ind,zlim_list,snap_ind)
            %%plot the snapshot of salinity and generate video if any
            if nargin<3 || isempty(video_ind)
               video_ind=1; 
            end
            
            z_mesh=obj.z_list*ones(1,length(obj.x_list));
            if strcmp(variable_name,'S_tot')
                obj.(variable_name)=h5read_complex(obj.h5_name,['/tasks/','S'])...
                    +obj.dy_S_mean*z_mesh;
                
                %This is for fixed flux on S. 
                if obj.flux_S
                    for t_ind=1:length(obj.t_list)
                        obj.S=h5read_complex(obj.h5_name,['/tasks/','S']);
                        obj.(variable_name)(:,:,t_ind)=obj.S(:,:,t_ind)...
                            +mean(mean(obj.dy_S_mean_q(:,:,t_ind)))*z_mesh;
                    end
                end
                if obj.dy_S_mean==-1
                    obj.(variable_name)=obj.(variable_name)+1;
                end
                plot_config.colormap='jet';%bluewhitered
                plot_config.zlim_list=[1,0,1];
                plot_config.ztick_list=[1,0,0.2,0.4,0.6,0.8,1];
            elseif strcmp(variable_name,'T_tot')
                obj.(variable_name)=h5read_complex(obj.h5_name,['/tasks/','T'])...
                    +obj.dy_T_mean*z_mesh;
                
                if obj.flux_T
                    obj.T=h5read_complex(obj.h5_name,['/tasks/','T']);
                    for t_ind=1:length(obj.t_list)
                        obj.(variable_name)(:,:,t_ind)=obj.T(:,:,t_ind)...
                            +mean(mean(obj.dy_T_mean_q(:,:,t_ind)))*z_mesh;
                    end
                end
                
                if obj.dy_T_mean<0
                    obj.(variable_name)=obj.(variable_name)+1;
                end
                plot_config.colormap='jet';%bluewhitered
                plot_config.zlim_list=[1,0,1];
                plot_config.ztick_list=[1,0,0.2,0.4,0.6,0.8,1];
            else
                obj.(variable_name)=h5read_complex(obj.h5_name,['/tasks/',variable_name]);
                background=0*z_mesh;
                plot_config.colormap='bluewhitered';%bluewhitered
                variable_max=max(max(max(obj.(variable_name))));
                variable_min=min(min(min(obj.(variable_name))));
                plot_config.zlim_list=[1,variable_min,variable_max];
            end
            
            if strcmp(obj.z_basis_mode,'Fourier')
                plot_config.zlim_list=[1,min(min(min(obj.(variable_name)))),max(max(max(obj.(variable_name))))];
                plot_config.ztick_list=0;
            end
            
            %Update 2022/09/14, update the zlim_list option so I could turn
            %off this. 
            if nargin<4 || isempty(zlim_list)
                %do nothing
            else
               plot_config.zlim_list=zlim_list; 
            end
            
            if obj.video
                snapshot_ind=1;
                for t_ind=1:video_ind:length(obj.t_list)
                    data{1}.z=obj.(variable_name)(:,:,t_ind);
                    
                    data{1}.x=obj.x_list;
                    data{1}.y=obj.z_list;
                    plot_config.label_list={1,'$x$','$z$'};

                    plot_config.fontsize=40;
                    plot_config.ylim_list=[1,round(min(data{1}.y)),round(max(data{1}.y))];
                    plot_config.ytick_list=[1,0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2];
                    
                    if obj.title_time
                        plot_config.title_list={1,['$t=$',num2str(round(obj.t_list(t_ind),2))]};
                    else
                        plot_config.title_list={0};
                    end
                    
                    plot_config.print_size=[1,1200,1200];
                    plot_config.name=[obj.h5_name(1:end-3),'_snapshot_',variable_name,'_t_',num2str(round(obj.t_list(t_ind))),'.png'];
                    plot_config.print=obj.print;
                    plot_config.visible=obj.visible;
                    snapshot(snapshot_ind)=plot_contour(data,plot_config);
                    snapshot_ind=snapshot_ind+1;
                    %plot_config.label_list={1,'$x$',''};
                    %plot_config.name=[obj.h5_name(1:end-3),'_snapshot_',variable_name,'_t_',num2str(round(obj.t_list(t_ind))),'_no_ylabel.png'];
                    %plot_contour(data,plot_config);
                end
               plot_config.name=[obj.h5_name(1:end-3),'_snapshot_',variable_name,'_t_video.mp4'];
               plot_video(snapshot,plot_config);
            end
            
            if nargin<5 || isempty(snap_ind)
                %do nothing
            else
                for t_ind=snap_ind
                    data{1}.z=obj.(variable_name)(:,:,t_ind);

                    data{1}.x=obj.x_list;
                    data{1}.y=obj.z_list;
                    plot_config.label_list={1,'$x$','$z$'};

%                     plot_config.fontsize=40;
                    plot_config.ylim_list=[1,round(min(data{1}.y)),round(max(data{1}.y))];
                    plot_config.ytick_list=[1,0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2.0];

                    if obj.title_time
                        plot_config.title_list={1,['$t=$',num2str(round(obj.t_list(t_ind),2))]};
                    else
                        plot_config.title_list={0};
                    end

                    plot_config.print_size=[1,1600,900];
                    plot_config.name=[obj.h5_name(1:end-3),'_snapshot_',variable_name,'_t_ind_',num2str(t_ind),'.png'];
                    plot_config.print=obj.print;
                    plot_config.visible=obj.visible;
                    plot_config.fontsize=38;
                    plot_contour(data,plot_config);
                    %plot_config.label_list={1,'$x$',''};
                    %plot_config.name=[obj.h5_name(1:end-3),'_snapshot_',variable_name,'_t_',num2str(t_ind),'_no_ylabel.png'];
                    %plot_contour(data,plot_config);
                end
            end
        end
        
        
        function obj=spectrum_t(obj,variable_name,z,x,t_range)
            %z is the vertical position, x is the first horizontal... 
            %direction
            %right now assume it as the 2D results... 
            
            if nargin<2 || isempty(variable_name)
               variable_name='S'; 
            end
            
            if nargin<3 || isempty(z)
                z=obj.Lz/2; %if not provided, just plot the middle plan x-t contour
            end
            if nargin<4 || isempty(x)
               x=obj.Lx/2;
            end
            if nargin<5 || isempty(t_range)
               t_range(1)=obj.t_list(1);
               t_range(2)=obj.t_list(end);
            elseif length(t_range)<2
                t_range(2)=obj.t_list(end);
            end
            
            if strcmp(x,'ave')
                x_ind_list=1:length(obj.x_list);
            else
                [val,x_ind_list]=min(abs(obj.x_list-x));
            end
            [val,z_ind]=min(abs(obj.z_list-z));
            [val,t_ind_begin]=min(abs(obj.t_list-t_range(1)));
            [val,t_ind_end]=min(abs(obj.t_list-t_range(2)));
            t_list_dedalus=obj.t_list(t_ind_begin:t_ind_end);
            dt=mean(diff(obj.t_list(t_ind_begin:t_ind_end)));
            t_list_uniform=obj.t_list(t_ind_begin):dt:obj.t_list(t_ind_end);
            
            %             [val,z_ind]=min(abs(z-obj.z_list));
            
%             data{1}.x=obj.t_list;
%             if strcmp(obj.flow(1:7),'IFSC_2D')
%                 data{1}.y=obj.z_list/(2*pi/obj.k_opt);
%                 plot_config.label_list={1,'$t$','$z/l_{opt}$'};
%             else
%                 data{1}.y=obj.x_list;
%                 plot_config.label_list={1,'$t$','$x$'};
%             end
            spec_t=0;
            period_t=0;
            
            for x_ind=x_ind_list
                switch variable_name
                    case {'u','v','w','S','T','p'}
                        obj.(variable_name)=h5read_complex(obj.h5_name,['/tasks/',variable_name]);
                        variable=squeeze(obj.(variable_name)(z_ind,x_ind,t_ind_begin:t_ind_end));
                        %plot_config.colormap='bluewhitered';
                    case {'u_x_ave'}
                        obj.(variable_name(1))=h5read_complex(obj.h5_name,['/tasks/',variable_name(1)]);
                        variable=squeeze(mean(obj.(variable_name(1))(z_ind,:,t_ind_begin:t_ind_end),2));
                    case 'U_0'    
                        obj.(variable_name)=h5read_complex(obj.h5_name,['/tasks/',variable_name]);
                        variable=obj.(variable_name)(z_ind,t_ind_begin:t_ind_end);
                    case {'S_tot','T_tot'}
                        obj.(variable_name)=h5read_complex(obj.h5_name,['/tasks/',variable_name(1)]);
                        variable=squeeze(obj.(variable_name)(z_ind,x_ind,t_ind_begin:t_ind_end))+obj.(['dy_',variable_name(1),'_mean'])*obj.z_list(z_ind);
                        %plot_config.colormap='jet';
                    case {'Nu_T_t'}
                        variable=obj.(variable_name)(t_ind_begin:t_ind_end);
                end
                variable_uniform=real(variable);
                
                ind_local_min=find(islocalmin(variable_uniform));
                if ~isempty(ind_local_min) && length(ind_local_min)>=2
                    variable_uniform=variable_uniform(ind_local_min(1):ind_local_min(end));
                    t_list_uniform=t_list_uniform(ind_local_min(1):ind_local_min(end));
                end
                Nt=length(t_list_uniform);
                Fs=1/dt;
                freq=Fs*(0:(Nt/2))/Nt;
                %variable_uniform=interp1(t_list_dedalus,variable,t_list_uniform,'linear');
                %variable_uniform=sin(2*pi*1/10*Fs*t_list_uniform); %This is to
                %test fft results using sinusoidal functino
                spec_tmp=abs(fft(variable_uniform)/Nt);
                spec_tmp=spec_tmp(1:Nt/2+1);
                spec_tmp(2:end-1)=2*spec_tmp(2:end-1);
                spec_t=spec_t+spec_tmp;
                
                ind_local_max=islocalmax(real(variable_uniform));
                period_t=period_t+mean(diff(t_list_uniform(ind_local_max)));
            end
            spec_t=spec_t/length(x_ind_list);
            period_t=period_t/length(x_ind_list);
            data{1}.x=freq;
            data{1}.y=spec_t;
            plot_config.label_list={1,'$f $','PSD'};
            
            if strcmp(x,'ave')
                plot_config.name=[obj.h5_name(1:end-3),'_',variable_name,'_spectrum_t_at_z=',num2str(round(z,2)),'_x=',x,'.png'];
            else
                plot_config.name=[obj.h5_name(1:end-3),'_',variable_name,'_spectrum_t_at_z=',num2str(round(z,2)),'_x=',num2str(round(x,2)),'.png'];
            end
            %plot_config.xlim_list=[1,0,10];
            plot_config.print=obj.print;
            plot_config.visible=obj.visible;
            plot_config.ylim_list=[1,0,1.1*max(spec_t(2:end))];
            plot_config.xlim_list=[1,0,0.4*max(freq)];
            plot_line(data,plot_config);
            
            ind_local_max=islocalmax(data{1}.y);
            obj.freq_local_max=data{1}.x(ind_local_max);
            
            data{1}.x=t_list_uniform;
            data{1}.y=variable_uniform;
            plot_config.label_list={1,'$t $','Var'};

            if strcmp(x,'ave')
                plot_config.name=[obj.h5_name(1:end-3),'_',variable_name,'_spectrum_t_history_at_z=',num2str(round(z,2)),'_x=',x,'.png'];
            else
                plot_config.name=[obj.h5_name(1:end-3),'_',variable_name,'_spectrum_t_history_at_z=',num2str(round(z,2)),'_x=',num2str(round(x,2)),'.png'];
            end
            plot_config.print=obj.print;
            plot_config.visible=obj.visible;
            plot_config.xlim_list=0;
            plot_config.ylim_list=0;
            if strcmp(variable_name,'Nu_T_t')
                plot_config.label_list={1,'$t$','$nu(t)$'};
            end
            plot_line(data,plot_config);
            
            obj.freq=freq;
            obj.spec_t=spec_t;
            [spec_t_sort,ind]=sort(spec_t,'descend');
            obj.freq_sort=obj.freq(ind); %the corresponding frequency in descending way... 
        
            obj.period_t=period_t;
        end
        
        
        
        function obj=spectrum_snapshot(obj,variable_name)
            %%plot the spectrum of salnity as time varies, also generate
            %%video if any
            obj.([variable_name,'_coeff'])=h5read_complex(obj.h5_name,['/tasks/',variable_name,'_coeff']);
            %coeff.r+1i*coeff.i;
            if obj.video
                for t_ind=1:length(obj.t_list)
                    clear data plot_config;
                    
                    data{1}.x=obj.kx_list;
                    data{1}.y=obj.kz_list(1:obj.Nz/2);
                    plot_config.label_list={1,'$k_x$','$k_z$'};
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
                    
                    data{1}.x=obj.kx_list;
                    data{2}.x=obj.kz_list(1:obj.Nz/2);
                    plot_config.label_list={1,'$k_x$ or $k_z$',''};
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
               plot_config.name=[obj.h5_name(1:end-3),'_spectrum_',variable_name,'_2D_t_video.mp4'];
               plot_video(frame_spectrum_2D,plot_config);
               plot_config.name=[obj.h5_name(1:end-3),'_spectrum_',variable_name,'_1D_t_video.mp4'];
               plot_video(frame_spectrum_1D,plot_config);
           end
        end
        
        function obj=spectrum_average(obj,variable_name,t_range)
            %%This function plot the 
            %%plot the overall spectrum averaged over time
            if nargin<3 || isempty(t_range)
               t_range(1)=obj.t_list(1);
               t_range(2)=obj.t_list(end); 
            end
           
            [val,t_ind_begin]=min(abs(obj.t_list-t_range(1)));
            [val,t_ind_end]=min(abs(obj.t_list-t_range(2)));
            
            %coeff.r+1i*coeff.i;
            
%             obj.(variable_name)=h5read_complex(obj.h5_name,['/tasks/',variable_name]);
%             for t_ind=1:length(obj.t_list)
%                 obj.(['E_',variable_name])(t_ind)=sum(sum(obj.(variable_name)(:,:,t_ind).^2))/obj.Nx/obj.Nz/2;
%             end
%             [val,max_ind]=max(obj.(['E_',variable_name]));
%            
            if strcmp(obj.z_basis_mode,'Chebyshev')
                %modify for the vertical as the bounded domain.. Chebyshev
                %basis instead of the Fourier.. so we do not have spectrum
                %int he vertical
                obj.([variable_name])=h5read_complex(obj.h5_name,['/tasks/',variable_name]);

                x_list=obj.x_list;
                dx=mean(diff(obj.x_list));
                Nx=length(x_list);
                Fs=1/dx;
                obj.kx_list=2*pi*Fs*(0:(Nx/2-1))/Nx; %Note that this wavenumber corresponding to circular frequency
                t_ave_num=0;
                spec_kx_z=zeros(length(obj.kx_list),length(obj.z_list));
                for t_ind=t_ind_begin:t_ind_end
                    for z_ind=1:length(obj.z_list)
                        spec_tmp=abs(fft(obj.([variable_name])(z_ind,:,t_ind))/Nx);
                        spec_tmp=spec_tmp(1:Nx/2);
                        %spec_tmp(2:end)=2*spec_tmp(2:end);
                        
                        spec_kx_z(:,z_ind)=spec_kx_z(:,z_ind)+(spec_tmp)';
                    end
                    t_ave_num=t_ave_num+1;
                end
                
                spec_kx_z=spec_kx_z/t_ave_num;
%                 spec_kx_z=abs(spec_kx_z);
                obj.(['spec_kx_z_',variable_name])=spec_kx_z;
                data{1}.x=obj.kx_list;
                data{1}.y=obj.z_list;
                data{1}.z=(spec_kx_z)';
                plot_config.label_list={1,'$k_x$','$z$'};

                %spectrum_average=mean(abs(obj.([variable_name,'_coeff'])(:,:,t_ind_begin:t_ind_end)).^2,3);
                %data{1}.z=log10(spectrum_average);
                plot_config.loglog=[0,0];
                plot_config.print_size=[1,800,900];
                plot_config.colormap='white_zero';
                plot_config.name=[obj.h5_name(1:end-3),'_spectrum_',variable_name,'_2D_time_average.png'];
                plot_config.print=obj.print;
                plot_config.ytick_list=[1,0,0.2,0.4,0.6,0.8,1];
                plot_config.xtick_list=[1,0,10,20,30];
                plot_config.xlim_list=[1,0,20];
                plot_config.ylim_list=[1,0,1];
                
                %This range is for the plotting of the bounded salt-finger
                plot_config.zlim_list=[1,0,0.16];
                plot_config.ztick_list=[1,0,0.04,0.08,0.12,0.16]
                
                plot_config.zlim_list=0;
                plot_config.ztick_list=0;
                plot_config.xlim_list=[1,0,50];
                
                plot_config.visible=1;%obj.visible;
                plot_config.fontsize=32;
                plot_contour(data,plot_config);

            elseif strcmp(obj.z_basis_mode,'Fourier')
                try 
                    obj.([variable_name,'_coeff'])=h5read_complex(obj.h5_name,['/tasks/',variable_name,'_coeff']);
                catch
                    obj.([variable_name])=h5read_complex(obj.h5_name,['/tasks/',variable_name]);
                    for t_ind=1:length(obj.t_list)
                        spectrum_tmp=fft2(obj.(variable_name)(:,:,t_ind));
                    end
                    obj.([variable_name,'_coeff'])=spectrum_tmp(1:obj.Nz-1,1:obj.Nx/2)/obj.Nx/obj.Nz;
                end
                data{1}.x=obj.kx_list;
                data{1}.y=obj.kz_list(1:obj.Nz/2);
                plot_config.label_list={1,'$k_x$','$k_z$'};

                spectrum_average=mean(abs(obj.([variable_name,'_coeff'])(1:obj.Nz/2,:,1:end)).^2,3);
                data{1}.z=log10(spectrum_average);
                plot_config.loglog=[0,0];
                plot_config.print_size=[1,1100,900];
                plot_config.colormap='white_zero';
                plot_config.name=[obj.h5_name(1:end-3),'_spectrum_',variable_name,'_2D_time_average.png'];
                plot_config.print=obj.print;
                plot_config.ytick_list=[0,10^(-8),10^(-7),10^(-6),10^(-5),10^(-4),...
                    0.001,0.01,0.1,1,10,100,1000];
                plot_config.xtick_list=[0,10^(-8),10^(-7),10^(-6),10^(-5),10^(-4),...
                    0.001,0.01,0.1,1,10,100,1000];
                plot_config.visible=obj.visible;
                plot_contour(data,plot_config);

                dx=diff(obj.kx_list); dx=dx(1);
                dz=diff(obj.kz_list); dz=dz(1);

                data{1}.x=obj.kx_list;
                data{1}.y=2*dz*sum(spectrum_average,1);
                data{2}.x=obj.kz_list(1:obj.Nz/2);
                data{2}.y=2*dx*sum(spectrum_average,2);

                plot_config.loglog=[1,1];

                plot_config.label_list={1,'$k_x$ or $k_z$',''};
                plot_config.legend_list={1,['$\int E_',variable_name,'(k_x,k_z)dk_z$'],['$\int E_',variable_name,'(k_x,k_z)d k_x$']};
                plot_config.name=[obj.h5_name(1:end-3),'_spectrum_',variable_name,'_1D_time_average.png'];
                plot_config.print=obj.print;
                plot_config.visible=obj.visible;
                plot_line(data,plot_config);
            end
        end
        
        
        function obj=spectrum_TKE_average(obj)
            %%This function plot the kinetic energy, (deduct the background laminar flow)
            %%This is average over time... plot the turbulence kinetic
            %%energy spectrum as kx and kz
            
            %%This is the post-processing for the TKE in 2D... 
            obj.w_coeff=h5read_complex(obj.h5_name,'/tasks/w_coeff');
            %=w_coeff.r+1i*w_coeff.i;
            obj.u_coeff=h5read_complex(obj.h5_name,'/tasks/u_coeff');
            %obj.u_coeff=u_coeff.r+1i*u_coeff.i;
            for t_ind=1:length(obj.t_list)
                obj.spectrum_TKE(:,:,t_ind)=abs(obj.u_coeff(:,:,t_ind)).^2+abs(obj.w_coeff(:,:,t_ind)).^2;
            end

%             obj.u=h5read_complex(obj.h5_name,'/tasks/u');
%             obj.w=h5read_complex(obj.h5_name,'/tasks/w');
% 
%             for t_ind=1:length(obj.t_list)
%                 obj.TKE_time(t_ind)=sum(sum(obj.u(:,:,t_ind).^2+obj.w(:,:,t_ind).^2))/obj.Nx/obj.Nz/2;
%             end
% 
%             [val,max_ind]=max(obj.TKE_time);
            spectrum_TKE_average=mean(abs(obj.spectrum_TKE(1:obj.Nz/2,:,1:end)).^2,3);%max(3*max_ind,30)

            data{1}.z=log10(spectrum_TKE_average);
            data{1}.x=obj.kx_list;
            data{1}.y=obj.kz_list(1:obj.Nz/2);
            
            plot_config.loglog=[0,0];
            plot_config.print_size=[1,1100,900];
            plot_config.label_list={1,'$k_x$','$k_z$'};
            plot_config.colormap='white_zero';
            plot_config.name=[obj.h5_name(1:end-3),'_spectrum_TKE_2D_time_average.png'];
            plot_config.print=obj.print;
            plot_config.visible=obj.visible;
            plot_contour(data,plot_config);
            
            dx=diff(obj.kx_list); dx=dx(1);
            dz=diff(obj.kz_list); dz=dz(1);

            data{1}.x=obj.kx_list;
            data{1}.y=2*dz*sum(spectrum_TKE_average,1);
            data{2}.x=obj.kz_list(1:obj.Nz/2);
            data{2}.y=2*dx*sum(spectrum_TKE_average,2);
            
            plot_config.loglog=[1,1];
            plot_config.ytick_list=[0,0.001,0.01,0.1,1,10,100,1000];
            plot_config.label_list={1,'$k_x$ or $k_z$',''};
            plot_config.legend_list={1,'$\int E_u(k_x,k_z)dk_z$','$\int E_u(k_x,k_z)d k_x$'};
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
            obj.(variable_name)=h5read_complex(obj.h5_name,['/tasks/',variable_name]);
            
            for t_ind=1:length(obj.t_list)
                if strcmp(obj.z_basis_mode,'Fourier')
                    obj.(['E_',variable_name])(t_ind)=sum(sum(obj.(variable_name)(:,:,t_ind).^2))/obj.Nx/obj.Nz/2;
                else
                    obj.(['E_',variable_name])(t_ind)=sum(sum(obj.(variable_name)(:,:,t_ind).^2))/obj.Nx/obj.Nz/2;
                    warning('Note that energy is assuming uniform grid!!!');
                end
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
            if obj.flux_T
               derivative_factor=min(2*pi./[obj.Lx,obj.Lz]);
               data{2}.x=obj.t_list;
               obj.E_T_upper_bound=1/8/derivative_factor^2;
               data{2}.y=obj.E_T_upper_bound*ones(size(obj.t_list));
               plot_config.legend_list={1,'DNS','Upper bound'};
               plot_config.user_color_style_marker_list={'k-','b--'};

            end
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
               plot_config.name=[obj.h5_name(1:end-3),'_E_',variable_name,'_t_video.mp4'];
               plot_video(E_time,plot_config);
            end
            
        end
                
        function obj=E_TKE_time(obj)
            %%Plot the turbulence kinetic energy as a function over time
            obj=obj.u_fluctuation_read;
            obj.w=h5read_complex(obj.h5_name,'/tasks/w');
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
  
            data{1}.y=obj.z_list;
            data{2}.y=obj.z_list;
            data{3}.y=obj.z_list;
            data{4}.y=obj.z_list;
            plot_config.label_list={1,'','$z$'};
            
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
            plot_config.ylim_list=[1,round(min(data{1}.y)),round(max(data{1}.y))];
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
            obj.u=h5read_complex(obj.h5_name,'/tasks/u');
            obj.u_fluctuation=obj.u;
            for z_ind=1:length(u_laminar_num)
                obj.u_fluctuation(z_ind,:,:)=obj.u(z_ind,:,:)-u_laminar_num(z_ind);
            end
        end
        
        function obj=u_fluctuation_x_ave(obj)
            %%plot the streamwise averaged u fluctuations...
            %%as a function of z (vertical axis) and time
            obj=obj.u_fluctuation_read();

            data{1}.x=obj.t_list;
%             if strcmp(obj.flow(1:7),'IFSC_2D')
%                 data{1}.y=obj.z_list/(2*pi/obj.k_opt);
%                 plot_config.label_list={1,'$t$','$z/l_{opt}$'};
%             else
                data{1}.y=obj.z_list;
                plot_config.label_list={1,'$t$','$z$'};
%             end
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

        
        
        function obj=x_ave(obj,variable_name,t_ind)
            %%plot the streamwise averaged salnity
            %%as a function of z (vertical axis) and time
            if nargin<3 || isempty(t_ind)
                %The default option, just average over the second half of
                %data...
                t_ind_begin=1;
                t_ind_end=length(obj.t_list);
            else
                t_ind_begin=t_ind(1);
                if length(t_ind)==1
                    t_ind_end=length(obj.t_list);
                else
                    t_ind_end=t_ind(2);
                end
            end
            
            data{1}.x=obj.t_list(t_ind_begin:t_ind_end);
%             if strcmp(obj.flow(1:7),'IFSC_2D')
%                 data{1}.y=obj.z_list/(2*pi/obj.k_opt);
%                 plot_config.label_list={1,'$t$','$z/l_{opt}$'};
%             else
            data{1}.y=obj.z_list;
            plot_config.label_list={1,'$t$','$z$'};
%             end
            if strcmp(obj.flow,'HB_benard_shear_periodic')
                switch variable_name
                    case {'u','w','T'}
                        obj.(variable_name)=h5read_complex(obj.h5_name,['/tasks/',upper(variable_name),'_0']);
                        data{1}.z=real(obj.(variable_name)(:,t_ind_begin:t_ind_end));
                    case 'dy_T_mean_q'
                        obj.(variable_name)=h5read_complex(obj.h5_name,['/tasks/',variable_name]);
                        data{1}.z=obj.(variable_name)(:,t_ind_begin:t_ind_end);
                end
            else
                switch variable_name
                    case {'u','v','w','S','T','p','dy_T_mean_q','dy_S_mean_q'}
                        obj.(variable_name)=h5read_complex(obj.h5_name,['/tasks/',variable_name]);
                    case {'uS','wS','uT','wT','uw','ww'}%%
                        var_1=variable_name(1);
                        var_2=variable_name(2);
                        if strcmp(var_1,'u') %%this require in default, the u is always in the first variable....
                           obj=obj.u_fluctuation_read();
                           var_1_data=obj.u_fluctuation;
                        else
                           var_1_data=h5read_complex(obj.h5_name,['/tasks/',var_1]);
                        end
                        var_2_data=h5read_complex(obj.h5_name,['/tasks/',var_2]);
                        obj.(variable_name)=var_1_data.*var_2_data;
                end

                data{1}.z=squeeze(mean(obj.(variable_name)(:,:,t_ind_begin:t_ind_end),2));

            end
            %             plot_config.label_list={1,'$t$','$z/l_{opt}$'};
            plot_config.colormap='bluewhitered';
            plot_config.print_size=[1,1600,900];
%             plot_config.ztick_list=[1,-0.001,0.001];
            plot_config.print=obj.print;
            plot_config.visible=obj.visible;
            plot_config.name=[obj.h5_name(1:end-3),'_',variable_name,'_x_ave.png'];
            plot_config.ylim_list=[1,round(min(data{1}.y),1),round(max(data{1}.y),1)];
            if round(min(data{1}.x),1)==round(max(data{1}.x),1)
                plot_config.xlim_list=[1,min(data{1}.x),max(data{1}.x)];
            else
                plot_config.xlim_list=[1,round(min(data{1}.x),1),round(max(data{1}.x),1)];
            end
%             plot_config.ytick_list=[1,0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2];
            %plot_config.fontsize=28;
            %plot_config.xlim_list=[1,150,300];
            %plot_config.xtick_list=[1,5,10,15];
            plot_contour(data,plot_config);
            
            if strcmp(variable_name,'u')
                for t_ind=1:length(data{1}.x)
                    [val,ind]=mink(abs(data{1}.z(:,t_ind)),5);
                    
                    Nz=obj.Nz;
                    if abs(mod(ind(2),Nz)-mod(ind(1),Nz))>1
                        z_max(t_ind)=(data{1}.y(ind(1))+data{1}.y(ind(2)))/2;
                    elseif abs(mod(ind(3),Nz)-mod(ind(1),Nz))>2
                        z_max(t_ind)=(data{1}.y(ind(1))+data{1}.y(ind(3)))/2;
                    else
                        z_max(t_ind)=(data{1}.y(ind(1))+data{1}.y(ind(4)))/2;
                    end
                end
                
                data_z{1}.x=data{1}.x;
                data_z{1}.y=z_max';
                plot_config.label_list={1,'$t$','$z$'};
                plot_config.name=[obj.h5_name(1:end-3),'_',variable_name,'_TW_c_z.png'];
                plot_config.ylim_list=0;
                plot_line(data_z,plot_config);
                
                for t_ind=1:length(z_max)-5
                    if abs(z_max(t_ind+1)-z_max(t_ind))>0.2
                         if all(abs(z_max(t_ind+2:t_ind+5)-z_max(t_ind+1))<0.2)
                             z_max(t_ind+1:end)=z_max(t_ind+1:end)+0.5;
                         else
                             z_max(t_ind+1)=z_max(t_ind+1)+0.5;
                         end
%                         if abs(z_max(t_ind+2)-z_max(t_ind+1))>0.45
%                             z_max(t_ind+1)=z_max(t_ind+1)+0.5;
%                         else
%                             z_max(t_ind+1:end)=z_max(t_ind+1:end)+0.5;
%                         end
                    end
                end
                data_z{1}.y=z_max';
                plot_config.label_list={1,'$t$','$z$'};
                plot_config.name=[obj.h5_name(1:end-3),'_',variable_name,'_TW_c_z_processed.png'];
                plot_config.ytick_list=0;
                plot_line(data_z,plot_config);
                
                coeff=[ones(length(data_z{1}.x),1),data_z{1}.x]\z_max';
                obj.phase_c_z=coeff(2);
                
            end
            
        end
        
        
        
        function obj=phase_diagram(obj,variable_name1,variable_name2,t_ind,option)
            if nargin<2 || isempty(variable_name1)
               variable_name1='u'; 
            end
            
            if nargin<3 || isempty(variable_name2)
               variable_name2='w'; 
            end
            
            if nargin<4 || isempty(t_ind)
                %The default option, just average over the second half of
                %data...
                t_ind_begin=1;
                t_ind_end=length(obj.t_list);
            else
                t_ind_begin=t_ind(1);
                if length(t_ind)==1
                    t_ind_end=length(obj.t_list);
                else
                    t_ind_end=t_ind(2);
                end
            end
            
            if nargin<5 || isempty(option)
                option='max_z';
            end
            %large-scale shear
            obj.(variable_name1)=h5read_complex(obj.h5_name,['/tasks/',variable_name1]);
                
            %mean temperature
            obj.(variable_name2)=h5read_complex(obj.h5_name,['/tasks/',variable_name2]);
            if strcmp(option(1:5),'max_z')
                if strcmp(option(6),'2')
                    x_ave=squeeze(mean(obj.(variable_name2)(:,:,t_ind_begin:t_ind_end),2));
                elseif strcmp(option(6),'1')
                    x_ave=squeeze(mean(obj.(variable_name1)(:,:,t_ind_begin:t_ind_end),2));
                end
                [val,z_ind]=max(x_ave(:,end));
                obj.z_phase_diagram=obj.z_list(z_ind);
                for t_ind=1:length(obj.t_list(t_ind_begin:t_ind_end))
                    data{1}.x(t_ind)=squeeze(mean(obj.(variable_name1)(z_ind,:,t_ind)));
                    data{1}.y(t_ind)=squeeze(mean(obj.(variable_name2)(z_ind,:,t_ind)));
                end
            end
            plot_config.name=[obj.h5_name(1:end-3),'_phase_diagram',variable_name1,'_',variable_name2,'_',option,'.png'];
            plot_config.label_list={1,['$\langle ',variable_name1,'\rangle_h(z_p,t)$'],['$\langle ',variable_name2,'\rangle_h(z_p,t)$']};
            plot_config.Markerindex=3;
            plot_config.user_color_style_marker_list={'k*'};
            plot_line(data,plot_config);
        end
        
        function obj=x_ave_z_max(obj,variable_name,ylim_list,t_ind)
            %%plot the streamwise averaged salnity
            %%as a function of z (vertical axis) and time
            if nargin<4 || isempty(t_ind)
                %The default option, just average over the second half of
                %data...
                t_ind_begin=1;
                t_ind_end=length(obj.t_list);
            else
                t_ind_begin=t_ind(1);
                if length(t_ind)==1
                    t_ind_end=length(obj.t_list);
                else
                    t_ind_end=t_ind(2);
                end
            end
            data{1}.x=obj.t_list(t_ind_begin:t_ind_end);
            switch variable_name
                case {'u','v','w','S','T','p','dy_T_mean_q','dy_S_mean_q'}
                    obj.(variable_name)=h5read_complex(obj.h5_name,['/tasks/',variable_name]);
                case {'uS','wS','uT','wT','uw','ww'}%%
                    var_1=variable_name(1);
                    var_2=variable_name(2);
                    if strcmp(var_1,'u') %%this require in default, the u is always in the first variable....
                       obj=obj.u_fluctuation_read();
                       var_1_data=obj.u_fluctuation;
                    else
                       var_1_data=h5read_complex(obj.h5_name,['/tasks/',var_1]);
                    end
                    var_2_data=h5read_complex(obj.h5_name,['/tasks/',var_2]);
                    obj.(variable_name)=var_1_data.*var_2_data;
            end
            data_tmp=squeeze(mean(obj.(variable_name)(:,:,t_ind_begin:t_ind_end),2));
            data{1}.y=max(abs(data_tmp));

%             plot_config.label_list={1,'$t$',['$\underset{z}{max}|\langle ',variable_name,'\rangle_{h}(z,t)|$']};
            plot_config.label_list={1,'$t$',''};
%             plot_config.label_list={1,'$t$','$z/l_{opt}$'};
            plot_config.print_size=[1,1600,1600];
%             plot_config.ztick_list=[1,-0.001,0.001];
            plot_config.print=obj.print;
            plot_config.name=[obj.h5_name(1:end-3),'_',variable_name,'_x_ave_z_max.png'];
            
            if nargin<3 || isempty(ylim_list)
                plot_config.ylim_list=[1,round(min(data{1}.y),1),round(max(data{1}.y),1)];
            else
                plot_config.ylim_list=[1,ylim_list];
            end
            plot_config.xlim_list=[1,round(min(data{1}.x),1),round(max(data{1}.x),1)]
            %plot_config.fontsize=28;
            %plot_config.ytick_list=[1,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2];
            plot_line(data,plot_config);
            
            
        end
        
        
        function obj=z_slice(obj,variable_name,z_list,t_inter)
            if nargin<2 || isempty(variable_name)
               variable_name='S_tot'; 
            end
            if nargin<3 || isempty(z_list)
                z_list=obj.Lz/2; %if not provided, just plot the middle plan x-t contour
            end
            if nargin<4 || isempty(t_inter)
                t_begin=1; 
                t_end=length(obj.t_list);
            else
                t_begin=t_inter(1);
                t_end=t_inter(2);
            end
            z_list_string=[];
            for z_plot_ind=1:length(z_list)
                z=z_list(z_plot_ind);
                [val,z_ind]=min(abs(z-obj.z_list));

                data{1}.x=obj.t_list(t_begin:t_end);
    %             if strcmp(obj.flow(1:7),'IFSC_2D')
    %                 data{1}.y=obj.z_list/(2*pi/obj.k_opt);
    %                 plot_config.label_list={1,'$t$','$z/l_{opt}$'};
    %             else
                    data{1}.y=obj.x_list;
                    plot_config.label_list={1,'$t$','$x$'};
    %             end
                switch variable_name
                    case {'u','v','w','S','T','p'}
                        obj.(variable_name)=h5read_complex(obj.h5_name,['/tasks/',variable_name]);
                        data{1}.z=squeeze(obj.(variable_name)(z_ind,:,t_begin:t_end));
                        plot_config.colormap='bluewhitered';
                        label=['$',variable_name,'$'];
                    case {'S_tot','T_tot'}
                        obj.(variable_name)=h5read_complex(obj.h5_name,['/tasks/',variable_name(1)]);
                        data{1}.z=squeeze(obj.(variable_name)(z_ind,:,t_begin:t_end))+obj.(['dy_',variable_name(1),'_mean'])*obj.z_list(z_ind);
                        if obj.flux_T
                            %old version use the instanteous mean
                            %temperature gradient
                            %data{1}.z=squeeze(obj.(variable_name)(z_ind,:,t_begin:t_end))+ones(obj.Nx,1)*transpose(squeeze(obj.(['dy_',variable_name(1),'_mean_q'])(1,:,t_begin:t_end)))*obj.z_list(z_ind);
                            
                            %new version that use time-averaged mean
                            %temperature gradient
                            data{1}.z=squeeze(obj.(variable_name)(z_ind,:,t_begin:t_end))+ones(obj.Nx,1)*ones(1,length(obj.t_list(t_begin:t_end)))*mean(obj.(['dy_',variable_name(1),'_mean_q'])(1,:,t_begin:t_end))*obj.z_list(z_ind);

                        end
                        if obj.(['dy_',variable_name(1),'_mean'])<0
                           data{1}.z=data{1}.z+1; 
                        end
                        plot_config.colormap='jet';
                        if obj.(['dy_',variable_name(1),'_mean'])==1
                            label=['$z+',variable_name(1),'$'];
                        elseif obj.(['dy_',variable_name(1),'_mean'])==-1
                            label=['$1-z+',variable_name(1),'$'];
                        end
                        
                        if obj.flux_T
                           label=['$1+\bar{\mathcal{T}}_{z,q}z+',variable_name(1),'$'];
                        end
                end
    %             data{1}.z=squeeze(mean(obj.(variable_name),2));
    %             plot_config.label_list={1,'$t$','$z/l_{opt}$'};
    
                plot_config.print_size=[1,1600,900];
                plot_config.print=obj.print;
                plot_config.name=[obj.h5_name(1:end-3),'_',variable_name,'_z_slice_at_z=',num2str(z),'.png'];
                plot_contour(data,plot_config);
                
                %also plot the version that put time in the vertical axis.
                data{1}.x=obj.x_list;
                data{1}.y=obj.t_list;
                data{1}.z=data{1}.z';
                plot_config.label_list={1,'$x$','$t$'};
                plot_config.name=[obj.h5_name(1:end-3),'_',variable_name,'_z_slice_at_z=',num2str(z),'_time_vertical.png'];
                plot_contour(data,plot_config);                
                
                %flip back the dimension of data{1}.z
                data{1}.z=data{1}.z';
                
                %add the data for the line plotting
                data_line{z_plot_ind}.y=data{1}.z(:,end);
                data_line{z_plot_ind}.x=obj.x_list;
                plot_config_line.legend_list{z_plot_ind+1}=['$z=',num2str(z),'$'];
                z_list_string=[z_list_string,num2str(z),'_'];
            end
            plot_config_line.legend_list{1}=1;
            plot_config_line.label_list={1,'$x$',label};
            plot_config_line.linewidth=5; 
            plot_config_line.name=[obj.h5_name(1:end-3),'_',variable_name,'_x_variation_at_z=',z_list_string(1:end-1),'.png'];
            plot_config_line.fontsize_legend=28;
            plot_line(data_line,plot_config_line);
            
        end
        
        function obj=rms_xt(obj,variable_name)
            switch variable_name
                case {'u','v','w','S','T','p','dy_T_mean_q','dy_S_mean_q'}
                    obj.(variable_name)=h5read_complex(obj.h5_name,['/tasks/',variable_name]);
                    obj.([variable_name,'_rms_xt'])=sqrt(mean(mean(obj.(variable_name).^2,2),3));
                    %plot_config.colormap='bluewhitered';
                case {'S_tot','T_tot'}
                    obj.(variable_name(1))=h5read_complex(obj.h5_name,['/tasks/',variable_name(1)]);
                    obj.([variable_name,'rms'])=sqrt(mean(mean(obj.(variable_name).^2,2),3))+obj.(['dy_',variable_name(1),'_mean'])*obj.z_list;
                    if obj.dy_T_mean==-1
                        T_mean_var=1+obj.dy_T_mean*obj.z_list;
                        T_mean_sign='+1-z';
                        dy_T_mean_sign='-1';
                    elseif obj.dy_T_mean==1
                        T_mean_var=obj.dy_T_mean*obj.z_list;
                        T_mean_sign='+z';
                        dy_T_mean_sign='+1';
                    elseif obj.dy_T_mean==0
                        T_mean_var=0*obj.z_list;
                        T_mean_sign='';
                        dy_T_mean_sign='';
                    end
                    
                    %Add fixed flux. 
                    if obj.flux_T
                        T_mean_var=mean(mean(obj.dy_T_mean_q))*obj.z_list;
                        T_mean_sign='+$\bar{\mathcal{T}}_{z,q}z$';
                        dy_T_mean_sign='+$\bar{\mathcal{T}}_{z,q}$';
                    end

                    if obj.dy_S_mean==-1
                        S_mean_var=1+obj.dy_S_mean*obj.z_list;
                        S_mean_sign='+1-z';
                        dy_S_mean_sign='-1';
                    elseif obj.dy_S_mean==1
                        S_mean_var=obj.dy_S_mean*obj.z_list;
                        S_mean_sign='+z';
                        dy_S_mean_sign='+1';
                    elseif obj.dy_S_mean==0
                        S_mean_var=0*obj.z_list;
                        S_mean_sign='';
                        dy_S_mean_sign='';
                    end
                    
                    %Add fixed flux. 
                    if obj.flux_S
                        S_mean_var=mean(mean(obj.dy_T_mean_q))*obj.z_list;
                        S_mean_sign='+$\bar{\mathcal{S}}_{z,q}z$';
                        dy_S_mean_sign='+$\bar{\mathcal{S}}_{z,q}$';
                    end
                    
                    %data{1}.z=squeeze(obj.(variable_name)(z_ind,:,:))+obj.(['dy_',variable_name(1),'_mean'])*obj.z_list(z_ind);
                    %plot_config.colormap='jet';
%                     plot_config.label_list={1,[variable_name(1),+],'$z$'};
            end
            data{1}.x=obj.([variable_name,'_rms_xt']);
            data{1}.y=obj.z_list;
            plot_config.print_size=[1,500,900];
            plot_config.label_list={1,['$',variable_name,'$'],'$z$'};
            plot_config.name=[obj.h5_name(1:end-3),'_',variable_name,'_rms_xt.png'];
            plot_line(data,plot_config);
            
        end
        
        
        function obj=total_xt_ave(obj,variable_name,x_ind,t_ind)
%             data{1}.y=obj.z_list/(2*pi/obj.k_opt);
%             dz=diff(obj.z_list); dz=dz(1);
%             z_list_full=[obj.z_list;obj.z_list(end)+dz];
            x_len=length(obj.x_list);
            time_len=length(obj.t_list);
            
            if nargin<4 || isempty(x_ind)
                %The default option, just average over the second half of
                %data...
                x_ind_begin=1;
                x_ind_end=x_len;
            else
                x_ind_begin=x_ind(1);
                if length(x_ind)==1
                    x_ind_end=x_len;
                else
                    x_ind_end=x_ind(2);
                end
            end
            
            
            if nargin<4 || isempty(t_ind)
                %The default option, just average over the second half of
                %data...
                t_ind_begin=time_len/2;
                t_ind_end=time_len;
            else
                t_ind_begin=t_ind(1);
                if length(t_ind)==1
                    t_ind_end=time_len;
                else
                    t_ind_end=t_ind(2);
                end
            end

            switch variable_name
                case {'T','S'}
                    if obj.(['flux_',variable_name])
                        if strcmp(obj.flow,'HB_benard_shear_periodic')
                            variable_data=h5read_complex(obj.h5_name,['/tasks/',variable_name,'_0']);
                            data{2}.x=mean(obj.(['dy_',variable_name,'_mean_q'])(1,t_ind_begin:t_ind_end),2)*obj.z_list;
                            data{1}.x=obj.(['Pe_',variable_name])*squeeze(mean(variable_data(:,t_ind_begin:t_ind_end),2))+mean(obj.(['dy_',variable_name,'_mean_q'])(1,t_ind_begin:t_ind_end),2)*obj.z_list;
                        else
                            variable_data=h5read_complex(obj.h5_name,['/tasks/',variable_name]);
                            data{2}.x=mean(mean(obj.(['dy_',variable_name,'_mean_q'])(1,1,t_ind_begin:t_ind_end),2),3)*obj.z_list;
                            data{1}.x=obj.(['Pe_',variable_name])*squeeze(mean(mean(variable_data(:,:,t_ind_begin:t_ind_end),2),3))+mean(mean(obj.(['dy_',variable_name,'_mean_q'])(1,1,t_ind_begin:t_ind_end),2),3)*obj.z_list;
                        end
                    else
                        variable_data=h5read_complex(obj.h5_name,['/tasks/',variable_name]);
                        data{2}.x=obj.(['dy_',variable_name,'_mean'])*obj.z_list;
                        data{1}.x=obj.(['Pe_',variable_name])*squeeze(mean(mean(variable_data(:,x_ind_begin:x_ind_end,t_ind_begin:t_ind_end),2),3))+obj.(['dy_',variable_name,'_mean'])*obj.z_list;
                    
                    end
                    
                    if obj.(['dy_',variable_name,'_mean'])>0
                        plot_config.legend_list={1,['$ z+\langle ',variable_name,'\rangle_{h,t}$'],['$z$']};
                    else obj.(['dy_',variable_name,'_mean'])<0
                        data{2}.x=data{2}.x+1;
                        data{1}.x=data{1}.x+1;
                        plot_config.legend_list={1,['$1-z+\langle ',variable_name,'\rangle_{h,t}$'],['$1-z$']};
                        if obj.(['flux_',variable_name])
                            plot_config.legend_list={1,['$1+\bar{\mathcal{T}}_{z,q} z+\langle ',variable_name,'\rangle_{h,t}(z)$'],['$1+\bar{\mathcal{T}}_{z,q} z$']};
                        end
                    end
                    plot_config.label_list={1,'','$z$'};

                case {'rho'}
                    %error('not ready');
                    variable_data_T=h5read_complex(obj.h5_name,['/tasks/T']);
                    variable_data_S=h5read_complex(obj.h5_name,['/tasks/S']);

                    R_rho_T2S=obj.Ra_T/obj.Ra_S2T;
                    data{1}.x=-obj.dy_T_mean*obj.z_list+1/R_rho_T2S*obj.dy_S_mean*obj.z_list;
                    data{2}.x=-(obj.Pe_T*squeeze(mean(mean(variable_data_T,2),3))+obj.dy_T_mean*obj.z_list)...
                        +1/R_rho_T2S*(obj.Pe_S*squeeze(mean(mean(variable_data_S,2),3))+obj.dy_S_mean*obj.z_list);
                    plot_config.legend_list={1,['$-(\bar{\mathcal{T}}_z z+Pe_T \langle ','T','\rangle_h)+R_\rho^{-1}(\bar{\mathcal{S}}_z z+Pe_S \langle ','S','\rangle_h)$'],['$-\bar{\mathcal{T}}_z z+R_\rho^{-1}\bar{\mathcal{S}}_z z$']};
                    plot_config.fontsize_legend=24;
                    plot_config.label_list={1,'','$z$'};

                case 'u'
                    u=h5read_complex(obj.h5_name,'/tasks/u');
                    syms z;
                    u_laminar=obj.F_sin/obj.ks^2*sin(obj.ks*z)...
                              +obj.F_sin_2ks/(2*obj.ks)^2*sin(2*obj.ks*z+obj.phase_2ks)...
                              +obj.F_sin_3ks/(3*obj.ks)^2*sin(3*obj.ks*z+obj.phase_3ks)...
                              +obj.F_sin_4ks/(4*obj.ks)^2*sin(4*obj.ks*z+obj.phase_4ks);
                    u_laminar_num=double(subs(u_laminar,z,obj.z_list));
                    data{1}.x=u_laminar_num;
                    data{2}.x=squeeze(mean(mean(u,2),3));
                    plot_config.legend_list={1,['$\bar{',variable_name,'}$'],['$\bar{',variable_name,'} +''\langle ',variable_name,'''\rangle_h$']};
                    plot_config.label_list={1,'','$z$'};

                case 'ww'
                    w=h5read_complex(obj.h5_name,'/tasks/w');
                    w=w(:,:,t_ind_begin:t_ind_end);
                    ww=w.*w;
                    data{1}.x=sqrt(mean(mean(ww,2),3));
                    data{2}.x=NaN*ones(size(obj.z_list));
                    %plot_config.label_list={1,'$\sqrt{\langle ww\rangle_{h,t}}$','$z$'};
                    plot_config.label_list={1,'','$z$'};
            end
%             if strcmp(obj.flow(1:7),'IFSC_2D')
%                 data{1}.y=obj.z_list/(2*pi/obj.k_opt);
%                 data{2}.y=obj.z_list/(2*pi/obj.k_opt);
%                 plot_config.label_list={1,'','$z/l_{opt}$'};
%             else
                data{1}.y=obj.z_list;
                data{2}.y=obj.z_list;
%             end
            plot_config.ytick_list=[1,0,0.2,0.4,0.6,0.8,1];
            plot_config.ylim_list=[1,round(min(data{1}.y)),round(max(data{1}.y))];
%             plot_config.label_list={1,'','$z/l_{opt}$'};
            plot_config.print_size=[1,1600,1600];
            plot_config.print=obj.print;
            plot_config.name=[obj.h5_name(1:end-3),'_',variable_name,'_total_xt_ave.png'];
            plot_config.linewidth=3;
            plot_line(data,plot_config);
            
            %data{1}.x=NaN; data{1}.y=NaN;
            %plot_config.label_list={1,plot_config.legend_list{3},'$z$'};
            %plot_config.legend_list={0};
            plot_config.print_size=[1,500,900];
            plot_config.name=[obj.h5_name(1:end-3),'_',variable_name,'_total_xt_ave_profile_only.png'];
            plot_config.fontsize_legend=16;
            plot_config.linewidth=3;
            plot_config.legend_list={0};
            plot_config.ytick_list=[1,0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2];
            plot_line(data,plot_config);
            
            plot_config.print=obj.print;
            plot_config.visible=0;
            plot_config.print_size=[1,500,900];

            switch variable_name
                case {'T','S'}
                if obj.video

                    for t_ind=1:length(obj.t_list)
                        data{1}.x=squeeze(mean(variable_data(:,:,t_ind),2))+obj.z_list;
                        data{2}.x=obj.z_list;
                        plot_config.legend_list={1,['$ z+\langle ',variable_name,'\rangle_{h}$'],['$z$']};
                        if obj.title_time
                            plot_config.title_list={1,['$t=$',num2str(round(obj.t_list(t_ind)))]};
                        else
                            plot_config.title_list={0};
                        end
                        if obj.no_ylabel
                           plot_config.label_list={1,'',''}; 
                        end
                        plot_config.name=[obj.h5_name(1:end-3),'_',variable_name,'_total_x_ave_',num2str(round(obj.t_list(t_ind))),'.png'];
                        snapshot(t_ind)=plot_line(data,plot_config);
                    end
                   plot_config.name=[obj.h5_name(1:end-3),'_',variable_name,'_total_xt_ave_video.mp4'];
                   plot_video(snapshot,plot_config);
                end
                otherwise
                    warning('No video for this case');
            end
            
        end

        function obj=flux_check(obj)
            w=h5read_complex(obj.h5_name,['/tasks/w']);
            T=h5read_complex(obj.h5_name,['/tasks/T']);
            obj.t_list=h5read_complex(obj.h5_name,'/scales/sim_time');


        end
        
        
        function obj=get_Nu(obj,variable_name,t_ind)
            obj.t_list=h5read_complex(obj.h5_name,'/scales/sim_time');
            time_len=length(obj.t_list);
            if nargin<3 || isempty(t_ind)
                %The default option, just average over the second half of
                %data...
                t_ind_begin=time_len/2;
                t_ind_end=time_len;
            else
                t_ind_begin=t_ind(1);
                if length(t_ind)==1
                    t_ind_end=time_len;
                else
                    t_ind_end=t_ind(2);
                end
            end
            
            if obj.flux_T~=0 && obj.flux_S~=0
                if strcmp(obj.store_variable,'all')
                    d_variable_data=h5read_complex(obj.h5_name,['/tasks/d_',variable_name]);
                else
                    error('Not ready, this interpolation is wrong!!');
                    z_list_cheb=obj.z_list*2-1;
                    variable_data=h5read_complex(obj.h5_name,['/tasks/',variable_name]);
                    [x,DM]=chebdif(length(z_list_cheb),1);
                    D1=DM(:,:,1);
                    for x_ind=1:size(variable_data,2)
                        for t_ind=1:size(variable_data,3)
                            variable_data_y=squeeze(variable_data(:,x_ind,t_ind));
                            variable_data_y_int=chebint(variable_data_y,z_list_cheb);
                            d_variable_data(:,x_ind,t_ind)=D1*variable_data_y_int;
                        end
                    end
                end

                d_variable_data_total_xt_ave=obj.(['Pe_',variable_name])*squeeze(mean(mean(d_variable_data(:,:,t_ind_begin:t_ind_end),2),3))+obj.(['dy_',variable_name,'_mean']);
                obj.(['d_',variable_name])=d_variable_data;
                switch variable_name
                    case 'T'
                        obj.Nu=d_variable_data_total_xt_ave;
                    case 'S'
                        obj.Nu_S=d_variable_data_total_xt_ave;
                end
            
            else
                %Update 2023/04/07, take the average over the first and the
                %end of the local maximum. 
                if strcmp(obj.flow,'HB_benard_shear_periodic')
                    obj.Nu_T_t=(-1)./(squeeze(obj.(['dy_',variable_name,'_mean_q'])(1,t_ind_begin:t_ind_end)));
                    dy_T_mean_q_tmp=(squeeze(obj.(['dy_',variable_name,'_mean_q'])(1,t_ind_begin:t_ind_end)));
                    t_ind_local_min=obj.t_list(t_ind_begin:t_ind_end);
                    
                    dy_T_mean_q_mid=(max(dy_T_mean_q_tmp)+min(dy_T_mean_q_tmp))/2;
                    ind_local_min=find(islocalmin(dy_T_mean_q_tmp).*(dy_T_mean_q_tmp<dy_T_mean_q_mid));
                    
                    t_ind_local_min=obj.t_list(t_ind_begin:t_ind_end);
                    if (~isempty(ind_local_min)) && length(ind_local_min)>2
                        obj.Nu_T_t=obj.Nu_T_t(ind_local_min(1):ind_local_min(end));
                        t_ind_local_min=t_ind_local_min(ind_local_min(1):ind_local_min(end));
                        dy_T_mean_q_tmp=dy_T_mean_q_tmp(ind_local_min(1):ind_local_min(end));

                    end
                    %d_variable_data_total_xt_ave=mean((-1)./(squeeze(obj.(['dy_',variable_name,'_mean_q'])(:,t_ind_begin:t_ind_end))));
                   d_variable_data_total_xt_ave=mean(obj.Nu_T_t);
                else
                    obj.Nu_T_t=(-1)./(squeeze(obj.(['dy_',variable_name,'_mean_q'])(1,1,t_ind_begin:t_ind_end)));
                    dy_T_mean_q_tmp=(squeeze(obj.(['dy_',variable_name,'_mean_q'])(1,1,t_ind_begin:t_ind_end)));
                    t_ind_local_min=obj.t_list(t_ind_begin:t_ind_end);
                    
                    %get the full time-history of the nu(t)
                    obj.Nu_T_t_full=obj.Nu_T_t;
                    
                    %get the index of several local minimum, use
                    %dy_T_mean_q to identify that is better than Nu. 
                    dy_T_mean_q_mid=(max(dy_T_mean_q_tmp)+min(dy_T_mean_q_tmp))/2;
                    ind_local_min=find(islocalmin(dy_T_mean_q_tmp).*(dy_T_mean_q_tmp<dy_T_mean_q_mid));
                    
                    if (~isempty(ind_local_min)) && length(ind_local_min)>=2
                        obj.Nu_T_t=obj.Nu_T_t(ind_local_min(1):ind_local_min(end));
                        t_ind_local_min=t_ind_local_min(ind_local_min(1):ind_local_min(end));
                        
                        %compute the Nusselt number by firstly computing
                        %time averaged mean temperature gradient and then
                        %take the recipocal. 
                        dy_T_mean_q_tmp=dy_T_mean_q_tmp(ind_local_min(1):ind_local_min(end));
                    end
                    %d_variable_data_total_xt_ave=mean((-1)./(squeeze(obj.(['dy_',variable_name,'_mean_q'])(:,:,t_ind_begin:t_ind_end))));
                    d_variable_data_total_xt_ave=mean(obj.Nu_T_t);
                end
                 switch variable_name
                    case 'T'
                        obj.Nu_average_after=d_variable_data_total_xt_ave;
                        obj.Nu=-1/mean(dy_T_mean_q_tmp);

                        Ra_T_q=obj.Ra_T;
                        obj.Ra_T_no_q=Ra_T_q/mean(obj.Nu);
                    case 'S'
                        obj.Nu_S=d_variable_data_total_xt_ave;
                        Ra_S2T_q=obj.Ra_S2T;
                        obj.Ra_S2T_no_q=Ra_S2T_q/mean(obj.Nu_S);
                 end
                 
            end
            
            %Update 2023/04/07, take from one local max to the end of local
            %maximum. 
            data{1}.x=t_ind_local_min; 
            data{1}.y=obj.Nu_T_t;
            
            %some markers for Pr=7, Ra_{T,q}=300000 case
%             data{2}.x=data{1}.x(362);
%             data{2}.y=data{1}.y(362);
%             data{3}.x=data{1}.x(366);
%             data{3}.y=data{1}.y(366);
%             data{4}.x=data{1}.x(367);
%             data{4}.y=data{1}.y(367);
%             data{5}.x=data{1}.x(368);
%             data{5}.y=data{1}.y(368);
%             plot_config.xlim_list=[1,3,4];


            %some markers for Pr=0.1, Ra_{T,q}=10800 case
%             [~,ind160]=min(abs(data{1}.x-160));
%             [~,ind166]=min(abs(data{1}.x-166));
%             [~,ind166p5]=min(abs(data{1}.x-166.5));
%             [~,ind167]=min(abs(data{1}.x-167));
%             data{2}.x=data{1}.x(ind160);
%             data{2}.y=data{1}.y(ind160);
%             data{3}.x=data{1}.x(ind166);
%             data{3}.y=data{1}.y(ind166);
%             data{4}.x=data{1}.x(ind166p5);
%             data{4}.y=data{1}.y(ind166p5);
%             data{5}.x=data{1}.x(ind167);
%             data{5}.y=data{1}.y(ind167);
%             plot_config.xlim_list=[1,150,300];
%             
            
            plot_config.user_color_style_marker_list={'k-','k*','bo','msquare','rx'};
            plot_config.Markerindex=3;
            plot_config.markersize=26;
            plot_config.label_list={1,'$t$','$nu(t)$'};
            plot_config.name=[obj.h5_name(1:end-3),'Nu_T_t.png'];
            plot_config.print=obj.print;
            plot_config.visible=obj.visible;
            plot_line(data,plot_config);
            
            data{1}.y=dy_T_mean_q_tmp;
            plot_config.label_list={1,'$t$','$\bar{\mathcal{T}}_{z,q}(t)$'};
            plot_config.name=[obj.h5_name(1:end-3),'dy_T_mean_q_T_t.png'];
            plot_line(data,plot_config);

            
        end
        
        function obj=elevator_growing(obj,variable_name)
            if strcmp(obj.flow,'HB_benard_shear_periodic')
                variable_name=[variable_name,'_hat'];
                obj.(variable_name)=h5read_complex(obj.h5_name,['/tasks/',variable_name]);
                data{1}.y=squeeze(max(obj.(variable_name)));
            else
                obj.(variable_name)=h5read_complex(obj.h5_name,['/tasks/',variable_name]);
                data{1}.y=squeeze(max(max(obj.(variable_name))));
            end
            data{1}.x=obj.t_list;
            plot_config.loglog=[0,1];
            switch variable_name
                case 'w_hat'
                    ylabel=['$max_z \;\hat{w}$'];
                case 'T_hat'
                    ylabel=['$max_z \; \hat{T}$'];
                case 'T'
                    ylabel=['$max_{z,x}\;T$'];
                case 'w'
                    ylabel=['$max_{z,x}\;w$'];
            end
            plot_config.label_list={1,'$t$',ylabel};
            plot_config.ytick_list=[1,0.1,1,10,10^2,10^3,10^4,10^5,10^6,10^7,10^8,10^9,10^10];
            plot_config.name=[obj.h5_name(1:end-3),variable_name,'_elevator_growing.png'];
            plot_line(data,plot_config);
            X=[data{1}.x(2:end),ones(size(data{1}.x(2:end)))];
            if size(data{1}.y,2)>1
                Y=log(data{1}.y(2:end))';
            else
                Y=log(data{1}.y(2:end));
            end
            fit=X\Y;
            obj.growth_rate=fit(1);
        end
        
        function obj=get_G(obj)
            %This is checking the residue
            
            if strcmp(obj.flow, 'HB_benard') ...
               & strcmp(obj.z_bc_w_left,'dirichlet')...
               & strcmp(obj.z_bc_w_right,'dirichlet')...
               & strcmp(obj.z_bc_u_v_left,'dirichlet')...
               & strcmp(obj.z_bc_u_v_right,'dirichlet')...
               & strcmp(obj.z_bc_T_left,'dirichlet')...
               & strcmp(obj.z_bc_T_right,'dirichlet')...
               & strcmp(obj.z_bc_S_left,'dirichlet')...
               & strcmp(obj.z_bc_S_right,'dirichlet')
           
               %[z_cheb, DM] = chebdif(obj.Nz+1, 2);
%                [D1,z_cheb]=cheb(512);

               %%-------------------------
               %% Note that the chebyshev grid point generated is not the same as the grid point from dedalus.
               %% this will lead to error... 
               %z_cheb=z_cheb(2:2:end);
               %D1=D1(2:2:end,2:2:end);
               %Lz_cheb=2;
               %D1=Lz_cheb/obj.Lz*fliplr(flipud(D1));
               %D2=Lz_cheb/obj.Lz*fliplr(flipud(DM(:,:,2)));
               %z_cheb_rescale=(obj.Lz/Lz_cheb*flipud(z_cheb)+0.5)
               %obj.z_list_res=obj.z_list-z_cheb_rescale;
               %%-------------------------
               
               %%Then we will type the residue of the following dedalus
               %%code
%                problem.add_equation('dz(u_tilde)-d_u_tilde=0')
%                problem.add_equation('dz(d_u_tilde)-(kx*p_hat+(kx*kx+ky*ky)*u_tilde)=0')
%                problem.add_equation('dz(v_tilde)-d_v_tilde=0')
%                problem.add_equation('dz(d_v_tilde)-(ky*p_hat+(kx*kx+ky*ky)*v_tilde)=0')
%                problem.add_equation('dz(w_hat)-(kx*u_tilde+ky*v_tilde)=0')
%                problem.add_equation('dz(p_hat)-(kx*d_u_tilde+ky*d_v_tilde-(kx*kx+ky*ky)*w_hat+Ra_T*T_hat-Ra_S2T*S_hat)=0')
%                problem.add_equation('dz(T_hat)-d_T_hat=0')
%                problem.add_equation('dz(S_hat)-d_S_hat=0')
%                problem.add_equation('dz(T_0)-d_T_0=0')
%                problem.add_equation('dz(S_0)-d_S_0=0')
%                problem.add_equation('dz(d_T_hat)-w_hat*dy_T_mean-(kx*kx+ky*ky)*T_hat=Pe_T*w_hat*d_T_0')
%                problem.add_equation('dz(d_T_0)=Pe_T*(2*kx*u_tilde*T_hat+2*ky*v_tilde*T_hat+2*w_hat*d_T_hat)')
%                problem.add_equation('dz(d_S_hat)-1/tau*w_hat*dy_S_mean-(kx*kx+ky*ky)*S_hat=Pe_S/tau*(w_hat*d_S_0)')   
%                problem.add_equation('dz(d_S_0)=Pe_S/tau*(2*kx*u_tilde*S_hat+2*ky*v_tilde*S_hat+2*w_hat*d_S_hat)')
%              

               %% Use poldif to generate the differential matrix corresponding to the grid point from dedalus
               D1=poldif(obj.z_list,1);
               
               %%get the local variable, and only select the first one just
               %%in case I am solving an IVP
               u_tilde=obj.u_tilde(:,1);
               d_u_tilde=obj.d_u_tilde(:,1);
               v_tilde=obj.v_tilde(:,1);
               d_v_tilde=obj.d_v_tilde(:,1);
               w_hat=obj.w_hat(:,1);
               p_hat=obj.p_hat(:,1);
               T_hat=obj.T_hat(:,1);
               d_T_hat=obj.d_T_hat(:,1);
               S_hat=obj.S_hat(:,1);
               d_S_hat=obj.d_S_hat(:,1);
               T_0=obj.T_0(:,1);
               d_T_0=obj.d_T_0(:,1);
               S_0=obj.S_0(:,1);
               d_S_0=obj.d_S_0(:,1);
               kx=obj.kx;
               ky=obj.ky;
               Pr=obj.Pr;
               Ra_T=obj.Ra_T;
               Ra_S2T=obj.Ra_S2T;
               tau=obj.tau;
               Pe_T=obj.Pe_T;
               Pe_S=obj.Pe_S;
               dy_T_mean=obj.dy_T_mean;
               dy_S_mean=obj.dy_S_mean;
               
               obj.G_num=[D1*(u_tilde)-d_u_tilde;
                       D1*(d_u_tilde)-(kx*p_hat+(kx*kx+ky*ky)*u_tilde);
                       D1*v_tilde-d_v_tilde;
                       D1*(d_v_tilde)-(ky*p_hat+(kx*kx+ky*ky)*v_tilde);
                       D1*(w_hat)-(kx*u_tilde+ky*v_tilde);
                       D1*(p_hat)-(kx*d_u_tilde+ky*d_v_tilde-(kx*kx+ky*ky)*w_hat+Ra_T*T_hat-Ra_S2T*S_hat);
                       D1*(T_hat)-d_T_hat;
                       D1*(S_hat)-d_S_hat;
                       D1*(T_0)-d_T_0;
                       D1*(S_0)-d_S_0;
                       D1*(d_T_hat)-w_hat*dy_T_mean-(kx*kx+ky*ky)*T_hat-Pe_T*diag(w_hat)*d_T_0;
                       D1*(d_T_0)-Pe_T*(2*kx*diag(u_tilde)*T_hat+2*ky*diag(v_tilde)*T_hat+2*diag(w_hat)*d_T_hat);
                       D1*(d_S_hat)-1/tau*w_hat*dy_S_mean-(kx*kx+ky*ky)*S_hat-Pe_S/tau*(diag(w_hat)*d_S_0);   
                       D1*(d_S_0)-Pe_S/tau*(2*kx*diag(u_tilde)*S_hat+2*ky*diag(v_tilde)*S_hat+2*diag(w_hat)*d_S_hat)
                       ];
               syms u_tilde d_u_tilde v_tilde d_v_tilde w_hat p_hat S_hat d_S_hat T_hat d_T_hat S_0 d_S_0 T_0 d_T_0;
               syms D1 
               %kx ky Pr Pe_S Pe_T Ra_T Ra_S2T tau dy_T_mean dy_S_mean;
               obj.var_list=[u_tilde d_u_tilde v_tilde d_v_tilde w_hat p_hat S_hat d_S_hat T_hat d_T_hat S_0 d_S_0 T_0 d_T_0];
               obj.G_sym=[D1*(u_tilde)-d_u_tilde;
                       D1*(d_u_tilde)-(kx*p_hat+(kx*kx+ky*ky)*u_tilde);
                       D1*v_tilde-d_v_tilde;
                       D1*(d_v_tilde)-(ky*p_hat+(kx*kx+ky*ky)*v_tilde);
                       D1*(w_hat)-(kx*u_tilde+ky*v_tilde);
                       D1*(p_hat)-(kx*d_u_tilde+ky*d_v_tilde-(kx*kx+ky*ky)*w_hat+Ra_T*T_hat-Ra_S2T*S_hat);
                       D1*(T_hat)-d_T_hat;
                       D1*(S_hat)-d_S_hat;
                       D1*(T_0)-d_T_0;
                       D1*(S_0)-d_S_0;
                       D1*(d_T_hat)-w_hat*dy_T_mean-(kx*kx+ky*ky)*T_hat-Pe_T*diag(w_hat)*d_T_0;
                       D1*(d_T_0)-Pe_T*(2*kx*diag(u_tilde)*T_hat+2*ky*diag(v_tilde)*T_hat+2*diag(w_hat)*d_T_hat);
                       D1*(d_S_hat)-1/tau*w_hat*dy_S_mean-(kx*kx+ky*ky)*S_hat-Pe_S/tau*(diag(w_hat)*d_S_0);   
                       D1*(d_S_0)-Pe_S/tau*(2*kx*diag(u_tilde)*S_hat+2*ky*diag(v_tilde)*S_hat+2*diag(w_hat)*d_S_hat)
                       ];
                   
            elseif strcmp(obj.flow, 'HB_porous') ...
               & strcmp(obj.z_bc_w_left,'dirichlet')...
               & strcmp(obj.z_bc_w_right,'dirichlet')...
               & strcmp(obj.z_bc_T_left,'dirichlet')...
               & strcmp(obj.z_bc_T_right,'dirichlet')...
               & strcmp(obj.z_bc_S_left,'dirichlet')...
               & strcmp(obj.z_bc_S_right,'dirichlet')
              
                
            else
                error('sG residue is not supported for this flow configuration')
            end
        end
        
        function obj=get_Gjac(obj)
            %This is computing the Jacobian matrix for the stability
            %analysis
            %make sure I get the G
            if obj.G_sym==0
                obj=get_G(obj);
            end   
            
            if strcmp(obj.flow, 'HB_benard') ...
               & strcmp(obj.z_bc_w_left,'dirichlet')...
               & strcmp(obj.z_bc_w_right,'dirichlet')...
               & strcmp(obj.z_bc_u_v_left,'dirichlet')...
               & strcmp(obj.z_bc_u_v_right,'dirichlet')...
               & strcmp(obj.z_bc_T_left,'dirichlet')...
               & strcmp(obj.z_bc_T_right,'dirichlet')...
               & strcmp(obj.z_bc_S_left,'dirichlet')...
               & strcmp(obj.z_bc_S_right,'dirichlet')
               
               syms D1;
               obj.Gjac_sym=jacobian(obj.G_sym,obj.var_list);
           
            %% Use poldif to generate the differential matrix corresponding to the grid point from dedalus
                
               D1_num=poldif(obj.z_list,1);
               I=eye(obj.Nz,obj.Nz);
%                O=zeros(obj.Nz,obj.Nz);
               %%get the local variable, and only select the first one just
               %%in case I am solving an IVP
               
               u_tilde=obj.u_tilde(:,1);
               d_u_tilde=obj.d_u_tilde(:,1);
               v_tilde=obj.v_tilde(:,1);
               d_v_tilde=obj.d_v_tilde(:,1);
               w_hat=obj.w_hat(:,1);
               p_hat=obj.p_hat(:,1);
               T_hat=obj.T_hat(:,1);
               d_T_hat=obj.d_T_hat(:,1);
               S_hat=obj.S_hat(:,1);
               d_S_hat=obj.d_S_hat(:,1);
               T_0=obj.T_0(:,1);
               d_T_0=obj.d_T_0(:,1);
               S_0=obj.S_0(:,1);
               d_S_0=obj.d_S_0(:,1);
               
               %This loop generate the numerical Jacobian matrix.
               for m =1:size(obj.Gjac_sym,1)
                   for n=1:size(obj.Gjac_sym,2)
                       m_block=(1+(m-1)*obj.Nz:m*obj.Nz);
                       n_block=(1+(n-1)*obj.Nz:n*obj.Nz);
                       try 
%                            double(obj.Gjac_sym(m,n));
                           obj.Gjac_num(m_block,n_block)=I*double(obj.Gjac_sym(m,n));
                       catch
                           if obj.Gjac_sym(m,n)==D1
                               obj.Gjac_num(m_block,n_block)=D1_num;
                           else
                               Gjac_num_function=matlabFunction(obj.Gjac_sym(m,n),'Vars',obj.var_list);
                               obj.Gjac_num(m_block,n_block)=...
                                   diag(Gjac_num_function((u_tilde),(d_u_tilde),...
                                                     (v_tilde),(d_v_tilde),...
                                                     (w_hat), (p_hat),...
                                                     (S_hat), (d_S_hat),...
                                                     (T_hat), (d_T_hat),...
                                                     (S_0), (d_S_0),...
                                                     (T_0), (d_T_0)));
                           end
                       end
                           
                   end
               end
               
%                
%                kx=obj.kx;
%                ky=obj.ky;
%                Pr=obj.Pr;
%                Ra_T=obj.Ra_T;
%                Ra_S2T=obj.Ra_S2T;
%                Pe_T=obj.Pe_T;
%                Pe_S=obj.Pe_S;
%                dy_T_mean=obj.dy_T_mean;
%                dy_S_mean=obj.dy_S_mean;
%                
               tau=obj.tau;
               Pr=obj.Pr;
               
               obj.Gjac_B=blkdiag(0*I,1/Pr*I,0*I,1/Pr*I,0*I,1/Pr*I,0*I,I,0*I,1/tau*I,0*I,I,0*I,1/tau*I);
               bc_ind=[1,obj.Nz,... %%b.c. index for u_tilde
                        2*obj.Nz+1,3*obj.Nz,...%b.c. index for v_tilde
                        4*obj.Nz+1,5*obj.Nz,...%b.c. index for w_hat
                        6*obj.Nz+1,7*obj.Nz,...%b.c. index for T_hat
                        8*obj.Nz+1,9*obj.Nz,...%b.c. index for S_hat
                        10*obj.Nz+1,11*obj.Nz,...%b.c. index for T_0
                        12*obj.Nz+1,13*obj.Nz,...%b.c. index for S_0
                        ];
                    
               obj.Gjac_num_bc=obj.Gjac_num;
               obj.Gjac_num_bc(bc_ind,:)=[];
               obj.Gjac_num_bc(:,bc_ind)=[];
               
               obj.Gjac_B_bc=obj.Gjac_B;
               obj.Gjac_B_bc(bc_ind,:)=[];
               obj.Gjac_B_bc(:,bc_ind)=[];
               
               obj.Gjac_eig=eig(obj.Gjac_num_bc,obj.Gjac_B_bc);
               %obj.Gjac_eig=eigs(obj.Gjac_num,obj.Gjac_B,1,'largestreal');
               %obj.Gjac_eig=sort(arrayfun(@(sigma)eigs(obj.Gjac_num,obj.Gjac_B,1,sigma),linspace(-10,10,10).'),'descend')
               toc
               %I=eye(obj.Nz,obj.Nz);
               %O=zeros(obj.Nz,obj.Nz);
%                obj.sG_jac=[D1,-I, O, O, O, O, O, O, O, O, O, O, O, O;
%                           -(kx*kx+ky*ky)*I, D1, O,O,O,-kx*I,O,O,O,O,O,O,O,O;
%                           O,O,D1,-I,O,O,O,O,O,O,O,O,O,O;
%                           O,O,-(kx*kx+ky*ky)*I,D1,O,-ky*I,O,O,O,O,O,O,O,O;
%                           -kx*I, O, -ky*I, O, D1, O,O,O,O,O,O,O,O,O;
%                           ];
%                
%                obj.sG=[D1*(u_tilde)-d_u_tilde;
%                        D1*(d_u_tilde)-(kx*p_hat+(kx*kx+ky*ky)*u_tilde);
%                        D1*v_tilde-d_v_tilde;
%                        D1*(d_v_tilde)-(ky*p_hat+(kx*kx+ky*ky)*v_tilde);
%                        D1*(w_hat)-(kx*u_tilde+ky*v_tilde);
%                        D1*(p_hat)-(kx*d_u_tilde+ky*d_v_tilde-(kx*kx+ky*ky)*w_hat+Ra_T*T_hat-Ra_S2T*S_hat);
%                        D1*(T_hat)-d_T_hat;
%                        D1*(S_hat)-d_S_hat;
%                        D1*(T_0)-d_T_0;
%                        D1*(S_0)-d_S_0;
%                        D1*(d_T_hat)-w_hat*dy_T_mean-(kx*kx+ky*ky)*T_hat-Pe_T*diag(w_hat)*d_T_0;
%                        D1*(d_T_0)-Pe_T*(2*kx*diag(u_tilde)*T_hat+2*ky*diag(v_tilde)*T_hat+2*diag(w_hat)*d_T_hat);
%                        D1*(d_S_hat)-1/tau*w_hat*dy_S_mean-(kx*kx+ky*ky)*S_hat-Pe_S/tau*(diag(w_hat)*d_S_0);   
%                        D1*(d_S_0)-Pe_S/tau*(2*kx*diag(u_tilde)*S_hat+2*ky*diag(v_tilde)*S_hat+2*diag(w_hat)*d_S_hat)
%                        ];
%                
           
           
            elseif strcmp(obj.flow, 'HB_porous') ...
               & strcmp(obj.z_bc_w_left,'dirichlet')...
               & strcmp(obj.z_bc_w_right,'dirichlet')...
               & strcmp(obj.z_bc_T_left,'dirichlet')...
               & strcmp(obj.z_bc_T_right,'dirichlet')...
               & strcmp(obj.z_bc_S_left,'dirichlet')...
               & strcmp(obj.z_bc_S_right,'dirichlet')
              
                
            else
                error('sG residue is not supported for this flow configuration')
            end
            
            
            
        end
        
        
        function Nu_kx=get_Nu_kx_Toomre(obj)
            Nu_kx=[0.759282	5.46973
                0.868806	5.84515
                1.14156	6.47067
                1.41393	6.98488
                1.65953	7.58263
                1.98612	8.12456
                2.42145	8.80539
                2.82923	9.34715
                3.1282	9.72217
                3.53561	10.1526
                4.13284	10.694
                4.72984	11.1658
                5.89592	11.859
                7.11565	12.413
                8.11822	12.7726
                9.49968	13.1315
                11.2325	13.3505
                12.6401	13.431
                13.3436	13.4017
                14.534	13.2879
                15.3725	13.1471
                16.4272	12.8944
                17.1301	12.6703
                18.1842	12.2786
                18.9409	11.9431
                19.9405	11.4401];
                        
        end
        
        
        function Nu_Ra=get_Nu_Ra_Toomre(obj)
            Nu_Ra=[1.02E+05	2.48817
                2.02E+05	3.10341
                4.11E+05	3.82865
                8.08E+05	4.7231
                1.48E+06	5.63688
                2.76E+06	6.65424
                5.26E+06	7.85553
                1.10E+07	9.58605
                2.52E+07	11.6997
                5.01E+07	13.6621
                1.31E+08	16.863
                3.03E+08	20.5816
                7.30E+08	25.3996
                1.47E+09	29.988
                3.53E+09	37.008
                9.36E+09	46.6941];

        end
        
        function porous_hewitt_2D=get_porous_hewitt_2D(obj)
            
            %shape of z over T_0, Ra=10000, from 
            %figure 3(a) of Hewitt DR, Neufeld JA, Lister JR. Ultimate regime of high Rayleigh number convection in a porous medium. Physical Review Letters. 2012 May 30;108(22):224503.
            porous_hewitt_2D.Ra_10000_z_T_0=[0.0124888	0.549645
                        0.0124888	0.546625
                        0.0133809	0.541652
                        0.0133809	0.538455
                        0.014273	0.534192
                        0.0160571	0.528508
                        0.0178412	0.524245
                        0.0223015	0.521403
                        0.0321142	0.5246
                        0.0392507	0.526909
                        0.0481713	0.527975
                        0.0633363	0.526554
                        0.0802855	0.5246
                        0.0990187	0.52389
                        0.134701	0.52389
                        0.161463	0.523179
                        0.187333	0.521581
                        0.204282	0.520515
                        0.227475	0.518739
                        0.255129	0.516607
                        0.280107	0.514476
                        0.29884	0.512877
                        0.321142	0.511101
                        0.345227	0.509325
                        0.366637	0.507726
                        0.392507	0.50595
                        0.417484	0.504352
                        0.44603	0.502575
                        0.464764	0.501332
                        0.482605	0.500266
                        0.521855	0.497602
                        0.546833	0.495826
                        0.577163	0.493694
                        0.600357	0.492096
                        0.628011	0.490142
                        0.654773	0.488366
                        0.689563	0.485702
                        0.723461	0.483393
                        0.752899	0.481439
                        0.809099	0.476821
                        0.859054	0.473446
                        0.896521	0.472202
                        0.915254	0.471314
                        0.933095	0.469538
                        0.947368	0.468472
                        0.956289	0.469361
                        0.963426	0.471137
                        0.976806	0.475933
                        0.981267	0.472913
                        0.983051	0.468828
                        0.984835	0.464565
                        0.984835	0.46119
                        0.986619	0.456927
                        0.987511	0.45
                        ];
                    
            %shape of z over T_0, Ra=20000, from 
            %figure 3(a) of Hewitt DR, Neufeld JA, Lister JR. Ultimate regime of high Rayleigh number convection in a porous medium. Physical Review Letters. 2012 May 30;108(22):224503.
            porous_hewitt_2D.Ra_20000_z_T_0=[0.00624442	0.549822
                    0.00624442	0.546625
                    0.00802855	0.540941
                    0.00802855	0.535613
                    0.00892061	0.531705
                    0.0115968	0.529218
                    0.0178412	0.533837
                    0.0214095	0.53579
                    0.0312221	0.534725
                    0.0374665	0.531883
                    0.0454951	0.52833
                    0.0517395	0.526554
                    0.0624442	0.524778
                    0.0749331	0.523179
                    0.0892061	0.522114
                    0.105263	0.521048
                    0.118644	0.520337
                    0.132025	0.519627
                    0.148082	0.518561
                    0.162355	0.517318
                    0.17752	0.516607
                    0.190901	0.515542
                    0.205174	0.514654
                    0.219447	0.513766
                    0.234612	0.512877
                    0.247993	0.512167
                    0.263158	0.511279
                    0.276539	0.510568
                    0.291704	0.50968
                    0.306869	0.508792
                    0.321142	0.507904
                    0.335415	0.507371
                    0.349688	0.506661
                    0.364853	0.505773
                    0.379126	0.505062
                    0.394291	0.504174
                    0.409456	0.503641
                    0.422837	0.502931
                    0.438002	0.50222
                    0.451383	0.50151
                    0.468332	0.500622
                    0.481713	0.499911
                    0.495986	0.499378
                    0.511151	0.49849
                    0.525424	0.497602
                    0.540589	0.497069
                    0.555754	0.496181
                    0.569135	0.495471
                    0.585192	0.49476
                    0.599465	0.493872
                    0.613738	0.493162
                    0.628903	0.492274
                    0.64496	0.491563
                    0.657449	0.49103
                    0.672614	0.489964
                    0.686887	0.489432
                    0.70116	0.488544
                    0.715433	0.487655
                    0.729706	0.486945
                    0.744871	0.485879
                    0.760036	0.485346
                    0.773417	0.484458
                    0.789474	0.48357
                    0.803747	0.482682
                    0.81802	0.481794
                    0.833185	0.480728
                    0.847458	0.480018
                    0.860839	0.47913
                    0.876004	0.478419
                    0.890277	0.477531
                    0.90455	0.476643
                    0.918822	0.475577
                    0.933095	0.474156
                    0.944692	0.472735
                    0.954505	0.470249
                    0.961641	0.467584
                    0.966994	0.46492
                    0.97413	0.462611
                    0.983051	0.465453
                    0.986619	0.468472
                    0.990187	0.468472
                    0.991079	0.46492
                    0.992864	0.459414
                    0.993756	0.454973
                    0.993756	0.450533
                    ];
            
            %shape of z over T_0, Ra=40000, from 
            %figure 3(a) of Hewitt DR, Neufeld JA, Lister JR. Ultimate regime of high Rayleigh number convection in a porous medium. Physical Review Letters. 2012 May 30;108(22):224503.       
            porous_hewitt_2D.Ra_40000_z_T_0=[0.00356824	0.549822
                    0.0044603	0.545204
                    0.0044603	0.540586
                    0.0044603	0.537211
                    0.00624442	0.533304
                    0.00981267	0.53952
                    0.0124888	0.542362
                    0.0178412	0.53952
                    0.0249777	0.533659
                    0.0330062	0.527975
                    0.0401427	0.524778
                    0.0463872	0.522114
                    0.0633363	0.517318
                    0.0749331	0.515364
                    0.088314	0.513766
                    0.101695	0.5127
                    0.115968	0.511989
                    0.131133	0.511279
                    0.146298	0.510746
                    0.176628	0.509858
                    0.206066	0.508792
                    0.235504	0.507904
                    0.265834	0.506838
                    0.295272	0.506128
                    0.325602	0.50524
                    0.355932	0.504352
                    0.384478	0.503464
                    0.4157	0.502575
                    0.445138	0.50151
                    0.475468	0.500622
                    0.504906	0.499556
                    0.536128	0.49849
                    0.563782	0.497602
                    0.595004	0.496536
                    0.638715	0.494938
                    0.669045	0.49405
                    0.714541	0.492451
                    0.759144	0.49103
                    0.803747	0.489787
                    0.849242	0.488366
                    0.892061	0.486767
                    0.907226	0.485879
                    0.920607	0.484458
                    0.933095	0.482504
                    0.950937	0.478064
                    0.964318	0.472558
                    0.973238	0.467229
                    0.981267	0.461012
                    0.986619	0.457282
                    0.991971	0.462256
                    0.99554	0.465986
                    0.996432	0.462078
                    0.996432	0.457105
                    0.996432	0.451066
                    ];
                
                %%Ra over Nu, 
                %from figure 2 of %shape Hewitt DR, Neufeld JA, Lister JR. Ultimate regime of high Rayleigh number convection in a porous medium. Physical Review Letters. 2012 May 30;108(22):224503.
                porous_hewitt_2D.Ra_Nu=[34.5828	1.00241
                    39.4326	0.984924
                    46.0482	1.28769
                    50.059	1.43349
                    57.7641	1.71409
                    69.9133	2.04981
                    89.8189	2.45157
                    110.015	2.82899
                    134.752	3.15002
                    159.247	3.44519
                    183.759	3.70119
                    219.77	3.97649
                    269.186	4.27246
                    310.619	4.58994
                    371.491	4.93134
                    482.99	5.79366
                    570.79	6.45068
                    724.609	8.1392
                    898.189	10.4542
                    1061.47	12.0627
                    1284.72	14.4253
                    1380.05	12.0689
                    1554.93	13.4363
                    1751.96	14.694
                    1997.64	16.6538
                    2505.9	20.2756
                    2996.98	23.8173
                    3990.58	30.5957
                    5005.9	37.2497
                    6058.76	44.5454
                    6991.33	51.3972
                    7971.76	58.2523
                    8981.89	64.8523
                    10000	72.1983
                    11001.5	78.9527
                    11959.7	84.8096
                    13001.4	92.7417
                    14133.8	99.6215
                    15924.7	112.906
                    16903.6	119.131
                    17942.6	125.698
                    19045.5	135.017
                    19976.4	142.457
                    21977	155.784
                    25059	176.562
                    27899.5	196.562
                    32193.8	222.784
                    35843	252.486
                    39905.8	281.086
                    ];
                
                
                %%Ra over T_rms, 
                %from figure 3(b) of Hewitt DR, Neufeld JA, Lister JR. Ultimate regime of high Rayleigh number convection in a porous medium. Physical Review Letters. 2012 May 30;108(22):224503.
                porous_hewitt_2D.Ra_T_rms=[1998.81	0.0998958
                    2488.21	0.0964592
                    4018.95	0.0937523
                    5002.98	0.0911571
                    7052.27	0.0892338
                    7985.71	0.0886005
                    9096.39	0.0867591
                    10059.4	0.0867459
                    11058.6	0.0867336
                    12085.4	0.0861226
                    13052.1	0.0861126
                    14096.1	0.0861026
                    16056.7	0.0854908
                    17136.9	0.0848916
                    18074.6	0.0848848
                    20106.8	0.0848711
                    32094.2	0.0836431
                    36342.2	0.0836274
                    40189.5	0.0836148
                    ];
                
                %%Ra over w_rms, 
                %from figure 3(b) of Hewitt DR, Neufeld JA, Lister JR. Ultimate regime of high Rayleigh number convection in a porous medium. Physical Review Letters. 2012 May 30;108(22):224503.
                porous_hewitt_2D.Ra_w_rms=[1998.81	0.0839935
                    2502.98	0.0839651
                    3995.24	0.0821782
                    5002.98	0.0821504
                    7010.65	0.0809776
                    8033.12	0.0815244
                    9042.7	0.0815099
                    10059.4	0.0814968
                    10993.4	0.0809227
                    12014.1	0.0809119
                    13052.1	0.0820318
                    14096.1	0.0814554
                    16056.7	0.0814395
                    17136.9	0.0808687
                    18074.6	0.0819916
                    19063.6	0.0814184
                    20106.8	0.0814119
                    32094.2	0.0819208
                    36127.7	0.0813401
                    40189.5	0.081893
                    ];
                
                
                %%Ra over u_rms, 
                %from figure 3(b) of Hewitt DR, Neufeld JA, Lister JR. Ultimate regime of high Rayleigh number convection in a porous medium. Physical Review Letters. 2012 May 30;108(22):224503.
                porous_hewitt_2D.Ra_u_rms=[1991.7	0.0257095
                    2493.85	0.0238291
                    3977.98	0.0245167
                    3979.46	0.0211938
                    4981.96	0.0209089
                    4982.76	0.0196437
                    6978.67	0.0203472
                    6980.03	0.0188526
                    7996.63	0.0188565
                    7997.77	0.0178387
                    8951.51	0.015639
                    9957.45	0.0151084
                    9954.96	0.0166489
                    10013.1	0.0172366
                    10012.2	0.0178448
                    11006.2	0.017359
                    11007.6	0.0165364
                    11009.4	0.0155357
                    11958.7	0.0158643
                    11960.6	0.0149044
                    12994	0.0142987
                    13949.7	0.0141032
                    15979.5	0.0148078
                    15981.4	0.0141061
                    15984.9	0.0129796
                    16954.1	0.0143044
                    17996.2	0.0116184
                    18975.2	0.0127157
                    19893.1	0.0129839
                    31932.4	0.0114683
                    31938.1	0.0106999
                    35949.1	0.010409
                    39986.8	0.0102672
                    ];
                
        end
        
        function porous_hewitt_3D=get_porous_hewitt_3D(obj)
            
            %from figure 7 of Hewitt DR, Neufeld JA, Lister JR. High Rayleigh number convection in a three-dimensional porous medium. Journal of fluid mechanics. 2014 Jun;748:879-95.
            porous_hewitt_3D.Ra_4000_z_T_0=[0.0153551	0.573563
                        0.0172745	0.563793
                        0.0211132	0.552299
                        0.0326296	0.537931
                        0.049904	0.537356
                        0.074856	0.53046
                        0.0959693	0.521264
                        0.12476	0.511494
                        0.149712	0.506322
                        0.1881	0.502874
                        0.24952	0.501724
                        0.351248	0.501149
                        0.43762	0.499425
                        0.59309	0.497126
                        0.71977	0.495977
                        0.808061	0.494828
                        0.838772	0.493103
                        0.869482	0.488506
                        0.884837	0.484483
                        0.907869	0.474713
                        0.930902	0.464943
                        0.948177	0.45977
                        0.955854	0.459195
                        0.96737	0.458621
                        0.976967	0.445977
                        0.980806	0.434483
                        0.982726	0.425862
                        ];
                    
            porous_hewitt_3D.Ra_8000_z_T_0=[0.00959693	0.573563
                        0.0115163	0.558046
                        0.0134357	0.553448
                        0.0307102	0.548851
                        0.0403071	0.543678
                        0.0518234	0.535057
                        0.059501	0.529885
                        0.0729367	0.522989
                        0.0902111	0.51954
                        0.117083	0.517816
                        0.138196	0.517816
                        0.165067	0.517241
                        0.18618	0.516667
                        0.213052	0.515517
                        0.234165	0.514368
                        0.262956	0.513218
                        0.28215	0.512069
                        0.309021	0.50977
                        0.330134	0.508621
                        0.355086	0.506322
                        0.37428	0.504598
                        0.403071	0.502299
                        0.422265	0.501149
                        0.451056	0.498851
                        0.47025	0.497701
                        0.495202	0.495977
                        0.516315	0.494253
                        0.545106	0.492529
                        0.56238	0.491379
                        0.591171	0.489655
                        0.610365	0.488506
                        0.639155	0.486782
                        0.658349	0.485632
                        0.68714	0.484483
                        0.706334	0.483333
                        0.733205	0.481609
                        0.754319	0.481034
                        0.783109	0.479885
                        0.802303	0.47931
                        0.829175	0.479885
                        0.850288	0.47931
                        0.879079	0.479885
                        0.898273	0.479885
                        0.923225	0.476437
                        0.93666	0.471264
                        0.946257	0.463793
                        0.955854	0.457471
                        0.96737	0.45
                        0.975048	0.447701
                        0.982726	0.447126
                        0.986564	0.438506
                        0.988484	0.432184
                        ];
                    
            porous_hewitt_3D.Ra_16000_z_T_0=[0.00575816	0.574138
                        0.0115163	0.566667
                        0.0191939	0.556322
                        0.024952	0.547126
                        0.0326296	0.536782
                        0.0422265	0.52931
                        0.049904	0.527011
                        0.0729367	0.525862
                        0.0978887	0.527011
                        0.134357	0.527011
                        0.161228	0.526437
                        0.197697	0.524138
                        0.222649	0.522989
                        0.259117	0.520115
                        0.28215	0.518391
                        0.318618	0.515517
                        0.357006	0.512069
                        0.403071	0.507471
                        0.43762	0.504598
                        0.462572	0.502299
                        0.497121	0.498276
                        0.522073	0.495402
                        0.556622	0.491954
                        0.59309	0.488506
                        0.641075	0.483908
                        0.698656	0.47931
                        0.760077	0.474138
                        0.821497	0.470115
                        0.857965	0.468966
                        0.882917	0.468966
                        0.919386	0.47069
                        0.946257	0.471264
                        0.965451	0.462644
                        0.973129	0.453448
                        0.978887	0.444828
                        0.982726	0.436207
                        0.990403	0.425862
                        ];

            %from figure 8(b) of Hewitt DR, Neufeld JA, Lister JR. High Rayleigh number convection in a three-dimensional porous medium. Journal of fluid mechanics. 2014 Jun;748:879-95.
            porous_hewitt_3D.Ra_4000_z_T_rms=[0   0
                        0.00718563	0.101401
                        0.011976	0.132399
                        0.0155689	0.144221
                        0.0227545	0.148424
                        0.0371257	0.141856
                        0.057485	0.131611
                        0.0778443	0.126883
                        0.134132	0.120315
                        0.182036	0.115587
                        0.227545	0.112434
                        0.282635	0.11007
                        0.379641	0.107968
                        0.445509	0.107443
                        0.508982	0.10718
                        0.582036	0.107443
                        0.652695	0.108231
                        0.741317	0.110858
                        0.814371	0.115324
                        0.861078	0.120053
                        0.902994	0.124518
                        0.932934	0.12951
                        0.947305	0.133713
                        0.956886	0.138441
                        0.967665	0.145271
                        0.977246	0.148424
                        0.984431	0.14317
                        0.989222	0.123468
                        0.992814	0.0916813
                        1	0
                        ];

            porous_hewitt_3D.Ra_8000_z_T_rms=[0	0
                    0.00479042	0.109545
                    0.00598802	0.124256
                    0.00718563	0.137391
                    0.00958084	0.138967
                    0.0167665	0.137653
                    0.0311377	0.125832
                    0.0431138	0.121366
                    0.0658683	0.117426
                    0.0862275	0.114011
                    0.111377	0.11007
                    0.132934	0.107706
                    0.160479	0.105867
                    0.184431	0.104553
                    0.215569	0.103765
                    0.239521	0.10324
                    0.269461	0.102715
                    0.293413	0.102452
                    0.324551	0.102189
                    0.349701	0.102189
                    0.378443	0.101926
                    0.403593	0.101926
                    0.433533	0.101926
                    0.458683	0.101664
                    0.488623	0.101664
                    0.513772	0.101664
                    0.542515	0.101664
                    0.567665	0.101664
                    0.598802	0.101664
                    0.622754	0.101926
                    0.652695	0.101926
                    0.676647	0.102189
                    0.707784	0.102452
                    0.731737	0.102715
                    0.762874	0.10324
                    0.785629	0.103765
                    0.815569	0.104553
                    0.839521	0.105867
                    0.867066	0.107968
                    0.88982	0.110333
                    0.91497	0.114273
                    0.932934	0.117688
                    0.955689	0.121891
                    0.97006	0.12662
                    0.976048	0.132662
                    0.982036	0.138179
                    0.991617	0.139492
                    0.991617	0.133975
                    0.994012	0.121891
                    0.99521	0.11007
                    0.99521	0.097986
                    1	0
                    ];
                
            porous_hewitt_3D.Ra_16000_z_T_rms=[0	0
                    0.00239521	0.106655
                    0.00359281	0.121366
                    0.00479042	0.141331
                    0.0131737	0.127671
                    0.0203593	0.122942
                    0.0335329	0.119002
                    0.0467066	0.114799
                    0.0610778	0.110858
                    0.0790419	0.107443
                    0.0982036	0.105079
                    0.11976	0.103503
                    0.142515	0.102452
                    0.164072	0.101926
                    0.188024	0.101664
                    0.210778	0.101401
                    0.25509	0.101138
                    0.301796	0.100876
                    0.346108	0.10035
                    0.392814	0.10035
                    0.437126	0.10035
                    0.483832	0.100613
                    0.552096	0.100613
                    0.620359	0.10035
                    0.665868	0.100613
                    0.711377	0.100613
                    0.756886	0.100876
                    0.802395	0.101401
                    0.847904	0.102189
                    0.870659	0.102977
                    0.892216	0.104291
                    0.912575	0.106392
                    0.931737	0.109282
                    0.946108	0.113222
                    0.959281	0.117426
                    0.972455	0.121366
                    0.983234	0.125832
                    0.989222	0.135552
                    0.994012	0.140543
                    0.99521	0.139229
                    0.996407	0.134238
                    0.996407	0.124518
                    0.996407	0.114273
                    0.997605	0.0943082
                    1	0
                    ];

            %from figure 8(b) of Hewitt DR, Neufeld JA, Lister JR. High Rayleigh number convection in a three-dimensional porous medium. Journal of fluid mechanics. 2014 Jun;748:879-95.
            porous_hewitt_3D.Ra_4000_z_w_rms=[0	0
                    0.00239521	0.0197023
                    0.00838323	0.0467601
                    0.011976	0.058056
                    0.0155689	0.0675131
                    0.0239521	0.0795972
                    0.0431138	0.0851138
                    0.0730539	0.0877408
                    0.111377	0.0898424
                    0.158084	0.0929947
                    0.197605	0.0956217
                    0.244311	0.0977233
                    0.31018	0.0990368
                    0.364072	0.100088
                    0.415569	0.100613
                    0.483832	0.100613
                    0.553293	0.100613
                    0.631138	0.0998249
                    0.692216	0.0987741
                    0.754491	0.0971979
                    0.8	0.0950963
                    0.850299	0.0916813
                    0.886228	0.089317
                    0.929341	0.0869527
                    0.962874	0.0835377
                    0.977246	0.0759194
                    0.985629	0.0635727
                    0.990419	0.0488616
                    0.99521	0.0296848
                    0.997605	0.0126095
                    1	0
                    ];

            porous_hewitt_3D.Ra_8000_z_w_rms=[0	0
                    0.00479042	0.0491243
                    0.0107784	0.0730298
                    0.0263473	0.0843257
                    0.0479042	0.08669
                    0.0742515	0.089317
                    0.097006	0.0922067
                    0.124551	0.0945709
                    0.148503	0.0958844
                    0.177246	0.0971979
                    0.232335	0.0985114
                    0.287425	0.0990368
                    0.342515	0.0992995
                    0.39521	0.0992995
                    0.45509	0.0995622
                    0.508982	0.0995622
                    0.567665	0.0995622
                    0.620359	0.0992995
                    0.673054	0.0990368
                    0.725749	0.0987741
                    0.782036	0.0977233
                    0.832335	0.0961471
                    0.857485	0.0948336
                    0.88503	0.092732
                    0.935329	0.0872154
                    0.979641	0.0819615
                    0.986826	0.0764448
                    0.991617	0.0661996
                    0.994012	0.0530648
                    0.994012	0.040718
                    1	0
                    ];
                
            porous_hewitt_3D.Ra_16000_z_w_rms=[0	0
                    0.00479042	0.0661996
                    0.00598802	0.0761821
                    0.0215569	0.0848511
                    0.0407186	0.0877408
                    0.057485	0.0911559
                    0.0754491	0.0937828
                    0.0958084	0.0956217
                    0.11976	0.0971979
                    0.162874	0.0982487
                    0.209581	0.0987741
                    0.25509	0.0990368
                    0.277844	0.0992995
                    0.318563	0.0995622
                    0.391617	0.0995622
                    0.437126	0.0995622
                    0.482635	0.0995622
                    0.550898	0.0995622
                    0.71018	0.0992995
                    0.732934	0.0992995
                    0.756886	0.0990368
                    0.802395	0.0987741
                    0.847904	0.0982487
                    0.891018	0.0958844
                    0.912575	0.0943082
                    0.932934	0.0914186
                    0.949701	0.0882662
                    0.967665	0.0853765
                    0.986826	0.0827496
                    0.991617	0.078021
                    0.99521	0.0683012
                    0.996407	0.058056
                    1	0
                    ];

            %from figure 8(b) of Hewitt DR, Neufeld JA, Lister JR. High Rayleigh number convection in a three-dimensional porous medium. Journal of fluid mechanics. 2014 Jun;748:879-95.
            porous_hewitt_3D.Ra_4000_z_u_rms=[0.00359281	0.0845884
                    0.0179641	0.0661996
                    0.0323353	0.0549037
                    0.0622754	0.0438704
                    0.097006	0.0357268
                    0.132934	0.028634
                    0.189222	0.0204904
                    0.232335	0.0168126
                    0.31018	0.0126095
                    0.392814	0.0105079
                    0.479042	0.00998249
                    0.562874	0.00998249
                    0.638323	0.0115587
                    0.711377	0.0136602
                    0.77006	0.0168126
                    0.814371	0.0207531
                    0.846707	0.0249562
                    0.88503	0.032049
                    0.922156	0.0401926
                    0.955689	0.0501751
                    0.980838	0.0656743
                    0.990419	0.0772329
                    0.997605	0.0851138
                    ];

            porous_hewitt_3D.Ra_8000_z_u_rms=[0.00598802	0.0701401
                    0.0131737	0.0559545
                    0.0239521	0.0464974
                    0.045509	0.0359895
                    0.0790419	0.0254816
                    0.117365	0.0176007
                    0.167665	0.0131349
                    0.219162	0.0110333
                    0.275449	0.00971979
                    0.329341	0.0089317
                    0.384431	0.0084063
                    0.438323	0.00788091
                    0.493413	0.00788091
                    0.548503	0.00761821
                    0.626347	0.00814361
                    0.681437	0.0089317
                    0.736527	0.00971979
                    0.789222	0.0110333
                    0.843114	0.0136602
                    0.892216	0.0189142
                    0.928144	0.0275832
                    0.958084	0.0378284
                    0.978443	0.0488616
                    0.988024	0.060683
                    0.994012	0.0722417
                    ];

            porous_hewitt_3D.Ra_16000_z_u_rms=[0.00239521	0.0743433
                    0.00479042	0.0643608
                    0.00598802	0.0593695
                    0.011976	0.0496497
                    0.0167665	0.0446585
                    0.0227545	0.0399299
                    0.0287425	0.0352014
                    0.0371257	0.0304729
                    0.045509	0.0257443
                    0.0562874	0.0215412
                    0.0706587	0.0176007
                    0.0874251	0.0144483
                    0.107784	0.0118214
                    0.129341	0.0105079
                    0.152096	0.00971979
                    0.196407	0.00788091
                    0.243114	0.00683012
                    0.287425	0.00630473
                    0.332934	0.00604203
                    0.379641	0.00604203
                    0.423952	0.00630473
                    0.469461	0.00630473
                    0.51497	0.00630473
                    0.560479	0.00604203
                    0.607186	0.00577933
                    0.652695	0.00577933
                    0.698204	0.00604203
                    0.742515	0.00656743
                    0.788024	0.00761821
                    0.833533	0.0084063
                    0.879042	0.0102452
                    0.900599	0.0123468
                    0.918563	0.0149737
                    0.934132	0.0183888
                    0.948503	0.0225919
                    0.956886	0.0270578
                    0.965269	0.0315236
                    0.972455	0.0365149
                    0.978443	0.0412434
                    0.984431	0.0462347
                    0.988024	0.0512259
                    0.992814	0.0559545
                    0.992814	0.0609457
                    0.996407	0.0661996
                    ];

            %from figure 2 of Hewitt DR, Neufeld JA, Lister JR. High Rayleigh number convection in a three-dimensional porous medium. Journal of fluid mechanics. 2014 Jun;748:879-95.
            porous_hewitt_3D.Ra_Nu=[45.3231	1.24819
                    50.1955	1.43731
                    64.135	1.86782
                    71.0299	2.04514
                    90.293	2.4767
                    100	2.60468
                    129.082	3.06039
                    141.506	3.28406
                    180.803	3.74373
                    201.266	3.93719
                    254.545	5.1683
                    281.91	5.43538
                    301.256	5.6022
                    362.041	6.13403
                    403.015	5.43538
                    400.963	6.45101
                    501.955	6.58235
                    600.165	8.29922
                    699.503	9.46085
                    802.891	9.94974
                    802.891	10.7851
                    898.332	11.0046
                    902.93	11.9285
                    1005.12	13.0609
                    1090.67	14.1574
                    1251.88	14.889
                    1496.81	18.0309
                    1753.49	20.5547
                    2002.41	23.1967
                    2506.76	28.0917
                    2997.22	33.3408
                    4009.63	42.8929
                    5019.55	51.9441
                    6001.65	61.0321
                    7030.83	70.991
                    7988.02	80.9275
                    9029.3	89.5076
                    ];
            
            %from figure 8(a) of Hewitt DR, Neufeld JA, Lister JR. High Rayleigh number convection in a three-dimensional porous medium. Journal of fluid mechanics. 2014 Jun;748:879-95.
            porous_hewitt_3D.Ra_T_rms=[98.2752	0.180517
                    196.627	0.161938
                    294.15	0.149605
                    388.998	0.148205
                    491.96	0.145341
                    688.189	0.135603
                    787.125	0.130318
                    890.148	0.129055
                    995.418	0.129086
                    1076.85	0.119171
                    1204.14	0.120399
                    1472.57	0.119252
                    1741.62	0.115766
                    1991.71	0.114646
                    2463.21	0.112425
                    3012.75	0.10806
                    3984.4	0.105982
                    4982.77	0.104977
                    5959.12	0.102936
                    6968.89	0.101946
                    7880.63	0.101973
                    8911.65	0.102
                    9965.56	0.102025
                    11917.7	0.101048
                    13781.5	0.10108
                    15937.7	0.100105
                    19931.2	0.0991556
                    ];

            porous_hewitt_3D.Ra_w_rms=[98.4043	0.137762
                196.857	0.127351
                294.464	0.120031
                394.212	0.0954047
                493.014	0.0935589
                696.838	0.103488
                788.005	0.103516
                891.1	0.103543
                985.311	0.105661
                1089.8	0.101535
                1191.68	0.102577
                1232.2	0.104659
                1473.5	0.1047
                1723.18	0.103692
                1992.68	0.103725
                2464.4	0.101715
                2979.98	0.102781
                3984.98	0.102846
                4983.49	0.101871
                5893.45	0.100893
                6892.09	0.0999223
                7970.37	0.0989582
                8912.52	0.0999781
                9966.52	0.100002
                11918.3	0.100041
                13782.2	0.100073
                15760.5	0.100102
                19931.2	0.0991556
                ];

            porous_hewitt_3D.Ra_u_rms=[99.0927	0.00763127
                200.629	0.0139358
                300.34	0.0212384
                401.155	0.0192272
                501.256	0.0182974
                699.44	0.0138345
                799.542	0.0137007
                904.757	0.0167423
                1000.17	0.016414
                1090.9	0.010052
                1192.22	0.00928013
                1248.05	0.0118017
                1508.18	0.0112301
                1761.71	0.0096675
                1991.8	0.00986565
                1986.69	0.00580348
                2485.44	0.00709343
                2487.37	0.00832573
                2489.3	0.00977211
                2975.21	0.00958211
                4020.41	0.00958838
                4017.68	0.00833441
                4967.21	0.00885446
                4965.76	0.00833825
                4964.08	0.00777392
                6001.05	0.00801424
                6003.96	0.00885811
                6938.84	0.00842828
                6935.14	0.0075494
                7924.22	0.00683219
                7929.22	0.00778184
                8965.26	0.00802124
                9911.67	0.00802299
                9907.83	0.00740548
                11974.6	0.00726167
                13837.9	0.00677233
                15810.6	0.00606789
                19978.4	0.00583267
                ];
            
            porous_hewitt_3D.Ra_kpp=[826.858	29.2733
                1000	31.3059
                1372.82	45.5924
                1787.7	48.1078
                2727.59	62.0873
                3705.16	79.0604
                4824.88	90.4208
                6023.05	107.664
                7284.26	118.273
                8625.41	137.096
                10000	152.642
                13022	186.697
                16255.8	207.866
                19659.7	247.506
                26989.4	290.775
                ];

        end
        
        function porous_otero_bound=get_porous_otero_bound(obj)
            porous_otero_bound.Ra_Nu_bound=[8.19727	1
                    13.3407	1
                    24.7062	1
                    33.2891	1
                    39.0272	1
                    42.6792	1.14941
                    48.5658	1.40239
                    56.374	1.71104
                    68.7717	2.19408
                    81.4313	2.65049
                    92.663	3.0465
                    105.444	3.39874
                    135.187	4.10574
                    155.37	4.58044
                    182.151	5.31743
                    228.935	6.61812
                    276.52	7.9948
                    365.25	10.6679
                    512.097	14.9605
                    867.22	25.3447
                    1240.3	36.6197
                    1454.09	42.9368
                    2142.57	63.2843
                    2802.09	82.7802
                    3889.81	114.941
                    5618.76	166.074
                    7799.85	230.596
                    10000	295.694];

            
        end
        
        function porous_wen_2D=get_porous_wen_2D(obj)
            porous_wen_2D.Ra_Nu_DNS=[1959.9	16.1351
                    2431.05	19.8296
                    3015.36	23.9539
                    3842.38	30.4685
                    4830.39	36.8042
                    6072.65	45.2299
                    7634.64	56.5502
                    9470.28	70.7063
                    18818.9	138.193
                    37396.1	270.092
                    47015	337.692
                    58319	422.225
                    73319.7	527.902
                    92178.9	660.027
                    ];
            
            porous_wen_2D.Ra_Nu_steady=[2430.73	18.5093
                    3014.55	20.8701
                    3789.31	23.9402
                    4763.35	27.9391
                    6068.38	31.5006
                    7526.17	36.1357
                    9333.84	40.7448
                    11891.4	46.7369
                    14947.1	52.6965
                    18788.7	60.4483
                    23617.5	69.3405
                    29687.3	79.5407
                    37318.4	92.8268
                    ];
        end
        
        function porous_pirozzoli_3D=get_porous_pirozzoli_3D(obj)
            porous_pirozzoli_3D.Ra_Nu_DNS=[10000,99.84;
                                           20000,193.71;
                                           30000,281.14;
                                           40000,370.17;
                                           80000,709];
            
        end
        
        function porous_trevisan_tau=get_porous_trevisan_tau(obj)
            porous_trevisan_tau.Sh_Le_Ra_50=[0.2,1.02;
                                            0.4,1.08;
                                            1,1.42;
                                            2,2.12;
                                            4,3.19;
                                            10,5.19;
                                            20,7.49;
                                            40,10.62;
                                            100,14.71];
            porous_trevisan_tau.Sh_Le_Ra_100=[0.2,1.15;
                                            0.4,1.50;
                                            1,2.68;
                                            2,4.04;
                                            4,5.84;
                                            10,9.55;
                                            20,13.73;
                                            40,18.50];
            porous_trevisan_tau.Sh_Le_Ra_200=[0.1,1.09;
                                            0.2,1.32;
                                            0.4,2.01;
                                            1,4.06;
                                            2,6.17;
                                            4,8.93;
                                            10,14.71;
                                            20,21.52;
                                            40,30.11];
            porous_trevisan_tau.Sh_Le_Ra_400=[0.1,1.16;
                                            0.2,1.58;
                                            0.4,2.79;
                                            1,6.26;
                                            2,9.83;
                                            4,14.47;
                                            10,24.03;
                                            20,32.76];
            porous_trevisan_tau.Sh_Le_Ra_1000=[0.02,1.01;
                                            0.04,1.04;
                                            0.1,1.24;
                                            0.2,1.87;
                                            0.4,3.87;
                                            1,10.62;
                                            2,18.32;
                                            4,27.84;
                                            10,46.29];
        end
        
        function porous_rosenberg_tau=get_porous_rosenberg_tau(obj)
            porous_rosenberg_tau.Nu_Le_Ra_100=[10.1601	2.3531
                                                30.1779	2.43848
                                                50.1957	2.44032
                                                100.16	2.47005
                                                ];
            porous_rosenberg_tau.Nu_Le_Ra_150=[10	2.87112
                                                19.9288	2.94728
                                                30.0178	2.98991
                                                40.1068	2.98242
                                                50.1957	3.00834
                                                100	3.0716
                                                ];
            porous_rosenberg_tau.Nu_Le_Ra_300=[10	4.55847
                                                19.9288	4.61793
                                                29.8577	4.61055
                                                40.1068	4.51941
                                                50.0356	4.56216
                                                100.16	4.57506
                                                ];
            porous_rosenberg_tau.Nu_Le_Ra_600=[10	7.28162
                                                30.0178	7.26676
                                                50.0356	7.28531
                                                100.32	7.3148
                                                ];

            porous_rosenberg_tau.Sh_Le_Ra_100=[30.0065	14.1058
                                                50.0044	17.4086
                                                100.165	25.2837
                                                ];
            porous_rosenberg_tau.Sh_Le_Ra_150=[10.0012	10.6767
                                                20.0329	15.2547
                                                30.021	18.5302
                                                40.1235	20.7513
                                                50.0315	23.63
                                                100.624	33.3952
                                                ];
            porous_rosenberg_tau.Sh_Le_Ra_300=[9.96697	15.4732
                                                20.1263	21.0477
                                                30.0367	24.8794
                                                40.1468	28.7889
                                                50.261	31.2109
                                                100.679	45.3298
                                                ];
            porous_rosenberg_tau.Sh_Le_Ra_600=[9.93587	26.5573
                                                30.1923	45.8387
                                                50.3155	57.5053
                                                100.774	77.3766
                                                ];
            
        end
        
        function porous_rosenberg_R_rho_S2T=get_porous_rosenberg_R_rho_S2T(obj)
            porous_rosenberg_R_rho_S2T.Nu_R_rho_Ra_100=[0.00E+00	2.61741
                                                    0.0997811	2.46966
                                                    0.19979	2.36412
                                                    0.301229	2.13193
                                                    0.40056	1.68865
                                                    ];
            porous_rosenberg_R_rho_S2T.Nu_R_rho_Ra_150=[-2.88E-04	3.25066
                                                    0.100434	3.1029
                                                    0.201156	2.95515
                                                    0.301177	2.72296
                                                    0.399781	2.46966
                                                    ];
            porous_rosenberg_R_rho_S2T.Nu_R_rho_Ra_300=[2.64E-04	5.02375
                                                    0.100283	4.81266
                                                    0.20102	4.49604
                                                    0.301046	4.20053
                                                    0.400363	3.90501
                                                    ];
            porous_rosenberg_R_rho_S2T.Nu_R_rho_Ra_600=[5.61E-06	7.93668
                                                    0.100036	7.59894
                                                    0.251861	7.00792
                                                    0.400148	6.33245
                                                    ];

            porous_rosenberg_R_rho_S2T.Sh_R_rho_Ra_100=[0	12.0533
                                                0.0998217	11.3867
                                                0.20107	10.613
                                                0.300178	9.41329
                                                0.4	6.82667
                                                ];
            porous_rosenberg_R_rho_S2T.Sh_R_rho_Ra_150=[0	20.48
                                                0.100535	19.4932
                                                0.199643	18.4001
                                                0.300178	15.92
                                                0.400713	13.3331
                                                ];
            porous_rosenberg_R_rho_S2T.Sh_R_rho_Ra_300=[0	23.4667
                                                0.0998217	22.2667
                                                0.200357	20.5332
                                                0.300178	18.48
                                                0.4	16
                                                ];
            porous_rosenberg_R_rho_S2T.Sh_R_rho_Ra_600=[7.13E-04	41.4931
                                                0.100535	39.7599
                                                0.250267	36.1999
                                                0.4	31.68
                                                ];
            
        end
        
        function finger_width_yang=get_finger_width_yang(obj)
            finger_width_yang=[3.37E+05	0.672211
                6.39E+05	0.584655
                3.37E+06	0.384664
                3.06E+07	0.224896
                3.48E+08	0.124615
                1.47E+09	0.087443
                5.81E+09	0.0626907
                2.37E+10	0.0449449
                6.19E+10	0.0358739
                1.56E+11	0.0283279
                3.96E+11	0.0226106
                1.03E+12	0.0180472
                2.61E+12	0.0144048
                4.64E+12	0.0123948
                ];
        end
        
        function Nu_salt_finger_yang=get_Nu_salt_finger_yang(obj)
           
           %Data from Yang Y, Van Der Poel EP, Ostilla-Mónico R, Sun C, Verzicco R, Grossmann S, Lohse D. Salinity transfer in bounded double diffusive convection. Journal of fluid mechanics. 2015 Apr;768:476-91. 
           %THe first column is the density ratio
           Nu_salt_finger_yang.Ra_T_1e5_Nu=[10,1.0052;
                                            5,1.0125;
                                            2,1.035;
                                            1,1.0775;
                                            0.5,1.1706;
                                            0.2,1.4265;
                                            0.1,1.8826];
           Nu_salt_finger_yang.Ra_T_1e5_Nu_S=[10,8.6347;
                                            5,11.064;
                                            2,15.05;
                                            1,17.854;
                                            0.5,22.107;
                                            0.2,29.259;
                                            0.1,35.342]; 
           Nu_salt_finger_yang.Ra_T_1e5_Re=[10,0.1107;
                                            5,0.1814;
                                            2,0.3521;
                                            1,0.5254;
                                            0.5,0.8275;
                                            0.2,1.4652;
                                            0.1,2.3496]; 
           Nu_salt_finger_yang.Ra_T_1e6_Nu=[10,1.0116;
                                            5,1.0277;
                                            2,1.0789;
                                            1,1.1791;
                                            0.5,1.3929;
                                            0.2,2.0197;
                                            0.1,3.0231];
           Nu_salt_finger_yang.Ra_T_1e6_Nu_S=[10,17.352;
                                            5,22.037;
                                            2,29.542;
                                            1,35.516;
                                            0.5,42.5;
                                            0.2,56.184;
                                            0.1,68.098]; 
           Nu_salt_finger_yang.Ra_T_1e6_Re=[10,0.2773;
                                            5,0.4584;
                                            2,0.8727;
                                            1,1.3349;
                                            0.5,2.0749;
                                            0.2,3.8484;
                                            0.1,6.2142]; 

        end
        
        function porous_hewitt_2_layer=get_porous_hewitt_2_layer(obj)
            %Here is the data digitized from Hewitt DR, Neufeld JA, Lister JR. High Rayleigh number convection in a porous medium containing a thin low-permeability layer. Journal of fluid mechanics. 2014 Oct;756:844-69.
            
            %Figure 8(a) of hewitt (2014) et al. 
            porous_hewitt_2_layer.Omega_Nu=[0.00250719	37.4737
                    0.00497608	37.4737
                    0.00628603	37.2632
                    0.0100312	37.4737
                    0.012476	37.2632
                    0.0202217	37.4737
                    0.0251502	37.6842
                    0.0303202	37.6842
                    0.0401345	37.6842
                    0.0499161	38.3158
                    0.0601772	38.1053
                    0.0809066	39.3684
                    0.100625	38.1053
                    0.12131	38.9474
                    0.160577	40.4211
                    0.199713	41.0526
                    0.240767	41.2632
                    0.402598	38.7368
                    0.642463	36.2105
                    0.799045	36.6316
                    0.799045	33.0526
                    1.27511	29.2632
                    1.61078	30.7368
                    1.61078	25.8947
                    2.57048	25.6842
                    3.19695	23.3684
                    5.18178	13.0526
                    6.44468	13.2632
                    10.1254	12.4211
                    20.0961	12.8421
                    30.1318	12.4211
                    40.5115	12.4211
                    50.3849	12.4211
                    100	12.2105
                    ];
                
            %Figure 8(b) of hewitt (2014) et al., T_rms, w_rms and u_rms 
            porous_hewitt_2_layer.Omega_T_rms=[0.00254335	0.0936508
                    0.00504316	0.0952381
                    0.00627043	0.0944444
                    0.01	0.0968254
                    0.0126285	0.0992063
                    0.0201397	0.1
                    0.0250408	0.101587
                    0.0297148	0.103175
                    0.0399348	0.107937
                    0.0504316	0.111905
                    0.0598449	0.114286
                    0.0804277	0.119841
                    0.1	0.127778
                    0.120526	0.134127
                    0.159479	0.138095
                    0.201397	0.146032
                    0.242737	0.153175
                    0.405609	0.165873
                    0.646861	0.196032
                    0.804277	0.2
                    1.28265	0.203175
                    1.59479	0.218254
                    1.59479	0.19127
                    2.54335	0.196032
                    3.21186	0.169841
                    5.12223	0.047619
                    6.46861	0.0507937
                    10	0.0365079
                    20.1397	0.0412698
                    29.7148	0.0380952
                    39.9348	0.0380952
                    49.6531	0.0373016
                    101.568	0.0333333
                    ];
            porous_hewitt_2_layer.Omega_w_rms=[0.00254335	0.081746
                    0.00504316	0.0801587
                    0.00627043	0.0801587
                    0.01	0.0785714
                    0.0128265	0.0761905
                    0.0198288	0.0746032
                    0.0254335	0.0738095
                    0.0301807	0.0722222
                    0.0405609	0.068254
                    0.0504316	0.065873
                    0.0598449	0.0634921
                    0.0804277	0.0634921
                    0.1	0.0547619
                    0.120526	0.0539683
                    0.161979	0.0563492
                    0.201397	0.0531746
                    0.238989	0.0507937
                    0.405609	0.0444444
                    0.646861	0.0309524
                    0.804277	0.0246032
                    0.804277	0.031746
                    1.28265	0.0214286
                    1.61979	0.0166667
                    1.59479	0.0230159
                    2.58322	0.015873
                    3.16228	0.0126984
                    5.20255	0.0031746
                    6.46861	0.0031746
                    10	0.0031746
                    19.8288	0.0031746
                    30.1807	0.00238095
                    39.9348	0.00238095
                    49.6531	0.00238095
                    100	0
                    ];
            porous_hewitt_2_layer.Omega_u_rms=[0.00250408	0.0190476
                    0.00504316	0.0222222
                    0.00627043	0.0198413
                    0.01	0.0222222
                    0.0126285	0.0253968
                    0.0204555	0.0253968
                    0.0254335	0.0269841
                    0.0301807	0.0293651
                    0.0399348	0.0365079
                    0.0504316	0.0396825
                    0.0598449	0.0404762
                    0.0804277	0.0412698
                    0.1	0.0507937
                    0.120526	0.052381
                    0.159479	0.05
                    0.201397	0.0531746
                    0.238989	0.0579365
                    0.399348	0.0555556
                    0.636875	0.0746032
                    0.804277	0.0730159
                    0.804277	0.0769841
                    1.28265	0.0746032
                    1.59479	0.0706349
                    1.59479	0.0753968
                    2.54335	0.0777778
                    3.21186	0.0793651
                    5.12223	0.0531746
                    6.36875	0.0547619
                    10	0.0484127
                    19.8288	0.0555556
                    29.7148	0.052381
                    39.9348	0.0515873
                    49.6531	0.0484127
                    100	0.0492063
                    ];
                
            %temperature and RMS profile over z direction...
            %from figure 7(a)
            porous_hewitt_2_layer.z_T_0_Omega_0p04=[0.25	0.991736
                    0.392857	0.983471
                    0.470238	0.975207
                    0.494048	0.966942
                    0.488095	0.929752
                    0.488095	0.859504
                    0.482143	0.789256
                    0.47619	0.706612
                    0.470238	0.619835
                    0.470238	0.541322
                    0.488095	0.508264
                    0.5	0.5
                    0.517857	0.483471
                    0.529762	0.450413
                    0.529762	0.404959
                    0.52381	0.334711
                    0.517857	0.239669
                    0.511905	0.152893
                    0.505952	0.0661157
                    0.505952	0.0413223
                    0.505952	0.0330579
                    0.541667	0.0206612
                    0.589286	0.0165289
                    0.684524	0.00826446
                    0.791667	0.00826446
                    ];
            porous_hewitt_2_layer.z_T_rms_Omega_0p04=[0	1
                    0.072619	0.991736
                    0.1083334	0.983471
                    0.122619	0.966942
                    0.120238	0.942149
                    0.1130952	0.917355
                    0.1035714	0.859504
                    0.0964286	0.760331
                    0.0940476	0.68595
                    0.0940476	0.603306
                    0.1011904	0.528926
                    0.1071428	0.5
                    0.0988096	0.450413
                    0.095238	0.413223
                    0.0940476	0.347107
                    0.0964286	0.247934
                    0.1035714	0.140496
                    0.1142858	0.0743802
                    0.1214286	0.0371901
                    0.1035714	0.0123967
                    0.075	0.00413223
                    0	0
                     ];
            porous_hewitt_2_layer.z_w_rms_Omega_0p04=[0	1
                    0.0345238	0.991736
                    0.0535714	0.983471
                    0.0630952	0.979339
                    0.070238	0.96281
                    0.0761904	0.900826
                    0.079762	0.826446
                    0.0821428	0.764463
                    0.0833334	0.698347
                    0.0833334	0.623967
                    0.0821428	0.566116
                    0.0738096	0.516529
                    0.070238	0.5
                    0.077381	0.46281
                    0.0833334	0.396694
                    0.0833334	0.334711
                    0.0833334	0.272727
                    0.0809524	0.202479
                    0.0785714	0.140496
                    0.0738096	0.0785124
                    0.0642858	0.0206612
                    0.0380952	0.00826446
                    0	0
                    ];
            porous_hewitt_2_layer.z_u_rms_Omega_0p04=[0.079762	0.987603
                    0.0595238	0.96281
                    0.0404762	0.921488
                    0.027381	0.85124
                    0.0214286	0.77686
                    0.01904762	0.690083
                    0.01785714	0.586777
                    0.0238096	0.533058
                    0.0333334	0.495868
                    0.0214286	0.454545
                    0.01785714	0.429752
                    0.01785714	0.384298
                    0.01904762	0.297521
                    0.0214286	0.231405
                    0.0261904	0.152893
                    0.0392858	0.0785124
                    0.0583334	0.0371901
                    0.0833334	0.00826446
                    ];

                     
            %temperature and RMS profile over z direction...
            %from figure 7(b)
            porous_hewitt_2_layer.z_T_0_Omega_0p25=[0.261905	0.991736
                0.404762	0.987603
                0.464286	0.971074
                0.470238	0.933884
                0.482143	0.904959
                0.482143	0.855372
                0.464286	0.780992
                0.428571	0.644628
                0.428571	0.53719
                0.488095	0.504132
                0.553571	0.475207
                0.583333	0.404959
                0.553571	0.301653
                0.529762	0.198347
                0.517857	0.119835
                0.529762	0.0330579
                0.654762	0.0123967
                0.77381	0.00826446
                ];
            porous_hewitt_2_layer.z_T_rms_Omega_0p25=[0.039759	1
                0.1144578	0.991701
                0.1493976	0.979253
                0.1578314	0.966805
                0.1385542	0.929461
                0.1156626	0.86722
                0.1036144	0.780083
                0.0987952	0.680498
                0.1024096	0.605809
                0.1156626	0.556017
                0.1349398	0.526971
                0.1542168	0.502075
                0.133735	0.473029
                0.1120482	0.435685
                0.1024096	0.394191
                0.0987952	0.323651
                0.1024096	0.228216
                0.113253	0.145228
                0.1289156	0.0954357
                0.146988	0.0622407
                0.1578314	0.0373444
                0.139759	0.0165975
                0.0891566	0.00414938
                ];
            porous_hewitt_2_layer.z_w_rms_Omega_0p25=[0.0240964	1
                0.0542168	0.987552
                0.0626506	0.983402
                0.0626506	0.937759
                0.0698796	0.883817
                0.0783132	0.804979
                0.0819278	0.73444
                0.0819278	0.66805
                0.0819278	0.626556
                0.073494	0.568465
                0.0518072	0.514523
                0.0518072	0.489627
                0.0674698	0.452282
                0.079518	0.39834
                0.0819278	0.33195
                0.0819278	0.26971
                0.079518	0.20332
                0.0722892	0.136929
                0.0626506	0.0705394
                0.0638554	0.0248963
                0.0421686	0.00414938
                ];
            porous_hewitt_2_layer.z_u_rms_Omega_0p25=[0.0915662	0.983402
                0.0554216	0.93361
                0.0361446	0.887967
                0.0253012	0.842324
                0.01566266	0.763485
                0.01445784	0.680498
                0.01686746	0.60166
                0.026506	0.560166
                0.039759	0.53112
                0.0566266	0.502075
                0.0325302	0.460581
                0.01807228	0.414938
                0.01445784	0.323651
                0.01566266	0.248963
                0.0216868	0.182573
                0.0373494	0.112033
                0.060241	0.0622407
                0.0819278	0.033195
                0.1024096	0.0124481
                ];
     
            %temperature and RMS profile over z direction...
            %from figure 7(c)
            porous_hewitt_2_layer.z_T_0_Omega_1p28=[0.255952	0.991736
                0.345238	0.983471
                0.392857	0.971074
                0.392857	0.942149
                0.392857	0.896694
                0.39881	0.855372
                0.392857	0.797521
                0.375	0.752066
                0.357143	0.702479
                0.333333	0.652893
                0.309524	0.57438
                0.327381	0.545455
                0.392857	0.520661
                0.47619	0.504132
                0.60119	0.487603
                0.678571	0.454545
                0.696429	0.42562
                0.684524	0.380165
                0.666667	0.338843
                0.630952	0.260331
                0.607143	0.144628
                0.613095	0.0950413
                0.613095	0.0619835
                0.613095	0.0330579
                0.672619	0.0165289
                0.761905	0.00826446
                0.880952	0.00826446
                ];
            porous_hewitt_2_layer.z_T_rms_Omega_1p28=[0.0464286	0.995868
                0.1011904	0.991736
                0.1285714	0.983471
                0.145238	0.966942
                0.1357142	0.929752
                0.1214286	0.884298
                0.1130952	0.822314
                0.1071428	0.764463
                0.102381	0.694215
                0.1011904	0.644628
                0.1035714	0.607438
                0.120238	0.553719
                0.1488096	0.528926
                0.1892858	0.512397
                0.1988096	0.504132
                0.197619	0.495868
                0.1809524	0.487603
                0.154762	0.475207
                0.1309524	0.458678
                0.1166666	0.442149
                0.1083334	0.421488
                0.1011904	0.380165
                0.1011904	0.347107
                0.1035714	0.285124
                0.1083334	0.22314
                0.1166666	0.157025
                0.1285714	0.0909091
                0.1369048	0.0661157
                0.145238	0.0371901
                0.127381	0.0165289
                0.1	0.00826446
                0.0666666	0.00826446
                ];
            porous_hewitt_2_layer.z_w_rms_Omega_1p28=[0.0285714	0.991736
                0.0511904	0.979339
                0.0559524	0.946281
                0.054762	0.880165
                0.0571428	0.809917
                0.0607142	0.752066
                0.0630952	0.68595
                0.0619048	0.619835
                0.054762	0.557851
                0.0357142	0.520661
                0.022619	0.5
                0.0404762	0.475207
                0.054762	0.442149
                0.0630952	0.380165
                0.0619048	0.289256
                0.0583334	0.22314
                0.0559524	0.157025
                0.054762	0.0785124
                0.0571428	0.0371901
                0.0428572	0.0123967
                0.020238	0.00413223
                ];
            porous_hewitt_2_layer.z_u_rms_Omega_1p28=[0.0880952	0.983471
                0.0666666	0.946281
                0.0488096	0.896694
                0.0333334	0.834711
                0.020238	0.764463
                0.01785714	0.681818
                0.0261904	0.61157
                0.0428572	0.553719
                0.0571428	0.524793
                0.072619	0.504132
                0.054762	0.471074
                0.0357142	0.429752
                0.0238096	0.371901
                0.01666666	0.297521
                0.0214286	0.22314
                0.0380952	0.144628
                0.05	0.0991736
                0.070238	0.0495868
                0.095238	0.0123967
                ];
            
                 
            %temperature and RMS profile over z direction...
            %from figure 7(d)
            porous_hewitt_2_layer.z_T_0_Omega_10=[0.130952	0.991736
                0.244048	0.975207
                0.285714	0.958678
                0.285714	0.904959
                0.27381	0.809917
                0.27381	0.632231
                0.267857	0.557851
                0.279762	0.53719
                0.339286	0.516529
                0.5	0.5
                0.619048	0.491736
                0.72619	0.466942
                0.732143	0.42562
                0.72619	0.293388
                0.72619	0.194215
                0.720238	0.0867769
                0.720238	0.0454545
                0.797619	0.0206612
                0.928571	0.00826446
                ];
            porous_hewitt_2_layer.z_T_rms_Omega_10=[0.0202381	0.991736
                0.047619	0.983471
                0.0583333	0.971074
                0.0642857	0.958678
                0.0619048	0.917355
                0.0571429	0.876033
                0.0547619	0.822314
                0.0547619	0.768595
                0.0559524	0.665289
                0.0595238	0.599174
                0.0654762	0.541322
                0.0571429	0.516529
                0.0392857	0.5
                0.0571429	0.483471
                0.0666667	0.46281
                0.0630952	0.421488
                0.0583333	0.376033
                0.0559524	0.309917
                0.0559524	0.214876
                0.0571429	0.144628
                0.0607143	0.0909091
                0.0654762	0.0413223
                0.0535714	0.0206612
                0.0297619	0.00826446
                ];
            porous_hewitt_2_layer.z_w_rms_Omega_10=[0.0130952	0.991736
                0.0321429	0.971074
                0.0416667	0.946281
                0.0464286	0.880165
                0.047619	0.793388
                0.047619	0.694215
                0.0464286	0.628099
                0.0428571	0.561983
                0.027381	0.516529
                0.00357143	0.5
                0.0345238	0.471074
                0.0416667	0.446281
                0.0452381	0.376033
                0.047619	0.289256
                0.047619	0.181818
                0.0452381	0.0991736
                0.0404762	0.0495868
                0.0214286	0.0123967
                ];
            porous_hewitt_2_layer.z_u_rms_Omega_10=[0.0345238	0.975207
                0.0166667	0.917355
                0.0119048	0.838843
                0.0107143	0.756198
                0.0119048	0.669421
                0.0154762	0.595041
                0.022619	0.549587
                0.0357143	0.520661
                0.0511905	0.5
                0.0309524	0.466942
                0.0166667	0.409091
                0.0119048	0.322314
                0.0119048	0.227273
                0.0119048	0.140496
                0.0202381	0.0661157
                0.0369048	0.0206612
                ];

            %Here, reorder the data for the following cases..
            %make the first column as z just make consistent with other
            %dataset from Hewitt also encoded here...
            struct_name={'z_T_0_Omega_0p04','z_T_rms_Omega_0p04','z_w_rms_Omega_0p04','z_u_rms_Omega_0p04',...
                'z_T_0_Omega_0p25','z_T_rms_Omega_0p25','z_w_rms_Omega_0p25','z_u_rms_Omega_0p25',...
                'z_T_0_Omega_1p28','z_T_rms_Omega_1p28','z_w_rms_Omega_1p28','z_u_rms_Omega_1p28',...
                'z_T_0_Omega_10','z_T_rms_Omega_10','z_w_rms_Omega_10','z_u_rms_Omega_10'}
            for struct_ind=1:length(struct_name)
               tmp= porous_hewitt_2_layer.(struct_name{struct_ind})(:,1);
               porous_hewitt_2_layer.(struct_name{struct_ind})(:,1)=porous_hewitt_2_layer.(struct_name{struct_ind})(:,2);
               porous_hewitt_2_layer.(struct_name{struct_ind})(:,2)=tmp;
            end
            
        end
    end
end
%%-------------Old code that are repeated....

          
%             if strcmp(obj.flow(1:7),'IFSC_2D')
%                 data{1}.y=obj.z_list/(2*pi/obj.k_opt);
%                 data{2}.y=obj.z_list/(2*pi/obj.k_opt);
%                 data{3}.y=obj.z_list/(2*pi/obj.k_opt);
%                 data{4}.y=obj.z_list/(2*pi/obj.k_opt);
%                 plot_config.label_list={1,'','$z/l_{opt}$'};
%             else


%             syms z;
%             u_laminar=obj.F_sin/obj.ks^2*sin(obj.ks*z)...
%                       +obj.F_sin_2ks/(2*obj.ks)^2*sin(2*obj.ks*z+obj.phase_2ks)...
%                       +obj.F_sin_3ks/(3*obj.ks)^2*sin(3*obj.ks*z+obj.phase_3ks)...
%                       +obj.F_sin_4ks/(4*obj.ks)^2*sin(4*obj.ks*z+obj.phase_4ks);
%             u_laminar_num=double(subs(u_laminar,z,obj.z_list));
%             obj.u=h5read_complex(obj.h5_name,'/tasks/u');
%             obj.u_fluctuation=obj.u;
%             for z_ind=1:length(u_laminar_num)
%                 obj.u_fluctuation(z_ind,:,:)=obj.u(z_ind,:,:)-u_laminar_num(z_ind);
%             end
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
%             S=h5read_complex(obj.h5_name,'/tasks/S');
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
%             S=h5read_complex(obj.h5_name,'/tasks/S');
%             w=h5read_complex(obj.h5_name,'/tasks/w');
%             obj.wS=w.*S;
%             data{1}.z=squeeze(mean(obj.wS,2));
% %             plot_config.label_list={1,'$t$','$z/l_{opt}$'};
%             plot_config.colormap='bluewhitered';
%             plot_config.print_size=[1,1200,1200];
%             plot_config.print=obj.print;
%             plot_config.name=[obj.h5_name(1:end-3),'_wS_x_ave.png'];
%             plot_contour(data,plot_config);
%         end



%         function obj=T_total_xt_ave(obj)
% %             data{1}.y=obj.z_list/(2*pi/obj.k_opt);
%             T=h5read_complex(obj.h5_name,'/tasks/T');
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
%             u=h5read_complex(obj.h5_name,'/tasks/u');
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
%             obj.w=h5read_complex(obj.h5_name,'/tasks/w');
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
        
        



%         function obj=E_T_time(obj,elevator_growth_rate)
%             %%Plot the salinity potential energy as a function over time
%             
%             
%             if nargin<2 || isempty(elevator_growth_rate)
%                 elevator_growth_rate=0;    
%                 %flag.mean='laminar_cou';   %%default value of flag_mean if not given, just set the laminar  flow.
%                     %error('The flag_mean is missing.')
%             end
%             obj.T=h5read_complex(obj.h5_name,'/tasks/T');
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


