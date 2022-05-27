classdef my_default
    %MY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        no_ylabel=0;
        folder_name=0;
        plot_config=struct;
        Nu=0;
        Nu_S=0;
        obj=struct;
        par_num=4;
%         obj.point_list=10000;
        bif_type={'bpt'};
        depth=1;
        
        %%These are parameter for the Kolmogorov shear flow
        F_sin=0;
        ks=2*pi;
        
        z_basis_mode='Chebyshev';
        
        %background large scale shear flow and their derivative
        U_bg=0;
        d_U_bg=0;
        dd_U_bg=0;
        
        root_folder_name=0;
        bpt_name=0;
        store_folder_name_suffix='';

        %In default, do not detech the bifurcaiton point and do not
        %compute the spectrum... They are quite computational
        %expensive..
        bifcheck=0;
        spcalc=0;
        tol=1e-8;
        tol_first=1e-8;
        tol_list=1e-8;
        ds=0.01;
        
        ilam_phase=10; %The index for the phase speed.. used for flow with phase condition
        ilam_phase_U_r=14; %The index for the phase condition of transition in the vertical direction
        ilam_phase_T_r=15;
        ilam_phase_S_r=16;
        ilam_phase_z=17;
        
        phase_tr=0;
        
        point_list=0; %The point list that you want specific point (lam)
        dsmax=1;
        dsmin=10^(-10);
        ilam=1;
        lammax=10^6;
        nfloq=0;
        neig=40;
        eigref=0;
        foldcheck=0;
        no_phase_step=1;
        sec_bif_stop=0;
        rotation_z=1;
        rotation_y=1;
        variable_version='complex_zonal_with_omega_z';
%             {'real_minimal',... only have five  variable, w_hat, T_hat, S_hat, T_0, S_0
%             'complex_minimal',... only have 8 variable, add imag of the hat, w_hat, T_hat, S_hat, T_0, S_0, w_hat_imag, T_hat_imag, S_hat_imag
%             'complex_zonal',... add one variable for the zonal flow, that should be the 9th variable
%             'complex_zonal_with_omega_z'... add two more variable, that should be 11 variable, (most comprehensive one), w_hat, T_hat, S_hat, T_0, S_0, w_hat_imag, T_hat_imag, S_hat_imag, U_0, omega_z_hat, omega_z_hat_imag
%             'real_with_omega_z',...A short version to test the hexagon planform, %This should be size variable, w_hat, T_hat, S_hat, T_0, S_0, omega_z_hat
%             'complex_with_omega_z'...A short version to test the hexagon planform, but  with the complex one. This should have 10 variable, w_hat, T_hat, S_hat, T_0, S_0, w_hat_imag, T_hat_
%             };
            %kx='scaling';%%{'scaling','equal','value',}
            %ky='scaling';%%{'scaling','equal','value',}
        kx_square=100;
        ilam_new=[];
        lam_new=[];
        bpt1_pt0=0.1;
        jac=1;
        resfac=10^(-3);
        
        z_bc_T_left='dirichlet'
        z_bc_T_right='dirichlet'
        z_bc_S_left='dirichlet'
        z_bc_S_right='dirichlet'
        z_bc_w_left='dirichlet'
        z_bc_w_right='dirichlet'
        z_bc_u_v_left='dirichlet'
        z_bc_u_v_right='dirichlet'
    
        grid=[0,1,64];%The grid points, that should be a n*3 matrix
        %for a compound grid, we can have something like
        %grid=[0,1/3,16;
%                 1/3,1/2,16;
%                 1/2,2/3,16;
%                 2/3,1,16]; This will have the benefit to cluster some
%                 grid points near the sharp interface... 
        mat=[];
        mass_matrix_Laplacian=0;
        compound=0;
        
        bpt_list=1; %the bpt point that you want to do the bifurcation analysis... 
        bpt_root_list={'tr/bpt1'}; %the root of secondary bifurcation
        bifcheck_list=0;
        spcalc_list=0;
        
        parameters;
        
        intol=0;
    end
    
    
    
    methods
        function obj = my_default(obj)
            obj.no_ylabel=0;
            obj.folder_name=0;
            obj.plot_config=struct;
            obj.plot_config.visible=0;
            obj.plot_config.print=0;
            obj.plot_config.post=0;
            obj.Nu=0;
            obj.Nu_S=0;
            obj.obj=struct;
            obj.par_num=4;
            obj.depth=1;
            obj.bif_type={'bpt'};
            %obj.point_list=10000;
            
            %MY Construct an instance of this class
            %This is the flag used for switching branch.
            obj.root_folder_name=0;
            obj.bpt_name=0;
            obj.store_folder_name_suffix='';
            
            %In default, do not detech the bifurcaiton point and do not
            %compute the spectrum... They are quite computational
            %expensive..
            obj.bifcheck=0;
            obj.spcalc=0;
            obj.tol=1e-8;
            obj.tol_first=1e-8;
            obj.ds=0.01;
            obj.ilam_phase=10; %The index for the phase speed.. used for flow with phase condition
            obj.ilam_phase_U_r=14; %The index for the phase condition of transition in the vertical direction
            obj.ilam_phase_T_r=15;
            obj.ilam_phase_S_r=16;
            obj.ilam_phase_z=17;

            obj.phase_tr=0;
            
            obj.point_list=0; %The point list that you want specific point (lam)
            obj.dsmax=1;
            obj.ilam=1;
            obj.bif_type='steady'; %This can be steady, hopf, fold or tw (traveling wave)
            obj.lammax=10^6;
            obj.nfloq=0;
            
            obj.neig=40;
            obj.eigref=0;
            obj.foldcheck=0;
            obj.no_phase_step=1;
            obj.sec_bif_stop=0;
            obj.variable_version='complex_zonal_with_omega_z';
            obj.bpt1_pt0=0.1;
            obj.jac=1;
            obj.resfac=10^(-3);
            
            obj.z_bc_T_left='dirichlet';
            obj.z_bc_T_right='dirichlet';
            obj.z_bc_S_left='dirichlet';
            obj.z_bc_S_right='dirichlet';
            obj.z_bc_w_left='dirichlet';
            obj.z_bc_w_right='dirichlet';
            obj.z_bc_u_v_left='dirichlet';
            obj.z_bc_u_v_right='dirichlet';
            obj.mass_matrix_Laplacian=0;
            obj.compound=0;

            obj.intol=0;
        end
        
        function obj = salt_finger_init_tr(obj,folder_name)

                Ta_sqrt_z=0;
                Ta_sqrt_y=0;
                A_U_bg=0;
                C=0;
                dy_T_mean=1;
                dy_S_mean=1;
                c_phase_kx=0;%phase velocity for horizontal translation invariance
                U_r=0; %reference large-scale shear
                T_r=0; %refernce temperature
                S_r=0; %reference salinity
                c_phase_z=0;%phase velocity for vertical translation invariance (periodic B.C. in vertical)
                obj.mass_matrix_Laplacian=0;
                Lx=1;
                switch folder_name
                    case {'salt_finger_Ra_S2T_3D_no_slip'}
                        kx=-10^5; ky=-10^5; 
                        Pr=7; nx=512; tau=0.01; Ra_T=10^5; Ra_S2T=1090;
                        obj.ilam=4; obj.lammax=Ra_T; Ta_sqrt_z=0; Ta_sqrt_y=0;
                        obj.point_list=[1/90,1/85,1/80,1/75,1/70,1/65,1/60,...
                              1/55,1/50,1/45,1/40,1/35,1/30,1/25,1/20,1/15,...
                              0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1/1.2,0.9,1,...
                              1.2,2,3,4,5,6,7,8,9,10]*Ra_T;
                        obj.dsmax=Ra_T;
                    case {'salt_finger_Ra_S2T_3D_stress_free'}
                        kx=-10^5; ky=-10^5; 
                        Pr=7; nx=512; tau=0.01; Ra_T=10^5; Ra_S2T=1050;
                        obj.ilam=4; obj.lammax=Ra_T; Ta_sqrt_z=0; Ta_sqrt_y=0;
                        obj.point_list=[1/90,1/85,1/80,1/75,1/70,1/65,1/60,...
                              1/55,1/50,1/45,1/40,1/35,1/30,1/25,1/20,1/15,...
                              0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1/1.2,0.9,1,...
                              1.2,2,3,4,5,6,7,8,9,10]*Ra_T;
                        obj.dsmax=Ra_T;
                        obj.z_bc_u_v_left='neumann';
                        obj.z_bc_u_v_right='neumann';
                        %obj.point_list=[0.1:0.1:1,1.2,2:10]*Ra_T;
                     case {'salt_finger_Ra_S2T_3D_stress_free'}
                        kx=-10^5; ky=-10^5; 
                        Pr=7; nx=512; tau=0.01; Ra_T=10^5; Ra_S2T=900;
                        obj.ilam=4; obj.lammax=Ra_T; Ta_sqrt_z=0; Ta_sqrt_y=0;
                        obj.point_list=[1/95,1/90,1/85,1/80,1/75,1/70,1/65,1/60,...
                              1/55,1/50,1/45,1/40,1/35,1/30,1/25,1/20,1/15,...
                              0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1/1.2,0.9,1,...
                              1.2,2,3,4,5,6,7,8,9,10]*Ra_T;
                        obj.dsmax=Ra_T;
                        %obj.point_list=[0.1:0.1:1,1.2,2:10]*Ra_T;
                        obj.z_bc_u_v_left='neumann';
                        obj.z_bc_u_v_right='neumann';
                    case 'salt_finger_Ra_S2T_Nu'
                        kx=-10^5; ky=-10^5; 
                        Pr=7; nx=512; tau=0.01; Ra_T=10^5; Ra_S2T=900;
                        obj.ilam=4; obj.lammax=Ra_T;
                        obj.point_list=10*Ra_T; Ta_sqrt_z=0; Ta_sqrt_y=0;
                        obj.dsmax=Ra_T;
                    case 'salt_finger_Ra_S2T_IC'
                        kx=-10^5; ky=-10^5; 
                        Pr=7; nx=512; tau=0.01; Ra_T=10^5; Ra_S2T=900; Ta_sqrt_z=0; Ta_sqrt_y=0;
                        obj.ilam=4; obj.lammax=Ra_T; 
                        obj.point_list=[1300:100:4000];
                        obj.dsmax=Ra_T;

                    case {'salt_finger_Ra_S2T_2D'}
                        kx=-10^5; ky=0; 
                        Pr=7; nx=512, tau=0.01; Ra_T=10^5; Ra_S2T=900; Ta_sqrt_z=0; Ta_sqrt_y=0;
                        obj.point_list=[1/95,1/90,1/85,1/80,1/75,1/70,1/65,1/60,...
                              1/55,1/50,1/45,1/40,1/35,1/30,1/25,1/20,1/15,...
                              0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1/1.2,0.9,1,...
                              1.2,2,3,4,5,6,7,8,9,10]*Ra_T;
                        %obj.point_list=[0.1:0.1:1,1.2,2:10]*Ra_T;
                        obj.ilam=4; obj.lammax=Ra_T; obj.dsmax=Ra_T;
                    case 'salt_finger_Ra_S2T_2D_Iw'
                        kx=-10^5; ky=0;
                        Pr=7; nx=512, tau=0.01; Ra_T=10^5; Ra_S2T=900; Ta_sqrt_z=0; Ta_sqrt_y=0;
                        obj.ilam=4; obj.lammax=Ra_T;
                    case 'salt_finger_Ra_S2T_low_Pr'
                        kx=-10^5; ky=-10^5;
                        Pr=0.03; nx=512; tau=0.01; Ra_T=10^5; Ra_S2T=900; Ta_sqrt_z=0; Ta_sqrt_y=0;
                        obj.ilam=4; obj.lammax=Ra_T; obj.dsmax=Ra_T;
                    case 'salt_finger_Ra_S2T_high_Pr'
                        kx=-10^5; ky=-10^5;
                        Pr=Inf; nx=512; tau=0.01; Ra_T=10^5; Ra_S2T=900; Ta_sqrt_z=0; Ta_sqrt_y=0;
                        obj.ilam=4; obj.lammax=Ra_T;
                    case {'salt_finger_Ra_S2T_hopf'}
                        kx=-10^5; ky=-10^5; 
                        Pr=7; nx=128; tau=1/3; Ra_T=10^5; Ra_S2T=Ra_T*(tau-0.005); Ta_sqrt_z=0; Ta_sqrt_y=0;
                        obj.ilam=4; obj.lammax=Ra_T; obj.dsmax=Ra_T;
                        obj.point_list=(0.4:0.1:1)*Ra_T; 
                    case 'salt_finger_Ra_S2T_tw'
                        kx=-10^5; ky=-10^5; 
                        Pr=7; nx=128; tau=1/3; Ra_T=10^5; Ra_S2T=Ra_T*(tau-0.005); Ta_sqrt_z=0; Ta_sqrt_y=0;
                        obj.ilam=4; obj.lammax=Ra_T; obj.dsmax=Ra_T;
                        obj.point_list=(0.4:0.1:1)*Ra_T; 
                    case {'salt_finger_Ra_S2T_hopf_2D'}
                        kx=-10^5; ky=0; 
                        Pr=7; nx=128; tau=1/3; Ra_T=10^5; Ra_S2T=Ra_T*(tau-0.005); Ta_sqrt_z=0; Ta_sqrt_y=0;
                        obj.ilam=4; obj.lammax=Ra_T; obj.dsmax=Ra_T;
                        obj.point_list=(0.4:0.1:1)*Ra_T;
                    case 'salt_finger_Ra_S2T_hopf_2D_Iw'
                        kx=-10^5; ky=0; 
                        Pr=7; nx=128; tau=1/3; Ra_T=10^5; Ra_S2T=Ra_T*(tau-0.005); Ta_sqrt_z=0; Ta_sqrt_y=0;
                        obj.ilam=4; obj.lammax=Ra_T;
                        obj.dsmax=Ra_T;
                        obj.point_list=(0.4:0.1:1)*Ra_T;
                    case 'salt_finger_Ra_S2T_hopf_tau0p01' 
                        kx=-10^5; ky=-10^5; 
                        Pr=7; nx=128; tau=0.01; Ra_T=10^5; Ra_S2T=Ra_T*(tau-0.005); Ta_sqrt_z=0; Ta_sqrt_y=0;
                        obj.ilam=4; obj.lammax=Ra_T; obj.dsmax=Ra_T;
                        obj.point_list=(0.4:0.1:1)*Ra_T; 
                    case 'salt_finger_Ra_S2T_high_Ra'
                        kx=-10^5; ky=-10^5; 
                        Pr=7; nx=512; tau=1/3; Ra_T=3*10^10; Ra_S2T=1.007*10^10; Ta_sqrt_z=0; Ta_sqrt_y=0;
                        obj.ilam=4; obj.lammax=Ra_T; obj.dsmax=Ra_T;
                        obj.point_list=Ra_T;
                    case {'salt_finger_tau'}
                        kx=-10^5; ky=-10^5; 
                        Pr=7; R_rho_T2S=1.2; Ra_T=10^5; nx=512; tau=1; Ra_S2T=Ra_T/R_rho_T2S; Ta_sqrt_z=0; Ta_sqrt_y=0;
                        obj.point_list=1./[0.5,1/3,0.33:-0.02:0.01];
                        obj.ilam=5; obj.lammax=100; obj.dsmax=10;
                    case 'salt_finger_tau_Nu'
                        kx=-10^5; ky=-10^5; 
                        Pr=7; R_rho_T2S=1.2; Ra_T=10^5; nx=512; tau=1; Ra_S2T=Ra_T/R_rho_T2S; Ta_sqrt_z=0; Ta_sqrt_y=0;
                        obj.point_list=100;
                        obj.ilam=5; obj.lammax=100; obj.dsmax=10;
                    case 'salt_finger_tau_low_Ra_S2T_Nu'
                        kx=-10^5; ky=-10^5; 
                        Pr=7; R_rho_T2S=20; Ra_T=10^5; nx=512; tau=1; Ra_S2T=Ra_T/R_rho_T2S; Ta_sqrt_z=0; Ta_sqrt_y=0;
                        obj.point_list=100;
                        obj.ilam=5; obj.lammax=100; obj.dsmax=10;
                    case 'salt_finger_tau_high_Ra_S2T'
                        kx=-10^5; ky=-10^5; 
                        Pr=7; R_rho_T2S=1.2; Ra_T=3*10^10; nx=512; tau=1; Ra_S2T=Ra_T/R_rho_T2S; Ta_sqrt_z=0; Ta_sqrt_y=0;
                        obj.point_list=3;
                        obj.ilam=5; obj.lammax=3; obj.dsmax=1;
                    case 'salt_finger_tau_low_Pr'
                        kx=-10^5; ky=-10^5; 
                        Pr=0.03; R_rho_T2S=1.2; Ra_T=10^5; nx=512; tau=1; Ra_S2T=Ra_T/R_rho_T2S; Ta_sqrt_z=0; Ta_sqrt_y=0;
                        obj.point_list=1./[0.5,1/3,0.33:-0.02:0.01];
                        obj.ilam=5; obj.lammax=100; obj.dsmax=10;
                    case 'salt_finger_tau_high_Pr'
                        kx=-10^5; ky=-10^5; 
                        Pr=100; R_rho_T2S=1.2; Ra_T=10^5; nx=512; tau=1; Ra_S2T=Ra_T/R_rho_T2S; Ta_sqrt_z=0; Ta_sqrt_y=0;
                        obj.point_list=1./[0.5,1/3,0.33:-0.02:0.01];
                        obj.ilam=5; obj.lammax=100; obj.dsmax=10;
                   case 'salt_finger_kx_low_Ra_S2T_2D_no_slip'
                        kx=-20; ky=0; %ky corresponding to 3D setup
                        Pr=7; tau=0.01; 
                        Ra_T=10^5; 
                        Ra_S2T=Ra_T/40; nx=128; Ta_sqrt_z=0; Ta_sqrt_y=0;
                        obj.ilam=1; obj.lammax=-0.01; obj.dsmax=3;
                        obj.ds=0.01;
                        obj.dsmax=1; obj.dsmin=1e-10;
                        obj.point_list=[-19,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6.873,-6,-5,-4,-3,-2,-1,-0.5,-0.01];
                    case 'salt_finger_kx_low_Ra_S2T_3D_no_slip'
                        kx=-15; ky=-10^4; %ky corresponding to 3D setup
                        Pr=7; tau=0.01; 
                        Ra_T=10^5; 
                        Ra_S2T=Ra_T/40; nx=128; Ta_sqrt_z=0; Ta_sqrt_y=0;
                        obj.ilam=1; obj.lammax=-0.01; obj.dsmax=3;
                        obj.ds=0.01; obj.dsmax=1; obj.dsmin=1e-10;
                        %obj.point_list=[-19,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6.873,-6,-5,-4,-3,-2,-1,-0.01];
                        obj.point_list=[-13,-12,-11,-10,-9,-8,-6.873,-6,-5,-4,-3,-2,-1,-0.5,-0.01];

                    case 'salt_finger_kx_low_Ra_S2T_low_Pr_2D_no_slip'
                        kx=-20; ky=0; %ky corresponding to 3D setup
                        Pr=0.05; tau=0.01; 
                        Ra_T=10^5; Ta_sqrt_z=0; Ta_sqrt_y=0;
                        Ra_S2T=Ra_T/40; nx=128; 
                        obj.ilam=1; obj.lammax=-0.01; obj.dsmax=3;
                        obj.ds=0.01; obj.dsmax=1; obj.dsmin=1e-10;
                        obj.point_list=[-19,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-6.873,-6,-5,-4,-3,-2,-1,-0.5,-0.01];
                        %obj.point_list=[-fliplr(obj.point_list),25];
                        %                 obj.point_list=-0.01;
                        %obj.point_list=[-13,-12,-11,-10,-9,-8,-7,-6.873,-6,-5,-4,-3,-2,-1,-0.01];
                        %obj.point_list=[obj.point_list,sqrt(2)*obj.point_list];
                        %obj.point_list=sort(obj.point_list,'ascend');
                    case 'salt_finger_kx_low_Ra_S2T_low_Pr_3D_no_slip'
                        kx=-15; ky=-10^4; %ky corresponding to 3D setup
                        Pr=0.05; tau=0.01; 
                        Ra_T=10^5; 
                        Ra_S2T=Ra_T/40; nx=128; Ta_sqrt_z=0; Ta_sqrt_y=0;
                        obj.ilam=1; obj.lammax=-0.01; obj.dsmax=3;
                        obj.ds=0.01; obj.dsmax=1; obj.dsmin=1e-10;
                        obj.point_list=[-13,-12,-11,-10,-9,-8,-6.873,-6,-5,-4,-3,-2,-1,-0.5,-0.01];
                    
                    case 'salt_finger_kx_low_Ra_S2T_high_Pr_2D_no_slip'
                        kx=-20; ky=0; %ky corresponding to 3D setup
                        Pr=Inf; tau=0.01; 
                        Ra_T=10^5; Ta_sqrt_z=0; Ta_sqrt_y=0;
                        Ra_S2T=Ra_T/40; nx=128; 
                        obj.ilam=1; obj.lammax=-0.01; obj.dsmax=3;
                        obj.ds=0.01; obj.dsmax=1; obj.dsmin=1e-10;
                        obj.point_list=[-19,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-6.873,-6,-5,-4,-3,-2,-1,-0.5,-0.01];
                       
                    case 'salt_finger_kx_low_Ra_S2T_high_Pr_3D_no_slip'
                        kx=-15; ky=-10^4; %ky corresponding to 3D setup
                        Pr=Inf; tau=0.01; 
                        Ra_T=10^5; 
                        Ra_S2T=Ra_T/40; nx=128; Ta_sqrt_z=0; Ta_sqrt_y=0;
                        obj.ilam=1; obj.lammax=-0.01; obj.dsmax=3;
                        obj.ds=0.01; obj.dsmax=1; obj.dsmin=1e-10;
                        obj.point_list=[-13,-12,-11,-10,-9,-8,-6.873,-6,-5,-4,-3,-2,-1,-0.5,-0.01];
                    
                    case 'salt_finger_kx_mid_Ra_S2T_2D_no_slip'
                        kx=-50; ky=0; %ky corresponding to 3D setup
                        Pr=7; tau=0.01; 
                        Ra_T=10^5; 
                        Ra_S2T=Ra_T/2; nx=512; Ta_sqrt_z=0; Ta_sqrt_y=0;
                        obj.ilam=1; obj.lammax=-0.01; obj.dsmax=1;
                        obj.ds=0.01; obj.dsmax=1; obj.dsmin=1e-10;
                        obj.point_list=[-40,-37.5,-35,-32.5,-30,-27.5,-25,-22.5,-20,-17.5,-15,-12.5,-9.7815,-5,-3,-0.01];
                    case 'salt_finger_kx_mid_Ra_S2T_3D_no_slip'
                        kx=-35; ky=-10^4; %ky corresponding to 3D setup
                        Pr=7; tau=0.01; 
                        Ra_T=10^5; 
                        Ra_S2T=Ra_T/2; nx=512; Ta_sqrt_z=0; Ta_sqrt_y=0;
                        obj.ilam=1; obj.lammax=-0.01; obj.dsmax=1;
                        obj.ds=0.01; obj.dsmax=1; obj.dsmin=1e-10;
                        obj.point_list=[-25,-22.5,-20,-17.5,-15,-12.5,-9.7815,-5,-3,-0.01];
                    case 'salt_finger_kx_mid_Ra_S2T_3D_rotation_no_slip'
                        kx=-35; ky=-10^4; %ky corresponding to 3D setup
                        Pr=7; tau=0.01;
                        Ra_T=10^5; 
                        Ra_S2T=Ra_T/2; nx=512; Ta_sqrt_z=-10000; Ta_sqrt_y=0;
                        obj.ilam=1; obj.lammax=-0.01; obj.dsmax=1;
                        obj.ds=0.01; obj.dsmax=1; obj.dsmin=1e-10;
                        obj.point_list=[-25,-20,-15,-9.7815,-5,-0.01];


                    case 'salt_finger_kx_low_Ra_S2T_zero_Pr_3D_no_slip'
                        kx=-15; ky=-10^4; %ky corresponding to 3D setup
                        Pr=0.01; tau=0.01; 
                        Ra_T=10^5; 
                        Ra_S2T=10^5/40; nx=128; Ta_sqrt_z=0; Ta_sqrt_y=0; 
                        obj.rotation_z=1; obj.rotation_y='z^2';
                        obj.F_sin=1; A_U_bg=0;
                        obj.ks=2*pi;
                        obj.ilam=1; obj.lammax=-0.01; obj.dsmax=10^6;
                        obj.ds=0.01; obj.dsmax=1; obj.dsmin=1e-10;
                        C=0;                
        %                 obj.point_list=[-13,-10,-6.873,-3,-0.01];
                        obj.point_list=[-6.873,-3,-0.01];
                    case 'salt_finger_kx_low_Ra_S2T_3D_vorticity'
                        kx=-15; ky=-15; %ky corresponding to 3D setup
                        Pr=0.03; tau=0.01; 
                        Ra_T=10^5; 
                        Ra_S2T=10^5/40; nx=128; Ta_sqrt_z=0; Ta_sqrt_y=0; 
                        obj.rotation_z=1; obj.rotation_y='z^2';
                        obj.F_sin=1; A_U_bg=0;
                        obj.ks=2*pi;
                        obj.ilam=1; obj.lammax=-0.01; obj.dsmax=10^6;
                        obj.ds=0.01; obj.dsmax=1; obj.dsmin=1e-10;
                        C=0;
        %                 C=1/sqrt(6);

        %                 obj.variable_version='real_with_omega_z';
        %                 obj.point_list=-0.01;
                        obj.point_list=[-13,-12,-11,-10,-9,-8,-7,-6.873,-6,-5,-4,-3,-2,-1,-0.5,-0.01];
        %                 obj.point_list=[-13,-10,-6.873,-3,-0.01];
        %                 obj.point_list=8000;

                    case 'salt_finger_kx_low_Ra_S2T_2D_stress_free'
                        kx=-20; ky=0; %ky corresponding to 3D setup
                        Pr=7; tau=0.01; 
                        Ra_T=10^5; 
                        Ra_S2T=Ra_T/40; nx=128; Ta_sqrt_z=0; Ta_sqrt_y=0;
                        obj.ilam=1; obj.lammax=-0.01; obj.dsmax=3;
                        obj.ds=0.01;
                        obj.dsmax=1; obj.dsmin=1e-10;
                        obj.point_list=[-19,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6.873,-6,-5,-4,-3,-2,-1,-0.5,-0.01];
                        obj.z_bc_u_v_left='neumann';
                        obj.z_bc_u_v_right='neumann';
                    case 'salt_finger_kx_low_Ra_S2T_3D_stress_free'
                        kx=-15; ky=-10^4; %ky corresponding to 3D setup
                        Pr=7; tau=0.01; 
                        Ra_T=10^5; 
                        Ra_S2T=Ra_T/40; nx=128; Ta_sqrt_z=0; Ta_sqrt_y=0;
                        obj.ilam=1; obj.lammax=-0.01; obj.dsmax=3;
                        obj.ds=0.01; obj.dsmax=1; obj.dsmin=1e-10;
                        %obj.point_list=[-19,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6.873,-6,-5,-4,-3,-2,-1,-0.01];
                        obj.point_list=[-13,-12,-11,-10,-9,-8,-6.873,-6,-5,-4,-3,-2,-1,-0.5,-0.01];
                        obj.z_bc_u_v_left='neumann';
                        obj.z_bc_u_v_right='neumann';
                    case 'salt_finger_kx_low_Ra_S2T_low_Pr_2D_stress_free'
                        kx=-25; ky=0; %ky corresponding to 3D setup
                        Pr=0.05; tau=0.01; 
                        Ra_T=10^5; Ta_sqrt_z=0; Ta_sqrt_y=0;
                        Ra_S2T=Ra_T/40; nx=128; 
                        obj.ilam=1; obj.lammax=-0.01; obj.dsmax=3;
                        obj.ds=0.01; obj.dsmax=1; obj.dsmin=1e-10;
                        obj.point_list=[-19,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-6.873,-6,-5,-4,-3,-2,-1,-0.5,-0.01];
                        obj.z_bc_u_v_left='neumann';
                        obj.z_bc_u_v_right='neumann';
                    case 'salt_finger_kx_low_Ra_S2T_low_Pr_3D_stress_free'
                        kx=-15; ky=-10^4; %ky corresponding to 3D setup
                        Pr=0.05; tau=0.01; 
                        Ra_T=10^5; 
                        Ra_S2T=Ra_T/40; nx=128; Ta_sqrt_z=0; Ta_sqrt_y=0;
                        obj.ilam=1; obj.lammax=-0.01; obj.dsmax=3;
                        obj.ds=0.01; obj.dsmax=1; obj.dsmin=1e-10;
                        obj.point_list=[-13,-12,-11,-10,-9,-8,-6.873,-6,-5,-4,-3,-2,-1,-0.5,-0.01];
                        obj.z_bc_u_v_left='neumann';
                        obj.z_bc_u_v_right='neumann';

                    case 'salt_finger_kx_mid_Ra_S2T_2D_stress_free'
                        kx=-50; ky=0; %ky corresponding to 3D setup
                        Pr=7; tau=0.01; 
                        Ra_T=10^5; 
                        Ra_S2T=Ra_T/2; nx=512; Ta_sqrt_z=0; Ta_sqrt_y=0;
                        obj.ilam=1; obj.lammax=-0.01; obj.dsmax=1;
                        obj.ds=0.01; obj.dsmax=1; obj.dsmin=1e-10;
                        obj.point_list=[-40,-37.5,-35,-32.5,-30,-27.5,-25,-22.5,-20,-17.5,-15,-12.5,-9.7815,-5,-3,-0.01];
                        obj.z_bc_u_v_left='neumann';
                        obj.z_bc_u_v_right='neumann';
                    case 'salt_finger_kx_mid_Ra_S2T_3D_stress_free'
                        kx=-35; ky=-10^4; %ky corresponding to 3D setup
                        Pr=7; tau=0.01; 
                        Ra_T=10^5; 
                        Ra_S2T=Ra_T/2; nx=512; Ta_sqrt_z=0; Ta_sqrt_y=0;
                        obj.ilam=1; obj.lammax=-0.01; obj.dsmax=1;
                        obj.ds=0.01; obj.dsmax=1; obj.dsmin=1e-10;
                        obj.point_list=[-25,-22.5,-20,-17.5,-15,-12.5,-9.7815,-5,-3,-0.01];
                        obj.z_bc_u_v_left='neumann';
                        obj.z_bc_u_v_right='neumann';
                        
                    case 'salt_finger_kx_low_Ra_S2T_2D_periodic'
                        kx=-20; ky=0; %ky corresponding to 3D setup
                        Pr=7; tau=0.01; 
                        Ra_T=10^5; Lx=1; 
                        Ra_S2T=Ra_T/40; nx=128; Ta_sqrt_z=0; Ta_sqrt_y=0;
                        obj.ilam=1; obj.lammax=-0.01; obj.dsmax=3;
                        obj.ds=0.01; obj.intol=0;
                        obj.dsmax=0.3; obj.dsmin=1e-10;
                        obj.point_list=[-19,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6.873,-6,-5,-4,-3,-2,-1,-0.5,-0.01];
                        obj.z_basis_mode='Fourier';
                        obj.z_bc_u_v_left='periodic'; obj.z_bc_u_v_right='periodic';
                        obj.z_bc_w_left='periodic'; obj.z_bc_w_right='periodic';
                        obj.z_bc_T_left='periodic'; obj.z_bc_T_right='periodic';
                        obj.z_bc_S_left='periodic'; obj.z_bc_S_right='periodic';
                        obj.ks=2*pi/Lx; obj.F_sin=1;
                        A_U_bg=0;
                    case 'salt_finger_kx_low_Ra_S2T_3D_periodic'
                        kx=-15; ky=-10^4; %ky corresponding to 3D setup
                        Pr=7; tau=0.01; 
                        Ra_T=10^5; 
                        Ra_S2T=Ra_T/40; nx=128; Ta_sqrt_z=0; Ta_sqrt_y=0;
                        obj.ilam=1; obj.lammax=-0.01; obj.dsmax=3;
                        obj.ds=0.01; obj.dsmax=1; obj.dsmin=1e-10;
                        %obj.point_list=[-19,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6.873,-6,-5,-4,-3,-2,-1,-0.01];
                        obj.point_list=[-13,-12,-11,-10,-9,-8,-6.873,-6,-5,-4,-3,-2,-1,-0.5-0.01];
                        obj.z_basis_mode='Fourier';
                        obj.z_bc_u_v_left='periodic'; obj.z_bc_u_v_right='periodic';
                        obj.z_bc_w_left='periodic'; obj.z_bc_w_right='periodic';
                        obj.z_bc_T_left='periodic'; obj.z_bc_T_right='periodic';
                        obj.z_bc_S_left='periodic'; obj.z_bc_S_right='periodic';
                        obj.ks=2*pi; obj.F_sin=1;
                        A_U_bg=0.0001;
                    case 'salt_finger_kx_low_Ra_S2T_low_Pr_2D_periodic'
                        kx=-25; ky=0; %ky corresponding to 3D setup
                        Pr=0.05; tau=0.01; 
                        Ra_T=10^5; Ta_sqrt_z=0; Ta_sqrt_y=0;
                        Ra_S2T=Ra_T/40; nx=128; 
                        obj.ilam=1; obj.lammax=-0.01; obj.dsmax=3;
                        obj.ds=0.01; obj.dsmax=1; obj.dsmin=1e-10;
                        obj.point_list=[-19,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-6.873,-6,-5,-4,-3,-2,-1,-0.5,-0.01];
                        obj.z_basis_mode='Fourier';
                        obj.z_bc_u_v_left='periodic'; obj.z_bc_u_v_right='periodic';
                        obj.z_bc_w_left='periodic'; obj.z_bc_w_right='periodic';
                        obj.z_bc_T_left='periodic'; obj.z_bc_T_right='periodic';
                        obj.z_bc_S_left='periodic'; obj.z_bc_S_right='periodic';
                        obj.ks=2*pi; obj.F_sin=1;
                        A_U_bg=0.0001;
                    case 'salt_finger_kx_low_Ra_S2T_low_Pr_3D_periodic'
                        kx=-15; ky=-10^4; %ky corresponding to 3D setup
                        Pr=0.05; tau=0.01; 
                        Ra_T=10^5; 
                        Ra_S2T=Ra_T/40; nx=128; Ta_sqrt_z=0; Ta_sqrt_y=0;
                        obj.ilam=1; obj.lammax=-0.01; obj.dsmax=3;
                        obj.ds=0.01; obj.dsmax=1; obj.dsmin=1e-10;
                        obj.point_list=[-13,-12,-11,-10,-9,-8,-6.873,-6,-5,-4,-3,-2,-1,-0.5,-0.01];
                        obj.z_basis_mode='Fourier';
                        obj.z_bc_u_v_left='periodic'; obj.z_bc_u_v_right='periodic';
                        obj.z_bc_w_left='periodic'; obj.z_bc_w_right='periodic';
                        obj.z_bc_T_left='periodic'; obj.z_bc_T_right='periodic';
                        obj.z_bc_S_left='periodic'; obj.z_bc_S_right='periodic';
                        obj.ks=2*pi; obj.F_sin=1;
                        A_U_bg=0.0001;
                    case 'salt_finger_kx_mid_Ra_S2T_3D_rotation'
                        kx=-35; ky=-10^4; %ky corresponding to 3D setup
                        Pr=7; tau=0.01;
                        Ra_T=10^5; 
                        Ra_S2T=Ra_T/2; nx=512; Ta_sqrt_z=-10000; Ta_sqrt_y=0;
                        obj.ilam=1; obj.lammax=-0.01; obj.dsmax=1;
                        obj.ds=0.01; obj.dsmax=1; obj.dsmin=1e-10;
                        obj.point_list=[-25,-20,-15,-9.7815,-5,-0.01];

                    case 'salt_finger_kx_high_Ra_S2T_3D'
                        kx=-500; ky=-10^4; %ky corresponding to 3D setup
                        Pr=7; tau=1/3; 
                        Ra_T=3.6*10^10; 
                        Ra_S2T=Ra_T/1.2; nx=512; Ta_sqrt_z=0; Ta_sqrt_y=0;
                        obj.ilam=1; obj.lammax=-0.01; obj.dsmax=30;
                        obj.ds=0.01; obj.dsmax=1; obj.dsmin=1e-10;
                        obj.point_list=-0.01;
                    otherwise
                        error('wrong folder name');
                end
                obj.grid=[0,Lx,nx/2];
%                 obj.grid=[-1,1,nx/4;
%                         1,3,nx/4];
            %initialize the continuation
                obj.parameters=[kx ky Ra_T Ra_S2T 1/tau Pr ...
                    dy_T_mean dy_S_mean C c_phase_kx ...
                    Ta_sqrt_z Ta_sqrt_y A_U_bg ...
                    U_r T_r S_r c_phase_z];

        end
        
        function obj=salt_finger_init_primary(obj,folder_name)
            
            %primary bifurcation
            switch folder_name
                case {'salt_finger_kx_low_Ra_S2T_2D_no_slip',...
                        'salt_finger_kx_low_Ra_S2T_3D_no_slip',...
                        'salt_finger_kx_low_Ra_S2T_2D_stress_free',...
                        'salt_finger_kx_low_Ra_S2T_3D_stress_free',...
                        'salt_finger_kx_low_Ra_S2T_high_Pr_2D_no_slip',...
                       'salt_finger_kx_low_Ra_S2T_high_Pr_3D_no_slip'
                        }
                    obj.bpt_list={'bpt1','bpt2','bpt3'};
                    obj.bpt_root_list={'tr'};
                    obj.bifcheck_list=[2,2,2];
                    obj.spcalc_list=[1,1,1];
                    obj.tol_list=[1e-8,1e-8,1e-8];
                    obj.bif_type='steady';
                case {'salt_finger_kx_mid_Ra_S2T_2D_no_slip',...
                        'salt_finger_kx_mid_Ra_S2T_3D_no_slip',...
                        'salt_finger_kx_mid_Ra_S2T_2D_stress_free',...
                        'salt_finger_kx_mid_Ra_S2T_3D_stress_free'}
                    obj.bpt_list={'bpt1','bpt2','bpt3'};
                    obj.bpt_root_list={'tr'};                    
                    obj.bifcheck_list=[2,0,0];
                    obj.spcalc_list=[1,0,0];
                    obj.tol_list=[1e-8,1e-3,1e-3];
                    obj.bif_type='steady';
                    %p.my.bifcheck=2;
                case {'salt_finger_kx_low_Ra_S2T_low_Pr_2D_no_slip',...
                        'salt_finger_kx_low_Ra_S2T_low_Pr_3D_no_slip',...
                        'salt_finger_kx_low_Ra_S2T_low_Pr_2D_stress_free',...
                        'salt_finger_kx_low_Ra_S2T_low_Pr_3D_stress_free'}
                    obj.bpt_list={'bpt1','bpt2','bpt3'};
                    obj.bpt_root_list={'tr'};                    
                    obj.bifcheck_list=[2,2,2];
                    obj.spcalc_list=[1,1,1];
                    obj.tol_list=[1e-8,1e-8,1e-8];
                    obj.bif_type='steady';
                case {'salt_finger_Ra_S2T_3D_no_slip','salt_finger_Ra_S2T_3D_stress_free'}
                    %obj.bpt_list={'bpt1','bpt2','bpt3'};
                    obj.bpt_list={'bpt1'};
                    obj.bpt_root_list={'tr'};
                    obj.bifcheck_list=[2,0,0];
                    obj.spcalc_list=[1,0,0];
                    obj.tol_list=[1e-8,1e-3,1e-3];
                    obj.bif_type='steady';
                case {'salt_finger_kx_low_Ra_S2T_2D_periodic',...
                        'salt_finger_kx_low_Ra_S2T_3D_periodic',...
                        'salt_finger_kx_low_Ra_S2T_low_Pr_2D_periodic',...
                        'salt_finger_kx_low_Ra_S2T_low_Pr_3D_periodic'}
                    %Add the case for the periodic boundary conditions. not
                    %merge to other one just for debug..
                    obj.bpt_list={'bpt1','bpt2','bpt3'};
                    obj.bpt_root_list={'tr'};
                    obj.bifcheck_list=[2,2,2];
                    obj.spcalc_list=[1,1,1];
                    obj.tol_list=[1e-8,1e-8,1e-8];
                    obj.no_phase_step=10;
                    obj.bif_type='steady';
                otherwise
                    error('Wrong folder name');
            end
        end
        
        function obj=salt_finger_init_secondary(obj,folder_name)

            switch folder_name
                case {'salt_finger_kx_low_Ra_S2T_2D_no_slip',...
                        'salt_finger_kx_low_Ra_S2T_3D_no_slip',...
                        'salt_finger_kx_low_Ra_S2T_3D_rotation_no_slip',...
                        'salt_finger_kx_low_Ra_S2T_2D_stress_free',...
                        'salt_finger_kx_low_Ra_S2T_3D_stress_free',...
                        'salt_finger_kx_low_Ra_S2T_high_Pr_2D_no_slip',...
                        'salt_finger_kx_low_Ra_S2T_high_Pr_3D_no_slip'
                        }
                    %secondary bifurcation
                    %for each bpt1, bpt2, bpt3 branch, perform a secondary
                    %bifurcation
                    obj.bpt_root_list={'tr/bpt1','tr/bpt2','tr/bpt3'};
                    obj.bpt_list={'bpt1';
                               'bpt1';
                               'bpt1'};
                    obj.bifcheck_list=[0;0;0];
                    obj.spcalc_list=[1;1;1];
                    obj.bif_type='steady';
                    obj.tol_list=[1e-8; 1e-8; 1e-8];

                case {'salt_finger_kx_low_Ra_S2T_low_Pr_2D_no_slip',...
                        'salt_finger_kx_low_Ra_S2T_low_Pr_3D_no_slip',...
                        'salt_finger_kx_low_Ra_S2T_zero_Pr_3D_no_slip',...
                        'salt_finger_kx_low_Ra_S2T_3D_vorticity_no_slip',...
                        'salt_finger_kx_low_Ra_S2T_zero_Pr_3D_stress_free',...
                        'salt_finger_kx_low_Ra_S2T_3D_vorticity_stress_free'}
                    %secondary bifurcation, just do the bpt1, and first and
                    %second one...
                    obj.bpt_root_list={'tr/bpt1'};
                    obj.bpt_list={'bpt1','bpt2'};
                    obj.bifcheck_list=[2,2];
                    obj.spcalc_list=[1,1];
                    obj.tol_list=[1e-8,1e-8];
                    obj.bif_type='steady';

    %             case {'salt_finger_kx_mid_Ra_S2T_2D'}
    %                 bpt_root_list={'tr/bpt1'};
    %                 bpt_list={'bpt1','bpt2','bpt3'};
                
                case {'salt_finger_kx_low_Ra_S2T_low_Pr_2D_stress_free',...
                        'salt_finger_kx_low_Ra_S2T_low_Pr_3D_stress_free'}
                    obj.bpt_root_list={'tr/bpt1'};
                    obj.bpt_list={'bpt1','bpt2'};
                    obj.bifcheck_list=[2,0];
                    obj.spcalc_list=[1,1];
                    obj.tol_list=[1e-8,1e-8];
                    obj.bif_type='steady';
                case {'salt_finger_kx_mid_Ra_S2T_2D_no_slip',...
                        'salt_finger_kx_mid_Ra_S2T_3D_no_slip',...
                        'salt_finger_kx_mid_Ra_S2T_2D_stress_free',...
                        'salt_finger_kx_mid_Ra_S2T_3D_stress_free'
                        }
                    obj.bpt_root_list={'tr/bpt1'};
                    obj.bpt_list={'bpt1','bpt2'};
                    obj.bifcheck_list=[0,0];
                    obj.spcalc_list=[1,1];
                    obj.tol_list=[1e-6,1e-6];
                    obj.bif_type='steady';

                case {'salt_finger_Ra_S2T_3D_no_slip',...
                        'salt_finger_Ra_S2T_3D_stress_free'}
                    obj.bpt_root_list={'tr/bpt1'};
                    obj.bpt_list={'bpt1'};
                    obj.bifcheck_list=[0];
                    obj.spcalc_list=[1];
                    obj.tol_list=1e-6;
                    obj.bif_type='steady';
                case {'salt_finger_kx_low_Ra_S2T_2D_periodic',...
                        'salt_finger_kx_low_Ra_S2T_3D_periodic',...
                        'salt_finger_kx_low_Ra_S2T_low_Pr_2D_periodic',...
                        'salt_finger_kx_low_Ra_S2T_low_Pr_3D_periodic'}
                    obj.bpt_root_list={'tr/bpt1','tr/bpt2','tr/bpt3'};
                    obj.bpt_list={'bpt1';
                               'bpt1';
                               'bpt1'};
                    obj.bifcheck_list=[0;0;0];
                    obj.spcalc_list=[1;1;1];
                    obj.bif_type='steady';
                    obj.tol_list=[1e-8; 1e-8; 1e-8];
                case 'salt_finger_kx_mid_Ra_S2T_3D_rotation'
                    %I can stop here... 
                    error('1')
                otherwise
                    error('Wrong folder name');
            end 
        end
        
        function obj=salt_finger_init_third(obj,folder_name)
            switch folder_name
                case 'salt_finger_kx_low_Ra_S2T_low_Pr_2D_no_slip'
                    obj.bpt_root_list={'tr/bpt1/bpt1'};
                    obj.bpt_list={'hpt1','hpt2'};
                    obj.bifcheck_list=[0,0];
                    obj.spcalc_list=[0,0];
                    obj.tol_list=[1e-6,1e-6];
                    obj.bif_type='hopf';
                case 'salt_finger_kx_low_Ra_S2T_low_Pr_2D_stress_free'
                    obj.bpt_root_list={'tr/bpt1/bpt1'};
                    obj.bpt_list={'hpt1'};
                    obj.bifcheck_list=[0];
                    obj.spcalc_list=[0];
                    obj.tol_list=[1e-6];
                    obj.bif_type='hopf';
                case 'salt_finger_kx_low_Ra_S2T_2D_no_slip'
                    obj.bpt_root_list={'tr/bpt2'};
                    obj.bpt_list={'hpt1'};
                    obj.bifcheck_list=[0];
                    obj.spcalc_list=[0];
                    obj.tol_list=1e-6;
                    obj.bif_type='hopf';
                otherwise 
                    error('Wrong folder name');
            end
        end
        
        
        function p=swibra_cont(obj,p)
            for bpt_root_ind=1:length(obj.bpt_root_list)
               for bpt_ind=1:length(obj.bpt_list(bpt_root_ind,:))
                   p.my.root_folder_name=[obj.folder_name,obj.bpt_root_list{bpt_root_ind}];
                   p.my.bpt_name=obj.bpt_list{bpt_root_ind,bpt_ind};
                   p.my.store_folder_name_suffix='';
                   p.my.bif_type=obj.bif_type;
                   %p.my.root_foblder_name=[p.my.folder_name,'tr/bpt',num2str(bpt_root_ind)];
                   %p.my.bpt_name=['bpt',num2str(bpt_ind)];
                   if strcmp(obj.bpt_root_list{1},'tr')
                       %from trivial branch, the first step need to enforce
                       % have no phase condition
                       p.my.no_phase_step=obj.no_phase_step;
                       p.my.tol_first=obj.tol_first;
                       p.my.sec_bif_stop=0;
                   else
                       if p.nc.ilam(1)==1 %secondary bifurcation, also only for wavenumber continuation... 
                            p.my.sec_bif_stop=1;
                       end
                   end
                   p.my.bifcheck=obj.bifcheck_list(bpt_root_ind,bpt_ind); 
                   p.my.spcalc=obj.spcalc_list(bpt_root_ind,bpt_ind);
                   p.my.point_list=obj.point_list;
    %                p.my.dsmin=1e-10;
    %                p.my.dsmax=1;
    %                p.my.dsmin=1e-10;
                   p.my.foldcheck=0;
                   p.my.ds=obj.ds;
                   p.my.dsmin=obj.dsmin;
                   p.my.tol=obj.tol_list(bpt_root_ind,bpt_ind);
%                    p.my.tol_first=obj.tol_first;
    %                p.my.par_num=4;
    %                p.my.jac=0;
                   p=my_swibra_cont(p);
               end
           end
        end
        
        function p=tr_cont(obj,p,folder_name)
            %map typical parameter setting from my to p.
            p.pm.mst=p.my.par_num;%number of parallel predictor
            p.pm.resfac=p.my.resfac;
            p.nc.ilam=p.my.ilam; %continuation parameter, the fourth one corresponding to Ra_S2T
            p.nc.lammax=p.my.lammax;
            p.nc.ds=p.my.ds;
            p.nc.neig=p.my.neig;
            p.nc.dsmax=p.my.dsmax;
            p.nc.eigref=p.my.eigref;
            p.nc.intol=p.my.intol;
            %continuation of trivial branch         %starting from the trivial branch.
            switch folder_name
                case {'salt_finger_kx_low_Ra_S2T_2D_no_slip',...
                        'salt_finger_kx_low_Ra_S2T_3D_no_slip',...
                        'salt_finger_kx_low_Ra_S2T_3D_rotation_no_slip',...
                        'salt_finger_kx_low_Ra_S2T_2D_stress_free',...
                        'salt_finger_kx_low_Ra_S2T_3D_stress_free',...
                        'salt_finger_kx_low_Ra_S2T_zero_Pr_3D_stress_free',...
                        'salt_finger_kx_low_Ra_S2T_3D_vorticity_stress_free',...
                        'salt_finger_kx_low_Ra_S2T_low_Pr_2D_stress_free',...
                        'salt_finger_kx_low_Ra_S2T_low_Pr_3D_stress_free',...
                        'salt_finger_kx_low_Ra_S2T_low_Pr_2D_no_slip',...
                        'salt_finger_kx_low_Ra_S2T_low_Pr_3D_no_slip',...
                        'salt_finger_kx_low_Ra_S2T_high_Pr_2D_no_slip',...
                        'salt_finger_kx_low_Ra_S2T_high_Pr_3D_no_slip',...
                        'salt_finger_kx_low_Ra_S2T_zero_Pr_3D_no_slip',...
                        'salt_finger_kx_low_Ra_S2T_3D_vorticity_no_slip',...
                        }
                    p.nc.lammax=-1;
                    p=pmcont(p);
                    p.sol.ds=0.02;
                    p.nc.dsmax=0.3;
                    p.nc.lammax=-0.01;
                    p=pmcont(p);
                case {'salt_finger_kx_low_Ra_S2T_2D_periodic',...
                        'salt_finger_kx_low_Ra_S2T_3D_periodic',...
                        'salt_finger_kx_low_Ra_S2T_low_Pr_2D_periodic',...
                        'salt_finger_kx_low_Ra_S2T_low_Pr_3D_periodic'}
%                     p.nc.ilam=[p.my.ilam,...
%                         p.my.ilam_phase_U_r,...
%                         p.my.ilam_phase_T_r,...
%                         p.my.ilam_phase_S_r];
%                     p.nc.nq=3;
%                     p.my.phase_tr=1;
                    %p.fuha.qf=@qf; p.fuha.qfder=@qfder;
                    
                    p.nc.lammax=-1;
                    p=pmcont(p);
                    p.sol.ds=0.02;
                    p.nc.dsmax=0.3;
                    p.nc.lammax=-0.01;
                    p=pmcont(p);

                case {'salt_finger_kx_mid_Ra_S2T_2D_no_slip',...
                      'salt_finger_kx_mid_Ra_S2T_2D_stress_free'}
                    p.nc.dsmax=0.3;
                    p.nc.lammax=-40;
                    p=pmcont(p);
                case {'salt_finger_kx_mid_Ra_S2T_3D_no_slip',...
                        'salt_finger_kx_mid_Ra_S2T_3D_stress_free'}
                    p.nc.dsmax=0.3;
                    %p.nc.ntot=60;
                    p.nc.lammax=-30;
                    p=pmcont(p);
                case {'salt_finger_Ra_S2T_3D_no_slip',...
                        'salt_finger_Ra_S2T_3D_stress_free'}
                    p.nc.lammax=1100;
                    p=pmcont(p);
                    
                otherwise
                    error('wrong folder name');
            end
            
%                 if p.my.ilam==1 && my.parameters(4)< 10^4
%                     p.nc.lammax=-1;
%                     p=pmcont(p);
%                     p.sol.ds=0.02;
%                     p.nc.dsmax=0.3;
%                     p.nc.lammax=-0.01;
%                     p=pmcont(p);
%                 else
%                     if p.my.ilam==1
%                         p.nc.dsmax=0.3;
%                         %p.nc.ntot=60;
%                         p.nc.lammax=-30;
%                         p=pmcont(p);
%                     elseif p.my.ilam==4
%                         %p.nc.ntot=300;
%                         p.nc.lammax=1100;
%                         p=pmcont(p);
%                     end
%                 end
        end
        
    end
end

