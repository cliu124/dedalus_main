classdef DDC_LST
    %LST_DOUBLE_DIFFUSIVE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        %%parameter for primitive
        Pr=1;%Prandtl number
        tau=1;%diffusivity ratio
        R_rho_T2S=1;%density ratio, temperature over salinity
        
        %%parameter for IFSC
        Ra_ratio=1;%Rayleigh ratio
        
        %%additional parameter for MRBC
        Sc=1;%Schmidt number
        
        %parameter for shear as Radko (2016)
        Pe=1;%real Peclet number
        Ri=1;%Richardson number
        
        %Parameter for unified shear flow parameter
        Re=1;%Reynolds number
        Pe_T=1;%Peclet number appearing in temperature equation
        Pe_S=1;%Peclet number appearing in salinity equation
        Ra_T=1;%Rayleigh nubmer defined as Ra_T=g \alpha dTdz L^4/(\nu\kappa_T) 
        Ra_S2T=1;%Fixed Rayleigh number based on salinity over temperature as Ra_{S,T}=g \beta dSdz L^4/(\nu\kappa_T) 
        
        %%This in default is equal to one... but sometimes will be removed
        %%to study the mechanism...This coefficient can be only 1 or
        %%zero... This is just used to compute results without thermal
        %%diffusivity...
        kappa_T_elevator=1; %%coefficient for the thermal diffusivity..
        kappa_S_elevator=1;
        
        tau_scaling=1;%scaling parameter for diffusivity ratio
    
        %%Note that although here we are doing the double-diffusive
        %%convection. the vertical direction is remain as the y instead of
        %%z to be consistent with shear flow code...
        kx_list=1; %streamwise wavenumber list.
        kz_list=0; %%spanwise wavenumber. This is in default setting up as zero....
    
        kx=1;%streamwise wavenumber
        kz=0;%%THis is the spanwise direction...
        
        eig_val_max;%the eigenvalue associated with the largest real part
        eig_vec_max;%the eigenvector associated with the largest real part
        
        eig_val_max_A_AT_time_dependent;%store eigenvalue of A+A* only involving time-dependent part... 
        eig_vec_max_A_AT_time_dependent;
        
        eig_val_max_A_AT_all;%store the eigenvalue of A+A* for full matrix...
        eig_vec_max_A_AT_all;%eigenvector 
        
        eig_val_max_list={};%putting different kx kz results of eig_val_max into a list
        eig_vec_max_list={};%putting different kx kz results of eig_vec_max into a cell.
        
        eig_val_max_A_AT_time_dependent_list={};%cell storing results for different kx and kz pair
        eig_vec_max_A_AT_time_dependent_list={};%cell storing results for different kx and kz pair
        
        eig_val_max_A_AT_all_list={};%same as agove
        eig_vec_max_A_AT_all_list={};
        
        eig_vec_kx_kz=struct;%eigenvector associated with the largest growth rate among different kx kz..
        
        %dy_T_mean=dy_S_mean=1 for salt finger, and dy_T_mean=dy_S_mean=-1
        %for diffusive regime
        dy_T_mean=1; %background temperature gradient
        dy_S_mean=1; %background salinty gradient
        
        
        mean='no';% background shear flow
        mean_kolmogorov=zeros(4,3);%%parameter to set up the kolmogorov flow...
        Ny_full=62;%the number of grid points in the shearwise direction, denoted as y here
        
        %%Grid points and the differential matrix...
        y_list_full;%grid points in y
        D1_full;%first to 4th order differential matrix
        D2_full;
        D3_full;
        D4_full;
        
        U_bar_full;%background mean velocity on the y_list_full
        d_U_bar_full;%dU/dy on the y_list_full
        dd_U_bar_full;%ddU/dyy on the y_list_full
            
        %The background gradient of the temperature and salinity...
        %computed from elevator mode...
        d_T_bar_full;
        d_S_bar_full;
        
        grid_diff={'fourier',[0,2*pi]}; %%This is in default using the fourier collocation method with domain [0,2*pi].
%         mean_elevator_W_list Note that the domain in Radko (2016) is up
%         to 1...
        
        %The sub flow within the unified formulation of double diffusive
        flow_sub_double_diffusive_shear_2D='primitive_Radko2013'
        
        %the reduced equations within the shear flow formulation of Radko
        %(2016)
        shear_Radko2016_reduced='primitive';
        
        %the flag for the post processing of eigenvalue comparison with
        %Radko (2016)
        post_eig_Radko2016=1; %%flag for post-processing. whether compare with their results
        post_eigenvector_Radko2016_unit=1;%convert the unit of eigenvector to be consistent with Radko 2016. This is in general useful that move density ratio to salinty equation
        post_eig_Holyer1984=1; %compare with the secondary instability analysis of Holyer 1984
        
        %path to store the figure and data.
        path_fig='C:/Figure/DDC_LST/';
        path_data='C:/Data/DDC_LST/';
        
        mean_elevator_amp_list={'T',0};%The amplitude of the elevator mode. The first symbol select the variable that you want to specify and the scond is the amplitude value        
        mean_elevator_kx=1;%This is the wavenumber index provided by the user
        mean_elevator_kx_max=-1;%%default value meaning they have not been computed
        mean_elevator_kx_steady=-1;%default value meaning they have not been computed
        mean_elevator_kx_local=0; %%This is the local wavenumber for elevator mode that will be modifed in the code for computing...
        
        elevator_lambda_max=0;%The largest growth rate of the elevator mode, maximized over kx.
        elevator_eig_vec_max=0;%The eigenvector associated with the largest growth rate of the elevator mode
        
        %These parameter elevator_lambda_balance are used to do the balance
        %between primary and secondary instability of elevator mode using
        %frozen assumption...
        elevator_lambda_balance_C=4;%%The C parameter used to do the growth rate balance as Radko & Smith (2012)
        elevator_lambda_balance_bisection=[0,30,0.001]; %%The minimal bound for bisection, the maximal bound for bisetion, and the residue error
        elevator_lambda_balance_error=1; %The error of the 
        
        %results name
        result_name='test';
        debug='';
        %the solver I would like to call
        solve='finished';
        
        %Flux of elevator mode, based on equation (3.12) of Radko &
        %Smith(2012)
        F_T=0;
        F_S=0;
        
        A;%%operator for the linear stability analysis
        A_time_dependent;%%operator that introduced by elevator mode only..
        
        M;
        
        operator='v_omega_y'; %%set the flag for operator..

        darcy=0; %%if 1, then we need to change the viscosity term based on the darcy law... Note that this is only for momentum equation
    
        flux_T=0;%whether this is for fixed flux results... 
    end
    
    methods
        %%------------
        %%First part of the methods... core solver... 
        function obj = DDC_LST(flow_sub_double_diffusive_shear_2D)
            %LST_DOUBLE_DIFFUSIVE Construct an instance of this class
            %   Detailed explanation goes here
            %only need to give the name of which sub flow within the shear
            %flow unit..
            obj.flow_sub_double_diffusive_shear_2D=flow_sub_double_diffusive_shear_2D;
        end
        
        function obj = convert_shear(obj)
            %%convert the parameter into the shear flow unit...
            
            switch obj.flow_sub_double_diffusive_shear_2D
                case 'primitive_Radko2013'
                    %This is convert the formulation into the double
                    %diffusive convection of Radko (2013). Note the density
                    %ratio is on the buoyancy term in the momentum. 
                    obj.Re=1/obj.Pr;
                    obj.Pe_T=1;
                    obj.Pe_S=1;
                    obj.tau=obj.tau; %#Set this as zero if remove salinity diffusivity
                    obj.Ra_T=1;
                    obj.Ra_S2T=1/obj.R_rho_T2S;
                    obj.tau_scaling=obj.tau;

                case 'primitive_IFSC_unit_tuS'
%                 ##parameter for primitive equations
              
              %  #map to the extended parameter in primitive_IFSC_unit_tus
                    obj.Re=obj.tau/obj.Pr;
                    obj.Ra_S2T=1/obj.R_rho_T2S/obj.tau;
                    obj.Pe_T=obj.tau;
                    obj.Pe_S=1;
                    obj.tau=1; %#Set this as zero if remove salinity diffusivity
                    obj.Ra_T=1;
                    obj.tau_scaling=1;
                case 'IFSC'
                    %#map to the extended parameter in double_diffusive_shear_2D
                    obj.Re=0;%no inertial term in momentum
                    obj.Pe_T=0;%no inertial term in the temperature
                    obj.Pe_S=1;
                    obj.tau=1; %#Set this as zero if remove salinity diffusivity
                    obj.Ra_T=1;
                    obj.Ra_S2T=obj.Ra_ratio;
                    obj.tau_scaling=1;
                case 'MRBC'
                
                    %#map to the extended parameter in double_diffusive_shear_2D
                    obj.Re=1/obj.Sc;
                    obj.Pe_T=0;%no inertial term in the temperature
                    obj.Pe_S=1;
                    obj.tau=1; %#Set this as zero if remove salinity diffusivity
                    obj.Ra_T=1;
                    obj.Ra_S2T=obj.Ra_ratio;
                    obj.tau_scaling=1;
                case 'Stokes'
                    %#map to the extended parameter in double_diffusive_shear_2D
                    obj.Re=0;%no inertial term in the momentum
                    obj.Pe_T=obj.tau;
                    obj.Pe_S=1;
                    obj.tau=1; %#Set this as zero if remove salinity diffusivity
                    obj.Ra_T=1;
                    obj.Ra_S2T=obj.Ra_ratio;
                    obj.tau_scaling=1;
                case 'unified'
                    %%need to do nothing...
                    
                case 'shear_Radko2016'
               
                    %#map to the extended parameter in double_diffusive_shear_2D
                    switch obj.shear_Radko2016_reduced
                        %set different term in the unified flow
                        %configuration as zero to obtain the reduced
                        %formulation
                        case 'primitive'%full primitive
                            obj.Re=obj.Pe/obj.Pr;
                            obj.Pe_T=obj.Pe;
                        case 'IFSC'%inertial free salt-finger, no inertial term in momentum and temperature
                            obj.Re=0;
                            obj.Pe_T=0;
                        case 'MRBC'%modified Rayleigh Benard convection, no inertial term in temperature
                            obj.Re=obj.Pe/obj.Pr;
                            obj.Pe_T=0;
                        case 'Stokes'%s
                            obj.Re=0;
                            obj.Pe_T=obj.Pe;
                        otherwise
                            error('Wrong obj.shear_Radko2016_reduced');
                    end
                    obj.Pe_S=obj.Pe;
                    obj.tau=obj.tau;
                    obj.Ra_T=4*pi^2*obj.Ri/(1/obj.R_rho_T2S-1)*obj.Pe^2/obj.Pr;
                    obj.Ra_S2T=obj.Ra_T/obj.R_rho_T2S;
                    
                    obj.tau_scaling=1;
                otherwise 
                    error('Wrong obj.flow_sub_double_diffusive_shear_2D');
            end
        end
        
        function obj=operator_A_M(obj)
            %perform the linear stability analysis for one kx and kz...
            %initilization.. remove these A and M field.
%             obj=rmfield(obj,'A');
%             obj=rmfield(obj,'M');
%             obj.A=[];
%             obj.M=[];
%             obj.A_time_dependent=[];
%             
            zi=sqrt(-1); %imaginary unit
            
            U_bar=obj.U_bar_full; 
            d_U_bar=obj.d_U_bar_full;
            dd_U_bar=obj.dd_U_bar_full;
            
            I_bc=eye(obj.Ny_full,obj.Ny_full);
            zero_bc=zeros(obj.Ny_full,obj.Ny_full);
            D4_bc=obj.D4_full;
            D2_bc=obj.D2_full;
            D1_bc=obj.D1_full;
            K2= obj.kx^2+obj.kz^2; % Total wave number, in convienent for calculation
            
            switch obj.operator
                case 'v_omega_y' %based on Orr-Sommerfeld and Squire operator...
                    for mean_elevator_W_ind=1:length(obj.mean_elevator_amp_list{2})
                        if K2>0.000001 
                            switch obj.mean
                                case {'kolmogorov','no'}%primary instability
                                    %construct the A matrix for standard shear flow. Note that
                                    %Re is associated with the shear term
                                    A11=(D4_bc-2*K2*D2_bc+K2^2*I_bc)+obj.Re*(diag(dd_U_bar)*zi*obj.kx*I_bc-zi*obj.kx*diag(U_bar)*(D2_bc-K2*I_bc)); %%Orr-Sommerfeld operator
                                    A21= -zi*obj.kz*diag(d_U_bar)*I_bc*obj.Re; %Coulping operator
                                    A22= -zi*obj.kx*diag(U_bar)*I_bc*obj.Re+(D2_bc-K2*I_bc); %Squire operator  
                                    %inv_lap=inv([D2_bc-K2*I_bc, zero_bc; zero_bc, I_bc]);
                                    inv_lap=inv(D2_bc-K2*I_bc);
                                    A_shear= [inv_lap*A11, zero_bc; A21, A22];

                                    %%standard shear flow model couple with T and S
                                    obj.A(:,:,mean_elevator_W_ind)=[A_shear, obj.Ra_T*[inv_lap*(-K2*I_bc); zero_bc], -obj.Ra_S2T*[inv_lap*(-K2*I_bc); zero_bc];
                                      -obj.dy_T_mean*I_bc,zero_bc,-zi*obj.kx*diag(U_bar)*I_bc*obj.Pe_T+(D2_bc-K2*I_bc),zero_bc;
                                      -obj.dy_S_mean*I_bc,zero_bc,zero_bc,-zi*obj.kx*diag(U_bar)*I_bc*obj.Pe_S+obj.tau*(D2_bc-K2*I_bc)];
                                    
                                case 'elevator'%secondary instability of elevator mode

                                    A11=(D4_bc-2*K2*D2_bc+K2^2*I_bc)+obj.Re*(diag(dd_U_bar(:,mean_elevator_W_ind))*zi*obj.kx*I_bc-zi*obj.kx*diag(U_bar(:,mean_elevator_W_ind))*(D2_bc-K2*I_bc)); %%Orr-Sommerfeld operator
                                    A21= -zi*obj.kz*diag(d_U_bar(:,mean_elevator_W_ind))*I_bc*obj.Re; %Coulping operator
                                    A22= -zi*obj.kx*diag(U_bar(:,mean_elevator_W_ind))*I_bc*obj.Re+(D2_bc-K2*I_bc); %Squire operator  
                                    inv_lap=inv(D2_bc-K2*I_bc);%, zero_bc; zero_bc, I_bc]);
                                    A_shear= [inv_lap*A11, zero_bc; A21, A22];

                                    Cu=[zi*obj.kx*D1_bc, -zi*obj.kz*I_bc]/K2;
                                    Cv=[I_bc,zero_bc];
                                    Bu=[inv_lap*(-zi*obj.kx*D1_bc); zi*obj.kz*I_bc];
                                    obj.A(:,:,mean_elevator_W_ind)=[A_shear, obj.Ra_T*Bu, -obj.Ra_S2T*Bu;
                                                    -diag(obj.d_T_bar_full(:,mean_elevator_W_ind))*Cv*obj.Pe_T-diag(obj.dy_T_mean)*Cu,-zi*obj.kx*diag(U_bar(:,mean_elevator_W_ind))*I_bc*obj.Pe_T+obj.kappa_T_elevator*(D2_bc-K2*I_bc),zero_bc;
                                                    -diag(obj.d_S_bar_full(:,mean_elevator_W_ind))*Cv*obj.Pe_S-diag(obj.dy_S_mean)*Cu,zero_bc,-zi*obj.kx*diag(U_bar(:,mean_elevator_W_ind))*I_bc*obj.Pe_S+obj.tau*obj.kappa_S_elevator*(D2_bc-K2*I_bc)];
                                    
                                              
                                                
                                otherwise 
                                    error('obj.mean not supported in this case');
                            end

                        elseif K2<=0.000001 %For small kx and kz, avoide singularity... use the formulation directly set up kx=kz=0.
                            %%Also, we need to avoide the inverse of Laplacian as there
                            %%is a zero eigenvalue of D2_bc, and it will lead to
                            %%singularity when computing the inverse.
                            switch obj.mean
                                case {'kolmogorov','no'}
                                    A_shear= [D2_bc, zero_bc; zero_bc, D2_bc]; 
                                    obj.A(:,:,mean_elevator_W_ind)=[A_shear, obj.Ra_T*[zero_bc; zero_bc], -obj.Ra_S2T*[zero_bc; zero_bc];
                                    -obj.dy_T_mean*I_bc,zero_bc,(D2_bc-K2*I_bc),zero_bc;
                                    -obj.dy_S_mean*I_bc,zero_bc,zero_bc,obj.tau*(D2_bc-K2*I_bc)];
                                case 'elevator'

                                    %%Note that this is different...
                                    %%for the elevator mode,,, the kx is actually
                                    %%kz (vertical wavenumber).
                                    %%Thus, kx=kz=0 will get rid of u velocity and
                                    %%it still remain the buoyancy force here.
                                    A_shear= [D2_bc, zero_bc; zero_bc, D2_bc]; 
                                    obj.A(:,:,mean_elevator_W_ind)=[A_shear, obj.Ra_T*[I_bc; zero_bc], -obj.Ra_S2T*[I_bc; zero_bc];
                                    -obj.dy_T_mean*I_bc,zero_bc,obj.kappa_T_elevator*D2_bc,zero_bc;
                                    -obj.dy_S_mean*I_bc,zero_bc,zero_bc,obj.tau*obj.kappa_S_elevator*D2_bc];
                                    %A11=(D4_bc-2*K2*D2_bc+K2^2*I_bc)+obj.Re*(diag(dd_U_bar(:,mean_elevator_W_ind))*zi*obj.kx*I_bc-zi*obj.kx*diag(U_bar(:,mean_elevator_W_ind))*(D2_bc-K2*I_bc)); %%Orr-Sommerfeld operator
                                    %A21= -zi*obj.kz*diag(d_U_bar(:,mean_elevator_W_ind))*I_bc*obj.Re; %Coulping operator
                                    %A22= -zi*obj.kx*diag(U_bar(:,mean_elevator_W_ind))*I_bc*obj.Re+(D2_bc-K2*I_bc); %Squire operator  
                                    %inv_lap=inv(D2_bc-K2*I_bc);%, zero_bc; zero_bc, I_bc]);
                                    %A_shear= [inv_lap*A11, zero_bc; A21, A22];
        %                             A_shear= [D2_bc, zero_bc; zero_bc, D2_bc]; 
        % 
        %                             Cu=[zi*obj.kx*D1_bc, -zi*obj.kz*I_bc]/K2;
        %                             Cv=[I_bc,zero_bc];
        %                             Bu=[inv_lap*(-zi*obj.kx*D1_bc); zi*obj.kz*I_bc];
        %                             A(:,:,mean_elevator_W_ind)=[A_shear, obj.Ra_T*Bu, -obj.Ra_S2T*Bu;
        %                                             -diag(obj.d_T_bar_full(:,mean_elevator_W_ind))*Cv*obj.Pe_T-diag(obj.dy_T_mean)*Cu,-zi*obj.kx*diag(U_bar(:,mean_elevator_W_ind))*I_bc*obj.Pe_T+(D2_bc-K2*I_bc),zero_bc;
        %                                             -diag(obj.d_S_bar_full(:,mean_elevator_W_ind))*Cv*obj.Pe_S-diag(obj.dy_S_mean)*Cu,zero_bc,-zi*obj.kx*diag(U_bar(:,mean_elevator_W_ind))*I_bc*obj.Pe_S+obj.tau*(D2_bc-K2*I_bc)];
        %                             
                                otherwise
                                    error('obj.mean not supported in this case');
                            end
                        end
                        
                        
                        if strcmp(obj.mean,'elevator')
                            A_shear_elevator=[obj.Re*(diag(dd_U_bar(:,mean_elevator_W_ind))*zi*obj.kx*I_bc-zi*obj.kx*diag(U_bar(:,mean_elevator_W_ind))*(D2_bc-K2*I_bc)), zero_bc;
                                -zi*obj.kz*diag(d_U_bar(:,mean_elevator_W_ind))*I_bc*obj.Re,-zi*obj.kx*diag(U_bar(:,mean_elevator_W_ind))*I_bc*obj.Re];
                            obj.A_time_dependent(:,:,mean_elevator_W_ind)=[A_shear_elevator,[zero_bc;zero_bc], [zero_bc;zero_bc];
                                            -diag(obj.d_T_bar_full(:,mean_elevator_W_ind))*Cv*obj.Pe_T,-zi*obj.kx*diag(U_bar(:,mean_elevator_W_ind))*I_bc*obj.Pe_T,zero_bc;
                                            -diag(obj.d_S_bar_full(:,mean_elevator_W_ind))*Cv*obj.Pe_S,zero_bc,-zi*obj.kx*diag(U_bar(:,mean_elevator_W_ind))*I_bc*obj.Pe_S];
                        end
                        
                    end

                    %%the matrix for the genearalized eigenvalue problem
                    obj.M=[obj.Re*I_bc, zero_bc,zero_bc,zero_bc;
                       zero_bc, obj.Re*I_bc, zero_bc, zero_bc;
                       zero_bc, zero_bc, obj.Pe_T*I_bc, zero_bc;
                       zero_bc, zero_bc, zero_bc, obj.Pe_S*I_bc];
                case 'uvwpTS'%formulation using primitive, with uvw p and T, S
                    
                    for mean_elevator_W_ind=1:length(obj.mean_elevator_amp_list{2})
                        if obj.darcy %%if darcy,, modifiy the advection diffusion operator based on Darcy law...
                            Diagterm=-1i*obj.kx*diag(U_bar(:,mean_elevator_W_ind))*obj.Re-I_bc;
                        else
                            Diagterm=-1i*obj.kx*diag(U_bar(:,mean_elevator_W_ind))*obj.Re+(D2_bc-K2*I_bc);
                        end
                                    
                        if K2>0.000001
                            switch obj.mean
                                case {'kolmogorov','no'}
                                    obj.A(:,:,mean_elevator_W_ind)=[Diagterm, -diag(d_U_bar(:,mean_elevator_W_ind))*obj.Re, zero_bc, -1i*obj.kx*I_bc, zero_bc, zero_bc;...
                                          zero_bc, Diagterm, zero_bc,     -D1_bc, obj.Ra_T*I_bc, -obj.Ra_S2T*I_bc;...
                                          zero_bc, zero_bc, Diagterm,     -1i*obj.kz*I_bc, zero_bc, zero_bc; ...
                                          1i*obj.kx*I_bc,   D1_bc,   1i*obj.kz*I_bc,    zero_bc, zero_bc, zero_bc;
                                          zero_bc, -obj.dy_T_mean*I_bc, zero_bc, zero_bc,-1i*obj.kx*diag(U_bar(:,mean_elevator_W_ind))*obj.Pe_T+(D2_bc-K2*I_bc), zero_bc;
                                          zero_bc, -obj.dy_S_mean*I_bc, zero_bc, zero_bc, zero_bc, -1i*obj.kx*diag(U_bar(:,mean_elevator_W_ind))*obj.Pe_S+obj.tau*(D2_bc-K2*I_bc)]; 
                                case 'elevator'
                                    obj.A(:,:,mean_elevator_W_ind)=[Diagterm, -diag(d_U_bar(:,mean_elevator_W_ind))*obj.Re, zero_bc, -1i*obj.kx*I_bc, obj.Ra_T*I_bc, -obj.Ra_S2T*I_bc;...
                                          zero_bc, Diagterm, zero_bc,     -D1_bc, zero_bc,zero_bc;...
                                          zero_bc, zero_bc, Diagterm,     -1i*obj.kz*I_bc, zero_bc, zero_bc; ...
                                          1i*obj.kx*I_bc,D1_bc,1i*obj.kz*I_bc,zero_bc, zero_bc, zero_bc;
                                          -obj.dy_T_mean*I_bc, -diag(obj.d_T_bar_full(:,mean_elevator_W_ind))*obj.Pe_T, zero_bc, zero_bc,-1i*obj.kx*diag(U_bar(:,mean_elevator_W_ind))*obj.Pe_T+obj.kappa_T_elevator*(D2_bc-K2*I_bc), zero_bc;
                                          -obj.dy_S_ean*I_bc, -diag(obj.d_S_bar_full(:,mean_elevator_W_ind))*obj.Pe_S, zero_bc, zero_bc, zero_bc, -1i*obj.kx*diag(U_bar(:,mean_elevator_W_ind))*obj.Pe_S+obj.tau*obj.kappa_S_elevator*(D2_bc-K2*I_bc)]; 
                                  
                                otherwise
                                    error('obj.mean not supported in this case');
                            end
                            obj.M=blkdiag(obj.Re*I_bc,obj.Re*I_bc,obj.Re*I_bc,zero_bc,obj.Pe_T*I_bc,obj.Pe_S*I_bc);

                        elseif K2<=0.000001
                            error('This branch has not been fully debug.. set a larger kx')
                            switch obj.mean
                                case {'kolmogorov','no'}
                                    obj.A(:,:,mean_elevator_W_ind)=[Diagterm, -diag(d_U_bar(:,mean_elevator_W_ind)), zero_bc, zero_bc, zero_bc, zero_bc;...
                                          zero_bc, Diagterm, zero_bc,     -zero_bc, obj.Ra_T*zero_bc, -obj.Ra_S2T*zero_bc;...
                                          zero_bc, zero_bc, Diagterm,     -zero_bc, zero_bc, zero_bc; ...
                                          zero_bc,zero_bc, zero_bc,zero_bc, zero_bc, zero_bc;
                                          zero_bc, -obj.dy_T_mean*I_bc, zero_bc, zero_bc,(D2_bc-K2*I_bc), zero_bc;
                                          zero_bc, -obj.dy_S_mean*I_bc, zero_bc, zero_bc, zero_bc,obj.tau*(D2_bc-K2*I_bc)]; 
                                case 'elevator'
                                     obj.A(:,:,mean_elevator_W_ind)=[Diagterm, -diag(d_U_bar(:,mean_elevator_W_ind)), zero_bc, zero_bc, zero_bc, zero_bc;...
                                          zero_bc, Diagterm, zero_bc,     -zero_bc, obj.Ra_T*I_bc, -obj.Ra_S2T*I_bc;...
                                          zero_bc, zero_bc, Diagterm,     -zero_bc, zero_bc, zero_bc; ...
                                          zero_bc, zero_bc, zero_bc,I_bc, zero_bc, zero_bc;
                                          zero_bc, -obj.dy_T_mean*I_bc, zero_bc, zero_bc,-1i*obj.kx*diag(U_bar(:,mean_elevator_W_ind))*obj.Pe_T+obj.kappa_T_elevator*(D2_bc-K2*I_bc), zero_bc;
                                          zero_bc, -obj.dy_S_mean*I_bc, zero_bc, zero_bc, zero_bc, -1i*obj.kx*diag(U_bar(:,mean_elevator_W_ind))*obj.Pe_S+obj.tau*obj.kappa_S_elevator*(D2_bc-K2*I_bc)]; 
                                otherwise
                                    error('obj.mean not supported in this case');
                            end
                            obj.M=blkdiag(obj.Re*I_bc,obj.Re*I_bc,obj.Re*I_bc,zero_bc,obj.Pe_T*I_bc,obj.Pe_S*I_bc);

                        end
                        
                        if strcmp(obj.mean,'elevator')
                            %%This is the term only induced by the elevator
                            %%mode...Update 2021/10/15
                            Diagterm_elevator=-1i*obj.kx*diag(U_bar(:,mean_elevator_W_ind))*obj.Re;
                            obj.A_time_dependent(:,:,mean_elevator_W_ind)=[Diagterm_elevator, -diag(d_U_bar(:,mean_elevator_W_ind))*obj.Re, zero_bc, zero_bc, zero_bc, zero_bc;...
                                  zero_bc, Diagterm_elevator, zero_bc,     zero_bc, zero_bc,zero_bc;...
                                  zero_bc, zero_bc, Diagterm_elevator,     zero_bc, zero_bc, zero_bc; ...
                                  zero_bc, zero_bc, zero_bc,zero_bc, zero_bc, zero_bc;
                                  zero_bc, -diag(obj.d_T_bar_full(:,mean_elevator_W_ind))*obj.Pe_T, zero_bc, zero_bc,-1i*obj.kx*diag(U_bar(:,mean_elevator_W_ind))*obj.Pe_T, zero_bc;
                                  zero_bc, -diag(obj.d_S_bar_full(:,mean_elevator_W_ind))*obj.Pe_S, zero_bc, zero_bc, zero_bc, -1i*obj.kx*diag(U_bar(:,mean_elevator_W_ind))*obj.Pe_S];
                        end
                        
                    end
                    
                    %six block diagonal matrix to perform the 
                    %[obj.Re*I_bc, zero_bc,zero_bc,zero_bc;
                    %   zero_bc, obj.Re*I_bc, zero_bc, zero_bc;
                    %   zero_bc, zero_bc, obj.Pe_T*I_bc, zero_bc;
                    %   zero_bc, zero_bc, zero_bc, obj.Pe_S*I_bc];
                
                otherwise
                    error('Wrong obj.operator');
            end
        end
        
        function obj=LST(obj)
            if strcmp(obj.mean,'elevator') && obj.grid_diff{2}(2) == 2*pi 
                %%This require the mean flow is the elevator mode and we do
                %%not specify any user defined domain size.. it remain
                %%2pi...
                %%make the grid adaptive to the elevator mode...
                %get 2 times the wavelength of elevator mode
                obj.grid_diff{2}(2)=2*pi/obj.mean_elevator_kx_local;
            end
            obj=grid_diff_fourier(obj); %%setup the fourier grid
            obj=mean_profile(obj); %%set up the mean profile..
            
            obj=obj.operator_A_M();
            
            for mean_elevator_W_ind=1:length(obj.mean_elevator_amp_list{2})

                %%This is a list of three different type of eigenvalue
                %%problem we want to solve....
                if strcmp(obj.mean,'elevator') %strcmp(obj.operator,'v_omega_y') && 
                    eig_name_list={'max','max_A_AT_all'};%,'max_A_AT_time_dependent'
                    %eig_name_list={'max'};
                else
                    eig_name_list={'max'};
                end
                
                %Update 2021/10/21. just go through the eig_name_list to
                %compute and store the eigenvalue and eigenvector... such
                %that I do not need to copy the code...
                for eig_name_ind=1:length(eig_name_list)
                    switch eig_name_list{eig_name_ind}
                        case 'max'
                            A=obj.A(:,:,mean_elevator_W_ind);
                            M=obj.M;
                        case {'max_A_AT_all','max_A_AT_time_dependent'}
                            if strcmp(eig_name_list{eig_name_ind},'max_A_AT_time_dependent')
                                A=obj.A_time_dependent(:,:,mean_elevator_W_ind);
                            elseif strcmp(eig_name_list{eig_name_ind},'max_A_AT_all')
                                A=obj.A(:,:,mean_elevator_W_ind);
                            end
                            zero_ind=find(diag(obj.M)==0);
                            non_zero_ind=find(diag(obj.M)~=0);
                            if isempty(zero_ind) %%The M matrix is not singular
                                A=inv(obj.M)*A;%+obj.A_time_dependent(:,:,mean_elevator_W_ind)';
                            else %%This step will eliminate the algebraic constraints...
                                A=A(non_zero_ind,non_zero_ind)-...
                                    A(non_zero_ind,zero_ind)*inv(A(zero_ind,zero_ind))*A(zero_ind,non_zero_ind);
                                A=inv(obj.M(non_zero_ind,non_zero_ind))*A;
                            end
                            A=A+A';
                            M=eye(size(A));
%                         case 'max_A_AT_all'
%                             A=inv(obj.M)*obj.A(:,:,mean_elevator_W_ind);
%                             A=A+A';%+obj.A(:,:,mean_elevator_W_ind)';
%                             M=eye(size(A));
                        otherwise
                            error('Wrong eig_name_list{eig_name_ind}');
                    end
                    
                    [eig_vec,eig_val_mat]=eig(A,M);
                    %Update 2021/10/13, also remove the eigenvalue that has
                    %growth rate much larger than 10^6...
                    eig_val_mat(isinf(eig_val_mat)|isnan(eig_val_mat)|real(eig_val_mat)>10^6) = -Inf;
                    [~,eig_val_max_ind]=max(real(diag(eig_val_mat))); %#compute the index of the eigenvalue
                    obj.(['eig_val_',eig_name_list{eig_name_ind}])(1,mean_elevator_W_ind)=eig_val_mat(eig_val_max_ind,eig_val_max_ind);
                    obj.(['eig_vec_',eig_name_list{eig_name_ind}])(:,mean_elevator_W_ind)=eig_vec(:,eig_val_max_ind); %#get the corresponding eigen vector

                end
%                 if strcmp(obj.operator,'v_omega_y') && strcmp(obj.mean,'elevator')
% 
% %                     %                     disp('Also compute the eigenvalue of A+A^T');
%                     A=inv(obj.M)*squeeze(obj.A_time_dependent(:,:,mean_elevator_W_ind));
%                     [eig_vec_A_AT_time_dependent,eig_val_mat_A_AT_time_dependent]=...
%                         max(eig(A+A'));
%                     eig_val_mat_A_AT_time_dependent(isinf(eig_val_mat_A_AT_time_dependent)|isnan(eig_val_mat_A_AT_time_dependent)|real(eig_val_mat_A_AT_time_dependent)>10^6) = -Inf;
%                     [~,eig_val_max_ind]=max(real(diag(eig_val_mat_A_AT_time_dependent))); %#compute the index of the eigenvalue
%                     obj.eig_val_max_A_AT_time_dependent(1,mean_elevator_W_ind)=eig_val_mat(eig_val_max_ind,eig_val_max_ind);
%                     obj.eig_val_max_A_AT_time_dependent(:,mean_elevator_W_ind)=eig_vec_A_AT_time_dependent(:,eig_val_max_ind); %#get the corresponding eigen vector
%                 
%                     A_all=inv(obj.M)*squeeze(obj.A(:,:,mean_elevator_W_ind));
%                     [eig_vec_A_AT_all,eig_val_mat_A_AT_all]=...
%                         max(eig(A_all+A_all'));
%                     eig_val_mat_A_AT_all(isinf(eig_val_mat_A_AT_all)|isnan(eig_val_mat_A_AT_all)|real(eig_val_mat_A_AT_all)>10^6) = -Inf;
%                     [~,eig_val_max_ind]=max(real(diag(eig_val_mat_A_AT_all))); %#compute the index of the eigenvalue
%                     obj.eig_val_max_A_AT_all(1,mean_elevator_W_ind)=eig_val_mat(eig_val_max_ind,eig_val_max_ind);
%                     obj.eig_val_max_A_AT_all(:,mean_elevator_W_ind)=eig_vec_A_AT_all(:,eig_val_max_ind); %#get the corresponding eigen vector
%                 end
                
            end
        end
        
%         

%         function obj=LST(obj)
%             
%             %perform the linear stability analysis for one kx and kz...
%             zi=sqrt(-1); %imaginary unit
%             if strcmp(obj.mean,'elevator')
%                 %%make the grid adaptive to the elevator mode...
%                 obj.grid_diff{2}(2)=2*pi/obj.mean_elevator_kx_local;
%             end
%             obj=grid_diff_fourier(obj); %%setup the fourier grid
%             obj=mean_profile(obj); %%set up the mean profile..
%             
%             U_bar=obj.U_bar_full; 
%             d_U_bar=obj.d_U_bar_full;
%             dd_U_bar=obj.dd_U_bar_full;
%             
%             I_bc=eye(obj.Ny_full,obj.Ny_full);
%             zero_bc=zeros(obj.Ny_full,obj.Ny_full);
%             D4_bc=obj.D4_full;
%             D2_bc=obj.D2_full;
%             D1_bc=obj.D1_full;
%             K2= obj.kx^2+obj.kz^2; % Total wave number, in convienent for calculation
%             for mean_elevator_W_ind=1:length(obj.mean_elevator_amp_list{2})
%                 if K2>0.000001 
%                     switch obj.mean
%                         case {'kolmogorov','no'}
%                             %construct the A matrix for standard shear flow. Note that
%                             %Re is associated with the shear term
%                             A11=(D4_bc-2*K2*D2_bc+K2^2*I_bc)+obj.Re*(diag(dd_U_bar)*zi*obj.kx*I_bc-zi*obj.kx*diag(U_bar)*(D2_bc-K2*I_bc)); %%Orr-Sommerfeld operator
%                             A21= -zi*obj.kz*diag(d_U_bar)*I_bc*obj.Re; %Coulping operator
%                             A22= -zi*obj.kx*diag(U_bar)*I_bc*obj.Re+(D2_bc-K2*I_bc); %Squire operator  
%                             %inv_lap=inv([D2_bc-K2*I_bc, zero_bc; zero_bc, I_bc]);
%                             inv_lap=inv(D2_bc-K2*I_bc);
%                             A_shear= [inv_lap*A11, zero_bc; A21, A22];
% 
%                             %%standard shear flow model couple with T and S
%                             A(:,:,mean_elevator_W_ind)=[A_shear, obj.Ra_T*[inv_lap*(-K2*I_bc); zero_bc], -obj.Ra_S2T*[inv_lap*(-K2*I_bc); zero_bc];
%                               -obj.dy_T_mean*I_bc,zero_bc,-zi*obj.kx*diag(U_bar)*I_bc*obj.Pe_T+(D2_bc-K2*I_bc),zero_bc;
%                               -obj.dy_S_mean*I_bc,zero_bc,zero_bc,-zi*obj.kx*diag(U_bar)*I_bc*obj.Pe_S+obj.tau*(D2_bc-K2*I_bc)];
% 
%                         case 'elevator'
% 
%                             A11=(D4_bc-2*K2*D2_bc+K2^2*I_bc)+obj.Re*(diag(dd_U_bar(:,mean_elevator_W_ind))*zi*obj.kx*I_bc-zi*obj.kx*diag(U_bar(:,mean_elevator_W_ind))*(D2_bc-K2*I_bc)); %%Orr-Sommerfeld operator
%                             A21= -zi*obj.kz*diag(d_U_bar(:,mean_elevator_W_ind))*I_bc*obj.Re; %Coulping operator
%                             A22= -zi*obj.kx*diag(U_bar(:,mean_elevator_W_ind))*I_bc*obj.Re+(D2_bc-K2*I_bc); %Squire operator  
%                             inv_lap=inv(D2_bc-K2*I_bc);%, zero_bc; zero_bc, I_bc]);
%                             A_shear= [inv_lap*A11, zero_bc; A21, A22];
% 
%                             Cu=[zi*obj.kx*D1_bc, -zi*obj.kz*I_bc]/K2;
%                             Cv=[I_bc,zero_bc];
%                             Bu=[inv_lap*(-zi*obj.kx*D1_bc); zi*obj.kz*I_bc];
%                             A(:,:,mean_elevator_W_ind)=[A_shear, obj.Ra_T*Bu, -obj.Ra_S2T*Bu;
%                                             -diag(obj.d_T_bar_full(:,mean_elevator_W_ind))*Cv*obj.Pe_T-diag(obj.dy_T_mean)*Cu,-zi*obj.kx*diag(U_bar(:,mean_elevator_W_ind))*I_bc*obj.Pe_T+(D2_bc-K2*I_bc),zero_bc;
%                                             -diag(obj.d_S_bar_full(:,mean_elevator_W_ind))*Cv*obj.Pe_S-diag(obj.dy_S_mean)*Cu,zero_bc,-zi*obj.kx*diag(U_bar(:,mean_elevator_W_ind))*I_bc*obj.Pe_S+obj.tau*(D2_bc-K2*I_bc)];
%                         otherwise 
%                             error('obj.mean not supported in this case');
%                     end
% 
%                 elseif K2<=0.000001 %For small kx and kz, avoide singularity... use the formulation directly set up kx=kz=0.
%                     %%Also, we need to avoide the inverse of Laplacian as there
%                     %%is a zero eigenvalue of D2_bc, and it will lead to
%                     %%singularity when computing the inverse.
%                     switch obj.mean
%                         case {'kolmogorov','no'}
%                             A_shear= [D2_bc, zero_bc; zero_bc, D2_bc]; 
%                             A(:,:,mean_elevator_W_ind)=[A_shear, obj.Ra_T*[zero_bc; zero_bc], -obj.Ra_S2T*[zero_bc; zero_bc];
%                             -obj.dy_T_mean*I_bc,zero_bc,(D2_bc-K2*I_bc),zero_bc;
%                             -obj.dy_S_mean*I_bc,zero_bc,zero_bc,obj.tau*(D2_bc-K2*I_bc)];
%                         case 'elevator'
%                             
%                             %%Note that this is different...
%                             %%for the elevator mode,,, the kx is actually
%                             %%kz (vertical wavenumber).
%                             %%Thus, kx=kz=0 will get rid of u velocity and
%                             %%it still remain the buoyancy force here.
%                             A_shear= [D2_bc, zero_bc; zero_bc, D2_bc]; 
%                             A(:,:,mean_elevator_W_ind)=[A_shear, obj.Ra_T*[I_bc; zero_bc], -obj.Ra_S2T*[I_bc; zero_bc];
%                             -obj.dy_T_mean*I_bc,zero_bc,D2_bc,zero_bc;
%                             -obj.dy_S_mean*I_bc,zero_bc,zero_bc,obj.tau*D2_bc];
%                             %A11=(D4_bc-2*K2*D2_bc+K2^2*I_bc)+obj.Re*(diag(dd_U_bar(:,mean_elevator_W_ind))*zi*obj.kx*I_bc-zi*obj.kx*diag(U_bar(:,mean_elevator_W_ind))*(D2_bc-K2*I_bc)); %%Orr-Sommerfeld operator
%                             %A21= -zi*obj.kz*diag(d_U_bar(:,mean_elevator_W_ind))*I_bc*obj.Re; %Coulping operator
%                             %A22= -zi*obj.kx*diag(U_bar(:,mean_elevator_W_ind))*I_bc*obj.Re+(D2_bc-K2*I_bc); %Squire operator  
%                             %inv_lap=inv(D2_bc-K2*I_bc);%, zero_bc; zero_bc, I_bc]);
%                             %A_shear= [inv_lap*A11, zero_bc; A21, A22];
% %                             A_shear= [D2_bc, zero_bc; zero_bc, D2_bc]; 
% % 
% %                             Cu=[zi*obj.kx*D1_bc, -zi*obj.kz*I_bc]/K2;
% %                             Cv=[I_bc,zero_bc];
% %                             Bu=[inv_lap*(-zi*obj.kx*D1_bc); zi*obj.kz*I_bc];
% %                             A(:,:,mean_elevator_W_ind)=[A_shear, obj.Ra_T*Bu, -obj.Ra_S2T*Bu;
% %                                             -diag(obj.d_T_bar_full(:,mean_elevator_W_ind))*Cv*obj.Pe_T-diag(obj.dy_T_mean)*Cu,-zi*obj.kx*diag(U_bar(:,mean_elevator_W_ind))*I_bc*obj.Pe_T+(D2_bc-K2*I_bc),zero_bc;
% %                                             -diag(obj.d_S_bar_full(:,mean_elevator_W_ind))*Cv*obj.Pe_S-diag(obj.dy_S_mean)*Cu,zero_bc,-zi*obj.kx*diag(U_bar(:,mean_elevator_W_ind))*I_bc*obj.Pe_S+obj.tau*(D2_bc-K2*I_bc)];
% %                             
%                         otherwise
%                             error('obj.mean not supported in this case');
%                     end
%                 end
%             end
% 
%             %%the matrix for the genearalized eigenvalue problem
%             M=[obj.Re*I_bc, zero_bc,zero_bc,zero_bc;
%                zero_bc, obj.Re*I_bc, zero_bc, zero_bc;
%                zero_bc, zero_bc, obj.Pe_T*I_bc, zero_bc;
%                zero_bc, zero_bc, zero_bc, obj.Pe_S*I_bc];
%             
%             
%             for mean_elevator_W_ind=1:length(obj.mean_elevator_amp_list{2})
% 
%                 [eig_vec,eig_val_mat]=eig(A(:,:,mean_elevator_W_ind),M);% #use linear algebra package to compute eigenvalue
%                 eig_val_mat(isinf(eig_val_mat)|isnan(eig_val_mat)) = -Inf;
% 
%                 [~,eig_val_max_ind]=max(real(diag(eig_val_mat))); %#compute the index of the eigenvalue
%                 obj.eig_val_max(1,mean_elevator_W_ind)=eig_val_mat(eig_val_max_ind,eig_val_max_ind);
%                 obj.eig_vec_max(:,mean_elevator_W_ind)=eig_vec(:,eig_val_max_ind); %#get the corresponding eigen vector
%             end
%         end

            %             eig_val[np.isinf(eig_val)]=-np.inf
            
            %get the eigenvalue with the largest real part and associatd
            %eigenvector
            
%             
%             A=[laplacian^2, 0, obj.Ra_T*laplacian_horizontal, -obj.Ra_S2T*laplacian_horizontal;
%                 0,laplacian,0,0;
%                -obj.dy_T_mean,0, laplacian, 0;
%                -obj.dy_S_mean,0, 0, obj.tau*laplacian];
            
            %%old version with analytical in the vertical..
%             laplacian=-(obj.kx^2+obj.ky^2+obj.kz^2);
%             laplacian_horizontal=-(obj.kx^2+obj.ky^2);
%             A=[laplacian^2, 0, obj.Ra_T*laplacian_horizontal, -obj.Ra_S2T*laplacian_horizontal;
%                 0,laplacian,0,0;
%                -obj.dy_T_mean,0, laplacian, 0;
%                -obj.dy_S_mean,0, 0, obj.tau*laplacian];
%             B=[obj.Re*laplacian, 0,0,0;
%                0, obj.Re, 0, 0;
%                0,0, obj.Pe_T, 0;
%                0,0,0, obj.Pe_S];
%             [eig_vec,eig_val_mat]=eig(A,B);% #use linear algebra package to compute eigenvalue
%             eig_val_mat(isinf(eig_val_mat)|isnan(eig_val_mat)) = -Inf;
% %             eig_val[np.isinf(eig_val)]=-np.inf
%             [~,eig_val_max_ind]=max(real(diag(eig_val_mat))); %#compute the index of the eigenvalue
%             obj.eig_val_max=eig_val_mat(eig_val_max_ind,eig_val_max_ind);
%             obj.eig_vec_max=eig_vec(:,eig_val_max_ind); %#get the corresponding eigen vector
%         

%old version for running in series...           
        function obj=solve_kxkz(obj)
            %%set up result name..
            obj.result_name=['DDC_LST_Re=',num2str(round(obj.Re,2)),...
                    '_Pe_T=',num2str(round(obj.Pe_T,2)),'_Pe_S=',num2str(round(obj.Pe_S,2)),...
                    '_tau=',num2str(obj.tau),'_R_rho_T2S=',num2str(round(obj.R_rho_T2S,2)),...
                    '_Pe=',num2str(round(obj.Pe,2)),'_Ri=',num2str(round(obj.Ri,3)),...
                    '_Ra_T=',num2str(round(obj.Ra_T,2)),'_Ra_S2T',num2str(round(round(obj.Ra_S2T,2))),...
                    '_dy_T_mean=',num2str(obj.dy_T_mean),...
                    '_dy_S_mean=',num2str(obj.dy_S_mean),...
                    '_',obj.flow_sub_double_diffusive_shear_2D,...
                    '_reduced_',obj.shear_Radko2016_reduced,'_mean=',obj.mean,...
                    '_mean_elevator_',obj.mean_elevator_amp_list{1},'=',num2str(round(obj.mean_elevator_amp_list{2}(1),2)) , obj.debug];
                
            switch obj.solve
                case 'finished' %%if finished, just load the data.
                    %%read the data and overwrite the old data obj file
                    %%variable...
                    load([obj.path_data,obj.result_name,'.mat'],'obj');
                    obj.solve='finished'; %This flag from the load data will be other stuff.. move this back to 'finished'
                case 'LST'
                    %go through linear stability analysis throughout the kx
                    %and kz
                    for kx_ind=1:length(obj.kx_list)
                        for kz_ind=1:length(obj.kz_list)
                            obj.kx=obj.kx_list(kx_ind);
                            obj.kz=obj.kz_list(kz_ind);

                            obj=LST(obj);
                            obj.eig_val_max_list{kx_ind,kz_ind}=obj.eig_val_max;
                            obj.eig_vec_max_list{kx_ind,kz_ind}=obj.eig_vec_max;
                            
                            %store the eigenvalue of A+A^T...
                            if  strcmp(obj.mean,'elevator') %strcmp(obj.operator,'v_omega_y') &&
                                obj.eig_val_max_A_AT_time_dependent_list{kx_ind,kz_ind}=obj.eig_val_max_A_AT_time_dependent;
                                obj.eig_vec_max_A_AT_time_dependent_list{kx_ind,kz_ind}=obj.eig_vec_max_A_AT_time_dependent;

                                obj.eig_val_max_A_AT_all_list{kx_ind,kz_ind}=obj.eig_val_max_A_AT_all;
                                obj.eig_vec_max_A_AT_all_list{kx_ind,kz_ind}=obj.eig_vec_max_A_AT_all;
                            end
                        end
                    end
                    save([obj.path_data,obj.result_name,'.mat'],'obj');

                otherwise 
                    error('Wrong flag for solve');
            end
            
        end

        function obj=elevator_lambda_balance_error_compute(obj)
            %%compute the 
            C=obj.elevator_lambda_balance_C;%%The C parameter used to do the growth rate balance as Radko & Smith (2012)
            obj.elevator_lambda_balance_error=max(real(cell2mat(obj.eig_val_max_list)))-C*real(obj.elevator_lambda_max);

        end
        
        function obj=lambda_balance_bisection(obj)
            residue_error=obj.elevator_lambda_balance_bisection(3);
            min_bound=obj.elevator_lambda_balance_bisection(1);
            max_bound=obj.elevator_lambda_balance_bisection(2);
            
            obj_lambda_balance_min=obj; %A local copy...
            obj_lambda_balance_max=obj;
            obj_lambda_balance_mid=obj;

            obj_lambda_balance_min.mean_elevator_amp_list{2}=min_bound;
            obj_lambda_balance_min=obj_lambda_balance_min.solve_kxkz();
            obj_lambda_balance_min=obj_lambda_balance_min.elevator_lambda_balance_error_compute();
            
            obj_lambda_balance_max.mean_elevator_amp_list{2}=max_bound;
            obj_lambda_balance_max=obj_lambda_balance_max.solve_kxkz();
            obj_lambda_balance_max=obj_lambda_balance_max.elevator_lambda_balance_error_compute();
            
            error=1; %%initial error
            bisection_ind=1;
            while abs(error)>residue_error
                %%get the middle object, compute the corresponding growth
                %%rate of the secondary instability and compute the residue
                obj_lambda_balance_mid.mean_elevator_amp_list{2}=(min_bound+max_bound)/2;
                obj_lambda_balance_mid=obj_lambda_balance_mid.solve_kxkz();
                obj_lambda_balance_mid=obj_lambda_balance_mid.elevator_lambda_balance_error_compute();
            
                %%If the middle error is positive, and the min one is
                %%negative, set the max as the middle one...
                %This is the basic algorithm of bisection...
                if obj_lambda_balance_mid.elevator_lambda_balance_error >0 & obj_lambda_balance_min.elevator_lambda_balance_error <0
                    obj_lambda_balance_max=obj_lambda_balance_mid;
                    max_bound=(min_bound+max_bound)/2;
                elseif obj_lambda_balance_mid.elevator_lambda_balance_error <0 & obj_lambda_balance_max.elevator_lambda_balance_error >0
                    obj_lambda_balance_min=obj_lambda_balance_mid;
                    min_bound=(min_bound+max_bound)/2;
                end
                error=(obj_lambda_balance_mid.elevator_lambda_balance_error)
                bisection_ind=bisection_ind+1
            end
            obj=obj_lambda_balance_mid;
        
        end
%         function obj=solve_kxkz(obj)
%             %%set up result name..
%             obj.result_name=['DDC_LST_Re=',num2str(obj.Re),...
%                     '_Pe_T=',num2str(obj.Pe_T),'_Pe_S=',num2str(obj.Pe_S),...
%                     '_tau=',num2str(obj.tau),'_R_rho_T2S=',num2str(obj.R_rho_T2S),...
%                     '_Pe=',num2str(obj.Pe),'_Ri=',num2str(obj.Ri),...
%                     '_Ra_T=',num2str(round(obj.Ra_T)),'_Ra_S2T',num2str(round(obj.Ra_S2T)),...
%                     '_dy_T_mean=',num2str(obj.dy_T_mean),...
%                     '_dy_S_mean=',num2str(obj.dy_S_mean),...
%                     '_',obj.flow_sub_double_diffusive_shear_2D,...
%                     '_reduced_',obj.shear_Radko2016_reduced,'_mean=',obj.mean];
%                 
%             switch obj.solve
%                 case 'finished' %%if finished, just load the data.
%                     %%read the data and overwrite the old data obj file
%                     %%variable...
%                     load([obj.path_data,obj.result_name,'.mat'],'obj');
%                 case 'LST'
%                     
%                     for kx_ind=1:length(obj.kx_list)
%                         
%                     end
%                     save([obj.path_data,obj.result_name,'.mat'],'obj');
% 
%                 otherwise 
%                     error('Wrong flag for solve');
%             end
%             
%         end        
% 
% 
%         function obj=solve_kz(obj)
%             for kz_ind=1:length(obj.kz_list)
%                 obj.kx=obj.kx_list(kx_ind);
%                 obj.kz=obj.kz_list(kz_ind);
% 
%                 obj=LST(obj);
%                 obj.eig_val_max_list(kx_ind,kz_ind,:)=obj.eig_val_max;
%                 obj.eig_vec_max_list{kx_ind,kz_ind}=obj.eig_vec_max;
%             end
%         end
        
%%--------------------
%%auxilliary function that are used for the core solver... 
%%These function should be private and not called by the user...
%%These functions only called by the first part the core solver..
        function obj=grid_diff_fourier(obj)
            %setup the chebyshev grid
            if strcmp(obj.grid_diff{1},'fourier')
                Ny_full=obj.Ny_full;
                Lz=obj.grid_diff{2}(2)-obj.grid_diff{2}(1);
                [y_list_full,D1_full]=fourdif(Ny_full,1);
                [~,D2_full]=fourdif(Ny_full,2);
                [~,D3_full]=fourdif(Ny_full,3);
                [~,D4_full]=fourdif(Ny_full,4);

                %stretch the grid to designed domain. The Fourier grid is
                %uniform, so there is no special trick on this.
                obj.y_list_full=y_list_full/(2*pi)*Lz;
                obj.D1_full=D1_full*2*pi/Lz;
                obj.D2_full=D2_full*(2*pi/Lz)^2;
                obj.D3_full=D3_full*(2*pi/Lz)^3;
                obj.D4_full=D4_full*(2*pi/Lz)^4;
            else
                error('The grid type is not supported right now!');
            end
        end
        
        function obj=mean_profile_elevator_kx(obj)

            obj_elevator=obj;
            obj_elevator.mean='no';
            obj_elevator.kz=0;
            obj_elevator.mean_elevator_amp_list{2}=0;

            if isnumeric(obj.mean_elevator_kx)
                 %setup the wavenumber as the wavenumber provided
                 %by user
                obj_elevator.kx=obj.mean_elevator_kx;
                obj.mean_elevator_kx_local=obj_elevator.kx;
                obj_elevator=obj_elevator.LST();
                obj.elevator_lambda_max=max(real(obj_elevator.eig_val_max));
%                 obj.elevator_eig_vec_max=obj_elevator.eig_vec_max; %%copy this to the obj which will be used later ...
            elseif strcmp(obj.mean_elevator_kx,'max')
                %%For this case, solve the relation in 
                obj_elevator.kx_list=linspace(0.01,2,50);
                obj_elevator=obj_elevator.solve_kxkz();
                [obj.elevator_lambda_max,kx_max_ind]=max(real(cell2mat(obj_elevator.eig_val_max_list)));
                obj_elevator.kx=obj_elevator.kx_list(kx_max_ind);
                obj.mean_elevator_kx_local=obj_elevator.kx;
                obj.mean_elevator_kx_max=obj_elevator.kx;
            elseif strcmp(obj.mean_elevator_kx,'steady')
                %%For this case, solve the relation in 
                obj_elevator.kx_list=linspace(0.01,2,50);
                obj_elevator=obj_elevator.solve_kxkz();
                kx_03_ind=min(find(obj_elevator.kx_list>0.3));
                [obj.elevator_lambda_max,kx_steady_ind]=min(abs(real(cell2mat(obj_elevator.eig_val_max_list(kx_03_ind:end)))));
                obj_elevator.kx=obj_elevator.kx_list(kx_03_ind+kx_steady_ind);
                obj.mean_elevator_kx_local=obj_elevator.kx;
                obj.mean_elevator_kx_steady=obj_elevator.kx;
            end
            obj_elevator=obj_elevator.LST();
            obj.elevator_eig_vec_max=obj_elevator.eig_vec_max; %%put these into the global variable for other functions
            
        end
        
        function obj=mean_profile(obj)
            %setup the mean profile.
            switch obj.mean
                case 'no'
                    %%set the mean flow as zero..
                    obj.U_bar_full=zeros(size(obj.y_list_full));
                    obj.d_U_bar_full=zeros(size(obj.y_list_full));
                    obj.dd_U_bar_full=zeros(size(obj.y_list_full));
                  
                case 'kolmogorov'
                    %%set the mean flow as kolmogorov flow
                    %obj.mean_kolmogorov contain the amplitude, wavenumber,
                    %and phase of the background shear... up to four
                    %different wavenumber
                    syms x

                    A1=obj.mean_kolmogorov(1,1);
                    k1=obj.mean_kolmogorov(1,2);
                    phase_k1=obj.mean_kolmogorov(1,3);

                    A2=obj.mean_kolmogorov(2,1);
                    k2=obj.mean_kolmogorov(2,2);
                    phase_k2=obj.mean_kolmogorov(2,3);

                    A3=obj.mean_kolmogorov(3,1);
                    k3=obj.mean_kolmogorov(3,2);
                    phase_k3=obj.mean_kolmogorov(3,3);

                    A4=obj.mean_kolmogorov(4,1);
                    k4=obj.mean_kolmogorov(4,2);
                    phase_k4=obj.mean_kolmogorov(4,3);

                    U_sym=A1*sin(k1*x+phase_k1)...
                            +A2*sin(k2*x+phase_k2)...
                            +A3*sin(k3*x+phase_k3)...
                            +A4*sin(k4*x+phase_k4);
                    dU_sym=diff(U_sym);
                    ddU_sym=diff(dU_sym);

                    obj.U_bar_full=double(subs(U_sym,obj.y_list_full));
                    obj.d_U_bar_full=double(subs(dU_sym,obj.y_list_full));
                    obj.dd_U_bar_full=double(subs(ddU_sym,obj.y_list_full));      

                case 'elevator'
                    %setup the background state as the elevator mode.
                    %make a copy of this object and then compute the
                    %eigenvalue analysis without shear...
                    %I need to make the following copy to make sure it is
                    %computing the double diffusive convection without
                    %shear...
                    
%                     v_amp=mean(obj_elevator.eig_vec_max(1:obj_elevator.Ny_full)); 
%                     T_amp=mean(obj_elevator.eig_vec_max(2*obj_elevator.Ny_full+1:3*obj_elevator.Ny_full));
%                     S_amp=mean(obj_elevator.eig_vec_max(3*obj_elevator.Ny_full+1:end));
                    
                    for mean_elevator_W_ind=1:length(obj.mean_elevator_amp_list{2})
                        %%change the amplitude of mean_elevator_W_ind;
                        mean_elevator_amp=obj.mean_elevator_amp_list{2}(mean_elevator_W_ind);
                        
                        switch obj.operator
                            case 'v_omega_y'
                                switch obj.mean_elevator_amp_list{1}
                                    %%note that the flag of
                                    %%obj.mean_elevator_amp_list{1} is
                                    %%using the physical coordinate for
                                    %%double-diffusive convection...
                                    %%Inside of the code, the coordinate is
                                    %%still shear flow type..
                                    case 'W'
                                         eig_vec_max_rescale=obj.elevator_eig_vec_max/obj.elevator_eig_vec_max(1)*mean_elevator_amp;
                                    case 'T'
                                         eig_vec_max_rescale=obj.elevator_eig_vec_max/obj.elevator_eig_vec_max(2*obj.Ny_full+1)*mean_elevator_amp;
                                    case 'S'
                                         eig_vec_max_rescale=obj.elevator_eig_vec_max/obj.elevator_eig_vec_max(3*obj.Ny_full+1)*mean_elevator_amp;
                                    otherwise 
                                        error('obj.mean_elevator_amp_list{1} is not supported');
                                end

                                v_amp=mean(eig_vec_max_rescale(1:obj.Ny_full)); 
                                T_amp=mean(eig_vec_max_rescale(2*obj.Ny_full+1:3*obj.Ny_full));
                                S_amp=mean(eig_vec_max_rescale(3*obj.Ny_full+1:4*obj.Ny_full));

                            case 'uvwpTS'
                                switch obj.mean_elevator_amp_list{1}
                                    case 'W'
                                         eig_vec_max_rescale=obj.elevator_eig_vec_max/obj.elevator_eig_vec_max(obj.Ny_full+1)*mean_elevator_amp;
                                    case 'T'
                                         eig_vec_max_rescale=obj.elevator_eig_vec_max/obj.elevator_eig_vec_max(4*obj.Ny_full+1)*mean_elevator_amp;
                                    case 'S'
                                         eig_vec_max_rescale=obj.elevator_eig_vec_max/obj.elevator_eig_vec_max(5*obj.Ny_full+1)*mean_elevator_amp;
                                    otherwise 
                                        error('obj.mean_elevator_amp_list{1} is not supported');
                                end
                                %%Note that this is the v in the code
                                %%coordinate that should be vertical
                                %%motion. W in the obj flag is for user...
                                v_amp=mean(eig_vec_max_rescale(obj.Ny_full+1:2*obj.Ny_full)); 
                                T_amp=mean(eig_vec_max_rescale(4*obj.Ny_full+1:5*obj.Ny_full));
                                S_amp=mean(eig_vec_max_rescale(5*obj.Ny_full+1:6*obj.Ny_full));
                                
                            otherwise
                                error('Wrong obj.operator');
                        end
                        %set up the numerical values of these
                        %mean_profile...
                        syms y 
                        U_sym=real(v_amp)*sin(obj.mean_elevator_kx_local*y);
                        dU_sym=diff(U_sym);
                        ddU_sym=diff(dU_sym);

                        T_sym=real(T_amp)*sin(obj.mean_elevator_kx_local*y);
                        dT_sym=diff(T_sym);

                        S_sym=real(S_amp)*sin(obj.mean_elevator_kx_local*y);
                        dS_sym=diff(S_sym);

                        obj.U_bar_full(:,mean_elevator_W_ind)=double(subs(U_sym,obj.y_list_full));
                        obj.d_U_bar_full(:,mean_elevator_W_ind)=double(subs(dU_sym,obj.y_list_full));
                        obj.dd_U_bar_full(:,mean_elevator_W_ind)=double(subs(ddU_sym,obj.y_list_full));      

                        %obj.T_bar_full(:,mean_elevator_W_ind)=double(subs(T_sym,y_list_full));
                        obj.d_T_bar_full(:,mean_elevator_W_ind)=double(subs(dT_sym,obj.y_list_full));

                        %obj.S_bar_full(:,mean_elevator_W_ind)=double(subs(S_sym,y_list_full));
                        obj.d_S_bar_full(:,mean_elevator_W_ind)=double(subs(dS_sym,obj.y_list_full));

                    end
                    
                    
                otherwise
                    error('The mean flow is not supported right now!');
            end
        end
        
        
%%--------------------------------
%%Third part, the function for the post-processing...
%%These post processing code can be called by the user..
        function obj=convert_IFSC_unit_tuS(obj)
            %%Right now just convert the result of the growth rate... some
            %%other eigenvector of velocity and salinity needs to be done..
            
            switch obj.flow_sub_double_diffusive_shear_2D
                case 'primitive_Radko2013'
                    obj.eig_val_max=obj.eig_val_max/obj.tau_scaling;
                    obj.eig_val_max_list=obj.eig_val_max_list/obj.tau_scaling;
                case {'primitive_IFSC_unit_tuS','IFSC','MRBC','Stokes','shear_Radko2016'}
                    %%need to do nothing 
                otherwise
                    error('Wrong obj.flow_sub_double_diffusive_shear_2D');
            end
            
        end
                


        function obj=post_eig_kxkz_contour(obj)
            %%validation against figure 7 of Radko (2016) Radko T. Thermohaline layering in dynamically and diffusively stable shear flows. Journal of Fluid Mechanics. 2016 Oct;805:147-70.
                
                %%if not comparing,, just plot the contour...
                data{1}.x=obj.kx_list;
                data{1}.y=obj.kz_list;
                data{1}.z=real(obj.eig_val_max_list)';
                if all(all(real(data{1}.z)<0))
                    %there is no instability in this model. do nothing...
                    plot_config.zlim_list=[0,min(min(data{1}.z)),max(max(data{1}.z))];
                else
                    data{1}.z(find(data{1}.z<0))=NaN;
                    data{1}.z(find(data{1}.z>10^4))=NaN;
                    plot_config.zlim_list=[1,0,max(max(data{1}.z))];
                end

                plot_config.label_list={1,'$k_x$','$k_y$'}; %%note spanwise is ky in plotting
                plot_config.Markerindex=3;
                plot_config.user_color_style_marker_list={'k^'};
                
                %add data from Radko 2016
                if obj.post_eig_Radko2016
                    %call the function to get the wavenumber
                    wavenumber_neutral=obj.get_wavenumber_neutral_Radko2016();

                %select different set up the comparing data as digitized
                %from Radko (2016) and the tick range as his paper...
                    if obj.Pe== 10^(4) & obj.Re~=0 %%This is figure 7(a) case
                        data{2}.x=wavenumber_neutral.figure7a(:,1);
                        data{2}.y=wavenumber_neutral.figure7a(:,2);
%                         Pe=10000;
                        plot_config.ztick_list=[1,0,10^(-4),2*10^(-4),3*10^(-4),4*10^(-4),5*10^(-4),6*10^(-4)];
                        
                    elseif obj.Pe==100 & obj.Ri==10  %%This is figure 7(c) case
                        data{2}.x=wavenumber_neutral.figure7c(:,1);
                        data{2}.y=wavenumber_neutral.figure7c(:,2);
%                         Pe=100;
                        plot_config.ztick_list=[1,0,0.005,0.01,0.015,0.02,0.025];
                        
                    elseif obj.Pe==100 & obj.Ri==1 %%This is figure 7(b) case
                        data{2}.x=wavenumber_neutral.figure7b(:,1);
                        data{2}.y=wavenumber_neutral.figure7b(:,2);
%                         Pe=100;
                        plot_config.ztick_list=[1,0,0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05];
                %         plot_config.ztick_list=[1,0,];
                    end
                    obj.result_name=[obj.result_name,'_Kolmogorov_Radko2016'];
                end
                plot_config.name=[obj.path_fig,obj.result_name,'.png'];
                plot_config.print_size=[1,1300,1100];
                plot_config.xlim_list=[1,min(obj.kx_list),max(obj.kx_list)];
                plot_config.ylim_list=[1,min(obj.kz_list),max(obj.kz_list)];
                plot_contour(data,plot_config);

        end
        
        function obj=elevator_flux(obj)
            %%use equation (3.12) of Radko & Smith (2012) to compute the
            %%temperature and salinity flux...
            
            if strcmp(obj.mean_elevator_amp_list{1},'T')
                obj.F_T=1/2*obj.mean_elevator_amp_list{2}.^2*(obj.elevator_lambda_max+obj.mean_elevator_kx_max^2);
                obj.F_S=1/2/obj.R_rho_T2S*obj.mean_elevator_amp_list{2}.^2*(obj.elevator_lambda_max+obj.mean_elevator_kx_max^2)^2/(obj.elevator_lambda_max+obj.tau*obj.mean_elevator_kx_max^2);
            else
                error('The flux computation is based on amplitude of temperature');
            end
            
        end
        
        
        function obj=post_eig_kx(obj)
             %%if not comparing,, just plot the contour...
            data{1}.x=obj.kx_list;
            data{1}.y=real(cell2mat(obj.eig_val_max_list));
            plot_config.Markerindex=3;
            plot_config.user_color_style_marker_list={'k-','b^'};
            if strcmp(obj.mean,'elevator')
                plot_config.label_list={1,'$k_z$','$\lambda$'}; %%note spanwise is ky in plotting
            else
                plot_config.label_list={1,'$k_x$','$\lambda$'}; %%note spanwise is ky in plotting
            end
            %add data from Radko 2016
            if obj.post_eig_Holyer1984
                %call the function to get the wavenumber
                elevator_Holyer1984=get_elevator_Holyer1984(obj);
            %select different set up the comparing data as digitized
            %from Radko (2016) and the tick range as his paper...
                if obj.Pr==10 & obj.tau==0.01 & obj.mean_elevator_amp_list{2}==4 %%This is figure 7(a) case
                    data{2}.x=elevator_Holyer1984.figure3(:,1);
                    data{2}.y=elevator_Holyer1984.figure3(:,2);
%                         Pe=10000;
                    plot_config.ytick_list=[1,0,0.1,0.2,0.3,0.4,0.5];
                    plot_config.ylim_list=[1,0,0.5];
                elseif obj.Pr==1000 & obj.tau==1/3 & obj.mean_elevator_amp_list{2}==100  %%This is figure 7(c) case
                    data{2}.x=elevator_Holyer1984.figure6(:,1);
                    data{2}.y=elevator_Holyer1984.figure6(:,2);
%                         Pe=100;
                    plot_config.ytick_list=[1,0,1,2,3,4];
                    plot_config.ylim_list=[1,0,4];
                end
                obj.result_name=[obj.result_name,'_elevator_Holyer1984'];
            end
            plot_config.name=[obj.path_fig,obj.result_name,'.png'];
            plot_config.print_size=[1,1300,1100];
            plot_config.xlim_list=[1,min(obj.kx_list),max(obj.kx_list)];
%                 plot_config.ylim_list=[1,min(obj.kz_list),max(obj.kz_list)];
            plot_line(data,plot_config);

        end
        
        
        function obj=post_eigenvector(obj)
                %%Then plot the eigenvector of the most unstable mode...
                %remove the negative kx part...
                data{1}.z=real(cell2mat(obj.eig_val_max_list))';
                data{1}.z(find(data{1}.z>10^4))=NaN;

                data{1}.z(:,find(obj.kx_list<0))=-Inf; %%set the negative wavenumber part as zero..
%                 kz0_ind=find(obj.kz_list==0);
                [val,kx_ind]=max(max(data{1}.z)); %%only concern with the zero kz...
                [kz_ind,kx_ind]=find(data{1}.z==val); %%Note that this value has been transposed
                kx_max=obj.kx_list(kx_ind);%get kx and kz that associated with the largest growth rate
                kz_max=obj.kz_list(kz_ind); 
                
                Ny=length(obj.y_list_full);

                %%Note that Pe is because this is in the non-dimensional system
                %%different from Radko (2016)
                %get the eigenvector and also compute eigenvector of u
                %based on the output matrix. This formulation should be
                %general for 3D.
                switch obj.operator
                    case 'v_omega_y'
                        
                        obj.eig_vec_kx_kz.u=1i*kx_max/(kx_max^2+kz_max^2)*obj.D1_full*squeeze(obj.eig_vec_max_list{kx_ind,kz_ind}(1:Ny))...
                            -1i*kz_max/(kx_max^2+kz_max^2)*squeeze(obj.eig_vec_max_list{kx_ind,kz_ind}(Ny+1:2*Ny));
                        %%Note that here I am using w as the vertical velocity.
                        obj.eig_vec_kx_kz.v=squeeze(obj.eig_vec_max_list{kx_ind,kz_ind}(1:Ny));
                        obj.eig_vec_kx_kz.T=squeeze(obj.eig_vec_max_list{kx_ind,kz_ind}(2*Ny+1:3*Ny));
                        obj.eig_vec_kx_kz.S=squeeze(obj.eig_vec_max_list{kx_ind,kz_ind}(3*Ny+1:4*Ny)); %%This 2 is due to the density ratio... that is 1/2 in my system...
                        obj.eig_vec_kx_kz.omega_z=1i*kx_max*obj.eig_vec_kx_kz.v-obj.D1_full*obj.eig_vec_kx_kz.u;
                    case 'uvwpTS'
                        obj.eig_vec_kx_kz.u=squeeze(obj.eig_vec_max_list{kx_ind,kz_ind}(1:Ny));
                        obj.eig_vec_kx_kz.v=squeeze(obj.eig_vec_max_list{kx_ind,kz_ind}(Ny+1:2*Ny));
                        obj.eig_vec_kx_kz.T=squeeze(obj.eig_vec_max_list{kx_ind,kz_ind}(4*Ny+1:5*Ny));
                        obj.eig_vec_kx_kz.S=squeeze(obj.eig_vec_max_list{kx_ind,kz_ind}(5*Ny+1:6*Ny));
                        obj.eig_vec_kx_kz.omega_z=1i*kx_max*obj.eig_vec_kx_kz.v-obj.D1_full*obj.eig_vec_kx_kz.u;

                    otherwise
                        error('obj.operator is not supported');
                        
                end
                
                %scale the eigenvector in the same unit as Radko... Right
                %now I just need to rescale the temperature and salinity...
                if obj.post_eigenvector_Radko2016_unit
                    obj.eig_vec_kx_kz.u=obj.eig_vec_kx_kz.u;
                    obj.eig_vec_kx_kz.v=obj.eig_vec_kx_kz.v;
                    obj.eig_vec_kx_kz.T=obj.eig_vec_kx_kz.T*obj.Pe;
                    obj.eig_vec_kx_kz.S=obj.eig_vec_kx_kz.S*obj.Pe/obj.R_rho_T2S;

                    if obj.Ra_T> 10^(8) %%This is figure 7(a) case
                        eig_vec_tick_list.u=[1,-2,-1,0,1,2];
                        eig_vec_tick_list.v=[1,-6,-4,-2,0,2,4,6];
                        eig_vec_tick_list.T=[1,-1,-0.5,0,0.5,1];
                        eig_vec_tick_list.S=[1,-1.5,-1,-0.5,0,0.5,1,1.5];

                    elseif obj.Ra_T>10^5 %%This is figure 7(c) case
                        eig_vec_tick_list.u=[1,-1,-0.5,0,0.5,1];
                        eig_vec_tick_list.v=[1,-0.04,-0.02,0,0.02,0.04];
                        eig_vec_tick_list.T=[1,-1,-0.5,0,0.5,1];
                        eig_vec_tick_list.S=[1,-1.5,-1,-0.5,0,0.5,1,1.5];

                    else
                        eig_vec_tick_list.u=0;
                        eig_vec_tick_list.v=0;
                        eig_vec_tick_list.T=0;
                        eig_vec_tick_list.S=0;
                    end
                    eig_vec_tick_list.omega_z=0;
                end
                
                
                T_mag=max(abs(obj.eig_vec_kx_kz.T));

                variable_list={'u','v','T','S','omega_z'};
                plot_config.Markerindex=3;
                plot_config.user_color_style_marker_list={'k-','b-.','r--'};
                x_list=linspace(0,2*pi,100);
                
                %plot the eigenvector of these four components...
                for variable_ind=1:length(variable_list)
                   clear data plot_config;
                   switch obj.mean
                       case {'kolmogorov','no'}
                           data{1}.x=x_list;
                           data{1}.y=obj.y_list_full;
                           data{1}.z=real(obj.eig_vec_kx_kz.(variable_list{variable_ind})*exp(1i*x_list))/T_mag;
                           if strcmp(variable_list{variable_ind},'omega_z')
                               data{2}.x=x_list;
                               data{2}.y=obj.y_list_full;
                               data{2}.u=real(obj.eig_vec_kx_kz.u*exp(1i*x_list))';
                               data{2}.v=real(obj.eig_vec_kx_kz.v*exp(1i*x_list))';
                               plot_config.panel_num=2;
                           end
                           
                           plot_config.label_list={1,'$k_x x$','$z$'};
                           plot_config.ylim_list=[1,min(obj.y_list_full),max(obj.y_list_full)];
                           plot_config.xlim_list=[1,0,2*pi];
                           plot_config.xtick_list=[1,0,pi/2,pi,3*pi/2,2*pi];
                           plot_config.xticklabels_list={1,'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'};

                       case 'elevator'
                           data{1}.x=obj.y_list_full;
                           data{1}.y=x_list;
                           data{1}.z=real(obj.eig_vec_kx_kz.(variable_list{variable_ind})*exp(1i*x_list))'/T_mag;
                           if strcmp(variable_list{variable_ind},'omega_z')
                               %%Note that the coordinate of elevator mode
                               %%is differernt...
                               data{2}.x=obj.y_list_full;
                               data{2}.y=x_list;
                               data{2}.u=real(obj.eig_vec_kx_kz.v*exp(1i*x_list))'/T_mag;
                               data{2}.v=real(obj.eig_vec_kx_kz.u*exp(1i*x_list))'/T_mag;
                               plot_config.panel_num=2;
                           end
                           plot_config.label_list={1,'$x$','$k_z z$'};
                           plot_config.xlim_list=[1,min(obj.y_list_full),max(obj.y_list_full)];
                           plot_config.ylim_list=[1,0,2*pi];
                           plot_config.ytick_list=[1,0,pi/2,pi,3*pi/2,2*pi];
                           plot_config.yticklabels_list={1,'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'};

                   end
                   plot_config.colormap='jet';
                   plot_config.ztick_list=eig_vec_tick_list.(variable_list{variable_ind});
                   plot_config.name=[obj.path_fig,obj.result_name,'_primitive_eigenvector_contour_',variable_list{variable_ind},'.png'];
                   plot_config.print_size=[1,1400,1000];
                   plot_contour(data,plot_config);
                end    
            
        end
        
        function obj=post_super_exponential(obj)
            %%Compute the super-exponential results in 
            %%
            lnln_super_exponential=get_lnln_super_exponential(obj);
            
            if obj.Pr==7 && obj.tau==1/24
                switch obj.R_rho_T2S
                    case 2
                        %%This is the 2D results, figure 
                        data{1}.x=lnln_super_exponential.figure3_2D_lower(:,1);
                        data{1}.y=lnln_super_exponential.figure3_2D_lower(:,2);
                        data{2}.x=lnln_super_exponential.figure3_2D_upper(:,1);
                        data{2}.y=lnln_super_exponential.figure3_2D_upper(:,2);
                        lambda=0.251;
                        
                        data{3}.x=data{1}.x;
                        data{3}.y=data{1}.y(1)+(data{3}.x-data{3}.x(1))*lambda;
                        data{4}.x=data{2}.x;
                        data{4}.y=data{2}.y(1)+(data{4}.x-data{4}.x(1))*lambda;
                        
                        plot_config.label_list={1,'$t$','$G_m(t)=lnln\frac{\|T(t)\|}{\|T(t-\delta t)\|}$'};
                        plot_config.Markerindex=3;
                        plot_config.user_color_style_marker_list={'k*','b^','k-','b-.'};
                        plot_config.name=[obj.path_fig,obj.result_name,'_stern2005_figure3_2D.png'];
                        plot_line(data,plot_config);
                    case 1.5
                        data{1}.x=lnln_super_exponential.figure4_3D_lower(:,1);
                        data{1}.y=lnln_super_exponential.figure4_3D_lower(:,2);
                        data{2}.x=lnln_super_exponential.figure4_3D_upper(:,1);
                        data{2}.y=lnln_super_exponential.figure4_3D_upper(:,2);
                        lambda=0.372;
                        data{3}.x=data{1}.x;
                        data{3}.y=data{1}.y(1)+(data{3}.x-data{3}.x(1))*lambda;
                        data{4}.x=data{2}.x;
                        data{4}.y=data{2}.y(1)+(data{4}.x-data{4}.x(1))*lambda;
                        plot_config.Markerindex=3;
                        plot_config.user_color_style_marker_list={'k*','b^','k-','b-.'};
                        plot_config.label_list={1,'$t$','$G_m(t)=lnln \frac{\|T(t)\|}{\|T(t-\delta t)\|}$'};
                        plot_config.Markerindex=3;
                        plot_config.name=[obj.path_fig,obj.result_name,'_stern2005_figure4_3D.png'];
                        plot_line(data,plot_config);

                    otherwise
                        error('R_rho_T2S is not supported');
                end
                
            else
                error('The Prandtl and tau are not correct');
            end
            
            
        end
        
        
 %%------------------------
 %%The Fourth part, storing the function to obtain the data digitized from
 %%other guys paper for validation. %%These function should be not be
 %%called by the user and they can only be called by the third part of the
 %%function...
        
        function wavenumber_neutral=get_wavenumber_neutral_Radko2016(obj)
            %Data digitized form Radko 2016, the zero growth rate boundary 
            %figure number is indicatd in varaible name.
            %Radko T. Thermohaline layering in dynamically and diffusively stable shear flows. Journal of Fluid Mechanics. 2016 Oct;805:147-70.
            wavenumber_neutral.figure7a=[
                    0.27972	0.0210896
                    0.545455	0.161687
                    0.993007	0.358524
                    1.32867	0.541301
                    1.55245	0.653779
                    1.77622	0.794376
                    2.06993	0.949033
                    2.43357	1.07557
                    2.55944	1.07557
                    2.85315	0.991213
                    3.02098	0.822496
                    3.14685	0.55536
                    3.23077	0.274165
                    3.25874	-0.0773286
                    3.14685	-0.541301
                    3.03497	-0.766257
                    2.83916	-0.977153
                    2.64336	-1.04745
                    2.40559	-1.04745
                    2.15385	-0.977153
                    1.88811	-0.836555
                    1.69231	-0.724077
                    1.3986	-0.56942
                    1.1049	-0.400703
                    0.867133	-0.288225
                    0.615385	-0.189807
                    0.433566	-0.105448
                    ];

                    wavenumber_neutral.figure7b=[0.00497278	0.145623
                    0.00821258	0.285979
                    0.0150834	0.342132
                    0.0254446	0.394781
                    0.0357753	0.464975
                    0.0496025	0.528157
                    0.0669079	0.594854
                    0.0877099	0.65454
                    0.108536	0.70019
                    0.16769	0.745908
                    0.206054	0.724923
                    0.23399	0.689884
                    0.261944	0.644319
                    0.279451	0.595227
                    0.29348	0.54262
                    0.311005	0.483002
                    0.325071	0.409342
                    0.339185	0.307613
                    0.356741	0.23045
                    0.363844	0.15327
                    0.367482	0.065557
                    0.367781	-0.106372
                    0.36102	-0.225682
                    0.343758	-0.316941
                    0.333427	-0.387134
                    0.305767	-0.50999
                    0.267641	-0.625846
                    0.23292	-0.696083
                    0.173766	-0.7418
                    0.114453	-0.69629
                    0.0899716	-0.643702
                    0.0654775	-0.584096
                    0.0444127	-0.492905
                    0.0268322	-0.401707
                    0.0197903	-0.359614
                    0.0127361	-0.310504
                    0.00916007	-0.257879
                    0.00556573	-0.194727
                    ];

                    wavenumber_neutral.figure7c=[0.00251255	0.0577977
                    0.00762002	0.126108
                    0.0179845	0.176914
                    0.0248982	0.20845
                    0.0335632	0.234735
                    0.0474547	0.261029
                    0.0700569	0.287338
                    0.0979161	0.296143
                    0.131035	0.285694
                    0.148481	0.271714
                    0.171169	0.248986
                    0.192124	0.221002
                    0.207849	0.19476
                    0.220093	0.16676
                    0.227098	0.145757
                    0.235855	0.119502
                    0.241121	0.0967445
                    0.246393	0.0704841
                    0.249929	0.0407179
                    0.251735	0.00394349
                    0.250069	-0.0398423
                    0.241453	-0.0941481
                    0.232803	-0.129189
                    0.222402	-0.15898
                    0.20852	-0.190528
                    0.194626	-0.215071
                    0.177247	-0.23962
                    0.159859	-0.258914
                    0.140722	-0.27471
                    0.119838	-0.287006
                    0.0971991	-0.292299
                    0.0658159	-0.278344
                    0.0483578	-0.257358
                    0.0396165	-0.239861
                    0.0326174	-0.22236
                    0.0238578	-0.194354
                    0.0168525	-0.173351
                    0.0115772	-0.145339
                    0.00805937	-0.12608
                    0.0045415	-0.106822
                    0.00273222	-0.0682964
                    ];
        end
        
        function elevator_Holyer1984=get_elevator_Holyer1984(obj)
            %Data digized from the following paper to validate the
            %secondary instability of a steady salt finger. 
            %Holyer JY. The stability of long, steady, two-dimensional salt fingers. Journal of Fluid Mechanics. 1984 Oct;147:169-85.
            %figure number is indicated in their variable name
            elevator_Holyer1984.figure3=[0.00557528	0.0141882
                    0.0204285	0.0558042
                    0.0408747	0.106883
                    0.0659865	0.164586
                    0.0920471	0.219456
                    0.111601	0.258245
                    0.146094	0.315966
                    0.175022	0.356664
                    0.206771	0.395476
                    0.245126	0.423904
                    0.285403	0.440048
                    0.307903	0.443871
                    0.352963	0.43545
                    0.381162	0.42038
                    0.406554	0.403415
                    0.42727	    0.38266
                    0.438569	0.371339
                    0.45932	    0.341133
                    0.474418	0.317532
                    0.489541	0.287315
                    0.502784	0.258039
                    ];
    
            elevator_Holyer1984.figure6=[0.00618316	0.390244
                    0.0147881	0.899729
                    0.0233974	1.34417
                    0.0309319	1.71274
                    0.0395456	2.09214
                    0.0470853	2.38482
                    0.0557055	2.66667
                    0.064598	2.91599
                    0.0748411	3.15447
                    0.0888631	3.40379
                    0.103967	3.6206
                    0.119886	3.77236
                    0.133919	3.85908
                    0.146066	3.8916
                    0.167398	3.8374
                    0.180632	3.75068
                    0.191978	3.65312
                    0.207106	3.5122
                    0.216294	3.39295
                    0.225752	3.26287
                    0.239263	3.07859
                    0.24764	2.95935
                    ];
    
        end
        
        function lambda_balance=get_lambda_balance_F_T_Radko_Smith2012(obj)
            %%Data digitized from Radko & Smith (2012)
            %Radko T, Smith DP. Equilibrium transport in double-diffusive convection. Journal of fluid mechanics. 2012 Feb;692:5-27.
            %Figure 5 of their paper regarding 2D growth rate (lambda)
            %balance based computation assuming frozen (quasi-steady)
            %elevator mode.
            lambda_balance.C1p5=[1.20503	8.82353
                1.31055	6.68449
                1.39799	5.34759
                1.50653	4.54545
                1.59698	4.0107
                1.70553	3.47594
                1.79598	3.20856
                1.89849	2.94118
                1.99799	2.6738
                2.1005	2.40642
                2.19698	2.13904
                2.31156	1.87166
                2.39598	1.87166
                ];
            
            lambda_balance.C2=[1.1206	25.4011
                1.2201	18.7166
                1.31658	14.4385
                1.40402	12.0321
                1.51859	9.89305
                1.6	8.55615
                1.71457	7.48663
                1.79899	6.68449
                1.90452	5.88235
                1.99799	5.34759
                2.11256	4.81283
                2.2	4.54545
                2.2995	4.27807
                2.40201	4.0107
                2.54372	3.47594
                ];
            
            lambda_balance.C3=[1.11156	64.7059
                1.16583	54.2781
                1.23216	44.385
                1.28945	37.7005
                1.36482	31.5508
                1.45528	26.4706
                1.55779	22.4599
                1.6603	19.2513
                1.76884	16.5775
                1.87739	14.4385
                2.01005	12.5668
                2.14271	10.9626
                2.30553	9.62567
                2.47136	8.55615
                ];
            
            lambda_balance.C4=[1.10553	118.984
                1.14774	103.743
                1.19296	90.107
                1.23819	78.6096
                1.29246	66.8449
                1.35879	56.9519
                1.4191	50
                1.50955	42.7807
                1.61809	36.0963
                1.70553	31.5508
                1.80503	27.8075
                1.9407	23.7968
                2.0794	20.5882
                2.20905	17.9144
                2.30854	16.5775
                2.50452	14.1711
                ];
            
        end
        
        function secondary_porous_media=get_secondary_porous_media(obj)
            secondary_porous_media.sigma=[2.03409	0.958915
                    3.93291	0.961373
                    7.99988	1.00811
                    15.4678	0.9665
                    31.4628	1.21196
                    63.9981	1.58924
                    123.74	2.08356
                    251.699	2.8571
                    511.978	3.74652
                    989.908	5.13644
                    2013.56	7.0434
                    3893.21	9.65643
                    7919.15	12.6625
                    16108.3	17.3636
                    32765.6	23.81
                    63352.3	32.6432
                    1.29E+05	44.7624
                    2.49E+05	61.3688
                    5.07E+05	84.1526
                    1.03E+06	115.395
                    2.10E+06	158.237
                    4.27E+06	216.984
                    8.25E+06	297.482
                    1.68E+07	407.926
                    3.24E+07	559.262
                    6.60E+07	766.894
                    1.28E+08	1005.43
                    2.59E+08	1441.75
                    5.28E+08	1890.57
                    1.07E+09	2711
                    2.08E+09	3554.24
                    4.22E+09	4873.78
                    8.16E+09	6681.91
                    ];
            
            secondary_porous_media.c=[31.4628	37.0363
                    60.8332	72.9029
                    123.74	143.504
                    239.251	282.475
                    511.978	556.03
                    989.908	1094.5
                    1913.98	2154.43
                    3893.21	4240.83
                    7919.15	9345.19
                    15311.7	18395.3
                    31145.3	40536.2
                    63352.3	71275.6
                    1.29E+05	1.40E+05
                    2.62E+05	3.09E+05
                    5.07E+05	6.09E+05
                    1.03E+06	1.20E+06
                    1.99E+06	2.36E+06
                    4.05E+06	5.20E+06
                    8.25E+06	1.02E+07
                    1.59E+07	2.01E+07
                    3.08E+07	3.96E+07
                    6.27E+07	8.73E+07
                    1.28E+08	1.72E+08
                    2.59E+08	3.38E+08
                    5.02E+08	7.46E+08
                    1.02E+09	1.47E+09
                    2.08E+09	2.58E+09
                    4.01E+09	5.69E+09
                    8.16E+09	1.12E+10
                    ];
                secondary_porous_media.kx=[34.9263	0.161188
                    68.7505	0.161282
                    135.271	0.158217
                    266.095	0.153685
                    523.091	0.144921
                    1083.53	0.138021
                    2130.01	0.13015
                    4411.1	0.122734
                    8671.37	0.115735
                    17038.6	0.107
                    33494.6	0.100898
                    65829.2	0.0942091
                    1.36E+05	0.0879673
                    2.68E+05	0.0813277
                    5.26E+05	0.0751892
                    1.09E+06	0.0709048
                    2.14E+06	0.0649084
                    4.21E+06	0.0606052
                    8.71E+06	0.0571518
                    1.71E+07	0.0523185
                    3.36E+07	0.0483696
                    6.96E+07	0.0451649
                    1.37E+08	0.0417559
                    2.69E+08	0.0386042
                    5.28E+08	0.0356904
                    1.09E+09	0.0329981
                    2.15E+09	0.0308104
                    4.45E+09	0.0287691
                    8.74E+09	0.0263361
                    ];
                
        end
        
        function lnln_super_exponential=get_lnln_super_exponential(obj)
            lnln_super_exponential.figure3_2D_lower=[9.35357	-0.52723
                9.6085	-0.458547
                10.0273	-0.366894
                10.61	-0.218051
                11.0653	-0.114957
                11.648	0.022479
                12.085	0.125556
                12.4856	0.217192
                12.8862	0.308827
                13.305	0.411887
                13.6328	0.48064
                13.997	0.57224
                14.3429	0.652417
                14.7071	0.732611
                14.9621	0.789888
                ];
            
            lnln_super_exponential.figure3_2D_upper=[7.82398	-0.722601
                8.20637	-0.608169
                9.0258	-0.379253
                9.51745	-0.241904
                9.89985	-0.138878
                10.4461	0.00992978
                10.8467	0.124379
                11.2473	0.227421
                11.6115	0.319022
                12.0303	0.433489
                12.3581	0.536462
                12.868	0.673829
                13.3961	0.82262
                13.8149	0.937086
                14.1426	1.01725
                14.5615	1.12031
                14.9621	1.24616
                ];
            lnln_super_exponential.figure4_3D_lower=[8.80888	0.059633
                9.11776	0.197248
                9.39768	0.334862
                9.55212	0.415138
                9.74517	0.506881
                9.92857	0.587156
                10.0637	0.655963
                10.1699	0.713303
                10.2664	0.759174
                10.4305	0.827982
                10.6332	0.931193
                10.7876	1
                ];
            
            lnln_super_exponential.figure4_3D_upper=[8.33591	-0.112385
                8.63514	0.0252294
                8.83784	0.12844
                9.09846	0.243119
                9.33977	0.369266
                9.53282	0.461009
                9.74517	0.56422
                9.94788	0.667431
                10.2085	0.78211
                10.4208	0.873853
                10.6718	1.01147
                10.778	1.06881
                ];
            
        end 
        
    end
end

