classdef DDC_LST
    %LST_DOUBLE_DIFFUSIVE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        %%parameter for primitive
        Pr=1;
        tau=1;
        R_rho_T2S=1;
        
        %%parameter for IFSC
        Ra_ratio=1;
        
        %%additional parameter for MRBC
        Sc=1;
        
        %parameter for shear as Radko (2016)
        Pe=1;%real Peclet number
        Ri=1;
        
        %Parameter for unified shear flow parameter
        Re=1;
        Pe_T=1;
        Pe_S=1;
        Ra_T=1;
        Ra_S2T=1;
        
        tau_scaling=1;
        
        kx_list=1;
        ky_list=0; %%spanwise wavenumber
        kz_list=0; %%vertical wavenumber. This is in default setting up as zero....
    
        kx=1;
        ky=0;
        kz=0;
        
        eig_val_max=0;
        eig_vec_max=0;
        
        eig_val_max_list=0;
        
        dy_T_mean=1;
        dy_S_mean=1;
        
        flow_sub_double_diffusive_shear_2D='primitive_Radko2013'
    end
    
    methods
        function obj = DDC_LST(flow_sub_double_diffusive_shear_2D)
            %LST_DOUBLE_DIFFUSIVE Construct an instance of this class
            %   Detailed explanation goes here
            obj.flow_sub_double_diffusive_shear_2D=flow_sub_double_diffusive_shear_2D;
        end
        
        function obj = convert_shear(obj)
            switch obj.flow_sub_double_diffusive_shear_2D
                case 'primitive_Radko2013'
                    obj.Re=1/obj.Pr;
                    obj.Pe_T=1;
                    obj.Pe_S=1;
                    obj.tau=obj.tau; %#Set this as zero if remove salinity diffusivity
                    obj.Ra_T=1;
                    obj.Ra_S2T=1/obj.R_rho_T2S;
                    obj.tau_scaling=obj.tau

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
                    obj.Re=0;
                    obj.Pe_T=0;
                    obj.Pe_S=1;
                    obj.tau=1; %#Set this as zero if remove salinity diffusivity
                    obj.Ra_T=1;
                    obj.Ra_S2T=obj.Ra_ratio;
                    obj.tau_scaling=1;
                case 'MRBC'
                
                    %#map to the extended parameter in double_diffusive_shear_2D
                    obj.Re=1/obj.Sc;
                    obj.Pe_T=0;
                    obj.Pe_S=1;
                    obj.tau=1; %#Set this as zero if remove salinity diffusivity
                    obj.Ra_T=1;
                    obj.Ra_S2T=obj.Ra_ratio;
                    obj.tau_scaling=1;
                case 'shear_Radko2016'
               
                    %#map to the extended parameter in double_diffusive_shear_2D
                    obj.Re=obj.Pe/obj.Pr;
                    obj.Pe_T=Pe;
                    obj.Pe_S=Pe;
                    obj.tau=tau;
                    obj.Ra_T=4*pi*obj.Ri/(1/obj.R_rho_T2S-1)*obj.Pr/obj.Pe^2;
                    obj.Ra_S2T=obj.Ra_T/obj.R_rho_T2S;
                otherwise 
                    error('Wrong obj.flow_sub_double_diffusive_shear_2D');
            end
        end
        
        function obj=LST(obj)

            laplacian=-(obj.kx^2+obj.ky^2+obj.kz^2);
            laplacian_horizontal=-(obj.kx^2+obj.ky^2);
            A=[laplacian^2, 0, obj.Ra_T*laplacian_horizontal, -obj.Ra_S2T*laplacian_horizontal;
                0,laplacian,0,0;
               -obj.dy_T_mean,0, laplacian, 0;
               -obj.dy_S_mean,0, 0, obj.tau*laplacian];
            B=[obj.Re*laplacian, 0,0,0;
               0, obj.Re, 0, 0;
               0,0, obj.Pe_T, 0;
               0,0,0, obj.Pe_S];
            [eig_vec,eig_val_mat]=eig(A,B);% #use linear algebra package to compute eigenvalue
            eig_val_mat(isinf(eig_val_mat)|isnan(eig_val_mat)) = -Inf;
%             eig_val[np.isinf(eig_val)]=-np.inf
            [~,eig_val_max_ind]=max(real(diag(eig_val_mat))); %#compute the index of the eigenvalue
            obj.eig_val_max=eig_val_mat(eig_val_max_ind,eig_val_max_ind);
            obj.eig_vec_max=eig_vec(:,eig_val_max_ind); %#get the corresponding eigen vector
        end

            
        function obj=LST_kxkykz(obj)
            for kx_ind=1:length(obj.kx_list)
                for ky_ind=1:length(obj.ky_list)
                    for kz_ind=1:length(obj.kz_list)
                        obj.kx=obj.kx_list(kx_ind);
                        obj.ky=obj.ky_list(ky_ind);
                        obj.kz=obj.kz_list(kz_ind);
                        
                        obj=LST(obj);
                        obj.eig_val_max_list(kx_ind,ky_ind,kz_ind)=obj.eig_val_max;
                    end
                end
            end
        end
        
        function obj=convert_IFSC_unit_tuS(obj)
            %%Right now just convert the result of the growth rate... some
            %%other eigenvector of velocity and salinity needs to be done..
            
            switch obj.flow_sub_double_diffusive_shear_2D
                case 'primitive_Radko2013'
                    obj.eig_val_max=obj.eig_val_max/obj.tau_scaling;
                    obj.eig_val_max_list=obj.eig_val_max_list/obj.tau_scaling;
                case {'primitive_IFSC_unit_tuS','IFSC','MRBC'}
                    %%need to do nothing 
                otherwise
                    error('Wrong obj.flow_sub_double_diffusive_shear_2D');
            end
            
        end
        
    end
end

