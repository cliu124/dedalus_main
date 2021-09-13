classdef convert
    %CONVERT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %%These are three parameter for the double-diffusive convection
        Pr;
        tau;
        R_rho_T2S;
        
        %%additional two parameter required to specify kolmogorov type
        %%shear, this is used in 
        %%
        %Radko T. Thermohaline layering in dynamically and diffusively stable shear flows. Journal of Fluid Mechanics. 2016 Oct;805:147-70.
        Pe;
        Ri;
        
        %%another way to specify the kolmogorov shear.. 
        %%THis is used in 
        %%Brown JM, Radko T. Initiation of diffusive layering by time-dependent shear. Journal of Fluid Mechanics. 2019 Jan;858:588-608.
        H;
        U0;
        
        %%The equivelant Rayleigh number
        Ra;
        
        %%the parameter that I need to store in the inertial free salt
        %%finger.
        Ra_ratio;
        u_L_IFSC;%%The velocity amplitude in 
        F_sin_IFSC;
        ks;
        
        time_IFSC_over_primitive;
        
        
        
        
    end
    
    methods
        function obj = convert(param)
            obj.Pr=param.Pr;
            obj.tau=param.tau;
            obj.R_rho_T2S=param.R_rho_T2S;
            obj.Pe=param.Pe;
            obj.Ri=param.Ri;
            
            obj.Ra=obj.Pe^2*4*pi*obj.Ri/obj.Pr/(1/obj.R_rho_T2S-1);
            obj.U0=obj.Ra^(-1/4)*obj.Pe;
            obj.H=obj.Ra^(1/4);
            obj.ks=2*pi/obj.H;
            obj.Ra_ratio=1/(obj.tau*obj.R_rho_T2S);
            obj.u_L_IFSC=obj.U0/obj.tau;
            obj.F_sin_IFSC=1/obj.tau*(2*pi)^2*obj.Pe*obj.Ra^(-3/4);
            obj.time_IFSC_over_primitive=obj.tau*obj.Ra^(1/2)/obj.Pe;
        end
        
    end
end

