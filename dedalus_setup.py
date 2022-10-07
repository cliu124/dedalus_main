import numpy as np
from dedalus import public as de
import time
import pathlib
from scipy import linalg
from dedalus.tools import post
from dedalus.extras import flow_tools
import shutil
import logging
import h5py


class flag(object):
    
    def __init__(self):
        self.Lx=np.pi
        self.Lz=np.pi
        self.Ly=np.pi
        self.Nx=0 ##these default values are zero.. this can avoid some error if I forget to setup the grid point numbers
        self.Ny=0
        self.Nz=0
        
        ##option for self.flow that has been developed:
        #delete the old IFSC_2D_with_shear and IFSC_2D_without_shear option.. these two has been unified together..
        #Update 2021/09/13: also add the option of double_diffusive_2D
        #
        
        #IFSC_2D
        #double_diffusive_2D
        self.flow='not_defined'#
        ##two additional flag that are added using the Harmonic balance formulation for single mode in horizontal
        ##['HB_porous','HB_benard'] 
        
        self.Ra_ratio=1.1# the parameter for IFSC... fix this default value as 1.1... although not used for double diffusive.. Ra_ratio=1 will cause division by zero error.
        
        self.post_store_dt=1
        self.stop_sim_time=1
        self.ks=2*np.pi# parameter for the large scale shear in IFSC with shear
        self.F_sin=0# amplitude for the large scale shear in IFSC with shear
        self.F_sin_2ks=0
        self.F_sin_3ks=0
        self.F_sin_4ks=0
        
        self.phase_2ks=0
        self.phase_3ks=0
        self.phase_4ks=0
        
        self.current_path='./'#This is the current folder path that might need to be specified is run on cluster
    
        #These two values as one corresponds to the salf finger
        #If both of them are set as -1, then we have the diffusive regime
        ##IFSC should also apply in that regime...
        self.dy_T_mean=1
        self.dy_S_mean=1
    
        ##These three parameters are used for the double diffusive convection in primitive variable...
        self.R_rho_T2S=1
        self.tau=1
        self.Pr=1
        
        ##These six parameters are used for the double_diffusive_shear_2D. This is the most general formulation that has six parameters
        self.Re=1 #The Reynolds number appearing in front of the inertial term in momentum
        self.Pe_T=1 #The Peclet number appearing in front of the inertial term in temperature
        self.Pe_S=1 #The Peclet number appearing in front of the inertial term in salinity
        self.Ra_T=1 #The Rayleigh number appearing in front of the temperature term, defined as Ra_T=g\alpha T_z L^4/\nu \kappa_T
        self.Ra_S2T=1 #The Rayleigh number appearing in front of the salinity term, this is defined based salintiy over temperature, thus Ra_T=g\beta S_z L^4/\nu \kappa_T
        #self.tau=1 #This tau is not necessary as it has been defined... #This is the diffusivity ratio, \kappa_S/\kappa_T 
        
        self.A_elevator=0
        self.A_elevator_imag=0
        self.k_elevator=0.5
        self.A_noise=0
        self.A_shear=0        
        
        self.A_secondary_T=0
        self.A_secondary_S=0
        self.A_secondary_w=0
        self.A_secondary_U0=0
        self.k_secondary=0
        
        self.flow_sub_double_diffusive_shear_2D='double_diffusive_2D'
        self.shear_Radko2016_reduced='primitive'
        
        self.lambda_elevator=0
             
        #This is the wavenumber pair for these...
        #Add the governing equations for SMHB... Single mode Harmonic balance...
        
        self.kx=1
        self.ky=1
        #These two parameter add the second harmonic into the govering equations
        self.kx_2=0 
        self.ky_2=0
        self.problem='IVP' #This can be IVP, BVP, EVP depends on the problem you want to solve
        self.bvp_tolerance=1e-11 #This is the tolerance for BVP.
        self.EVP_homogeneous_tolerance=1e-10
        #self.z_bc_T_S_w='dirichlet' #This can be also dirichlet
        #Update 2022/02/25, change the default setting as periodic B.C.
        self.z_bc_u_v_left='periodic' #This can be periodic, dirichlet, or neumann
        self.z_bc_T_left='periodic'
        self.z_bc_S_left='periodic'
        self.z_bc_w_left='periodic'
        self.z_bc_u_v_right='periodic' #This can be periodic, dirichlet, or neumann
        self.z_bc_T_right='periodic'
        self.z_bc_S_right='periodic'
        self.z_bc_w_right='periodic'
        self.z_basis_mode='Fourier'
        
        #flag for the time stepper.. default value
        self.timesteppers='RK443'
        self.analysis=0#Add analysis
        self.solver=0#add solver as a temporal variable
        
        self.initial_dt=0.01#initial time step.
        self.continuation=0 #if yes, use the existing data to continue the next computation
    
        self.continuation_asymmetric=0
    
        self.IBM_A=0
        self.IBM_z0=1/2
        self.IBM_sigma=0.0001
        self.HB_porous_shear_phi=0
    
        #impedance as Hewitt et al. (2014) for the two layer porous medium convection
        self.HB_porous_2_layer_Omega=0
        self.HB_porous_3_layer_Pi=1
        self.HB_porous_3_layer_h=0.1
        
        #wavenumber in the vertical 
        self.initial_kz=2*np.pi
        
        #the default EVP_trivial=1, when EVP_trivial=0, do the linearization of non-trivial state and compute secondary stability results. 
        self.EVP_secondary=0
        
        self.store_variable='all'
        self.nx_trunc_num=0
        
        self.flux_T=0
        self.flux_S=0
    def print_screen(self,logger):
        #print the flag onto the screen
        flag_attrs=vars(self)
        #print(', '.join("%s: %s, \n" % item for item in flag_attrs.items()))
        logger.info(', Attributes: Value,\n,')
        logger.info(', '.join("%s: %s, \n" % item for item in flag_attrs.items()))

    def print_file(self):
        #print the flag onto file
        flag_text=open(self.current_path+'/analysis'+'/flag.txt','w+')
        flag_attrs=vars(self)
        print(', Attributes: 123,\n ------\n-------\n------',file=flag_text)
        print(', test: 123,',file=flag_text)
        print(', '+', '.join("%s: %s, \n" % item for item in flag_attrs.items()),file=flag_text)
        flag_text.close()
        
    def build_domain(self):
        #build domain
        #For these, right now is the Fourier in the veritical
        if self.flow in ['IFSC_2D','double_diffusive_2D','double_diffusive_shear_2D','porous_media_2D']:
            x_basis = de.Fourier('x', self.Nx, interval=(0,self.Lx), dealias=3/2)
            if self.z_bc_u_v_left=='periodic' and self.z_bc_T_left=='periodic' and \
                self.z_bc_S_left=='periodic' and self.z_bc_w_left=='periodic' and \
                self.z_bc_u_v_right=='periodic' and self.z_bc_T_right=='periodic' and \
                self.z_bc_S_right=='periodic' and self.z_bc_w_right=='periodic':
                #if all B.C. are periodic, then just use Fourier mode in vertical
                self.z_basis_mode='Fourier'
                z_basis = de.Fourier('z', self.Nz, interval=(0,self.Lz), dealias=3/2)
                print('Fourier basis in the vertical z direction')
            else: 
                self.z_basis_mode='Chebyshev'
                z_basis = de.Chebyshev('z', self.Nz, interval=(0,self.Lz), dealias=3/2)
                print('Chebyshev basis in the vertical z direction')
            domain = de.Domain([x_basis, z_basis], grid_dtype=np.float64)
        
        #For Harmonic balance method. Chebyshev in the vertical
        #elif self.flow in ['HB_porous','HB_benard','test_periodic','HB_porous_2_layer','HB_porous_3_layer']:
            #if self.z_bc_w =='periodic' and self.z_bc_S =='periodic' and self.z_bc_T=='periodic' and self.z_bc_u_v == 'periodic':
            #    z_basis = de.Fourier('z', self.Nz, interval=(0,self.Lz), dealias=3/2)
            #    #domain = de.Domain([z_basis],grid_dtype=np.float64)
            #else:
            #    z_basis = de.Chebyshev('z', self.Nz, interval=(0, self.Lz), dealias=2)
            #z_basis = de.Chebyshev('z', self.Nz, interval=(0, self.Lz), dealias=2)
            #if self.F_sin==0:
            #    domain = de.Domain([z_basis],grid_dtype=np.float64) 
            #else:
            #    domain = de.Domain([z_basis],grid_dtype=np.complex128) 
        elif self.flow in ['HB_porous','HB_benard','test_periodic',\
                           'HB_porous_2_layer','HB_porous_3_layer',\
                           'HB_porous_shear','HB_benard_shear']:
            if self.problem == 'EVP':
                z_basis = de.Chebyshev('z', self.Nz, interval=(0, self.Lz), dealias=1)
                domain = de.Domain([z_basis],grid_dtype=np.complex128) 
            else:
                z_basis = de.Chebyshev('z', self.Nz, interval=(0, self.Lz), dealias=3/2)
                domain = de.Domain([z_basis],grid_dtype=np.float64) 

        return domain

    def governing_equation(self,domain,solver_in=0):
        #Update 2022/02/27, add one argument as solver_in=0, which in default is not useful, but will be used when solving the eigenvalue problem of a non-trivial state
        
        #This function setup the governing equations
        
        if self.flow in ['IFSC_2D']:
            #This is the 2D inertial free salt finger convection
            problem = de.IVP(domain,variables=['p','u','w','S','T'])
            problem.parameters['Ra_ratio']=self.Ra_ratio
            problem.parameters['dy_T_mean']=self.dy_T_mean
            problem.parameters['dy_S_mean']=self.dy_S_mean
            
            #Update 2021/09/12, change the language to specify the background shear
            #test whether these amplitude of shear is zero....
            #with additional shear needs to modify the x-momentum equation
            if self.F_sin == 0:
                #without shear
                print('without shear')
                problem.add_equation("- (dx(dx(u))+dz(dz(u)) ) +dx(p) = 0", condition="(nx!=0) or (nz!=0)")
            else:
                print('with shear')
                ##specify the background shear... this is kolmogorov type shear... 
                ##This is the amplitude and wavenumber of the fundamental frequency forcing
                problem.parameters['ks']=self.ks
                problem.parameters['F_sin']=self.F_sin
                
                ##Amplitude of the other frequency
                problem.parameters['F_sin_2ks']=self.F_sin_2ks
                problem.parameters['F_sin_3ks']=self.F_sin_3ks
                problem.parameters['F_sin_4ks']=self.F_sin_4ks
                
                ##phase of other frequency
                problem.parameters['phase_2ks']=self.phase_2ks
                problem.parameters['phase_3ks']=self.phase_3ks
                problem.parameters['phase_4ks']=self.phase_4ks
                problem.add_equation("- (dx(dx(u))+dz(dz(u)) ) +dx(p) = F_sin*sin(ks*z)+F_sin_2ks*sin(2*ks*z+phase_2ks)+F_sin_3ks*sin(3*ks*z+phase_3ks)+F_sin_4ks*sin(4*ks*z+phase_4ks)", condition="(nx!=0) or (nz!=0)")

            #These are other governing equations
            problem.add_equation(" - ( dx(dx(w)) + dz(dz(w)) ) + dz(p) -(T-S*Ra_ratio)  =0")            
            problem.add_equation("p=0",condition="(nx==0) and (nz==0)")
            problem.add_equation("u=0",condition="(nx==0) and (nz==0)")
            problem.add_equation("dx(u)+dz(w)=0",condition="(nx!=0) or (nz!=0)")
            problem.add_equation(" - ( dx(dx(T)) + dz(dz(T)) ) + dy_T_mean*w =0")
            problem.add_equation("dt(S) - (dx(dx(S)) + dz(dz(S))) + dy_S_mean*w =-u*dx(S)-w*dz(S) ")

            
        elif self.flow in ['double_diffusive_2D']:
            #double diffusive convection in 2D
            
            problem = de.IVP(domain,variables=['p','u','w','S','T'])
            problem.parameters['R_rho_T2S']=self.R_rho_T2S
            problem.parameters['tau']=self.tau
            problem.parameters['Pr']=self.Pr
            problem.parameters['dy_T_mean']=self.dy_T_mean
            problem.parameters['dy_S_mean']=self.dy_S_mean
            
            #Update 2021/09/12, change the language to specify the background shear
            #test whether these amplitude of shear is zero....
            
            ##Note that this is different from the IFSC,,, here I do not need to constraint that (nx!=0) or (nz!=0) because at nx=nz=0, it is just dt(u)=0, a valid equation.. 
            if self.F_sin == 0:
                print('without shear')
                problem.add_equation("dt(u)- Pr*(dx(dx(u))+dz(dz(u)) ) + Pr*dx(p) = -u*dx(u)-w*dz(u)")
            else:
                print('with shear')
                ##specify the background shear... this is kolmogorov type shear... 
                ##This is the amplitude and wavenumber of the fundamental frequency forcing
                problem.parameters['ks']=self.ks
                problem.parameters['F_sin']=self.F_sin
                
                ##Amplitude of the other frequency
                problem.parameters['F_sin_2ks']=self.F_sin_2ks
                problem.parameters['F_sin_3ks']=self.F_sin_3ks
                problem.parameters['F_sin_4ks']=self.F_sin_4ks
                
                ##phase of other frequency
                problem.parameters['phase_2ks']=self.phase_2ks
                problem.parameters['phase_3ks']=self.phase_3ks
                problem.parameters['phase_4ks']=self.phase_4ks
                problem.add_equation("dt(u) - Pr*(dx(dx(u))+dz(dz(u)) ) +Pr*dx(p) = -u*dx(u)-w*dz(u)+ Pr*(F_sin*sin(ks*z)+F_sin_2ks*sin(2*ks*z+phase_2ks)+F_sin_3ks*sin(3*ks*z+phase_3ks)+F_sin_4ks*sin(4*ks*z+phase_4ks))")

            #problem.add_equation("u=0",condition="(nx==0) and (nz==0)") #Note that for the primitive equation,,, this singularity for u momentum is not there...
            problem.add_equation("dx(u)+dz(w)=0",condition="(nx!=0) or (nz!=0)")
            problem.add_equation("p=0",condition="(nx==0) and (nz==0)")
            problem.add_equation(" dt(w) - Pr*( dx(dx(w)) + dz(dz(w)) ) + Pr*dz(p) -Pr*(T-S/R_rho_T2S)  =-u*dx(w)-w*dz(w)")
            problem.add_equation("dt(S) - tau*(dx(dx(S)) + dz(dz(S))) + dy_S_mean*w =-u*dx(S)-w*dz(S) ")
            problem.add_equation(" dt(T) - ( dx(dx(T)) + dz(dz(T)) ) + dy_T_mean*w =-u*dx(T)-w*dz(T)")
       
        elif self.flow in ['double_diffusive_shear_2D']:
            #This is using the unified formulation, where the velocity and length scale are arbitrary or determined by shear
            if self.flux_T:
                #problem = de.IVP(domain,variables=['p','u','w','S','T','d_u','d_w','d_S','d_T','dy_T_mean_q'])
                #Update 2022/10/07, get rid of dy_T_mean_q as a variable, just a post-processing variable
                problem = de.IVP(domain,variables=['p','u','w','S','T','d_u','d_w','d_S','d_T'])
            else:
                problem = de.IVP(domain,variables=['p','u','w','S','T','d_u','d_w','d_S','d_T'])
            
            problem.parameters['Re']=self.Re
            problem.parameters['Pe_T']=self.Pe_T
            problem.parameters['Pe_S']=self.Pe_S
            problem.parameters['Ra_T']=self.Ra_T
            problem.parameters['Ra_S2T']=self.Ra_S2T
            problem.parameters['tau']=self.tau
            
            problem.parameters['dy_T_mean']=self.dy_T_mean
            problem.parameters['dy_S_mean']=self.dy_S_mean
            
            problem.parameters['Lx']=self.Lx
            problem.parameters['Lz']=self.Lz
            
            #Update 2022/02/25, use the first order formulation. This gives the flexibility of doing Periodic/Dirichlet/Neumann B.C. in the vertical direction by just changing the vertical basis.
            problem.add_equation('dz(u)-d_u=0')
            problem.add_equation('dz(w)-d_w=0')
            problem.add_equation('dz(S)-d_S=0')
            problem.add_equation('dz(T)-d_T=0')
            
            #Update 2021/09/12, change the language to specify the background shear
            #test whether these amplitude of shear is zero....
            
            
            #get the index of the z_basis
            # if self.z_bc_u_v_left=='periodic' and self.z_bc_T_left=='periodic' and \
            #     self.z_bc_S_left=='periodic' and self.z_bc_w_left=='periodic' and \
            #     self.z_bc_u_v_right=='periodic' and self.z_bc_T_right=='periodic' and \
            #     self.z_bc_S_right=='periodic' and self.z_bc_w_right=='periodic':
            #     z_basis_mode='Fourier'
            # else:
            #     z_basis_mode='Chebyshev'
            
            if self.nx_trunc_num>0:
                #Update 2022/06/16, this version make the single mode truncation.
                ##Note that this is different from the IFSC,,, here I do not need to constraint that (nx!=0) or (nz!=0) because at nx=nz=0, it is just dt(u)=0, a valid equation.. 
                #firstly set up the x-momentum equation. if Re=0, then no inertial term
                #Also it needs to distinguish whether we have shear driven by body force or not
                
                nx_trunc_str=str(self.nx_trunc_num)
                if self.F_sin == 0:
                    print('without shear')
                    if self.Re == 0:
                        problem.add_equation("- (dx(dx(u))+dz(d_u) ) + dx(p) = 0",condition="(nx!=0) or (nz!=0)")
                        problem.add_equation("u=0",condition="(nx==0) and (nz==0)")
                    else:
                        problem.add_equation("Re*dt(u)- (dx(dx(u))+dz(d_u)) + dx(p) = Re*(-u*dx(u)-w*d_u)",condition="(nx<=" + nx_trunc_str + ")")
                        problem.add_equation("u=0",condition="(nx>" + nx_trunc_str + ")")
                else:
                    print('with shear')
                    ##specify the background shear... this is kolmogorov type shear... 
                    ##This is the amplitude and wavenumber of the fundamental frequency forcing
                    problem.parameters['ks']=self.ks
                    problem.parameters['F_sin']=self.F_sin
                    
                    ##Amplitude of the other frequency
                    problem.parameters['F_sin_2ks']=self.F_sin_2ks
                    problem.parameters['F_sin_3ks']=self.F_sin_3ks
                    problem.parameters['F_sin_4ks']=self.F_sin_4ks
                    
                    ##phase of other frequency
                    problem.parameters['phase_2ks']=self.phase_2ks
                    problem.parameters['phase_3ks']=self.phase_3ks
                    problem.parameters['phase_4ks']=self.phase_4ks
                    
                    if self.Re == 0:
                        problem.add_equation("- (dx(dx(u))+dz(d_u) ) +dx(p) =  (F_sin*sin(ks*z)+F_sin_2ks*sin(2*ks*z+phase_2ks)+F_sin_3ks*sin(3*ks*z+phase_3ks)+F_sin_4ks*sin(4*ks*z+phase_4ks))",condition="(nx!=0) or (nz!=0)")
                        problem.add_equation("u=0",condition="(nx==0) and (nz==0)")
                    else:
                        problem.add_equation("Re*dt(u) - (dx(dx(u))+dz(d_u) ) +dx(p) = Re*( -u*dx(u)-w*d_u )+ (F_sin*sin(ks*z)+F_sin_2ks*sin(2*ks*z+phase_2ks)+F_sin_3ks*sin(3*ks*z+phase_3ks)+F_sin_4ks*sin(4*ks*z+phase_4ks))")
    
                if self.Re ==0:
                    #no inertial term in the momentum
                    problem.add_equation("- ( dx(dx(w)) + dz(d_w) ) + dz(p) -(Ra_T*T-Ra_S2T*S)  =0")
                else:
                    problem.add_equation("Re*dt(w)- ( dx(dx(w)) + dz(d_w) ) + dz(p) -(Ra_T*T-Ra_S2T*S)  = Re*(-u*dx(w)-w*d_w)",condition="(nx<=" + nx_trunc_str + ")")
                    problem.add_equation("w=0",condition="(nx>" + nx_trunc_str + ")")

    
                #divergence free and pressure gauge
                if self.z_basis_mode=='Fourier':
                    problem.add_equation("dx(u)+d_w=0",condition="(((nx!=0) or (nz!=0)) and (nx<=" + nx_trunc_str + "))")
                    problem.add_equation("p=0",condition="(nx==0) and (nz==0)")
                    problem.add_equation("p=0",condition="(nx>" + nx_trunc_str + ")")
                elif self.z_basis_mode=='Chebyshev':
                    problem.add_equation("dx(u)+d_w=0")
                    #problem.add_equation("dx(u)+d_w=0",condition="(nx!=0)")
                    #problem.add_equation("p=0",condition="(nx==0)")
                
                
                if self.Pe_T == 0:
                    #no inertial term in the temperature
                    problem.add_equation(" - ( dx(dx(T)) + dz(d_T) ) + dy_T_mean*w =0")
                else:
                    problem.add_equation(" Pe_T*dt(T) - ( dx(dx(T)) + dz(d_T) ) + dy_T_mean*w =Pe_T*( -u*dx(T)-w*d_T )",condition="(nx<=" + nx_trunc_str + ")")
                    problem.add_equation(" T=0",condition="(nx>" + nx_trunc_str + ")")

    
                #Add salinity equation
                problem.add_equation("Pe_S*dt(S) - tau*(dx(dx(S)) + dz(d_S)) + dy_S_mean*w =Pe_S*( -u*dx(S)-w*d_S ) ",condition="(nx<=" + nx_trunc_str + ")")
                problem.add_equation("S=0",condition="(nx>" + nx_trunc_str + ")")

            else:
                #This is the branch that has the full equation
                
                ##Note that this is different from the IFSC,,, here I do not need to constraint that (nx!=0) or (nz!=0) because at nx=nz=0, it is just dt(u)=0, a valid equation.. 
                #firstly set up the x-momentum equation. if Re=0, then no inertial term
                #Also it needs to distinguish whether we have shear driven by body force or not
                if self.F_sin == 0:
                    print('without shear')
                    if self.Re == 0:
                        problem.add_equation("- (dx(dx(u))+dz(d_u) ) + dx(p) = 0",condition="(nx!=0) or (nz!=0)")
                        problem.add_equation("u=0",condition="(nx==0) and (nz==0)")
                    else:
                        problem.add_equation("Re*dt(u)- (dx(dx(u))+dz(d_u)) + dx(p) = Re*(-u*dx(u)-w*d_u)")
                else:
                    print('with shear')
                    ##specify the background shear... this is kolmogorov type shear... 
                    ##This is the amplitude and wavenumber of the fundamental frequency forcing
                    problem.parameters['ks']=self.ks
                    problem.parameters['F_sin']=self.F_sin
                    
                    ##Amplitude of the other frequency
                    problem.parameters['F_sin_2ks']=self.F_sin_2ks
                    problem.parameters['F_sin_3ks']=self.F_sin_3ks
                    problem.parameters['F_sin_4ks']=self.F_sin_4ks
                    
                    ##phase of other frequency
                    problem.parameters['phase_2ks']=self.phase_2ks
                    problem.parameters['phase_3ks']=self.phase_3ks
                    problem.parameters['phase_4ks']=self.phase_4ks
                    
                    if self.Re == 0:
                        problem.add_equation("- (dx(dx(u))+dz(d_u) ) +dx(p) =  (F_sin*sin(ks*z)+F_sin_2ks*sin(2*ks*z+phase_2ks)+F_sin_3ks*sin(3*ks*z+phase_3ks)+F_sin_4ks*sin(4*ks*z+phase_4ks))",condition="(nx!=0) or (nz!=0)")
                        problem.add_equation("u=0",condition="(nx==0) and (nz==0)")
                    else:
                        #
                        problem.add_equation("Re*dt(u) - (dx(dx(u))+dz(d_u) ) +dx(p) = Re*( -u*dx(u)-w*d_u )+ (F_sin*sin(ks*z)+F_sin_2ks*sin(2*ks*z+phase_2ks)+F_sin_3ks*sin(3*ks*z+phase_3ks)+F_sin_4ks*sin(4*ks*z+phase_4ks))")
                        
                        #Below is to enforce that large-scale shear has a zero mean value, but for simulation seems not necessary
                        #if self.z_bc_u_v_left=='neumann' and self.z_bc_u_v_left=='neumann':
                        #    #modify the large-scale shear 
                        #    problem.add_equation("Re*dt(u) - (dx(dx(u))+dz(d_u) ) +dx(p) = Re*( -u*dx(u)-w*d_u )+ (F_sin*sin(ks*z)+F_sin_2ks*sin(2*ks*z+phase_2ks)+F_sin_3ks*sin(3*ks*z+phase_3ks)+F_sin_4ks*sin(4*ks*z+phase_4ks))",condition="(nx!=0) or (nz!=0)")
                        #    problem.add_equation("u=0",condition="(nx==0) and (nz==0)")
                        #else:
                        #    problem.add_equation("Re*dt(u) - (dx(dx(u))+dz(d_u) ) +dx(p) = Re*( -u*dx(u)-w*d_u )+ (F_sin*sin(ks*z)+F_sin_2ks*sin(2*ks*z+phase_2ks)+F_sin_3ks*sin(3*ks*z+phase_3ks)+F_sin_4ks*sin(4*ks*z+phase_4ks))")
    
                if self.Re ==0:
                    #no inertial term in the momentum
                    problem.add_equation("- ( dx(dx(w)) + dz(d_w) ) + dz(p) -(Ra_T*T-Ra_S2T*S)  =0")
                else:
                    problem.add_equation("Re*dt(w)- ( dx(dx(w)) + dz(d_w) ) + dz(p) -(Ra_T*T-Ra_S2T*S)  = Re*(-u*dx(w)-w*d_w)")
    
                #divergence free and pressure gauge
                if self.z_basis_mode=='Fourier':
                    problem.add_equation("dx(u)+d_w=0",condition="(nx!=0) or (nz!=0)")
                    problem.add_equation("p=0",condition="(nx==0) and (nz==0)")
                elif self.z_basis_mode=='Chebyshev':
                    problem.add_equation("dx(u)+d_w=0")
                    #problem.add_equation("dx(u)+d_w=0",condition="(nx!=0)")
                    #problem.add_equation("p=0",condition="(nx==0)")
                #add flux feedback 
                if self.flux_T:
                    #problem.add_equation(" Pe_T*dt(T) - ( dx(dx(T)) + dz(d_T) )  =-dy_T_mean_q*w+Pe_T*( -u*dx(T)-w*d_T )")
                    #problem.add_equation(" Pe_T*dt(T) - ( dx(dx(T)) + dz(d_T) )  =-dy_T_mean_q*w+Pe_T*( -u*dx(T)-w*d_T )",condition="(nx!=0) or (nz!=0)")
                    
                    #Update 2022/10/05, try this version to have the time-varying effect of dy_T_mean_q on the right hand size
                    #problem.add_equation(" Pe_T*dt(T) - ( dx(dx(T)) + dz(d_T) )  =-dy_T_mean_q*w+Pe_T*( -u*dx(T)-w*d_T)",condition="(nx!=0) or (nz!=0)")
                    #problem.add_equation("T=0",condition="(nx==0) and (nz==0)")                     
                    #problem.add_equation("dy_T_mean_q=0",condition="(nx!=0) or (nz!=0)")
                    #problem.add_equation("-dy_T_mean_q=1-integ(w*T)/Lx/Lz",condition="(nx==0) and (nz==0)")
                
                    #Update 2022/10/07, substitute dy_T_mean_q into the T equation
                    problem.add_equation(" Pe_T*dt(T) - ( dx(dx(T)) + dz(d_T) )  =(1-integ(w*T)/Lx/Lz)*w+Pe_T*( -u*dx(T)-w*d_T)",condition="(nx!=0) or (nz!=0)")
                    problem.add_equation("T=0",condition="(nx==0) and (nz==0)")                     
                    
                
                else: 
                    if self.Pe_T == 0:
                        #no inertial term in the temperature
                        problem.add_equation(" - ( dx(dx(T)) + dz(d_T) ) + dy_T_mean*w =0")
                    else:
                        problem.add_equation(" Pe_T*dt(T) - ( dx(dx(T)) + dz(d_T) ) + dy_T_mean*w =Pe_T*( -u*dx(T)-w*d_T )")
        
        
    
                #Add salinity equation
                problem.add_equation("Pe_S*dt(S) - tau*(dx(dx(S)) + dz(d_S)) + dy_S_mean*w =Pe_S*( -u*dx(S)-w*d_S ) ")

            #Add B.C. conditions for non-periodic vertical domain
            
            if self.z_bc_w_left=='dirichlet':
                problem.add_bc("left(w)=0")
                print("Dirichlet for w left")
            elif self.z_bc_w_left=='neumann':
                problem.add_bc("left(d_w)=0")
                print("Neumann for w left")
                
            if self.z_bc_w_right=='dirichlet':
                if self.z_bc_w_left=='dirichlet':
                    problem.add_bc("right(w)=0",condition="(nx!=0)")
                    problem.add_equation("right(p)=0",condition="(nx==0)")
                else:
                    problem.add_bc("right(w)=0")
                print("Dirichlet for w right")
            elif self.z_bc_w_right=='neumann':
                problem.add_bc("right(d_w)=0")
                print("Neumann for w right")
            
            if self.z_bc_w_left=='periodic' and self.z_bc_w_right=='periodic':
                print('Periodic for w')
                #For DNS, the periodic domain in the vertical does not need to do anything but just change the vertical basis as the Fourier one.
                #problem.add_bc("left(w_hat)-right(w_hat)=0")
                #problem.add_bc("left(p_hat)-right(p_hat)=0")
                
            if self.z_bc_T_left=='dirichlet':
                problem.add_bc("left(T)=0")
                #problem.add_bc("left(d_T)=0")
                print("Dirichlet for T left")
            elif self.z_bc_T_left=='neumann':
                problem.add_bc("left(d_T)=0")
                #problem.add_bc("left(d_T_0)=0")
                print("Neumann for T left")
                
            if self.z_bc_T_right=='dirichlet':
                problem.add_bc("right(T)=0")
                #problem.add_bc("right(T_0)=0")
                print("Dirichlet for T right")
            elif self.z_bc_T_right=='neumann':
                problem.add_bc("right(d_T)=0")
                #problem.add_bc("right(d_T_0)=0")
                print("Neumann for T right")
            
            if self.z_bc_T_left=='periodic' and self.z_bc_T_right=='periodic':
                #problem.add_bc("left(T_hat)-right(T_hat)=0")
                #problem.add_bc("left(d_T_hat)-right(d_T_hat)=0")
                #problem.add_bc("left(T_0)=0")
                #problem.add_bc("right(T_0)=0")
                #problem.add_bc("left(d_T_0)-right(d_T_0)=0")
                print("Periodic for T")
               
            if self.z_bc_S_left=='dirichlet':
                problem.add_bc("left(S)=0")
                #problem.add_bc("left(S_0)=0")
                print("Dirichlet for S left")
            elif self.z_bc_S_left=='neumann':
                problem.add_bc("left(d_S)=0")
                #problem.add_bc("left(d_S_0)=0")
                print("Neumann for S left")
                
            if self.z_bc_S_right=='dirichlet':
                problem.add_bc("right(S)=0")
                #problem.add_bc("right(S_0)=0")
                print("Dirichlet for S right")
            elif self.z_bc_S_right=='neumann':
                problem.add_bc("right(d_S)=0")
                #problem.add_bc("right(d_S_0)=0")
                print("Neumann for S right")
            
            if self.z_bc_S_left=='periodic' and self.z_bc_S_right=='periodic':
                #problem.add_bc("left(S_hat)-right(S_hat)=0")
                #problem.add_bc("left(d_S_hat)-right(d_S_hat)=0")
                #problem.add_bc("left(S_0)=0")
                #problem.add_bc("right(S_0)=0")
                #problem.add_bc("left(d_S_0)-right(d_S_0)=0")
                print("Periodic for S")
           
            if self.z_bc_u_v_left=='dirichlet':
                problem.add_bc("left(u)=0")
                #problem.add_bc("left(v_tilde)=0")
            elif self.z_bc_u_v_left=='neumann':
                problem.add_bc("left(d_u)=0")
                #problem.add_bc("left(d_v_tilde)=0")
                
            if self.z_bc_u_v_right=='dirichlet':
                problem.add_bc("right(u)=0")
                #problem.add_bc("right(v_tilde)=0")
            elif self.z_bc_u_v_right=='neumann':
                problem.add_bc("right(d_u)=0")
                #problem.add_bc("right(d_v_tilde)=0")
                
            if self.z_bc_u_v_left=='periodic' and self.z_bc_u_v_right=='periodic':
                #problem.add_bc("left(u_tilde)-right(u_tilde)=0")
                #problem.add_bc("left(v_tilde)-right(v_tilde)=0")
                #problem.add_bc("left(d_u_tilde)-right(d_u_tilde)=0")
                #problem.add_bc("left(d_v_tilde)-right(d_v_tilde)=0")
                print("Periodic for u and v")


        elif self.flow =='porous_media_2D':
            #porous media in 2D!!!
            #This has not been fully benchmarked!!!
            problem = de.IVP(domain,variables=['p','u','w','T'])
            problem.parameters['Ra_T']=self.Ra_T
            problem.add_equation(" u + dx(p) = 0", condition="(nx!=0) or (nz!=0)")
            problem.add_equation(" w + dz(p) - Ra_T*T  =0")            
            problem.add_equation("p=0",condition="(nx==0) and (nz==0)")
            problem.add_equation("u=0",condition="(nx==0) and (nz==0)")
            problem.add_equation("dx(u)+dz(w)=0",condition="(nx!=0) or (nz!=0)")
            problem.add_equation("dt(T) - (dx(dx(T)) + dz(dz(T)))  - w =-u*dx(T)-w*dz(T) ")

        elif self.flow == "channel":
            #channel flow
            #This has not been fully benchmarked!!!

            problem.parameters['Re']=self.Re
            problem = de.IVP(domain,variables=['p','u','w','v'])

            problem.add_equation("dx(u) + dy(v)+wz=0")
            problem.add_equation("dt(u) - 1/Re*(dx(dx(u)) + dy(dy(u)) + dz(uz) ) + dx(p) = -u*dx(u) - v*dy(u) - w*uz + 2/Re")
            problem.add_equation("dt(v) - 1/Re*(dx(dx(v)) + dy(dy(v)) + dz(vz) ) + dy(p) = -u*dx(v) - v*dy(v) - w*vz")
            problem.add_equation("dt(w) - 1/Re*(dx(dx(w)) + dy(dy(w)) + dz(wz) ) + dz(p) = -u*dx(w) v*dy(w) - w*wz")
            problem.add_equation("uz - dz(u) =0")
            problem.add_equation("vz - dz(v) =0")
            problem.add_equation("wz - dz(w) =0")
            problem.add_bc("left(u) = 0")
            problem.add_bc("left(v) = 0")
            problem.add_bc("left(w) = 0")
            problem.add_bc("right(u) = 0")
            problem.add_bc("right(v) = 0")
            problem.add_bc("right(w) = 0", condition="(nx !=0) or (ny !=0)")
            problem.add_bc("right(p) = 0", condition="(nx == 0) and (ny == 0)")
        
        elif self.flow =='HB_porous':
            #harmonic balance for the porous media
            #For different problem, we need to claim different dedalus problem.
            if self.kx_2==0 and self.kx_2==0:
                #setup the problem, only one scale
                if self.problem =='BVP':
                    problem = de.NLBVP(domain, variables=[\
                        'w_hat','p_hat','T_hat','d_T_hat','S_hat','d_S_hat', \
                        'T_0','d_T_0','S_0','d_S_0'])
                elif self.problem == 'IVP':
                    problem = de.IVP(domain, variables=[\
                        'w_hat','p_hat','T_hat','d_T_hat','S_hat','d_S_hat', \
                        'T_0','d_T_0','S_0','d_S_0'])
                elif self.problem == 'EVP':
                    problem = de.EVP(domain, variables=[\
                        'w_hat','p_hat','T_hat','d_T_hat','S_hat','d_S_hat', \
                        'T_0','d_T_0','S_0','d_S_0'],eigenvalue='eig_val',tolerance=self.EVP_homogeneous_tolerance)
            else:   
                #setup the problem variable. This has two scales...
                if self.problem =='BVP':
                    problem = de.NLBVP(domain, variables=[\
                        'w_hat','p_hat','T_hat','d_T_hat','S_hat','d_S_hat', \
                        'w_hat_2','p_hat_2','T_hat_2','d_T_hat_2','S_hat_2','d_S_hat_2', \
                        'T_0','d_T_0','S_0','d_S_0'])
                elif self.problem == 'IVP':
                    problem = de.IVP(domain, variables=[\
                        'w_hat','p_hat','T_hat','d_T_hat','S_hat','d_S_hat', \
                        'w_hat_2','p_hat_2','T_hat_2','d_T_hat_2','S_hat_2','d_S_hat_2', \
                        'T_0','d_T_0','S_0','d_S_0'])
                elif self.problem == 'EVP':
                    problem = de.EVP(domain, variables=[\
                        'w_hat','p_hat','T_hat','d_T_hat','S_hat','d_S_hat', \
                        'w_hat_2','p_hat_2','T_hat_2','d_T_hat_2','S_hat_2','d_S_hat_2', \
                        'T_0','d_T_0','S_0','d_S_0'],eigenvalue='eig_val',tolerance=self.EVP_homogeneous_tolerance)
                
                #also get the kx_2 and ky_2 as parameters
                problem.parameters['kx_2']=self.kx_2
                problem.parameters['ky_2']=self.ky_2
            
            problem.parameters['Ra_T'] = self.Ra_T
            problem.parameters['Ra_S2T'] = self.Ra_S2T
            problem.parameters['tau']=self.tau
            problem.parameters['dy_T_mean']=self.dy_T_mean
            problem.parameters['dy_S_mean']=self.dy_S_mean
            problem.parameters['kx']=self.kx
            problem.parameters['ky']=self.ky
            #if self.kx_2!=0 or self.kx_2!=0:
                
            if self.problem =='BVP':
                #variable with _hat is the harmonic term
                #variable with _0 is the horizontal average term
                problem.add_equation('dz(w_hat)-(-(kx*kx+ky*ky)*p_hat)=0')
                problem.add_equation('dz(p_hat)-(-w_hat+Ra_T*T_hat-Ra_S2T*S_hat)=0')
                problem.add_equation('dz(T_hat)-d_T_hat=0')
                problem.add_equation('dz(d_T_hat)-(w_hat*dy_T_mean+(kx*kx+ky*ky)*T_hat)=w_hat*d_T_0')
                problem.add_equation('dz(S_hat)-d_S_hat=0')
                problem.add_equation('dz(d_S_hat)-1/tau*w_hat*dy_S_mean-(kx*kx+ky*ky)*S_hat=1/tau*(w_hat*d_S_0)')   
                problem.add_equation('dz(T_0)-d_T_0=0')
                problem.add_equation('dz(S_0)-d_S_0=0')
                if self.kx_2==0 and self.ky_2==0:
                    problem.add_equation('dz(d_T_0)=-2*(kx*kx+ky*ky)*p_hat*T_hat+2*w_hat*d_T_hat')
                    problem.add_equation('dz(d_S_0)=1/tau*(-2*(kx*kx+ky*ky)*p_hat*S_hat+2*w_hat*d_S_hat)')
                else:
                    problem.add_equation('dz(d_T_0)=-2*(kx*kx+ky*ky)*p_hat*T_hat+2*w_hat*d_T_hat  -2*(kx_2*kx_2+ky_2*ky_2)*p_hat_2*T_hat_2+2*w_hat_2*d_T_hat_2')
                    problem.add_equation('dz(d_S_0)=1/tau*(-2*(kx*kx+ky*ky)*p_hat*S_hat+2*w_hat*d_S_hat -2*(kx_2*kx_2+ky_2*ky_2)*p_hat_2*S_hat_2+2*w_hat_2*d_S_hat_2)')
                    #For the second harmonic
                    problem.add_equation('dz(w_hat_2)-(-(kx_2*kx_2+ky_2*ky_2)*p_hat_2)=0')
                    problem.add_equation('dz(p_hat_2)-(-w_hat_2+Ra_T*T_hat_2-Ra_S2T*S_hat_2)=0')
                    problem.add_equation('dz(T_hat_2)-d_T_hat_2=0')
                    problem.add_equation('dz(d_T_hat_2)-(w_hat_2*dy_T_mean+(kx_2*kx_2+ky_2*ky_2)*T_hat_2)=w_hat_2*d_T_0')
                    problem.add_equation('dz(S_hat_2)-d_S_hat_2=0')
                    problem.add_equation('dz(d_S_hat_2)-1/tau*w_hat_2*dy_S_mean-(kx_2*kx_2+ky_2*ky_2)*S_hat_2=1/tau*(w_hat_2*d_S_0)')   
                
            
            elif self.problem =='IVP':
                #This is the same version, but just add the time dependence term in the temperature, both harmonic and horizontal average term
                problem.add_equation('dz(w_hat)-(-(kx*kx+ky*ky)*p_hat)=0')
                problem.add_equation('dz(p_hat)-(-w_hat+Ra_T*T_hat-Ra_S2T*S_hat)=0')
                problem.add_equation('dz(T_hat)-d_T_hat=0')
                problem.add_equation('-dt(T_hat)+dz(d_T_hat)-(w_hat*dy_T_mean+(kx*kx+ky*ky)*T_hat)=w_hat*d_T_0')
                problem.add_equation('dz(S_hat)-d_S_hat=0')
                problem.add_equation('-1/tau*dt(S_hat)+dz(d_S_hat)-1/tau*w_hat*dy_S_mean-(kx*kx+ky*ky)*S_hat=1/tau*(w_hat*d_S_0)')   
                problem.add_equation('dz(T_0)-d_T_0=0')
                problem.add_equation('dz(S_0)-d_S_0=0')
            
                if self.kx_2==0 and self.ky_2==0:
                    problem.add_equation('-dt(T_0)+dz(d_T_0)=-2*(kx*kx+ky*ky)*p_hat*T_hat+2*w_hat*d_T_hat')
                    problem.add_equation('-1/tau*dt(S_0)+dz(d_S_0)=1/tau*(-2*(kx*kx+ky*ky)*p_hat*S_hat+2*w_hat*d_S_hat)')
                else:
                    #Modify the mean equations, also add the contribution from the second harmonic
                    problem.add_equation('-dt(T_0)+dz(d_T_0)=-2*(kx*kx+ky*ky)*p_hat*T_hat+2*w_hat*d_T_hat-2*(kx_2*kx_2+ky_2*ky_2)*p_hat_2*T_hat_2+2*w_hat_2*d_T_hat_2')
                    problem.add_equation('-1/tau*dt(S_0)+dz(d_S_0)=1/tau*(-2*(kx*kx+ky*ky)*p_hat*S_hat+2*w_hat*d_S_hat -2*(kx_2*kx_2+ky_2*ky_2)*p_hat_2*S_hat_2+2*w_hat_2*d_S_hat_2)')
                    #Add the dynamics of the second harmonic
                    problem.add_equation('dz(w_hat_2)-(-(kx_2*kx_2+ky_2*ky_2)*p_hat_2)=0')
                    problem.add_equation('dz(p_hat_2)-(-w_hat_2+Ra_T*T_hat_2-Ra_S2T*S_hat_2)=0')
                    problem.add_equation('dz(T_hat_2)-d_T_hat_2=0')
                    problem.add_equation('-dt(T_hat_2)+dz(d_T_hat_2)-(w_hat_2*dy_T_mean+(kx_2*kx_2+ky_2*ky_2)*T_hat_2)=w_hat_2*d_T_0_2')
                    problem.add_equation('dz(S_hat_2)-d_S_hat_2=0')
                    problem.add_equation('-1/tau*dt(S_hat_2)+dz(d_S_hat_2)-1/tau*w_hat_2*dy_S_mean-(kx_2*kx_2+ky_2*ky_2)*S_hat_2=1/tau*(w_hat_2*d_S_0)')   
            
            elif self.problem == 'EVP':    
                problem.substitutions['dt(A)'] = "eig_val*A"
                if self.EVP_secondary:
                    state=solver_in.state
                    if self.kx_2==0 and self.ky_2==0:
                        var_list=['w_hat','p_hat','T_hat','d_T_hat','S_hat','d_S_hat', \
                        'T_0','d_T_0','S_0','d_S_0']
                    else:
                        var_list=['w_hat','p_hat','T_hat','d_T_hat','S_hat','d_S_hat', \
                        'w_hat_2','p_hat_2','T_hat_2','d_T_hat_2','S_hat_2','d_S_hat_2', \
                        'T_0','d_T_0','S_0','d_S_0']
                        
                    for varname in var_list:
                        problem.substitutions['{0}_tot'.format(varname)]='{0}_base'.format(varname)+'+'+varname 
                        ncc = domain.new_field(name='{0}_base'.format(varname))
                        ncc['c'] = state[varname]['c']
                        problem.parameters['{0}_base'.format(varname)] = ncc
                    
                    #setup the problem
                    #put total variable on the RHS, but need to substract the linear part of the base state. This will guarantee that when perturbation=0 (without _tot or _base) , RHS is zero (homogeneous). Dedalus will check this
                    #The LHS is not tricky, just the original problem will be OK.
                    #print('Add equation of EVP for HB_porous')
                    problem.add_equation('dz(w_hat)-(-(kx*kx+ky*ky)*p_hat)=0')
                    problem.add_equation('dz(p_hat)-(-w_hat+Ra_T*T_hat-Ra_S2T*S_hat)=0')
                    problem.add_equation('dz(T_hat)-d_T_hat=0')
                    problem.add_equation('-dt(T_hat)+dz(d_T_hat)-(w_hat*dy_T_mean+(kx*kx+ky*ky)*T_hat)=w_hat_tot*d_T_0_tot  -dz(d_T_hat_base)+(w_hat_base*dy_T_mean+(kx*kx+ky*ky)*T_hat_base) ')
                    problem.add_equation('dz(S_hat)-d_S_hat=0')
                    problem.add_equation('-1/tau*dt(S_hat)+dz(d_S_hat)-1/tau*w_hat*dy_S_mean-(kx*kx+ky*ky)*S_hat=1/tau*(w_hat_tot*d_S_0_tot)  -dz(d_S_hat_base)+1/tau*w_hat_base*dy_S_mean+(kx*kx+ky*ky)*S_hat_base ')   
                    problem.add_equation('dz(T_0)-d_T_0=0')
                    problem.add_equation('dz(S_0)-d_S_0=0')
                
                    if self.kx_2==0 and self.ky_2==0:
                        problem.add_equation('-dt(T_0)+dz(d_T_0)=-2*(kx*kx+ky*ky)*p_hat_tot*T_hat_tot+2*w_hat_tot*d_T_hat_tot -dz(d_T_0_base)')
                        problem.add_equation('-1/tau*dt(S_0)+dz(d_S_0)=1/tau*(-2*(kx*kx+ky*ky)*p_hat_tot*S_hat_tot+2*w_hat_tot*d_S_hat_tot) -dz(d_S_0_base)')
                    else:
                        #Modify the mean equations, also add the contribution from the second harmonic
                        problem.add_equation('-dt(T_0)+dz(d_T_0)=-2*(kx*kx+ky*ky)*p_hat_tot*T_hat_tot+2*w_hat_tot*d_T_hat_tot-2*(kx_2*kx_2+ky_2*ky_2)*p_hat_2_tot*T_hat_2_tot+2*w_hat_2_tot*d_T_hat_2_tot -dz(d_T_0_base)')
                        problem.add_equation('-1/tau*dt(S_0)+dz(d_S_0)=1/tau*(-2*(kx*kx+ky*ky)*p_hat_tot*S_hat_tot+2*w_hat_tot*d_S_hat_tot -2*(kx_2*kx_2+ky_2*ky_2)*p_hat_2_tot*S_hat_2_tot+2*w_hat_2_tot*d_S_hat_2_tot) -dz(d_S_0_base)')
                        #Add the dynamics of the second harmonic
                        problem.add_equation('dz(w_hat_2)-(-(kx_2*kx_2+ky_2*ky_2)*p_hat_2)=0')
                        problem.add_equation('dz(p_hat_2)-(-w_hat_2+Ra_T*T_hat_2-Ra_S2T*S_hat_2)=0')
                        problem.add_equation('dz(T_hat_2)-d_T_hat_2=0')
                        problem.add_equation('-dt(T_hat_2)+dz(d_T_hat_2)-(w_hat_2*dy_T_mean+(kx_2*kx_2+ky_2*ky_2)*T_hat_2)=w_hat_2_tot*d_T_0_2_tot  -dz(d_T_hat_2_base)+(w_hat_2_base*dy_T_mean+(kx_2*kx_2+ky_2*ky_2)*T_hat_2_base)')
                        problem.add_equation('dz(S_hat_2)-d_S_hat_2=0')
                        problem.add_equation('-1/tau*dt(S_hat_2)+dz(d_S_hat_2)-1/tau*w_hat_2*dy_S_mean-(kx_2*kx_2+ky_2*ky_2)*S_hat_2=1/tau*(w_hat_2_tot*d_S_0_tot)  -dz(d_S_hat_2_base)+1/tau*w_hat_2_base*dy_S_mean+(kx_2*kx_2+ky_2*ky_2)*S_hat_2_base ')   
                    
                
                else:
                    problem.add_equation('dz(w_hat)-(-(kx*kx+ky*ky)*p_hat)=0')
                    problem.add_equation('dz(p_hat)-(-w_hat+Ra_T*T_hat-Ra_S2T*S_hat)=0')
                    problem.add_equation('dz(T_hat)-d_T_hat=0')
                    problem.add_equation('-dt(T_hat)+dz(d_T_hat)-(w_hat*dy_T_mean+(kx*kx+ky*ky)*T_hat)=w_hat*d_T_0')
                    problem.add_equation('dz(S_hat)-d_S_hat=0')
                    problem.add_equation('-1/tau*dt(S_hat)+dz(d_S_hat)-1/tau*w_hat*dy_S_mean-(kx*kx+ky*ky)*S_hat=1/tau*(w_hat*d_S_0)')   
                    problem.add_equation('dz(T_0)-d_T_0=0')
                    problem.add_equation('dz(S_0)-d_S_0=0')
                
                    if self.kx_2==0 and self.ky_2==0:
                        problem.add_equation('-dt(T_0)+dz(d_T_0)=-2*(kx*kx+ky*ky)*p_hat*T_hat+2*w_hat*d_T_hat')
                        problem.add_equation('-1/tau*dt(S_0)+dz(d_S_0)=1/tau*(-2*(kx*kx+ky*ky)*p_hat*S_hat+2*w_hat*d_S_hat)')
                    else:
                        #Modify the mean equations, also add the contribution from the second harmonic
                        problem.add_equation('-dt(T_0)+dz(d_T_0)=-2*(kx*kx+ky*ky)*p_hat*T_hat+2*w_hat*d_T_hat-2*(kx_2*kx_2+ky_2*ky_2)*p_hat_2*T_hat_2+2*w_hat_2*d_T_hat_2')
                        problem.add_equation('-1/tau*dt(S_0)+dz(d_S_0)=1/tau*(-2*(kx*kx+ky*ky)*p_hat*S_hat+2*w_hat*d_S_hat -2*(kx_2*kx_2+ky_2*ky_2)*p_hat_2*S_hat_2+2*w_hat_2*d_S_hat_2)')
                        #Add the dynamics of the second harmonic
                        problem.add_equation('dz(w_hat_2)-(-(kx_2*kx_2+ky_2*ky_2)*p_hat_2)=0')
                        problem.add_equation('dz(p_hat_2)-(-w_hat_2+Ra_T*T_hat_2-Ra_S2T*S_hat_2)=0')
                        problem.add_equation('dz(T_hat_2)-d_T_hat_2=0')
                        problem.add_equation('-dt(T_hat_2)+dz(d_T_hat_2)-(w_hat_2*dy_T_mean+(kx_2*kx_2+ky_2*ky_2)*T_hat_2)=w_hat_2*d_T_0_2')
                        problem.add_equation('dz(S_hat_2)-d_S_hat_2=0')
                        problem.add_equation('-1/tau*dt(S_hat_2)+dz(d_S_hat_2)-1/tau*w_hat_2*dy_S_mean-(kx_2*kx_2+ky_2*ky_2)*S_hat_2=1/tau*(w_hat_2*d_S_0)')   
                
            # if self.z_bc_w=='periodic' and self.z_bc_T=='periodic' and self.z_bc_S=='periodic':
            #     problem.add_equation('T_0=0',condition="(nz==0)")
            #     problem.add_equation('S_0=0',condition="(nz==0)")
            # else:
            #     problem.add_equation('dz(d_T_0)=-2*(kx*kx+ky*ky)*p_hat*T_hat+2*w_hat*d_T_hat',condition="(nz==0)")
            #     problem.add_equation('dz(d_S_0)=1/tau*(-2*(kx*kx+ky*ky)*p_hat*S_hat+2*w_hat*d_S_hat)',condition="(nz=0)")
            
            #Setup the B.C. update 2021/11/29...
            if self.z_bc_w_left=='dirichlet':
                problem.add_bc("left(w_hat)=0")
                print("Dirichlet B.C. for w left")
            elif self.z_bc_w_left=='neumann':
                problem.add_bc("left(p_hat)=0")
                print("Neumann B.C. for w left")
            #else:
                #raise TypeError('flag.z_bc_w_left is not supported yet') 

            if self.z_bc_w_right=='dirichlet':
                problem.add_bc("right(w_hat)=0")
                print("Dirichlet B.C. for w right")
            elif self.z_bc_w_right=='neumann':
                problem.add_bc("right(p_hat)=0")
                print("Neumann B.C. for w right")
            #else:
                #raise TypeError('flag.z_bc_w_left is not supported yet') 
            if self.z_bc_w_left=='periodic' and self.z_bc_w_right=='periodic':
                problem.add_bc("left(w_hat)-right(w_hat)=0")
                problem.add_bc("left(p_hat)-right(p_hat)=0")
                print("Periodic B.C. for w and p")

            if self.z_bc_T_left=='dirichlet':
                problem.add_bc("left(T_hat)=0")
                problem.add_bc("left(T_0)=0")
                print("Dirichlet B.C. for T left")
            elif self.z_bc_T_left=='neumann':
                problem.add_bc("left(d_T_hat)=0")
                problem.add_bc("left(d_T_0)=0")
                print("Neumann B.C. for T left")
            #else:
            #    raise TypeError('flag.z_bc_T_left is not supported yet') 

            if self.z_bc_T_right=='dirichlet':
                problem.add_bc("right(T_hat)=0")
                problem.add_bc("right(T_0)=0")
                print("Dirichlet B.C. for T right")
            elif self.z_bc_T_right=='neumann':
                problem.add_bc("right(d_T_hat)=0")
                if self.z_bc_T_left=='neumann' and self.z_bc_w_left=='dirichlet' and self.z_bc_w_right=='dirichlet':
                    ##If both sides are Neumann B.C., just set one side for the mean temperature...
                    problem.add_bc("right(T_0)=0")
                else:
                    problem.add_bc("right(d_T_0)=0")
                print("Neumann B.C. for T right")
            #else:
            #    raise TypeError('flag.z_bc_T_right is not supported yet') 

            if self.z_bc_T_left =='periodic' and self.z_bc_T_right=='periodic':
                problem.add_bc("left(T_hat)-right(T_hat)=0")
                problem.add_bc("left(d_T_hat)-right(d_T_hat)=0")
                problem.add_bc("left(T_0)=0")
                problem.add_bc("right(T_0)=0")
                #problem.add_bc("left(d_T_0)-right(d_T_0)=0")
                print("Periodic B.C. for T")
                
            if self.z_bc_S_left=='dirichlet':
                problem.add_bc("left(S_hat)=0")
                problem.add_bc("left(S_0)=0")
                print("Dirichlet B.C. for S left")
            elif self.z_bc_S_left=='neumann':
                problem.add_bc("left(d_S_hat)=0")
                problem.add_bc("left(d_S_0)=0")
                print("Neumann B.C. for S left")
            #else:
            #    raise TypeError('flag.z_bc_S_left is not supported yet') 

            if self.z_bc_S_right=='dirichlet':
                problem.add_bc("right(S_hat)=0")
                problem.add_bc("right(S_0)=0")
                print("Dirichlet B.C. for S right")
            elif self.z_bc_S_right=='neumann':
                problem.add_bc("right(d_S_hat)=0")
                if self.z_bc_S_left=='neumann' and self.z_bc_w_left=='dirichlet' and self.z_bc_w_right=='dirichlet':
                    ##If both sides are Neumann B.C., just set one side for the mean temperature...
                    problem.add_bc("right(S_0)=0")
                else:
                    problem.add_bc("right(d_S_0)=0")
                print("Neumann B.C. for S right")
            #else:
            #    raise TypeError('flag.z_bc_S_right is not supported yet') 

            if self.z_bc_S_left=='periodic' and self.z_bc_S_right=='periodic':
                problem.add_bc("left(S_hat)-right(S_hat)=0")
                problem.add_bc("left(d_S_hat)-right(d_S_hat)=0")
                problem.add_bc("left(S_0)=0")
                problem.add_bc("right(S_0)=0")
                #problem.add_bc("left(d_S_0)-right(d_S_0)=0")
                print("Periodic B.C. for S")
         
        elif self.flow =='HB_porous_2_layer':
            if self.problem =='BVP':
                problem = de.NLBVP(domain, variables=[\
                    'w_hat','p_hat','T_hat','d_T_hat','S_hat','d_S_hat', \
                    'T_0','d_T_0','S_0','d_S_0',\
                    'w_hat_top','p_hat_top','T_hat_top','d_T_hat_top','S_hat_top','d_S_hat_top', \
                    'T_0_top','d_T_0_top','S_0_top','d_S_0_top'])
            problem.parameters['Ra_T'] = self.Ra_T
            problem.parameters['Ra_S2T'] = self.Ra_S2T
            problem.parameters['tau']=self.tau
            problem.parameters['dy_T_mean']=self.dy_T_mean
            problem.parameters['dy_S_mean']=self.dy_S_mean
            problem.parameters['kx']=self.kx
            problem.parameters['ky']=self.ky
            problem.parameters['Omega']=self.HB_porous_2_layer_Omega
            
            #bottom layer equation z\in [0,0.5]
            problem.add_equation('dz(w_hat)-(-(kx*kx+ky*ky)*p_hat)=0')
            problem.add_equation('dz(p_hat)-(-w_hat+Ra_T*T_hat-Ra_S2T*S_hat)=0')
            problem.add_equation('dz(T_hat)-d_T_hat=0')
            problem.add_equation('dz(d_T_hat)-(w_hat*dy_T_mean+(kx*kx+ky*ky)*T_hat)=w_hat*d_T_0')
            problem.add_equation('dz(S_hat)-d_S_hat=0')
            problem.add_equation('dz(d_S_hat)-1/tau*w_hat*dy_S_mean-(kx*kx+ky*ky)*S_hat=1/tau*(w_hat*d_S_0)')   
            problem.add_equation('dz(T_0)-d_T_0=0')
            problem.add_equation('dz(S_0)-d_S_0=0')
            problem.add_equation('dz(d_T_0)=-2*(kx*kx+ky*ky)*p_hat*T_hat+2*w_hat*d_T_hat')
            problem.add_equation('dz(d_S_0)=1/tau*(-2*(kx*kx+ky*ky)*p_hat*S_hat+2*w_hat*d_S_hat)')
             
            #top layer equation z\in [0.5,1]
            problem.add_equation('dz(w_hat_top)-(-(kx*kx+ky*ky)*p_hat_top)=0')
            problem.add_equation('dz(p_hat_top)-(-w_hat_top+Ra_T*T_hat_top-Ra_S2T*S_hat_top)=0')
            problem.add_equation('dz(T_hat_top)-d_T_hat_top=0')
            problem.add_equation('dz(d_T_hat_top)-(w_hat_top*dy_T_mean+(kx*kx+ky*ky)*T_hat_top)=w_hat_top*d_T_0_top')
            problem.add_equation('dz(S_hat_top)-d_S_hat_top=0')
            problem.add_equation('dz(d_S_hat_top)-1/tau*w_hat_top*dy_S_mean-(kx*kx+ky*ky)*S_hat_top=1/tau*(w_hat_top*d_S_0_top)')   
            problem.add_equation('dz(T_0_top)-d_T_0_top=0')
            problem.add_equation('dz(S_0_top)-d_S_0_top=0')
            problem.add_equation('dz(d_T_0_top)=-2*(kx*kx+ky*ky)*p_hat_top*T_hat_top+2*w_hat_top*d_T_hat_top')
            problem.add_equation('dz(d_S_0_top)=1/tau*(-2*(kx*kx+ky*ky)*p_hat_top*S_hat_top+2*w_hat_top*d_S_hat_top)')
               
            #Add B.C. for z=0
            problem.add_bc('left(w_hat)=0')
            problem.add_bc('left(T_hat)=0')
            problem.add_bc('left(S_hat)=0')
            problem.add_bc('left(T_0)=0')
            problem.add_bc('left(S_0)=0')
            
            #Add B.C. at z=1, note that z=1 needs to be implemented at the variable with top
            problem.add_bc('right(w_hat_top)=0')
            problem.add_bc('right(T_hat_top)=0')
            problem.add_bc('right(S_hat_top)=0')
            problem.add_bc('right(T_0_top)=0')
            problem.add_bc('right(S_0_top)=0')
            
            #Add B.C. that they are continuous at the interface.
            problem.add_bc('right(w_hat)-left(w_hat_top)=0')
            problem.add_bc('right(p_hat)-left(p_hat_top)=Omega*right(w_hat)')
            problem.add_bc('right(T_hat)-left(T_hat_top)=0')
            problem.add_bc('right(d_T_hat)-left(d_T_hat_top)=0')
            #Try this B.C. .. This is equivalent that the second order derivative of pressure should be zero... for small Omega, it is also ok...
            #problem.add_bc('right(d_T_hat)-left(d_T_hat_top)=(-kx*kx-ky*ky)/Ra_T*Omega*right(w_hat)')
            
            problem.add_bc('right(S_hat)-left(S_hat_top)=0')
            problem.add_bc('right(d_S_hat)-left(d_S_hat_top)=0')
            problem.add_bc('right(T_0)-left(T_0_top)=0')
            problem.add_bc('right(d_T_0)-left(d_T_0_top)=0')
            problem.add_bc('right(S_0)-left(S_0_top)=0')
            problem.add_bc('right(d_S_0)-left(d_S_0_top)=0')

        elif self.flow=='HB_porous_3_layer':
            problem = de.NLBVP(domain, variables=[\
                    'w_hat','p_hat','T_hat','d_T_hat','S_hat','d_S_hat', \
                    'T_0','d_T_0','S_0','d_S_0',\
                    'w_hat_mid','p_hat_mid','T_hat_mid','d_T_hat_mid','S_hat_mid','d_S_hat_mid', \
                    'T_0_mid','d_T_0_mid','S_0_mid','d_S_0_mid',\
                    'w_hat_top','p_hat_top','T_hat_top','d_T_hat_top','S_hat_top','d_S_hat_top', \
                    'T_0_top','d_T_0_top','S_0_top','d_S_0_top'])
            problem.parameters['Ra_T'] = self.Ra_T
            problem.parameters['Ra_S2T'] = self.Ra_S2T
            problem.parameters['tau']=self.tau
            problem.parameters['dy_T_mean']=self.dy_T_mean
            problem.parameters['dy_S_mean']=self.dy_S_mean
            problem.parameters['kx']=self.kx
            problem.parameters['ky']=self.ky
            problem.parameters['Pi']=self.HB_porous_3_layer_Pi
            problem.parameters['h']=self.HB_porous_3_layer_h
            problem.parameters['z_scale']=2*self.Lz/(1-self.HB_porous_3_layer_h)
            problem.parameters['z_scale_mid']=self.Lz/self.HB_porous_3_layer_h
            problem.parameters['z_scale_top']=2*self.Lz/(1-self.HB_porous_3_layer_h)
            
            #Note that the main difference from the 1 layer and 2 layer model is that
            #we need have additional parameter h and Pi... They need to be multiplied in a couple of places.
            
            #bottom layer equation z\in [0,(1-h)/2]
            problem.add_equation('z_scale*dz(w_hat)-(-(kx*kx+ky*ky)*p_hat)=0')
            problem.add_equation('z_scale*dz(p_hat)-(-w_hat+Ra_T*T_hat-Ra_S2T*S_hat)=0')
            problem.add_equation('z_scale*dz(T_hat)-d_T_hat=0')
            problem.add_equation('z_scale*dz(d_T_hat)-(w_hat*dy_T_mean+(kx*kx+ky*ky)*T_hat)=w_hat*d_T_0')
            problem.add_equation('z_scale*dz(S_hat)-d_S_hat=0')
            problem.add_equation('z_scale*dz(d_S_hat)-1/tau*w_hat*dy_S_mean-(kx*kx+ky*ky)*S_hat=1/tau*(w_hat*d_S_0)')   
            problem.add_equation('z_scale*dz(T_0)-d_T_0=0')
            problem.add_equation('z_scale*dz(S_0)-d_S_0=0')
            problem.add_equation('z_scale*dz(d_T_0)=-2*(kx*kx+ky*ky)*p_hat*T_hat+2*w_hat*d_T_hat')
            problem.add_equation('z_scale*dz(d_S_0)=1/tau*(-2*(kx*kx+ky*ky)*p_hat*S_hat+2*w_hat*d_S_hat)')
             
            #middle layer equation z\in [(1-h)/2,(1+h)/2]
            problem.add_equation('z_scale_mid*dz(w_hat_mid)-(-(kx*kx+ky*ky)*p_hat_mid*Pi)=0')
            problem.add_equation('z_scale_mid*dz(p_hat_mid)-(-w_hat_mid/Pi+Ra_T*T_hat_mid-Ra_S2T*S_hat_mid)=0')
            problem.add_equation('z_scale_mid*dz(T_hat_mid)-d_T_hat_mid=0')
            problem.add_equation('z_scale_mid*dz(d_T_hat_mid)-(w_hat_mid*dy_T_mean+(kx*kx+ky*ky)*T_hat_mid)=w_hat_mid*d_T_0_mid')
            problem.add_equation('z_scale_mid*dz(S_hat_mid)-d_S_hat_mid=0')
            problem.add_equation('z_scale_mid*dz(d_S_hat_mid)-1/tau*w_hat_mid*dy_S_mean-(kx*kx+ky*ky)*S_hat_mid=1/tau*(w_hat_mid*d_S_0_mid)')   
            problem.add_equation('z_scale_mid*dz(T_0_mid)-d_T_0_mid=0')
            problem.add_equation('z_scale_mid*dz(S_0_mid)-d_S_0_mid=0')
            problem.add_equation('z_scale_mid*dz(d_T_0_mid)=-2*(kx*kx+ky*ky)*p_hat_mid*Pi*T_hat_mid+2*w_hat_mid*d_T_hat_mid')
            problem.add_equation('z_scale_mid*dz(d_S_0_mid)=1/tau*(-2*(kx*kx+ky*ky)*p_hat_mid*Pi*S_hat_mid+2*w_hat_mid*d_S_hat_mid)')
            
            #top layer equation z\in [(1+h)/2,1]
            problem.add_equation('z_scale_top*dz(w_hat_top)-(-(kx*kx+ky*ky)*p_hat_top)=0')
            problem.add_equation('z_scale_top*dz(p_hat_top)-(-w_hat_top+Ra_T*T_hat_top-Ra_S2T*S_hat_top)=0')
            problem.add_equation('z_scale_top*dz(T_hat_top)-d_T_hat_top=0')
            problem.add_equation('z_scale_top*dz(d_T_hat_top)-(w_hat_top*dy_T_mean+(kx*kx+ky*ky)*T_hat_top)=w_hat_top*d_T_0_top')
            problem.add_equation('z_scale_top*dz(S_hat_top)-d_S_hat_top=0')
            problem.add_equation('z_scale_top*dz(d_S_hat_top)-1/tau*w_hat_top*dy_S_mean-(kx*kx+ky*ky)*S_hat_top=1/tau*(w_hat_top*d_S_0_top)')   
            problem.add_equation('z_scale_top*dz(T_0_top)-d_T_0_top=0')
            problem.add_equation('z_scale_top*dz(S_0_top)-d_S_0_top=0')
            problem.add_equation('z_scale_top*dz(d_T_0_top)=-2*(kx*kx+ky*ky)*p_hat_top*T_hat_top+2*w_hat_top*d_T_hat_top')
            problem.add_equation('z_scale_top*dz(d_S_0_top)=1/tau*(-2*(kx*kx+ky*ky)*p_hat_top*S_hat_top+2*w_hat_top*d_S_hat_top)')
               
            #Add B.C. for z=0
            problem.add_bc('left(w_hat)=0')
            problem.add_bc('left(T_hat)=0')
            problem.add_bc('left(S_hat)=0')
            problem.add_bc('left(T_0)=0')
            problem.add_bc('left(S_0)=0')
            
            #Add B.C. at z=1, note that z=1 needs to be implemented at the variable with top
            problem.add_bc('right(w_hat_top)=0')
            problem.add_bc('right(T_hat_top)=0')
            problem.add_bc('right(S_hat_top)=0')
            problem.add_bc('right(T_0_top)=0')
            problem.add_bc('right(S_0_top)=0')
            
            #Add B.C. that they are continuous at the interface. This is at interface z=(1-h)/2
            problem.add_bc('right(w_hat)-left(w_hat_mid)=0')
            problem.add_bc('right(p_hat)-left(p_hat_mid)=0')
            problem.add_bc('right(T_hat)-left(T_hat_mid)=0')
            problem.add_bc('right(d_T_hat)-left(d_T_hat_mid)=0')
            problem.add_bc('right(S_hat)-left(S_hat_mid)=0')
            problem.add_bc('right(d_S_hat)-left(d_S_hat_mid)=0')
            problem.add_bc('right(T_0)-left(T_0_mid)=0')
            problem.add_bc('right(d_T_0)-left(d_T_0_mid)=0')
            problem.add_bc('right(S_0)-left(S_0_mid)=0')
            problem.add_bc('right(d_S_0)-left(d_S_0_mid)=0')

            #Add B.C. that they are continuous at the interface. This is at interface z=(1+h)/2
            problem.add_bc('right(w_hat_mid)-left(w_hat_top)=0')
            problem.add_bc('right(p_hat_mid)-left(p_hat_top)=0')
            problem.add_bc('right(T_hat_mid)-left(T_hat_top)=0')
            problem.add_bc('right(d_T_hat_mid)-left(d_T_hat_top)=0')
            problem.add_bc('right(S_hat_mid)-left(S_hat_top)=0')
            problem.add_bc('right(d_S_hat_mid)-left(d_S_hat_top)=0')
            problem.add_bc('right(T_0_mid)-left(T_0_top)=0')
            problem.add_bc('right(d_T_0_mid)-left(d_T_0_top)=0')
            problem.add_bc('right(S_0_mid)-left(S_0_top)=0')
            problem.add_bc('right(d_S_0_mid)-left(d_S_0_top)=0')

            
            
        elif self.flow =='HB_porous_shear':
            #harmonic balance for the porous media
            #For different problem, we need to claim different dedalus problem.
            if self.kx_2==0 and self.kx_2==0:
                #setup the problem, only one scale
                if self.problem =='BVP':
                    problem = de.NLBVP(domain, variables=[\
                        'w_hat_real','p_hat_real','T_hat_real','d_T_hat_real','S_hat_real','d_S_hat_real', \
                        'w_hat_imag','p_hat_imag','T_hat_imag','d_T_hat_imag','S_hat_imag','d_S_hat_imag', \
                            'T_0','d_T_0','S_0','d_S_0','eta','int_phase_cond'])
                    problem.meta['eta']['z']['constant'] = True

                elif self.problem == 'IVP':
                    problem = de.IVP(domain, variables=[\
                        'w_hat','p_hat','T_hat','d_T_hat','S_hat','d_S_hat', \
                        'T_0','d_T_0','S_0','d_S_0'])
            
            problem.parameters['Ra_T'] = self.Ra_T
            problem.parameters['Ra_S2T'] = self.Ra_S2T
            problem.parameters['tau']=self.tau
            problem.parameters['dy_T_mean']=self.dy_T_mean
            problem.parameters['dy_S_mean']=self.dy_S_mean
            problem.parameters['kx']=self.kx
            problem.parameters['ky']=self.ky
            problem.parameters['phi']=self.HB_porous_shear_phi
            #if self.kx_2!=0 or self.kx_2!=0:
                
            if self.problem =='BVP':
                #variable with _hat is the harmonic term
                #variable with _0 is the horizontal average term
                problem.substitutions["U"] = "Ra_T*sin(phi)*(1/2+dy_T_mean*z)-Ra_S2T*sin(phi)*(1/2+dy_S_mean*z)"
                problem.add_equation('dz(w_hat_real)-(-(kx*kx+ky*ky)*p_hat_real)=0')
                problem.add_equation('dz(p_hat_real)-(-w_hat_real+Ra_T*cos(phi)*T_hat_real-Ra_S2T*cos(phi)*S_hat_real)=0')
                problem.add_equation('dz(T_hat_real)-d_T_hat_real=0')
                problem.add_equation('dz(d_T_hat_real)-(w_hat_real*dy_T_mean+(kx*kx+ky*ky)*T_hat_real)+kx*U*T_hat_imag=w_hat_real*d_T_0-kx*eta*T_hat_imag')
                problem.add_equation('dz(S_hat_real)-d_S_hat_real=0')
                problem.add_equation('dz(d_S_hat_real)-1/tau*w_hat_real*dy_S_mean-(kx*kx+ky*ky)*S_hat_real+kx*U*S_hat_imag/tau=1/tau*(w_hat_real*d_S_0)-kx*eta/tau*S_hat_imag')   
                
                problem.add_equation('dz(w_hat_imag)-(-(kx*kx+ky*ky)*p_hat_imag)=0')
                problem.add_equation('dz(p_hat_imag)-(-w_hat_imag+Ra_T*cos(phi)*T_hat_imag-Ra_S2T*cos(phi)*S_hat_imag)=0')
                problem.add_equation('dz(T_hat_imag)-d_T_hat_imag=0')
                problem.add_equation('dz(d_T_hat_imag)-(w_hat_imag*dy_T_mean+(kx*kx+ky*ky)*T_hat_imag)-kx*U*T_hat_real=w_hat_imag*d_T_0+kx*eta*T_hat_real')
                problem.add_equation('dz(S_hat_imag)-d_S_hat_imag=0')
                problem.add_equation('dz(d_S_hat_imag)-1/tau*w_hat_imag*dy_S_mean-(kx*kx+ky*ky)*S_hat_imag-kx*U*S_hat_real/tau=1/tau*(w_hat_imag*d_S_0)+kx*eta/tau*S_hat_real')   
                
                
                problem.add_equation('dz(T_0)-d_T_0=0')
                problem.add_equation('dz(S_0)-d_S_0=0')
                problem.add_equation('dz(d_T_0)=-2*(kx*kx+ky*ky)*(p_hat_real*T_hat_real+p_hat_imag*T_hat_imag)+2*(w_hat_real*d_T_hat_real+w_hat_imag*d_T_hat_imag)')
                problem.add_equation('dz(d_S_0)=1/tau*(-2*(kx*kx+ky*ky)*(p_hat_real*S_hat_real+p_hat_imag*S_hat_imag)+2*(w_hat_real*d_S_hat_real+w_hat_imag*d_S_hat_imag) )')

                #problem.add_bc("left(p_hat_imag)=0")

                #integration of phase condition                
                problem.add_equation('dz(int_phase_cond)-T_hat_imag=0')
                problem.add_bc("left(int_phase_cond)=0")
                problem.add_bc("right(int_phase_cond)=0")
                
            #Setup the B.C. update 2021/11/29...
            if self.z_bc_w_left=='dirichlet':
                problem.add_bc("left(w_hat_real)=0")
                problem.add_bc("left(w_hat_imag)=0")
                print("Dirichlet B.C. for w left")
            elif self.z_bc_w_left=='neumann':
                problem.add_bc("left(p_hat_real)=0")
                problem.add_bc("left(p_hat_imag)=0")
                print("Neumann B.C. for w left")
            #else:
                #raise TypeError('flag.z_bc_w_left is not supported yet') 

            if self.z_bc_w_right=='dirichlet':
                problem.add_bc("right(w_hat_real)=0")
                problem.add_bc("right(w_hat_imag)=0")
                print("Dirichlet B.C. for w right")
            elif self.z_bc_w_right=='neumann':
                problem.add_bc("right(p_hat_real)=0")
                problem.add_bc("right(p_hat_imag)=0")
                print("Neumann B.C. for w right")
            #else:
                #raise TypeError('flag.z_bc_w_left is not supported yet') 
            if self.z_bc_w_left=='periodic' and self.z_bc_w_right=='periodic':
                problem.add_bc("left(w_hat_real)-right(w_hat_real)=0")
                problem.add_bc("left(p_hat_real)-right(p_hat_real)=0")
                problem.add_bc("left(w_hat_imag)-right(w_hat_imag)=0")
                problem.add_bc("left(p_hat_imag)-right(p_hat_imag)=0")
                print("Periodic B.C. for w and p")

            if self.z_bc_T_left=='dirichlet':
                problem.add_bc("left(T_hat_real)=0")
                problem.add_bc("left(T_hat_imag)=0")
                problem.add_bc("left(T_0)=0")
                
                print("Dirichlet B.C. for T left")
            elif self.z_bc_T_left=='neumann':
                problem.add_bc("left(d_T_hat_real)=0")
                problem.add_bc("left(d_T_hat_imag)=0")
                problem.add_bc("left(d_T_0)=0")
                print("Neumann B.C. for T left")
            #else:
            #    raise TypeError('flag.z_bc_T_left is not supported yet') 

            if self.z_bc_T_right=='dirichlet':
                problem.add_bc("right(T_hat_real)=0")
                problem.add_bc("right(T_hat_imag)=0")
                problem.add_bc("right(T_0)=0")
                print("Dirichlet B.C. for T right")
            elif self.z_bc_T_right=='neumann':
                problem.add_bc("right(d_T_hat_real)=0")
                problem.add_bc("right(d_T_hat_imag)=0")
                if self.z_bc_T_left=='neumann' and self.z_bc_w_left=='dirichlet' and self.z_bc_w_right=='dirichlet':
                    ##If both sides are Neumann B.C., just set one side for the mean temperature...
                    problem.add_bc("right(T_0)=0")
                else:
                    problem.add_bc("right(d_T_0)=0")
                print("Neumann B.C. for T right")
            #else:
            #    raise TypeError('flag.z_bc_T_right is not supported yet') 

            if self.z_bc_T_left =='periodic' and self.z_bc_T_right=='periodic':
                problem.add_bc("left(T_hat_real)-right(T_hat_real)=0")
                problem.add_bc("left(d_T_hat_real)-right(d_T_hat_real)=0")
                problem.add_bc("left(T_hat_imag)-right(T_hat_imag)=0")
                problem.add_bc("left(d_T_hat_imag)-right(d_T_hat_imag)=0")
                problem.add_bc("left(T_0)=0")
                problem.add_bc("right(T_0)=0")
                #problem.add_bc("left(d_T_0)-right(d_T_0)=0")
                print("Periodic B.C. for T")
                
            if self.z_bc_S_left=='dirichlet':
                problem.add_bc("left(S_hat_real)=0")
                problem.add_bc("left(S_hat_imag)=0")
                problem.add_bc("left(S_0)=0")
                print("Dirichlet B.C. for S left")
            elif self.z_bc_S_left=='neumann':
                problem.add_bc("left(d_S_hat_real)=0")
                problem.add_bc("left(d_S_hat_imag)=0")
                problem.add_bc("left(d_S_0)=0")
                print("Neumann B.C. for S left")
            #else:
            #    raise TypeError('flag.z_bc_S_left is not supported yet') 

            if self.z_bc_S_right=='dirichlet':
                problem.add_bc("right(S_hat_real)=0")
                problem.add_bc("right(S_hat_imag)=0")
                problem.add_bc("right(S_0)=0")
                print("Dirichlet B.C. for S right")
            elif self.z_bc_S_right=='neumann':
                problem.add_bc("right(d_S_hat_real)=0")
                problem.add_bc("right(d_S_hat_imag)=0")
                if self.z_bc_S_left=='neumann' and self.z_bc_w_left=='dirichlet' and self.z_bc_w_right=='dirichlet':
                    ##If both sides are Neumann B.C., just set one side for the mean temperature...
                    problem.add_bc("right(S_0)=0")
                else:
                    problem.add_bc("right(d_S_0)=0")
                print("Neumann B.C. for S right")
            #else:
            #    raise TypeError('flag.z_bc_S_right is not supported yet') 

            if self.z_bc_S_left=='periodic' and self.z_bc_S_right=='periodic':
                problem.add_bc("left(S_hat_real)-right(S_hat_real)=0")
                problem.add_bc("left(d_S_hat_real)-right(d_S_hat_real)=0")
                problem.add_bc("left(S_hat_imag)-right(S_hat_imag)=0")
                problem.add_bc("left(d_S_hat_imag)-right(d_S_hat_imag)=0")
            
                problem.add_bc("left(S_0)=0")
                problem.add_bc("right(S_0)=0")
                #problem.add_bc("left(d_S_0)-right(d_S_0)=0")
                print("Periodic B.C. for S")
        
            # if self.z_bc_w_left=='dirichlet' and self.z_bc_w_right=='dirichlet':
            #     problem.add_bc("left(p_hat_imag)=0")
            # else:
            #     problem.add_bc("left(w_hat_imag)=0")
        
            
            #The bottom block is the old version that try to implement the second harmonic and the periodic B.C.... but both of them seems not working very well.
            #Here, try to add the second harmonic and the periodic B.C., 2021/11/28
            """
            #setup different B.C. 
            #for w is always no penetration
            if self.z_bc_w == 'dirichlet':
                problem.add_bc("left(w_hat) = 0")
                problem.add_bc("right(w_hat) = 0")
                #Also add BC for the second harmonic
                if not (self.kx_2==0 and self.ky_2==0):
                    problem.add_bc("left(w_hat_2) = 0")
                    problem.add_bc("right(w_hat_2) = 0")
                print('Dirichlet B.C. for w')
            elif self.z_bc_w=='periodic':
                problem.add_bc("left(w_hat)-right(w_hat)=0")
                problem.add_bc("left(p_hat)-right(p_hat)=0")
                if not (self.kx_2==0 and self.ky_2==0):
                    problem.add_bc("left(w_hat_2)-right(w_hat_2) = 0")
                    problem.add_bc("left(p_hat_2)-right(p_hat_2) = 0")
                
                print('Periodic B.C. for w')
            else:
                raise TypeError('flag.z_bc_w is not supported yet') 


            #set up the B.C. for temperature
            if self.z_bc_T =='dirichlet':
                problem.add_bc("left(T_0) = 0")
                problem.add_bc("right(T_0) = 0")
                problem.add_bc("left(T_hat) = 0")
                problem.add_bc("right(T_hat) = 0")
                if not (self.kx_2==0 and self.ky_2==0):
                     problem.add_bc("left(T_hat_2) = 0")
                     problem.add_bc("right(T_hat_2) = 0")
                
                print('Dirichlet B.C. for T')

            elif self.z_bc_T =='neumann':
                
                #Neumann B.C. for horizontal average temperature, 
                #it only needs to be set up in one side and the other side will automatically satisfy the Neumann B.C. 
                problem.add_bc("left(d_T_0) = 0")
                problem.add_bc("right(T_0) = 0")
                problem.add_bc("left(d_T_hat) = 0")
                problem.add_bc("right(d_T_hat) = 0")
                if not (self.kx_2==0 and self.ky_2==0):
                    problem.add_bc("left(d_T_hat_2) = 0")
                    problem.add_bc("right(d_T_hat_2) = 0")
                
                print('Neumann B.C. for T')

            elif self.z_bc_T =='periodic':
                #not fully reddy!!!!
                #problem.add_bc("left(T_0)=0")
                problem.add_bc("left(T_0)-right(T_0)=0")
                problem.add_bc("left(d_T_0)-right(d_T_0)=0")
                problem.add_bc("left(T_hat)-right(T_hat)=0")
                problem.add_bc("left(d_T_hat)-right(d_T_hat)=0")
                if not (self.kx_2==0 and self.ky_2==0):
                    problem.add_bc("left(T_hat_2)-right(T_hat_2)=0")
                    problem.add_bc("left(d_T_hat_2)-right(d_T_hat_2)=0")
                
                #problem.add_bc("left(T_0)=0")
                print('Periodic B.C. for T')
            else:
                raise TypeError('flag.z_bc_T is not supported yet') 

            ##set up different B.C. for salinity 
            if self.z_bc_S=='dirichlet':
                problem.add_bc("left(S_0) = 0")
                problem.add_bc("right(S_0) = 0")
                problem.add_bc("left(S_hat) = 0")
                problem.add_bc("right(S_hat) = 0")
                if not (self.kx_2==0 and self.ky_2==0):
                    problem.add_bc("left(S_hat_2) = 0")
                    problem.add_bc("right(S_hat_2) = 0")
               
                print('Dirichlet B.C. for S')

            elif self.z_bc_S=='neumann':
                
                #Neumann B.C. for horizontal average temperature, 
                #it only needs to be set up in one side and the other side will automatically satisfy the Neumann B.C. 
                #Instead, we set up a gauge by Dirichlet B.C. 
                problem.add_bc("left(d_S_0) = 0")
                problem.add_bc("right(S_0) = 0") 
                problem.add_bc("left(d_S_hat) = 0")
                problem.add_bc("right(d_S_hat) = 0")
                #second harmonic
                if not (self.kx_2==0 and self.ky_2==0):
                    problem.add_bc("left(d_S_hat_2) = 0")
                    problem.add_bc("right(d_S_hat_2) = 0")
                
                print('Neumann B.C. for S')

            elif self.z_bc_S=='periodic':
                #not fully ready!!!
                print('Periodic B.C. for S')
                #problem.add_bc("left(S_0)=0")
                problem.add_bc("left(S_0)-right(S_0)=0")
                problem.add_bc("left(d_S_0)-right(d_S_0)=0")
                problem.add_bc("left(S_hat)-right(S_hat)=0")
                problem.add_bc("left(d_S_hat)-right(d_S_hat)=0")
                #second harmonic
                if not (self.kx_2==0 and self.ky_2==0):
                    problem.add_bc("left(S_hat_2)-right(S_hat_2)=0")
                    problem.add_bc("left(d_S_hat_2)-right(d_S_hat_2)=0")
                
            else:
                raise TypeError('flag.z_bc_S is not supported yet') 
            
            #elif self.z_bc_T_S_w == 'periodic':
                #need to do nothing for periodic BC but change the basis as Fourier at the beginning
            """
             
        elif self.flow =='HB_benard':
            #Harmonic balance for Benard problem at high Prandtl number
            #Not fully benchmarked!!!!
            
            if self.problem =='BVP':
                problem = de.NLBVP(domain, variables=\
                    ['u_tilde','d_u_tilde','v_tilde','d_v_tilde', \
                    'w_hat','p_hat','T_hat','d_T_hat', \
                    'S_hat','d_S_hat','T_0','d_T_0','S_0','d_S_0'])
            elif self.problem == 'IVP':
                problem = de.IVP(domain, variables=\
                    ['u_tilde','d_u_tilde','v_tilde','d_v_tilde', \
                    'w_hat','p_hat','T_hat','d_T_hat', \
                    'S_hat','d_S_hat','T_0','d_T_0','S_0','d_S_0'])
            elif self.problem =='EVP':
                problem = de.EVP(domain, variables=\
                    ['u_tilde','d_u_tilde','v_tilde','d_v_tilde', \
                    'w_hat','p_hat','T_hat','d_T_hat', \
                    'S_hat','d_S_hat','T_0','d_T_0','S_0','d_S_0'] \
                        ,eigenvalue='eig_val',tolerance=self.EVP_homogeneous_tolerance)
            
            problem.parameters['Pr'] = self.Pr 
            problem.parameters['Ra_T'] = self.Ra_T
            problem.parameters['Ra_S2T'] = self.Ra_S2T
            problem.parameters['tau']=self.tau
            problem.parameters['dy_T_mean']=self.dy_T_mean
            problem.parameters['dy_S_mean']=self.dy_S_mean
            problem.parameters['kx']=self.kx
            problem.parameters['ky']=self.ky
            problem.parameters['Pe_T']=self.Pe_T
            problem.parameters['Pe_S']=self.Pe_S
            problem.parameters['F_sin']=self.F_sin
            problem.parameters['ks']=self.ks
            problem.parameters['j']=1j
            if self.problem =='BVP':
                problem.add_equation('dz(u_tilde)-d_u_tilde=0')
                problem.add_equation('dz(d_u_tilde)-(kx*p_hat+(kx*kx+ky*ky)*u_tilde)=0')
                problem.add_equation('dz(v_tilde)-d_v_tilde=0')
                problem.add_equation('dz(d_v_tilde)-(ky*p_hat+(kx*kx+ky*ky)*v_tilde)=0')
                problem.add_equation('dz(w_hat)-(kx*u_tilde+ky*v_tilde)=0')
                problem.add_equation('dz(p_hat)-(kx*d_u_tilde+ky*d_v_tilde-(kx*kx+ky*ky)*w_hat+Ra_T*T_hat-Ra_S2T*S_hat)=0')
                problem.add_equation('dz(T_hat)-d_T_hat=0')
                problem.add_equation('dz(S_hat)-d_S_hat=0')
                problem.add_equation('dz(T_0)-d_T_0=0')
                problem.add_equation('dz(S_0)-d_S_0=0')

                if self.F_sin==0:
                    problem.add_equation('dz(d_T_hat)-w_hat*dy_T_mean-(kx*kx+ky*ky)*T_hat=Pe_T*w_hat*d_T_0')
                    problem.add_equation('dz(d_T_0)=Pe_T*(2*kx*u_tilde*T_hat+2*ky*v_tilde*T_hat+2*w_hat*d_T_hat)')
                    problem.add_equation('dz(d_S_hat)-1/tau*w_hat*dy_S_mean-(kx*kx+ky*ky)*S_hat=Pe_S/tau*(w_hat*d_S_0)')   
                    problem.add_equation('dz(d_S_0)=Pe_S/tau*(2*kx*u_tilde*S_hat+2*ky*v_tilde*S_hat+2*w_hat*d_S_hat)')
                else:
                    problem.add_equation('dz(d_T_hat)-w_hat*dy_T_mean-(kx*kx+ky*ky)*T_hat-Pe_T*j*kx*F_sin*sin(ks*z)*T_hat=Pe_T*w_hat*d_T_0')
                    #problem.add_equation('dz(d_T_0)=Pe_T*(kx*conj(u_tilde)*T_hat+kx*u_tilde*conj(T_hat)+ky*conj(v_tilde)*T_hat+ky*v_tilde*conj(T_hat)+conj(w_hat)*d_T_hat+w_hat*conj(d_T_hat))')
                    problem.add_equation('dz(d_T_0)=Pe_T*(kx*(u_tilde)*T_hat+kx*u_tilde*(T_hat)+ky*(v_tilde)*T_hat+ky*v_tilde*(T_hat)+(w_hat)*d_T_hat+w_hat*(d_T_hat))')

                    problem.add_equation('dz(d_S_hat)-1/tau*w_hat*dy_S_mean-(kx*kx+ky*ky)*S_hat-Pe_S/tau*j*kx*F_sin*sin(ks*z)*S_hat=Pe_S/tau*(w_hat*d_S_0)')   
                    #problem.add_equation('dz(d_S_0)=Pe_S/tau*(kx*conj(u_tilde)*S_hat+kx*u_tilde*conj(S_hat)+ky*conj(v_tilde)*S_hat+ky*v_tilde*conj(S_hat)+conj(w_hat)*d_S_hat+w_hat*conj(d_S_hat))')
                    problem.add_equation('dz(d_S_0)=Pe_S/tau*(kx*(u_tilde)*S_hat+kx*u_tilde*(S_hat)+ky*(v_tilde)*S_hat+ky*v_tilde*(S_hat)+(w_hat)*d_S_hat+w_hat*(d_S_hat))')

            elif self.problem=='IVP':
                problem.add_equation('dz(u_tilde)-d_u_tilde=0')
                problem.add_equation('-1/Pr*dt(u_tilde)+dz(d_u_tilde)-(kx*p_hat+(kx*kx+ky*ky)*u_tilde)=0')
                problem.add_equation('dz(v_tilde)-d_v_tilde=0')
                problem.add_equation('-1/Pr*dt(v_tilde)+dz(d_v_tilde)-(ky*p_hat+(kx*kx+ky*ky)*v_tilde)=0')
                problem.add_equation('dz(w_hat)-(kx*u_tilde+ky*v_tilde)=0')
                problem.add_equation('1/Pr*dt(w_hat)+dz(p_hat)-(kx*d_u_tilde+ky*d_v_tilde-(kx*kx+ky*ky)*w_hat+Ra_T*T_hat-Ra_S2T*S_hat)=0')
                problem.add_equation('dz(T_hat)-d_T_hat=0')
                problem.add_equation('-dt(T_hat)+dz(d_T_hat)-(w_hat*dy_T_mean+(kx*kx+ky*ky)*T_hat)=w_hat*d_T_0')
                problem.add_equation('dz(S_hat)-d_S_hat=0')
                problem.add_equation('-1/tau*dt(S_hat)+dz(d_S_hat)-1/tau*w_hat*dy_S_mean-(kx*kx+ky*ky)*S_hat=1/tau*(w_hat*d_S_0)')   
                problem.add_equation('dz(T_0)-d_T_0=0')
                problem.add_equation('-dt(T_0)+dz(d_T_0)=2*kx*u_tilde*T_hat+2*ky*v_tilde*T_hat+2*w_hat*d_T_hat')
                problem.add_equation('dz(S_0)-d_S_0=0')
                problem.add_equation('-1/tau*dt(S_0)+dz(d_S_0)=1/tau*(2*kx*u_tilde*S_hat+2*ky*v_tilde*S_hat+2*w_hat*d_S_hat)')
            
            elif self.problem =='EVP':
                problem.substitutions['dt(A)'] = "eig_val*A"
                if self.EVP_secondary:
                
                    state=solver_in.state
                    for varname in ['u_tilde','d_u_tilde','v_tilde','d_v_tilde', \
                                    'w_hat','p_hat','T_hat','d_T_hat', \
                                        'S_hat','d_S_hat','T_0','d_T_0','S_0','d_S_0']:
                        problem.substitutions['{0}_tot'.format(varname)]='{0}_base'.format(varname)+'+'+varname 
                        ncc = domain.new_field(name='{0}_base'.format(varname))
                        ncc['c'] = state[varname]['c']
                        problem.parameters['{0}_base'.format(varname)] = ncc
                               
                    #Make sure that the RHS is zero when variable in EVP is evaluated as zero.    
                    #put total variable on the RHS, but need to substract the linear part of the base state. This will guarantee that when perturbation=0 (without _tot or _base) , RHS is zero (homogeneous). Dedalus will check this
                    #The LHS is not tricky, just the original problem will be OK.
                    
                    problem.add_equation('dz(u_tilde)-d_u_tilde=0')
                    problem.add_equation('-1/Pr*dt(u_tilde)+dz(d_u_tilde)-(kx*p_hat+(kx*kx+ky*ky)*u_tilde)=0')
                    problem.add_equation('dz(v_tilde)-d_v_tilde=0')
                    problem.add_equation('-1/Pr*dt(v_tilde)+dz(d_v_tilde)-(ky*p_hat+(kx*kx+ky*ky)*v_tilde)=0')
                    problem.add_equation('dz(w_hat)-(kx*u_tilde+ky*v_tilde)=0')
                    problem.add_equation('1/Pr*dt(w_hat)+dz(p_hat)-(kx*d_u_tilde+ky*d_v_tilde-(kx*kx+ky*ky)*w_hat+Ra_T*T_hat-Ra_S2T*S_hat)=0')
                    problem.add_equation('dz(T_hat)-d_T_hat=0')
                    problem.add_equation('-dt(T_hat)+dz(d_T_hat)-(w_hat*dy_T_mean+(kx*kx+ky*ky)*T_hat)=w_hat_tot*d_T_0_tot -dz(d_T_hat_base)+(w_hat_base*dy_T_mean+(kx*kx+ky*ky)*T_hat_base) ')
                    problem.add_equation('dz(S_hat)-d_S_hat=0')
                    problem.add_equation('-1/tau*dt(S_hat)+dz(d_S_hat)-1/tau*w_hat*dy_S_mean-(kx*kx+ky*ky)*S_hat=1/tau*(w_hat_tot*d_S_0_tot)  -dz(d_S_hat_base)+1/tau*w_hat_base*dy_S_mean+(kx*kx+ky*ky)*S_hat_base ')   
                    problem.add_equation('dz(T_0)-d_T_0=0')
                    problem.add_equation('-dt(T_0)+dz(d_T_0)=2*kx*u_tilde_tot*T_hat_tot+2*ky*v_tilde_tot*T_hat_tot+2*w_hat_tot*d_T_hat_tot  -dz(d_T_0_base)')
                    problem.add_equation('dz(S_0)-d_S_0=0')
                    problem.add_equation('-1/tau*dt(S_0)+dz(d_S_0)=1/tau*(2*kx*u_tilde_tot*S_hat_tot+2*ky*v_tilde_tot*S_hat_tot+2*w_hat_tot*d_S_hat_tot) -dz(d_S_0_base)')
                
                else:
                    #I do not have any data to load, just solve the eigenvalue problem
                    problem.add_equation('dz(u_tilde)-d_u_tilde=0')
                    problem.add_equation('-1/Pr*dt(u_tilde)+dz(d_u_tilde)-(kx*p_hat+(kx*kx+ky*ky)*u_tilde)=0')
                    problem.add_equation('dz(v_tilde)-d_v_tilde=0')
                    problem.add_equation('-1/Pr*dt(v_tilde)+dz(d_v_tilde)-(ky*p_hat+(kx*kx+ky*ky)*v_tilde)=0')
                    problem.add_equation('dz(w_hat)-(kx*u_tilde+ky*v_tilde)=0')
                    problem.add_equation('1/Pr*dt(w_hat)+dz(p_hat)-(kx*d_u_tilde+ky*d_v_tilde-(kx*kx+ky*ky)*w_hat+Ra_T*T_hat-Ra_S2T*S_hat)=0')
                    problem.add_equation('dz(T_hat)-d_T_hat=0')
                    problem.add_equation('-dt(T_hat)+dz(d_T_hat)-(w_hat*dy_T_mean+(kx*kx+ky*ky)*T_hat)=w_hat*d_T_0')
                    problem.add_equation('dz(S_hat)-d_S_hat=0')
                    problem.add_equation('-1/tau*dt(S_hat)+dz(d_S_hat)-1/tau*w_hat*dy_S_mean-(kx*kx+ky*ky)*S_hat=1/tau*(w_hat*d_S_0)')   
                    problem.add_equation('dz(T_0)-d_T_0=0')
                    problem.add_equation('-dt(T_0)+dz(d_T_0)=2*kx*u_tilde*T_hat+2*ky*v_tilde*T_hat+2*w_hat*d_T_hat')
                    problem.add_equation('dz(S_0)-d_S_0=0')
                    problem.add_equation('-1/tau*dt(S_0)+dz(d_S_0)=1/tau*(2*kx*u_tilde*S_hat+2*ky*v_tilde*S_hat+2*w_hat*d_S_hat)')
                
                    
            if self.z_bc_w_left=='dirichlet':
                problem.add_bc("left(w_hat)=0")
                print("Dirichlet for w left")
            if self.z_bc_w_right=='dirichlet':
                problem.add_bc("right(w_hat)=0")
                print("Dirichlet for w right")
            
            if self.z_bc_w_left=='periodic' and self.z_bc_w_right=='periodic':
                problem.add_bc("left(w_hat)-right(w_hat)=0")
                problem.add_bc("left(p_hat)-right(p_hat)=0")
                
            if self.z_bc_T_left=='dirichlet':
                problem.add_bc("left(T_hat)=0")
                problem.add_bc("left(T_0)=0")
                print("Dirichlet for T left")
            elif self.z_bc_T_left=='neumann':
                problem.add_bc("left(d_T_hat)=0")
                problem.add_bc("left(d_T_0)=0")
                print("Neumann for T left")
                
            if self.z_bc_T_right=='dirichlet':
                problem.add_bc("right(T_hat)=0")
                problem.add_bc("right(T_0)=0")
                print("Dirichlet for T right")
            elif self.z_bc_T_right=='neumann':
                problem.add_bc("right(d_T_hat)=0")
                problem.add_bc("right(d_T_0)=0")
                print("Neumann for T right")
            
            if self.z_bc_T_left=='periodic' and self.z_bc_T_right=='periodic':
                problem.add_bc("left(T_hat)-right(T_hat)=0")
                problem.add_bc("left(d_T_hat)-right(d_T_hat)=0")
                problem.add_bc("left(T_0)=0")
                problem.add_bc("right(T_0)=0")
                #problem.add_bc("left(d_T_0)-right(d_T_0)=0")
                print("Periodic for T")
               
            if self.z_bc_S_left=='dirichlet':
                problem.add_bc("left(S_hat)=0")
                problem.add_bc("left(S_0)=0")
                print("Dirichlet for S left")
            elif self.z_bc_S_left=='neumann':
                problem.add_bc("left(d_S_hat)=0")
                problem.add_bc("left(d_S_0)=0")
                print("Neumann for S left")
                
            if self.z_bc_S_right=='dirichlet':
                problem.add_bc("right(S_hat)=0")
                problem.add_bc("right(S_0)=0")
                print("Dirichlet for S right")
            elif self.z_bc_S_right=='neumann':
                problem.add_bc("right(d_S_hat)=0")
                problem.add_bc("right(d_S_0)=0")
                print("Neumann for S right")
            
            if self.z_bc_S_left=='periodic' and self.z_bc_S_right=='periodic':
                problem.add_bc("left(S_hat)-right(S_hat)=0")
                problem.add_bc("left(d_S_hat)-right(d_S_hat)=0")
                problem.add_bc("left(S_0)=0")
                problem.add_bc("right(S_0)=0")
                #problem.add_bc("left(d_S_0)-right(d_S_0)=0")
                print("Periodic for S")
           
            if self.z_bc_u_v_left=='dirichlet':
                problem.add_bc("left(u_tilde)=0")
                problem.add_bc("left(v_tilde)=0")
                print("Dirichet for u,v left")
            elif self.z_bc_u_v_left=='neumann':
                problem.add_bc("left(d_u_tilde)=0")
                problem.add_bc("left(d_v_tilde)=0")
                print('Neumann for u,v left')
                
            if self.z_bc_u_v_right=='dirichlet':
                problem.add_bc("right(u_tilde)=0")
                problem.add_bc("right(v_tilde)=0")
                print('Dirichlet for u,v, right')
            elif self.z_bc_u_v_right=='neumann':
                problem.add_bc("right(d_u_tilde)=0")
                problem.add_bc("right(d_v_tilde)=0")
                print('Neumann for u,v right')
                
            if self.z_bc_u_v_left=='periodic' and self.z_bc_u_v_right=='periodic':
                problem.add_bc("left(u_tilde)-right(u_tilde)=0")
                problem.add_bc("left(v_tilde)-right(v_tilde)=0")
                problem.add_bc("left(d_u_tilde)-right(d_u_tilde)=0")
                problem.add_bc("left(d_v_tilde)-right(d_v_tilde)=0")
                print('Periodic for u,v')
            # if self.z_bc_u_v =='dirichlet':
            #     problem.add_bc("left(u_tilde) = 0")
            #     problem.add_bc("right(u_tilde) = 0")
            #     problem.add_bc("left(v_tilde) = 0")
            #     problem.add_bc("right(v_tilde) = 0")
                
            # elif self.z_bc_u_v =='neumann':
            #     problem.add_bc("left(d_u_tilde) = 0")
            #     problem.add_bc("right(d_u_tilde) = 0")
            #     problem.add_bc("left(d_v_tilde) = 0")
            #     problem.add_bc("right(d_v_tilde) = 0")
             
             
            # if self.z_bc_w == 'dirichlet':
            #     problem.add_bc("left(w_hat) = 0")
            #     problem.add_bc("right(w_hat) = 0")
            # elif self.z_bc_W=='periodic':
            #     print('Periodic B.C. for w')
            # else:
            #     raise TypeError('flag.z_bc_w is not supported yet') 
    
    
            # if self.z_bc_S=='dirichlet':
            #     problem.add_bc("left(S_hat) = 0")
            #     problem.add_bc("right(S_hat) = 0")
            #     problem.add_bc("left(S_0) = 0")
            #     problem.add_bc("right(S_0) = 0")
            # elif self.z_bc_S=='neumann':
            #     problem.add_bc("left(d_S_hat) = 0")
            #     problem.add_bc("right(d_S_hat) = 0")
            #     problem.add_bc("left(d_S_0) = 0")
            #     problem.add_bc("right(S_0) = 0")  
            # elif self.z_bc_S=='periodic':
            #     print('Periodic B.C. for S')
            # else:
            #     raise TypeError('flag.z_bc_S is not supported yet') 

        elif self.flow =='HB_benard_2_layer':
            #Harmonic balance for Benard problem at high Prandtl number
            #Not fully benchmarked!!!!
            
            if self.problem =='BVP':
                problem = de.NLBVP(domain, variables=\
                    ['u_tilde','d_u_tilde','v_tilde','d_v_tilde', \
                    'w_hat','p_hat','T_hat','d_T_hat', \
                    'S_hat','d_S_hat','T_0','d_T_0','S_0','d_S_0', \
                    'u_tilde_top','d_u_tilde_top','v_tilde_top','d_v_tilde_top',\
                    'w_hat_top','p_hat_top','T_hat_top','d_T_hat_top',\
                    'S_hat_top','d_S_hat_top','T_0_top','d_T_0_top','S_0_top','d_S_0_top'])
            
            elif self.problem == 'IVP':
                problem = de.IVP(domain, variables=\
                    ['u_tilde','d_u_tilde','v_tilde','d_v_tilde', \
                    'w_hat','p_hat','T_hat','d_T_hat', \
                    'S_hat','d_S_hat','T_0','d_T_0','S_0','d_S_0',\
                    'u_tilde_top','d_u_tilde_top','v_tilde_top','d_v_tilde_top', \
                    'w_hat_top','p_hat_top','T_hat_top','d_T_hat_top', \
                    'S_hat_top','d_S_hat_top','T_0_top','d_T_0_top','S_0_top','d_S_0_top'])
                    
            problem.parameters['Pr'] = self.Pr 
            problem.parameters['Ra_T'] = self.Ra_T
            problem.parameters['Ra_S2T'] = self.Ra_S2T
            problem.parameters['tau']=self.tau
            problem.parameters['dy_T_mean']=self.dy_T_mean
            problem.parameters['dy_S_mean']=self.dy_S_mean
            problem.parameters['kx']=self.kx
            problem.parameters['ky']=self.ky
            problem.parameters['Pe_T']=self.Pe_T
            problem.parameters['Pe_S']=self.Pe_S
            problem.parameters['F_sin']=self.F_sin
            problem.parameters['ks']=self.ks
            problem.parameters['j']=1j
            if self.problem =='BVP':
                problem.add_equation('dz(u_tilde)-d_u_tilde=0')
                problem.add_equation('dz(d_u_tilde)-(kx*p_hat+(kx*kx+ky*ky)*u_tilde)=0')
                problem.add_equation('dz(v_tilde)-d_v_tilde=0')
                problem.add_equation('dz(d_v_tilde)-(ky*p_hat+(kx*kx+ky*ky)*v_tilde)=0')
                problem.add_equation('dz(w_hat)-(kx*u_tilde+ky*v_tilde)=0')
                problem.add_equation('dz(p_hat)-(kx*d_u_tilde+ky*d_v_tilde-(kx*kx+ky*ky)*w_hat+Ra_T*T_hat-Ra_S2T*S_hat)=0')
                problem.add_equation('dz(T_hat)-d_T_hat=0')
                problem.add_equation('dz(S_hat)-d_S_hat=0')
                problem.add_equation('dz(T_0)-d_T_0=0')
                problem.add_equation('dz(S_0)-d_S_0=0')

                problem.add_equation('dz(u_tilde_top)-d_u_tilde_top=0')
                problem.add_equation('dz(d_u_tilde_top)-(kx*p_hat_top+(kx*kx+ky*ky)*u_tilde_top)=0')
                problem.add_equation('dz(v_tilde_top)-d_v_tilde_top=0')
                problem.add_equation('dz(d_v_tilde_top)-(ky*p_hat_top+(kx*kx+ky*ky)*v_tilde_top)=0')
                problem.add_equation('dz(w_hat_top)-(kx*u_tilde_top+ky*v_tilde_top)=0')
                problem.add_equation('dz(p_hat_top)-(kx*d_u_tilde_top+ky*d_v_tilde_top-(kx*kx+ky*ky)*w_hat_top+Ra_T*T_hat_top-Ra_S2T*S_hat_top)=0')
                problem.add_equation('dz(T_hat_top)-d_T_hat_top=0')
                problem.add_equation('dz(S_hat_top)-d_S_hat_top=0')
                problem.add_equation('dz(T_0_top)-d_T_0_top=0')
                problem.add_equation('dz(S_0_top)-d_S_0_top=0')


                if self.F_sin==0:
                    problem.add_equation('dz(d_T_hat)-w_hat*dy_T_mean-(kx*kx+ky*ky)*T_hat=Pe_T*w_hat*d_T_0')
                    problem.add_equation('dz(d_T_0)=Pe_T*(2*kx*u_tilde*T_hat+2*ky*v_tilde*T_hat+2*w_hat*d_T_hat)')
                    problem.add_equation('dz(d_S_hat)-1/tau*w_hat*dy_S_mean-(kx*kx+ky*ky)*S_hat=Pe_S/tau*(w_hat*d_S_0)')   
                    problem.add_equation('dz(d_S_0)=Pe_S/tau*(2*kx*u_tilde*S_hat+2*ky*v_tilde*S_hat+2*w_hat*d_S_hat)')
                    
                    problem.add_equation('dz(d_T_hat_top)-w_hat_top*dy_T_mean-(kx*kx+ky*ky)*T_hat_top=Pe_T*w_hat_top*d_T_0_top')
                    problem.add_equation('dz(d_T_0_top)=Pe_T*(2*kx*u_tilde_top*T_hat_top+2*ky*v_tilde_top*T_hat_top+2*w_hat_top*d_T_hat_top)')
                    problem.add_equation('dz(d_S_hat_top)-1/tau*w_hat_top*dy_S_mean-(kx*kx+ky*ky)*S_hat_top=Pe_S/tau*(w_hat_top*d_S_0_top)')   
                    problem.add_equation('dz(d_S_0_top)=Pe_S/tau*(2*kx*u_tilde_top*S_hat_top+2*ky*v_tilde_top*S_hat_top+2*w_hat_top*d_S_hat_top)')
                
                else:
                    problem.add_equation('dz(d_T_hat)-w_hat*dy_T_mean-(kx*kx+ky*ky)*T_hat-Pe_T*j*kx*F_sin*sin(ks*z)*T_hat=Pe_T*w_hat*d_T_0')
                    #problem.add_equation('dz(d_T_0)=Pe_T*(kx*conj(u_tilde)*T_hat+kx*u_tilde*conj(T_hat)+ky*conj(v_tilde)*T_hat+ky*v_tilde*conj(T_hat)+conj(w_hat)*d_T_hat+w_hat*conj(d_T_hat))')
                    problem.add_equation('dz(d_T_0)=Pe_T*(kx*(u_tilde)*T_hat+kx*u_tilde*(T_hat)+ky*(v_tilde)*T_hat+ky*v_tilde*(T_hat)+(w_hat)*d_T_hat+w_hat*(d_T_hat))')

                    problem.add_equation('dz(d_S_hat)-1/tau*w_hat*dy_S_mean-(kx*kx+ky*ky)*S_hat-Pe_S/tau*j*kx*F_sin*sin(ks*z)*S_hat=Pe_S/tau*(w_hat*d_S_0)')   
                    #problem.add_equation('dz(d_S_0)=Pe_S/tau*(kx*conj(u_tilde)*S_hat+kx*u_tilde*conj(S_hat)+ky*conj(v_tilde)*S_hat+ky*v_tilde*conj(S_hat)+conj(w_hat)*d_S_hat+w_hat*conj(d_S_hat))')
                    problem.add_equation('dz(d_S_0)=Pe_S/tau*(kx*(u_tilde)*S_hat+kx*u_tilde*(S_hat)+ky*(v_tilde)*S_hat+ky*v_tilde*(S_hat)+(w_hat)*d_S_hat+w_hat*(d_S_hat))')

                    problem.add_equation('dz(d_T_hat_top)-w_hat_top*dy_T_mean-(kx*kx+ky*ky)*T_hat_top-Pe_T*j*kx*F_sin*sin(ks*z)*T_hat_top=Pe_T*w_hat_top*d_T_0_top')
                    #problem.add_equation('dz(d_T_0)=Pe_T*(kx*conj(u_tilde)*T_hat+kx*u_tilde*conj(T_hat)+ky*conj(v_tilde)*T_hat+ky*v_tilde*conj(T_hat)+conj(w_hat)*d_T_hat+w_hat*conj(d_T_hat))')
                    problem.add_equation('dz(d_T_0_top)=Pe_T*(kx*(u_tilde_top)*T_hat_top+kx*u_tilde_top*(T_hat_top)+ky*(v_tilde_top)*T_hat_top+ky*v_tilde_top*(T_hat_top)+(w_hat_top)*d_T_hat_top+w_hat_top*(d_T_hat_top))')

                    problem.add_equation('dz(d_S_hat_top)-1/tau*w_hat_top*dy_S_mean-(kx*kx+ky*ky)*S_hat_top-Pe_S/tau*j*kx*F_sin*sin(ks*z)*S_hat_top=Pe_S/tau*(w_hat_top*d_S_0_top)')   
                    #problem.add_equation('dz(d_S_0)=Pe_S/tau*(kx*conj(u_tilde)*S_hat+kx*u_tilde*conj(S_hat)+ky*conj(v_tilde)*S_hat+ky*v_tilde*conj(S_hat)+conj(w_hat)*d_S_hat+w_hat*conj(d_S_hat))')
                    problem.add_equation('dz(d_S_0_top)=Pe_S/tau*(kx*(u_tilde_top)*S_hat_top+kx*u_tilde_top*(S_hat_top)+ky*(v_tilde_top)*S_hat_top+ky*v_tilde_top*(S_hat_top)+(w_hat_top)*d_S_hat_top+w_hat_top*(d_S_hat_top))')


            elif self.problem=='IVP':
                problem.add_equation('dz(u_tilde)-d_u_tilde=0')
                problem.add_equation('-1/Pr*dt(u_tilde)+dz(d_u_tilde)-(kx*p_hat+(kx*kx+ky*ky)*u_tilde)=0')
                problem.add_equation('dz(v_tilde)-d_v_tilde=0')
                problem.add_equation('-1/Pr*dt(v_tilde)+dz(d_v_tilde)-(ky*p_hat+(kx*kx+ky*ky)*v_tilde)=0')
                problem.add_equation('dz(w_hat)-(kx*u_tilde+ky*v_tilde)=0')
                problem.add_equation('1/Pr*dt(w_hat)+dz(p_hat)-(kx*d_u_tilde+ky*d_v_tilde-(kx*kx+ky*ky)*w_hat+Ra_T*T_hat-Ra_S2T*S_hat)=0')
                problem.add_equation('dz(T_hat)-d_T_hat=0')
                problem.add_equation('-dt(T_hat)+dz(d_T_hat)-(w_hat*dy_T_mean+(kx*kx+ky*ky)*T_hat)=w_hat*d_T_0')
                problem.add_equation('dz(S_hat)-d_S_hat=0')
                problem.add_equation('-1/tau*dt(S_hat)+dz(d_S_hat)-1/tau*w_hat*dy_S_mean-(kx*kx+ky*ky)*S_hat=1/tau*(w_hat*d_S_0)')   
                problem.add_equation('dz(T_0)-d_T_0=0')
                problem.add_equation('-dt(T_0)+dz(d_T_0)=2*kx*u_tilde*T_hat+2*ky*v_tilde*T_hat+2*w_hat*d_T_hat')
                problem.add_equation('dz(S_0)-d_S_0=0')
                problem.add_equation('-1/tau*dt(S_0)+dz(d_S_0)=1/tau*(2*kx*u_tilde*S_hat+2*ky*v_tilde*S_hat+2*w_hat*d_S_hat)')
            
                problem.add_equation('dz(u_tilde_top)-d_u_tilde_top=0')
                problem.add_equation('-1/Pr*dt(u_tilde_top)+dz(d_u_tilde_top)-(kx*p_hat_top+(kx*kx+ky*ky)*u_tilde_top)=0')
                problem.add_equation('dz(v_tilde_top)-d_v_tilde_top=0')
                problem.add_equation('-1/Pr*dt(v_tilde_top)+dz(d_v_tilde_top)-(ky*p_hat_top+(kx*kx+ky*ky)*v_tilde_top)=0')
                problem.add_equation('dz(w_hat_top)-(kx*u_tilde_top+ky*v_tilde_top)=0')
                problem.add_equation('1/Pr*dt(w_hat_top)+dz(p_hat_top)-(kx*d_u_tilde_top+ky*d_v_tilde_top-(kx*kx+ky*ky)*w_hat_top+Ra_T*T_hat_top-Ra_S2T*S_hat_top)=0')
                problem.add_equation('dz(T_hat_top)-d_T_hat_top=0')
                problem.add_equation('-dt(T_hat_top)+dz(d_T_hat_top)-(w_hat_top*dy_T_mean+(kx*kx+ky*ky)*T_hat_top)=w_hat_top*d_T_0_top')
                problem.add_equation('dz(S_hat_top)-d_S_hat_top=0')
                problem.add_equation('-1/tau*dt(S_hat_top)+dz(d_S_hat_top)-1/tau*w_hat_top*dy_S_mean-(kx*kx+ky*ky)*S_hat_top=1/tau*(w_hat_top*d_S_0_top)')   
                problem.add_equation('dz(T_0_top)-d_T_0_top=0')
                problem.add_equation('-dt(T_0_top)+dz(d_T_0_top)=2*kx*u_tilde_top*T_hat_top+2*ky*v_tilde_top*T_hat_top+2*w_hat_top*d_T_hat_top')
                problem.add_equation('dz(S_0_top)-d_S_0_top=0')
                problem.add_equation('-1/tau*dt(S_0_top)+dz(d_S_0_top)=1/tau*(2*kx*u_tilde_top*S_hat_top+2*ky*v_tilde_top*S_hat_top+2*w_hat_top*d_S_hat_top)')
                        
            if self.z_bc_w_left=='dirichlet':
                problem.add_bc("left(w_hat)=0")
                print("Dirichlet for w left")
            if self.z_bc_w_right=='dirichlet':
                problem.add_bc("right(w_hat_top)=0")
                print("Dirichlet for w right")
            
            if self.z_bc_w_left=='periodic' and self.z_bc_w_right=='periodic':
                problem.add_bc("left(w_hat)-right(w_hat_top)=0")
                problem.add_bc("left(p_hat)-right(p_hat_top)=0")
                
            if self.z_bc_T_left=='dirichlet':
                problem.add_bc("left(T_hat)=0")
                problem.add_bc("left(T_0)=0")
                print("Dirichlet for T left")
            elif self.z_bc_T_left=='neumann':
                problem.add_bc("left(d_T_hat)=0")
                problem.add_bc("left(d_T_0)=0")
                print("Neumann for T left")
                
            if self.z_bc_T_right=='dirichlet':
                problem.add_bc("right(T_hat_top)=0")
                problem.add_bc("right(T_0_top)=0")
                print("Dirichlet for T right")
            elif self.z_bc_T_right=='neumann':
                problem.add_bc("right(d_T_hat_top)=0")
                problem.add_bc("right(d_T_0_top)=0")
                print("Neumann for T right")
            
            if self.z_bc_T_left=='periodic' and self.z_bc_T_right=='periodic':
                problem.add_bc("left(T_hat)-right(T_hat_top)=0")
                problem.add_bc("left(d_T_hat)-right(d_T_hat_top)=0")
                problem.add_bc("left(T_0)=0")
                problem.add_bc("right(T_0_top)=0")
                #problem.add_bc("left(d_T_0)-right(d_T_0)=0")
                print("Periodic for T")
               
            if self.z_bc_S_left=='dirichlet':
                problem.add_bc("left(S_hat)=0")
                problem.add_bc("left(S_0)=0")
                print("Dirichlet for S left")
            elif self.z_bc_S_left=='neumann':
                problem.add_bc("left(d_S_hat)=0")
                problem.add_bc("left(d_S_0)=0")
                print("Neumann for S left")
                
            if self.z_bc_S_right=='dirichlet':
                problem.add_bc("right(S_hat_top)=0")
                problem.add_bc("right(S_0_top)=0")
                print("Dirichlet for S right")
            elif self.z_bc_S_right=='neumann':
                problem.add_bc("right(d_S_hat_top)=0")
                problem.add_bc("right(d_S_0_top)=0")
                print("Neumann for S right")
            
            if self.z_bc_S_left=='periodic' and self.z_bc_S_right=='periodic':
                problem.add_bc("left(S_hat)-right(S_hat_top)=0")
                problem.add_bc("left(d_S_hat)-right(d_S_hat_top)=0")
                problem.add_bc("left(S_0)=0")
                problem.add_bc("right(S_0_top)=0")
                #problem.add_bc("left(d_S_0)-right(d_S_0)=0")
                print("Periodic for S")
           
            if self.z_bc_u_v_left=='dirichlet':
                problem.add_bc("left(u_tilde)=0")
                problem.add_bc("left(v_tilde)=0")
            elif self.z_bc_u_v_left=='neumann':
                problem.add_bc("left(d_u_tilde)=0")
                problem.add_bc("left(d_v_tilde)=0")
                
            if self.z_bc_u_v_right=='dirichlet':
                problem.add_bc("right(u_tilde_top)=0")
                problem.add_bc("right(v_tilde_top)=0")
            elif self.z_bc_u_v_right=='neumann':
                problem.add_bc("right(d_u_tilde_top)=0")
                problem.add_bc("right(d_v_tilde_top)=0")
                
            if self.z_bc_u_v_left=='periodic' and self.z_bc_u_v_right=='periodic':
                problem.add_bc("left(u_tilde)-right(u_tilde_top)=0")
                problem.add_bc("left(v_tilde)-right(v_tilde_top)=0")
                problem.add_bc("left(d_u_tilde)-right(d_u_tilde_top)=0")
                problem.add_bc("left(d_v_tilde)-right(d_v_tilde_top)=0")
            
            #Add B.C. that they are continuous at the interface.
            problem.add_bc('right(u_tilde)-left(u_tilde_top)=0')
            problem.add_bc('right(v_tilde)-left(v_tilde_top)=0')
            problem.add_bc('right(d_u_tilde)-left(d_u_tilde_top)=0')
            problem.add_bc('right(d_v_tilde)-left(d_v_tilde_top)=0')
            problem.add_bc('right(w_hat)-left(w_hat_top)=0')
            problem.add_bc('right(p_hat)-left(p_hat_top)=0')
            problem.add_bc('right(T_hat)-left(T_hat_top)=0')
            problem.add_bc('right(d_T_hat)-left(d_T_hat_top)=0')
            
            problem.add_bc('right(S_hat)-left(S_hat_top)=0')
            problem.add_bc('right(d_S_hat)-left(d_S_hat_top)=0')
            problem.add_bc('right(T_0)-left(T_0_top)=0')
            problem.add_bc('right(d_T_0)-left(d_T_0_top)=0')
            problem.add_bc('right(S_0)-left(S_0_top)=0')
            problem.add_bc('right(d_S_0)-left(d_S_0_top)=0')
    
            



        elif self.flow =='HB_benard_shear':
            """
            if self.problem =='BVP':
                problem = de.NLBVP(domain, variables=\
                    ['u_tilde_real','d_u_tilde_real','v_tilde_real','d_v_tilde_real', \
                     'w_hat_real','p_hat_real','T_hat_real','d_T_hat_real', \
                     'T_hat_imag','d_T_hat_imag', \
                    'S_hat_real','d_S_hat_real','S_hat_imag','d_S_hat_imag', \
                        'T_0','d_T_0','S_0','d_S_0'])
                    #'u_tilde_imag','d_u_tilde_imag','v_tilde_imag','d_v_tilde_imag', \
                    #'w_hat_imag','p_hat_imag',
            elif self.problem == 'IVP':
                problem = de.IVP(domain, variables=\
                    ['u_tilde','d_u_tilde','v_tilde','d_v_tilde', \
                    'w_hat','p_hat','T_hat','d_T_hat', \
                    'S_hat','d_S_hat','T_0','d_T_0','S_0','d_S_0'])
                    
            problem.parameters['Pr'] = self.Pr 
            problem.parameters['Ra_T'] = self.Ra_T
            problem.parameters['Ra_S2T'] = self.Ra_S2T
            problem.parameters['tau']=self.tau
            problem.parameters['dy_T_mean']=self.dy_T_mean
            problem.parameters['dy_S_mean']=self.dy_S_mean
            problem.parameters['kx']=self.kx
            problem.parameters['ky']=self.ky
            problem.parameters['Pe_T']=self.Pe_T
            problem.parameters['Pe_S']=self.Pe_S
            problem.parameters['ks']=self.ks
            problem.parameters['j']=1j
            if self.F_sin=='z':
                print('Couette shear')
            else:
                problem.parameters['F_sin']=self.F_sin
            
            if self.problem =='BVP':
                #real
                problem.add_equation('dz(u_tilde_real)-d_u_tilde_real=0')
                problem.add_equation('dz(d_u_tilde_real)-(kx*p_hat_real+(kx*kx+ky*ky)*u_tilde_real)=0')
                problem.add_equation('dz(v_tilde_real)-d_v_tilde_real=0')
                problem.add_equation('dz(d_v_tilde_real)-(ky*p_hat_real+(kx*kx+ky*ky)*v_tilde_real)=0')
                problem.add_equation('dz(w_hat_real)-(kx*u_tilde_real+ky*v_tilde_real)=0')
                problem.add_equation('dz(p_hat_real)-(kx*d_u_tilde_real+ky*d_v_tilde_real-(kx*kx+ky*ky)*w_hat_real+Ra_T*T_hat_real-Ra_S2T*S_hat_real)=0')
                problem.add_equation('dz(T_hat_real)-d_T_hat_real=0')
                problem.add_equation('dz(S_hat_real)-d_S_hat_real=0')
                
                #imag
                #problem.add_equation('dz(u_tilde_imag)-d_u_tilde_imag=0')
                #problem.add_equation('dz(d_u_tilde_imag)-(kx*p_hat_imag+(kx*kx+ky*ky)*u_tilde_imag)=0')
                #problem.add_equation('dz(v_tilde_imag)-d_v_tilde_imag=0')
                #problem.add_equation('dz(d_v_tilde_imag)-(ky*p_hat_imag+(kx*kx+ky*ky)*v_tilde_imag)=0')
                #problem.add_equation('dz(w_hat_imag)-(kx*u_tilde_imag+ky*v_tilde_imag)=0')
                #problem.add_equation('dz(p_hat_imag)-(kx*d_u_tilde_imag+ky*d_v_tilde_imag-(kx*kx+ky*ky)*w_hat_imag+Ra_T*T_hat_imag-Ra_S2T*S_hat_imag)=0')
                #problem.add_equation('Ra_T*T_hat_imag-Ra_S2T*S_hat_imag=0')
                problem.add_equation('dz(T_hat_imag)-d_T_hat_imag=0')
                problem.add_equation('dz(S_hat_imag)-d_S_hat_imag=0')
                #problem.add_equation('dz(T_0_imag)-d_T_0_imag=0')
                #problem.add_equation('dz(S_0_imag)-d_S_0_imag=0')

                #mean temperature and salinity
                problem.add_equation('dz(T_0)-d_T_0=0')
                problem.add_equation('dz(S_0)-d_S_0=0')
                #problem.add_equation('dz(d_T_0)=Pe_T*(2*kx*u_tilde_real*T_hat_real+2*kx*u_tilde_imag*T_hat_imag+2*ky*v_tilde_real*T_hat_real+2*ky*v_tilde_imag*T_hat_imag+2*w_hat_real*d_T_hat_real+2*w_hat_imag*d_T_hat_imag)')
                #problem.add_equation('dz(d_S_0)=Pe_S/tau*(2*kx*u_tilde_real*S_hat_real+2*kx*u_tilde_imag*S_hat_imag+2*ky*v_tilde_real*S_hat_real+2*ky*v_tilde_imag*S_hat_imag+2*w_hat_real*d_S_hat_real+2*w_hat_imag*d_S_hat_imag)')
                problem.add_equation('dz(d_T_0)=Pe_T*(2*kx*u_tilde_real*T_hat_real+2*ky*v_tilde_real*T_hat_real+2*w_hat_real*d_T_hat_real)')
                problem.add_equation('dz(d_S_0)=Pe_S/tau*(2*kx*u_tilde_real*S_hat_real+2*ky*v_tilde_real*S_hat_real+2*w_hat_real*d_S_hat_real)')

                #coupling between real and imag due to shear
                if self.F_sin=='z':
                    #problem.add_equation('dz(d_T_hat_real)-w_hat_real*dy_T_mean-(kx*kx+ky*ky)*T_hat_real+Pe_T*kx*(z-1/2)*T_hat_imag=Pe_T*w_hat_real*d_T_0')
                    #problem.add_equation('dz(d_T_hat_imag)-w_hat_imag*dy_T_mean-(kx*kx+ky*ky)*T_hat_imag-Pe_T*kx*(z-1/2)*T_hat_real=Pe_T*w_hat_imag*d_T_0')
                    #problem.add_equation('dz(d_S_hat_real)-1/tau*w_hat_real*dy_S_mean-(kx*kx+ky*ky)*S_hat_real+Pe_S/tau*kx*(z-1/2)*S_hat_imag=Pe_S/tau*(w_hat_real*d_S_0)')   
                    #problem.add_equation('dz(d_S_hat_imag)-1/tau*w_hat_imag*dy_S_mean-(kx*kx+ky*ky)*S_hat_imag-Pe_S/tau*kx*(z-1/2)*S_hat_real=Pe_S/tau*(w_hat_imag*d_S_0)')   
                    problem.add_equation('dz(d_T_hat_real)-w_hat_real*dy_T_mean-(kx*kx+ky*ky)*T_hat_real+Pe_T*kx*(z-1/2)*T_hat_imag=Pe_T*w_hat_real*d_T_0')
                    problem.add_equation('dz(d_T_hat_imag)-(kx*kx+ky*ky)*T_hat_imag-Pe_T*kx*(z-1/2)*T_hat_real=0')
                    problem.add_equation('dz(d_S_hat_real)-1/tau*w_hat_real*dy_S_mean-(kx*kx+ky*ky)*S_hat_real+Pe_S/tau*kx*(z-1/2)*S_hat_imag=Pe_S/tau*(w_hat_real*d_S_0)')   
                    problem.add_equation('dz(d_S_hat_imag)-(kx*kx+ky*ky)*S_hat_imag-Pe_S/tau*kx*(z-1/2)*S_hat_real=0')   
                
                
                else:
                    print(self.F_sin)
                    problem.add_equation('dz(d_T_hat_real)-w_hat_real*dy_T_mean-(kx*kx+ky*ky)*T_hat_real+Pe_T*kx*F_sin*sin(ks*z)*T_hat_imag=Pe_T*w_hat_real*d_T_0')
                    problem.add_equation('dz(d_T_hat_imag)-(kx*kx+ky*ky)*T_hat_imag-Pe_T*kx*F_sin*sin(ks*z)*T_hat_real=0')
                    problem.add_equation('dz(d_S_hat_real)-1/tau*w_hat_real*dy_S_mean-(kx*kx+ky*ky)*S_hat_real+Pe_S/tau*kx*F_sin*sin(ks*z)*S_hat_imag=Pe_S/tau*(w_hat_real*d_S_0)')   
                    problem.add_equation('dz(d_S_hat_imag)-(kx*kx+ky*ky)*S_hat_imag-Pe_S/tau*kx*F_sin*sin(ks*z)*S_hat_real=0')   
                #else:
                #    print("Wrong flag of F_sin.")

            elif self.problem=='IVP':
                problem.add_equation('dz(u_tilde)-d_u_tilde=0')
                problem.add_equation('-1/Pr*dt(u_tilde)+dz(d_u_tilde)-(kx*p_hat+(kx*kx+ky*ky)*u_tilde)=0')
                problem.add_equation('dz(v_tilde)-d_v_tilde=0')
                problem.add_equation('-1/Pr*dt(v_tilde)+dz(d_v_tilde)-(ky*p_hat+(kx*kx+ky*ky)*v_tilde)=0')
                problem.add_equation('dz(w_hat)-(kx*u_tilde+ky*v_tilde)=0')
                problem.add_equation('1/Pr*dt(w_hat)+dz(p_hat)-(kx*d_u_tilde+ky*d_v_tilde-(kx*kx+ky*ky)*w_hat+Ra_T*T_hat-Ra_S2T*S_hat)=0')
                problem.add_equation('dz(T_hat)-d_T_hat=0')
                problem.add_equation('-dt(T_hat)+dz(d_T_hat)-(w_hat*dy_T_mean+(kx*kx+ky*ky)*T_hat)=w_hat*d_T_0')
                problem.add_equation('dz(S_hat)-d_S_hat=0')
                problem.add_equation('-1/tau*dt(S_hat)+dz(d_S_hat)-1/tau*w_hat*dy_S_mean-(kx*kx+ky*ky)*S_hat=1/tau*(w_hat*d_S_0)')   
                problem.add_equation('dz(T_0)-d_T_0=0')
                problem.add_equation('-dt(T_0)+dz(d_T_0)=2*kx*u_tilde*T_hat+2*ky*v_tilde*T_hat+2*w_hat*d_T_hat')
                problem.add_equation('dz(S_0)-d_S_0=0')
                problem.add_equation('-1/tau*dt(S_0)+dz(d_S_0)=1/tau*(2*kx*u_tilde*S_hat+2*ky*v_tilde*S_hat+2*w_hat*d_S_hat)')
            
            if self.z_bc_w_left=='periodic' and self.z_bc_w_right=='periodic':
                problem.add_bc("left(w_hat_real)-right(w_hat_real)=0")
                problem.add_bc("left(p_hat_real)-right(p_hat_real)=0")
                #problem.add_bc("left(w_hat_imag)-right(w_hat_imag)=0")
                #problem.add_bc("left(p_hat_imag)-right(p_hat_imag)=0")
           
            if self.z_bc_w_left=='dirichlet':
                problem.add_bc("left(w_hat_real)=0")
                #problem.add_bc("left(w_hat_imag)=0")
                print("Dirichlet for w left")
                
            if self.z_bc_w_right=='dirichlet':
                problem.add_bc("right(w_hat_real)=0")
                #problem.add_bc("right(w_hat_imag)=0")
                print("Dirichlet for w right")
           
            if self.z_bc_T_left=='periodic' and self.z_bc_T_right=='periodic':
                problem.add_bc("left(T_hat_real)-right(T_hat_real)=0")
                problem.add_bc("left(d_T_hat_real)-right(d_T_hat_real)=0")
                problem.add_bc("left(T_hat_imag)-right(T_hat_imag)=0")
                problem.add_bc("left(d_T_hat_imag)-right(d_T_hat_imag)=0")
                problem.add_bc("left(T_0)=0")
                problem.add_bc("right(T_0)=0")
                #problem.add_bc("left(d_T_0)-right(d_T_0)=0")
                print("Periodic for T")
              
            if self.z_bc_T_left=='dirichlet':
                problem.add_bc("left(T_hat_real)=0")
                problem.add_bc("left(T_hat_imag)=0")
                problem.add_bc("left(T_0)=0")
                print("Dirichlet for T left")
            elif self.z_bc_T_left=='neumann':
                problem.add_bc("left(d_T_hat_real)=0")
                problem.add_bc("left(d_T_hat_imag)=0")
                problem.add_bc("left(d_T_0)=0")
                print("Neumann for T left")
                
            if self.z_bc_T_right=='dirichlet':
                problem.add_bc("right(T_hat_real)=0")
                problem.add_bc("right(T_hat_imag)=0")
                problem.add_bc("right(T_0)=0")
                print("Dirichlet for T right")
            elif self.z_bc_T_right=='neumann':
                problem.add_bc("right(d_T_hat_real)=0")
                problem.add_bc("right(d_T_hat_imag)=0")
                problem.add_bc("right(d_T_0)=0")
                print("Neumann for T right")
                
            if self.z_bc_S_left=='periodic' and self.z_bc_S_right=='periodic':
                problem.add_bc("left(S_hat_real)-right(S_hat_real)=0")
                problem.add_bc("left(d_S_hat_real)-right(d_S_hat_real)=0")
                problem.add_bc("left(S_hat_imag)-right(S_hat_imag)=0")
                problem.add_bc("left(d_S_hat_imag)-right(d_S_hat_imag)=0")
                problem.add_bc("left(S_0)=0")
                problem.add_bc("right(S_0)=0")
                #problem.add_bc("left(d_S_0)-right(d_S_0)=0")
                print("Periodic for S")
           
            if self.z_bc_S_left=='dirichlet':
                problem.add_bc("left(S_hat_real)=0")
                problem.add_bc("left(S_hat_imag)=0")
                problem.add_bc("left(S_0)=0")
                print("Dirichlet for S left")
            elif self.z_bc_S_left=='neumann':
                problem.add_bc("left(d_S_hat_real)=0")
                problem.add_bc("left(d_S_hat_imag)=0")
                problem.add_bc("left(d_S_0)=0")
                print("Neumann for S left")
                
            if self.z_bc_S_right=='dirichlet':
                problem.add_bc("right(S_hat_real)=0")
                problem.add_bc("right(S_hat_imag)=0")
                problem.add_bc("right(S_0)=0")
                print("Dirichlet for S right")
            elif self.z_bc_S_right=='neumann':
                problem.add_bc("right(d_S_hat_real)=0")
                problem.add_bc("right(d_S_hat_imag)=0")
                problem.add_bc("right(d_S_0)=0")
                print("Neumann for S right")
            
            if self.z_bc_u_v_left=='periodic' and self.z_bc_u_v_right=='periodic':
                problem.add_bc("left(u_tilde_real)-right(u_tilde_real)=0")
                problem.add_bc("left(v_tilde_real)-right(v_tilde_real)=0")
                problem.add_bc("left(d_u_tilde_real)-right(d_u_tilde_real)=0")
                problem.add_bc("left(d_v_tilde_real)-right(d_v_tilde_real)=0")
                #problem.add_bc("left(u_tilde_imag)-right(u_tilde_imag)=0")
                #problem.add_bc("left(v_tilde_imag)-right(v_tilde_imag)=0")
                #problem.add_bc("left(d_u_tilde_imag)-right(d_u_tilde_imag)=0")
                #problem.add_bc("left(d_v_tilde_imag)-right(d_v_tilde_imag)=0")
                print("Periodic for u,v")

            if self.z_bc_u_v_left=='dirichlet':
                problem.add_bc("left(u_tilde_real)=0")
                #problem.add_bc("left(u_tilde_imag)=0")
                problem.add_bc("left(v_tilde_real)=0")
                #problem.add_bc("left(v_tilde_imag)=0")
                print("Dirichlet for u,v left")

            elif self.z_bc_u_v_left=='neumann':
                problem.add_bc("left(d_u_tilde_real)=0")
                #problem.add_bc("left(d_u_tilde_imag)=0")
                problem.add_bc("left(d_v_tilde_real)=0")
                #problem.add_bc("left(d_v_tilde_imag)=0")
                print("Neumann for u,v left")

            if self.z_bc_u_v_right=='dirichlet':
                problem.add_bc("right(u_tilde_real)=0")
                #problem.add_bc("right(u_tilde_imag)=0")
                problem.add_bc("right(v_tilde_real)=0")
                #problem.add_bc("right(v_tilde_imag)=0")
                print("Dirichlet for u,v right")

            elif self.z_bc_u_v_right=='neumann':
                problem.add_bc("right(d_u_tilde_real)=0")
                #problem.add_bc("right(d_u_tilde_imag)=0")
                problem.add_bc("right(d_v_tilde_real)=0")
                #problem.add_bc("right(d_v_tilde_imag)=0")
                print("Neumann for u,v right")
            """
            
            
            ##This version have one issue that cannot fix the phase, for any solution, multiply by e^{i\theta} is again a solution 
            ##Thus, the Newton iteration cannot converge.
            ##
            if self.problem =='BVP':
                problem = de.NLBVP(domain, variables=\
                    ['u_tilde_real','d_u_tilde_real','v_tilde_real','d_v_tilde_real', \
                    'u_tilde_imag','d_u_tilde_imag','v_tilde_imag','d_v_tilde_imag', \
                     'w_hat_real','p_hat_real','T_hat_real','d_T_hat_real', \
                     'w_hat_imag','p_hat_imag','T_hat_imag','d_T_hat_imag', \
                    'S_hat_real','d_S_hat_real','S_hat_imag','d_S_hat_imag', \
                        'T_0','d_T_0','S_0','d_S_0','eta',\
                            'U_0','d_U_0']) # add large scale shear, 2022/05/04
                problem.meta['eta']['z']['constant'] = True

            elif self.problem == 'IVP':
                problem = de.IVP(domain, variables=\
                    ['u_tilde_real','d_u_tilde_real','v_tilde_real','d_v_tilde_real', \
                    'u_tilde_imag','d_u_tilde_imag','v_tilde_imag','d_v_tilde_imag', \
                     'w_hat_real','p_hat_real','T_hat_real','d_T_hat_real', \
                     'w_hat_imag','p_hat_imag','T_hat_imag','d_T_hat_imag', \
                    'S_hat_real','d_S_hat_real','S_hat_imag','d_S_hat_imag', \
                        'T_0','d_T_0','S_0','d_S_0',\
                            'U_0','d_U_0']) # add large scale shear, 2022/05/04
            elif self.problem == 'EVP': 
                problem = de.EVP(domain, variables=\
                    ['u_tilde','d_u_tilde','v_tilde','d_v_tilde', \
                    'w_hat','p_hat','T_hat','d_T_hat', \
                    'S_hat','d_S_hat'],eigenvalue='eig_val',tolerance=self.EVP_homogeneous_tolerance)
            
                
            problem.parameters['Pr'] = self.Pr 
            problem.parameters['Ra_T'] = self.Ra_T
            problem.parameters['Ra_S2T'] = self.Ra_S2T
            problem.parameters['tau']=self.tau
            problem.parameters['dy_T_mean']=self.dy_T_mean
            problem.parameters['dy_S_mean']=self.dy_S_mean
            problem.parameters['kx']=self.kx
            problem.parameters['ky']=self.ky
            problem.parameters['Pe_T']=self.Pe_T
            problem.parameters['Pe_S']=self.Pe_S
            problem.parameters['ks']=self.ks
            #problem.parameters['j']=1j
            #problem.parameters['z0']=self.IBM_z0
            #problem.parameters['sigma']=self.IBM_sigma
            #problem.parameters['A']=self.IBM_A
            if self.F_sin=='z':
                print('Couette shear')
            else:
                problem.parameters['F_sin']=self.F_sin
            
            if self.problem =='BVP':
                #real
                problem.add_equation('dz(u_tilde_real)-d_u_tilde_real=0')
                problem.add_equation('dz(d_u_tilde_real)-(kx*p_hat_real+(kx*kx+ky*ky)*u_tilde_real)=0')
                problem.add_equation('dz(v_tilde_real)-d_v_tilde_real=0')
                problem.add_equation('dz(d_v_tilde_real)-(ky*p_hat_real+(kx*kx+ky*ky)*v_tilde_real)=0')
                problem.add_equation('dz(w_hat_real)-(kx*u_tilde_real+ky*v_tilde_real)=0')
                problem.add_equation('dz(p_hat_real)-(kx*d_u_tilde_real+ky*d_v_tilde_real-(kx*kx+ky*ky)*w_hat_real+Ra_T*T_hat_real-Ra_S2T*S_hat_real)=0')
                problem.add_equation('dz(T_hat_real)-d_T_hat_real=0')
                problem.add_equation('dz(S_hat_real)-d_S_hat_real=0')
                
                #IBM: +A*exp(-(z-z0)**2/sigma**2)*w_hat_imag
                #imag
                problem.add_equation('dz(u_tilde_imag)-d_u_tilde_imag=0')
                problem.add_equation('dz(d_u_tilde_imag)-(kx*p_hat_imag+(kx*kx+ky*ky)*u_tilde_imag)=0')
                problem.add_equation('dz(v_tilde_imag)-d_v_tilde_imag=0')
                problem.add_equation('dz(d_v_tilde_imag)-(ky*p_hat_imag+(kx*kx+ky*ky)*v_tilde_imag)=0')
                problem.add_equation('dz(w_hat_imag)-(kx*u_tilde_imag+ky*v_tilde_imag)=0')
                problem.add_equation('dz(p_hat_imag)-(kx*d_u_tilde_imag+ky*d_v_tilde_imag-(kx*kx+ky*ky)*w_hat_imag+Ra_T*T_hat_imag-Ra_S2T*S_hat_imag)=0')
                problem.add_equation('dz(T_hat_imag)-d_T_hat_imag=0')
                problem.add_equation('dz(S_hat_imag)-d_S_hat_imag=0')
                #problem.add_equation('dz(T_0_imag)-d_T_0_imag=0')
                #problem.add_equation('dz(S_0_imag)-d_S_0_imag=0')

                #mean temperature and salinity
                problem.add_equation('dz(T_0)-d_T_0=0')
                problem.add_equation('dz(S_0)-d_S_0=0')
                problem.add_equation('dz(d_T_0)=Pe_T*(2*kx*u_tilde_real*T_hat_real+2*kx*u_tilde_imag*T_hat_imag+2*ky*v_tilde_real*T_hat_real+2*ky*v_tilde_imag*T_hat_imag+2*w_hat_real*d_T_hat_real+2*w_hat_imag*d_T_hat_imag)')
                problem.add_equation('dz(d_S_0)=Pe_S/tau*(2*kx*u_tilde_real*S_hat_real+2*kx*u_tilde_imag*S_hat_imag+2*ky*v_tilde_real*S_hat_real+2*ky*v_tilde_imag*S_hat_imag+2*w_hat_real*d_S_hat_real+2*w_hat_imag*d_S_hat_imag)')

                #coupling between real and imag due to shear
                if self.F_sin=='z':
                    problem.add_equation('dz(d_T_hat_real)-w_hat_real*dy_T_mean-(kx*kx+ky*ky)*T_hat_real+Pe_T*kx*(z-1/2)*T_hat_imag=Pe_T*w_hat_real*d_T_0-Pe_T*kx*eta*T_hat_imag')
                    problem.add_equation('dz(d_T_hat_imag)-w_hat_imag*dy_T_mean-(kx*kx+ky*ky)*T_hat_imag-Pe_T*kx*(z-1/2)*T_hat_real=Pe_T*w_hat_imag*d_T_0+Pe_T*kx*eta*T_hat_real')
                    problem.add_equation('dz(d_S_hat_real)-1/tau*w_hat_real*dy_S_mean-(kx*kx+ky*ky)*S_hat_real+Pe_S/tau*kx*(z-1/2)*S_hat_imag=Pe_S/tau*(w_hat_real*d_S_0)-Pe_S*kx*eta/tau*S_hat_imag')   
                    problem.add_equation('dz(d_S_hat_imag)-1/tau*w_hat_imag*dy_S_mean-(kx*kx+ky*ky)*S_hat_imag-Pe_S/tau*kx*(z-1/2)*S_hat_real=Pe_S/tau*(w_hat_imag*d_S_0)+Pe_S*kx*eta/tau*S_hat_real')   
                else:
                    print(self.F_sin)
                    problem.add_equation('dz(d_T_hat_real)-w_hat_real*dy_T_mean-(kx*kx+ky*ky)*T_hat_real+Pe_T*kx*F_sin*sin(ks*z)*T_hat_imag=Pe_T*w_hat_real*d_T_0-Pe_T*kx*eta*T_hat_imag')
                    problem.add_equation('dz(d_T_hat_imag)-w_hat_imag*dy_T_mean-(kx*kx+ky*ky)*T_hat_imag-Pe_T*kx*F_sin*sin(ks*z)*T_hat_real=Pe_T*w_hat_imag*d_T_0+Pe_T*kx*eta*T_hat_real')
                    problem.add_equation('dz(d_S_hat_real)-1/tau*w_hat_real*dy_S_mean-(kx*kx+ky*ky)*S_hat_real+Pe_S/tau*kx*F_sin*sin(ks*z)*S_hat_imag=Pe_S/tau*(w_hat_real*d_S_0)-Pe_S*kx*eta/tau*S_hat_imag')   
                    problem.add_equation('dz(d_S_hat_imag)-1/tau*w_hat_imag*dy_S_mean-(kx*kx+ky*ky)*S_hat_imag-Pe_S/tau*kx*F_sin*sin(ks*z)*S_hat_real=Pe_S/tau*(w_hat_imag*d_S_0)+Pe_S*kx*eta/tau*S_hat_real')   
                #else:
                #    print("Wrong flag of F_sin.")

                # #integral of the phase condition
                # problem.add_equation('dz(int_phase_cond)-w_hat_imag=0')
                # problem.add_bc('left(int_phase_cond)=0')
                # problem.add_bc('right(int_phase_cond)=0')
                
                #This B.C. works for the periodic B.C. and they can converge smoothly...
                #problem.add_bc("left(T_hat_imag)=0")
    
                problem.add_bc("left(p_hat_imag)=0")
    
            elif self.problem =='IVP':
                #real
                #with the large scale shear
                #Update 2022/05/04
                #
                if self.F_sin=='z':
                    problem.substitutions['U_bg'] = "z-1/2"
                    problem.substitutions['d_U_bg']="1"
                
                else:
                    problem.substitutions['U_bg'] = "F_sin*sin(ks*z)"
                    problem.substitutions['d_U_bg'] = "ks*F_sin*cos(ks*z)"
                 
                problem.add_equation('dz(u_tilde_real)-d_u_tilde_real=0')
                problem.add_equation('-1/Pr*dt(u_tilde_real)+dz(d_u_tilde_real)-(kx*p_hat_real+(kx*kx+ky*ky)*u_tilde_real)=1/Pr*(-kx*(U_bg+U_0)*u_tilde_imag+(d_U_bg+d_U_0)*w_hat_imag)')
                problem.add_equation('dz(v_tilde_real)-d_v_tilde_real=0')
                problem.add_equation('-1/Pr*dt(v_tilde_real)+dz(d_v_tilde_real)-(ky*p_hat_real+(kx*kx+ky*ky)*v_tilde_real)=1/Pr*(-kx*(U_bg+U_0)*v_tilde_imag)')
                problem.add_equation('dz(w_hat_real)-(kx*u_tilde_real+ky*v_tilde_real)=0')
                problem.add_equation('-1/Pr*dt(w_hat_real)-dz(p_hat_real)+(kx*d_u_tilde_real+ky*d_v_tilde_real-(kx*kx+ky*ky)*w_hat_real+Ra_T*T_hat_real-Ra_S2T*S_hat_real)=1/Pr*(-kx*(U_bg+U_0)*w_hat_imag)')
                problem.add_equation('dz(T_hat_real)-d_T_hat_real=0')
                problem.add_equation('dz(S_hat_real)-d_S_hat_real=0')
                
                #IBM: +A*exp(-(z-z0)**2/sigma**2)*w_hat_imag
                #imag
                problem.add_equation('dz(u_tilde_imag)-d_u_tilde_imag=0')
                problem.add_equation('-1/Pr*dt(u_tilde_imag)+dz(d_u_tilde_imag)-(kx*p_hat_imag+(kx*kx+ky*ky)*u_tilde_imag)=1/Pr*(kx*(U_bg+U_0)*u_tilde_real-(d_U_bg+d_U_0)*w_hat_real)')
                problem.add_equation('dz(v_tilde_imag)-d_v_tilde_imag=0')
                problem.add_equation('-1/Pr*dt(v_tilde_imag)+dz(d_v_tilde_imag)-(ky*p_hat_imag+(kx*kx+ky*ky)*v_tilde_imag)=1/Pr*(kx*(U_bg+U_0)*v_tilde_real)')
                problem.add_equation('dz(w_hat_imag)-(kx*u_tilde_imag+ky*v_tilde_imag)=0')
                problem.add_equation('-1/Pr*dt(w_hat_imag)-dz(p_hat_imag)+(kx*d_u_tilde_imag+ky*d_v_tilde_imag-(kx*kx+ky*ky)*w_hat_imag+Ra_T*T_hat_imag-Ra_S2T*S_hat_imag)=1/Pr*(kx*(U_bg+U_0)*w_hat_real)')
                problem.add_equation('dz(T_hat_imag)-d_T_hat_imag=0')
                problem.add_equation('dz(S_hat_imag)-d_S_hat_imag=0')
                #problem.add_equation('dz(T_0_imag)-d_T_0_imag=0')
                #problem.add_equation('dz(S_0_imag)-d_S_0_imag=0')

                #harmonnic of the temperature and salinity
                problem.add_equation('-dt(T_hat_real)+dz(d_T_hat_real)-w_hat_real*dy_T_mean-(kx*kx+ky*ky)*T_hat_real=Pe_T*w_hat_real*d_T_0-Pe_T*kx*(U_bg+U_0)*T_hat_imag')
                problem.add_equation('-dt(T_hat_imag)+dz(d_T_hat_imag)-w_hat_imag*dy_T_mean-(kx*kx+ky*ky)*T_hat_imag=Pe_T*w_hat_imag*d_T_0+Pe_T*kx*(U_bg+U_0)*T_hat_real')
                problem.add_equation('-1/tau*dt(S_hat_real)+dz(d_S_hat_real)-1/tau*w_hat_real*dy_S_mean-(kx*kx+ky*ky)*S_hat_real=Pe_S/tau*(w_hat_real*d_S_0)-Pe_S/tau*kx*(U_bg+U_0)*S_hat_imag')   
                problem.add_equation('-1/tau*dt(S_hat_imag)+dz(d_S_hat_imag)-1/tau*w_hat_imag*dy_S_mean-(kx*kx+ky*ky)*S_hat_imag=Pe_S/tau*(w_hat_imag*d_S_0)+Pe_S/tau*kx*(U_bg+U_0)*S_hat_real')   

                #mean temperature and salinity
                problem.add_equation('dz(T_0)-d_T_0=0')
                problem.add_equation('dz(S_0)-d_S_0=0')
                problem.add_equation('-dt(T_0)+dz(d_T_0)=(2*kx*u_tilde_real*T_hat_real+2*kx*u_tilde_imag*T_hat_imag+2*ky*v_tilde_real*T_hat_real+2*ky*v_tilde_imag*T_hat_imag+2*w_hat_real*d_T_hat_real+2*w_hat_imag*d_T_hat_imag)')
                problem.add_equation('-1/tau*dt(S_0)+dz(d_S_0)=1/tau*(2*kx*u_tilde_real*S_hat_real+2*kx*u_tilde_imag*S_hat_imag+2*ky*v_tilde_real*S_hat_real+2*ky*v_tilde_imag*S_hat_imag+2*w_hat_real*d_S_hat_real+2*w_hat_imag*d_S_hat_imag)')

                #large scale shear U_0
                problem.add_equation('dz(U_0)-d_U_0=0')
                problem.add_equation('-1/Pr*dt(U_0)+dz(d_U_0)=1/Pr*(2*kx*u_tilde_real*(-u_tilde_imag)+2*kx*u_tilde_imag*u_tilde_real+2*ky*v_tilde_real*(-u_tilde_imag)+2*ky*v_tilde_imag*u_tilde_real+2*w_hat_real*(-d_u_tilde_imag)+2*w_hat_imag*d_u_tilde_real)')
                """
                #coupling between real and imag due to shear
                if self.F_sin=='z':
                    problem.add_equation('-dt(T_hat_real)+dz(d_T_hat_real)-w_hat_real*dy_T_mean-(kx*kx+ky*ky)*T_hat_real+Pe_T*kx*(z-1/2)*T_hat_imag=Pe_T*w_hat_real*d_T_0')
                    problem.add_equation('-dt(T_hat_imag)+dz(d_T_hat_imag)-w_hat_imag*dy_T_mean-(kx*kx+ky*ky)*T_hat_imag-Pe_T*kx*(z-1/2)*T_hat_real=Pe_T*w_hat_imag*d_T_0')
                    problem.add_equation('-1/tau*dt(S_hat_real)+dz(d_S_hat_real)-1/tau*w_hat_real*dy_S_mean-(kx*kx+ky*ky)*S_hat_real+Pe_S/tau*kx*(z-1/2)*S_hat_imag=Pe_S/tau*(w_hat_real*d_S_0)')   
                    problem.add_equation('-1/tau*dt(S_hat_imag)+dz(d_S_hat_imag)-1/tau*w_hat_imag*dy_S_mean-(kx*kx+ky*ky)*S_hat_imag-Pe_S/tau*kx*(z-1/2)*S_hat_real=Pe_S/tau*(w_hat_imag*d_S_0)')   
                else:
                    print(self.F_sin)
                    problem.add_equation('-dt(T_hat_real)+dz(d_T_hat_real)-w_hat_real*dy_T_mean-(kx*kx+ky*ky)*T_hat_real+Pe_T*kx*F_sin*sin(ks*z)*T_hat_imag=Pe_T*w_hat_real*d_T_0')
                    problem.add_equation('-dt(T_hat_imag)+dz(d_T_hat_imag)-w_hat_imag*dy_T_mean-(kx*kx+ky*ky)*T_hat_imag-Pe_T*kx*F_sin*sin(ks*z)*T_hat_real=Pe_T*w_hat_imag*d_T_0')
                    problem.add_equation('-1/tau*dt(S_hat_real)+dz(d_S_hat_real)-1/tau*w_hat_real*dy_S_mean-(kx*kx+ky*ky)*S_hat_real+Pe_S/tau*kx*F_sin*sin(ks*z)*S_hat_imag=Pe_S/tau*(w_hat_real*d_S_0)')   
                    problem.add_equation('-1/tau*dt(S_hat_imag)+dz(d_S_hat_imag)-1/tau*w_hat_imag*dy_S_mean-(kx*kx+ky*ky)*S_hat_imag-Pe_S/tau*kx*F_sin*sin(ks*z)*S_hat_real=Pe_S/tau*(w_hat_imag*d_S_0)')   
                """
                
                
                
                """
                #real
                problem.add_equation('dz(u_tilde_real)-d_u_tilde_real=0')
                problem.add_equation('-1/Pr*dt(u_tilde_real)+dz(d_u_tilde_real)-(kx*p_hat_real+(kx*kx+ky*ky)*u_tilde_real)=0')
                problem.add_equation('dz(v_tilde_real)-d_v_tilde_real=0')
                problem.add_equation('-1/Pr*dt(v_tilde_real)+dz(d_v_tilde_real)-(ky*p_hat_real+(kx*kx+ky*ky)*v_tilde_real)=0')
                problem.add_equation('dz(w_hat_real)-(kx*u_tilde_real+ky*v_tilde_real)=0')
                problem.add_equation('1/Pr*dt(w_hat_real)+dz(p_hat_real)-(kx*d_u_tilde_real+ky*d_v_tilde_real-(kx*kx+ky*ky)*w_hat_real+Ra_T*T_hat_real-Ra_S2T*S_hat_real)=0')
                problem.add_equation('dz(T_hat_real)-d_T_hat_real=0')
                problem.add_equation('dz(S_hat_real)-d_S_hat_real=0')
                
                #IBM: +A*exp(-(z-z0)**2/sigma**2)*w_hat_imag
                #imag
                problem.add_equation('dz(u_tilde_imag)-d_u_tilde_imag=0')
                problem.add_equation('-1/Pr*dt(u_tilde_imag)+dz(d_u_tilde_imag)-(kx*p_hat_imag+(kx*kx+ky*ky)*u_tilde_imag)=0')
                problem.add_equation('dz(v_tilde_imag)-d_v_tilde_imag=0')
                problem.add_equation('-1/Pr*dt(v_tilde_imag)+dz(d_v_tilde_imag)-(ky*p_hat_imag+(kx*kx+ky*ky)*v_tilde_imag)=0')
                problem.add_equation('dz(w_hat_imag)-(kx*u_tilde_imag+ky*v_tilde_imag)=0')
                problem.add_equation('1/Pr*dt(w_hat_imag)+dz(p_hat_imag)-(kx*d_u_tilde_imag+ky*d_v_tilde_imag-(kx*kx+ky*ky)*w_hat_imag+Ra_T*T_hat_imag-Ra_S2T*S_hat_imag)=0')
                problem.add_equation('dz(T_hat_imag)-d_T_hat_imag=0')
                problem.add_equation('dz(S_hat_imag)-d_S_hat_imag=0')
                #problem.add_equation('dz(T_0_imag)-d_T_0_imag=0')
                #problem.add_equation('dz(S_0_imag)-d_S_0_imag=0')

                #mean temperature and salinity
                problem.add_equation('dz(T_0)-d_T_0=0')
                problem.add_equation('dz(S_0)-d_S_0=0')
                problem.add_equation('-dt(T_0)+dz(d_T_0)=Pe_T*(2*kx*u_tilde_real*T_hat_real+2*kx*u_tilde_imag*T_hat_imag+2*ky*v_tilde_real*T_hat_real+2*ky*v_tilde_imag*T_hat_imag+2*w_hat_real*d_T_hat_real+2*w_hat_imag*d_T_hat_imag)')
                problem.add_equation('-1/tau*dt(S_0)+dz(d_S_0)=Pe_S/tau*(2*kx*u_tilde_real*S_hat_real+2*kx*u_tilde_imag*S_hat_imag+2*ky*v_tilde_real*S_hat_real+2*ky*v_tilde_imag*S_hat_imag+2*w_hat_real*d_S_hat_real+2*w_hat_imag*d_S_hat_imag)')

                #coupling between real and imag due to shear
                if self.F_sin=='z':
                    problem.add_equation('-dt(T_hat_real)+dz(d_T_hat_real)-w_hat_real*dy_T_mean-(kx*kx+ky*ky)*T_hat_real+Pe_T*kx*(z-1/2)*T_hat_imag=Pe_T*w_hat_real*d_T_0')
                    problem.add_equation('-dt(T_hat_imag)+dz(d_T_hat_imag)-w_hat_imag*dy_T_mean-(kx*kx+ky*ky)*T_hat_imag-Pe_T*kx*(z-1/2)*T_hat_real=Pe_T*w_hat_imag*d_T_0')
                    problem.add_equation('-1/tau*dt(S_hat_real)+dz(d_S_hat_real)-1/tau*w_hat_real*dy_S_mean-(kx*kx+ky*ky)*S_hat_real+Pe_S/tau*kx*(z-1/2)*S_hat_imag=Pe_S/tau*(w_hat_real*d_S_0)')   
                    problem.add_equation('-1/tau*dt(S_hat_imag)+dz(d_S_hat_imag)-1/tau*w_hat_imag*dy_S_mean-(kx*kx+ky*ky)*S_hat_imag-Pe_S/tau*kx*(z-1/2)*S_hat_real=Pe_S/tau*(w_hat_imag*d_S_0)')   
                else:
                    print(self.F_sin)
                    problem.add_equation('-dt(T_hat_real)+dz(d_T_hat_real)-w_hat_real*dy_T_mean-(kx*kx+ky*ky)*T_hat_real+Pe_T*kx*F_sin*sin(ks*z)*T_hat_imag=Pe_T*w_hat_real*d_T_0')
                    problem.add_equation('-dt(T_hat_imag)+dz(d_T_hat_imag)-w_hat_imag*dy_T_mean-(kx*kx+ky*ky)*T_hat_imag-Pe_T*kx*F_sin*sin(ks*z)*T_hat_real=Pe_T*w_hat_imag*d_T_0')
                    problem.add_equation('-1/tau*dt(S_hat_real)+dz(d_S_hat_real)-1/tau*w_hat_real*dy_S_mean-(kx*kx+ky*ky)*S_hat_real+Pe_S/tau*kx*F_sin*sin(ks*z)*S_hat_imag=Pe_S/tau*(w_hat_real*d_S_0)')   
                    problem.add_equation('-1/tau*dt(S_hat_imag)+dz(d_S_hat_imag)-1/tau*w_hat_imag*dy_S_mean-(kx*kx+ky*ky)*S_hat_imag-Pe_S/tau*kx*F_sin*sin(ks*z)*S_hat_real=Pe_S/tau*(w_hat_imag*d_S_0)')   
                
                
                """
                
                
                # ##old time dependent version
                # problem.add_equation('dz(u_tilde)-d_u_tilde=0')
                # problem.add_equation('-1/Pr*dt(u_tilde)+dz(d_u_tilde)-(kx*p_hat+(kx*kx+ky*ky)*u_tilde)=0')
                # problem.add_equation('dz(v_tilde)-d_v_tilde=0')
                # problem.add_equation('-1/Pr*dt(v_tilde)+dz(d_v_tilde)-(ky*p_hat+(kx*kx+ky*ky)*v_tilde)=0')
                # problem.add_equation('dz(w_hat)-(kx*u_tilde+ky*v_tilde)=0')
                # problem.add_equation('1/Pr*dt(w_hat)+dz(p_hat)-(kx*d_u_tilde+ky*d_v_tilde-(kx*kx+ky*ky)*w_hat+Ra_T*T_hat-Ra_S2T*S_hat)=0')
                # problem.add_equation('dz(T_hat)-d_T_hat=0')
                # problem.add_equation('-dt(T_hat)+dz(d_T_hat)-(w_hat*dy_T_mean+(kx*kx+ky*ky)*T_hat)=w_hat*d_T_0')
                # problem.add_equation('dz(S_hat)-d_S_hat=0')
                # problem.add_equation('-1/tau*dt(S_hat)+dz(d_S_hat)-1/tau*w_hat*dy_S_mean-(kx*kx+ky*ky)*S_hat=1/tau*(w_hat*d_S_0)')   
                # problem.add_equation('dz(T_0)-d_T_0=0')
                # problem.add_equation('-dt(T_0)+dz(d_T_0)=2*kx*u_tilde*T_hat+2*ky*v_tilde*T_hat+2*w_hat*d_T_hat')
                # problem.add_equation('dz(S_0)-d_S_0=0')
                # problem.add_equation('-1/tau*dt(S_0)+dz(d_S_0)=1/tau*(2*kx*u_tilde*S_hat+2*ky*v_tilde*S_hat+2*w_hat*d_S_hat)')
            
            if self.problem in ['IVP','BVP']:
                #Add BC for IVP and BVP. 
                if self.z_bc_w_left=='periodic' and self.z_bc_w_right=='periodic':
                    problem.add_bc("left(w_hat_real)-right(w_hat_real)=0")
                    problem.add_bc("left(p_hat_real)-right(p_hat_real)=0")
                    problem.add_bc("left(w_hat_imag)-right(w_hat_imag)=0")
                    problem.add_bc("left(p_hat_imag)-right(p_hat_imag)=0")
    
                    #problem.add_bc("left(w_hat_imag)=0")
                    #problem.add_bc("right(w_hat_imag)=0")
               
                if self.z_bc_w_left=='dirichlet':
                    problem.add_bc("left(w_hat_real)=0")
                    problem.add_bc("left(w_hat_imag)=0")
                    print("Dirichlet for w left")
                    
                if self.z_bc_w_right=='dirichlet':
                    problem.add_bc("right(w_hat_real)=0")
                    problem.add_bc("right(w_hat_imag)=0")
                    print("Dirichlet for w right")
               
                if self.z_bc_T_left=='periodic' and self.z_bc_T_right=='periodic':
                    problem.add_bc("left(T_hat_real)-right(T_hat_real)=0")
                    problem.add_bc("left(d_T_hat_real)-right(d_T_hat_real)=0")
                    problem.add_bc("left(T_hat_imag)-right(T_hat_imag)=0")
                    problem.add_bc("left(d_T_hat_imag)-right(d_T_hat_imag)=0")
                    problem.add_bc("left(T_0)=0")
                    problem.add_bc("right(T_0)=0")
                    #problem.add_bc("left(d_T_0)-right(d_T_0)=0")
                    print("Periodic for T")
                  
                if self.z_bc_T_left=='dirichlet':
                    problem.add_bc("left(T_hat_real)=0")
                    problem.add_bc("left(T_hat_imag)=0")
                    problem.add_bc("left(T_0)=0")
                    print("Dirichlet for T left")
                elif self.z_bc_T_left=='neumann':
                    problem.add_bc("left(d_T_hat_real)=0")
                    problem.add_bc("left(d_T_hat_imag)=0")
                    problem.add_bc("left(d_T_0)=0")
                    print("Neumann for T left")
                    
                if self.z_bc_T_right=='dirichlet':
                    problem.add_bc("right(T_hat_real)=0")
                    problem.add_bc("right(T_hat_imag)=0")
                    problem.add_bc("right(T_0)=0")
                    print("Dirichlet for T right")
                elif self.z_bc_T_right=='neumann':
                    problem.add_bc("right(d_T_hat_real)=0")
                    problem.add_bc("right(d_T_hat_imag)=0")
                    problem.add_bc("right(d_T_0)=0")
                    print("Neumann for T right")
                    
                if self.z_bc_S_left=='periodic' and self.z_bc_S_right=='periodic':
                    problem.add_bc("left(S_hat_real)-right(S_hat_real)=0")
                    problem.add_bc("left(d_S_hat_real)-right(d_S_hat_real)=0")
                    problem.add_bc("left(S_hat_imag)-right(S_hat_imag)=0")
                    problem.add_bc("left(d_S_hat_imag)-right(d_S_hat_imag)=0")
                    problem.add_bc("left(S_0)=0")
                    problem.add_bc("right(S_0)=0")
                    #problem.add_bc("left(d_S_0)-right(d_S_0)=0")
                    print("Periodic for S")
               
                if self.z_bc_S_left=='dirichlet':
                    problem.add_bc("left(S_hat_real)=0")
                    problem.add_bc("left(S_hat_imag)=0")
                    problem.add_bc("left(S_0)=0")
                    print("Dirichlet for S left")
                elif self.z_bc_S_left=='neumann':
                    problem.add_bc("left(d_S_hat_real)=0")
                    problem.add_bc("left(d_S_hat_imag)=0")
                    problem.add_bc("left(d_S_0)=0")
                    print("Neumann for S left")
                    
                if self.z_bc_S_right=='dirichlet':
                    problem.add_bc("right(S_hat_real)=0")
                    problem.add_bc("right(S_hat_imag)=0")
                    problem.add_bc("right(S_0)=0")
                    print("Dirichlet for S right")
                elif self.z_bc_S_right=='neumann':
                    problem.add_bc("right(d_S_hat_real)=0")
                    problem.add_bc("right(d_S_hat_imag)=0")
                    problem.add_bc("right(d_S_0)=0")
                    print("Neumann for S right")
                
                if self.z_bc_u_v_left=='periodic' and self.z_bc_u_v_right=='periodic':
                    problem.add_bc("left(u_tilde_real)-right(u_tilde_real)=0")
                    problem.add_bc("left(v_tilde_real)-right(v_tilde_real)=0")
                    problem.add_bc("left(d_u_tilde_real)-right(d_u_tilde_real)=0")
                    problem.add_bc("left(d_v_tilde_real)-right(d_v_tilde_real)=0")
                    problem.add_bc("left(u_tilde_imag)-right(u_tilde_imag)=0")
                    problem.add_bc("left(v_tilde_imag)-right(v_tilde_imag)=0")
                    problem.add_bc("left(d_u_tilde_imag)-right(d_u_tilde_imag)=0")
                    problem.add_bc("left(d_v_tilde_imag)-right(d_v_tilde_imag)=0")
                    print("Periodic for u,v")
                    problem.add_bc("left(U_0)-right(U_0)=0")
                    problem.add_bc("left(d_U_0)-right(d_U_0)=0")
                    print("Periodic for U_0")
                    
                if self.z_bc_u_v_left=='dirichlet':
                    problem.add_bc("left(u_tilde_real)=0")
                    problem.add_bc("left(u_tilde_imag)=0")
                    problem.add_bc("left(v_tilde_real)=0")
                    problem.add_bc("left(v_tilde_imag)=0")
                    print("Dirichlet for u,v left")
                    problem.add_bc("left(U_0)=0")
                    print("Dirichlet for U_0 left")
                elif self.z_bc_u_v_left=='neumann':
                    problem.add_bc("left(d_u_tilde_real)=0")
                    problem.add_bc("left(d_u_tilde_imag)=0")
                    problem.add_bc("left(d_v_tilde_real)=0")
                    problem.add_bc("left(d_v_tilde_imag)=0")
                    print("Neumann for u,v left")
                    problem.add_bc("left(d_U_0)=0")
                    print("Neumann for U_0 left")
                    
                if self.z_bc_u_v_right=='dirichlet':
                    problem.add_bc("right(u_tilde_real)=0")
                    problem.add_bc("right(u_tilde_imag)=0")
                    problem.add_bc("right(v_tilde_real)=0")
                    problem.add_bc("right(v_tilde_imag)=0")
                    print("Dirichlet for u,v right")
                    problem.add_bc("right(U_0)=0")
                    print("Dirichlet for U_0 right")
                    
                elif self.z_bc_u_v_right=='neumann':
                    problem.add_bc("right(d_u_tilde_real)=0")
                    problem.add_bc("right(d_u_tilde_imag)=0")
                    problem.add_bc("right(d_v_tilde_real)=0")
                    problem.add_bc("right(d_v_tilde_imag)=0")
                    print("Neumann for u,v right")
                    problem.add_bc("right(d_U_0)=0")
                    print("Neumann for U_0 right")


            # elif self.problem=='IVP':
            #     problem.add_equation('dz(u_tilde)-d_u_tilde=0')
            #     problem.add_equation('-1/Pr*dt(u_tilde)+dz(d_u_tilde)-(kx*p_hat+(kx*kx+ky*ky)*u_tilde)=0')
            #     problem.add_equation('dz(v_tilde)-d_v_tilde=0')
            #     problem.add_equation('-1/Pr*dt(v_tilde)+dz(d_v_tilde)-(ky*p_hat+(kx*kx+ky*ky)*v_tilde)=0')
            #     problem.add_equation('dz(w_hat)-(kx*u_tilde+ky*v_tilde)=0')
            #     problem.add_equation('1/Pr*dt(w_hat)+dz(p_hat)-(kx*d_u_tilde+ky*d_v_tilde-(kx*kx+ky*ky)*w_hat+Ra_T*T_hat-Ra_S2T*S_hat)=0')
            #     problem.add_equation('dz(T_hat)-d_T_hat=0')
            #     problem.add_equation('-dt(T_hat)+dz(d_T_hat)-(w_hat*dy_T_mean+(kx*kx+ky*ky)*T_hat)=w_hat*d_T_0')
            #     problem.add_equation('dz(S_hat)-d_S_hat=0')
            #     problem.add_equation('-1/tau*dt(S_hat)+dz(d_S_hat)-1/tau*w_hat*dy_S_mean-(kx*kx+ky*ky)*S_hat=1/tau*(w_hat*d_S_0)')   
            #     problem.add_equation('dz(T_0)-d_T_0=0')
            #     problem.add_equation('-dt(T_0)+dz(d_T_0)=2*kx*u_tilde*T_hat+2*ky*v_tilde*T_hat+2*w_hat*d_T_hat')
            #     problem.add_equation('dz(S_0)-d_S_0=0')
            #     problem.add_equation('-1/tau*dt(S_0)+dz(d_S_0)=1/tau*(2*kx*u_tilde*S_hat+2*ky*v_tilde*S_hat+2*w_hat*d_S_hat)')
            
            elif self.problem =='EVP':
                #real
                problem.add_equation('dz(u_tilde)-d_u_tilde=0')
                problem.add_equation('dz(d_u_tilde)-(kx*p_hat+(kx*kx+ky*ky)*u_tilde)=0')
                problem.add_equation('dz(v_tilde)-d_v_tilde=0')
                problem.add_equation('dz(d_v_tilde)-(ky*p_hat+(kx*kx+ky*ky)*v_tilde)=0')
                problem.add_equation('dz(w_hat)-(kx*u_tilde+ky*v_tilde)=0')
                problem.add_equation('dz(p_hat)-(kx*d_u_tilde+ky*d_v_tilde-(kx*kx+ky*ky)*w_hat+Ra_T*T_hat-Ra_S2T*S_hat)=0')
                problem.add_equation('dz(T_hat)-d_T_hat=0')
                problem.add_equation('dz(S_hat)-d_S_hat=0')
                
                #coupling between real and imag due to shear
                if self.F_sin=='z':
                    problem.add_equation('dz(d_T_hat)-w_hat*dy_T_mean-(kx*kx+ky*ky)*T_hat-1j*Pe_T*kx*(z-1/2)*T_hat-Pe_T*eig_val*T_hat=0')
                    problem.add_equation('dz(d_S_hat)-1/tau*w_hat*dy_S_mean-(kx*kx+ky*ky)*S_hat-1j*Pe_S/tau*kx*(z-1/2)*S_hat-Pe_S*eig_val*S_hat=0')   
                else:
                    print(self.F_sin)
                    problem.add_equation('dz(d_T_hat)-w_hat*dy_T_mean-(kx*kx+ky*ky)*T_hat-1j*Pe_T*kx*F_sin*sin(ks*z)*T_hat-Pe_T*eig_val*T_hat=0')
                    problem.add_equation('dz(d_S_hat)-1/tau*w_hat*dy_S_mean-(kx*kx+ky*ky)*S_hat-1j*Pe_S/tau*kx*F_sin*sin(ks*z)*S_hat-Pe_S/tau*eig_val*S_hat=0')   
                
                # if self.F_sin=='z':
                #     problem.add_equation('dz(d_T_hat)/Pe_T-w_hat*dy_T_mean/Pe_T-(kx*kx+ky*ky)*T_hat/Pe_T-1j*kx*(z-1/2)*T_hat-eig_val*T_hat=0')
                #     problem.add_equation('dz(d_S_hat)/Pe_S-1/tau*w_hat*dy_S_mean/Pe_S-(kx*kx+ky*ky)*S_hat/Pe_S-1j/tau*kx*(z-1/2)*S_hat-eig_val*S_hat=0')   
                # else:
                #     print(self.F_sin)
                #     problem.add_equation('dz(d_T_hat)/Pe_T-w_hat*dy_T_mean/Pe_T-(kx*kx+ky*ky)*T_hat/Pe_T-1j*kx*F_sin*sin(ks*z)*T_hat-eig_val*T_hat=0')
                #     problem.add_equation('dz(d_S_hat)/Pe_S-1/tau*w_hat*dy_S_mean/Pe_S-(kx*kx+ky*ky)*S_hat/Pe_S-1j/tau*kx*F_sin*sin(ks*z)*S_hat-eig_val*S_hat=0')   
                
                
                if self.z_bc_w_left=='periodic' and self.z_bc_w_right=='periodic':
                    problem.add_bc("left(w_hat)-right(w_hat)=0")
                    problem.add_bc("left(p_hat)-right(p_hat)=0")
    
                    #problem.add_bc("left(w_hat_imag)=0")
                    #problem.add_bc("right(w_hat_imag)=0")
               
                if self.z_bc_w_left=='dirichlet':
                    problem.add_bc("left(w_hat)=0")
                    print("Dirichlet for w left")
                    
                if self.z_bc_w_right=='dirichlet':
                    problem.add_bc("right(w_hat)=0")
                    print("Dirichlet for w right")
               
                if self.z_bc_T_left=='periodic' and self.z_bc_T_right=='periodic':
                    problem.add_bc("left(T_hat)-right(T_hat)=0")
                    problem.add_bc("left(d_T_hat)-right(d_T_hat)=0")
                    
                    #problem.add_bc("left(d_T_0)-right(d_T_0)=0")
                    print("Periodic for T")
                  
                if self.z_bc_T_left=='dirichlet':
                    problem.add_bc("left(T_hat)=0")
                    print("Dirichlet for T left")
                elif self.z_bc_T_left=='neumann':
                    problem.add_bc("left(d_T_hat)=0")
                    print("Neumann for T left")
                    
                if self.z_bc_T_right=='dirichlet':
                    problem.add_bc("right(T_hat)=0")
                    print("Dirichlet for T right")
                elif self.z_bc_T_right=='neumann':
                    problem.add_bc("right(d_T_hat)=0")
                    print("Neumann for T right")
                    
                if self.z_bc_S_left=='periodic' and self.z_bc_S_right=='periodic':
                    problem.add_bc("left(S_hat)-right(S_hat)=0")
                    problem.add_bc("left(d_S_hat)-right(d_S_hat)=0")
                    #problem.add_bc("left(d_S_0)-right(d_S_0)=0")
                    print("Periodic for S")
               
                if self.z_bc_S_left=='dirichlet':
                    problem.add_bc("left(S_hat)=0")
                    print("Dirichlet for S left")
                elif self.z_bc_S_left=='neumann':
                    problem.add_bc("left(d_S_hat)=0")
                    print("Neumann for S left")
                    
                if self.z_bc_S_right=='dirichlet':
                    problem.add_bc("right(S_hat)=0")
                    print("Dirichlet for S right")
                elif self.z_bc_S_right=='neumann':
                    problem.add_bc("right(d_S_hat)=0")
                    print("Neumann for S right")
                
                if self.z_bc_u_v_left=='periodic' and self.z_bc_u_v_right=='periodic':
                    problem.add_bc("left(u_tilde)-right(u_tilde)=0")
                    problem.add_bc("left(v_tilde)-right(v_tilde)=0")
                    problem.add_bc("left(d_u_tilde)-right(d_u_tilde)=0")
                    problem.add_bc("left(d_v_tilde)-right(d_v_tilde)=0")
                    print("Periodic for u,v")
    
                if self.z_bc_u_v_left=='dirichlet':
                    problem.add_bc("left(u_tilde)=0")
                    problem.add_bc("left(v_tilde)=0")
                    print("Dirichlet for u,v left")
    
                elif self.z_bc_u_v_left=='neumann':
                    problem.add_bc("left(d_u_tilde)=0")
                    problem.add_bc("left(d_v_tilde)=0")
                    print("Neumann for u,v left")
    
                if self.z_bc_u_v_right=='dirichlet':
                    problem.add_bc("right(u_tilde)=0")
                    problem.add_bc("right(v_tilde)=0")
                    print("Dirichlet for u,v right")
    
                elif self.z_bc_u_v_right=='neumann':
                    problem.add_bc("right(d_u_tilde)=0")
                    problem.add_bc("right(d_v_tilde)=0")
                    print("Neumann for u,v right")
                              
            

            # if self.z_bc_w_left=='dirichlet' and self.z_bc_w_right=='dirichlet':
            #     problem.add_bc("left(p_hat_imag)=0")
            # else:
            #     problem.add_bc("left(w_hat_imag)=0")
            

            #elif self.z_bc_u_v =='periodic':
                #need to to nothing for periodic BC. but change the basis as Fourier at the beginning    
        elif self.flow=='test_periodic':
            #problem = de.NLBVP(domain, variables=['T_hat','d_T_hat','w_hat','p_hat'])
            #problem = de.NLBVP(domain,variables=['T_hat','d_T_hat'])
            problem = de.NLBVP(domain, variables=['T_hat','d_T_hat','w_hat','p_hat','u_tilde','d_u_tilde','v_tilde','d_v_tilde','T_0','d_T_0'])
            problem.parameters['F_sin']=self.F_sin
            problem.parameters['ks']=self.ks
            problem.parameters['kx']=self.kx
            problem.parameters['ky']=self.ky
            problem.parameters['Ra_T']=self.Ra_T
            #problem.parameters['Ra_S2T']=self.Ra_S2T
            problem.parameters['dy_T_mean']=self.dy_T_mean
            #
            problem.add_equation('dz(T_hat)-d_T_hat=0')
            problem.add_equation('dz(d_T_hat)-w_hat*dy_T_mean-(kx*kx+ky*ky)*T_hat=F_sin*sin(ks*z)*T_hat+w_hat*d_T_0')
            #problem.add_equation('dz(w_hat)-p_hat=0')
            #problem.add_equation('dz(p_hat)+w_hat-T_hat=0')
            problem.add_equation('dz(u_tilde)-d_u_tilde=0')
            problem.add_equation('dz(d_u_tilde)-(kx*p_hat+(kx*kx+ky*ky)*u_tilde)=0')
            problem.add_equation('dz(v_tilde)-d_v_tilde=0')
            problem.add_equation('dz(d_v_tilde)-(ky*p_hat+(kx*kx+ky*ky)*v_tilde)=0')
            problem.add_equation('dz(w_hat)-(kx*u_tilde+ky*v_tilde)=0')
            problem.add_equation('dz(p_hat)-(kx*d_u_tilde+ky*d_v_tilde-(kx*kx+ky*ky)*w_hat+Ra_T*T_hat)=0')
            #problem.add_equation('dz(T_hat)-d_T_hat=0')
            #problem.add_equation('dz(d_T_hat)-w_hat*dy_T_mean-(kx*kx+ky*ky)*T_hat-Pe_T*j*kx*F_sin*sin(ks*z)*T_hat=Pe_T*w_hat*d_T_0')
            #problem.add_equation('dz(S_hat)-d_S_hat=0')
            #problem.add_equation('dz(d_S_hat)-1/tau*w_hat*dy_S_mean-(kx*kx+ky*ky)*S_hat-Pe_S/tau*j*kx*F_sin*sin(ks*z)*S_hat=Pe_S/tau*(w_hat*d_S_0)')   
            problem.add_equation('dz(T_0)-d_T_0=0')
            problem.add_equation('dz(d_T_0)=(2*kx*u_tilde*T_hat+2*ky*v_tilde*T_hat+2*w_hat*d_T_hat)')
            #problem.add_equation('dz(S_0)-d_S_0=0')
            #problem.add_equation('dz(d_S_0)=Pe_S/tau*(2*kx*u_tilde*S_hat+2*ky*v_tilde*S_hat+2*w_hat*d_S_hat)')
            #problem.add_bc('left(T_0)-right(T_0)=0')
            #problem.add_bc('left(d_T_0)-right(d_T_0)=0')
            problem.add_bc('left(T_0)=0')
            problem.add_bc('right(T_0)=0')
            problem.add_bc('left(T_hat)-right(T_hat)=0')
            problem.add_bc('left(d_T_hat)-right(d_T_hat)=0')
            problem.add_bc('left(w_hat)-right(w_hat)=0')
            problem.add_bc('left(p_hat)-right(p_hat)=0')
            problem.add_bc('left(u_tilde)-right(u_tilde)=0')
            problem.add_bc('left(d_u_tilde)-right(d_u_tilde)=0')
            problem.add_bc('left(v_tilde)-right(v_tilde)=0')
            problem.add_bc('left(d_v_tilde)-right(d_v_tilde)=0')
        else:
            raise TypeError('flag.flow is not defined yet') 
        
        if self.problem =='IVP':
            if self.timesteppers == 'RK443':
                ts = de.timesteppers.RK443
            elif self.timesteppers == 'RK222':
                ts = de.timesteppers.RK222
            elif self.timesteppers == 'RK111':
                ts = de.timesteppers.RK111
            elif self.timesteppers == 'RKGFY':
                ts = de.timesteppers.RKGFY
            elif self.timesteppers == 'RKSMR':
                ts = de.timestepers.RKSMR
            elif self.timesteppers == 'SBDF1':
                ts = de.timesteppers.SBDF1
            elif self.timesteppers == 'SBDF2':
                ts = de.timesteppers.SBDF2
            elif self.timesteppers == 'SBDF3':
                ts = de.timesteppers.SBDF3
            elif self.timesteppers == 'SBDF4':
                ts = de.timesteppers.SBDF4
                
            
            else:
                raise TypeError('flag.timesteppers is not defined yet') 
            solver =  problem.build_solver(ts)
            
        elif self.problem in ['BVP','EVP']:
            solver =  problem.build_solver()
            
        return solver

    def initial_condition(self,domain,solver):
        #This function setup the initial condition for IVP and intitial guess for BVP
        #if self.problem == 'EVP':
        #    print('Solving eigenvalue problem, no need for initial condition or initial guess')
        #else:  
        if not pathlib.Path('restart.h5').exists():
            print('setup initial condition')
            #This initial condition also need to be modified
            if self.flow in ['IFSC_2D','double_diffusive_2D','double_diffusive_shear_2D']:
        
                #if self.flux_T:
                #    dy_T_mean_q=solver.state['dy_T_mean_q']
                #    dy_T_mean_q['g']=self.dy_T_mean
        
                x = domain.grid(0)
                z = domain.grid(1)
                u = solver.state['u']
                w = solver.state['w']
                S = solver.state['S']
                p = solver.state['p']
                T = solver.state['T']
                
                gshape = domain.dist.grid_layout.global_shape(scales=1)
                slices = domain.dist.grid_layout.slices(scales=1)
                #print(gshape)
                
                rand = np.random.RandomState(seed=23)
                noise = rand.standard_normal(gshape)[slices]
                
                ##Add the random noise
                u0=self.A_noise*noise
                w0=self.A_noise*noise
                S0=self.A_noise*noise
                T0=self.A_noise*noise
                p0=self.A_noise*noise
                #Add the background shear
                u0 = u0+ self.A_shear*self.F_sin/self.ks**2*np.sin(self.ks*z)\
                    + self.A_shear*self.F_sin_2ks/(2*self.ks)**2*np.sin(2*self.ks*z+self.phase_2ks) \
                    + self.A_shear*self.F_sin_3ks/(3*self.ks)**2*np.sin(3*self.ks*z+self.phase_3ks) \
                    + self.A_shear*self.F_sin_4ks/(4*self.ks)**2*np.sin(4*self.ks*z+self.phase_4ks)
                
                #Compute the eigenvalue problem of elevator mode and add the elevator mode into initial condition
                if self.flow=='IFSC_2D':
                    k_opt=(1/2*(-2-self.Ra_ratio+np.sqrt(self.Ra_ratio**2+8*self.Ra_ratio)))**(1/4)
        
                    w0 =w0 +self.A_elevator*np.sin(k_opt*x)
          
                    S0 =S0 -1/self.Ra_ratio*(k_opt**2+1/k_opt**2)*self.A_elevator*np.sin(k_opt*x)
                    T0 =T0 -1/(k_opt**2)*self.A_elevator*np.sin(k_opt*x)
                elif self.flow in ['double_diffusive_2D']:
                    k2=self.k_elevator**2
                    
                    #This A is the linear system for the eigenvalue problem... We already take the elevator mode, so that k_y=0 (vertical wavenumber)
                    #This should work for both finger regime and diffusive regime..
                    A=[[-k2*self.Pr, self.Pr, -self.Pr/self.R_rho_T2S],
                        [-self.dy_T_mean, -k2, 0],
                        [-self.dy_S_mean, 0, -self.tau*k2]];
                    #Compute eigenvalue and eigenvector of 
                    
                    eig_val,eig_vec=linalg.eig(A) #use linear algebra package to compute eigenvalue
                    eig_val_max_ind=np.argmax(np.real(eig_val)) #compute the index of the eigenvalue
                    eig_vec_max=eig_vec[:,eig_val_max_ind] #get the corresponding eigen vector
                    
                    self.lambda_elevator=eig_val[eig_val_max_ind]
                    w0 =w0 + self.A_elevator*np.real(np.exp(1j*self.k_elevator*x))*np.real(eig_vec_max[0]) #set the results weighted by the corresponding eigenvector 
                    T0 =T0 + self.A_elevator*np.real(np.exp(1j*self.k_elevator*x))*np.real(eig_vec_max[1])
                    S0 =S0 + self.A_elevator*np.real(np.exp(1j*self.k_elevator*x))*np.real(eig_vec_max[2])
                    #print(w0)
                    ##This is sample code to setup the initial condition as whatever I want...
                    #S0=np.array([1,2,3,4],[5,6,7,8],[9,10,11,12],[13,14,15,16])
                    #S['g']=S0[slices]
                elif self.flow in ['double_diffusive_shear_2D']:
                    k2=self.k_elevator**2
                    
                    #This A is the linear system for the eigenvalue problem... We already take the elevator mode, so that k_y=0 (vertical wavenumber)
                    #This should work for both finger regime and diffusive regime..
                    #This is the eigenvalue problem only for no shear case....
                    A=[[-k2, self.Ra_T, -self.Ra_S2T],
                        [-self.dy_T_mean, -k2, 0],
                        [-self.dy_S_mean, 0, -self.tau*k2]]
                    B=[[self.Re, 0,0],
                       [0,self.Pe_T, 0],
                       [0,0,self.Pe_S]];
                    #Compute eigenvalue and eigenvector of 
                    
                    #Update 2021/10/05, use the scipy package of the linalg. This can also solve the generalized eigenvalue problem
                    eig_val,eig_vec=linalg.eig(A,B) #use linear algebra package to compute eigenvalue
                    eig_val[np.isinf(eig_val)]=-np.inf
                    eig_val_max_ind=np.argmax(np.real(eig_val)) #compute the index of the eigenvalue
                    eig_vec_max=eig_vec[:,eig_val_max_ind] #get the corresponding eigen vector
                    
                    self.lambda_elevator=eig_val[eig_val_max_ind]
                    w0 =w0 + self.A_elevator*np.real(np.exp(1j*self.k_elevator*x))*np.real(eig_vec_max[0]) #set the results weighted by the corresponding eigenvector 
                    T0 =T0 + self.A_elevator*np.real(np.exp(1j*self.k_elevator*x))*np.real(eig_vec_max[1])
                    S0 =S0 + self.A_elevator*np.real(np.exp(1j*self.k_elevator*x))*np.real(eig_vec_max[2])
                    
                    #Update the elevator mode analytically
                    if self.flux_T:
                        w0=self.A_elevator*np.real(np.exp(1j*self.k_elevator*x))
                        T0=self.k_elevator**2/self.Ra_T*self.A_elevator*np.real(np.exp(1j*self.k_elevator*x))
                        
                    
                    #print(w0)
                    ##This is sample code to setup the initial condition as whatever I want...
                    #S0=np.array([1,2,3,4],[5,6,7,8],[9,10,11,12],[13,14,15,16])
                    #S['g']=S0[slices]
                
                #superpose the field of secondary mode that has variation in vertical direction
                #Right now just superpose a simple one...
                
                #u0 =u0 + self.A_secondary*np.real(np.exp(1j*self.k_secondary*z))
                #w0 =w0 + self.A_secondary*np.real(np.exp(1j*self.k_secondary*z)) #set the results weighted by the corresponding eigenvector 
                u0 = u0+self.A_secondary_U0*np.real(np.exp(1j*self.k_secondary*z))
                w0 = w0+self.A_secondary_w*np.real(np.exp(1j*self.k_secondary*z))
                T0 = T0 + self.A_secondary_T*np.real(np.exp(1j*self.k_secondary*z))
                #S0 =S0 + self.A_secondary_S*np.real(np.exp(1j*self.k_secondary*z))
                S0 = S0 + self.A_secondary_S*np.sin(self.k_secondary*z)    
                
                u['g']=u0
                w['g']=w0
                T['g']=T0
                S['g']=S0
                p['g']=p0
            
            elif self.flow in ['porous_media_2D']:
                #not fully benchmarked!!
                x = domain.grid(0)
                z = domain.grid(1)
                u = solver.state['u']
                w = solver.state['w']
                p = solver.state['p']
                T = solver.state['T']
                
                gshape = domain.dist.grid_layout.global_shape(scales=1)
                slices = domain.dist.grid_layout.slices(scales=1)
                rand = np.random.RandomState(seed=23)
                noise = rand.standard_normal(gshape)[slices]
                
                ##Add the random noise
                u0=self.A_noise*noise
                w0=self.A_noise*noise
                T0=self.A_noise*noise
                p0=self.A_noise*noise
                
                k2=self.k_elevator**2
                    
                #The matrix to compute the eigenvalue of the porous media convection...
                A=[[-1, self.Ra_T],
                    [1, -k2] ]
                B=[[0, 0],
                   [0, 1]];
                
                eig_val,eig_vec=linalg.eig(A,B) #use linear algebra package to compute eigenvalue
                eig_val[np.isinf(eig_val)]=-np.inf
                eig_val_max_ind=np.argmax(np.real(eig_val)) #compute the index of the eigenvalue
                eig_vec_max=eig_vec[:,eig_val_max_ind] #get the corresponding eigen vector
                
                self.lambda_elevator=eig_val[eig_val_max_ind]
                w0 =w0 + self.A_elevator*np.real(np.exp(1j*self.k_elevator*x))*np.real(eig_vec_max[0]) #set the results weighted by the corresponding eigenvector 
                T0 =T0 + self.A_elevator*np.real(np.exp(1j*self.k_elevator*x))*np.real(eig_vec_max[1])
                
                T0 =T0 + self.A_secondary_T*np.real(np.exp(1j*self.k_secondary*z))

                u['g']=u0
                w['g']=w0
                T['g']=T0
                p['g']=p0
                
            elif self.flow =='HB_porous':
                #initial guess
                z = domain.grid(0)

                #initial guess for the HB_porous, harmonic balance method for double-diffusive convection within porous media
                w_hat = solver.state['w_hat']
                p_hat = solver.state['p_hat']
                T_hat = solver.state['T_hat']
                d_T_hat = solver.state['d_T_hat']
                S_hat = solver.state['S_hat']
                d_S_hat = solver.state['d_S_hat']
                T_0 = solver.state['T_0']
                d_T_0 = solver.state['d_T_0']
                S_0 = solver.state['S_0']
                d_S_0 = solver.state['d_S_0']
                
                
                W0=self.A_elevator;
                gshape = domain.dist.grid_layout.global_shape(scales=1)
                slices = domain.dist.grid_layout.slices(scales=1)
                rand = np.random.RandomState(seed=23)
                noise = rand.standard_normal(gshape)[slices]
                if self.z_bc_T_left=='periodic' and self.z_bc_S_left=='periodic' and self.z_bc_w_left=='periodic':
                    #periodic B.C.
                    w_hat['g'] = W0 +self.A_noise*noise
                    p_hat['g'] = W0/(-(self.kx*self.kx+self.ky*self.ky))+self.A_noise*noise
                    T_hat['g'] = 1/(-(self.kx*self.kx+self.ky*self.ky))*self.dy_T_mean*W0+self.A_noise*noise
                    d_T_hat['g'] =self.A_noise*noise
                    S_hat['g'] = 1/(-(self.kx*self.kx+self.ky*self.ky))*self.dy_S_mean/self.tau*W0+self.A_noise*noise
                    d_S_hat['g'] =self.A_noise*noise
                    T_0['g'] = self.A_noise*noise
                    d_T_0['g'] = self.A_noise*noise
                    S_0['g'] = self.A_noise*noise
                    d_S_0['g'] = self.A_noise*noise
                    
                else:
                    #This is for the other B.C. like the 
                    w_hat['g'] = W0*np.sin(np.pi*z) +self.A_noise*noise
                    p_hat['g'] = W0*np.pi*np.cos(np.pi*z)/(-(self.kx*self.kx+self.ky*self.ky))+self.A_noise*noise
                    T_hat['g'] = 1/(-np.pi**2-(self.kx*self.kx+self.ky*self.ky))*self.dy_T_mean*W0*np.sin(np.pi*z)+self.A_noise*noise
                    d_T_hat['g'] =1/(-np.pi**2-(self.kx*self.kx+self.ky*self.ky))*self.dy_T_mean* W0*np.pi*np.cos(np.pi*z)+self.A_noise*noise
                    S_hat['g'] = 1/(-np.pi**2-(self.kx*self.kx+self.ky*self.ky))*self.dy_S_mean/self.tau*W0*np.sin(np.pi*z)+self.A_noise*noise
                    d_S_hat['g'] =1/(-np.pi**2-(self.kx*self.kx+self.ky*self.ky))*self.dy_S_mean/self.tau* W0*np.pi*np.cos(np.pi*z)+self.A_noise*noise
                    T_0['g'] = self.A_noise*noise
                    d_T_0['g'] = self.A_noise*noise
                    S_0['g'] = self.A_noise*noise
                    d_S_0['g'] = self.A_noise*noise
                    
                    if not (self.kx_2==0 and self.ky_2==0):
                        w_hat_2 = solver.state['w_hat_2']
                        p_hat_2 = solver.state['p_hat_2']
                        T_hat_2 = solver.state['T_hat_2']
                        d_T_hat_2 = solver.state['d_T_hat_2']
                        S_hat_2 = solver.state['S_hat_2']
                        d_S_hat_2 = solver.state['d_S_hat_2']
                        w_hat_2['g'] = W0*np.sin(np.pi*z) +self.A_noise*noise
                        p_hat_2['g'] = W0*np.pi*np.cos(np.pi*z)/(-(self.kx_2*self.kx_2+self.ky_2*self.ky_2))+self.A_noise*noise
                        T_hat_2['g'] = 1/(-np.pi**2-(self.kx_2*self.kx_2+self.ky_2*self.ky_2))*self.dy_T_mean*W0*np.sin(np.pi*z)+self.A_noise*noise
                        d_T_hat_2['g'] =1/(-np.pi**2-(self.kx_2*self.kx_2+self.ky_2*self.ky_2))*self.dy_T_mean* W0*np.pi*np.cos(np.pi*z)+self.A_noise*noise
                        S_hat_2['g'] = 1/(-np.pi**2-(self.kx_2*self.kx_2+self.ky_2*self.ky_2))*self.dy_S_mean/self.tau*W0*np.sin(np.pi*z)+self.A_noise*noise
                        d_S_hat_2['g'] =1/(-np.pi**2-(self.kx_2*self.kx_2+self.ky_2*self.ky_2))*self.dy_S_mean/self.tau* W0*np.pi*np.cos(np.pi*z)+self.A_noise*noise
                   
            elif self.flow =='HB_porous_2_layer':
                #initial guess
                z = domain.grid(0)

                #initial guess for the HB_porous, harmonic balance method for double-diffusive convection within porous media
                w_hat = solver.state['w_hat']
                p_hat = solver.state['p_hat']
                T_hat = solver.state['T_hat']
                d_T_hat = solver.state['d_T_hat']
                S_hat = solver.state['S_hat']
                d_S_hat = solver.state['d_S_hat']
                T_0 = solver.state['T_0']
                d_T_0 = solver.state['d_T_0']
                S_0 = solver.state['S_0']
                d_S_0 = solver.state['d_S_0']
                
                w_hat_top = solver.state['w_hat_top']
                p_hat_top = solver.state['p_hat_top']
                T_hat_top = solver.state['T_hat_top']
                d_T_hat_top = solver.state['d_T_hat_top']
                S_hat_top = solver.state['S_hat_top']
                d_S_hat_top = solver.state['d_S_hat_top']
                T_0_top = solver.state['T_0_top']
                d_T_0_top = solver.state['d_T_0_top']
                S_0_top = solver.state['S_0_top']
                d_S_0_top = solver.state['d_S_0_top']
                
                
                W0=self.A_elevator;
                gshape = domain.dist.grid_layout.global_shape(scales=1)
                slices = domain.dist.grid_layout.slices(scales=1)
                rand = np.random.RandomState(seed=23)
                noise = rand.standard_normal(gshape)[slices]
            
                #This is for the other B.C. like the 
                w_hat['g'] = W0*np.sin(np.pi*z) +self.A_noise*noise
                p_hat['g'] = W0*np.pi*np.cos(np.pi*z)/(-(self.kx*self.kx+self.ky*self.ky))+self.A_noise*noise
                T_hat['g'] = 1/(-np.pi**2-(self.kx*self.kx+self.ky*self.ky))*self.dy_T_mean*W0*np.sin(np.pi*z)+self.A_noise*noise
                d_T_hat['g'] =1/(-np.pi**2-(self.kx*self.kx+self.ky*self.ky))*self.dy_T_mean* W0*np.pi*np.cos(np.pi*z)+self.A_noise*noise
                S_hat['g'] = 1/(-np.pi**2-(self.kx*self.kx+self.ky*self.ky))*self.dy_S_mean/self.tau*W0*np.sin(np.pi*z)+self.A_noise*noise
                d_S_hat['g'] =1/(-np.pi**2-(self.kx*self.kx+self.ky*self.ky))*self.dy_S_mean/self.tau* W0*np.pi*np.cos(np.pi*z)+self.A_noise*noise
                T_0['g'] = self.A_noise*noise
                d_T_0['g'] = self.A_noise*noise
                S_0['g'] = self.A_noise*noise
                d_S_0['g'] = self.A_noise*noise
                
                  
                #This is for the other B.C. like the 
                z_top=z+self.Lz
                w_hat_top['g'] = W0*np.sin(np.pi*z_top) +self.A_noise*noise
                p_hat_top['g'] = W0*np.pi*np.cos(np.pi*z_top)/(-(self.kx*self.kx+self.ky*self.ky))+self.A_noise*noise
                T_hat_top['g'] = 1/(-np.pi**2-(self.kx*self.kx+self.ky*self.ky))*self.dy_T_mean*W0*np.sin(np.pi*z_top)+self.A_noise*noise
                d_T_hat_top['g'] =1/(-np.pi**2-(self.kx*self.kx+self.ky*self.ky))*self.dy_T_mean* W0*np.pi*np.cos(np.pi*z_top)+self.A_noise*noise
                S_hat_top['g'] = 1/(-np.pi**2-(self.kx*self.kx+self.ky*self.ky))*self.dy_S_mean/self.tau*W0*np.sin(np.pi*z_top)+self.A_noise*noise
                d_S_hat_top['g'] =1/(-np.pi**2-(self.kx*self.kx+self.ky*self.ky))*self.dy_S_mean/self.tau* W0*np.pi*np.cos(np.pi*z_top)+self.A_noise*noise
                T_0_top['g'] = self.A_noise*noise
                d_T_0_top['g'] = self.A_noise*noise
                S_0_top['g'] = self.A_noise*noise
                d_S_0_top['g'] = self.A_noise*noise
                
            elif self.flow=='HB_porous_3_layer':
                #initial guess
                z = domain.grid(0)

                #initial guess for the HB_porous, harmonic balance method for double-diffusive convection within porous media
                w_hat = solver.state['w_hat']
                p_hat = solver.state['p_hat']
                T_hat = solver.state['T_hat']
                d_T_hat = solver.state['d_T_hat']
                S_hat = solver.state['S_hat']
                d_S_hat = solver.state['d_S_hat']
                T_0 = solver.state['T_0']
                d_T_0 = solver.state['d_T_0']
                S_0 = solver.state['S_0']
                d_S_0 = solver.state['d_S_0']
                
                w_hat_mid = solver.state['w_hat_mid']
                p_hat_mid = solver.state['p_hat_mid']
                T_hat_mid = solver.state['T_hat_mid']
                d_T_hat_mid = solver.state['d_T_hat_mid']
                S_hat_mid = solver.state['S_hat_mid']
                d_S_hat_mid = solver.state['d_S_hat_mid']
                T_0_mid = solver.state['T_0_mid']
                d_T_0_mid = solver.state['d_T_0_mid']
                S_0_mid = solver.state['S_0_mid']
                d_S_0_mid = solver.state['d_S_0_mid']
                
                w_hat_top = solver.state['w_hat_top']
                p_hat_top = solver.state['p_hat_top']
                T_hat_top = solver.state['T_hat_top']
                d_T_hat_top = solver.state['d_T_hat_top']
                S_hat_top = solver.state['S_hat_top']
                d_S_hat_top = solver.state['d_S_hat_top']
                T_0_top = solver.state['T_0_top']
                d_T_0_top = solver.state['d_T_0_top']
                S_0_top = solver.state['S_0_top']
                d_S_0_top = solver.state['d_S_0_top']
                
                
                W0=self.A_elevator;
                gshape = domain.dist.grid_layout.global_shape(scales=1)
                slices = domain.dist.grid_layout.slices(scales=1)
                rand = np.random.RandomState(seed=23)
                noise = rand.standard_normal(gshape)[slices]
            
                h=self.HB_porous_3_layer_h
                z_bot=(1-h)/2/self.Lz*z
                #This is for the other B.C. like the 
                w_hat['g'] = W0*np.sin(np.pi*z_bot) +self.A_noise*noise
                p_hat['g'] = W0*np.pi*np.cos(np.pi*z_bot)/(-(self.kx*self.kx+self.ky*self.ky))+self.A_noise*noise
                T_hat['g'] = 1/(-np.pi**2-(self.kx*self.kx+self.ky*self.ky))*self.dy_T_mean*W0*np.sin(np.pi*z_bot)+self.A_noise*noise
                d_T_hat['g'] =1/(-np.pi**2-(self.kx*self.kx+self.ky*self.ky))*self.dy_T_mean* W0*np.pi*np.cos(np.pi*z_bot)+self.A_noise*noise
                S_hat['g'] = 1/(-np.pi**2-(self.kx*self.kx+self.ky*self.ky))*self.dy_S_mean/self.tau*W0*np.sin(np.pi*z_bot)+self.A_noise*noise
                d_S_hat['g'] =1/(-np.pi**2-(self.kx*self.kx+self.ky*self.ky))*self.dy_S_mean/self.tau* W0*np.pi*np.cos(np.pi*z_bot)+self.A_noise*noise
                T_0['g'] = self.A_noise*noise
                d_T_0['g'] = self.A_noise*noise
                S_0['g'] = self.A_noise*noise
                d_S_0['g'] = self.A_noise*noise
                
                #This is for the other B.C. like the 
                z_mid=(1-h)/2+h/self.Lz*z
                w_hat_mid['g'] = W0*np.sin(np.pi*z_mid) +self.A_noise*noise
                p_hat_mid['g'] = W0*np.pi*np.cos(np.pi*z_mid)/(-(self.kx*self.kx+self.ky*self.ky))+self.A_noise*noise
                T_hat_mid['g'] = 1/(-np.pi**2-(self.kx*self.kx+self.ky*self.ky))*self.dy_T_mean*W0*np.sin(np.pi*z_mid)+self.A_noise*noise
                d_T_hat_mid['g'] =1/(-np.pi**2-(self.kx*self.kx+self.ky*self.ky))*self.dy_T_mean* W0*np.pi*np.cos(np.pi*z_mid)+self.A_noise*noise
                S_hat_mid['g'] = 1/(-np.pi**2-(self.kx*self.kx+self.ky*self.ky))*self.dy_S_mean/self.tau*W0*np.sin(np.pi*z_mid)+self.A_noise*noise
                d_S_hat_mid['g'] =1/(-np.pi**2-(self.kx*self.kx+self.ky*self.ky))*self.dy_S_mean/self.tau* W0*np.pi*np.cos(np.pi*z_mid)+self.A_noise*noise
                T_0_mid['g'] = self.A_noise*noise
                d_T_0_mid['g'] = self.A_noise*noise
                S_0_mid['g'] = self.A_noise*noise
                d_S_0_mid['g'] = self.A_noise*noise
                
                  
                #This is for the other B.C. like the 
                z_top=(1+h)/2+(1-h)/2/self.Lz*z
                w_hat_top['g'] = W0*np.sin(np.pi*z_top) +self.A_noise*noise
                p_hat_top['g'] = W0*np.pi*np.cos(np.pi*z_top)/(-(self.kx*self.kx+self.ky*self.ky))+self.A_noise*noise
                T_hat_top['g'] = 1/(-np.pi**2-(self.kx*self.kx+self.ky*self.ky))*self.dy_T_mean*W0*np.sin(np.pi*z_top)+self.A_noise*noise
                d_T_hat_top['g'] =1/(-np.pi**2-(self.kx*self.kx+self.ky*self.ky))*self.dy_T_mean* W0*np.pi*np.cos(np.pi*z_top)+self.A_noise*noise
                S_hat_top['g'] = 1/(-np.pi**2-(self.kx*self.kx+self.ky*self.ky))*self.dy_S_mean/self.tau*W0*np.sin(np.pi*z_top)+self.A_noise*noise
                d_S_hat_top['g'] =1/(-np.pi**2-(self.kx*self.kx+self.ky*self.ky))*self.dy_S_mean/self.tau* W0*np.pi*np.cos(np.pi*z_top)+self.A_noise*noise
                T_0_top['g'] = self.A_noise*noise
                d_T_0_top['g'] = self.A_noise*noise
                S_0_top['g'] = self.A_noise*noise
                d_S_0_top['g'] = self.A_noise*noise
                
                
                
            elif self.flow =='HB_porous_shear':
                #initial guess
                z = domain.grid(0)

                #initial guess for the HB_porous, harmonic balance method for double-diffusive convection within porous media
                w_hat_real = solver.state['w_hat_real']
                p_hat_real = solver.state['p_hat_real']
                T_hat_real = solver.state['T_hat_real']
                d_T_hat_real = solver.state['d_T_hat_real']
                S_hat_real = solver.state['S_hat_real']
                d_S_hat_real = solver.state['d_S_hat_real']
                #T_0_real = solver.state['T_0_real']
                #d_T_0_real = solver.state['d_T_0_real']
                #S_0_real = solver.state['S_0_real']
                #d_S_0_real = solver.state['d_S_0_real']
                
                
                W0=self.A_elevator;
                gshape = domain.dist.grid_layout.global_shape(scales=1)
                slices = domain.dist.grid_layout.slices(scales=1)
                rand = np.random.RandomState(seed=23)
                noise = rand.standard_normal(gshape)[slices]
                  
                #This is for the other B.C. like the 
                w_hat_real['g'] = W0*np.sin(np.pi*z) +self.A_noise*noise
                p_hat_real['g'] = W0*np.pi*np.cos(np.pi*z)/(-(self.kx*self.kx+self.ky*self.ky))+self.A_noise*noise
                T_hat_real['g'] = 1/(-np.pi**2-(self.kx*self.kx+self.ky*self.ky))*self.dy_T_mean*W0*np.sin(np.pi*z)+self.A_noise*noise
                d_T_hat_real['g'] =1/(-np.pi**2-(self.kx*self.kx+self.ky*self.ky))*self.dy_T_mean* W0*np.pi*np.cos(np.pi*z)+self.A_noise*noise
                S_hat_real['g'] = 1/(-np.pi**2-(self.kx*self.kx+self.ky*self.ky))*self.dy_S_mean/self.tau*W0*np.sin(np.pi*z)+self.A_noise*noise
                d_S_hat_real['g'] =1/(-np.pi**2-(self.kx*self.kx+self.ky*self.ky))*self.dy_S_mean/self.tau* W0*np.pi*np.cos(np.pi*z)+self.A_noise*noise
                #T_0['g'] = self.A_noise*noise
                #d_T_0['g'] = self.A_noise*noise
                #S_0['g'] = self.A_noise*noise
                #d_S_0['g'] = self.A_noise*noise
                            
            elif self.flow =='HB_benard':
                z = domain.grid(0)

                #initial guess for the HB_porous, harmonic balance method for double-diffusive convection within porous media
                u_tilde = solver.state['u_tilde']
                d_u_tilde = solver.state['d_u_tilde']
                v_tilde = solver.state['v_tilde']
                d_v_tilde = solver.state['d_v_tilde']
                w_hat = solver.state['w_hat']
                p_hat = solver.state['p_hat']
                T_hat = solver.state['T_hat']
                d_T_hat = solver.state['d_T_hat']
                S_hat = solver.state['S_hat']
                d_S_hat = solver.state['d_S_hat']
                T_0 = solver.state['T_0']
                d_T_0 = solver.state['d_T_0']
                S_0 = solver.state['S_0']
                d_S_0 = solver.state['d_S_0']
                
                W0=self.A_elevator;
                gshape = domain.dist.grid_layout.global_shape(scales=1)
                slices = domain.dist.grid_layout.slices(scales=1)
                rand = np.random.RandomState(seed=23)
                noise = rand.standard_normal(gshape)[slices]
                
                if self.dy_T_mean == 1 and self.dy_S_mean ==1:
                    #This is for the salt finger regime. This can work to find the three layer and asymmetric solutions.
                    
                    #Set up the initial guess such that the
                    if self.F_sin==0:
                        #without shear...
                        u_tilde['g'] = self.kx*W0*np.sin(np.pi*z)/((self.kx*self.kx+self.ky*self.ky))+self.A_noise*noise
                        d_u_tilde['g'] = self.kx*np.pi*W0*np.cos(np.pi*z)/((self.kx*self.kx+self.ky*self.ky))+self.A_noise*noise
                        v_tilde['g'] = self.ky*W0*np.sin(np.pi*z)/((self.kx*self.kx+self.ky*self.ky))+self.A_noise*noise
                        d_v_tilde['g'] = self.ky*W0*np.cos(np.pi*z)/((self.kx*self.kx+self.ky*self.ky))+self.A_noise*noise
                        w_hat['g'] = W0*np.sin(np.pi*z) +self.A_noise*noise
                        p_hat['g'] = (-np.pi*np.pi-self.kx*self.kx-self.ky*self.ky)*W0*np.sin(np.pi*z)/((self.kx*self.kx+self.ky*self.ky))+self.A_noise*noise
                        T_hat['g'] = 1/(-np.pi**2-(self.kx*self.kx+self.ky*self.ky))*self.dy_T_mean*W0*np.sin(np.pi*z)+self.A_noise*noise
                        d_T_hat['g'] =1/(-np.pi**2-(self.kx*self.kx+self.ky*self.ky))*self.dy_T_mean* W0*np.pi*np.cos(np.pi*z)+self.A_noise*noise
                        S_hat['g'] = 1/(-np.pi**2-(self.kx*self.kx+self.ky*self.ky))*self.dy_S_mean/self.tau*W0*np.sin(np.pi*z)+self.A_noise*noise
                        d_S_hat['g'] =1/(-np.pi**2-(self.kx*self.kx+self.ky*self.ky))*self.dy_S_mean/self.tau* W0*np.pi*np.cos(np.pi*z)+self.A_noise*noise
                        T_0['g'] = self.A_noise*noise
                        d_T_0['g'] = self.A_noise*noise
                        S_0['g'] = self.A_noise*noise
                        d_S_0['g'] = self.A_noise*noise
                    else:
                        #with shear...
                        u_tilde['g'] = self.kx*W0*np.sin(np.pi*z)/((self.kx*self.kx+self.ky*self.ky))+self.A_noise*noise
                        d_u_tilde['g'] = self.kx*np.pi*W0*np.cos(np.pi*z)/((self.kx*self.kx+self.ky*self.ky))+self.A_noise*noise
                        v_tilde['g'] = self.ky*W0*np.sin(np.pi*z)/((self.kx*self.kx+self.ky*self.ky))+self.A_noise*noise
                        d_v_tilde['g'] = self.ky*W0*np.cos(np.pi*z)/((self.kx*self.kx+self.ky*self.ky))+self.A_noise*noise
                        w_hat['g'] = W0*np.sin(np.pi*z) +self.A_noise*noise
                        p_hat['g'] = (-np.pi*np.pi-self.kx*self.kx-self.ky*self.ky)*W0*np.sin(np.pi*z)/((self.kx*self.kx+self.ky*self.ky))+self.A_noise*noise
                        T_hat['g'] = 1/(-np.pi**2-(self.kx*self.kx+self.ky*self.ky))*self.dy_T_mean*W0*np.sin(np.pi*z)+self.A_noise*noise
                        d_T_hat['g'] =1/(-np.pi**2-(self.kx*self.kx+self.ky*self.ky))*self.dy_T_mean* W0*np.pi*np.cos(np.pi*z)+self.A_noise*noise
                        S_hat['g'] = 1/(-np.pi**2-(self.kx*self.kx+self.ky*self.ky))*self.dy_S_mean/self.tau*W0*np.sin(np.pi*z)+self.A_noise*noise
                        d_S_hat['g'] =1/(-np.pi**2-(self.kx*self.kx+self.ky*self.ky))*self.dy_S_mean/self.tau* W0*np.pi*np.cos(np.pi*z)+self.A_noise*noise
                        T_0['g'] = self.A_noise*noise
                        d_T_0['g'] = self.A_noise*noise
                        S_0['g'] = self.A_noise*noise
                        d_S_0['g'] = self.A_noise*noise
                        
                elif self.dy_T_mean==-1 and self.dy_S_mean==-1:  
                    #This is the initial guess try to find layer like Gough & Toomre (1982) for diffusive regime
                    if self.F_sin==0:
                        #without shear...
                        u_tilde['g'] = self.kx*W0*np.sin(np.pi*z)/((self.kx*self.kx+self.ky*self.ky))+self.A_noise*noise
                        d_u_tilde['g'] = self.kx*np.pi*W0*np.cos(np.pi*z)/((self.kx*self.kx+self.ky*self.ky))+self.A_noise*noise
                        v_tilde['g'] = self.ky*W0*np.sin(np.pi*z)/((self.kx*self.kx+self.ky*self.ky))+self.A_noise*noise
                        d_v_tilde['g'] = self.ky*W0*np.cos(np.pi*z)/((self.kx*self.kx+self.ky*self.ky))+self.A_noise*noise
                        w_hat['g'] = W0*np.sin(np.pi*z) +self.A_noise*noise
                        p_hat['g'] = (-np.pi*np.pi-self.kx*self.kx-self.ky*self.ky)*W0*np.sin(np.pi*z)/((self.kx*self.kx+self.ky*self.ky))+self.A_noise*noise
                        T_hat['g'] = 1/(-np.pi**2-(self.kx*self.kx+self.ky*self.ky))*self.dy_T_mean*W0*np.sin(np.pi*z)+self.A_noise*noise
                        d_T_hat['g'] =1/(-np.pi**2-(self.kx*self.kx+self.ky*self.ky))*self.dy_T_mean* W0*np.pi*np.cos(np.pi*z)+self.A_noise*noise
                        S_hat['g'] = 1/(-np.pi**2-(self.kx*self.kx+self.ky*self.ky))*self.dy_S_mean/self.tau*W0*np.sin(np.pi*z)+self.A_noise*noise
                        d_S_hat['g'] =1/(-np.pi**2-(self.kx*self.kx+self.ky*self.ky))*self.dy_S_mean/self.tau* W0*np.pi*np.cos(np.pi*z)+self.A_noise*noise
                        
                        T_0['g'] = W0*np.sin(self.initial_kz*z)+ self.A_noise*noise
                        d_T_0['g'] = W0*self.initial_kz*np.cos(self.initial_kz*z)+self.A_noise*noise
                        S_0['g'] = W0*np.sin(self.initial_kz*z)+self.A_noise*noise
                        d_S_0['g'] = W0*self.initial_kz*np.cos(self.initial_kz*z)+self.A_noise*noise
                    
                    
                    
            elif self.flow =='test_periodic':
                z = domain.grid(0)
                w_hat =solver.state['w_hat']
                p_hat = solver.state['p_hat']
                T_hat = solver.state['T_hat']
                d_T_hat = solver.state['d_T_hat']
                w_hat['g']=self.F_sin*np.sin(self.ks*z)
                p_hat['g']=self.F_sin*np.cos(self.ks*z)
                T_hat['g']=self.F_sin*np.sin(self.ks*z)
                d_T_hat['g']=self.F_sin*np.cos(self.ks*z)
            elif self.flow=='HB_benard_shear':
                z = domain.grid(0)

                #initial guess for the HB_porous, harmonic balance method for double-diffusive convection within porous media
                u_tilde_real = solver.state['u_tilde_real']
                d_u_tilde_real = solver.state['d_u_tilde_real']
                v_tilde_real = solver.state['v_tilde_real']
                d_v_tilde_real = solver.state['d_v_tilde_real']
                w_hat_real = solver.state['w_hat_real']
                p_hat_real = solver.state['p_hat_real']
                T_hat_real = solver.state['T_hat_real']
                d_T_hat_real = solver.state['d_T_hat_real']
                S_hat_real = solver.state['S_hat_real']
                d_S_hat_real = solver.state['d_S_hat_real']
                
                u_tilde_imag = solver.state['u_tilde_imag']
                d_u_tilde_imag = solver.state['d_u_tilde_imag']
                v_tilde_imag = solver.state['v_tilde_imag']
                d_v_tilde_imag = solver.state['d_v_tilde_imag']
                w_hat_imag = solver.state['w_hat_imag']
                p_hat_imag = solver.state['p_hat_imag']
                T_hat_imag = solver.state['T_hat_imag']
                d_T_hat_imag = solver.state['d_T_hat_imag']
                S_hat_imag = solver.state['S_hat_imag']
                d_S_hat_imag = solver.state['d_S_hat_imag']
                
                T_0 = solver.state['T_0']
                d_T_0 = solver.state['d_T_0']
                S_0 = solver.state['S_0']
                d_S_0 = solver.state['d_S_0']
                
                W0=self.A_elevator;
                gshape = domain.dist.grid_layout.global_shape(scales=1)
                slices = domain.dist.grid_layout.slices(scales=1)
                rand = np.random.RandomState(seed=23)
                noise = rand.standard_normal(gshape)[slices]
                
                # #with shear...
                # T_hat_elevator=1/(-(self.kx*self.kx+self.ky*self.ky))*self.dy_T_mean*W0
                # S_hat_elevator=1/(-(self.kx*self.kx+self.ky*self.ky))*self.dy_S_mean/self.tau*W0
                # p_hat_elevator=0#-(self.kx*self.kx+self.ky*self.ky)*W0+self.Ra_T*T_hat_elevator-self.Ra_S2T*S_hat_elevator
                # u_tilde_elevator=0#-self.kx/(self.kx*self.kx+self.ky*self.ky)*p_hat_elevator
                # v_tilde_elevator=0#-self.ky/(self.kx*self.kx+self.ky*self.ky)*p_hat_elevator
                
                # u_tilde_real['g'] = u_tilde_elevator+self.A_noise*noise
                # #d_u_tilde_real['g'] = self.kx*np.pi*W0*np.cos(self.ks*z)/((self.kx*self.kx+self.ky*self.ky))+self.A_noise*noise
                # v_tilde_real['g'] = v_tilde_elevator+self.A_noise*noise
                # #d_v_tilde_real['g'] = self.ky*W0*np.cos(self.ks*z)/((self.kx*self.kx+self.ky*self.ky))+self.A_noise*noise
                # w_hat_real['g'] = W0 +self.A_noise*noise
                # p_hat_real['g'] = p_hat_elevator+self.A_noise*noise
                # T_hat_real['g'] = T_hat_elevator+self.A_noise*noise
                # #d_T_hat_real['g'] =1/(-self.ks*self.ks-(self.kx*self.kx+self.ky*self.ky))*self.dy_T_mean* W0*np.pi*np.cos(self.ks*z)+self.A_noise*noise
                # S_hat_real['g'] = S_hat_elevator+self.A_noise*noise
                # #d_S_hat_real['g'] =1/(-self.ks*self.ks-(self.kx*self.kx+self.ky*self.ky))*self.dy_S_mean/self.tau* W0*np.pi*np.cos(self.ks*z)+self.A_noise*noise
                
                #without shear...
                u_tilde_real['g'] = self.kx*W0*np.sin(np.pi*z)/((self.kx*self.kx+self.ky*self.ky))+self.A_noise*noise
                d_u_tilde_real['g'] = self.kx*np.pi*W0*np.cos(np.pi*z)/((self.kx*self.kx+self.ky*self.ky))+self.A_noise*noise
                v_tilde_real['g'] = self.ky*W0*np.sin(np.pi*z)/((self.kx*self.kx+self.ky*self.ky))+self.A_noise*noise
                d_v_tilde_real['g'] = self.ky*W0*np.cos(np.pi*z)/((self.kx*self.kx+self.ky*self.ky))+self.A_noise*noise
                w_hat_real['g'] = W0*np.sin(np.pi*z) +self.A_noise*noise
                p_hat_real['g'] = (-np.pi*np.pi-self.kx*self.kx-self.ky*self.ky)*W0*np.sin(np.pi*z)/((self.kx*self.kx+self.ky*self.ky))+self.A_noise*noise
                T_hat_real['g'] = 1/(-np.pi**2-(self.kx*self.kx+self.ky*self.ky))*self.dy_T_mean*W0*np.sin(np.pi*z)+self.A_noise*noise
                d_T_hat_real['g'] =1/(-np.pi**2-(self.kx*self.kx+self.ky*self.ky))*self.dy_T_mean* W0*np.pi*np.cos(np.pi*z)+self.A_noise*noise
                S_hat_real['g'] = 1/(-np.pi**2-(self.kx*self.kx+self.ky*self.ky))*self.dy_S_mean/self.tau*W0*np.sin(np.pi*z)+self.A_noise*noise
                d_S_hat_real['g'] =1/(-np.pi**2-(self.kx*self.kx+self.ky*self.ky))*self.dy_S_mean/self.tau* W0*np.pi*np.cos(np.pi*z)+self.A_noise*noise
                
                W0_imag=self.A_elevator_imag
                u_tilde_imag['g'] = self.kx*W0_imag*np.sin(np.pi*z)/((self.kx*self.kx+self.ky*self.ky))+self.A_noise*noise
                d_u_tilde_imag['g'] = self.kx*np.pi*W0_imag*np.cos(np.pi*z)/((self.kx*self.kx+self.ky*self.ky))+self.A_noise*noise
                v_tilde_imag['g'] = self.ky*W0_imag*np.sin(np.pi*z)/((self.kx*self.kx+self.ky*self.ky))+self.A_noise*noise
                d_v_tilde_imag['g'] = self.ky*W0_imag*np.cos(np.pi*z)/((self.kx*self.kx+self.ky*self.ky))+self.A_noise*noise
                w_hat_imag['g'] = W0_imag*np.sin(np.pi*z) +self.A_noise*noise
                p_hat_imag['g'] = (-np.pi*np.pi-self.kx*self.kx-self.ky*self.ky)*W0_imag*np.sin(np.pi*z)/((self.kx*self.kx+self.ky*self.ky))+self.A_noise*noise
                T_hat_imag['g'] = 1/(-np.pi**2-(self.kx*self.kx+self.ky*self.ky))*self.dy_T_mean*W0_imag*np.sin(np.pi*z)+self.A_noise*noise
                d_T_hat_imag['g'] =1/(-np.pi**2-(self.kx*self.kx+self.ky*self.ky))*self.dy_T_mean* W0_imag*np.pi*np.cos(np.pi*z)+self.A_noise*noise
                S_hat_imag['g'] = 1/(-np.pi**2-(self.kx*self.kx+self.ky*self.ky))*self.dy_S_mean/self.tau*W0_imag*np.sin(np.pi*z)+self.A_noise*noise
                d_S_hat_imag['g'] =1/(-np.pi**2-(self.kx*self.kx+self.ky*self.ky))*self.dy_S_mean/self.tau* W0_imag*np.pi*np.cos(np.pi*z)+self.A_noise*noise
                
                
                # T_0['g'] = self.A_noise*noise
                # d_T_0['g'] = self.A_noise*noise
                # S_0['g'] = self.A_noise*noise
                # d_S_0['g'] = self.A_noise*noise
           
            
                # u_tilde_imag['g'] = self.kx*W0*np.sin(self.ks*z)/((self.kx*self.kx+self.ky*self.ky))+self.A_noise*noise
                # d_u_tilde_imag['g'] = self.kx*self.ks*W0*np.cos(self.ks*z)/((self.kx*self.kx+self.ky*self.ky))+self.A_noise*noise
                # v_tilde_imag['g'] = self.ky*W0*np.sin(self.ks*z)/((self.kx*self.kx+self.ky*self.ky))+self.A_noise*noise
                # d_v_tilde_imag['g'] = self.ky*W0*np.cos(self.ks*z)/((self.kx*self.kx+self.ky*self.ky))+self.A_noise*noise
                # w_hat_imag['g'] = W0*np.sin(self.ks*z) +self.A_noise*noise
                # p_hat_imag['g'] = (-self.ks*self.ks-self.kx*self.kx-self.ky*self.ky)*W0*np.sin(self.ks*z)/((self.kx*self.kx+self.ky*self.ky))+self.A_noise*noise
                # T_hat_imag['g'] = 1/(-self.ks*self.ks-(self.kx*self.kx+self.ky*self.ky))*self.dy_T_mean*W0*np.sin(self.ks*z)+self.A_noise*noise
                # d_T_hat_imag['g'] =1/(-self.ks*self.ks-(self.kx*self.kx+self.ky*self.ky))*self.dy_T_mean* W0*np.pi*np.cos(self.ks*z)+self.A_noise*noise
                # S_hat_imag['g'] = 1/(-self.ks*self.ks-(self.kx*self.kx+self.ky*self.ky))*self.dy_S_mean/self.tau*W0*np.sin(self.ks*z)+self.A_noise*noise
                # d_S_hat_imag['g'] =1/(-self.ks*self.ks-(self.kx*self.kx+self.ky*self.ky))*self.dy_S_mean/self.tau* W0*np.pi*np.cos(self.ks*z)+self.A_noise*noise
                
                
                T_0['g'] =-np.sin(2*np.pi*z)+ self.A_noise*noise
                d_T_0['g'] =-2*np.pi*np.cos(2*np.pi*z)+ self.A_noise*noise
                S_0['g'] =-np.sin(2*np.pi*z)+ self.A_noise*noise
                d_S_0['g'] =-2*np.pi*np.cos(2*np.pi*z)+ self.A_noise*noise
            
                
        else:
            #Restart
            print('restart')
            #self.EVP_trivial=0
            write, last_dt = solver.load_state('restart.h5', -1)
            if self.A_noise !=0 and self.flow in ['HB_benard']:
                gshape = domain.dist.grid_layout.global_shape(scales=1.5)
                slices = domain.dist.grid_layout.slices(scales=1.5)
                rand = np.random.RandomState(seed=23)
                noise = rand.standard_normal(gshape)[slices]
                
                #print(gshape)
                #print(slices)
                z = domain.grid(0)

                #initial guess for the HB_porous, harmonic balance method for double-diffusive convection within porous media
                u_tilde = solver.state['u_tilde']
                d_u_tilde = solver.state['d_u_tilde']
                v_tilde = solver.state['v_tilde']
                d_v_tilde = solver.state['d_v_tilde']
                w_hat = solver.state['w_hat']
                p_hat = solver.state['p_hat']
                T_hat = solver.state['T_hat']
                d_T_hat = solver.state['d_T_hat']
                S_hat = solver.state['S_hat']
                d_S_hat = solver.state['d_S_hat']
                T_0 = solver.state['T_0']
                d_T_0 = solver.state['d_T_0']
                S_0 = solver.state['S_0']
                d_S_0 = solver.state['d_S_0']
                
                #print(np.size(u_tilde['g']))
                #print(np.size(noise))
                u_tilde['g'] = u_tilde['g']+self.A_noise*noise
                #d_u_tilde['g'] = d_u_tilde['g']+self.A_noise*noise
                v_tilde['g'] = v_tilde['g']+self.A_noise*noise
                #d_v_tilde['g'] = d_v_tilde['g']+self.A_noise*noise
                w_hat['g'] = w_hat['g']+self.A_noise*noise
                p_hat['g'] = p_hat['g']+self.A_noise*noise
                T_hat['g'] = T_hat['g']+self.A_noise*noise
                #d_T_hat['g'] = d_T_hat['g']+self.A_noise*noise
                S_hat['g'] = S_hat['g']+self.A_noise*noise
                #d_S_hat['g'] = d_S_hat['g']+self.A_noise*noise
                T_0['g'] = T_0['g']+self.A_noise*noise
                #d_T_0['g'] = d_T_0['g']+self.A_noise*noise
                S_0['g'] = S_0['g']+self.A_noise*noise
                #d_S_0['g'] = d_S_0['g']+self.A_noise*noise
            elif self.A_noise!=0 and self.flow in ['HB_benard_shear']:
                gshape = domain.dist.grid_layout.global_shape(scales=1.5)
                slices = domain.dist.grid_layout.slices(scales=1.5)
                rand = np.random.RandomState(seed=23)
                noise = rand.standard_normal(gshape)[slices]
                
                #print(gshape)
                #print(slices)
                z = domain.grid(0)

                #initial guess for the HB_porous, harmonic balance method for double-diffusive convection within porous media
                u_tilde_real = solver.state['u_tilde_real']
                u_tilde_imag = solver.state['u_tilde_imag']
                #d_u_tilde = solver.state['d_u_tilde']
                v_tilde_real = solver.state['v_tilde_real']
                v_tilde_imag = solver.state['v_tilde_imag']
                #d_v_tilde = solver.state['d_v_tilde']
                w_hat_real = solver.state['w_hat_real']
                w_hat_imag = solver.state['w_hat_imag']
                p_hat_real = solver.state['p_hat_real']
                p_hat_imag = solver.state['p_hat_imag']
                T_hat_real = solver.state['T_hat_real']
                T_hat_imag = solver.state['T_hat_imag']
                #d_T_hat = solver.state['d_T_hat']
                S_hat_real = solver.state['S_hat_real']
                S_hat_imag = solver.state['S_hat_imag']
                #d_S_hat = solver.state['d_S_hat']
                T_0 = solver.state['T_0']
                #d_T_0 = solver.state['d_T_0']
                S_0 = solver.state['S_0']
                #d_S_0 = solver.state['d_S_0']
                U_0 = solver.state['U_0']
                
                u_tilde_real['g'] = u_tilde_real['g']+self.A_noise*noise
                u_tilde_imag['g'] = u_tilde_imag['g']+self.A_noise*noise
                v_tilde_real['g'] = v_tilde_real['g']+self.A_noise*noise
                v_tilde_imag['g'] = v_tilde_imag['g']+self.A_noise*noise
                w_hat_real['g'] = w_hat_real['g']+self.A_noise*noise
                w_hat_imag['g'] = w_hat_imag['g']+self.A_noise*noise
                p_hat_real['g'] = p_hat_real['g']+self.A_noise*noise
                p_hat_imag['g'] = p_hat_imag['g']+self.A_noise*noise
                T_hat_real['g'] = T_hat_real['g']+self.A_noise*noise
                T_hat_imag['g'] = T_hat_imag['g']+self.A_noise*noise
                S_hat_real['g'] = S_hat_real['g']+self.A_noise*noise
                S_hat_imag['g'] = S_hat_imag['g']+self.A_noise*noise

                T_0['g'] = T_0['g']+self.A_noise*noise
                S_0['g'] = S_0['g']+self.A_noise*noise
                U_0['g'] = U_0['g']+self.A_noise*noise
              
            if self.A_noise !=0 and self.flow in ['double_diffusive_shear_2D']:
                x = domain.grid(0)
                z = domain.grid(1)
                u = solver.state['u']
                w = solver.state['w']
                S = solver.state['S']
                p = solver.state['p']
                T = solver.state['T']
                
                #Update 2022/05/05, we should use the same dealizing scales for all direction. the scales here should also be the same, thus 1.5
                gshape = domain.dist.grid_layout.global_shape(scales=1.5)
                slices = domain.dist.grid_layout.slices(scales=1.5)
                rand = np.random.RandomState(seed=23)
                noise = rand.standard_normal(gshape)[slices]
                
                ##Add the random noise
                u['g']=u['g']+self.A_noise*noise
                w['g']=w['g']+self.A_noise*noise
                S['g']=S['g']+self.A_noise*noise
                T['g']=S['g']+self.A_noise*noise
                p['g']=p['g']+self.A_noise*noise
                  
              
            if self.continuation_asymmetric ==1 and self.flow =='HB_benard':
                z = domain.grid(0)

                #initial guess for the HB_porous, harmonic balance method for double-diffusive convection within porous media
                u_tilde = solver.state['u_tilde']
                d_u_tilde = solver.state['d_u_tilde']
                v_tilde = solver.state['v_tilde']
                d_v_tilde = solver.state['d_v_tilde']
                w_hat = solver.state['w_hat']
                p_hat = solver.state['p_hat']
                T_hat = solver.state['T_hat']
                d_T_hat = solver.state['d_T_hat']
                S_hat = solver.state['S_hat']
                d_S_hat = solver.state['d_S_hat']
                T_0 = solver.state['T_0']
                d_T_0 = solver.state['d_T_0']
                S_0 = solver.state['S_0']
                d_S_0 = solver.state['d_S_0']
                
                u_tilde['g'] = np.flip(u_tilde['g'])
                d_u_tilde['g'] = -np.flip(d_u_tilde['g'])
                v_tilde['g'] = np.flip(v_tilde['g'])
                d_v_tilde['g'] = -np.flip(d_v_tilde['g'])
                w_hat['g'] = -np.flip(w_hat['g'])
                p_hat['g'] = np.flip(p_hat['g'])
                T_hat['g'] = -np.flip(T_hat['g'])
                d_T_hat['g'] = np.flip(d_T_hat['g'])
                S_hat['g'] = -np.flip(S_hat['g'])
                d_S_hat['g'] = np.flip(d_S_hat['g'])
                T_0['g'] = -np.flip(T_0['g'])
                d_T_0['g'] = np.flip(d_T_0['g'])
                S_0['g'] = -np.flip(S_0['g'])
                d_S_0['g'] = np.flip(d_S_0['g'])
              
        #If set the continuation... then just load the existing data...
        if self.continuation != 0:
            #firstly make a copy of the old data
            shutil.copytree('analysis','analysis'+str(self.continuation))
            
            #Then load_state as the initial condition of IVP/initial guess of BVP
            write, last_dt = solver.load_state('./analysis/analysis_s1.h5', -1)
              
                    
        
    def run(self,solver,domain,logger):
        ##This CFL condition need to be modified for different simulation configuration.
        ##Note this line is very important and it needs to be added!!!!!
        if self.problem =='IVP':
            solver.stop_sim_time = self.stop_sim_time
    
            logger.info('Starting loop')
            solver.stop_wall_time = np.inf
            solver.stop_iteration = np.inf
            start_time=time.time()
            dt = self.initial_dt
            
            if self.flow in ['IFSC_2D','double_diffusive_2D','double_diffusive_shear_2D','porous_media_2D']:
                #cfl is required for the IVP
                cfl = flow_tools.CFL(solver,self.initial_dt,safety=0.8,max_change=1,cadence=8)
                cfl.add_velocities(('u','w'))
            #elif self.flow =='HB_porous':
            #    cfl.add_velocities('w_hat')
            #elif self.flow =='HB_benard':
            #    cfl.add_velocities(('w_hat','u_tilde','v_tilde'))
                while solver.ok:
                    dt = cfl.compute_dt()    
                    solver.step(dt)
                    if solver.iteration % 100 == 0:
                        logger.info('Iteration: %i, Time: %e, dt: %e' %(solver.iteration, solver.sim_time, dt))
                        if self.flux_T:
                            dy_T_mean_q=-(1-np.sum(np.sum(solver.state['w']['g']*solver.state['T']['g']))/self.Lx/self.Lz)
                            logger.info('dy_T_mean: {}'.format(dy_T_mean_q))
                            logger.info('Nu: {}'.format(-1/dy_T_mean_q))

            elif self.flow in ['HB_porous','HB_benard','HB_porous_shear','HB_benard_shear']:
                #This harmonic balance is 1D simulation and thus no CFL condition is required.... 
                while solver.ok:
                    solver.step(dt)
                    if solver.iteration % 100 == 0:
                        logger.info('Iteration: %i, Time: %e, dt: %e' %(solver.iteration, solver.sim_time, dt))
                        logger.info('Nu: {}'.format(-(solver.state['d_T_0']['g'][0]+self.dy_T_mean)))
                        logger.info('Nu_S: {}'.format(-(solver.state['d_S_0']['g'][0]+self.dy_S_mean)))
                
            
            end_time = time.time()
            # Print statistics
            logger.info('Run time: %f' %(end_time-start_time))
            logger.info('Iterations: %i' %solver.iteration)
            logger.info('Run time: %f cpu-hr' %((end_time-start_time)/60/60*domain.dist.comm_cart.size))
        elif self.problem=='BVP':
            # Newton iterations for BVP
            pert = solver.perturbations.data
            pert.fill(1+self.bvp_tolerance)
            start_time = time.time()
            while np.sum(np.abs(pert)) > self.bvp_tolerance:
                solver.newton_iteration()
                logger.info('Perturbation norm: {}'.format(np.sum(np.abs(pert))))
                #Update 2022/03/01: also print the Nusselt number for T and S on the screen. 
                logger.info('Nu: {}'.format(-(solver.state['d_T_0']['g'][0]+self.dy_S_mean)))
                logger.info('Nu_S: {}'.format(-(solver.state['d_S_0']['g'][0]+self.dy_T_mean)))
                logger.info('T_0 norm: {}'.format(np.sum(np.abs(solver.state['T_0']['g']))))
                logger.info('S_0 norm: {}'.format(np.sum(np.abs(solver.state['S_0']['g']))))
                if self.flow in ['HB_porous','HB_benard']:
                    logger.info('w_hat norm: {}'.format(np.sum(np.abs(solver.state['w_hat']['g']))))
                    logger.info('T_hat norm: {}'.format(np.sum(np.abs(solver.state['T_hat']['g']))))
                    logger.info('S_hat norm: {}'.format(np.sum(np.abs(solver.state['S_hat']['g']))))
                elif self.flow in ['HB_porous_2_layer']:
                    logger.info('w_hat norm: {}'.format(np.sum(np.abs(solver.state['w_hat']['g']))))
                    logger.info('T_hat norm: {}'.format(np.sum(np.abs(solver.state['T_hat']['g']))))
                    logger.info('S_hat norm: {}'.format(np.sum(np.abs(solver.state['S_hat']['g']))))
                    logger.info('T_0_top norm: {}'.format(np.sum(np.abs(solver.state['T_0_top']['g']))))
                    logger.info('S_0_top norm: {}'.format(np.sum(np.abs(solver.state['S_0_top']['g']))))
                    logger.info('w_hat_top norm: {}'.format(np.sum(np.abs(solver.state['w_hat_top']['g']))))
                    logger.info('T_hat_top norm: {}'.format(np.sum(np.abs(solver.state['T_hat_top']['g']))))
                    logger.info('S_hat_top norm: {}'.format(np.sum(np.abs(solver.state['S_hat_top']['g']))))
                elif self.flow in ['HB_porous_3_layer']:
                    logger.info('w_hat norm: {}'.format(np.sum(np.abs(solver.state['w_hat']['g']))))
                    logger.info('T_hat norm: {}'.format(np.sum(np.abs(solver.state['T_hat']['g']))))
                    logger.info('S_hat norm: {}'.format(np.sum(np.abs(solver.state['S_hat']['g']))))
                    
                    logger.info('T_0_mid norm: {}'.format(np.sum(np.abs(solver.state['T_0_mid']['g']))))
                    logger.info('S_0_mid norm: {}'.format(np.sum(np.abs(solver.state['S_0_mid']['g']))))
                    logger.info('w_hat_mid norm: {}'.format(np.sum(np.abs(solver.state['w_hat_mid']['g']))))
                    logger.info('T_hat_mid norm: {}'.format(np.sum(np.abs(solver.state['T_hat_mid']['g']))))
                    logger.info('S_hat_mid norm: {}'.format(np.sum(np.abs(solver.state['S_hat_mid']['g']))))
                
                    logger.info('T_0_top norm: {}'.format(np.sum(np.abs(solver.state['T_0_top']['g']))))
                    logger.info('S_0_top norm: {}'.format(np.sum(np.abs(solver.state['S_0_top']['g']))))
                    logger.info('w_hat_top norm: {}'.format(np.sum(np.abs(solver.state['w_hat_top']['g']))))
                    logger.info('T_hat_top norm: {}'.format(np.sum(np.abs(solver.state['T_hat_top']['g']))))
                    logger.info('S_hat_top norm: {}'.format(np.sum(np.abs(solver.state['S_hat_top']['g']))))
                    
                elif self.flow in ['HB_benard_shear','HB_porous_shear']:
                    logger.info('w_hat_real norm: {}'.format(np.sum(np.abs(solver.state['w_hat_real']['g']))))
                    logger.info('w_hat_imag norm: {}'.format(np.sum(np.abs(solver.state['w_hat_imag']['g']))))
                    logger.info('T_hat_real norm: {}'.format(np.sum(np.abs(solver.state['T_hat_real']['g']))))
                    logger.info('T_hat_imag norm: {}'.format(np.sum(np.abs(solver.state['T_hat_imag']['g']))))
                    logger.info('S_hat_real norm: {}'.format(np.sum(np.abs(solver.state['S_hat_real']['g']))))
                    logger.info('S_hat_imag norm: {}'.format(np.sum(np.abs(solver.state['S_hat_imag']['g']))))
                    logger.info('eta: {}'.format(solver.state['eta']['g']))
                #logger.info('R iterate: {}'.format(R['g'][0]))
            end_time = time.time()
        elif self.problem=='EVP':
            logger.info('Computing max growth rate')
            solver.solve_dense(solver.pencils[0])
            
            # Filter infinite/nan eigenmodes
            finite = np.isfinite(solver.eigenvalues)
            #print(solver.eigenvectors.shape)
            solver.eigenvalues = solver.eigenvalues[finite]
            solver.eigenvectors = solver.eigenvectors[:, finite]
            #print(solver.eigenvectors.shape)   
            # Sort eigenmodes by eigenvalue
            # Update 2022/02/28, order is based on the real part of the eigenvalues and also sort from the largest one to the smallest one.
            # which also require a flip to flip the order. 
            order = np.flip(np.argsort(np.real(solver.eigenvalues)))
            solver.eigenvalues = solver.eigenvalues[order]
            #print(solver.eigenvectors.shape)
            #logger.info(order)
            solver.eigenvectors = solver.eigenvectors[:, order]
            
            #Update 2022/03/01, only store the eigenvector corresponding to the largest real part to save space. Also note that python start from 0 index
            solver.eigenvectors = solver.eigenvectors[:,0]
            logger.info('The maximum real part of the all eigenvalues are \n')
            logger.info(np.max(np.real(solver.eigenvalues)))
            logger.info('Corresponding complex eigenvalue is:\n')
            logger.info((solver.eigenvalues[0]))
            #logger.info(solver.eigenvalues)
            
            #store the eigenvalue and eigenvectors into a new field
            #Update 2022/02/28
            #eigenvalues=domain.new_field()
            #eigenvectors=domain.new_field()
            #eigenvalues['g']=solver.eigenvalues
            #eigenvectors['g']=solver.eigenvectors
            
    def post_store(self,solver):
        #This post-processing variable need to be modified for different flow configuration
        if self.flow in ['IFSC_2D','double_diffusive_2D','double_diffusive_shear_2D']:
            analysis = solver.evaluator.add_file_handler('analysis',sim_dt=self.post_store_dt)
            #max_write=400... 
            #Update 2022/02/26, add the system 
            
            if self.store_variable =='all':
                analysis.add_system(solver.state)
    
                #analysis.add_task('S',layout='g',name='S')
                #analysis.add_task('T',layout='g',name='T')
                #analysis.add_task('u',layout='g',name='u')
                #analysis.add_task('w',layout='g',name='w')
                #analysis.add_task('p',layout='g',name='p')
                
                analysis.add_task("S",layout='c',name='S_coeff')
                analysis.add_task("T",layout='c',name='T_coeff')
                analysis.add_task("u",layout='c',name='u_coeff')
                analysis.add_task("w",layout='c',name='w_coeff')
                analysis.add_task("p",layout='c',name='p_coeff')
                
            elif self.store_variable =='S_u_w':
                analysis.add_task('S',layout='g',name='S')
                analysis.add_task('u',layout='g',name='u')
                analysis.add_task('w',layout='g',name='w')
            elif self.store_variable =='S_T_u_w':
                analysis.add_task('S',layout='g',name='S')
                analysis.add_task('T',layout='g',name='T')
                analysis.add_task('u',layout='g',name='u')
                analysis.add_task('w',layout='g',name='w')
            elif self.store_variable =='T_u_w':
                analysis.add_task('T',layout='g',name='T')
                analysis.add_task('u',layout='g',name='u')
                analysis.add_task('w',layout='g',name='w')
              
            if (self.flux_T):
                analysis.add_task('-(1-integ(w*T)/Lx/Lz)',layout='g',name='dy_T_mean_q')
                
        elif self.flow in ['porous_media_2D']:
            analysis = solver.evaluator.add_file_handler('analysis',sim_dt=self.post_store_dt)
            analysis.add_system(solver.state)

            #analysis.add_task('T',layout='g',name='T')
            #analysis.add_task('u',layout='g',name='u')
            #analysis.add_task('w',layout='g',name='w')
            #analysis.add_task('p',layout='g',name='p')
            
            analysis.add_task("T",layout='c',name='T_coeff')
            analysis.add_task("u",layout='c',name='u_coeff')
            analysis.add_task("w",layout='c',name='w_coeff')
            analysis.add_task("p",layout='c',name='p_coeff')
        elif self.flow in ['HB_porous','HB_porous_shear','HB_benard','test_periodic','HB_benard_shear','HB_porous_2_layer','HB_porous_3_layer']:
            #For IVP and BVP, they have some small difference. IVP can also set the dt for storage.
            if self.problem == 'IVP':
                analysis = solver.evaluator.add_file_handler('analysis',sim_dt=self.post_store_dt)
                analysis.add_system(solver.state)
            elif self.problem in ['BVP']:
                self.analysis = solver.evaluator.add_file_handler('analysis')
                self.analysis.add_system(solver.state)
                #For BVP problem, here need to do nothing
            #elif self.problem in ['EVP']:
            #    self.analysis = solver.evaluator.add_file_handler('analysis')
            #    self.analysis.add_system(solver.state)

    def post_store_after_run(self,solver):
        #merge step for the IVP and BVP...
        if self.problem == 'IVP':
            post.merge_process_files('analysis',cleanup=True)
        elif self.problem in ['BVP']:
            #Here is the place to output the post-processing of BVP
            print('post store after run')
            solver.evaluator.evaluate_handlers([self.analysis], world_time=0, wall_time=0, sim_time=0, timestep=0, iteration=0)
            post.merge_process_files('analysis',cleanup=True)
            
            if self.EVP_secondary:
                
                #if I would like to compute the stability of the secondary state, 
                #Then, continuation number +1, and linearize around this state to redo the eigenvalue problem.
                #self.continuation=self.continuation+1
                self_EVP_secondary=self
                self_EVP_secondary.problem='EVP'
                domain_EVP_secondary=self_EVP_secondary.build_domain()
                solver_EVP_secondary=self_EVP_secondary.governing_equation(domain_EVP_secondary,solver)
                logger = logging.getLogger(__name__)
                #self_EVP_secondary.print_screen(logger)
                #shutil.copytree('analysis','analysis'+str(self.continuation))
                #self_EVP_secondary.initial_condition(domain_EVP_secondary,solver)
                #self_EVP_secondary.post_store(solver_EVP_secondary)
                #self_EVP_secondary.print_file() #move print file to here.
                self_EVP_secondary.run(solver_EVP_secondary,domain_EVP_secondary,logger)
                self_EVP_secondary.post_store_after_run(solver_EVP_secondary)
                
                #print('self.problem='+self.problem)
                #print('self_EVP_secondary.problem='+self_EVP_secondary.problem)
                #self.problem='BVP'
                #Update 2022/02/28, note that this local version of self works more like a pointer...
                #self_EVP_secondary.problem='EVP' will also change self.problem='EVP', so here, I need to change it back
                #Otherwise, the next continuation of the BVP will not work...
                self_EVP_secondary.problem='BVP'
                
                
        #write the eigenvalues and eigenvectors into the code...
        elif self.problem =='EVP':
            with h5py.File('./analysis/analysis_s1.h5', 'a') as f:
                eigenvalues = f.create_dataset("eigenvalues", data=solver.eigenvalues)
                eigenvectors = f.create_dataset("eigenvectors", data=solver.eigenvectors)
                print('Append eigenvalues and eigenvectors into ./analysis/analysis_s1.h5')
             
    def get_HB_porous_2_layer_Omega_k(self):
        ##This is the data from figure 10(b) at Ra=5000 from Hewitt DR, Neufeld JA, Lister JR. High Rayleigh number convection in a porous medium containing a thin low-permeability layer. Journal of fluid mechanics. 2014 Oct;756:844-69.
        Omega_list=[0.00249159,0.00498315,0.00631384,0.00996625,0.0126276,0.0199324,0.0252551, \
               0.0299069,0.0398646,0.0496633,0.0598134,0.0797288,0.099326,0.119626,0.159457, \
                0.198651, 0.239251, 0.3973, 0.637821, 0.808144, 0.794597, 1.27564, 1.58919, 1.61628, \
                2.55126, 3.17836, 5.01696, 10.0339, 19.7312, 50.0002,100]
        kx_list_2pi=[2.42908,2.05674,2.5, 2.30496, 2.07447, 2.28723, 2.16312, 2.21631, 2.2695, 2.02128, \
                    2.00355, 1.66667, 2.00355, 1.75532,1.29433, 1.25887, 1.25887, 0.762411, 1.01064,1.01064, \
                        0.744681, 0.744681, 0.744681, 0.514184, 0.496454, 0.496454, 2.42908, 2.76596, \
                        2.9078, 2.69504, 2.92553]
        kx_list=[kx*2*np.pi for kx in kx_list_2pi]
        return Omega_list, kx_list