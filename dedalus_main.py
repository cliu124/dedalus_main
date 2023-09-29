import numpy as np
from mpi4py import MPI
import pathlib
#from IPython import display
from dedalus import public as de
from dedalus.extras import flow_tools
import logging

##These lib are used for post-processing
import h5py
import pathlib
import subprocess
import shutil
import dedalus_setup
shutil.rmtree('analysis',ignore_errors=True)

logger = logging.getLogger(__name__)
flag=dedalus_setup.flag()



#------------select the flow configuration and special parameters for each
#flag.flow='HB_porous_3_layer'
#flag.flow='HB_porous'
#flag.flow='HB_benard_shear'
#flag.flow='test_periodic'

#This is runing 2D DNS general formulation
flag.flow='double_diffusive_shear_2D'#['IFSC_2D','double_diffusive_2D','double_diffusive_shear_2D','porous_media_2D']
#flag.flow='HB_benard_shear_periodic'
#flag.flow='porous_media_2D'
#flag.flow_sub_double_diffusive_shear_2D='primitive_dirichlet_salt_finger'
#flag.flow_sub_double_diffusive_shear_2D='primitive_stress_free_salt_finger'
#flag.flow_sub_double_diffusive_shear_2D='primitive_periodic_salt_finger'
flag.flow_sub_double_diffusive_shear_2D='primitive_periodic_RBC'


#flag.flow_sub_double_diffusive_shear_2D='double_diffusive'
#flag.flow_sub_double_diffusive_shear_2D='IFSC'
#flag.flow_sub_double_diffusive_shear_2D='MRBC'
#flag.flow_sub_double_diffusive_shear_2D='Stokes'
#flag.flow_sub_double_diffusive_shear_2D='primitive_IFSC_unit_tuS'
#flag.flow_sub_double_diffusive_shear_2D='shear_Radko2016'
#flag.shear_Radko2016_reduced='primitive'


if flag.flow=='HB_porous':
    flag.Nz=128
    flag.Lz=1
    flag.tau=1/5
    flag.dy_T_mean=-1
    flag.dy_S_mean=-1
    flag.Ra_T=55
    flag.Ra_S2T=5.5
    flag.kx=np.pi
    #Ra_T_list=[10000,20000,40000]
    
    #wavenumber-Ra scaling from Hewitt
    #flag.Ra_T=10000
    #flag.kx=0.48*flag.Ra_T**0.4 #2D
    #flag.Ra_S2T=0

    #wavenumber-Ra scaling from Trevisan Bejan (1987), salinity is passive scalar
    #flag.Ra_T=50#[50,100,200,400,1000]
    #flag.kx=np.pi*1#2*np.pi*[1,1.25,2,3,5.83]
    #flag.Ra_S2T=0

    #wavenumber-Ra scaling from Rosenberg & Spera (1991), salinity is also active
    #flag.Ra_T=600 #[100,150,300,600]
    #flag.kx=np.pi*1
    #R_rho_S2T=0.5
    #flag.Ra_S2T=R_rho_S2T*flag.Ra_T
    #flag.tau=1/20
    
    #wavenumber-Ra for the case of Mamou
    #flag.Ra_T=55
    #flag.Ra_S2T=0.1*flag.Ra_T
    #flag.kx=np.pi*1#This is aspect ratio A=1, which leads to a wavelength as 2...
    #flag.tau=1/5
    
    
    flag.continuation=0
    flag.ky=0
    flag.z_bc_w_left='dirichlet'
    flag.z_bc_w_right='dirichlet'  
    flag.z_bc_T_left='dirichlet'
    flag.z_bc_T_right='dirichlet'
    flag.z_bc_S_left='dirichlet'#BC for salinity
    flag.z_bc_S_right='dirichlet'
    flag.bvp_tolerance=1e-10
    
    flag.A_elevator=0.1
    flag.problem='IVP'
    flag.EVP_secondary=1
    if flag.problem =='IVP':
        flag.initial_dt=0.01 #This is the time step for double-diffusive convection in porous medium, Rosenberg case 

elif flag.flow == 'HB_porous_2_layer':
    flag.Nz=512
    flag.Lz=0.5
    flag.tau=0.01
    flag.dy_T_mean=-1
    flag.dy_S_mean=-1
    flag.Ra_T=5000
    flag.Ra_S2T=0
    flag.continuation=0
    flag.kx=0#0.48*flag.Ra_T**0.4 #2D
    flag.HB_porous_2_layer_Omega=0
    flag.problem='BVP'
    flag.A_elevator=1/10*flag.Ra_T
    flag.bvp_tolerance=1e-8
    flag.z_bc_T_left='dirichlet'
    flag.z_bc_T_right='dirichlet'
    flag.z_bc_S_left='dirichlet'
    flag.z_bc_S_right='dirichlet'
    flag.z_bc_w_left='dirichlet'
    flag.z_bc_w_right='dirichlet'
    
elif flag.flow == 'HB_porous_3_layer':
    flag.Nz=256
    flag.Lz=0.5
    flag.tau=0.01
    flag.dy_T_mean=-1
    flag.dy_S_mean=-1
    flag.Ra_T=5000
    flag.Ra_S2T=0
    flag.continuation=0
    flag.kx=2*np.pi/(4/9)#0.48*flag.Ra_T**0.4 #2D
    flag.HB_porous_3_layer_Pi=1#0.04
    flag.HB_porous_3_layer_h=0.005
    flag.problem='BVP'
    flag.A_elevator=1/10*flag.Ra_T
    flag.bvp_tolerance=1e-8
    flag.z_bc_T_left='dirichlet'
    flag.z_bc_T_right='dirichlet'
    flag.z_bc_S_left='dirichlet'
    flag.z_bc_S_right='dirichlet'
    flag.z_bc_w_left='dirichlet'
    flag.z_bc_w_right='dirichlet'
    
elif flag.flow == 'HB_porous_shear':
    flag.Nz=1024
    flag.Lz=1
    flag.tau=0.01
    flag.dy_T_mean=-1
    flag.dy_S_mean=-1
    #Ra_T_list=[10000,20000,40000]
    flag.Ra_T=10000
    flag.Ra_S2T=0
    flag.continuation=0
    flag.kx=0.48*flag.Ra_T**0.4 #2D
    flag.ky=0
    flag.problem='BVP'
    flag.z_bc_w_left='dirichlet'
    flag.z_bc_w_right='dirichlet'  
    flag.z_bc_T_left='dirichlet'
    flag.z_bc_T_right='dirichlet'
    flag.z_bc_S_left='dirichlet'#BC for salinity
    flag.z_bc_S_right='dirichlet'
    flag.bvp_tolerance=1e-10
    
    flag.A_elevator=1/10*flag.Ra_T
    
    
elif flag.flow=='HB_benard':
    flag.Nz=128
    flag.Lz=1
    #flag.tau=0.01
    #flag.Ra_T=1.25*100000
    #R_rho_T2S=1/4

    flag.tau=0.01
    flag.Ra_T=100000
    R_rho_T2S=2
    flag.continuation_asymmetric=0
    #flag.Ra_T=4*np.pi**2*Ri/(1/R_rho_T2S-1)*Pe*Pe/Pr
    flag.Ra_S2T=flag.Ra_T/R_rho_T2S
    
    flag.Pr=7
    flag.F_sin=0
    flag.ks=2*np.pi
    flag.dy_T_mean=1
    flag.dy_S_mean=1
    flag.bvp_tolerance=1e-10
    flag.EVP_homogeneous_tolerance=1e5
    #flag.kx=0.48*flag.Ra_T**0.4
    #flag.kx=2*np.pi/0.5
    flag.ky=0
    flag.EVP_secondary=1
    flag.problem='BVP'
    flag.z_bc_T_left='dirichlet'
    flag.z_bc_T_right='dirichlet'
    flag.z_bc_S_left='dirichlet'
    flag.z_bc_S_right='dirichlet'
    flag.z_bc_w_left='dirichlet'
    flag.z_bc_w_right='dirichlet'
    flag.z_bc_u_v_left='dirichlet'
    flag.z_bc_u_v_right='dirichlet'
    
    flag.A_elevator=1/100000*flag.Ra_T
    #For Nz=128, A_elevator=1 give 1-layer solution, A_elevator=650 give 2 layer solution, A_elevator=10 gives asymmetric solution
    #For Nz=1024, A_elevator=1 give 1-layer solution, A_elevator=10-100 gives 2-layer solution, A_elevator=1000 gives asymmetric solution.
    
    flag.A_noise=0.01
    if flag.problem =='IVP':
        flag.initial_dt=10**3/flag.Ra_T #This is the time step for double-diffusive convection in porous medium, Rosenberg case 

elif flag.flow in ['HB_benard_shear']:
    flag.Nz=128
    flag.Lz=1
    flag.tau=0.01
    
    # #Radko (2016) parameter
    #Pr=10
    #R_rho_T2S=0.5
    #Pe=100
    #Ri=10
    #flag.Ra_T=4*np.pi**2*Ri/(1/R_rho_T2S-1)*Pe*Pe/Pr
    #flag.Ra_S2T=0
    
    #flag.Ra_S2T=flag.Ra_T/R_rho_T2S
    #flag.kx=2*np.pi/64
    #flag.F_sin=1
    # flag.ks=2*np.pi
    # flag.HB_porous_shear_phi=0

    flag.Pr=0.05
    R_rho_T2S=40
    Pe=1
    flag.Ra_T=100000
    flag.Ra_S2T=flag.Ra_T/R_rho_T2S
    flag.F_sin=0
    
    flag.Pe_T=Pe
    flag.Pe_S=Pe
    flag.dy_T_mean=1
    flag.dy_S_mean=1
    flag.bvp_tolerance=1e-5
    #flag.kx=0.48*flag.Ra_T**0.4
    flag.kx=12
    flag.ky=0
    flag.problem='IVP'
    flag.z_bc_T_left='dirichlet'
    flag.z_bc_T_right='dirichlet'
    flag.z_bc_S_left='dirichlet'
    flag.z_bc_S_right='dirichlet'
    flag.z_bc_w_left='dirichlet'
    flag.z_bc_w_right='dirichlet'
    flag.z_bc_u_v_left='dirichlet'
    flag.z_bc_u_v_right='dirichlet'
    #flag.A_elevator=1/10000*flag.Ra_T
    #flag.A_elevator_imag=0#flag.A_elevator
    if flag.problem =='IVP':
        flag.initial_dt=0.01
        flag.post_store_dt=1
        flag.stop_sim_time=200
    
    flag.A_noise=0.001  
    """
    flag.Nz=1024
    flag.Lz=1
    flag.tau=0.01
    
    # #Radko (2016) parameter
    #Pr=10
    #R_rho_T2S=0.5
    #Pe=100
    #Ri=10
    #flag.Ra_T=4*np.pi**2*Ri/(1/R_rho_T2S-1)*Pe*Pe/Pr
    #flag.Ra_S2T=0
    
    #flag.Ra_S2T=flag.Ra_T/R_rho_T2S
    #flag.kx=2*np.pi/64
    #flag.F_sin=1
    # flag.ks=2*np.pi
    # flag.HB_porous_shear_phi=0

    Pr=10
    # Yang et al. (2021) parameter
    R_rho_T2S=0.5
    Pe=1
    Ri=1
    flag.Ra_T=100000
    flag.Ra_S2T=flag.Ra_T/R_rho_T2S
    flag.F_sin='z'
    
    flag.Pe_T=Pe
    flag.Pe_S=Pe
    flag.dy_T_mean=-1
    flag.dy_S_mean=-1
    flag.bvp_tolerance=1e-5
    #flag.kx=0.48*flag.Ra_T**0.4
    #flag.kx=1
    flag.ky=0
    flag.problem='BVP'
    flag.z_bc_T_left='dirichlet'
    flag.z_bc_T_right='dirichlet'
    flag.z_bc_S_left='dirichlet'
    flag.z_bc_S_right='dirichlet'
    flag.z_bc_w_left='dirichlet'
    flag.z_bc_w_right='dirichlet'
    flag.z_bc_u_v_left='neumann'
    flag.z_bc_u_v_right='neumann'
    flag.A_elevator=1/10000*flag.Ra_T
    flag.A_elevator_imag=0#flag.A_elevator
    if flag.problem =='IVP':
        flag.initial_dt=0.000001/flag.Ra_T
    """

elif flag.flow=='test_periodic':
    flag.F_sin=1
    flag.ks=2*np.pi
    flag.Nz=1024
    flag.Lz=1
    flag.problem='BVP'
elif flag.flow == 'IFSC_2D':
    #setup basic parameter for inertial free salt finger
    flag.Ra_ratio=2 ##This is the special parameter for the Rayleigh ratio
    flag.dy_T_mean=1
    flag.dy_S_mean=1
    
    #-------setup the grid points
    ##These are general setup
    k_opt=(1/2*(-2-flag.Ra_ratio+np.sqrt(flag.Ra_ratio**2+8*flag.Ra_ratio)))**(1/4)
    Lx2d = 8
    Lz2d = 24
    grid_l_opt=8
    flag.Lx, flag.Lz = (Lx2d*2*np.pi/k_opt, Lz2d*2*np.pi/k_opt)
    flag.Nx, flag.Nz = (grid_l_opt*Lx2d,grid_l_opt*Lz2d)
    
    ##These are setup for testing the layering based on Radko (2016)
    #flag.Lz=2*np.pi/0.3337
    #flag.Lx=flag.Lz
    #flag.Nz=32
    #flag.Nx=32

    #setup the uL and the ks, F_sin and dt...
    #u_L=9444.9
    u_L=50
    flag.ks=2*np.pi/flag.Lz
    flag.F_sin=u_L*flag.ks**2
    flag.initial_dt=np.min([flag.Lx/flag.Nx/(flag.F_sin/flag.ks**2)/flag.Ra_ratio,flag.Lx/flag.Nx/flag.Ra_ratio])

    #-----------------parameter for initial condition
    flag.A_elevator=1
    flag.A_noise=0.01
    flag.A_shear=0
    
elif flag.flow == 'double_diffusive_2D':
    #setup basic parameter for inertial free salt finger
    flag.tau=0.01
    flag.Pr=10
    flag.R_rho_T2S=0.5
    flag.dy_T_mean=-1#-------------These values as 1 corresponds to salt finger and -1 corresponds to diffusive regime
    flag.dy_S_mean=-1
    
    ##These are setup for testing the layering based on Radko (2016)
    Lx2d = 1
    Lz2d = 2
    flag.ks=0.0593
    #u_L=0.0593
    u_L=9444.9*flag.tau
    flag.F_sin=u_L*flag.ks**2

    flag.Lz=Lx2d*2*np.pi/flag.ks
    flag.Lx=flag.Lz
    flag.Nz=Lx2d*128
    flag.Nx=Lz2d*128
    
    #setup the uL and the ks, F_sin and dt...
    #u_L=531.126*flag.tau
    #u_L=0
    
    #Here, use the np.divide so divide by zero will give Inf...

    flag.initial_dt=np.min([np.divide(flag.Lx/flag.Nx,u_L),flag.Lx/flag.Nx])
    #print(initial_dt)
    #u_L=9444.9*flag.tau
    
    #-----------------parameter for initial condition
    #flag.A_elevator=0.1
    #flag.k_elevator=0.5
    flag.A_noise=1
    flag.A_shear=1

elif flag.flow == 'porous_media_2D':
    Lx2d = 4
    Lz2d = 4
    flag.Lx=Lx2d*2*np.pi
    flag.Lz=Lz2d*2*np.pi
    flag.Nx=Lx2d*16
    flag.Nz=Lz2d*16
    flag.initial_dt=np.min([flag.Lx/flag.Nx])
    flag.Ra_T=1
    flag.A_elevator=2**8
    flag.k_elevator=1
    flag.A_secondary_T=1
    flag.k_secondary=0.1795
    
    
elif flag.flow == 'HB_benard_shear_periodic':   
    flag.Nz=128
    flag.Lz=1
    flag.tau=1
    flag.Pr=1
    flag.flux_T=1
    flag.Pe_S=1
    
    #Thermal diffusivity unit
    #flag.Re=1/flag.Pr
    #flag.Pe_T=1
    
    #Thermal diffusivity unit, infinite Prandtl number
    #flag.Re=0
    #flag.Pe_T=1
    
    #Viscous unit, zero Prandtl number
    flag.Re=1
    flag.Pe_T=0
    
    
    flag.Ra_T=4*10**4
    flag.Ra_S2T=0#flag.Ra_T/R_rho_T2S
    
    #flag.kx=0.48*flag.Ra_T**0.4
    flag.kx=10
    flag.ky=0
    flag.problem='IVP'
    flag.z_bc_T_left='periodic'
    flag.z_bc_T_right='periodic'
    flag.z_bc_S_left='periodic'
    flag.z_bc_S_right='periodic'
    flag.z_bc_w_left='periodic'
    flag.z_bc_w_right='periodic'
    flag.z_bc_u_v_left='periodic'
    flag.z_bc_u_v_right='periodic'
    flag.z_basis_mode='Fourier'
    if flag.problem =='IVP':
        flag.initial_dt=1e-3
        flag.post_store_dt=0.01
        flag.stop_sim_time=20
    
    kx_2D=np.sqrt(flag.kx*flag.kx+flag.ky*flag.ky)
    Lx2d=1
    n_elevator=1
     
    flag.k_elevator=n_elevator*flag.kx
    flag.A_secondary_phase=0#phase for the second half domain
    #w_hat=np.sqrt((1-flag.k_elevator**4/flag.Ra_T)/(2*flag.k_elevator**2/flag.Ra_T))
    w_hat=np.sqrt((1-flag.kx**4/flag.Ra_T)/(2*flag.kx**2/flag.Ra_T))
    #flag.A_elevator=2*w_hat
    flag.A_w_hat=w_hat
    #flag.dy_T_mean=-flag.k_elevator**4/flag.Ra_T

    flag.A_noise=0
    flag.store_variable='T_u_w'#only store S and u variable
    flag.S_active=0
    flag.A_w_mean=0 #This is mean vertical velocity
    flag.A_u_mean=0 #This is the mean horizontal velocity
    flag.timesteppers ='RK443'


elif flag.flow == 'double_diffusive_shear_2D':
    ##These are setup for testing the layering based on Radko (2016)
    #Lx2d = 16
    #Lz2d = 16
    #flag.ks=0.5
    #flag.Lx=Lx2d*2*np.pi/flag.ks
    #flag.Lz=Lz2d*2*np.pi/flag.ks
    #flag.Nx=Lx2d*16
    #flag.Nz=Lz2d*16
    #u_L=0
    #flag.initial_dt=np.min([np.divide(flag.Lx/flag.Nx,u_L),flag.Lx/flag.Nx])
    #flag.initial_dt=0.01
    if flag.flow_sub_double_diffusive_shear_2D == 'primitive_Radko2013':
        ##parameter for Radko (2013) type
        Pr=10
        tau=0.01
        R_rho_T2S=2
        
        #map to the extended parameter in double_diffusive_shear_2D
        flag.Re=1/Pr
        flag.Pe_T=1
        flag.Pe_S=1
        flag.tau=tau #Set this as zero if remove salinity diffusivity
        flag.Ra_T=1
        flag.Ra_S2T=1/R_rho_T2S
        
        flag.dy_T_mean=1
        flag.dy_S_mean=1
        
        flag.A_elevator=1
        flag.k_elevator=1
        
        #flag.
        #flag.A_secondary_T=0.1
        
    elif flag.flow_sub_double_diffusive_shear_2D == 'primitive_IFSC_unit_tuS':
        ##parameter for primitive equations
        Pr=10
        tau=0.01
        R_rho_T2S=50
        
        #map to the extended parameter in primitive_IFSC_unit_tus
        flag.Re=tau/Pr
        flag.Pe_T=tau
        flag.Pe_S=1
        flag.tau=1 #Set this as zero if remove salinity diffusivity
        flag.Ra_T=1
        flag.Ra_S2T=1/R_rho_T2S/tau
        
        flag.dy_T_mean=1
        flag.dy_S_mean=1
        
        flag.A_elevator=0.1
        flag.k_elevator=0.5
        
    elif flag.flow_sub_double_diffusive_shear_2D == 'IFSC':
    ##parameter for 
        Ra_ratio=2 ##parameter of IFSC
        
        #map to the extended parameter in double_diffusive_shear_2D
        flag.Re=0
        flag.Pe_T=0
        flag.Pe_S=1
        flag.tau=1 #Set this as zero if remove salinity diffusivity
        flag.Ra_T=1
        flag.Ra_S2T=Ra_ratio
        
        flag.dy_T_mean=1
        flag.dy_S_mean=1
        
        flag.A_elevator=0.1
        flag.k_elevator=0.5
    elif flag.flow_sub_double_diffusive_shear_2D == 'MRBC':
        Ra_ratio=2
        Sc=1000
        
        #map to the extended parameter in double_diffusive_shear_2D
        flag.Re=1/Sc
        flag.Pe_T=0
        flag.Pe_S=1
        flag.tau=1 #Set this as zero if remove salinity diffusivity
        flag.Ra_T=1
        flag.Ra_S2T=Ra_ratio
    
        flag.dy_T_mean=1
        flag.dy_S_mean=1
        
        flag.A_elevator=0.1
        flag.k_elevator=0.5
    elif flag.flow_sub_double_diffusive_shear_2D == 'Stokes':
        #map to the extended parameter in double_diffusive_shear_2D
        Ra_ratio=2
        tau=0.01
        
        flag.Re=0
        flag.Pe_T=tau
        flag.Pe_S=1
        flag.tau=1 #Set this as zero if remove salinity diffusivity
        flag.Ra_T=1
        flag.Ra_S2T=Ra_ratio
    
        flag.dy_T_mean=1
        flag.dy_S_mean=1
        
        flag.A_elevator=0.1
        flag.k_elevator=0.5
        
        
    elif flag.flow_sub_double_diffusive_shear_2D == 'shear_Radko2016':
        flag.ks=2*np.pi
        flag.Lx=64
        flag.Lz=1
        flag.Nz=24
        flag.Nx=512
        u_L=1
        flag.F_sin=u_L*flag.ks*flag.ks
        flag.initial_dt=np.min([np.divide(flag.Lx/flag.Nx,u_L)/100,flag.Lx/flag.Nx])
        
        Pr=10
        tau=0.01
        R_rho_T2S=0.5
        Pe=100
        Ri=10
        
        #map to the extended parameter in double_diffusive_shear_2D
        if flag.shear_Radko2016_reduced =='primitive':
            flag.Re=Pe/Pr
            flag.Pe_T=Pe
        elif flag.shear_Radko2016_reduced =='IFSC':
            flag.Re=0
            flag.Pe_T=0
        elif flag.shear_Radko2016_reduced =='MRBC':
            flag.Re=Pe/Pr
            flag.Pe_T=0
        elif flag.shear_Radko2016_reduced =='Stokes':
            flag.Re=0
            flag.Pe_T=Pe

        flag.Pe_S=Pe
        flag.tau=tau
        flag.Ra_T=4*np.pi**2*Ri/(1/R_rho_T2S-1)*Pe*Pe/Pr
        flag.Ra_S2T=flag.Ra_T/R_rho_T2S
        
        flag.dy_T_mean=-1
        flag.dy_S_mean=-1
        
        flag.A_noise=0.001
        flag.A_shear=1
        #flag.A_elevator=0.1
        #flag.k_elevator=0.5
    elif flag.flow_sub_double_diffusive_shear_2D=='primitive_dirichlet_salt_finger':
        ##parameter for Radko (2013) type
        flag.Pr=0.05
        #flag.tau=0.02944
        #R_rho_T2S=20
        flag.tau=0.01
        #R_rho_T2S=40
        flag.initial_dt=0.01
        
        
        #map to the extended parameter in double_diffusive_shear_2D
        flag.Re=1/flag.Pr
        flag.Pe_T=1
        flag.Pe_S=1
        #flag.tau=tau #Set this as zero if remove salinity diffusivity
        flag.Ra_T=10**5
        flag.Ra_S2T=2500#flag.Ra_T#flag.Ra_T/R_rho_T2S
        R_rho_T2S=flag.Ra_T/flag.Ra_S2T
        #I need to overwrite these domain setup here
        Ra_S=flag.Ra_S2T/flag.tau
        flag.kx=6#2*np.pi/(2*14.8211*Ra_S**(-0.2428)/R_rho_T2S**(0.25/2))
        flag.ky=0
        kx_2D=np.sqrt(flag.kx*flag.kx+flag.ky*flag.ky)
        Lx2d=1
        flag.Lx=Lx2d*2*np.pi/kx_2D
        flag.Lz=1
        flag.Nx=128
        flag.Nz=128
         
        flag.dy_T_mean=1
        flag.dy_S_mean=1
        
        flag.z_bc_u_v_left='dirichlet' #This can be periodic, dirichlet, or neumann
        flag.z_bc_T_left='dirichlet'
        flag.z_bc_S_left='dirichlet'
        flag.z_bc_w_left='dirichlet'
        flag.z_bc_u_v_right='dirichlet' #This can be periodic, dirichlet, or neumann
        flag.z_bc_T_right='dirichlet'
        flag.z_bc_S_right='dirichlet'
        flag.z_bc_w_right='dirichlet'
        
        flag.A_elevator=0
        flag.k_elevator=1
        
        flag.A_noise=0
        flag.A_secondary_S=0
        flag.k_secondary=2*np.pi #4*np.pi, or 6*np.pi, will give 2 or 3 staircase
        flag.store_variable='S_u_w'#only store S and u variable
        
        
    elif flag.flow_sub_double_diffusive_shear_2D=='primitive_stress_free_salt_finger':
        ##parameter for Radko (2013) type
        flag.Pr=7
        #flag.tau=0.02944
        #R_rho_T2S=20
        flag.tau=0.01
        #R_rho_T2S=40
        flag.initial_dt=0.001
        
        
        #map to the extended parameter in double_diffusive_shear_2D
        flag.Re=1/flag.Pr
        flag.Pe_T=1
        flag.Pe_S=1
        #flag.tau=tau #Set this as zero if remove salinity diffusivity
        flag.Ra_T=10**5
        flag.Ra_S2T=2500#flag.Ra_T#flag.Ra_T/R_rho_T2S
        R_rho_T2S=flag.Ra_T/flag.Ra_S2T
        #I need to overwrite these domain setup here
        Ra_S=flag.Ra_S2T/flag.tau
        flag.kx=14#2*np.pi/(2*14.8211*Ra_S**(-0.2428)/R_rho_T2S**(0.25/2))
        flag.ky=0
        kx_2D=np.sqrt(flag.kx*flag.kx+flag.ky*flag.ky)
        Lx2d=1
        flag.Lx=Lx2d*2*np.pi/kx_2D
        flag.Lz=1
        flag.Nx=128
        flag.Nz=128
         
        flag.dy_T_mean=1
        flag.dy_S_mean=1
        
        flag.z_bc_u_v_left='neumann' #This can be periodic, dirichlet, or neumann
        flag.z_bc_T_left='dirichlet'
        flag.z_bc_S_left='dirichlet'
        flag.z_bc_w_left='dirichlet'
        flag.z_bc_u_v_right='neumann' #This can be periodic, dirichlet, or neumann
        flag.z_bc_T_right='dirichlet'
        flag.z_bc_S_right='dirichlet'
        flag.z_bc_w_right='dirichlet'
        
        flag.A_elevator=0
        flag.k_elevator=1
        
        flag.A_noise=0
        flag.A_secondary_S=0
        flag.k_secondary=2*np.pi #4*np.pi, or 6*np.pi, will give 2 or 3 staircase
        flag.store_variable='all'#only store S and u variable
        
        
    elif flag.flow_sub_double_diffusive_shear_2D=='primitive_periodic_salt_finger':
        ##parameter for Radko (2013) type
        flag.Pr=7
        #flag.tau=0.02944
        #R_rho_T2S=20
        flag.tau=0.01
        #R_rho_T2S=40
        flag.initial_dt=0.001
        
        
        #map to the extended parameter in double_diffusive_shear_2D
        flag.Re=1/flag.Pr
        flag.Pe_T=1
        flag.Pe_S=1
        #flag.tau=tau #Set this as zero if remove salinity diffusivity
        flag.Ra_T=10**5
        flag.Ra_S2T=flag.Ra_T/40#flag.Ra_T#flag.Ra_T/R_rho_T2S
        R_rho_T2S=flag.Ra_T/flag.Ra_S2T
        #I need to overwrite these domain setup here
        Ra_S=flag.Ra_S2T/flag.tau
        flag.kx=18#2*np.pi/(2*14.8211*Ra_S**(-0.2428)/R_rho_T2S**(0.25/2))
        flag.ky=0
        kx_2D=np.sqrt(flag.kx*flag.kx+flag.ky*flag.ky)
        Lx2d=1
        flag.Lx=Lx2d*2*np.pi/kx_2D
        flag.Lz=1
        flag.Nx=128
        flag.Nz=128
         
        flag.dy_T_mean=1
        flag.dy_S_mean=1
                
        flag.z_bc_u_v_left='periodic' #This can be periodic, dirichlet, or neumann
        flag.z_bc_T_left='periodic'
        flag.z_bc_S_left='periodic'
        flag.z_bc_w_left='periodic'
        flag.z_bc_u_v_right='periodic' #This can be periodic, dirichlet, or neumann
        flag.z_bc_T_right='periodic'
        flag.z_bc_S_right='periodic'
        flag.z_bc_w_right='periodic'
        
        flag.A_elevator=0
        
        flag.A_noise=0
        #flag.store_variable='S_u_w'#only store S and u variable
        
    elif flag.flow_sub_double_diffusive_shear_2D=='primitive_periodic_RBC':
        ##parameter for Radko (2013) type
        #flag.damping_1_beta=10
        flag.Pr=1
        flag.tau=1
        #R_rho_T2S=40
        flag.flux_T=1
        flag.Pe_S=1

        #finite Prandtl number, map to the extended parameter in double_diffusive_shear_2D
        flag.Re=1/flag.Pr
        flag.Pe_T=1
        #flag.Ta_sqrt_z=10 #Taylor number for rotation in 

        #zero Prandtl number
        #flag.Re=1
        #flag.Pe_T=0
        
        #infinite Prandtl number
        #flag.Re=0
        #flag.Pe_T=1
        
        #flag.tau=tau #Set this as zero if remove salinity diffusivity
        #flag.Ra_T=2*10**4
        #flag.Ra_T=6*10**4
        flag.Ra_T=10**8
        
        #flag.initial_dt=10**(-3)
        flag.initial_dt=10**(-6)
        #flag.initial_dt=50/flag.Ra_T

        flag.Ra_S2T=0#flag.Ra_T#flag.Ra_T/R_rho_T2S
        #Ra_S=flag.Ra_S2T/flag.tau
        flag.kx=10#2*np.pi#2*np.pi/(2*14.8211*Ra_S**(-0.2428)/R_rho_T2S**(0.25/2))
        flag.ky=0
        kx_2D=np.sqrt(flag.kx*flag.kx+flag.ky*flag.ky)
        Lx2d=1
        flag.Lx=Lx2d*2*np.pi/kx_2D
        flag.Lz=1
        flag.Nx=256
        flag.Nz=256
        n_elevator=1
         
        flag.dy_T_mean=-flag.kx**4/flag.Ra_T
        flag.dy_S_mean=0
                
        flag.z_bc_u_v_left='periodic' #This can be periodic, dirichlet, or neumann
        flag.z_bc_T_left='periodic'
        flag.z_bc_S_left='periodic'
        flag.z_bc_w_left='periodic'
        flag.z_bc_u_v_right='periodic' #This can be periodic, dirichlet, or neumann
        flag.z_bc_T_right='periodic'
        flag.z_bc_S_right='periodic'
        flag.z_bc_w_right='periodic'
        
        flag.k_elevator=n_elevator*flag.kx
        flag.A_secondary_phase=0#phase for the second half domain
        #w_hat=np.sqrt((1-flag.k_elevator**4/flag.Ra_T)/(2*flag.k_elevator**2/flag.Ra_T))
        w_hat=np.sqrt((1-flag.kx**4/flag.Ra_T)/(2*flag.kx**2/flag.Ra_T))
        flag.A_elevator=2*w_hat
        flag.dy_T_mean=-flag.k_elevator**4/flag.Ra_T
        
        flag.A_secondary_w=0#10
        flag.A_secondary_U0=0#10
        flag.A_secondary_T=0#1
        flag.k_secondary=0#2*np.pi
        flag.A_noise=0
        flag.store_variable='T_u_v_w'#only store S and u variable
        flag.S_active=0
        flag.A_w_mean=0 #This is mean vertical velocity
        
        flag.nx_trunc_num=0
        flag.nz_trunc_num=0
        flag.timesteppers ='RK443'
        
    else:
        raise TypeError('flag.flow_sub_double_diffusive_shear_2D is not found')
#--------------setup the background shear

#u_L_2ks=0
#u_L_3ks=0
#u_L_4ks=0


#flag.F_sin_2ks=u_L_2ks*(2*flag.ks)**2
#flag.F_sin_3ks=u_L_3ks*(3*flag.ks)**2
#flag.F_sin_4ks=u_L_4ks*(4*flag.ks)**2

#flag.phase_2ks=0
#flag.phase_3ks=0
#flag.phase_4ks=0


#-----------------setup storing for post-processing
if flag.flow in ['HB_benard_shear']:
    print('1')
    #flag.post_store_dt=0.01/flag.Ra_T
    #flag.stop_sim_time=10/flag.Ra_T
elif flag.flow in ['HB_benard']:
    flag.post_store_dt=10**3/flag.Ra_T
    flag.stop_sim_time=10**7/flag.Ra_T
elif flag.flow in ['HB_porous']:
    flag.post_store_dt=0.01
    flag.stop_sim_time=40
elif flag.flow == 'double_diffusive_shear_2D':
    if flag.flow_sub_double_diffusive_shear_2D=='primitive_dirichlet_salt_finger':
        flag.post_store_dt=1
        flag.stop_sim_time=4000
    elif flag.flow_sub_double_diffusive_shear_2D=='primitive_stress_free_salt_finger':
        flag.post_store_dt=0.001
        flag.stop_sim_time=0.01
    elif flag.flow_sub_double_diffusive_shear_2D=='primitive_periodic_salt_finger':
        flag.post_store_dt=0.3
        flag.stop_sim_time=30
    elif flag.flow_sub_double_diffusive_shear_2D=='primitive_periodic_RBC':
        #low Ra, fixed flux RBC
        #flag.post_store_dt=0.01
        #flag.stop_sim_time=10
        
        flag.post_store_dt=0.001
        flag.stop_sim_time=1
        
        #high Ra, fixed flux RBC
        #may need several run
        #flag.post_store_dt=0.001
        #flag.stop_sim_time=0.1
        
        #test
        #flag.post_store_dt=flag.initial_dt
        #flag.stop_sim_time=100*flag.initial_dt
        
        #flag.post_store_dt=10**4/flag.Ra_T
        #flag.stop_sim_time=10**7/flag.Ra_T
        
        
        #flag.post_store_dt=50/flag.Ra_T
        #flag.stop_sim_time=100/flag.Ra_T 
else:
    print('1')
    #flag.post_store_dt=0.000001/flag.Ra_T;
    #flag.stop_sim_time=0.00001/flag.Ra_T;

#for flag.Ra_T in [40000,35000,33000,32260,32255,32250,32000,31500,31000,30500,30000]:
#for flag.Ra_T in reversed([40000,35000,33000,32260,32255,32250,32000,31500,31000,30500,30000]):
domain=flag.build_domain()
solver=flag.governing_equation(domain)
flag.print_screen(logger)
flag.initial_condition(domain,solver)
flag.post_store(solver)
flag.print_file() #move print file to here.
flag.run(solver,domain,logger)
flag.post_store_after_run(solver)
    #flag.continuation=flag.continuation+1

#------------ print these parameters in the screen
#flag.ky=flag.kx                        
#Ra_T_list=[10000,20000,40000]
#np.linspace(0,20000,21)
#Ra_T_list=np.logspace(2,6,100)

#####This is the for loop for different Ra
#for flag.Ra_T in Ra_T_list:
#Lz_list=np.linspace(1,64,65)
#Pi_list=[1,0.5,0.2,0.1,0.04]
#for flag.HB_porous_3_layer_Pi in Pi_list:
#    flag.Pe_T=Pe
#    flag.Pe_S=Pe
#F_sin_list=np.linspace(0.1,1,19)
#####This is the for loop for different shear amplitude
#for flag.F_sin in F_sin_list:
#    flag.kx=0.48*flag.Ra_T**0.4 #2D
    #flag.kx=0.17*flag.Ra_T**0.52 #3D
#####This is for loop for the HB_porous_2_layer with different Omega    
#Omega_list,kx_list = flag.get_HB_porous_2_layer_Omega_k()
#for index, flag.HB_porous_2_layer_Omega in enumerate(Omega_list):
#    flag.kx=kx_list[index]
#for flag.Ra_S2T in np.linspace(0,20000,21):
#for flag.tau in np.divide(1,[0.02,0.04,0.1,0.2,0.4,1,2,4,10,20,40,100]):
#for R_rho_S2T in [0,0.1,0.2,0.3,0.4]:
    #flag.Ra_S2T=R_rho_S2T*flag.Ra_T
    
  
#--------Loop for solving the boundary value problem    
#This is the loop for the Yang (2015) comparison of bounded salt finger    
#for R_rho_T2S in [10,5,2,1,0.5,0.2,0.1]:    
"""
for R_rho_T2S in [2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10]:#[10,5,2,1,0.5,0.2,0.1]
#for R_rho_T2S in [2,1.8,1.6,1.4,1.3,1,0.8,0.6,0.4,0.2,0.1]:
    flag.Ra_S2T=flag.Ra_T/R_rho_T2S #10^6, 2*10^6, 5*10^6, 10^7, 2*10^7, 5*10^7, 10^8
    Ra_S=flag.Ra_S2T/flag.tau
    kx_final=2*np.pi/(2*14.8211*Ra_S**(-0.2428)/R_rho_T2S**(0.25/2))
    flag.kx=kx_final
#This is just try to study the wavenumber 
#for flag.kx in np.linspace(1,kx_final,10):
    #R_rho_T2S=2
    #flag.Ra_S2T=flag.Ra_T/R_rho_T2S
    
#This is just for nothing need to loop    
#for flag.ky in [0]:    
    flag.ky=0
    flag.kx_2=0
    flag.ky_2=flag.kx_2
    domain=flag.build_domain()
    solver=flag.governing_equation(domain)
    flag.print_screen(logger)
    flag.initial_condition(domain,solver)
    flag.post_store(solver)
    flag.print_file() #move print file to here.
    flag.run(solver,domain,logger)
    flag.post_store_after_run(solver)
    flag.continuation=flag.continuation+1
"""

##Old one try different B.C. and second harmonic
#---------main loop to run the dedalus 
# for bc in ['dirichlet']:
#     if bc in ['dirichlet','periodic']:
#         flag.z_bc_T=bc
#         flag.z_bc_S=bc
#         flag.z_bc_w=bc
#         flag.z_bc_u_v='neumann'
#     elif bc == 'neumann':
#         flag.z_bc_T=bc
#         flag.z_bc_S=bc
#         flag.z_bc_w='dirichlet'
#         flag.z_bc_u_v='dirichlet'
        
#     #for flag.Ra_T in Ra_T_list:
#     #flag.kx=0.48*flag.Ra_T**0.4 #2D
#     flag.kx=0.17*flag.Ra_T**0.52 #3D
#     flag.ky=flag.kx
#     for flag.kx_2 in [0,10,20,30,40,50,60,70,80]:
#         flag.ky_2=flag.kx_2
#         domain=flag.build_domain()
#         solver=flag.governing_equation(domain)
#         flag.print_screen(logger)
#         flag.initial_condition(domain,solver)
#         flag.post_store(solver)
#         flag.print_file() #move print file to here.
#         flag.run(solver,domain,logger)
#         flag.post_store_after_run(solver)
#         flag.continuation=flag.continuation+1



# flag.Ra_S2T=11000
# flag.continuation=1
# domain=flag.build_domain()
# solver=flag.governing_equation(domain)
# flag.initial_condition(domain,solver)
# flag.post_store(solver)
# flag.print_file() #move print file to here.
# flag.run(solver,domain,logger)
# flag.post_store_after_run(solver)

#-----------merge process data


#if flag.problem == 'IVP':
#    domain=flag.build_domain()
#    solver=flag.governing_equation(domain)
    #ts = de.timesteppers.RK443
    #solver =  problem.build_solver(ts)
#    flag.initial_condition(domain,solver)
#    flag.post_store(solver)
#    flag.print_file() #move print file to here.
#    flag.run(solver,domain,logger)
#    flag.post_store_after_run(solver)
#elif flag.problem =='BVP':