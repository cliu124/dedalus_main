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
flag.flow='HB_porous'
#flag.flow='double_diffusive_shear_2D'#['IFSC_2D','double_diffusive_2D','double_diffusive_shear_2D','porous_media_2D']
#flag.flow='porous_media_2D'
flag.flow_sub_double_diffusive_shear_2D='double_diffusive'
flag.flow_sub_double_diffusive_shear_2D='IFSC'
flag.flow_sub_double_diffusive_shear_2D='MRBC'
flag.flow_sub_double_diffusive_shear_2D='Stokes'
flag.flow_sub_double_diffusive_shear_2D='primitive_IFSC_unit_tuS'
flag.flow_sub_double_diffusive_shear_2D='shear_Radko2016'
flag.shear_Radko2016_reduced='primitive'


if flag.flow=='HB_porous':
    flag.Nz=1024
    flag.Lz=1
    flag.tau=0.01
    flag.dy_T_mean=-1
    flag.dy_S_mean=-1
    flag.Ra_T=10000
    #Ra_T_list=[10000,20000,40000]
    flag.continuation=0
    flag.Ra_S2T=0
    flag.ky=0
    flag.problem='BVP'
    flag.z_bc_T='dirichlet'
    flag.z_bc_S='dirichlet'
    flag.z_bc_w='dirichlet'
    flag.z_bc_u_v='dirichlet'
    flag.A_elevator=1/10*flag.Ra_T
    flag.A_noise=0
    flag.initial_dt=0.0001
elif flag.flow=='HB_benard':
    flag.Nz=1028
    flag.Lz=1
    flag.tau=0.01
    flag.dy_T_mean=-1
    flag.dy_S_mean=-1
    flag.Ra_T=10000
    flag.Ra_S2T=0
    flag.kx=0.48*flag.Ra_T**0.4
    flag.ky=0
    flag.problem='BVP'
    flag.z_bc_T_S_w='dirichlet'
    flag.z_bc_u_v='dirichlet'
    flag.A_elevator=1/10*flag.Ra_T
    flag.A_noise=0
    flag.initial_dt=0.0001
    
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
    

elif flag.flow == 'double_diffusive_shear_2D':
    ##These are setup for testing the layering based on Radko (2016)
    Lx2d = 16
    Lz2d = 16
    flag.ks=0.5
    flag.Lx=Lx2d*2*np.pi/flag.ks
    flag.Lz=Lz2d*2*np.pi/flag.ks
    flag.Nx=Lx2d*16
    flag.Nz=Lz2d*16
    u_L=0
    flag.initial_dt=np.min([np.divide(flag.Lx/flag.Nx,u_L),flag.Lx/flag.Nx])

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
        
        flag.A_secondary_T=0.1
        flag.k_secondary_T=0.5
        
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
flag.post_store_dt=0.5;
flag.stop_sim_time=20;

#------------ print these parameters in the screen
flag.print_screen(logger)


#---------main loop to run the dedalus 
for bc in ['dirichlet','neumann']:
    flag.z_bc_T=bc
    flag.z_bc_S=bc
    #for flag.Ra_T in Ra_T_list:
    flag.kx=0.48*flag.Ra_T**0.4
    domain=flag.build_domain()
    solver=flag.governing_equation(domain)
    flag.initial_condition(domain,solver)
    flag.post_store(solver)
    flag.print_file() #move print file to here.
    flag.run(solver,domain,logger)
    flag.post_store_after_run(solver)
    flag.continuation=flag.continuation+1



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