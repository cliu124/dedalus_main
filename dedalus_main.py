import numpy as np
from mpi4py import MPI
import pathlib
#from IPython import display
from dedalus import public as de
from dedalus.extras import flow_tools
from dedalus.tools import post
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


#-------setup the grid points
##These are general setup
#k_opt=(1/2*(-2-flag.Ra_ratio+np.sqrt(flag.Ra_ratio**2+8*flag.Ra_ratio)))**(1/4)
#Lx2d = 8
#Lz2d = 24
#grid_l_opt=8
#flag.Lx, flag.Lz = (Lx2d*2*np.pi/k_opt, Lz2d*2*np.pi/k_opt)
#flag.Nx, flag.Nz = (grid_l_opt*Lx2d,grid_l_opt*Lz2d)

##These are setup for testing the layering based on Radko (2016)
flag.Lz=2*np.pi/0.3337
flag.Lx=flag.Lz
flag.Nz=32
flag.Nx=32

#------------select the flow configuration and special parameters for each
flag.flow='IFSC_2D'

if flag.flow == 'IFSC_2D':
    #setup basic parameter for inertial free salt finger
    flag.Ra_ratio=200 ##This is the special parameter for the Rayleigh ratio
    flag.dy_T_mean=-1
    flag.dy_S_mean=-1
    #setup the uL and the ks, F_sin and dt...
    #u_L=9444.9
    u_L=531.126
    flag.ks=2*np.pi/flag.Lz
    flag.F_sin=u_L*flag.ks**2
    initial_dt=np.min([flag.Lx/flag.Nx/(flag.F_sin/flag.ks**2)/flag.Ra_ratio,flag.Lx/flag.Nx/flag.Ra_ratio])

    
elif flag.flow == 'double_diffusive_2D':
    #setup basic parameter for inertial free salt finger
    flag.tau=0.01
    flag.Pr=10
    flag.R_rho_T2S=0.5
    flag.dy_T_mean=-1#-------------These values as 1 corresponds to salt finger and -1 corresponds to diffusive regime
    flag.dy_S_mean=-1
    #setup the uL and the ks, F_sin and dt...
    #u_L=9444.9*flag.tau
    u_L=531.126*flag.tau
    flag.ks=2*np.pi/flag.Lz
    flag.F_sin=u_L*flag.ks**2
    initial_dt=np.min([flag.Lx/flag.Nx/(flag.F_sin/flag.ks**2),flag.Lx/flag.Nx])

    #u_L=9444.9*flag.tau

#--------------setup the background shear

u_L_2ks=0
u_L_3ks=0
u_L_4ks=0


flag.F_sin_2ks=u_L_2ks*(2*flag.ks)**2
flag.F_sin_3ks=u_L_3ks*(3*flag.ks)**2
flag.F_sin_4ks=u_L_4ks*(4*flag.ks)**2

flag.phase_2ks=0
flag.phase_3ks=0
flag.phase_4ks=0



#-----------------parameter for initial condition
flag.A_elevator=0
flag.A_noise=0.01
flag.A_shear=1

#-----------------setup storing for post-processing
flag.post_store_dt=0.001;
flag.stop_sim_time=0.01;

#------------ print these parameters in the screen
flag.print_screen(logger)


#---------main loop to run the dedalus 
domain=flag.build_domain()
problem=flag.governing_equation(domain)
ts = de.timesteppers.RK443
solver =  problem.build_solver(ts)
flag.initial_condition(domain,solver)
solver.stop_sim_time = flag.stop_sim_time
cfl = flow_tools.CFL(solver,initial_dt,safety=0.8,max_change=1,cadence=8)
flag.post_store(solver)
flag.run(solver,cfl,domain,logger)

#-----------merge process data
post.merge_process_files('analysis',cleanup=True)
flag.print_file()


# Create bases and domain
#x_basis = de.Fourier('x', nx, interval=(0, Lx), dealias=3/2)
#z_basis = de.Fourier('z',nz, interval=(0,Lz), dealias=3/2)
#domain = de.Domain([x_basis,  z_basis], grid_dtype=np.float64)

#problem = de.IVP(domain, variables=['p','u','w','S','T'])
#problem.parameters['Ra_ratio']=flag.Ra_ratio
#problem.parameters['pi']=np.pi
#problem.parameters['R_rho']=1/flag.Ra_ratio/tau
#problem.parameters['Pr']=Pr
#problem.parameters['tau']=tau
#problem.parameters['ks']=ks
#F_sin=ks/tau*np.sqrt(Pr*(1-1/R_rho)/Ri)
#logger.info('The shear forcing amplitude is %f' %(F_sin))

#problem.parameters['F_sin']=F_sin
#problem=dedalus_setup.governing_equation(problem,'IFSC_2D_with_shear')

#F_sin=1
#Ra_ratio=1
##These two equations are dealing with the singularity of the following two equations when kx=kx=0
#problem.add_equation("p=0",condition="(nx==0) and (nz==0)")
#problem.add_equation("u=0",condition="(nx==0) and (nz==0)")

#problem.add_equation("dx(u) + dz(w) = 0",condition="(nx!=0) or (nz!=0)") ##divergence free constraint is singular at kx=kz=0
#problem.add_equation("dt(S) - (dx(dx(S)) + dz(dz(S)))  + w                  = - u*dx(S) - w*dz(S)")

##with shear
#problem.parameters['Lz']=Lz
#problem.parameters['pi']=np.pi
#problem.add_equation(" -( dx(dx(u)) + dz(dz(u)) ) + dx(p) =F_sin*sin(ks*z)",condition="(nx!=0) or (nz!=0)")


##without shear
#problem.add_equation(" -( dx(dx(u)) + dz(dz(u)) ) + dx(p) =0",condition="(nx!=0) or (nz!=0)")
#This is also singular at kx=kz=0, this is unique for the IFSC, because the momentum is reduced as the algebraic constraint
#problem.add_equation(" -( dx(dx(w)) + dz(dz(w)) ) + dz(p) -(T-S*Ra_ratio)   =0")
#problem.add_equation(" -( dx(dx(T)) + dz(dz(T)) ) + w                       =0")



#flag.print_file()
#print(', '.join("%s: %s, \n" % item for item in flag_attrs.items()))

#Ra_ratio=1.1

#merge time snapshot
#set_paths = list(pathlib.Path("IFSC_2D_without_shear").glob("IFSC_2D_without_shear_s*.h5"))
#post.merge_sets("IFSC_2D_without_shear/IFSC_2D_without_shear.h5",set_paths,cleanup=True) 

#print the information for post-processing
#print(subprocess.check_output("find IFSC_2D_without_shear", shell=True).decode())


#tau=0.01
#Pr=7

#R_rho=1/tau/flag.Ra_ratio


#Ri=1/1000*Pr*(1-1/R_rho)/tau**2/ks**2