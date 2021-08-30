import numpy as np
from mpi4py import MPI
import time
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
flag.Ra_ratio=1.1
flag.name='IFSC_2D_without_shear'
flag.path='./'

flag.print_screen()
flag.current_path='/projects/chli3324/dedalus/'
flag.flow='IFSC_2D_without_shear'

#flag.print_file()
#print(', '.join("%s: %s, \n" % item for item in flag_attrs.items()))

#Ra_ratio=1.1
tau=0.01
Pr=7

R_rho=1/tau/flag.Ra_ratio

k_opt=(1/2*(-2-flag.Ra_ratio+np.sqrt(flag.Ra_ratio**2+8*flag.Ra_ratio)))**(1/4)
logger.info('The optimal horizontal wavenumber is %f for Rayleigh ratio %f' %(k_opt,flag.Ra_ratio))

Lx2d = 32
Lz2d = 32
grid_l_opt=8
Lx, Lz = (Lx2d*2*np.pi/k_opt, Lz2d*2*np.pi/k_opt)
nx, nz = (grid_l_opt*Lx2d,grid_l_opt*Lz2d)

ks=2*np.pi/Lz

#Ri=1/1000*Pr*(1-1/R_rho)/tau**2/ks**2
Ri=1
logger.info('The ks=%f, Ri=%f ' %(ks, Ri))
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
problem, domain=dedalus_setup.governing_equation(flag)
F_sin=1
Ra_ratio=1
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


ts = de.timesteppers.RK443

solver =  problem.build_solver(ts)
x = domain.grid(0)
z = domain.grid(1)
u = solver.state['u']
w = solver.state['w']
S = solver.state['S']
p = solver.state['p']
T = solver.state['T']

gshape = domain.dist.grid_layout.global_shape(scales=1)
slices = domain.dist.grid_layout.slices(scales=1)
rand = np.random.RandomState(seed=23)
noise = rand.standard_normal(gshape)[slices]

A_elevator=1
A_noise=0
A_shear=1
w['g'] = A_elevator*np.sin(k_opt*x) + A_noise*noise
u['g'] = A_noise*noise + A_shear*F_sin/ks**2*np.sin(ks*z)
S['g'] = -1/Ra_ratio*(k_opt**2+1/k_opt**2)*A_elevator*np.sin(k_opt*x) + A_noise*noise
T['g'] = -1/(k_opt**2)*A_elevator*np.sin(k_opt*x) + A_noise*noise
p['g'] = A_noise*noise

solver.stop_sim_time = 1
solver.stop_wall_time = np.inf
solver.stop_iteration = np.inf

#initial_dt = 0.02*Lx/nx
initial_dt=0.2*Lx/nx/(F_sin/ks**2)
cfl = flow_tools.CFL(solver,initial_dt,safety=0.8,max_change=1,cadence=8)
cfl.add_velocities(('u','w'))

analysis = solver.evaluator.add_file_handler(flag.name,sim_dt=0.1)
analysis.add_task('S')
analysis.add_task('u')
analysis.add_task('w')
analysis.add_task('T')

logger.info('Starting loop')
start_time=time.time()
while solver.ok:
    dt = cfl.compute_dt()    
    solver.step(dt)
    if solver.iteration % 10 == 0:
        logger.info('Iteration: %i, Time: %e, dt: %e' %(solver.iteration, solver.sim_time, dt))

end_time = time.time()


# Print statistics
logger.info('Run time: %f' %(end_time-start_time))
logger.info('Iterations: %i' %solver.iteration)
#logger.info('Sim end time: %f' %solver.sim_time)

logger.info('Run time: %f cpu-hr' %((end_time-start_time)/60/60*domain.dist.comm_cart.size))

#merge process data
post.merge_process_files(flag.name,cleanup=True)
flag.print_file()

#merge time snapshot
#set_paths = list(pathlib.Path("IFSC_2D_without_shear").glob("IFSC_2D_without_shear_s*.h5"))
#post.merge_sets("IFSC_2D_without_shear/IFSC_2D_without_shear.h5",set_paths,cleanup=True) 

#print the information for post-processing
#print(subprocess.check_output("find IFSC_2D_without_shear", shell=True).decode())
