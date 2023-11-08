import numpy as np
from mpi4py import MPI
import time
import pathlib
from dedalus import public as de
from dedalus.extras import flow_tools
from dedalus.tools import post

import logging
logger = logging.getLogger(__name__)

class flag:
    pass
# Parameters
flag=flag()

# Parameters
flag.Lx, flag.Ly, flag.Lz = (2*np.pi/10, 2*np.pi/10, 1.) #domain size
flag.Ra_T = 6*1e4 #Rayleigh number
flag.Pr=1
flag.Nx, flag.Ny, flag.Nz = (128,128,128) 

#parameter to control simulation and storage time
flag.initial_dt=0.001 #the initial time step
flag.stop_sim_time=10 #The simulation time to stop
flag.post_store_dt=0.01 #The time step to store the data
flag.k_elevator=10
flag.k_elevator_y=10

#paramter for the initial guess
flag.A_noise=0 #random noise magnitude in the initial condition
flag.restart_t0=1 #if 1, the simulation time will start from zero. Otherwise, will continue the previous one 
k_perp_square=flag.k_elevator**2+flag.k_elevator_y**2
flag.A_elevator=np.sqrt(flag.Ra_T/2/k_perp_square-k_perp_square/2)

# Create bases and domain
x_basis = de.Fourier('x', flag.Nx, interval=(0, flag.Lx), dealias=3/2)
y_basis = de.Fourier('y', flag.Ny, interval=(0, flag.Ly), dealias=3/2)
z_basis = de.Fourier('z', flag.Nz, interval=(0, flag.Lz), dealias=3/2)
domain = de.Domain([x_basis, y_basis, z_basis], grid_dtype=np.float64)

# 3D Boussinesq hydrodynamics
problem = de.IVP(domain, variables=['p','T','u','v','w'])
problem.parameters['Ra'] = flag.Ra_T
problem.parameters['Pr'] = flag.Pr
problem.parameters['vol'] = flag.Lx*flag.Ly*flag.Lz
problem.add_equation("dt(u) + dx(p) - Pr*(dx(dx(u))+dy(dy(u))+dz(dz(u)))= -u*dx(u)-v*dy(u)-w*dz(u)")
problem.add_equation("dt(v) + dy(p) - Pr*(dx(dx(v))+dy(dy(v))+dz(dz(v)))= -u*dx(v)-v*dy(v)-w*dz(v)")
problem.add_equation("dt(w) + dz(p) - Pr*(dx(dx(w))+dy(dy(w))+dz(dz(w))) - Pr*Ra*T= -u*dx(w)-v*dy(w)-w*dz(w)")
problem.add_equation("dx(u) + dy(v) + dz(w) = 0", condition="(nx!=0) or (ny!=0) or (nz!=0)")
problem.add_equation("p=0", condition="(nx==0) and (ny==0) and (nz==0)")
problem.add_equation("dt(T) - (dx(dx(T))+dy(dy(T))+dz(dz(T)))-w  = - w*integ(w*T)/vol -u*dx(T)-v*dy(T)-w*dz(T)",condition="(nx!=0) or (ny!=0) or (nz!=0)")
problem.add_equation("T=0", condition="(nx==0) and (ny==0) and (nz==0)")

# Build solver
solver = problem.build_solver(de.timesteppers.RK443)
logger.info('Solver built')

if not pathlib.Path('restart.h5').exists():

    print('Set up initial condition!')
    # Initial conditions setting up analytical elevator mode
    
    # Random perturbations, initialized globally for same results in parallel
    gshape = domain.dist.grid_layout.global_shape(scales=1)
    slices = domain.dist.grid_layout.slices(scales=1)
    rand = np.random.RandomState(seed=42)
    noise = rand.standard_normal(gshape)[slices]

    u0=flag.A_noise*noise
    v0=flag.A_noise*noise
    w0=flag.A_noise*noise
    T0=flag.A_noise*noise
    p0=flag.A_noise*noise

    x, y, z = domain.all_grids()
    w0=w0 + flag.A_elevator*np.real(np.exp(1j*(flag.k_elevator*x+flag.k_elevator_y*y)))                           
    T0=T0 + k_perp_square/flag.Ra_T*flag.A_elevator*np.real(np.exp(1j*(flag.k_elevator*x+flag.k_elevator_y*y)))
    
    u = solver.state['u']
    v = solver.state['v']
    w = solver.state['w']
    p = solver.state['p']
    T = solver.state['T']
    
    u['g'] = u0
    v['g'] = v0
    w['g'] = w0
    p['g'] = p0
    T['g'] = T0
    
    # Timestepping and output
    dt = flag.initial_dt
    stop_sim_time = flag.stop_sim_time
    fh_mode = 'overwrite'

else:
    # Restart
    print('Restart')
    write, last_dt = solver.load_state('restart.h5', -1)

    # Timestepping and output
    dt = last_dt
    stop_sim_time = flag.stop_sim_time
    fh_mode = 'append'
    if flag.restart_t0:
        solver.sim_time=0
        fh_mode='overwrite'

# Integration parameters
solver.stop_sim_time = flag.stop_sim_time

# Analysis
analysis = solver.evaluator.add_file_handler('analysis', sim_dt=flag.post_store_dt)
analysis.add_system(solver.state)

# CFL
CFL = flow_tools.CFL(solver, initial_dt=flag.initial_dt, cadence=10, safety=0.5,
                     max_change=1.5, min_change=0.5, max_dt=0.125, threshold=0.05)
CFL.add_velocities(('u', 'v', 'w'))

# Flow properties
flow = flow_tools.GlobalFlowProperty(solver, cadence=10)
flow.add_property("integ(w*T)/vol-1", name='dy_T_mean_q')
#flow_out.add_property('w*b',name='wb')
                            

           
def print_screen(flag,logger):
    #print the flag onto the screen
    flag_attrs=vars(flag)
    #print(', '.join("%s: %s, \n" % item for item in flag_attrs.items()))
    logger.info(', Attributes: Value,\n,')
    logger.info(', '.join("%s: %s, \n" % item for item in flag_attrs.items()))

def print_file(flag):
    #print the flag onto file
    flag_text=open('./analysis'+'/flag.txt','w+')
    flag_attrs=vars(flag)
    print(', Attributes: 123,\n ------\n-------\n------',file=flag_text)
    print(', test: 123,',file=flag_text)
    print(', '+', '.join("%s: %s, \n" % item for item in flag_attrs.items()),file=flag_text)
    flag_text.close()
    
        
# Main loop
try:
    logger.info('Starting loop')
    print_screen(flag,logger)
    print_file(flag)
    while solver.proceed:
        dt = CFL.compute_dt()
        dt = solver.step(dt)
        if (solver.iteration-1) % 1000 == 0:
            logger.info('Iteration: %i, Time: %e, dt: %e' %(solver.iteration, solver.sim_time, dt))
            logger.info('dy_T_mean_q = %f' %flow.max('dy_T_mean_q'))
            logger.info('Nu = %f' %-1/flow.max('dy_T_mean_q'))

    #add check point, only store one snapshot
    checkpoint=solver.evaluator.add_file_handler('checkpoint')
    checkpoint.add_system(solver.state)
    end_world_time = solver.get_world_time()
    end_wall_time = end_world_time - solver.start_time
    solver.evaluator.evaluate_handlers([checkpoint], timestep = flag.initial_dt, sim_time = solver.sim_time, world_time=end_world_time, wall_time=end_wall_time, iteration=solver.iteration)
    post.merge_process_files('checkpoint',cleanup=True)

except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    solver.log_stats()