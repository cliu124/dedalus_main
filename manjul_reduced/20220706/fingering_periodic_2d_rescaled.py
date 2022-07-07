'''
Dedalus script to simulate rescaled salt-finger equations in 2D periodic domain.
um, wm, Tm, Sm: <u>, <w>, <T>, <S> (<.> is horizontally averaged variable)
up, wp, Tp, Sp: u', w', T', S'

Parameters: Prandtl, tau, Ra_SH, R_rho
dTz = Background T profile
dSz = Background S profile
'''
import numpy as np
from mpi4py import MPI
import time

from dedalus import public as de
from dedalus.extras import flow_tools

#from filter_field import filter_field

import logging
logger = logging.getLogger(__name__)

# Parameters
Lx, Lz = (1., 1.)
Prandtl = 7.
#Ra_S = 1e8
tau = 1/100.
R_rho = 40.
R_rho_inv = 1./R_rho
tau_inv = 1./tau
Pr_tau = Prandtl/tau
#Ra_SH = Ra_S*tau**4
Ra_SH = 1e5
#dTz = 0.
#dSz = 0.
k_opt = 1.
beta = 1.#/tau           # For vertical diffusion terms
alpha = 1.              # To switch on/off non-linear term in mean salinity (Sm)
A = 2.*np.pi/16.
A_inv = 1./A

# Create bases and domain
x_basis = de.Fourier('x', 256, interval=(0, Lx), dealias=3/2)
z_basis = de.Fourier('z', 256, interval=(0, Lz), dealias=3/2)
domain = de.Domain([x_basis, z_basis], grid_dtype=np.float64)

# 2D re-scaled salt-fingers. Suffix 'p' for perturbation (or prime) and 'm' for mean (horizontal mean)
#problem = de.IVP(domain, variables=['p','Sp','up','wp','Tp','Sm','um','wm','Tm'])
problem = de.IVP(domain, variables=['p','Sp','up','wp','Tp','Sm','um','Tm'])
problem.parameters['Pr'] = Prandtl
problem.parameters['tau'] = tau
problem.parameters['tau_inv'] = tau_inv
problem.parameters['Pr_tau'] = Pr_tau
problem.parameters['Ra_SH'] = Ra_SH                      # Re-scaled Ra in the vertical momentum
problem.parameters['R_rho_inv'] = R_rho_inv
problem.parameters['R_rho'] = R_rho
problem.parameters['beta'] = beta
problem.parameters['alpha'] = alpha
problem.parameters['A'] = A
problem.parameters['A_inv'] = A_inv
#problem.parameters['dTz'] = dTz
#problem.parameters['dSz'] = dSz
problem.parameters['Lz'] = Lz
problem.parameters['Lx'] = Lx

problem.substitutions['vol_avg(A)']   = 'integ(A)/(Lx*Lz)'       # Volume average of A
problem.substitutions['h_mean(A)']   = "integ(A,'x')/Lx"         # Horizontal mean of A
problem.substitutions['u_full']   = "up + um"                    # u_total = u' + <u>
#problem.substitutions['w_full']   = "wp + wm"                    # w_total = w' + <w>
problem.substitutions['w_full']   = "wp"                    # w_total = w' + <w>
problem.substitutions['uw']   = "up*wp"                          # For correlation
problem.substitutions['ww']   = "wp*wp"                          # For correlation
problem.substitutions['wT']   = "wp*Tp"                          # For correlation
problem.substitutions['wS']   = "wp*Sp"                          # For correlation
problem.substitutions['wTm']  = "wp*Tm"
problem.substitutions['wSm']  = "wp*Sm"

#problem.substitutions['uuf']   = "u_full*up"                          # For correlation
#problem.substitutions['uwf']   = "u_full*wp"                          # For correlation
#problem.substitutions['uTf']   = "u_full*Tp"                          # For correlation
#problem.substitutions['uSf']   = "u_full*Sp"                          # For correlation
#problem.substitutions['wTf']   = "wp*Tp"                          # For correlation
#problem.substitutions['wSf']   = "wp*Sp"                          # For correlation

# u' equation in the non-conservative form, only the spatial correlation is in the conservative form
problem.add_equation("dt(up) + Pr*(A_inv*dx(p) - A_inv**2*dx(dx(up)) - dz(dz(up))) = -(A_inv*u_full*dx(up) + w_full*dz(up) + wp*dz(um)) + dz(h_mean(uw))", condition="(nx!=0)")
# w' equation in the non-conservative form
problem.add_equation("dt(wp)+Pr*(dz(p)-A_inv**2*dx(dx(wp))-dz(dz(wp))-Ra_SH*(Tp - R_rho_inv*Sp)) = -(A_inv*u_full*dx(wp)+w_full*dz(wp))+dz(h_mean(ww))",condition="(nx!=0)")
# <u> equation in the non-conservative form. Equation for <w> is not required
problem.add_equation("dt(um) - Pr*dz(dz(um)) = - dz(h_mean(uw))",condition="(nx ==0)")# and (nz==0)")
# perturbation continuity equation
problem.add_equation("dx(up) + dz(wp) = 0", condition="(nx != 0)")# or (nz != 0)")
# pressure gauge condition
problem.add_equation("p = 0", condition="(nx == 0)")# and (nz == 0)")
# T' equation in the non-conservative form
problem.add_equation("dt(Tp) - (A_inv**2*dx(dx(Tp)) + dz(dz(Tp)) - wp) = -(A_inv*u_full*dx(Tp) + w_full*dz(Tp) + wp*dz(Tm)) + dz(h_mean(wT))",condition="(nx!=0)")
# S' equation in the non-conservative form
problem.add_equation("dt(Sp) - tau*(A_inv**2*dx(dx(Sp)) - dz(dz(Sp))) + wp = -(A_inv*u_full*dx(Sp) + w_full*dz(Sp) + wp*dz(Sm)) + dz(h_mean(wS))",condition="(nx!=0)")
# <T> equation
problem.add_equation("dt(Tm) - dz(dz(Tm)) = - dz(h_mean(wT))", condition="(nx==0)")# and (nz==0)")
# <S> equation, multiplied by 'beta' to increase the diffusion
problem.add_equation("dt(Sm) - tau*dz(dz(Sm)) = - dz(h_mean(wS))",condition="(nx==0)")# and (nz==0)")

# Setting up fluctuations to have zero horizontal mean
problem.add_equation("Tp=0",condition="(nx==0)")# and (nz==0)")
problem.add_equation("Sp=0",condition="(nx==0)")# and (nz==0)")
problem.add_equation("up=0",condition="(nx==0)")# and (nz==0)")
problem.add_equation("wp=0",condition="(nx==0)")# and (nz==0)")
#problem.add_equation("wm=0",condition="(nx==0)")# and (nz==0)")

problem.add_equation("Tm=0",condition="(nx!=0)")# or (nz!=0)")
problem.add_equation("Sm=0",condition="(nx!=0)")# or (nz!=0)")
problem.add_equation("um=0",condition="(nx!=0)")# or (nz!=0)")
#problem.add_equation("wm=0",condition="(nx!=0)")# or (nz!=0)")

# Build solver
#solver = problem.build_solver(de.timesteppers.RK443)
solver = problem.build_solver(de.timesteppers.SBDF2)
#solver = problem.build_solver(de.timesteppers.RK222)
logger.info('Solver built')

# Initial conditions
x = domain.grid(0)
z = domain.grid(1)
up = solver.state['up']
um = solver.state['um']
wp = solver.state['wp']
#wm = solver.state['wm']
Tp = solver.state['Tp']
Tm = solver.state['Tm']
Sp = solver.state['Sp']
Sm = solver.state['Sm']
p = solver.state['p']

a = 0.05
sigma = 0.2
flow = 0.5
amp = -0.1

# Random perturbations, initialized globally for same results in parallel
gshape = domain.dist.grid_layout.global_shape(scales=1)
slices = domain.dist.grid_layout.slices(scales=1)
rand = np.random.RandomState(seed=42)
noise = rand.standard_normal(gshape)[slices]

pert1 =  1e-3 * noise#*np.sin(np.pi*z/Lz)
pert =  1e-2 * np.sin(k_opt*x)
#p['g'] = 0
up['g'] = pert1
wp['g'] = pert1
Tp['g'] = pert1
#T['g'] = -0.5*(1+np.tanh(z/a))
Sp['g'] = pert1

# Initial timestep
dt = 0.0001

# Integration parameters
solver.stop_sim_time = 4000
solver.stop_wall_time = np.inf
solver.stop_iteration = np.inf

# Analysis
snapshots = solver.evaluator.add_file_handler('snapshots', iter=1, max_writes=400)
snapshots.add_system(solver.state, scales=4)

profiles = solver.evaluator.add_file_handler('profiles', sim_dt=1., max_writes=np.inf)
profiles.add_task("h_mean(sqrt(Tm*Tm))",name='Tm_rms')
profiles.add_task("h_mean(sqrt(Sm*Sm))",name='Sm_rms')

timeseries = solver.evaluator.add_file_handler('timeseries',iter=1, max_writes=np.inf)
timeseries.add_task("vol_avg(sqrt(up*up))",name='urms')
timeseries.add_task("vol_avg(sqrt(wp*wp))",name='wrms')
timeseries.add_task("vol_avg(sqrt(Tp*Tp))",name='Trms') 
timeseries.add_task("vol_avg(sqrt(Sp*Sp))",name='Srms')
timeseries.add_task("h_mean(wp*Sp)",name='wSflux')
timeseries.add_task("h_mean(wp*Tp)",name='wTflux')

# CFL
CFL = flow_tools.CFL(solver, initial_dt=dt, cadence=1, safety=0.8,
                     max_change=1.5, min_change=0.5, max_dt=3.5, threshold=0.05)
#CFL.add_velocities(('u_full', 'w_full'))
CFL.add_velocities(('up', 'wp'))

# Flow properties
flow = flow_tools.GlobalFlowProperty(solver, cadence=1)
flow.add_property("0.5*(up*up + wp*wp)", name="Ekin")
flow.add_property("up", name="up")
flow.add_property("wp", name="wp")
flow.add_property("Tp", name="Tp")
flow.add_property("Sp", name="Sp")

# Main loop
try:
    logger.info('Starting loop')
    start_time = time.time()
    while solver.ok:
        dt = CFL.compute_dt()
        dt = solver.step(dt)#, trim=True)
        if (solver.iteration-1) % 10 == 0:
            logger.info('Iteration: %i, Time: %e, dt: %e, Kinetic Energy: %e, u max: %e, w max: %e' %(solver.iteration, solver.sim_time, dt, flow.volume_average('Ekin'), flow.max('up'), flow.max('wp')))
        #logger.info('urms = %f' %flow.max('urms'))
except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    end_time = time.time()
    logger.info('Iterations: %i' %solver.iteration)
    logger.info('Sim end time: %f' %solver.sim_time)
    logger.info('Run time: %.2f sec' %(end_time-start_time))
    logger.info('Run time: %f cpu-hr' %((end_time-start_time)/60/60*domain.dist.comm_cart.size))
