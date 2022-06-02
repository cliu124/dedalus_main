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

from filter_field import filter_field

import logging
logger = logging.getLogger(__name__)

# Parameters
Lx, Lz = (10., 10.)
Prandtl = 7.
Ra_S = 1e7
tau = 1/10.
R_rho = 2.
R_rho_inv = 1./R_rho
tau_inv = 1./tau
Pr_tau = Prandtl/tau
Ra_SH = Ra_S*tau**4
dTz = 1.
dSz = 1.

# Create bases and domain
x_basis = de.Fourier('x', 64, interval=(0, Lx), dealias=3/2)
z_basis = de.Fourier('z', 64, interval=(0, Lz), dealias=3/2)
domain = de.Domain([x_basis, z_basis], grid_dtype=np.float64)

# 2D re-scaled salt-fingers. Suffix 'p' for perturbation (or prime) and 'm' for mean (horizontal mean)
problem = de.IVP(domain, variables=['p','Sp','up','wp','Tp','Sm','um','wm','Tm'])
problem.parameters['Pr'] = Prandtl
problem.parameters['tau'] = tau
problem.parameters['tau_inv'] = tau_inv
problem.parameters['Pr_tau'] = Pr_tau
problem.parameters['Ra_SH'] = Ra_SH                      # Re-scaled Ra in the vertical momentum
problem.parameters['R_rho_inv'] = R_rho_inv
problem.parameters['R_rho'] = R_rho
problem.parameters['dTz'] = dTz
problem.parameters['dSz'] = dSz
problem.parameters['Lz'] = Lz
problem.parameters['Lx'] = Lx

problem.substitutions['vol_avg(A)']   = 'integ(A)/(Lx*Lz)'       # Volume average of A
problem.substitutions['h_mean(A)']   = "integ(A,'x')/Lx"         # Horizontal mean of A
problem.substitutions['u_full']   = "up + um"                    # u_total = u' + <u>
problem.substitutions['w_full']   = "wp + wm"                    # w_total = w' + <w>
problem.substitutions['uw']   = "up*wp"                          # For correlation
problem.substitutions['ww']   = "wp*wp"                          # For correlation
problem.substitutions['wT']   = "wp*Tp"                          # For correlation
problem.substitutions['wS']   = "wp*Sp"                          # For correlation

problem.add_equation("dt(up) + Pr_tau*(dx(p - h_mean(p)) - dx(dx(up)) + tau**2*dz(dz(up))) = -(u_full*dx(up) + w_full*dz(up)) + dz(h_mean(uw))")
problem.add_equation("dt(wp) + Pr_tau*(tau**2*dx(p - h_mean(p)) - dx(dx(wp)) + tau**2*dz(dz(wp))) - Ra_SH*(R_rho*Tm - Sm) = -(u_full*dx(wp) + w_full*dz(wp)) + dz(h_mean(ww))")
problem.add_equation("dt(um) - tau*Pr*dz(dz(um)) = - dz(h_mean(uw))")
problem.add_equation("Pr_tau*(dz(h_mean(p)) - Ra_SH*(R_rho*Tm - Sm)) = - dz(h_mean(ww))")
problem.add_equation("dx(u_full) + dz(w_full) = 0", condition="(nx != 0) and (nz != 0)")
problem.add_equation("dx(um) + dz(wm) = 0", condition="(nx != 0) and (nz != 0)")
problem.add_equation("p = 0", condition="(nx == 0) or (nz == 0)")
problem.add_equation("dt(Tp) - tau_inv*(dx(dx(Tp)) + tau**2*dz(dz(Tp)) + wp*(1 + dTz)) = -(u_full*dx(Tp) + w_full*dz(Tp)) + dz(h_mean(wT))")
problem.add_equation("dt(Sp) - dx(dx(Sp)) + tau**2*dz(dz(Sp)) + wp*(1 + dSz) = -(u_full*dx(Sp) + w_full*dz(Sp)) + dz(h_mean(wS))")
problem.add_equation("dt(Tm) - tau*dz(dz(Tm)) = - tau*dz(h_mean(wT))")
problem.add_equation("dt(Sm) - tau**2*dz(dz(Sm)) = - dz(h_mean(wS))")

# Build solver
solver = problem.build_solver(de.timesteppers.RK443)
logger.info('Solver built')

# Initial conditions
x = domain.grid(0)
z = domain.grid(1)
up = solver.state['up']
um = solver.state['um']
wp = solver.state['wp']
wm = solver.state['wm']
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

pert =  1e-3 * noise#*np.sin(np.pi*z/Lz)
#p['g'] = 0
up['g'] = pert
wp['g'] = pert
Tp['g'] = pert
#T['g'] = -0.5*(1+np.tanh(z/a))
Sp['g'] = pert

# Initial timestep
dt = 0.001

# Integration parameters
solver.stop_sim_time = 4000
solver.stop_wall_time = np.inf
solver.stop_iteration = np.inf

# Analysis
snapshots = solver.evaluator.add_file_handler('snapshots', iter=100, max_writes=400)
snapshots.add_system(solver.state, scales=4)

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
CFL.add_velocities(('u_full', 'w_full'))

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
            logger.info('Iteration: %i, Time: %e, dt: %e, Kinetic Energy: %e, u max: %e, w max: %e' %(solver.iteration, solver.sim_time, dt, flow.volume_average('Ekin'), flow.max('u'), flow.max('w')))
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
