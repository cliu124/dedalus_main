# -*- coding: utf-8 -*-
"""
Created on Sun Nov 14 20:30:28 2021

@author: Owner
"""



import time
import numpy as np
import matplotlib.pyplot as plt

from dedalus import public as de
import h5py
import logging
logger = logging.getLogger(__name__)



# Parameters
Nz = 1028
n = 3.0
ncc_cutoff = 1e-6
tolerance = 1e-12

#parameters
Ra_T=6000
kx=0.48*Ra_T**0.4
ky=kx
Ra_S2T=8000
tau=0.1
dy_T_mean=-1
dy_S_mean=-1

# Build domain
z_basis = de.Chebyshev('z', Nz, interval=(0, 1), dealias=2)
domain = de.Domain([z_basis], np.float64)

problem = de.NLBVP(domain, variables=[ 'w_hat','p_hat','T_hat','d_T_hat', \
        'S_hat','d_S_hat','T_0','d_T_0','S_0','d_S_0'])
problem.parameters['Ra_T'] = Ra_T
problem.parameters['Ra_S2T'] = Ra_S2T
problem.parameters['tau']=tau
problem.parameters['dy_T_mean']=dy_T_mean
problem.parameters['dy_S_mean']=dy_S_mean
problem.parameters['kx']=kx
problem.parameters['ky']=ky

problem.add_equation('dz(w_hat)-(-(kx*kx+ky*ky)*p_hat)=0')
problem.add_equation('dz(p_hat)-(-w_hat+Ra_T*T_hat-Ra_S2T*S_hat)=0')
problem.add_equation('dz(T_hat)-d_T_hat=0')
problem.add_equation('dz(d_T_hat)-(w_hat*dy_T_mean+(kx*kx+ky*ky)*T_hat)=w_hat*d_T_0')
problem.add_equation('dz(S_hat)-d_S_hat=0')
problem.add_equation('dz(d_S_hat)-1/tau*w_hat*dy_S_mean-(kx*kx+ky*ky)*S_hat=1/tau*(w_hat*d_S_0)')   
problem.add_equation('dz(T_0)-d_T_0=0')
problem.add_equation('dz(d_T_0)=-2*(kx*kx+ky*ky)*p_hat*T_hat+2*w_hat*d_T_hat')
problem.add_equation('dz(S_0)-d_S_0=0')
problem.add_equation('dz(d_S_0)=1/tau*(-2*(kx*kx+ky*ky)*p_hat*S_hat+2*w_hat*d_S_hat)')

problem.add_bc("left(w_hat) = 0")
problem.add_bc("right(w_hat) = 0")
problem.add_bc("left(T_hat) = 0")
problem.add_bc("right(T_hat) = 0")
problem.add_bc("left(S_hat) = 0")
problem.add_bc("right(S_hat) = 0")
problem.add_bc("left(T_0) = 0")
problem.add_bc("right(T_0) = 0")
problem.add_bc("left(S_0) = 0")
problem.add_bc("right(S_0) = 0")

# Setup initial guess
solver = problem.build_solver()
z = domain.grid(0)
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

W0=Ra_T;
w_hat['g'] = W0*np.sin(np.pi*z)
p_hat['g'] = W0*np.pi*np.cos(np.pi*z)/(-(kx*kx+ky*ky));
T_hat['g'] = 1/Ra_T*W0*np.sin(np.pi*z)
d_T_hat['g'] = 1/Ra_T*W0*np.pi*np.cos(np.pi*z)
S_hat['g'] = -1/Ra_S2T*W0*np.sin(np.pi*z)
d_S_hat['g'] = -1/Ra_S2T*W0*np.pi*np.cos(np.pi*z);
T_0['g'] = 0
d_T_0['g'] = 0
S_0['g'] = 0
d_S_0['g'] = 0




# Iterations
pert = solver.perturbations.data
pert.fill(1+tolerance)
start_time = time.time()
while np.sum(np.abs(pert)) > tolerance:
    solver.newton_iteration()
    logger.info('Perturbation norm: {}'.format(np.sum(np.abs(pert))))
    #logger.info('R iterate: {}'.format(R['g'][0]))
end_time = time.time()

hf = h5py.File('data.h5', 'w')
print(solver.state['w_hat']['g'])
print(solver.state['T_hat']['g'])
print(solver.state['S_hat']['g'])
print(solver.state['T_0']['g'])
print(solver.state['S_0']['g'])

hf.create_dataset('w_hat', data=solver.state['w_hat']['g'])
hf.create_dataset('p_hat', data=solver.state['p_hat']['g'])
hf.create_dataset('T_hat', data=solver.state['T_hat']['g'])
hf.create_dataset('d_T_hat', data=solver.state['d_T_hat']['g'])
hf.create_dataset('S_hat', data=solver.state['S_hat']['g'])
hf.create_dataset('d_S_hat', data=solver.state['d_S_hat']['g'])
hf.create_dataset('T_0', data=solver.state['T_0']['g'])
hf.create_dataset('d_T_0', data=solver.state['d_T_0']['g'])
hf.create_dataset('S_0', data=solver.state['S_0']['g'])
hf.create_dataset('d_S_0', data=solver.state['d_S_0']['g'])

hf.close()