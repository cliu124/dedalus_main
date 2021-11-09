"""
Dedalus script for the Lane-Emden equation.

This is a 1D script and should be ran serially.  It should converge within
roughly a dozen iterations, and should take under a minute to run.

In astrophysics, the Lane–Emden equation is a dimensionless form of Poisson's
equation for the gravitational potential of a Newtonian self-gravitating,
spherically symmetric, polytropic fluid [1].

It is usually written as:
    dr(dr(f)) + (2/r)*dr(f) + f**n = 0
    f(r=0) = 1
where n is the polytropic index, and the equation is solved over the interval
r=[0,R], where R is the n-dependent first zero of f(r). Although the equation
is second order, it is singular at r=0, and therefore only requires a single
outer boundary condition.

Following [2], we rescale the equation by defining r=R*x:
    dx(dx(f)) + (2/x)*dx(f) + (R**2)*(f**n) = 0
    f(x=0) = 1
    f(x=1) = 0
This is a nonlinear eigenvalue problem over the interval x=[0,1], with the
additional boundary condition fixing the eigenvalue R.

References:
    [1]: http://en.wikipedia.org/wiki/Lane–Emden_equation
    [2]: J. P. Boyd, "Chebyshev spectral methods and the Lane-Emden problem,"
         Numerical Mathematics Theory (2011).

"""

import time
import numpy as np
import matplotlib.pyplot as plt

from dedalus import public as de

import logging
logger = logging.getLogger(__name__)


# Parameters
Nx = 128
n = 3.0
ncc_cutoff = 1e-6
tolerance = 1e-12

# Build domain
x_basis = de.Chebyshev('x', Nx, interval=(0, 1), dealias=2)
domain = de.Domain([x_basis], np.float64)

# Setup problem
C=2
R2=10
problem = de.NLBVP(domain, variables=['theta_2_bar','d_theta_2_bar', 'theta_10','d_theta_10'], ncc_cutoff=ncc_cutoff)
#problem.meta['R']['x']['constant'] = True
problem.parameters['R2'] = R2
problem.parameters['C'] = C
#problem.parameters['n'] = n
problem.add_equation("d_theta_2_bar - dx(theta_2_bar) = 0")
problem.add_equation("dx(d_theta_2_bar) = theta_10*d_theta_10")
problem.add_equation('d_theta_10 - dx(theta_10) = 0')
problem.add_equation("dx(d_theta_10) + R2/C*theta_10 = 1/C*theta_10*d_theta_2_bar")
problem.add_bc("left(theta_2_bar) = 0")
problem.add_bc("right(theta_2_bar) = 0")
problem.add_bc("left(theta_10) = 0")
problem.add_bc("right(theta_10) = 0")

# Setup initial guess
solver = problem.build_solver()
x = domain.grid(0)
theta_10 = solver.state['theta_10']
d_theta_10 = solver.state['d_theta_10']
theta_2_bar = solver.state['theta_2_bar']
d_theta_2_bar = solver.state['d_theta_2_bar']

if C==3:
    ####setup initial guess. This is the high R_2 approximation in
    ####Blennerhassett & Bassom (1994)
    A2=1/24*R2**2-R2
    theta_10['g']=R2/2/np.sqrt(3)*(-1+np.tanh(1/12*R2*x)+np.tanh(1/12*R2*(1-x)))
    d_theta_2_bar['g']=(1/2*theta_10['g']**2-A2)
    theta_10.differentiate('x', out = d_theta_10)
    d_theta_2_bar.integrate('x', out = theta_2_bar)

elif C==2:
    ####%setup initial guess. This is the high R_2 approximation in
    ####%Lewis, Rees, Bassom (1997)
    A2=1/16*R2**2-R2
    theta_10['g']=R2/2/np.sqrt(2)*(-1+np.tanh(1/8*R2*x)+np.tanh(1/8*R2*(1-x)))
    d_theta_2_bar['g']=(1/2*theta_10['g']**2-A2)
    theta_10.differentiate('x', out = d_theta_10)
    d_theta_2_bar.integrate('x', out = theta_2_bar)


#f['g'] = np.cos(np.pi/2 * x)*0.9


# analysis = solver.evaluator.add_file_handler('analysis')
# analysis.add_task('theta_10',layout='g',name='theta_10')
# analysis.add_task('theta_2_bar',layout='g',name='theta_2_bar')
            
# Iterations
pert = solver.perturbations.data
pert.fill(1+tolerance)
start_time = time.time()
while np.sum(np.abs(pert)) > tolerance:
    solver.newton_iteration()
    logger.info('Perturbation norm: {}'.format(np.sum(np.abs(pert))))
    #logger.info('R iterate: {}'.format(R['g'][0]))
end_time = time.time()

print(theta_10['g'])
print(theta_2_bar['g'])

# # Compare to reference solutions from Boyd
# R_ref = {0.0: np.sqrt(6),
#          0.5: 2.752698054065,
#          1.0: np.pi,
#          1.5: 3.65375373621912608,
#          2.0: 4.3528745959461246769735700,
#          2.5: 5.355275459010779,
#          3.0: 6.896848619376960375454528,
#          3.25: 8.018937527,
#          3.5: 9.535805344244850444,
#          4.0: 14.971546348838095097611066,
#          4.5: 31.836463244694285264}
# logger.info('-'*20)
# logger.info('Iterations: {}'.format(solver.iteration))
# logger.info('Run time: %.2f sec' %(end_time-start_time))
# logger.info('Final R iteration: {}'.format(R['g'][0]))
# if n in R_ref:
#     logger.info('Error vs reference: {}'.format(R['g'][0]-R_ref[n]))


